///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2026  Washington State Department of Transportation
//                        Bridge and Structures Office
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the Alternate Route Open Source License as
// published by the Washington State Department of Transportation,
// Bridge and Structures Office.
//
// This program is distributed in the hope that it will be useful, but
// distribution is AS IS, WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
// the Alternate Route Open Source License for more details.
//
// You should have received a copy of the Alternate Route Open Source
// License along with this program; if not, write to the Washington
// State Department of Transportation, Bridge and Structures Office,
// P.O. Box  47340, Olympia, WA 98503, USA or e-mail
// Bridge_Support@wsdot.wa.gov
///////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "RoadwaySection.h"
#include "GirderSection.h"

#include <IFace\Project.h>
#include <IFace/Alignment.h>

#include <iomanip>
#include <map>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   const double inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);

   // The section of the roadway at a station. Offsets are from the alignment, positive to the right
   struct StationSection
   {
      Float64 station;
      std::vector<double> ridges;          // offsets of the ridge points, left to right; one is at the alignment (0)
      std::vector<double> slopes;          // rise to the right of the segments between the ridge points
      double left_slope, right_slope;      // rise to the right of the outer segments
      double elevation;                    // top of the deck at the alignment
      double from, to;                     // the offsets of the roadway measured

      size_t alignment_point() const { return std::find(ridges.begin(), ridges.end(), 0.0) - ridges.begin(); }

      // the elevation of the surface at offset o, relative to the alignment
      double surface(double o) const
      {
         size_t a = alignment_point();
         std::vector<double> z(ridges.size(), 0.0);
         for (size_t i = a; i + 1 < ridges.size(); i++)
            z[i + 1] = z[i] + slopes[i] * (ridges[i + 1] - ridges[i]);
         for (size_t i = a; 0 < i; i--)
            z[i - 1] = z[i] - slopes[i - 1] * (ridges[i] - ridges[i - 1]);
         if (o <= ridges.front())
            return z.front() + left_slope * (o - ridges.front());
         if (ridges.back() <= o)
            return z.back() + right_slope * (o - ridges.back());
         size_t i = 1;
         while (ridges[i] < o)
            i++;
         return z[i - 1] + slopes[i - 1] * (o - ridges[i - 1]);
      }

      // the section with a ridge point added at offset o (where the surface doesn't change slope), for a section whose
      // ridge point faded away (e.g. in a superelevation transition)
      StationSection with_ridge(double o) const
      {
         StationSection s(*this);
         s.ridges.insert(std::upper_bound(s.ridges.begin(), s.ridges.end(), o), o);
         s.slopes.clear();
         for (size_t i = 0; i + 1 < s.ridges.size(); i++)
            s.slopes.push_back((surface(s.ridges[i + 1]) - surface(s.ridges[i])) / (s.ridges[i + 1] - s.ridges[i]));
         return s;
      }

      // the same ridge points, and surfaces within 1/16 in of each other over the roadways measured
      bool same(const StationSection& other) const
      {
         if (ridges.size() != other.ridges.size() || alignment_point() != other.alignment_point())
            return false;
         std::vector<double> offsets{ from, to, other.from, other.to };
         offsets.insert(offsets.end(), ridges.begin(), ridges.end());
         offsets.insert(offsets.end(), other.ridges.begin(), other.ridges.end());
         return std::all_of(offsets.begin(), offsets.end(), [&](double o) {return IsEqual(surface(o), other.surface(o), inch / 16); });
      }
   };

   using Segments = std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>>;

   // The highest point of a cut at offset o (points are (offset, elevation)), or nullopt if the cut doesn't reach it
   std::optional<double> top_at(const Segments& segments, double o)
   {
      std::optional<double> top;
      for (const auto& [a, b] : segments)
      {
         if ((a.x() <= o) != (b.x() <= o))
            top = std::max(top.value_or(-DBL_MAX), a.y() + (o - a.x()) * (b.y() - a.y()) / (b.x() - a.x()));
      }
      return top;
   }

   // The polyline's elevation at o, extended beyond its ends along its end segments
   double elevation_at(const std::vector<Eigen::Vector2d>& line, double o)
   {
      size_t i = 1;
      while (i + 1 < line.size() && line[i].x() < o)
         i++;
      const auto& a = line[i - 1];
      const auto& b = line[i];
      return a.y() + (o - a.x()) * (b.y() - a.y()) / (b.x() - a.x());
   }

   std::optional<StationSection> measure(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, Float64 station)
   {
      auto cut = get_alignment_cut(pBroker, station);
      // the deck section in (offset, elevation), as the segments of the cut (the deck mesh may have haunches that touch
      // the slab, so it isn't joined into loops)
      Segments segments;
      for (const auto& [a, b] : CutMeshSegments(deck, cut.point, cut.axis))
         segments.emplace_back(Eigen::Vector2d(-a.x(), a.y()), Eigen::Vector2d(-b.x(), b.y()));
      if (segments.empty())
         return std::nullopt;
      double deck_from = DBL_MAX, deck_to = -DBL_MAX;
      for (const auto& [a, b] : segments)
      {
         deck_from = std::min({ deck_from, a.x(), b.x() });
         deck_to = std::max({ deck_to, a.x(), b.x() });
      }
      deck_from += inch; // clear of the edges
      deck_to -= inch;

      // the offsets the barriers cover (and an inch more): the deck there isn't roadway
      std::vector<std::pair<double, double>> covered;
      for (const auto& barrier : barriers)
      {
         for (const auto& barrier_loop : CutMesh(barrier.mesh, cut.point, cut.axis))
         {
            double from = DBL_MAX, to = -DBL_MAX;
            for (const auto& p : barrier_loop)
            {
               from = std::min(from, -p.x());
               to = std::max(to, -p.x());
            }
            covered.emplace_back(from - inch, to + inch);
         }
      }
      auto is_covered = [&covered](double o) {return std::any_of(covered.begin(), covered.end(), [o](const auto& c) {return c.first < o && o < c.second; }); };

      // the top of the deck at its corners, the edges of the barriers, and the ends
      std::vector<double> offsets{ deck_from, deck_to };
      for (const auto& [a, b] : segments)
      {
         offsets.push_back(a.x());
         offsets.push_back(b.x());
      }
      for (const auto& c : covered)
      {
         offsets.push_back(c.first);
         offsets.push_back(c.second);
      }
      std::sort(offsets.begin(), offsets.end());
      offsets.erase(std::unique(offsets.begin(), offsets.end(), [](double a, double b) {return IsEqual(a, b, 1e-6); }), offsets.end());

      // the uncovered stretches of the top of the deck; the roadway is the one at the alignment, else the longest
      std::vector<std::vector<Eigen::Vector2d>> stretches;
      std::vector<Eigen::Vector2d> stretch;
      for (size_t i = 0; i < offsets.size(); i++)
      {
         double o = offsets[i];
         bool in_deck = deck_from - 1e-9 <= o && o <= deck_to + 1e-9;
         auto top = in_deck ? top_at(segments, o) : std::nullopt;
         bool continues = !stretch.empty() && !is_covered((stretch.back().x() + o) / 2);
         if (!top || is_covered(o) || (!stretch.empty() && !continues))
         {
            if (1 < stretch.size())
               stretches.push_back(stretch);
            stretch.clear();
         }
         if (top && !is_covered(o))
            stretch.emplace_back(o, *top);
      }
      if (1 < stretch.size())
         stretches.push_back(stretch);
      if (stretches.empty())
         return std::nullopt;

      auto length = [](const std::vector<Eigen::Vector2d>& s) {return s.back().x() - s.front().x(); };
      auto roadway = std::max_element(stretches.begin(), stretches.end(), [&](const auto& a, const auto& b)
         {
            bool a_at = a.front().x() <= 0 && 0 <= a.back().x();
            bool b_at = b.front().x() <= 0 && 0 <= b.back().x();
            return a_at != b_at ? b_at : length(a) < length(b);
         });
      if (length(*roadway) < 36 * inch)
         return std::nullopt;

      // its corners are the ridge points, with one at the alignment (a corner within 1/2 in of it moves there)
      auto corners = SimplifyPolyline(*roadway, inch / 16);
      for (size_t i = 1; i + 1 < corners.size(); i++)
      {
         if (IsZero(corners[i].x(), inch / 2))
            corners[i].x() = 0.0;
      }
      StationSection section;
      section.station = station;
      section.from = corners.front().x();
      section.to = corners.back().x();
      section.ridges.push_back(0.0);
      for (size_t i = 1; i + 1 < corners.size(); i++)
      {
         if (corners[i].x() != 0.0)
            section.ridges.push_back(corners[i].x());
      }
      std::sort(section.ridges.begin(), section.ridges.end());

      const double h = 1e-6; // ridge points are corners or between corners, so the outer slopes are just outside them
      section.left_slope = (elevation_at(corners, section.ridges.front()) - elevation_at(corners, section.ridges.front() - h)) / h;
      section.right_slope = (elevation_at(corners, section.ridges.back() + h) - elevation_at(corners, section.ridges.back())) / h;
      for (size_t i = 0; i + 1 < section.ridges.size(); i++)
      {
         double a = section.ridges[i], b = section.ridges[i + 1];
         section.slopes.push_back((elevation_at(corners, b) - elevation_at(corners, a)) / (b - a));
      }
      section.elevation = elevation_at(corners, 0.0);
      return section;
   }

   std::string percent(double slope)
   {
      std::ostringstream os;
      os << std::fixed << std::setprecision(2) << 100 * slope << "%";
      return os.str();
   }

   std::string feet(double length)
   {
      std::ostringstream os;
      os << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(length, WBFL::Units::Measure::Feet) << " ft";
      return os.str();
   }

   std::string describe(const StationSection& s)
   {
      std::ostringstream os;
      os << "left " << percent(s.left_slope);
      for (size_t i = 0; i < s.ridges.size(); i++)
      {
         os << ", ridge at " << (s.ridges[i] == 0 ? std::string("the alignment") : feet(s.ridges[i]));
         if (i < s.slopes.size())
            os << ", " << percent(s.slopes[i]) << " for " << feet(s.ridges[i + 1] - s.ridges[i]);
      }
      os << ", right " << percent(s.right_slope) << " (rise to the right, measured from " << feet(s.from) << " to " << feet(s.to) << ")";
      return os.str();
   }

   // A profile fitted to elevations along the alignment: tangents and parabolic vertical curves with the grade
   // continuous (a quadratic spline). The pieces are between knots and are curved or straight
   struct ProfileFit
   {
      std::vector<double> knots; // the ends of the pieces, from the first station to the last
      std::vector<bool> curved;  // for each piece
      Eigen::VectorXd x;         // the elevation and grade at the first knot, then the rate of change of grade of each curved piece
      double max_residual = DBL_MAX;
      double sum_of_squares = DBL_MAX;

      // the elevation at s is row(s).dot(x); with grade_row, the grade
      Eigen::RowVectorXd row(double s, bool grade_row = false) const
      {
         Eigen::RowVectorXd r(x.size());
         r(0) = grade_row ? 0.0 : 1.0;
         r(1) = grade_row ? 1.0 : s - knots.front();
         Eigen::Index c = 2;
         for (size_t k = 0; k < curved.size(); k++)
         {
            if (!curved[k])
               continue;
            double a = knots[k], b = knots[k + 1];
            double u = std::clamp(s, a, b) - a; // distance into the piece
            double beyond = std::max(s - b, 0.0);
            r(c++) = grade_row ? u : u * u / 2 + (b - a) * beyond;
         }
         return r;
      }
      double elevation(double s) const { return row(s).dot(x); }
      double grade(double s) const { return row(s, true).dot(x); }
   };

   ProfileFit fit_pieces(const std::vector<Eigen::Vector2d>& samples, const std::vector<double>& knots, const std::vector<bool>& curved)
   {
      ProfileFit fit;
      fit.knots = knots;
      fit.curved = curved;
      fit.x = Eigen::VectorXd::Zero(2 + std::count(curved.begin(), curved.end(), true));
      if (samples.size() < (size_t)fit.x.size())
         return fit;
      Eigen::MatrixXd A(samples.size(), fit.x.size());
      Eigen::VectorXd z(samples.size());
      for (size_t i = 0; i < samples.size(); i++)
      {
         A.row(i) = fit.row(samples[i].x());
         z(i) = samples[i].y();
      }
      fit.x = A.colPivHouseholderQr().solve(z);
      Eigen::VectorXd residuals = A * fit.x - z;
      fit.max_residual = residuals.cwiseAbs().maxCoeff();
      fit.sum_of_squares = residuals.squaredNorm();
      return fit;
   }

   // The simplest profile within tolerance of the samples: a tangent, else one vertical curve, else vertical curves
   // between knots added at the samples one at a time (where the fit improves most); then pieces become tangents and
   // knots are removed where the fit stays within tolerance. The end pieces reach from and to the given stations
   ProfileFit fit_profile(const std::vector<Eigen::Vector2d>& samples, double from, double to, double tolerance)
   {
      double first = std::min(from, samples.front().x()), last = std::max(to, samples.back().x());
      auto best = fit_pieces(samples, { first, last }, { false });
      if (best.max_residual <= tolerance)
         return best;
      best = fit_pieces(samples, { first, last }, { true });

      while (tolerance < best.max_residual && best.knots.size() < samples.size() / 3 + 2)
      {
         ProfileFit next;
         for (const auto& p : samples)
         {
            double t = p.x();
            if (std::any_of(best.knots.begin(), best.knots.end(), [t](double k) {return IsEqual(k, t, 1e-6); }))
               continue;
            auto knots = best.knots;
            knots.insert(std::upper_bound(knots.begin(), knots.end(), t), t);
            auto f = fit_pieces(samples, knots, std::vector<bool>(knots.size() - 1, true));
            if (f.sum_of_squares < next.sum_of_squares)
               next = f;
         }
         if (best.sum_of_squares <= next.sum_of_squares)
            break;
         best = next;
      }

      double limit = std::max(tolerance, best.max_residual);
      for (bool changed = true; changed; )
      {
         changed = false;
         for (size_t k = 0; k < best.curved.size(); k++)
         {
            if (!best.curved[k])
               continue;
            auto curved = best.curved;
            curved[k] = false;
            auto f = fit_pieces(samples, best.knots, curved);
            if (f.max_residual <= limit)
            {
               best = f;
               changed = true;
            }
         }
         for (size_t j = 1; j + 1 < best.knots.size(); j++)
         {
            auto knots = best.knots;
            auto curved = best.curved;
            knots.erase(knots.begin() + j);
            curved[j - 1] = curved[j - 1] || curved[j];
            curved.erase(curved.begin() + j);
            auto f = fit_pieces(samples, knots, curved);
            if (f.max_residual <= limit)
            {
               best = f;
               changed = true;
               break;
            }
         }
      }
      return best;
   }

   // PGSuper's profile: from the first knot, with a symmetric vertical curve for each curved piece
   ProfileData2 to_profile_data(const ProfileFit& fit)
   {
      ProfileData2 profile;
      profile.Station = fit.knots.front();
      profile.Elevation = fit.elevation(profile.Station);
      profile.Grade = fit.grade(profile.Station);
      for (size_t k = 0; k < fit.curved.size(); k++)
      {
         if (!fit.curved[k])
            continue;
         VertCurveData curve;
         curve.PVIStation = (fit.knots[k] + fit.knots[k + 1]) / 2;
         curve.L1 = fit.knots[k + 1] - fit.knots[k];
         curve.L2 = 0; // symmetric
         curve.ExitGrade = fit.grade(fit.knots[k + 1]);
         profile.VertCurves.push_back(curve);
      }
      return profile;
   }

   std::string describe(const ProfileFit& fit)
   {
      std::ostringstream os;
      os << "elevation " << feet(fit.elevation(fit.knots.front())) << " at station " << feet(fit.knots.front());
      for (size_t k = 0; k < fit.curved.size(); k++)
      {
         double a = fit.knots[k], b = fit.knots[k + 1];
         if (fit.curved[k])
            os << ", a vertical curve from station " << feet(a) << " to " << feet(b) << " (" << percent(fit.grade(a)) << " to " << percent(fit.grade(b)) << ")";
         else
            os << ", a grade of " << percent(fit.grade(a)) << " from station " << feet(a) << " to " << feet(b);
      }
      return os.str();
   }
}

void import_roadway_section(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, const std::vector<Float64>& pier_stations)
{
   if (pier_stations.size() < 2)
      return;

   // the quarter points of each span
   std::vector<StationSection> sections;
   size_t not_measured = 0;
   for (size_t span = 0; span + 1 < pier_stations.size(); span++)
   {
      for (double fraction : { 0.25, 0.5, 0.75 })
      {
         if (auto section = measure(pBroker, deck, barriers, pier_stations[span] + fraction * (pier_stations[span + 1] - pier_stations[span])))
            sections.push_back(*section);
         else
            not_measured++;
      }
   }

   // and as near the ends of the bridge as the deck goes its full width (for superelevation transitions): at least 90%
   // of the median width measured at the quarter points (at a skewed end, a cut near it crosses part of the deck)
   std::vector<double> widths;
   for (const auto& section : sections)
      widths.push_back(section.to - section.from);
   if (!widths.empty())
   {
      std::nth_element(widths.begin(), widths.begin() + widths.size() / 2, widths.end());
      double full_width = 0.9 * widths[widths.size() / 2];
      auto near_end = [&](Float64 end, Float64 direction, Float64 span) -> std::optional<StationSection>
         {
            for (double distance : { 1.0, 2.0, 4.0, 8.0, 16.0 })
            {
               if (span / 4 <= distance * 12 * inch)
                  break;
               auto section = measure(pBroker, deck, barriers, end + direction * distance * 12 * inch);
               if (section && full_width <= section->to - section->from)
                  return section;
            }
            return std::nullopt;
         };
      if (auto section = near_end(pier_stations.front(), 1.0, pier_stations[1] - pier_stations[0]))
         sections.insert(sections.begin(), *section);
      if (auto section = near_end(pier_stations.back(), -1.0, pier_stations.back() - pier_stations[pier_stations.size() - 2]))
         sections.push_back(*section);
   }
   if (sections.empty())
   {
      WBFL::System::Logger::Info("The roadway section couldn't be measured from the top of the deck. The roadway section is not changed.");
      return;
   }

   // the number of ridge points and the alignment's must be the same at every station: the most common is used
   std::map<std::pair<size_t, size_t>, size_t> counts;
   for (const auto& s : sections)
      counts[{ s.ridges.size(), s.alignment_point() }]++;
   auto common = std::max_element(counts.begin(), counts.end(), [](const auto& a, const auto& b) {return a.second < b.second; })->first;
   auto is_common = [&common](const StationSection& s) {return std::make_pair(s.ridges.size(), s.alignment_point()) == common; };

   // A section with fewer ridge points gets the missing ones of the nearest common section (a ridge point more than a
   // foot from its own)
   for (auto& s : sections)
   {
      if (is_common(s) || common.first <= s.ridges.size())
         continue;
      auto nearest = std::min_element(sections.begin(), sections.end(), [&](const auto& a, const auto& b)
         {
            auto distance = [&](const StationSection& x) {return is_common(x) ? fabs(x.station - s.station) : DBL_MAX; };
            return distance(a) < distance(b);
         });
      auto augmented = s;
      for (double o : nearest->ridges)
      {
         if (std::none_of(augmented.ridges.begin(), augmented.ridges.end(), [o](double r) {return IsEqual(r, o, 12 * inch); }))
            augmented = augmented.with_ridge(o);
      }
      if (is_common(augmented))
         s = augmented;
   }

   std::vector<StationSection> different;
   std::vector<StationSection> used;
   for (const auto& s : sections)
      (is_common(s) ? used : different).push_back(s);

   // In a transition (the sections at the first two, or last two, stations differ), the first and last sections are
   // extended linearly to the ends of the bridge, as PGSuper interpolates between templates
   size_t extrapolated = 0;
   auto extrapolate = [&](const StationSection& a, const StationSection& b, Float64 station)
      {
         StationSection s(b);
         double f = (station - b.station) / (b.station - a.station);
         auto lerp = [f](double u, double v) {return v + f * (v - u); };
         s.station = station;
         s.left_slope = lerp(a.left_slope, b.left_slope);
         s.right_slope = lerp(a.right_slope, b.right_slope);
         for (size_t i = 0; i < s.ridges.size(); i++)
            s.ridges[i] = (b.ridges[i] == 0.0 ? 0.0 : lerp(a.ridges[i], b.ridges[i]));
         for (size_t i = 0; i < s.slopes.size(); i++)
            s.slopes[i] = lerp(a.slopes[i], b.slopes[i]);
         return s;
      };
   if (2 <= used.size() && pier_stations.front() < used.front().station && !used[0].same(used[1]))
   {
      used.insert(used.begin(), extrapolate(used[1], used[0], pier_stations.front()));
      extrapolated++;
   }
   if (2 <= used.size() && used.back().station < pier_stations.back() && !used[used.size() - 1].same(used[used.size() - 2]))
   {
      used.push_back(extrapolate(used[used.size() - 2], used.back(), pier_stations.back()));
      extrapolated++;
   }

   // A template where the section changes: at the first and last stations, and where the section isn't the one PGSuper
   // interpolates between the templates on either side (so a linear transition has a template at each end). One
   // template if the section is the same everywhere
   std::vector<const StationSection*> templates;
   bool all_same = std::all_of(used.begin(), used.end(), [&](const auto& s) {return s.same(used.front()); });
   for (size_t i = 0; i < used.size(); i++)
   {
      if (all_same && i != 0)
         break;
      if (i == 0 || i + 1 == used.size() || !extrapolate(*templates.back(), used[i + 1], used[i].station).same(used[i]))
         templates.push_back(&used[i]);
   }

   RoadwaySectionData data;
   data.slopeMeasure = RoadwaySectionData::RelativeToAlignmentPoint;
   data.NumberOfSegmentsPerSection = common.first + 1;
   data.AlignmentPointIdx = common.second + 1; // ridge points are numbered from 1
   data.ProfileGradePointIdx = data.AlignmentPointIdx;
   for (const auto* s : templates)
   {
      // slopes are measured away from the alignment point: the outer left segment and the segments left of the
      // alignment point are the negative of the rise to the right
      RoadwaySectionTemplate t;
      t.Station = s->station;
      t.LeftSlope = -s->left_slope;
      t.RightSlope = s->right_slope;
      for (size_t i = 0; i < s->slopes.size(); i++)
      {
         IndexType ridge_point = i + 1; // the segment starts at this ridge point
         RoadwaySegmentData segment;
         segment.Length = s->ridges[i + 1] - s->ridges[i];
         segment.Slope = (ridge_point < data.AlignmentPointIdx ? -1 : 1) * s->slopes[i];
         t.SegmentDataVec.push_back(segment);
      }
      data.RoadwaySectionTemplates.push_back(t);
   }

   GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
   pRoadwayData->SetRoadwaySectionData(data);

   std::ostringstream os;
   os << "Roadway section from the top of the deck at " << sections.size() << " stations (the quarter points of the spans and near the ends of the bridge): ";
   if (all_same)
      os << "the same at every station: " << describe(used.front()) << ".";
   else
   {
      os << templates.size() << " templates:";
      for (const auto* s : templates)
         os << " station " << feet(s->station) << ": " << describe(*s) << ";";
   }
   if (0 < extrapolated)
      os << " The sections are a transition at " << (extrapolated == 2 ? "both ends" : "one end") << " of the bridge, so the nearest section is extended linearly to the end pier.";
   if (0 < not_measured)
      os << " " << not_measured << " stations couldn't be measured.";
   for (const auto& s : different)
      os << " At station " << feet(s.station) << " the section has " << s.ridges.size() << " ridge points, not " << common.first << ", and is left out: " << describe(s) << ".";
   WBFL::System::Logger::Info(os.str());
}

void import_profile(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, const std::vector<Float64>& pier_stations, bool bReplace)
{
   if (pier_stations.size() < 2)
      return;

   // the roadway surface at the alignment (the top of the deck, or the roadway extended to the alignment where the deck
   // doesn't reach it) about every 10 ft, and at least 8 times in each span
   std::vector<Eigen::Vector2d> on_deck, extended;
   for (size_t span = 0; span + 1 < pier_stations.size(); span++)
   {
      double length = pier_stations[span + 1] - pier_stations[span];
      size_t n = std::max<size_t>(8, (size_t)ceil(length / (120 * inch)));
      for (size_t i = 0; i <= n; i++)
      {
         if (i == n && span + 2 < pier_stations.size())
            break; // the next span starts here
         double station = pier_stations[span] + length * i / n;
         // at the end piers, as near them as the deck goes (it may stop short of the CL pier)
         bool first = (span == 0 && i == 0);
         bool last = (span + 2 == pier_stations.size() && i == n);
         for (double distance : { 0.0, 1.0, 2.0, 4.0, 8.0 })
         {
            if (0 < distance && (!(first || last) || length / n <= distance * 12 * inch))
               break;
            double at = station + (first ? 1 : -1) * distance * 12 * inch;
            if (auto section = measure(pBroker, deck, barriers, at))
            {
               (section->from <= 0 && 0 <= section->to ? on_deck : extended).emplace_back(at, section->elevation);
               break;
            }
         }
      }
   }
   // the deck at the alignment where it's measured there at most stations (a skewed end can leave it off the deck)
   auto& samples = (extended.size() < on_deck.size() ? on_deck : extended);   if (samples.size() < 3)
   {
      WBFL::System::Logger::Info("The profile couldn't be measured from the top of the deck. The profile is not changed.");
      return;
   }

   // glitches in the deck (e.g. a cut that misses the top of the deck, which happens at a few stations of the skewed
   // PGSuper exports): a sample more than 1/2 in from the local profile is left out. The local profile is the parabola
   // through three of the eight samples nearest it (four on each side, or eight on one side at an end) that the most of
   // them agree with (least median of squares), so it holds when a few glitches are near each other
   std::vector<Eigen::Vector2d> glitches;
   std::vector<Eigen::Vector2d> kept;
   for (size_t i = 0; i < samples.size(); i++)
   {
      std::vector<Eigen::Vector2d> nearby;
      size_t from = (i < 4 ? 0 : i - 4);
      size_t to = std::min(from + 8, samples.size() - 1);
      from = (8 < to ? to - 8 : 0);
      for (size_t j = from; j <= to; j++)
      {
         if (j != i)
            nearby.push_back(samples[j] - samples[i]); // relative to the sample, which is at the origin
      }
      bool glitch = false;
      if (4 <= nearby.size())
      {
         size_t rank = nearby.size() / 2 + 1; // a majority of the nearby samples
         double best_rank_residual = DBL_MAX, best_at_sample = 0;
         std::vector<double> residuals(nearby.size());
         for (size_t a = 0; a < nearby.size(); a++)
         {
            for (size_t b = a + 1; b < nearby.size(); b++)
            {
               for (size_t c = b + 1; c < nearby.size(); c++)
               {
                  // the parabola through three samples (Lagrange form), at s
                  auto parabola = [&](double s)
                     {
                        const auto& p = nearby[a];
                        const auto& q = nearby[b];
                        const auto& r = nearby[c];
                        return p.y() * (s - q.x()) * (s - r.x()) / ((p.x() - q.x()) * (p.x() - r.x())) +
                               q.y() * (s - p.x()) * (s - r.x()) / ((q.x() - p.x()) * (q.x() - r.x())) +
                               r.y() * (s - p.x()) * (s - q.x()) / ((r.x() - p.x()) * (r.x() - q.x()));
                     };
                  for (size_t k = 0; k < nearby.size(); k++)
                     residuals[k] = fabs(nearby[k].y() - parabola(nearby[k].x()));
                  std::nth_element(residuals.begin(), residuals.begin() + rank - 1, residuals.end());
                  if (residuals[rank - 1] < best_rank_residual)
                  {
                     best_rank_residual = residuals[rank - 1];
                     best_at_sample = parabola(0.0);
                  }
               }
            }
         }
         glitch = !IsZero(best_at_sample, inch / 2);
      }
      (glitch ? glitches : kept).push_back(samples[i]);
   }
   if (kept.size() < 3)
   {
      WBFL::System::Logger::Info("The profile couldn't be measured from the top of the deck. The profile is not changed.");
      return;
   }

   // from end pier to end pier: the deck may stop short of them, and a vertical curve there continues to them
   auto fit = fit_profile(kept, pier_stations.front(), pier_stations.back(), inch / 8);

   std::ostringstream os;
   os << "Profile from the top of the deck at the alignment at " << kept.size() << " stations: " << describe(fit) << " (within " << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(fit.max_residual, WBFL::Units::Measure::Inch) << " in).";
   if (&samples == &extended)
      os << " The deck doesn't reach the alignment at most stations, so the roadway is extended to it.";
   for (const auto& g : glitches)
      os << " At station " << feet(g.x()) << " the top of the deck is out of line with its neighbors and is left out.";

   if (bReplace)
   {
      GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
      pRoadwayData->SetProfileData2(to_profile_data(fit));
   }
   else
   {
      // the alignment has a profile: the deck is a check of it (the profile grade is the roadway surface)
      GET_IFACE2(pBroker, IRoadway, pAlignment);
      double lowest = DBL_MAX, highest = -DBL_MAX;
      for (const auto& p : kept)
      {
         double difference = p.y() - pAlignment->GetElevation(p.x(), 0.0);
         lowest = std::min(lowest, difference);
         highest = std::max(highest, difference);
      }
      os << " The alignment has a profile, which is used.";
      if (inch / 8 < std::max(fabs(lowest), fabs(highest)))
         os << " The top of the deck at the alignment is " << feet(lowest) << " to " << feet(highest) << " above it.";
   }
   WBFL::System::Logger::Info(os.str());
}
