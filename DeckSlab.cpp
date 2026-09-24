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
#include "DeckSlab.h"
#include "Geometry.h"
#include "BeamSpacing.h"

#include <map>
#include <numeric>

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/element.h>
#include <ifcgeom/kernels/opencascade/opencascade_kernel.h>

#include <IFace\Project.h>
#include <IFace/Alignment.h>
#include <EAF/EAFProgress.h>
#include <EAF/AutoProgress.h>


#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file


// returns the id of the IfcSlab... assumes there is only one
int get_slab_id(ifcopenshell::file& file)
{
   auto slabs = file.instances_by_type<IfcSchema::IfcSlab>();
   auto it = std::find_if(slabs.begin(), slabs.end(), [](const auto& slab) {return slab.PredefinedType() == IfcSchema::IfcSlabTypeEnum::IfcSlabType_FLOOR; });
   if (it == slabs.end())
      return 0; // no slabs

   auto slab = *it;
   return slab.id();
}

using EdgeProfile = std::vector<std::pair<double, double>>; // (station, offset) along a deck edge, sorted by station

// offset of an edge at station, interpolated between points and held constant beyond the ends
double get_edge_offset(const EdgeProfile& edge, double station)
{
   if (station <= edge.front().first)
      return edge.front().second;
   if (edge.back().first <= station)
      return edge.back().second;

   auto next = std::lower_bound(edge.begin(), edge.end(), station, [](const auto& point, double s) {return point.first < s; });
   auto prev = std::prev(next);
   double ds = next->first - prev->first;
   return ds <= 0 ? prev->second : prev->second + (station - prev->first) * (next->second - prev->second) / ds;
}

// Douglas-Peucker simplification of an edge, measuring the deviation in offset.
// Keeps the fewest points that reproduce the edge within tolerance.
EdgeProfile simplify_edge(const EdgeProfile& edge, double tolerance)
{
   if (edge.size() < 3)
      return edge;

   std::vector<bool> keep(edge.size(), false);
   keep.front() = keep.back() = true;
   std::vector<std::pair<size_t, size_t>> ranges{ {0, edge.size() - 1} };
   while (!ranges.empty())
   {
      auto [first, last] = ranges.back();
      ranges.pop_back();

      const auto& [s1, o1] = edge[first];
      const auto& [s2, o2] = edge[last];
      double max_deviation = 0;
      size_t max_idx = first;
      for (size_t i = first + 1; i < last; i++)
      {
         double chord_offset = (s2 == s1) ? o1 : o1 + (edge[i].first - s1) * (o2 - o1) / (s2 - s1);
         double deviation = fabs(edge[i].second - chord_offset);
         if (max_deviation < deviation)
         {
            max_deviation = deviation;
            max_idx = i;
         }
      }

      if (tolerance < max_deviation)
      {
         keep[max_idx] = true;
         ranges.emplace_back(first, max_idx);
         ranges.emplace_back(max_idx, last);
      }
   }

   EdgeProfile result;
   for (size_t i = 0; i < edge.size(); i++)
   {
      if (keep[i])
         result.push_back(edge[i]);
   }
   return result;
}

// Station and offset of the vertices of a wire, sorted by station
EdgeProfile get_edge_profile(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Wire& wire)
{
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> point;
   point.CoCreateInstance(CLSID_Point2d);

   std::set<int> vertices;
   for (const auto& e : wire.edges)
   {
      vertices.insert(e[0]);
      vertices.insert(e[1]);
   }

   EdgeProfile edge;
   for (auto v : vertices)
   {
      point->Move(wire.verts[v].x(), wire.verts[v].y());
      Float64 station, offset;
      pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &station, &offset);
      edge.emplace_back(station, offset);
   }
   std::sort(edge.begin(), edge.end());
   return edge;
}

// Deck edge points that reproduce the left and right deck edges within tolerance.
// Where an edge doesn't change offset between points, it parallels the alignment.
std::vector<CDeckPoint> create_deck_points(const EdgeProfile& left, const EdgeProfile& right, double tolerance)
{
   // stations where either edge changes direction
   std::vector<double> stations;
   for (const auto* edge : { &left, &right })
   {
      for (const auto& [station, offset] : simplify_edge(*edge, tolerance))
         stations.push_back(station);
   }
   std::sort(stations.begin(), stations.end());
   stations.erase(std::unique(stations.begin(), stations.end(), [tolerance](double a, double b) {return fabs(a - b) <= tolerance; }), stations.end());

   // round to a micrometer to remove floating point noise from the geometry (e.g. 4.572000000003 m)
   const double accuracy = 1.0e-6;

   std::vector<CDeckPoint> deck_points;
   for (auto station : stations)
   {
      CDeckPoint deck_point;
      deck_point.Station = RoundOff(station, accuracy);
      deck_point.MeasurementType = pgsTypes::OffsetMeasurementType::omtAlignment; // offsets are computed relative to the alignment
      deck_point.LeftEdge = RoundOff(-get_edge_offset(left, station), accuracy); // offsets are positive to the right, LeftEdge is positive to the left
      deck_point.RightEdge = RoundOff(get_edge_offset(right, station), accuracy);
      deck_points.push_back(deck_point);
   }

   // an edge that doesn't change offset over a transition parallels the alignment
   for (size_t i = 0; i < deck_points.size(); i++)
   {
      bool bLast = (i == deck_points.size() - 1);
      deck_points[i].LeftTransitionType = (bLast || fabs(deck_points[i].LeftEdge - deck_points[i + 1].LeftEdge) <= tolerance) ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;
      deck_points[i].RightTransitionType = (bLast || fabs(deck_points[i].RightEdge - deck_points[i + 1].RightEdge) <= tolerance) ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;
   }

   // remove points in the middle of runs that parallel the alignment on both sides - they don't change the edges
   auto is_parallel = [](const CDeckPoint& p) {return p.LeftTransitionType == pgsTypes::DeckPointTransitionType::dptParallel && p.RightTransitionType == pgsTypes::DeckPointTransitionType::dptParallel; };
   for (size_t i = 1; i < deck_points.size(); )
   {
      bool bLast = (i == deck_points.size() - 1);
      if (is_parallel(deck_points[i - 1]) && (bLast || is_parallel(deck_points[i])))
         deck_points.erase(deck_points.begin() + i); // same edges as the previous point
      else
         i++;
   }

   return deck_points;
}

std::optional<Mesh> get_deck_mesh(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file)
{
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Processing deck geometry"));

   auto slab_id = get_slab_id(file);
   ASSERT(slab_id != 0); // should not be calling into this function if there isn't a deck slab

   ifcopenshell::geom::settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   // set up filter for the geometry iterator
   ifcopenshell::geom::instance_id_filter filter(true, false, { slab_id });
   std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

   std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
   ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters,1);
   if (!iterator.initialize())
   {
      WBFL::System::Logger::Info(_T("Unable to process deck slab geometry. Deck geometry will not be imported."));
      return std::nullopt;
   }

   auto element = iterator.get();
   auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
   auto geometry = triangulation->geometry_pointer();
   const auto& verts = geometry->verts();
   const auto& faces = geometry->faces();

   std::vector<Eigen::Vector3d> v;
   v.reserve(verts.size() / 3);
   for (auto i = 0; i < verts.size(); i += 3)
   {
      v.emplace_back(verts[i], verts[i + 1], verts[i + 2]);
   }

   std::vector<std::array<int, 3>> f;
   f.reserve(faces.size() / 3);
   for (auto i = 0; i < faces.size(); i += 3)
   {
      f.emplace_back(std::array<int, 3>({ faces[i],faces[i + 1],faces[i + 2] }));
   }

   return Mesh(v, f);
}

namespace
{
   bool contains_haunch(const std::string& text)
   {
      std::string lower(text);
      std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) {return (char)std::tolower(c); });
      return lower.find("haunch") != std::string::npos;
   }

   // Haunches that are modeled as elements separate from the deck slab: parts aggregated under the deck
   // and slabs or building element parts identified as haunches by their ObjectType or Name
   std::set<int> get_haunch_ids(ifcopenshell::file& file, int slab_id)
   {
      std::set<int> ids;

      auto deck = file.instance_by_id(slab_id).as<IfcSchema::IfcSlab>();
      for (auto& rel : deck.IsDecomposedBy())
      {
         for (auto& part : rel.RelatedObjects())
         {
            if (part.as<IfcSchema::IfcSlab>() || part.as<IfcSchema::IfcBuildingElementPart>())
               ids.insert(part.id());
         }
      }

      auto is_haunch = [](const auto& element) {return (element.ObjectType() && contains_haunch(*element.ObjectType())) || (element.Name() && contains_haunch(*element.Name())); };
      for (auto& slab : file.instances_by_type<IfcSchema::IfcSlab>())
      {
         if (slab.id() != slab_id && is_haunch(slab))
            ids.insert(slab.id());
      }
      for (auto& part : file.instances_by_type<IfcSchema::IfcBuildingElementPart>())
      {
         if (is_haunch(part))
            ids.insert(part.id());
      }

      return ids;
   }
}

std::optional<Mesh> get_deck_concrete_mesh(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file)
{
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Processing deck and haunch geometry"));

   auto slab_id = get_slab_id(file);
   ASSERT(slab_id != 0); // should not be calling into this function if there isn't a deck slab

   auto haunch_ids = get_haunch_ids(file, slab_id);
   if (haunch_ids.empty())
      return std::nullopt; // the deck slab is all of the deck concrete

   std::_tostringstream os;
   os << haunch_ids.size() << _T(" haunch elements found. They are measured with the deck slab.");
   WBFL::System::Logger::Info(os.str().c_str());

   std::set<int> ids(haunch_ids);
   ids.insert(slab_id);

   ifcopenshell::geom::settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   ifcopenshell::geom::instance_id_filter filter(true, false, ids);
   std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

   std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
   ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, 1);
   if (!iterator.initialize())
   {
      WBFL::System::Logger::Info(_T("Unable to process deck slab geometry. Deck geometry will not be imported."));
      return std::nullopt;
   }

   // all of the deck concrete in one mesh. The haunches are cast with the deck
   std::vector<Eigen::Vector3d> v;
   std::vector<std::array<int, 3>> f;
   do
   {
      auto element = iterator.get();
      auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
      auto geometry = triangulation->geometry_pointer();
      const auto& verts = geometry->verts();
      const auto& faces = geometry->faces();

      int base = (int)v.size();
      for (size_t i = 0; i < verts.size(); i += 3)
         v.emplace_back(verts[i], verts[i + 1], verts[i + 2]);
      for (size_t i = 0; i < faces.size(); i += 3)
         f.emplace_back(std::array<int, 3>({ base + faces[i], base + faces[i + 1], base + faces[i + 2] }));
   } while (iterator.next());

   return Mesh(v, f);
}

std::vector<CDeckPoint> get_deck_edge_points(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& m)
{

   // 1. Decompose mesh into smooth components
   std::vector<Mesh> decomposed = m.decompose();

   // 2. Select top component (max Z of vertex + avg normal)
   Mesh top_component = get_top_mesh(decomposed);

   // 3. Extract boundary wire
   Wire boundary = top_component.boundary();

   // 4. Decompose boundary into smooth segments
   std::vector<Wire> segments = boundary.decompose();
   if (segments.size() < 2)
   {
      WBFL::System::Logger::Info(_T("Unable to identify deck slab edges from the slab geometry."));
      return {};
   }

   // 5. Assume the two longest segments are the slab edges
   std::sort(segments.begin(), segments.end(), [](const Wire& a, const Wire& b) {return a.plan_length() < b.plan_length(); });
   EdgeProfile edge1 = get_edge_profile(pBroker, segments[segments.size() - 1]);
   EdgeProfile edge2 = get_edge_profile(pBroker, segments[segments.size() - 2]);

   // 6. The left edge has the smaller offsets (offsets are positive to the right)
   auto average_offset = [](const EdgeProfile& edge) {return std::accumulate(edge.begin(), edge.end(), 0.0, [](double sum, const auto& point) {return sum + point.second; }) / edge.size(); };
   bool bEdge1IsLeft = average_offset(edge1) < average_offset(edge2);
   const EdgeProfile& left = bEdge1IsLeft ? edge1 : edge2;
   const EdgeProfile& right = bEdge1IsLeft ? edge2 : edge1;

   // 7. Fewest deck points that reproduce the edges
   double tolerance = WBFL::Units::ConvertToSysUnits(0.125, WBFL::Units::Measure::Inch);
   return create_deck_points(left, right, tolerance);
}

namespace
{
   const Float64 inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);

   // A line across the deck, normal to the alignment at a station.
   // Points are located by offset from the alignment (positive to the right)
   class TransverseLine
   {
   public:
      TransverseLine(std::shared_ptr<WBFL::EAF::Broker> pBroker, Float64 station)
      {
         GET_IFACE2(pBroker, IRoadway, pAlignment);
         CComPtr<IPoint2d> point;
         pAlignment->GetPoint(station, 0.0, nullptr, pgsTypes::PlanCoordinateType::pcGlobal, &point);
         Float64 x, y;
         point->Location(&x, &y);
         m_Origin = Eigen::Vector2d(x, y);

         CComPtr<IDirection> normal;
         pAlignment->GetBearingNormal(station, &normal);
         Float64 angle;
         normal->get_Value(&angle);
         m_Direction = Eigen::Vector2d(cos(angle), sin(angle));

         // make the direction point to the right, the direction of positive offsets
         point->Move(m_Origin.x() + m_Direction.x(), m_Origin.y() + m_Direction.y());
         Float64 s, offset;
         pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &s, &offset);
         if (offset < 0)
            m_Direction = -m_Direction;
      }

      Eigen::Vector2d At(Float64 offset) const { return m_Origin + offset * m_Direction; }

   private:
      Eigen::Vector2d m_Origin;
      Eigen::Vector2d m_Direction;
   };

   // Deck top and bottom along a transverse line
   struct DeckProfile
   {
      Float64 left_edge, right_edge; // offsets of the deck edges
      std::function<std::optional<std::pair<double, double>>(Float64)> z_range; // bottom and top of deck at an offset
      Float64 thickness(Float64 offset) const { auto z = z_range(offset); return z ? z->second - z->first : 0.0; }
   };

   // Offsets of the left and right edges of a mesh along a transverse line, searching outward from an offset
   // that is on the mesh. Steps then bisects to the nearest 1/64 inch.
   std::optional<std::pair<Float64, Float64>> find_extents(const Mesh& mesh, const TransverseLine& line, Float64 start_offset, Float64 max_distance)
   {
      auto on_mesh = [&](Float64 offset) { auto p = line.At(offset); return mesh.z_range_at(p.x(), p.y()).has_value(); };
      if (!on_mesh(start_offset))
         return std::nullopt;

      auto search = [&](Float64 direction)
         {
            Float64 inside = start_offset;
            Float64 step = inch;
            Float64 outside = inside + direction * step;
            while (on_mesh(outside) && fabs(outside - start_offset) < max_distance)
            {
               inside = outside;
               outside += direction * step;
            }
            while (inch / 64 < fabs(outside - inside))
            {
               Float64 mid = (inside + outside) / 2;
               (on_mesh(mid) ? inside : outside) = mid;
            }
            return inside;
         };

      return std::make_pair(search(-1.0), search(1.0));
   }

   Float64 median(std::vector<Float64> values)
   {
      std::sort(values.begin(), values.end());
      size_t n = values.size();
      return n % 2 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2;
   }

   // soffit depth below the top of the deck at offset, from a least squares line through the soffit depths at the offsets.
   // if max_residual is given, it is the largest distance from a soffit depth to the line
   Float64 extrapolate_soffit(const DeckProfile& profile, const std::vector<Float64>& offsets, Float64 offset, Float64* max_residual = nullptr)
   {
      std::vector<Float64> depths;
      for (auto o : offsets)
         depths.push_back(profile.thickness(o));

      Float64 n = (Float64)offsets.size();
      Float64 sx = 0, sy = 0, sxx = 0, sxy = 0;
      for (size_t i = 0; i < offsets.size(); i++)
      {
         sx += offsets[i]; sy += depths[i]; sxx += offsets[i] * offsets[i]; sxy += offsets[i] * depths[i];
      }
      Float64 det = n * sxx - sx * sx;
      Float64 slope = fabs(det) < 1.0e-12 ? 0.0 : (n * sxy - sx * sy) / det;
      Float64 intercept = (sy - slope * sx) / n;

      if (max_residual)
      {
         *max_residual = 0;
         for (size_t i = 0; i < offsets.size(); i++)
            *max_residual = std::max(*max_residual, fabs(depths[i] - (intercept + slope * offsets[i])));
      }

      return intercept + slope * offset;
   }
}

DeckSectionData analyze_deck_section(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const GirderLayout& layout, const std::vector<Float64>& pier_stations, const CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Measuring deck cross sections"));

   DeckSectionData data;

   std::vector<Float64> gross_depths;
   std::array<std::vector<Float64>, 2> edge_depths;
   std::array<std::vector<pgsTypes::DeckOverhangTaper>, 2> tapers;
   std::vector<Float64> fillets; // width of the transition from the slab soffit to the bottom of the haunch
   bool bDeckHasHaunch = false;

   // Measure the deck at the quarter points of each span
   for (SpanIndexType spanIdx = 0; spanIdx < layout.size(); spanIdx++)
   {
      const auto& girders = layout[spanIdx];
      if (girders.size() < 2)
         continue;

      for (auto fraction : { 0.25, 0.50, 0.75 })
      {
         Float64 station = pier_stations[spanIdx] + fraction * (pier_stations[spanIdx + 1] - pier_stations[spanIdx]);
         TransverseLine line(pBroker, station);

         DeckProfile profile;
         profile.z_range = [&](Float64 offset) { auto p = line.At(offset); return deck.z_range_at(p.x(), p.y()); };

         // girder top flange extents along the line (offsets), left to right
         std::vector<std::pair<Float64, Float64>> flanges;
         for (const auto& girder : girders)
         {
            Float64 cl_offset = girder.start_offset + (station - girder.start_station) * (girder.end_offset - girder.start_offset) / (girder.end_station - girder.start_station);
            auto extents = find_extents(*girder.mesh, line, cl_offset, WBFL::Units::ConvertToSysUnits(10.0, WBFL::Units::Measure::Feet));
            if (extents)
               flanges.push_back(*extents);
         }
         if (flanges.size() != girders.size())
            continue; // couldn't find all the girders on this line

         // deck edges, searching outward from the first girder
         auto deck_extents = find_extents(deck, line, (flanges.front().first + flanges.front().second) / 2, WBFL::Units::ConvertToSysUnits(200.0, WBFL::Units::Measure::Feet));
         if (!deck_extents)
            continue;
         profile.left_edge = deck_extents->first;
         profile.right_edge = deck_extents->second;

         // Gross depth: thickness in the bays, 2 inches outside of the flanges (clear of fillets and of crown points in the middle of bays)
         Float64 gross_depth_here = 0;
         {
            std::vector<Float64> depths;
            for (size_t i = 0; i + 1 < flanges.size(); i++)
            {
               depths.push_back(profile.thickness(flanges[i].second + 2 * inch));
               depths.push_back(profile.thickness(flanges[i + 1].first - 2 * inch));
            }
            gross_depth_here = median(depths);
            gross_depths.push_back(gross_depth_here);
         }

         // Haunch in the deck solid? Thickness over the middle of the flanges exceeds the gross depth
         std::vector<Float64> over_flange;
         for (const auto& [left, right] : flanges)
            over_flange.push_back(profile.thickness((left + right) / 2));
         bool bHaunchHere = gross_depth_here + inch / 8 < median(over_flange);
         bDeckHasHaunch |= bHaunchHere;

         // Fillet: width of the transition at the flange edges of interior bays
         if (bHaunchHere)
         {
            for (size_t i = 0; i + 1 < flanges.size(); i++)
            {
               for (auto [tip, direction] : { std::make_pair(flanges[i].second, 1.0), std::make_pair(flanges[i + 1].first, -1.0) })
               {
                  // walk from 2 inches in the bay toward and over the flange in 1/32 inch steps
                  Float64 haunch_depth = profile.thickness(tip - direction * inch); // 1 inch over the flange
                  Float64 last_gross = tip + direction * 2 * inch;
                  Float64 first_haunch = tip - direction * inch;
                  for (Float64 d = 2 * inch; -inch <= d; d -= inch / 32)
                  {
                     Float64 t = profile.thickness(tip + direction * d);
                     if (t <= gross_depth_here + inch / 32)
                        last_gross = tip + direction * d;
                     if (haunch_depth - inch / 32 <= t)
                     {
                        first_haunch = tip + direction * d;
                        break;
                     }
                  }
                  fillets.push_back(fabs(first_haunch - last_gross));
               }
            }
         }

         // Edge depths and overhang tapers
         for (auto side : { pgsTypes::stLeft, pgsTypes::stRight })
         {
            Float64 deck_edge = (side == pgsTypes::stLeft ? profile.left_edge : profile.right_edge);
            Float64 tip = (side == pgsTypes::stLeft ? flanges.front().first : flanges.back().second);
            Float64 inward = (side == pgsTypes::stLeft ? 1.0 : -1.0);
            Float64 edge = deck_edge + inward * inch / 4; // clear of the edge of the mesh

            // soffit of the overhang, from 1 inch inside the edge to 3 inches outside the flange tip
            std::vector<Float64> soffit_offsets;
            Float64 overhang = fabs(tip - edge);
            if (overhang < 6 * inch)
            {
               edge_depths[side].push_back(profile.thickness(edge));
               continue; // too short to see a taper
            }

            for (Float64 d = inch; d <= overhang - 3 * inch; d += inch)
               soffit_offsets.push_back(edge + inward * d);

            // the soffit of a tapered overhang isn't parallel to the deck top, so where the soffit near the edge
            // is straight, the edge depth is the soffit extrapolated to the edge. otherwise (e.g. a chamfer or drip)
            // it's the depth just inside the edge
            std::vector<Float64> edge_offsets;
            for (Float64 d = 0; d <= 3 * inch; d += inch / 4)
               edge_offsets.push_back(edge + inward * d);
            Float64 max_residual;
            Float64 edge_depth = extrapolate_soffit(profile, edge_offsets, deck_edge, &max_residual);
            if (inch / 32 < max_residual)
               edge_depth = profile.thickness(edge);
            edge_depths[side].push_back(edge_depth);

            Float64 depth_at_tip = extrapolate_soffit(profile, soffit_offsets, tip);

            // depth to the bottom of the haunch (the top of the girder in the deck model), or to the top of the girder
            // if the haunch isn't part of the deck solid
            Float64 haunch_depth;
            if (bHaunchHere)
            {
               haunch_depth = profile.thickness(tip + inward * inch);
            }
            else
            {
               const auto& girder = (side == pgsTypes::stLeft ? girders.front() : girders.back());
               auto p = line.At(tip + inward * inch);
               auto deck_z = deck.z_range_at(p.x(), p.y());
               auto girder_z = girder.mesh->z_range_at(p.x(), p.y());
               if (!deck_z || !girder_z)
                  continue;
               haunch_depth = deck_z->second - girder_z->second;
            }

            const Float64 tolerance = inch / 4;
            pgsTypes::DeckOverhangTaper taper;
            if (fabs(depth_at_tip - edge_depth) <= tolerance)
               taper = pgsTypes::dotNone; // constant thickness overhang
            else if (fabs(depth_at_tip - haunch_depth) <= tolerance)
               taper = pgsTypes::dotTopTopFlange; // soffit runs to the top of the top flange
            else if (haunch_depth < depth_at_tip)
               taper = pgsTypes::dotBottomTopFlange; // soffit runs below the top of the top flange
            else
               taper = (fabs(depth_at_tip - edge_depth) < fabs(depth_at_tip - haunch_depth)) ? pgsTypes::dotNone : pgsTypes::dotTopTopFlange; // closest

            tapers[side].push_back(taper);
         }
      }
   }

   if (!gross_depths.empty())
      data.gross_depth = median(gross_depths);

   for (auto side : { pgsTypes::stLeft, pgsTypes::stRight })
   {
      if (!edge_depths[side].empty())
         data.edge_depth[side] = median(edge_depths[side]);

      if (!tapers[side].empty())
      {
         // most common taper
         std::map<pgsTypes::DeckOverhangTaper, int> count;
         for (auto t : tapers[side]) count[t]++;
         data.overhang_taper[side] = std::max_element(count.begin(), count.end(), [](const auto& a, const auto& b) {return a.second < b.second; })->first;
      }
   }

   if (bDeckHasHaunch && !fillets.empty())
   {
      Float64 fillet = median(fillets);
      if (fillet < inch / 8)
      {
         data.haunch_shape = pgsTypes::hsSquare;
      }
      else
      {
         data.haunch_shape = pgsTypes::hsFilleted;
         data.fillet = fillet;
      }
   }

   // Slab offset ("A" dimension): top of deck to top of girder at CL girder, CL bearing
   data.slab_offset.resize(layout.size());
   for (SpanIndexType spanIdx = 0; spanIdx < layout.size(); spanIdx++)
   {
      const auto* pSpan = bridge_desc.GetSpan(spanIdx);
      auto start_end_distance = pSpan->GetPrevPier()->GetGirderEndDistance(pgsTypes::Ahead).first;
      auto end_end_distance = pSpan->GetNextPier()->GetGirderEndDistance(pgsTypes::Back).first;

      for (const auto& girder : layout[spanIdx])
      {
         Eigen::Vector2d direction = (girder.end - girder.start).normalized();
         std::array<std::optional<Float64>, 2> slab_offset;
         Float64 girder_length = (girder.end - girder.start).norm();
         // CL bearing from the bearing geometry, otherwise located with the end distance
         Eigen::Vector2d start_bearing = girder.start_bearing.value_or(Eigen::Vector2d(girder.start + start_end_distance * direction));
         Eigen::Vector2d end_bearing = girder.end_bearing.value_or(Eigen::Vector2d(girder.end - end_end_distance * direction));
         for (auto [end, point, inward] : { std::make_tuple(pgsTypes::metStart, start_bearing, Eigen::Vector2d(direction)),
                                            std::make_tuple(pgsTypes::metEnd, end_bearing, Eigen::Vector2d(-direction)) })
         {
            // top of deck to top of girder at a distance from the CL bearing, toward the middle of the girder
            // The deck concrete must include the slab, not just a haunch that runs beyond the end of the slab
            auto slab_offset_at = [&](Float64 distance) -> std::optional<Float64>
               {
                  Eigen::Vector2d p = point + distance * inward;
                  auto deck_z = deck.z_range_at(p.x(), p.y());
                  auto girder_z = girder.mesh->z_range_at(p.x(), p.y());
                  if (!deck_z || !girder_z)
                     return std::nullopt;
                  if (data.gross_depth && deck_z->second - deck_z->first < *data.gross_depth - inch / 4)
                     return std::nullopt; // no slab here
                  return deck_z->second - girder_z->second;
               };

            slab_offset[end] = slab_offset_at(0.0);
            if (!slab_offset[end])
            {
               // The deck doesn't extend over the bearing (e.g. an integral abutment where the deck
               // ends at the abutment). Extrapolate from where the deck is over the girder.
               const Float64 step = 6 * inch;
               Float64 first = step;
               while (first < girder_length / 4 && !slab_offset_at(first))
                  first += step;

               std::vector<std::pair<Float64, Float64>> samples; // (distance, slab offset)
               for (int i = 0; i < 4; i++)
               {
                  Float64 d = first + i * WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Feet);
                  if (auto a = slab_offset_at(d))
                     samples.emplace_back(d, *a);
               }

               if (2 <= samples.size())
               {
                  Float64 n = (Float64)samples.size(), sx = 0, sy = 0, sxx = 0, sxy = 0;
                  for (auto [d, a] : samples) { sx += d; sy += a; sxx += d * d; sxy += d * a; }
                  Float64 slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
                  slab_offset[end] = (sy - slope * sx) / n; // at distance = 0
                  data.slab_offset_extrapolated = true;
               }
            }
         }
         data.slab_offset[spanIdx].push_back(slab_offset);
      }
   }

   return data;
}

bool create_alignment_from_deck(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file)
{
   // The bridge has a deck, but the model does not have an alignment.
   // We can't compute station and offset of the deck edge points without an alignment.
   // Assume the alignment is straight. Create the alignment using the center point of the deck ends to establish the start point and direction.
   // This is a best effort approach to get an alignment that is usable for stationing the deck geometry. 
   // It will not be correct, but it will at least allow the user to see the deck geometry in the correct location and station the deck geometry.

   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Progressing deck geometry to establish an alignment"));

   auto slab_id = get_slab_id(file);
   ASSERT(slab_id != 0); // should not be calling into this function if there isn't a deck slab

   ifcopenshell::geom::settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   // set up filter for the geometry iterator
   ifcopenshell::geom::instance_id_filter filter(true, false, { slab_id });
   std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

   std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
   ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, 1);
   if (!iterator.initialize())
   {
      WBFL::System::Logger::Info(_T("Unable to process deck slab geometry."));
      return false;
   }

   do
   {
      auto element = iterator.get();

      auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
      auto geometry = triangulation->geometry_pointer();
      const auto& verts = geometry->verts();
      const auto& faces = geometry->faces();

      std::vector<Eigen::Vector3d> v;
      v.reserve(verts.size() / 3);
      for (auto i = 0; i < verts.size(); i += 3)
      {
         v.emplace_back(verts[i], verts[i + 1], verts[i + 2]);
      }

      std::vector<std::array<int, 3>> f;
      f.reserve(faces.size() / 3);
      for (auto i = 0; i < faces.size(); i += 3)
      {
         // This line of code should be correct, however....
         f.emplace_back(std::array<int, 3>({ faces[i],faces[i + 1],faces[i + 2] }));

         // There appears to be a bug in the IfcOpenShell processing of IfcSectionedSolidHorizontal
         // The vertex indices of the triangular faces are in the wrong direction resulting in an
         // inward surface normal. Working around this problem by manually reversing the order
         //f.emplace_back(std::array<int, 3>({ faces[i+2],faces[i + 1],faces[i] }));
      }

      // 1. Build mesh
      Mesh m(v, f);

      // 2. Decompose mesh into smooth components
      std::vector<Mesh> decomposed = m.decompose();

      // 3. Select top component (max Z of vertex + avg normal)
      Mesh top_component = get_top_mesh(decomposed);

      // 4. Extract boundary wire
      Wire boundary = top_component.boundary();

      // 5. Decompose boundary into smooth segments
      std::vector<Wire> segments = boundary.decompose();

      // 6. Sort segments by length
      std::sort(segments.begin(), segments.end(),
         [](const Wire& a, const Wire& b)
         {
            return a.plan_length() < b.plan_length();
         });

      if (segments.size() < 2)
      {
         WBFL::System::Logger::Info(_T("Unable to identify deck slab ends from the slab geometry."));
         continue;
      }

      // 7. Assume two shortest segments are the ends of the slab
      Wire start_seg = segments[0];
      Wire end_seg = segments[1];

      // 8. Get center point of edge wire bounding-box
      Eigen::Vector3d c1 = start_seg.bbox_center();
      Eigen::Vector3d c2 = end_seg.bbox_center();

      // 9. Create alignment between those two points
      GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
      AlignmentData2 alignmentData = pRoadwayData->GetAlignmentData2();
      alignmentData.xRefPoint = c1.x();
      alignmentData.yRefPoint = c1.y();
      alignmentData.Direction = atan2(c2.y() - c1.y(), c2.x() - c1.x());
      pRoadwayData->SetAlignmentData2(alignmentData);

      return true;

   } while (iterator.next());

   return false;
}