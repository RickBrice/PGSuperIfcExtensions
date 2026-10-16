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
#include "GirderSection.h"
#include "Geometry.h"

#include <psgLib/GirderLibraryEntry.h>
#include <IFace/BeamFactory.h>
#include <WBFLGenericBridge.h>

#include <algorithm>
#include <functional>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <chrono>
#include <iomanip>
#include <map>
#include <numeric>
#include <sstream>

namespace
{
   using Loop = SectionOutline::Loop;

   double signed_area(const Loop& loop)
   {
      double a = 0;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& p = loop[i];
         const auto& q = loop[(i + 1) % loop.size()];
         a += p.x() * q.y() - q.x() * p.y();
      }
      return a / 2;
   }

   // area, first moment, and second moment about the x axis of a loop, counterclockwise positive
   void loop_moments(const Loop& loop, double& area, double& qx, double& ixx)
   {
      area = qx = ixx = 0;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& p = loop[i];
         const auto& q = loop[(i + 1) % loop.size()];
         double cross = p.x() * q.y() - q.x() * p.y();
         area += cross / 2;
         qx += cross * (p.y() + q.y()) / 6;
         ixx += cross * (p.y() * p.y() + p.y() * q.y() + q.y() * q.y()) / 12;
      }
   }

   // the length of a horizontal line at y inside a loop
   double loop_width(const Loop& loop, double y)
   {
      std::vector<double> xs;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& p = loop[i];
         const auto& q = loop[(i + 1) % loop.size()];
         if ((p.y() <= y) != (q.y() <= y))
            xs.push_back(p.x() + (y - p.y()) / (q.y() - p.y()) * (q.x() - p.x()));
      }
      std::sort(xs.begin(), xs.end());
      double w = 0;
      for (size_t i = 0; i + 1 < xs.size(); i += 2)
         w += xs[i + 1] - xs[i];
      return w;
   }

   // the area of a loop above a horizontal line (the loop clipped to y >= y0), counterclockwise positive
   double loop_area_above(const Loop& loop, double y0)
   {
      Loop clipped;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& p = loop[i];
         const auto& q = loop[(i + 1) % loop.size()];
         bool p_in = y0 <= p.y(), q_in = y0 <= q.y();
         if (p_in)
            clipped.push_back(p);
         if (p_in != q_in)
            clipped.push_back(p + (y0 - p.y()) / (q.y() - p.y()) * (q - p));
      }
      return clipped.size() < 3 ? 0.0 : signed_area(clipped);
   }

   // orients the loops (solids counterclockwise, voids clockwise) and moves them so the origin is at the top center
   SectionOutline normalize(SectionOutline outline)
   {
      for (auto& loop : outline.solids)
         if (signed_area(loop) < 0) std::reverse(loop.begin(), loop.end());
      for (auto& loop : outline.voids)
         if (0 < signed_area(loop)) std::reverse(loop.begin(), loop.end());

      double xmin = DBL_MAX, xmax = -DBL_MAX, ymax = -DBL_MAX;
      for (const auto& loop : outline.solids)
      {
         for (const auto& p : loop)
         {
            xmin = std::min(xmin, p.x());
            xmax = std::max(xmax, p.x());
            ymax = std::max(ymax, p.y());
         }
      }
      Eigen::Vector2d origin((xmin + xmax) / 2, ymax);
      for (auto* loops : { &outline.solids, &outline.voids })
         for (auto& loop : *loops)
            for (auto& p : loop)
               p -= origin;
      return outline;
   }

   double inches(double length) { return WBFL::Units::ConvertFromSysUnits(length, WBFL::Units::Measure::Inch); }
}

std::vector<SectionOutline::Loop> CutMesh(const Mesh& mesh, const Eigen::Vector2d& point2d, Eigen::Vector2d axis)
{
   if (axis.norm() <= 0)
      return {};
   axis.normalize();

   Eigen::Vector3d point(point2d.x(), point2d.y(), 0);
   Eigen::Vector3d normal(axis.x(), axis.y(), 0);
   Eigen::Vector3d left(-axis.y(), axis.x(), 0);

   // signed distance of each vertex from the plane. A vertex on the plane counts as in front of it
   std::vector<double> d(mesh.verts.size());
   for (size_t i = 0; i < mesh.verts.size(); i++)
      d[i] = (mesh.verts[i] - point).dot(normal);

   // Each face the plane cuts gives a segment between two cut edges. A cut edge is a node, identified by its vertices,
   // so the segments of neighboring faces join (the mesh has welded vertices)
   using Node = std::pair<int, int>;
   std::map<Node, Eigen::Vector2d> node_points;
   std::map<Node, std::vector<Node>> links;
   auto cut = [&](int a, int b) -> std::optional<Node>
      {
         if ((0 <= d[a]) == (0 <= d[b]))
            return std::nullopt;
         Node node(std::min(a, b), std::max(a, b));
         if (node_points.find(node) == node_points.end())
         {
            double t = d[node.first] / (d[node.first] - d[node.second]);
            Eigen::Vector3d p = mesh.verts[node.first] + t * (mesh.verts[node.second] - mesh.verts[node.first]);
            node_points.emplace(node, Eigen::Vector2d((p - point).dot(left), p.z()));
         }
         return node;
      };

   for (const auto& face : mesh.faces)
   {
      std::vector<Node> nodes;
      for (int i = 0; i < 3; i++)
      {
         if (auto node = cut(face[i], face[(i + 1) % 3]))
            nodes.push_back(*node);
      }
      if (nodes.size() == 2)
      {
         links[nodes[0]].push_back(nodes[1]);
         links[nodes[1]].push_back(nodes[0]);
      }
   }

   // walk the closed loops
   std::vector<Loop> loops;
   std::set<Node> visited;
   for (const auto& [first, neighbors] : links)
   {
      if (visited.count(first))
         continue;

      Loop loop;
      Node previous = first, current = first;
      bool closed = false;
      while (true)
      {
         visited.insert(current);
         loop.push_back(node_points[current]);
         const auto& next_nodes = links[current];
         auto next = std::find_if(next_nodes.begin(), next_nodes.end(), [&](const Node& n) {return n != previous && !visited.count(n); });
         if (next == next_nodes.end())
         {
            closed = std::find(next_nodes.begin(), next_nodes.end(), first) != next_nodes.end() && 3 <= loop.size();
            break;
         }
         previous = current;
         current = *next;
      }
      if (closed)
         loops.push_back(std::move(loop));
   }

   std::sort(loops.begin(), loops.end(), [](const Loop& a, const Loop& b) {return fabs(signed_area(a)) > fabs(signed_area(b)); });
   return loops;
}

std::optional<SectionOutline> CutGirderSection(const Mesh& mesh, const Eigen::Vector2d& start, const Eigen::Vector2d& end, double fraction)
{
   auto loops = CutMesh(mesh, start + fraction * (end - start), end - start);
   if (loops.empty())
      return std::nullopt;

   // the largest loop is the outside of the section, loops inside it are voids

   // A girder can be rotated about its axis (e.g. to the cross slope of the deck). The section is turned back to the
   // rotation, within 15 degrees, where it's the least high: girders have flat tops and bottoms, so that's upright, also
   // for sections that aren't symmetric. Less than a quarter degree is left alone (e.g. a slight crown on the top)
   {
      const auto& outer = loops.front();
      auto height = [&outer](double angle)
         {
            double s = sin(angle), c = cos(angle), ymin = DBL_MAX, ymax = -DBL_MAX;
            for (const auto& p : outer)
            {
               double y = s * p.x() + c * p.y();
               ymin = std::min(ymin, y);
               ymax = std::max(ymax, y);
            }
            return ymax - ymin;
         };
      const double step = M_PI / 3600; // 0.05 degree
      double best_angle = 0, best_height = height(0);
      for (double angle = -M_PI / 12; angle <= M_PI / 12; angle += step)
      {
         double h = height(angle);
         if (h < best_height)
         {
            best_height = h;
            best_angle = angle;
         }
      }
      if (M_PI / 720 < fabs(best_angle))
      {
         Eigen::Rotation2Dd rotation(best_angle);
         for (auto& loop : loops)
            for (auto& p : loop)
               p = rotation * p;
      }
   }
   SectionOutline outline;
   outline.solids.push_back(loops.front());
   for (size_t i = 1; i < loops.size(); i++)
   {
      // a loop is a void if its first point is inside the outside loop
      const auto& p = loops[i].front();
      if (0 < loop_width(loops.front(), p.y()))
         outline.voids.push_back(loops[i]);
   }

   return normalize(outline);
}

SectionOutline GetShapeOutline(IShape* pShape)
{
   auto points = [](IShape* shape)
      {
         Loop loop;
         CComPtr<IPoint2dCollection> poly_points;
         shape->get_PolyPoints(&poly_points);
         IndexType n = 0;
         poly_points->get_Count(&n);
         for (IndexType i = 0; i < n; i++)
         {
            CComPtr<IPoint2d> point;
            poly_points->get_Item(i, &point);
            Float64 x, y;
            point->Location(&x, &y);
            if (loop.empty() || (loop.back() - Eigen::Vector2d(x, y)).norm() > 1e-9)
               loop.emplace_back(x, y);
         }
         if (1 < loop.size() && (loop.back() - loop.front()).norm() < 1e-9)
            loop.pop_back(); // closed polygons repeat the first point
         return loop;
      };

   SectionOutline outline;
   CComQIPtr<ICompositeShape> composite(pShape);
   if (composite)
   {
      IndexType n = 0;
      composite->get_Count(&n);
      for (IndexType i = 0; i < n; i++)
      {
         CComPtr<ICompositeShapeItem> item;
         composite->get_Item(i, &item);
         CComPtr<IShape> shape;
         item->get_Shape(&shape);
         VARIANT_BOOL bVoid;
         item->get_Void(&bVoid);
         (bVoid == VARIANT_TRUE ? outline.voids : outline.solids).push_back(points(shape));
      }
   }
   else
   {
      outline.solids.push_back(points(pShape));
   }
   return normalize(outline);
}

std::optional<SectionOutline> GetLibraryEntryOutline(const GirderLibraryEntry* pEntry)
{
   try
   {
      auto factory = pEntry->GetBeamFactory();
      if (!factory)
         return std::nullopt;

      CComPtr<IGirderSection> section;
      factory->CreateGirderSection(nullptr, INVALID_ID, pEntry->GetDimensions(), -1, -1, &section);
      CComQIPtr<IShape> shape(section);
      if (!shape)
         return std::nullopt;

      return GetShapeOutline(shape);
   }
   catch (...)
   {
      return std::nullopt;
   }
}

SectionMeasures MeasureSection(const SectionOutline& outline)
{
   SectionMeasures m;
   m.outline = outline;
   double qx = 0, ixx = 0;
   double ymin = 0, xmin = DBL_MAX, xmax = -DBL_MAX;
   for (auto* loops : { &outline.solids, &outline.voids })
   {
      for (const auto& loop : *loops)
      {
         double a, q, i;
         loop_moments(loop, a, q, i);
         m.area += a;
         qx += q;
         ixx += i;
      }
   }
   for (const auto& loop : outline.solids)
   {
      for (const auto& p : loop)
      {
         ymin = std::min(ymin, p.y());
         xmin = std::min(xmin, p.x());
         xmax = std::max(xmax, p.x());
      }
   }

   m.height = -ymin;
   m.width = xmax - xmin;
   if (m.area <= 0 || !std::isfinite(m.area) || !std::isfinite(m.height) || !std::isfinite(m.width))
      return m;

   m.centroid_y = qx / m.area;
   m.ix = ixx - m.area * m.centroid_y * m.centroid_y;

   for (double y = -WIDTH_PROFILE_STEP / 2; -m.height < y; y -= WIDTH_PROFILE_STEP)
   {
      double w = 0;
      for (const auto& loop : outline.solids) w += loop_width(loop, y);
      for (const auto& loop : outline.voids) w -= loop_width(loop, y);
      m.width_profile.push_back(w);
   }

   for (double y = -AREA_PROFILE_STEP; -m.height < y; y -= AREA_PROFILE_STEP)
   {
      double a = 0;
      for (auto* loops : { &outline.solids, &outline.voids })
         for (const auto& loop : *loops)
            a += loop_area_above(loop, y);
      m.area_profile.push_back(a);
   }
   return m;
}

namespace
{
   struct Comparison
   {
      std::vector<std::string> differences; // not allowed
      std::vector<std::string> allowed;     // allowed, but worth reporting
      bool widths_differ = false;           // the last of the differences is the widths'
   };

   // Compares two sections. The widths may differ (by more than the length tolerance, at any depth within it) over a
   // total depth up to local_depth
   Comparison compare_sections(const SectionMeasures& model, const SectionMeasures& other, const SectionTolerances& tolerances, double local_depth)
   {
      Comparison result;
      auto length = [&](const char* name, double a, double b)
         {
            if (tolerances.length < fabs(a - b))
            {
               std::ostringstream os;
               os << std::fixed << std::setprecision(3) << name << " " << inches(a) << " in vs " << inches(b) << " in";
               result.differences.push_back(os.str());
            }
         };
      auto ratio = [&](const char* name, double a, double b, double unit, const char* unit_name)
         {
            if (b != 0 && tolerances.ratio < fabs(a - b) / fabs(b))
            {
               std::ostringstream os;
               os << std::fixed << std::setprecision(1) << name << " " << a / unit << " vs " << b / unit << " " << unit_name;
               result.differences.push_back(os.str());
            }
         };

      length("height", model.height, other.height);
      length("centroid (from the top)", -model.centroid_y, -other.centroid_y);
      ratio("area", model.area, other.area, WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch2), "in^2");
      ratio("I", model.ix, other.ix, WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch4), "in^4");

      // Width at each depth. A width may come from a depth within the tolerance, so a flange thickness that differs by
      // less than the tolerance isn't a difference
      int band = (int)ceil(tolerances.length / WIDTH_PROFILE_STEP - 1e-9);
      size_t n = std::min(model.width_profile.size(), other.width_profile.size());
      std::vector<size_t> mismatches;
      size_t largest = 0;
      for (size_t k = 0; k < n; k++)
      {
         bool match = false;
         for (int j = std::max(0, (int)k - band); !match && j <= std::min((int)n - 1, (int)k + band); j++)
            match = fabs(model.width_profile[k] - other.width_profile[j]) <= tolerances.length || fabs(model.width_profile[j] - other.width_profile[k]) <= tolerances.length;
         if (!match)
         {
            if (mismatches.empty() || fabs(model.width_profile[largest] - other.width_profile[largest]) < fabs(model.width_profile[k] - other.width_profile[k]))
               largest = k;
            mismatches.push_back(k);
         }
      }
      if (!mismatches.empty())
      {
         double depth = mismatches.size() * WIDTH_PROFILE_STEP;
         std::ostringstream os;
         os << std::fixed << std::setprecision(3) << "width differs over " << inches(depth) << " in of depth (the most "
            << inches((largest + 0.5) * WIDTH_PROFILE_STEP) << " in below the top: " << inches(model.width_profile[largest]) << " in vs "
            << inches(other.width_profile[largest]) << " in)";
         if (depth <= local_depth + 1e-9)
         {
            result.allowed.push_back(os.str());
         }
         else
         {
            result.differences.push_back(os.str());
            result.widths_differ = true;
         }
      }
      return result;
   }

   // As compare_sections, except that widths that differ over more than the local depth are allowed where the outlines
   // are within the length tolerance of each other everywhere (OutlineDistance). On a shallow taper (e.g. a flange's
   // underside) a small vertical offset is a large width difference: a model whose top and flange tips are a tenth of
   // an inch lower on one side (e.g. from a crown) differs in width over inches of depth
   Comparison compare_fitted_sections(const SectionMeasures& model, const SectionMeasures& other, const SectionTolerances& tolerances)
   {
      auto result = compare_sections(model, other, tolerances, tolerances.local_depth);
      if (result.widths_differ)
      {
         double distance = OutlineDistance(model.outline, other.outline);
         if (distance <= tolerances.length)
         {
            std::ostringstream os;
            os << result.differences.back() << std::fixed << std::setprecision(3) << ", but the outlines are within " << inches(distance) << " in of each other";
            result.differences.pop_back();
            result.allowed.push_back(os.str());
            result.widths_differ = false;
         }
      }
      return result;
   }
}

std::vector<std::string> ExactSectionComparison::Compare(const SectionMeasures& model, const SectionMeasures& other) const
{
   return compare_sections(model, other, m_Tolerances, 0.0).differences;
}

double ExactSectionComparison::Difference(const SectionMeasures& model, const SectionMeasures& other) const
{
   const auto& t = m_Tolerances;
   double d = fabs(model.height - other.height) / t.length;
   d += fabs(model.centroid_y - other.centroid_y) / t.length;
   if (other.area != 0) d += fabs(model.area - other.area) / fabs(other.area) / t.ratio;
   if (other.ix != 0) d += fabs(model.ix - other.ix) / fabs(other.ix) / t.ratio;

   // the area above each depth: a width difference of the tolerance over the whole height gives about 1
   size_t n = std::min(model.area_profile.size(), other.area_profile.size());
   double a = 0;
   for (size_t k = 0; k < n; k++)
      a += fabs(model.area_profile[k] - other.area_profile[k]) / ((k + 1) * AREA_PROFILE_STEP);
   if (0 < n)
      d += a / n / t.length;

   // the width at each depth, so a fit follows the outline (flange thicknesses, tapers, fillets) and doesn't trade the
   // width at one depth for the width at another: a width difference of the tolerance at every depth gives about 1
   size_t m = std::min(model.width_profile.size(), other.width_profile.size());
   double w = 0;
   for (size_t k = 0; k < m; k++)
      w += fabs(model.width_profile[k] - other.width_profile[k]);
   if (0 < m)
      d += w / m / t.length;
   return d;
}

std::vector<std::string> FittedSectionComparison::Compare(const SectionMeasures& model, const SectionMeasures& other) const
{
   return compare_fitted_sections(model, other, m_Tolerances).differences;
}

namespace
{
   namespace bg = boost::geometry;
   using BPoint = bg::model::d2::point_xy<double>;
   using BPolygon = bg::model::polygon<BPoint, false /*counterclockwise*/, false /*open*/>;
   using BMultiPolygon = bg::model::multi_polygon<BPolygon>;

   // the region of an outline: the union of its solids less its voids
   BMultiPolygon region(const SectionOutline& outline)
   {
      auto to_polygon = [](const SectionOutline::Loop& loop, bool ccw)
         {
            BPolygon polygon;
            if (ccw)
               for (const auto& p : loop) bg::append(polygon.outer(), BPoint(p.x(), p.y()));
            else
               for (auto it = loop.rbegin(); it != loop.rend(); ++it) bg::append(polygon.outer(), BPoint(it->x(), it->y()));
            bg::correct(polygon);
            return polygon;
         };

      BMultiPolygon result;
      for (const auto& loop : outline.solids)
      {
         BMultiPolygon merged;
         bg::union_(result, to_polygon(loop, true), merged);
         result = std::move(merged);
      }
      for (const auto& loop : outline.voids)
      {
         BMultiPolygon cut;
         bg::difference(result, to_polygon(loop, false), cut);
         result = std::move(cut);
      }
      return result;
   }
}

bool IsSimpleOutline(const SectionOutline& outline)
{
   try
   {
      for (const auto* loops : { &outline.solids, &outline.voids })
      {
         for (const auto& loop : *loops)
         {
            if (loop.size() < 3)
               return false;
            BPolygon polygon;
            for (const auto& p : loop)
            {
               if (!std::isfinite(p.x()) || !std::isfinite(p.y()))
                  return false;
               bg::append(polygon.outer(), BPoint(p.x(), p.y()));
            }
            if (bg::intersects(polygon))
               return false; // it crosses itself
         }
      }
      return true;
   }
   catch (...)
   {
      return false;
   }
}

double OutlineDifferenceArea(const SectionOutline& a, const SectionOutline& b)
{
   try
   {
      BMultiPolygon difference;
      bg::sym_difference(region(a), region(b), difference);
      return bg::area(difference);
   }
   catch (...)
   {
      return DBL_MAX;
   }
}

double FittedSectionComparison::Difference(const SectionMeasures& model, const SectionMeasures& other) const
{
   const auto& t = m_Tolerances;
   double d = fabs(model.height - other.height) / t.length;
   d += fabs(model.centroid_y - other.centroid_y) / t.length;
   if (other.area != 0) d += fabs(model.area - other.area) / fabs(other.area) / t.ratio;
   if (other.ix != 0) d += fabs(model.ix - other.ix) / fabs(other.ix) / t.ratio;

   // the area between the outlines: a strip as wide as the length tolerance down both sides gives about 1
   if (0 < model.height)
      d += OutlineDifferenceArea(model.outline, other.outline) / (2 * t.length * model.height);
   return d;
}

std::vector<std::string> FittedSectionComparison::Allowed(const SectionMeasures& model, const SectionMeasures& other) const
{
   return compare_fitted_sections(model, other, m_Tolerances).allowed;
}

std::string DescribeSection(const SectionMeasures& m)
{
   std::ostringstream os;
   os << std::fixed << std::setprecision(3) << "height " << inches(m.height) << " in, width " << inches(m.width) << " in, "
      << std::setprecision(1) << "area " << WBFL::Units::ConvertFromSysUnits(m.area, WBFL::Units::Measure::Inch2) << " in^2, "
      << std::setprecision(0) << "I " << WBFL::Units::ConvertFromSysUnits(m.ix, WBFL::Units::Measure::Inch4) << " in^4";
   return os.str();
}

namespace
{
   // The sections a beam factory makes for sets of dimensions, for fitting. It knows the factory's length dimensions,
   // rejects dimensions the factory doesn't accept or with a negative length, and keeps the fit within a budget
   // (evaluations and time), so a factory whose sections don't converge (or take long to make) can't stall an import
   class FactorySections
   {
   public:
      FactorySections(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, int max_evaluations = 20000, double max_seconds = 30.0) :
         m_Factory(factory), m_Start(start), m_MaxEvaluations(max_evaluations), m_MaxSeconds(max_seconds), m_StartTime(std::chrono::steady_clock::now())
      {
         const auto& names = factory->GetDimensionNames();
         const auto& units = factory->GetDimensionUnits(true);
         for (size_t i = 0; i < start.size(); i++)
         {
            auto found = std::find(names.begin(), names.end(), start[i].first);
            if (found == names.end())
               continue;
            auto unit = units[std::distance(names.begin(), found)];
            if (unit != (const WBFL::Units::Length*)BFDIMUNITSCALAR && unit != (const WBFL::Units::Length*)BFDIMUNITBOOLEAN)
               m_Lengths.push_back(i);
         }
      }

      const GirderDimensions& Start() const { return m_Start; }
      const std::vector<size_t>& Lengths() const { return m_Lengths; }

      bool OutOfBudget() const
      {
         return m_MaxEvaluations <= m_Evaluations || m_MaxSeconds < std::chrono::duration<double>(std::chrono::steady_clock::now() - m_StartTime).count();
      }

      // the section of a set of dimensions, or nullopt if they're rejected or the budget is spent
      std::optional<SectionMeasures> Section(const GirderDimensions& d)
      {
         if (OutOfBudget())
            return std::nullopt;
         m_Evaluations++;

         for (size_t i : m_Lengths)
            if (d[i].second < 0)
               return std::nullopt; // the factory may accept it, but a negative length isn't a girder

         std::_tstring message;
         if (!m_Factory->ValidateDimensions(d, true, &message))
            return std::nullopt;
         try
         {
            CComPtr<IGirderSection> section;
            m_Factory->CreateGirderSection(nullptr, INVALID_ID, d, -1, -1, &section);
            CComQIPtr<IShape> shape(section);
            if (!shape)
               return std::nullopt;
            auto outline = GetShapeOutline(shape);
            if (!IsSimpleOutline(outline))
               return std::nullopt; // the factory accepted the dimensions, but the section isn't a girder's
            return MeasureSection(outline);
         }
         catch (...)
         {
            return std::nullopt;
         }
      }

      // the dimensions with the active ones set from x
      GirderDimensions With(const GirderDimensions& base, const std::vector<size_t>& active, const std::vector<double>& x) const
      {
         auto d = base;
         for (size_t k = 0; k < active.size(); k++)
            d[active[k]].second = x[k];
         return d;
      }

      // The length dimensions that change the section (others, e.g. end blocks, keep their values), and a step size for each
      void ActiveDimensions(const GirderDimensions& base, const SectionMeasures& base_section, std::vector<size_t>& active, std::vector<double>& steps)
      {
         const double min_step = WBFL::Units::ConvertToSysUnits(0.25, WBFL::Units::Measure::Inch);
         for (size_t i : m_Lengths)
         {
            double h = std::max(0.1 * fabs(base[i].second), min_step);
            bool changes = false;
            for (double sign : { 1.0, -1.0 })
            {
               auto d = base;
               d[i].second += sign * h;
               if (auto s = Section(d); s && (1e-9 < fabs(s->area - base_section.area) || 1e-9 < fabs(s->height - base_section.height)))
                  changes = true;
            }
            if (changes)
            {
               active.push_back(i);
               steps.push_back(h);
            }
         }
      }

      // A dimension that changes a feature of the section (e.g. H or the web height for the height) solved for the
      // feature's value. The dimensions that change it the most are tried first, and of those that change it as much
      // the largest first (e.g. the height of an NU girder is the sum of its depths: the web's can reach a height the
      // flange's can't). Returns the dimensions with the first that reaches the value, or the closest, or the
      // dimensions unchanged if none gets the feature closer
      GirderDimensions SolveFeature(const GirderDimensions& base, const SectionMeasures& base_section, const std::function<double(const SectionMeasures&)>& feature, double target)
      {
         std::vector<size_t> active;
         std::vector<double> steps;
         ActiveDimensions(base, base_section, active, steps);

         double f_base = feature(base_section);
         auto feature_of = [&](size_t i, double value) -> std::optional<double>
            {
               auto d = base;
               d[i].second = value;
               auto s = Section(d);
               return s ? std::optional<double>(feature(*s)) : std::nullopt;
            };

         struct Option { size_t dimension; double slope; };
         std::vector<Option> options;
         for (size_t k = 0; k < active.size(); k++)
         {
            size_t i = active[k];
            auto f = feature_of(i, base[i].second + steps[k]);
            if (f && 1e-9 < fabs(*f - f_base))
               options.push_back({ i, fabs(*f - f_base) / steps[k] });
         }
         std::stable_sort(options.begin(), options.end(), [&](const Option& a, const Option& b)
            {
               if (1e-6 < fabs(a.slope - b.slope))
                  return a.slope > b.slope;
               return base[a.dimension].second > base[b.dimension].second;
            });

         GirderDimensions best = base;
         double best_error = fabs(target - f_base);
         for (const auto& option : options)
         {
            if (best_error <= 1e-7)
               break;

            // secant iterations on feature(value) = target
            size_t dimension = option.dimension;
            double x0 = base[dimension].second, f0 = f_base;
            double x1 = x0 + (target - f0) / option.slope; // the slope is about +/- slope; try both signs
            auto f1 = feature_of(dimension, x1);
            if (!f1 || fabs(target - *f1) > fabs(target - f0))
            {
               x1 = x0 - (target - f0) / option.slope;
               f1 = feature_of(dimension, x1);
            }
            for (int iteration = 0; f1 && iteration < 20 && 1e-7 < fabs(target - *f1) && 1e-12 < fabs(*f1 - f0); iteration++)
            {
               double x2 = x1 + (target - *f1) * (x1 - x0) / (*f1 - f0);
               x0 = x1; f0 = *f1;
               x1 = x2; f1 = feature_of(dimension, x1);
            }
            if (f1 && fabs(target - *f1) < best_error)
            {
               best = base;
               best[dimension].second = x1;
               best_error = fabs(target - *f1);
            }
         }
         return best;
      }

      GirderDimensions SolveHeight(const GirderDimensions& base, const SectionMeasures& base_section, double height)
      {
         return SolveFeature(base, base_section, [](const SectionMeasures& m) {return m.height; }, height);
      }

      // the fit of a set of dimensions, judged by the comparison
      std::optional<DimensionFit> Result(const GirderDimensions& d, const SectionMeasures& model, const ISectionComparison& comparison)
      {
         // the result is measured even if the budget is spent
         auto evaluations = m_Evaluations;
         auto start_time = m_StartTime;
         m_Evaluations = 0;
         m_MaxEvaluations = std::max(m_MaxEvaluations, 1);
         m_StartTime = std::chrono::steady_clock::now();
         auto section = Section(d);
         m_Evaluations = evaluations;
         m_StartTime = start_time;
         if (!section)
            return std::nullopt;

         DimensionFit fit;
         fit.dimensions = d;
         fit.section = *section;
         fit.difference = comparison.Difference(model, fit.section);
         fit.differences = comparison.Compare(model, fit.section);
         fit.allowed = comparison.Allowed(model, fit.section);
         return fit;
      }

   private:
      std::shared_ptr<PGS::Beams::BeamFactory> m_Factory;
      GirderDimensions m_Start;
      std::vector<size_t> m_Lengths;
      int m_MaxEvaluations;
      double m_MaxSeconds;
      int m_Evaluations = 0;
      std::chrono::steady_clock::time_point m_StartTime;
   };

   // the average width in each AREA_PROFILE_STEP band below the top, from the exact area profile (so it changes
   // continuously with the dimensions, also when a flange edge moves), and zero below the bottom
   double band_width(const SectionMeasures& m, size_t k)
   {
      double above = (k == 0 ? 0.0 : (k - 1 < m.area_profile.size() ? m.area_profile[k - 1] : m.area));
      double below = (k < m.area_profile.size() ? m.area_profile[k] : m.area);
      return (below - above) / AREA_PROFILE_STEP;
   }

   // the residuals of a section, or nullopt if they can't be computed. Their number must not change during a fit
   using Residuals = std::function<std::optional<Eigen::VectorXd>(const SectionMeasures&)>;

   std::vector<double> to_vector(const Eigen::VectorXd& v)
   {
      return std::vector<double>(v.data(), v.data() + v.size());
   }

   // Damped least squares (Levenberg-Marquardt) of the active dimensions, starting at x (their values). A step that
   // would make a length negative stops at zero, so a dimension that is zero (e.g. a chamfer or a flange's horizontal
   // part, which many factory defaults have) can stay zero or grow; without the bound every step from such a start is
   // rejected and the fit never moves. Returns the values of the active dimensions
   Eigen::VectorXd LeastSquares(FactorySections& sections, const GirderDimensions& dimensions, const std::vector<size_t>& active, Eigen::VectorXd x, const Residuals& residuals, int max_iterations)
   {
      size_t n = active.size();
      auto evaluate = [&](const Eigen::VectorXd& v) -> std::optional<Eigen::VectorXd>
         {
            auto s = sections.Section(sections.With(dimensions, active, to_vector(v)));
            return s ? residuals(*s) : std::nullopt;
         };

      auto r0 = evaluate(x);
      if (!r0)
         return x;
      Eigen::VectorXd r = *r0;
      double cost = r.squaredNorm();
      double lambda = 1e-2;
      for (int iteration = 0; iteration < max_iterations && !sections.OutOfBudget(); iteration++)
      {
         // the Jacobian by forward differences (backward at a bound; zero if neither works)
         Eigen::MatrixXd J(r.size(), n);
         for (size_t k = 0; k < n; k++)
         {
            double h = std::max(1e-3 * fabs(x(k)), 1e-5);
            Eigen::VectorXd xk = x;
            xk(k) += h;
            auto rk = evaluate(xk);
            if (rk && rk->size() == r.size())
            {
               J.col(k) = (*rk - r) / h;
               continue;
            }
            xk(k) = x(k) - h;
            rk = evaluate(xk);
            if (rk && rk->size() == r.size())
               J.col(k) = (r - *rk) / h;
            else
               J.col(k).setZero();
         }

         // damped Gauss-Newton steps, bounded at zero, more damped until one is better
         Eigen::MatrixXd A = J.transpose() * J;
         Eigen::VectorXd g = J.transpose() * r;
         bool improved = false;
         for (int attempt = 0; attempt < 12 && !improved && !sections.OutOfBudget(); attempt++)
         {
            Eigen::MatrixXd M = A;
            for (size_t k = 0; k < n; k++)
               M(k, k) += lambda * std::max(A(k, k), 1e-12);
            Eigen::VectorXd candidate = (x + M.ldlt().solve(-g)).cwiseMax(0.0);
            auto rc = evaluate(candidate);
            if (rc && rc->size() == r.size() && rc->squaredNorm() < cost)
            {
               double step = (candidate - x).norm();
               double previous = cost;
               x = candidate;
               r = *rc;
               cost = r.squaredNorm();
               lambda = std::max(lambda / 3, 1e-9);
               improved = true;
               if (previous - cost < 1e-10 * previous || step < 1e-7)
                  iteration = max_iterations; // converged
            }
            else
            {
               lambda *= 4;
            }
         }
         if (!improved)
            break;
      }
      return x;
   }

   // The edges of an outline on the left of its vertical axis (x <= 0), without the edges along the axis. Girder
   // sections are symmetric about the axis, so the left half is the section, and two halves can be compared without
   // pairing each edge with its mirror image
   using Edge = std::pair<Eigen::Vector2d, Eigen::Vector2d>;
   std::vector<Edge> LeftHalfEdges(const SectionOutline& outline)
   {
      const double on_axis = 1e-9;
      std::vector<Edge> edges;
      auto add = [&](const Loop& loop)
         {
            for (size_t i = 0; i < loop.size(); i++)
            {
               Eigen::Vector2d a = loop[i], b = loop[(i + 1) % loop.size()];
               if (on_axis < a.x() && on_axis < b.x())
                  continue;
               if (on_axis < a.x())
                  a = a + (b - a) * (-a.x() / (b.x() - a.x()));
               else if (on_axis < b.x())
                  b = a + (b - a) * (-a.x() / (b.x() - a.x()));
               if ((fabs(a.x()) <= on_axis && fabs(b.x()) <= on_axis) || (b - a).norm() < 1e-12)
                  continue;
               edges.emplace_back(a, b);
            }
         };
      for (const auto& loop : outline.solids)
         add(loop);
      for (const auto& loop : outline.voids)
         add(loop);
      return edges;
   }

   // n points spaced evenly along edges
   std::vector<Eigen::Vector2d> SampleEdges(const std::vector<Edge>& edges, size_t n)
   {
      std::vector<double> lengths{ 0.0 };
      for (const auto& [a, b] : edges)
         lengths.push_back(lengths.back() + (b - a).norm());
      std::vector<Eigen::Vector2d> points;
      if (edges.empty() || lengths.back() <= 0)
         return points;
      for (size_t k = 0; k < n; k++)
      {
         double s = lengths.back() * (k + 0.5) / n;
         size_t i = std::min((size_t)std::distance(lengths.begin(), std::upper_bound(lengths.begin(), lengths.end(), s)) - 1, edges.size() - 1);
         const auto& [a, b] = edges[i];
         double length = lengths[i + 1] - lengths[i];
         points.push_back(length <= 0 ? a : Eigen::Vector2d(a + (b - a) * ((s - lengths[i]) / length)));
      }
      return points;
   }

   double DistanceToEdges(const Eigen::Vector2d& p, const std::vector<Edge>& edges)
   {
      double distance = DBL_MAX;
      for (const auto& [a, b] : edges)
      {
         Eigen::Vector2d ab = b - a;
         double t = std::clamp((p - a).dot(ab) / std::max(ab.squaredNorm(), 1e-18), 0.0, 1.0);
         distance = std::min(distance, (p - (a + t * ab)).norm());
      }
      return distance;
   }

   // Residuals for the distance between the left halves of two outlines, each about 1 at the length tolerance: points
   // spaced evenly along the model's outline measured to the other outline, points along the other measured to the
   // model's, and the height. A point is measured to the nearest edge at any depth, so a fit isn't held where an edge
   // of one outline lines up with the wrong edge of the other (as it is when widths are compared depth by depth)
   class OutlineResiduals
   {
   public:
      OutlineResiduals(const SectionMeasures& model, double tolerance, size_t samples = 150) :
         m_Height(model.height), m_Edges(LeftHalfEdges(model.outline)), m_Samples(samples), m_Tolerance(tolerance)
      {
         m_Points = SampleEdges(m_Edges, samples);
      }

      std::optional<Eigen::VectorXd> operator()(const SectionMeasures& s) const
      {
         auto edges = LeftHalfEdges(s.outline);
         auto points = SampleEdges(edges, m_Samples);
         if (m_Points.size() != m_Samples || points.size() != m_Samples)
            return std::nullopt;
         double weight = 1.0 / sqrt((double)m_Samples) / m_Tolerance;
         Eigen::VectorXd r(1 + 2 * m_Samples);
         r(0) = (s.height - m_Height) / m_Tolerance;
         for (size_t k = 0; k < m_Samples; k++)
         {
            r(1 + k) = weight * DistanceToEdges(m_Points[k], edges);
            r(1 + m_Samples + k) = weight * DistanceToEdges(points[k], m_Edges);
         }
         return r;
      }

   private:
      double m_Height;
      std::vector<Edge> m_Edges;
      std::vector<Eigen::Vector2d> m_Points;
      size_t m_Samples;
      double m_Tolerance;
   };

   // The left half of a single-loop outline as a chain of points from the top of the axis down to the bottom of the
   // axis. nullopt if the outline isn't one loop (e.g. it has voids) or doesn't cross the axis
   std::optional<std::vector<Eigen::Vector2d>> LeftChain(const SectionOutline& outline)
   {
      if (outline.solids.size() != 1 || !outline.voids.empty())
         return std::nullopt;
      const auto& loop = outline.solids.front();

      // the loop clipped to x <= 0, without repeated points
      const double on_axis = 1e-7;
      std::vector<Eigen::Vector2d> clipped;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& a = loop[i];
         const auto& b = loop[(i + 1) % loop.size()];
         bool a_in = a.x() <= on_axis, b_in = b.x() <= on_axis;
         if (a_in)
            clipped.push_back(a);
         if (a_in != b_in)
            clipped.emplace_back(0.0, a.y() + (b.y() - a.y()) * (-a.x() / (b.x() - a.x())));
      }
      std::vector<Eigen::Vector2d> c;
      for (const auto& q : clipped)
         if (c.empty() || on_axis < (q - c.back()).norm())
            c.push_back(q);
      while (1 < c.size() && (c.front() - c.back()).norm() <= on_axis)
         c.pop_back();

      size_t n = c.size();
      size_t top = INVALID_INDEX, bottom = INVALID_INDEX;
      for (size_t i = 0; i < n; i++)
      {
         if (on_axis < fabs(c[i].x()))
            continue;
         if (top == INVALID_INDEX || c[top].y() < c[i].y())
            top = i;
         if (bottom == INVALID_INDEX || c[i].y() < c[bottom].y())
            bottom = i;
      }
      if (n < 3 || top == INVALID_INDEX || top == bottom)
         return std::nullopt;

      // from the top, the way that leaves the axis (not down the axis), to the bottom
      int direction = (c[(top + n - 1) % n].x() < c[(top + 1) % n].x() ? -1 : 1);
      std::vector<Eigen::Vector2d> chain;
      for (size_t k = 0, i = top; k < n; k++, i = (i + n + direction) % n)
      {
         chain.push_back(c[i]);
         if (i == bottom)
            return chain;
      }
      return std::nullopt;
   }

   // the chain without the points within the tolerance of a line between their neighbors (Douglas-Peucker): its corners
   std::vector<Eigen::Vector2d> Corners(const std::vector<Eigen::Vector2d>& chain, double tolerance)
   {
      if (chain.size() < 3)
         return chain;
      std::vector<bool> keep(chain.size(), false);
      keep.front() = keep.back() = true;
      std::vector<std::pair<size_t, size_t>> spans{ {0, chain.size() - 1} };
      while (!spans.empty())
      {
         auto [first, last] = spans.back();
         spans.pop_back();
         double farthest = -1;
         size_t index = first;
         for (size_t i = first + 1; i < last; i++)
         {
            double d = DistanceToEdges(chain[i], { {chain[first], chain[last]} });
            if (farthest < d)
            {
               farthest = d;
               index = i;
            }
         }
         if (tolerance < farthest)
         {
            keep[index] = true;
            spans.emplace_back(first, index);
            spans.emplace_back(index, last);
         }
      }
      std::vector<Eigen::Vector2d> corners;
      for (size_t i = 0; i < chain.size(); i++)
         if (keep[i])
            corners.push_back(chain[i]);
      return corners;
   }

   // Pairs of corners of two chains, in order along the chains (dynamic time warping): each corner has at least one
   // partner, and a corner can have several (e.g. a factory's chamfer that is zero, which is two corners at one place)
   std::vector<std::pair<size_t, size_t>> MatchCorners(const std::vector<Eigen::Vector2d>& a, const std::vector<Eigen::Vector2d>& b)
   {
      size_t n = a.size(), m = b.size();
      Eigen::MatrixXd D = Eigen::MatrixXd::Constant(n + 1, m + 1, DBL_MAX);
      D(0, 0) = 0;
      for (size_t i = 1; i <= n; i++)
         for (size_t j = 1; j <= m; j++)
            D(i, j) = (a[i - 1] - b[j - 1]).norm() + std::min({ D(i - 1, j - 1), D(i - 1, j), D(i, j - 1) });
      std::vector<std::pair<size_t, size_t>> pairs;
      size_t i = n, j = m;
      while (0 < i && 0 < j)
      {
         pairs.emplace_back(i - 1, j - 1);
         double diagonal = D(i - 1, j - 1), up = D(i - 1, j), left = D(i, j - 1);
         if (diagonal <= up && diagonal <= left) { i--; j--; }
         else if (up <= left) { i--; }
         else { j--; }
      }
      std::reverse(pairs.begin(), pairs.end());
      return pairs;
   }

   // Residuals for the corners of the left halves of two single-loop outlines: each corner of the model's outline to the
   // farthest of the other's corners paired with it (MatchCorners), about 1 at the length tolerance, with a part of the
   // outline residuals. The corners of a polygonal section (e.g. a WF girder's flange tips and tapers) fix its
   // dimensions directly, so they bring a fit from a factory's defaults to the right shape
   class CornerResiduals
   {
   public:
      CornerResiduals(const SectionMeasures& model, double tolerance, const OutlineResiduals& outline) :
         m_Tolerance(tolerance), m_Outline(outline)
      {
         if (auto chain = LeftChain(model.outline))
            m_Corners = Corners(*chain, tolerance);
      }

      bool Applies() const { return 2 < m_Corners.size(); }

      std::optional<Eigen::VectorXd> operator()(const SectionMeasures& s) const
      {
         auto chain = LeftChain(s.outline);
         auto outline = m_Outline(s);
         if (!Applies() || !chain || !outline)
            return std::nullopt;
         auto corners = Corners(*chain, 1e-9);
         std::vector<Eigen::Vector2d> offsets(m_Corners.size(), Eigen::Vector2d::Zero());
         for (auto [i, j] : MatchCorners(m_Corners, corners))
         {
            Eigen::Vector2d offset = corners[j] - m_Corners[i];
            if (offsets[i].norm() < offset.norm())
               offsets[i] = offset;
         }
         const double outline_weight = 0.3;
         Eigen::VectorXd r(2 * m_Corners.size() + outline->size());
         for (size_t i = 0; i < m_Corners.size(); i++)
         {
            r(2 * i) = offsets[i].x() / m_Tolerance;
            r(2 * i + 1) = offsets[i].y() / m_Tolerance;
         }
         r.tail(outline->size()) = outline_weight * *outline;
         return r;
      }

   private:
      double m_Tolerance;
      const OutlineResiduals& m_Outline;
      std::vector<Eigen::Vector2d> m_Corners;
   };
}

bool FlangesMatch(const SectionMeasures& model, const SectionMeasures& other, const SectionTolerances& tolerances, double fraction)
{
   // the top of each section, measured down from the top, and the bottom, measured up from the bottom
   size_t n = (size_t)(fraction * std::min(model.width_profile.size(), other.width_profile.size()));
   if (n == 0)
      return false;

   int band = (int)ceil(tolerances.length / WIDTH_PROFILE_STEP - 1e-9);
   auto matches = [&](auto width_model, auto width_other)
      {
         for (size_t k = 0; k < n; k++)
         {
            bool match = false;
            for (int j = std::max(0, (int)k - band); !match && j <= std::min((int)n - 1, (int)k + band); j++)
               match = fabs(width_model(k) - width_other(j)) <= tolerances.length || fabs(width_model(j) - width_other(k)) <= tolerances.length;
            if (!match)
               return false;
         }
         return true;
      };

   const auto& wm = model.width_profile;
   const auto& wo = other.width_profile;
   bool top = matches([&](size_t k) {return wm[k]; }, [&](size_t k) {return wo[k]; });
   bool bottom = matches([&](size_t k) {return wm[wm.size() - 1 - k]; }, [&](size_t k) {return wo[wo.size() - 1 - k]; });
   return top && bottom;
}

SectionOutline::Loop SimplifyLoop(const SectionOutline::Loop& loop, double tolerance)
{
   if (loop.size() < 4)
      return loop;
   auto chain = loop;
   chain.push_back(loop.front()); // closed: the first point is kept, and is the last too
   auto corners = Corners(chain, tolerance);
   corners.pop_back();
   return corners.size() < 3 ? loop : corners;
}

double OutlineDistance(const SectionOutline& a, const SectionOutline& b)
{
   auto edges = [](const SectionOutline& outline)
      {
         std::vector<Edge> edges;
         for (const auto* loops : { &outline.solids, &outline.voids })
            for (const auto& loop : *loops)
               for (size_t i = 0; i < loop.size(); i++)
                  edges.emplace_back(loop[i], loop[(i + 1) % loop.size()]);
         return edges;
      };
   auto points = [](const std::vector<Edge>& edges)
      {
         // the corners, and points spaced evenly between them
         auto points = SampleEdges(edges, 1000);
         for (const auto& edge : edges)
            points.push_back(edge.first);
         return points;
      };
   auto edges_a = edges(a), edges_b = edges(b);
   if (edges_a.empty() || edges_b.empty())
      return DBL_MAX;
   double distance = 0;
   for (const auto& p : points(edges_a))
      distance = std::max(distance, DistanceToEdges(p, edges_b));
   for (const auto& p : points(edges_b))
      distance = std::max(distance, DistanceToEdges(p, edges_a));
   return distance;
}

std::optional<DimensionFit> HeightDimensionFitter::Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const
{
   if (!factory)
      return std::nullopt;
   FactorySections sections(factory, start);
   auto start_section = sections.Section(start);
   if (!start_section)
      return std::nullopt;
   return sections.Result(sections.SolveHeight(start, *start_section, model.height), model, comparison);
}

std::optional<DimensionFit> LevenbergMarquardtDimensionFitter::Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const
{
   if (!factory)
      return std::nullopt;
   FactorySections sections(factory, start);
   auto start_section = sections.Section(start);
   if (!start_section)
      return std::nullopt;

   // start at the model's height
   GirderDimensions dimensions = sections.SolveHeight(start, *start_section, model.height);
   auto section = sections.Section(dimensions);
   if (!section)
      return std::nullopt;

   std::vector<size_t> active;
   std::vector<double> steps;
   sections.ActiveDimensions(dimensions, *section, active, steps);
   if (active.empty())
      return sections.Result(dimensions, model, comparison);

   // The residuals, each about 1 at the tolerances: the average width in each band of depth, then the height,
   // centroid, area, and I. The bands are weighted so all of them together count about as much as one of the others
   SectionTolerances t;
   size_t bands = model.area_profile.size() + 1;
   double band_weight = 1.0 / sqrt((double)bands);
   Residuals residuals = [&](const SectionMeasures& s) -> std::optional<Eigen::VectorXd>
      {
         Eigen::VectorXd r(bands + 4);
         for (size_t k = 0; k < bands; k++)
            r(k) = band_weight * (band_width(model, k) - band_width(s, k)) / t.length;
         r(bands) = (model.height - s.height) / t.length;
         r(bands + 1) = (model.centroid_y - s.centroid_y) / t.length;
         r(bands + 2) = s.area != 0 ? (model.area - s.area) / fabs(model.area) / t.ratio : 0;
         r(bands + 3) = s.ix != 0 ? (model.ix - s.ix) / fabs(model.ix) / t.ratio : 0;
         return r;
      };

   Eigen::VectorXd x(active.size());
   for (size_t k = 0; k < active.size(); k++)
      x(k) = dimensions[active[k]].second;
   x = LeastSquares(sections, dimensions, active, x, residuals, m_MaxIterations);
   return sections.Result(sections.With(dimensions, active, to_vector(x)), model, comparison);
}

std::optional<DimensionFit> OutlineDimensionFitter::Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const
{
   if (!factory)
      return std::nullopt;
   FactorySections sections(factory, start);
   auto start_section = sections.Section(start);
   if (!start_section)
      return std::nullopt;

   // start at the model's height
   GirderDimensions dimensions = sections.SolveHeight(start, *start_section, model.height);
   auto section = sections.Section(dimensions);
   if (!section)
      return std::nullopt;

   // The dimensions with every zero length (e.g. a chamfer, or a flange's horizontal part) at a small length where the
   // factory accepts it and the section changes, so the section has all of its features. A factory's outline leaves out
   // the corners of a zero dimension, and one zero dimension can hide another (e.g. a web's lower part with no
   // horizontal offset is just more web, so its depth changes nothing): the active dimensions are found with all of
   // the features, and the corner fit starts with them
   const double start_length = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
   GirderDimensions featured = dimensions;
   auto featured_section = section;
   for (size_t i : sections.Lengths())
   {
      if (featured[i].second != 0)
         continue;
      auto d = featured;
      d[i].second = start_length;
      auto s = sections.Section(d);
      if (s && (1e-9 < fabs(s->area - featured_section->area) || 1e-9 < fabs(s->height - featured_section->height)))
      {
         featured = d;
         featured_section = s;
      }
   }

   std::vector<size_t> active;
   std::vector<double> steps;
   sections.ActiveDimensions(featured, *featured_section, active, steps);
   if (active.empty())
      return sections.Result(dimensions, model, comparison);

   Eigen::VectorXd x0(active.size()), x_featured(active.size());
   for (size_t k = 0; k < active.size(); k++)
   {
      x0(k) = dimensions[active[k]].second;
      x_featured(k) = featured[active[k]].second;
   }

   std::optional<DimensionFit> best;
   auto keep = [&](const Eigen::VectorXd& x)
      {
         auto fit = sections.Result(sections.With(dimensions, active, to_vector(x)), model, comparison);
         if (fit && (!best || fit->difference < best->difference))
            best = fit;
      };

   SectionTolerances t;
   OutlineResiduals outline(model, t.length);
   Residuals outline_residuals = [&outline](const SectionMeasures& s) { return outline(s); };

   // the corners first, then the outline: for polygonal sections
   CornerResiduals corners(model, t.length, outline);
   if (corners.Applies() && LeftChain(section->outline))
   {
      // with all of the features: a corner that isn't there can't be paired with the model's
      Eigen::VectorXd x_corners = x_featured;
      Residuals corner_residuals = [&corners](const SectionMeasures& s) { return corners(s); };
      auto x = LeastSquares(sections, dimensions, active, x_corners, corner_residuals, m_MaxIterations);
      keep(LeastSquares(sections, dimensions, active, x, outline_residuals, m_MaxIterations));
   }

   // the outline only: for curved sections (e.g. fillets), whose corners are many short edges that pair poorly
   keep(LeastSquares(sections, dimensions, active, x0, outline_residuals, m_MaxIterations));
   return best;
}

std::optional<DimensionFit> SimplexDimensionFitter::Fit(std::shared_ptr<PGS::Beams::BeamFactory> factory, const GirderDimensions& start, const SectionMeasures& model, const ISectionComparison& comparison) const
{
   if (!factory)
      return std::nullopt;
   FactorySections sections(factory, start);
   auto start_section = sections.Section(start);
   if (!start_section)
      return std::nullopt;

   GirderDimensions dimensions = sections.SolveHeight(start, *start_section, model.height);
   auto section = sections.Section(dimensions);
   if (!section)
      return std::nullopt;

   std::vector<size_t> active;
   std::vector<double> steps;
   sections.ActiveDimensions(dimensions, *section, active, steps);

   const double rejected = 1e12;
   auto objective = [&](const std::vector<double>& x)
      {
         auto s = sections.Section(sections.With(dimensions, active, x));
         return s ? comparison.Difference(model, *s) : rejected;
      };

   // downhill simplex (Nelder-Mead), restarted from the best point
   size_t n = active.size();
   std::vector<double> best(n);
   for (size_t k = 0; k < n; k++)
      best[k] = dimensions[active[k]].second;
   double best_value = comparison.Difference(model, *section);

   for (int restart = 0; restart < 4 && 0 < n && !sections.OutOfBudget(); restart++)
   {
      std::vector<std::vector<double>> simplex(n + 1, best);
      std::vector<double> values(n + 1);
      for (size_t k = 0; k < n; k++)
         simplex[k + 1][k] += steps[k];
      for (size_t j = 0; j <= n; j++)
         values[j] = objective(simplex[j]);

      for (int iteration = 0; iteration < 4000 && !sections.OutOfBudget(); iteration++)
      {
         std::vector<size_t> order(n + 1);
         std::iota(order.begin(), order.end(), 0);
         std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {return values[a] < values[b]; });
         size_t lo = order.front(), hi = order.back(), next_hi = order[n - 1];
         if (fabs(values[hi] - values[lo]) < 1e-6)
            break;

         std::vector<double> centroid(n, 0.0);
         for (size_t j = 0; j <= n; j++)
            if (j != hi)
               for (size_t k = 0; k < n; k++)
                  centroid[k] += simplex[j][k] / n;

         auto along = [&](double t)
            {
               std::vector<double> x(n);
               for (size_t k = 0; k < n; k++)
                  x[k] = centroid[k] + t * (simplex[hi][k] - centroid[k]);
               return x;
            };

         auto reflected = along(-1.0);
         double fr = objective(reflected);
         if (fr < values[lo])
         {
            auto expanded = along(-2.0);
            double fe = objective(expanded);
            if (fe < fr) { simplex[hi] = expanded; values[hi] = fe; }
            else { simplex[hi] = reflected; values[hi] = fr; }
         }
         else if (fr < values[next_hi])
         {
            simplex[hi] = reflected; values[hi] = fr;
         }
         else
         {
            auto contracted = along(0.5);
            double fc = objective(contracted);
            if (fc < values[hi])
            {
               simplex[hi] = contracted; values[hi] = fc;
            }
            else
            {
               // shrink toward the best point
               for (size_t j = 0; j <= n; j++)
               {
                  if (j == lo)
                     continue;
                  for (size_t k = 0; k < n; k++)
                     simplex[j][k] = simplex[lo][k] + 0.5 * (simplex[j][k] - simplex[lo][k]);
                  values[j] = objective(simplex[j]);
               }
            }
         }
      }

      size_t lo = std::distance(values.begin(), std::min_element(values.begin(), values.end()));
      if (values[lo] < best_value)
      {
         best_value = values[lo];
         best = simplex[lo];
      }
   }

   // then a pattern search (Hooke-Jeeves), one dimension at a time with halving steps
   {
      std::vector<double> h(steps);
      const double min_h = 1e-5; // 0.01 mm
      for (int round = 0; round < 2000 && !sections.OutOfBudget() && 0 < n && min_h < *std::max_element(h.begin(), h.end()); round++)
      {
         bool improved = false;
         for (size_t k = 0; k < n; k++)
         {
            for (double sign : { 1.0, -1.0 })
            {
               auto x = best;
               x[k] += sign * h[k];
               double f = objective(x);
               if (f < best_value - 1e-9)
               {
                  best = x;
                  best_value = f;
                  improved = true;
                  break;
               }
            }
         }
         if (!improved)
            for (auto& hk : h)
               hk /= 2;
      }
   }

   return sections.Result(sections.With(dimensions, active, best), model, comparison);
}
