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
#include <iomanip>
#include <map>
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

std::optional<SectionOutline> CutGirderSection(const Mesh& mesh, const Eigen::Vector2d& start, const Eigen::Vector2d& end, double fraction)
{
   Eigen::Vector2d axis = end - start;
   if (axis.norm() <= 0)
      return std::nullopt;
   axis.normalize();

   Eigen::Vector2d point2d = start + fraction * (end - start);
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

   if (loops.empty())
      return std::nullopt;

   // the largest loop is the outside of the section, loops inside it are voids
   std::sort(loops.begin(), loops.end(), [](const Loop& a, const Loop& b) {return fabs(signed_area(a)) > fabs(signed_area(b)); });

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
   if (m.area <= 0)
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
   return m;
}

std::vector<std::string> CompareSections(const SectionMeasures& model, const SectionMeasures& library, const SectionTolerances& tolerances)
{
   std::vector<std::string> differences;
   auto length = [&](const char* name, double a, double b)
      {
         if (tolerances.length < fabs(a - b))
         {
            std::ostringstream os;
            os << std::fixed << std::setprecision(3) << name << " " << inches(a) << " in vs " << inches(b) << " in";
            differences.push_back(os.str());
         }
      };
   auto ratio = [&](const char* name, double a, double b, double unit, const char* unit_name)
      {
         if (b != 0 && tolerances.ratio < fabs(a - b) / fabs(b))
         {
            std::ostringstream os;
            os << std::fixed << std::setprecision(1) << name << " " << a / unit << " vs " << b / unit << " " << unit_name;
            differences.push_back(os.str());
         }
      };

   length("height", model.height, library.height);
   ratio("area", model.area, library.area, WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch2), "in^2");
   ratio("I", model.ix, library.ix, WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch4), "in^4");

   // Width at each depth. A width may come from a depth within the tolerance, so a flange thickness that differs by
   // less than the tolerance isn't a difference
   int band = (int)ceil(tolerances.length / WIDTH_PROFILE_STEP - 1e-9);
   size_t n = std::min(model.width_profile.size(), library.width_profile.size());
   for (size_t k = 0; k < n; k++)
   {
      bool match = false;
      for (int j = std::max(0, (int)k - band); !match && j <= std::min((int)n - 1, (int)k + band); j++)
         match = fabs(model.width_profile[k] - library.width_profile[j]) <= tolerances.length;
      if (!match)
      {
         std::ostringstream os;
         os << std::fixed << std::setprecision(3) << "width " << inches((k + 0.5) * WIDTH_PROFILE_STEP) << " in below the top "
            << inches(model.width_profile[k]) << " in vs " << inches(library.width_profile[k]) << " in";
         differences.push_back(os.str());
         break; // the first one says enough
      }
   }
   return differences;
}

double SectionDifference(const SectionMeasures& model, const SectionMeasures& library)
{
   SectionTolerances t;
   double d = fabs(model.height - library.height) / t.length;
   if (library.area != 0) d += fabs(model.area - library.area) / fabs(library.area) / t.ratio;
   if (library.ix != 0) d += fabs(model.ix - library.ix) / fabs(library.ix) / t.ratio;
   size_t n = std::min(model.width_profile.size(), library.width_profile.size());
   double w = 0;
   for (size_t k = 0; k < n; k++)
      w += fabs(model.width_profile[k] - library.width_profile[k]);
   if (0 < n)
      d += w / n / t.length;
   return d;
}

std::string DescribeSection(const SectionMeasures& m)
{
   std::ostringstream os;
   os << std::fixed << std::setprecision(3) << "height " << inches(m.height) << " in, width " << inches(m.width) << " in, "
      << std::setprecision(1) << "area " << WBFL::Units::ConvertFromSysUnits(m.area, WBFL::Units::Measure::Inch2) << " in^2, "
      << std::setprecision(0) << "I " << WBFL::Units::ConvertFromSysUnits(m.ix, WBFL::Units::Measure::Inch4) << " in^4";
   return os.str();
}
