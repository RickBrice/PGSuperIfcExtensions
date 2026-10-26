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
#include "GirderReinforcement.h"
#include "IfcImporter.h"
#include "IfcImportUnits.h"

#include <ifcgeom/abstract_mapping.h>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

const Float64 reinforcement_position_tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 8.0, WBFL::Units::Measure::Inch);
const Float64 end_distance_limit = WBFL::Units::ConvertToSysUnits(2.0, WBFL::Units::Measure::Feet);

namespace
{
   // The girder coordinate system of a girder: the origin is the top of the CL girder at the start, the x axis is along
   // the top line to the top of the CL girder at the end, and the y axis is up, normal to the x axis
   struct GirderFrame
   {
      Eigen::Vector3d origin, x_axis, left_axis, up_axis;
      Float64 plan_factor; // plan distance per distance along the top line
      Float64 length;      // in plan

      GirderFrame(const GirderCenterline& girder)
      {
         origin = Eigen::Vector3d(girder.start.x(), girder.start.y(), girder.start_top);
         Eigen::Vector3d end(girder.end.x(), girder.end.y(), girder.end_top);
         x_axis = (end - origin).normalized();
         left_axis = Eigen::Vector3d::UnitZ().cross(x_axis).normalized();
         up_axis = x_axis.cross(left_axis);
         length = (girder.end - girder.start).norm();
         plan_factor = length / (end - origin).norm();
      }

      GirderPoint ToGirder(const Eigen::Vector3d& p) const
      {
         Eigen::Vector3d d = p - origin;
         return { d.dot(x_axis) * plan_factor, -d.dot(left_axis), -d.dot(up_axis) };
      }
   };

   // The 4x4 transformation of a placement (object placement, axis placement, or cartesian transformation operator), in metres
   template <typename Placement>
   Eigen::Matrix4d placement_matrix(ifcopenshell::geom::abstract_mapping& mapping, Placement placement)
   {
      if (!placement)
         return Eigen::Matrix4d::Identity();
      auto matrix = ifcopenshell::geom::taxonomy::dcast<ifcopenshell::geom::taxonomy::matrix4>(mapping.map(placement));
      return matrix ? matrix->ccomponents() : Eigen::Matrix4d::Identity();
   }

   // A straight piece of a centerline (global coordinates, system units)
   struct Segment
   {
      Eigen::Vector3d a, b;
      Float64 radius = 0; // of the strand or bar, as drawn (system units)
   };

   Eigen::Vector3d transform(const Eigen::Matrix4d& m, const Eigen::Vector3d& p)
   {
      return (m * Eigen::Vector4d(p.x(), p.y(), p.z(), 1.0)).head<3>();
   }

   Float64 to_system_length(Float64 value)
   {
      return ConvertToSysUnits<IfcSchema::IfcLengthMeasure>(CIfcImporter::GetUnits(), value);
   }

   Eigen::Vector3d to_system_point(const std::vector<double>& c)
   {
      return Eigen::Vector3d(to_system_length(c[0]), to_system_length(1 < c.size() ? c[1] : 0.0), to_system_length(2 < c.size() ? c[2] : 0.0));
   }

   // The center of a circular profile (in profile coordinates, system units) and its radius, or nullopt if the profile isn't a
   // circle
   std::optional<Eigen::Vector2d> circle_center(IfcSchema::IfcProfileDef profile, Float64& radius)
   {
      auto center = [](IfcSchema::IfcAxis2Placement2D position) -> Eigen::Vector2d
         {
            if (!position)
               return Eigen::Vector2d::Zero();
            auto c = position.Location().as<IfcSchema::IfcCartesianPoint>().Coordinates();
            return Eigen::Vector2d(to_system_length(c[0]), to_system_length(c[1]));
         };

      if (auto circle = profile.as<IfcSchema::IfcCircleProfileDef>())
      {
         radius = to_system_length(circle.Radius());
         return center(circle.Position());
      }

      if (auto arbitrary = profile.as<IfcSchema::IfcArbitraryClosedProfileDef>())
      {
         auto conic = arbitrary.OuterCurve().as<IfcSchema::IfcConic>();
         auto ellipse = conic.as<IfcSchema::IfcEllipse>();
         auto circle = conic.as<IfcSchema::IfcCircle>();
         if (conic && (circle || (ellipse && IsEqual(ellipse.SemiAxis1(), ellipse.SemiAxis2()))))
         {
            radius = to_system_length(circle ? circle.Radius() : ellipse.SemiAxis1());
            return center(conic.Position().as<IfcSchema::IfcAxis2Placement2D>());
         }
      }
      return std::nullopt;
   }

   // The straight pieces of a swept disk's directrix: the lines of a polyline, or the line segments of an indexed poly curve
   // (its arcs are bends, left out)
   void collect_directrix_segments(IfcSchema::IfcCurve directrix, const Eigen::Matrix4d& m, Float64 radius, std::vector<Segment>& segments)
   {
      std::vector<Eigen::Vector3d> points;
      if (auto polyline = directrix.as<IfcSchema::IfcPolyline>())
      {
         for (const auto& point : polyline.Points())
            points.push_back(transform(m, to_system_point(point.Coordinates())));
         for (size_t i = 1; i < points.size(); i++)
            segments.push_back({ points[i - 1], points[i], radius });
      }
      else if (auto curve = directrix.as<IfcSchema::IfcIndexedPolyCurve>())
      {
         if (auto list3d = curve.Points().as<IfcSchema::IfcCartesianPointList3D>())
         {
            for (const auto& c : list3d.CoordList())
               points.push_back(transform(m, to_system_point(c)));
         }
         else if (auto list2d = curve.Points().as<IfcSchema::IfcCartesianPointList2D>())
         {
            for (const auto& c : list2d.CoordList())
               points.push_back(transform(m, to_system_point(c)));
         }

         if (auto curve_segments = curve.Segments())
         {
            for (const auto& segment : *curve_segments)
            {
               if (auto line = segment.as<IfcSchema::IfcLineIndex>())
               {
                  std::vector<int64_t> indices = line;
                  for (size_t i = 1; i < indices.size(); i++)
                  {
                     auto a = indices[i - 1] - 1, b = indices[i] - 1; // IFC indices are 1-based
                     if (0 <= a && 0 <= b && a < (int64_t)points.size() && b < (int64_t)points.size())
                        segments.push_back({ points[a], points[b], radius });
                  }
               }
            }
         }
         else
         {
            for (size_t i = 1; i < points.size(); i++)
               segments.push_back({ points[i - 1], points[i], radius });
         }
      }
   }

   // The straight pieces of a centerline in a representation item (m: from the item's coordinates to global): the axes of
   // circular profiles extruded along a straight line, and the straight pieces of swept disks. Revolved solids of circular
   // profiles are bends; their centerlines go in arcs, as short straight pieces. Mapped items are followed
   void collect_segments(ifcopenshell::geom::abstract_mapping& mapping, IfcSchema::IfcRepresentationItem item, const Eigen::Matrix4d& m, std::vector<Segment>& segments, std::vector<Segment>& arcs)
   {
      if (auto mapped = item.as<IfcSchema::IfcMappedItem>())
      {
         auto source = mapped.MappingSource();
         Eigen::Matrix4d to_global = m * placement_matrix(mapping, mapped.MappingTarget()) * placement_matrix(mapping, source.MappingOrigin());
         for (const auto& mapped_item : source.MappedRepresentation().Items())
            collect_segments(mapping, mapped_item, to_global, segments, arcs);
      }
      else if (auto revolution = item.as<IfcSchema::IfcRevolvedAreaSolid>())
      {
         Float64 radius = 0;
         auto center = circle_center(revolution.SweptArea(), radius);
         if (!center)
            return;
         Eigen::Matrix4d to_global = m * (revolution.Position() ? placement_matrix(mapping, revolution.Position()) : Eigen::Matrix4d::Identity());
         auto axis = revolution.Axis();
         auto location = axis.Location().as<IfcSchema::IfcCartesianPoint>().Coordinates();
         Eigen::Vector3d origin = to_system_point(location);
         Eigen::Vector3d direction = Eigen::Vector3d::UnitZ();
         if (auto axis_direction = axis.Axis())
         {
            auto ratios = axis_direction.DirectionRatios();
            direction = Eigen::Vector3d(ratios[0], ratios[1], 2 < ratios.size() ? ratios[2] : 0.0).normalized();
         }
         Float64 angle = ConvertToSysUnits<IfcSchema::IfcPlaneAngleMeasure>(CIfcImporter::GetUnits(), revolution.Angle());
         int n = std::max(1, (int)ceil(fabs(angle) / WBFL::Units::ConvertToSysUnits(5.0, WBFL::Units::Measure::Degree)));
         Eigen::Vector3d p(center->x(), center->y(), 0.0);
         Eigen::Vector3d previous = transform(to_global, p);
         for (int i = 1; i <= n; i++)
         {
            Eigen::Vector3d next = transform(to_global, origin + Eigen::AngleAxisd(angle * i / n, direction) * (p - origin));
            arcs.push_back({ previous, next, radius });
            previous = next;
         }
      }
      else if (auto extrusion = item.as<IfcSchema::IfcExtrudedAreaSolid>())
      {
         Float64 radius = 0;
         auto center = circle_center(extrusion.SweptArea(), radius);
         if (!center)
            return;
         Eigen::Matrix4d to_global = m * (extrusion.Position() ? placement_matrix(mapping, extrusion.Position()) : Eigen::Matrix4d::Identity());
         auto ratios = extrusion.ExtrudedDirection().DirectionRatios();
         Eigen::Vector3d direction = Eigen::Vector3d(ratios[0], ratios[1], ratios[2]).normalized();
         Eigen::Vector3d a(center->x(), center->y(), 0.0);
         segments.push_back({ transform(to_global, a), transform(to_global, a + direction * to_system_length(extrusion.Depth())), radius });
      }
      else if (auto swept_disk = item.as<IfcSchema::IfcSweptDiskSolid>())
      {
         collect_directrix_segments(swept_disk.Directrix(), m, to_system_length(swept_disk.Radius()), segments);
      }
   }

   // The straight pieces chained end to end into a polyline, starting at the open end farthest from the other pieces. Where
   // two pieces don't meet (a bend between them), the polyline goes through the middle of the gap. Empty if a gap is wider
   // than max_gap
   std::vector<Eigen::Vector3d> chain(std::vector<Segment> segments, Float64 max_gap)
   {
      std::vector<Eigen::Vector3d> points;
      if (segments.empty())
         return points;

      // the start: the segment end whose nearest other segment end is farthest away
      auto nearest_other = [&segments](size_t i, const Eigen::Vector3d& p)
         {
            Float64 d = std::numeric_limits<Float64>::max();
            for (size_t j = 0; j < segments.size(); j++)
               if (j != i)
                  d = std::min({ d, (segments[j].a - p).norm(), (segments[j].b - p).norm() });
            return d;
         };
      size_t first = 0;
      bool bReverse = false;
      Float64 farthest = -1;
      for (size_t i = 0; i < segments.size(); i++)
      {
         for (bool bEnd : {false, true})
         {
            Float64 d = nearest_other(i, bEnd ? segments[i].b : segments[i].a);
            if (farthest < d)
            {
               farthest = d;
               first = i;
               bReverse = bEnd;
            }
         }
      }
      if (bReverse)
         std::swap(segments[first].a, segments[first].b);
      points = { segments[first].a, segments[first].b };
      segments.erase(segments.begin() + first);

      while (!segments.empty())
      {
         const auto& end = points.back();
         size_t next = 0;
         bool bFlip = false;
         Float64 gap = std::numeric_limits<Float64>::max();
         for (size_t i = 0; i < segments.size(); i++)
         {
            for (bool bEnd : {false, true})
            {
               Float64 d = ((bEnd ? segments[i].b : segments[i].a) - end).norm();
               if (d < gap)
               {
                  gap = d;
                  next = i;
                  bFlip = bEnd;
               }
            }
         }
         if (max_gap < gap)
            return {};
         if (bFlip)
            std::swap(segments[next].a, segments[next].b);
         points.back() = (points.back() + segments[next].a) / 2;
         points.push_back(segments[next].b);
         segments.erase(segments.begin() + next);
      }
      return points;
   }

   // A centerline (global coordinates, system units) and the diameter as drawn
   struct Centerline
   {
      std::vector<Eigen::Vector3d> points;
      std::optional<Float64> diameter;
   };

   // the straight pieces and the arcs of the bends of a centerline
   struct Pieces
   {
      std::vector<Segment> segments, arcs;
   };

   // The centerline of the straight pieces chained through the gaps at the bends between them, else (a gap is too wide, e.g. a
   // bar curved to follow the girder) with the arcs of the bends
   Centerline make_centerline(const Pieces& pieces)
   {
      Centerline centerline;
      if (pieces.segments.empty())
         return centerline;
      auto radius = std::max_element(pieces.segments.begin(), pieces.segments.end(), [](const auto& a, const auto& b) {return a.radius < b.radius; })->radius;
      if (0 < radius)
         centerline.diameter = 2 * radius;
      const Float64 max_gap = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Feet);
      centerline.points = chain(pieces.segments, max_gap);
      if (centerline.points.empty() && !pieces.arcs.empty())
      {
         std::vector<Segment> all(pieces.segments);
         all.insert(all.end(), pieces.arcs.begin(), pieces.arcs.end());
         centerline.points = chain(all, max_gap);
      }
      return centerline;
   }

   // The centerlines of a strand or bar from the straight pieces of its body representation (see collect_segments): one for
   // each representation item whose pieces make a line when there's more than one (a group of bars, e.g. the bars of a type
   // written as one element with a mapped item for each bar), else one from all of the pieces. Empty if none
   std::vector<Centerline> get_centerlines(ifcopenshell::geom::abstract_mapping& mapping, IfcSchema::IfcElement element)
   {
      std::vector<Centerline> centerlines;
      if (!element.Representation())
         return centerlines;

      Eigen::Matrix4d to_global = placement_matrix(mapping, element.ObjectPlacement());
      for (const auto& representation : element.Representation().Representations())
      {
         auto identifier = representation.RepresentationIdentifier();
         if (identifier && *identifier != "Body")
            continue;
         std::vector<Pieces> item_pieces;
         Pieces all;
         for (const auto& item : representation.Items())
         {
            item_pieces.emplace_back();
            collect_segments(mapping, item, to_global, item_pieces.back().segments, item_pieces.back().arcs);
            all.segments.insert(all.segments.end(), item_pieces.back().segments.begin(), item_pieces.back().segments.end());
            all.arcs.insert(all.arcs.end(), item_pieces.back().arcs.begin(), item_pieces.back().arcs.end());
         }
         if (all.segments.empty())
            continue;

         // a group: several items that each make a line (side by side bars would chain as if bent into each other)
         std::vector<Centerline> item_centerlines;
         for (const auto& pieces : item_pieces)
         {
            auto centerline = make_centerline(pieces);
            if (2 <= centerline.points.size())
               item_centerlines.push_back(centerline);
         }
         if (1 < item_centerlines.size())
            return item_centerlines;

         auto centerline = make_centerline(all);
         if (2 <= centerline.points.size())
            centerlines.push_back(centerline);
         return centerlines;
      }
      return centerlines;
   }

   // The polyline with points added so they're no more than step apart (the vertices of a centerline are only at the ends of
   // its straight pieces, which may all be outside the girder, e.g. a strand projecting beyond both ends)
   std::vector<Eigen::Vector3d> densify(const std::vector<Eigen::Vector3d>& points, Float64 step)
   {
      std::vector<Eigen::Vector3d> dense;
      for (size_t i = 0; i < points.size(); i++)
      {
         if (0 < i)
         {
            auto n = (int)ceil((points[i] - points[i - 1]).norm() / step);
            for (int j = 1; j < n; j++)
               dense.push_back(points[i - 1] + (points[i] - points[i - 1]) * j / (Float64)n);
         }
         dense.push_back(points[i]);
      }
      return dense;
   }

   // the point halfway along a polyline
   Eigen::Vector3d halfway(const std::vector<Eigen::Vector3d>& points)
   {
      Float64 length = 0;
      for (size_t i = 1; i < points.size(); i++)
         length += (points[i] - points[i - 1]).norm();
      Float64 remaining = length / 2;
      for (size_t i = 1; i < points.size(); i++)
      {
         Float64 l = (points[i] - points[i - 1]).norm();
         if (remaining <= l && 0 < l)
            return points[i - 1] + (points[i] - points[i - 1]) * (remaining / l);
         remaining -= l;
      }
      return points.back();
   }

   // The part of a centerline (densified) inside the girder's plan footprint, from where it enters to where it leaves
   // (found by bisection), and whether it projects beyond the girder at its first and last point (by more than 1 in, e.g.
   // strands bent up into the diaphragms). Empty if no point is inside
   std::vector<Eigen::Vector3d> clip_to_girder(const std::vector<Eigen::Vector3d>& points, const Mesh& mesh, std::array<bool, 2>& projects)
   {
      auto inside = [&mesh](const Eigen::Vector3d& p) {return mesh.z_range_at(p.x(), p.y()).has_value(); };
      auto boundary = [&](Eigen::Vector3d in, Eigen::Vector3d out)
         {
            for (int i = 0; i < 30; i++)
            {
               Eigen::Vector3d mid = (in + out) / 2;
               (inside(mid) ? in : out) = mid;
            }
            return Eigen::Vector3d((in + out) / 2);
         };

      std::vector<bool> is_inside;
      for (const auto& p : points)
         is_inside.push_back(inside(p));
      auto first = std::find(is_inside.begin(), is_inside.end(), true);
      if (first == is_inside.end())
         return {};
      size_t i0 = std::distance(is_inside.begin(), first);
      size_t i1 = is_inside.size() - 1 - std::distance(is_inside.rbegin(), std::find(is_inside.rbegin(), is_inside.rend(), true));

      const Float64 projection_tolerance = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
      std::vector<Eigen::Vector3d> clipped;
      projects = { false, false };
      if (0 < i0)
      {
         auto b = boundary(points[i0], points[i0 - 1]);
         projects[0] = projection_tolerance < (points[i0 - 1] - b).norm() || 1 < i0;
         clipped.push_back(b);
      }
      clipped.insert(clipped.end(), points.begin() + i0, points.begin() + i1 + 1);
      if (i1 + 1 < points.size())
      {
         auto b = boundary(points[i1], points[i1 + 1]);
         projects[1] = projection_tolerance < (points[i1 + 1] - b).norm() || i1 + 2 < points.size();
         clipped.push_back(b);
      }
      return clipped;
   }

   // The girder (IfcBeam id in beam_ids) that the object is a part of, through its assemblies
   std::optional<int> get_assembly_girder(IfcSchema::IfcObjectDefinition object, const std::set<int>& beam_ids, int depth = 0)
   {
      if (!object || 5 < depth)
         return std::nullopt;
      for (const auto& rel : object.Decomposes())
      {
         auto whole = rel.RelatingObject();
         if (beam_ids.count(whole.id()))
            return whole.id();
         if (auto id = get_assembly_girder(whole, beam_ids, depth + 1))
            return id;
      }
      return std::nullopt;
   }

   // the path without points within the position tolerance of the previous point, and without interior points on the line
   // through their neighbors, so a path has its ends and its bends (a path may go across the girder, e.g. a stirrup)
   std::vector<GirderPoint> simplify(const std::vector<GirderPoint>& path)
   {
      const auto tolerance = reinforcement_position_tolerance;
      std::vector<GirderPoint> distinct;
      for (const auto& p : path)
      {
         if (distinct.empty() || tolerance < Eigen::Vector3d(p.x - distinct.back().x, p.z - distinct.back().z, p.y - distinct.back().y).norm())
            distinct.push_back(p);
      }
      if (distinct.size() < 3)
         return distinct;

      auto vector = [](const GirderPoint& p) {return Eigen::Vector3d(p.x, p.z, p.y); };
      std::vector<GirderPoint> simplified{ distinct.front() };
      for (size_t i = 1; i + 1 < distinct.size(); i++)
      {
         Eigen::Vector3d a = vector(simplified.back());
         Eigen::Vector3d ab = vector(distinct[i + 1]) - a;
         Eigen::Vector3d ap = vector(distinct[i]) - a;
         Float64 t = IsZero(ab.squaredNorm()) ? 0.0 : ap.dot(ab) / ab.squaredNorm();
         if (tolerance < (ap - t * ab).norm())
            simplified.push_back(distinct[i]);
      }
      simplified.push_back(distinct.back());
      return simplified;
   }
}

GirderPoint point_at(const std::vector<GirderPoint>& path, Float64 x)
{
   if (path.size() == 1)
      return path.front();
   size_t i = 1;
   while (i + 1 < path.size() && path[i].x < x)
      i++;
   const auto& a = path[i - 1];
   const auto& b = path[i];
   Float64 t = IsEqual(a.x, b.x) ? 1.0 : (x - a.x) / (b.x - a.x);
   return { x, a.z + t * (b.z - a.z), a.y + t * (b.y - a.y) };
}

// The elevations of the top and bottom of a girder along its length: the highest and lowest points of the girder's geometry
// across the girder, sampled at stations along the CL girder and interpolated (extrapolated beyond the first and last
// station). Depths below the top follow a girder that isn't straight, e.g. one drawn along the profile or cambered, so
// reinforcement parallel to its top has a constant depth. The top of an open-top section (U-beam) is at its webs
class GirderReinforcementReader::Profile
{
public:
   Profile(const GirderCenterline& girder)
   {
      if (!girder.mesh)
         return;
      GirderFrame frame(girder);
      Eigen::Vector2d left = frame.left_axis.head<2>().normalized();
      const Float64 half_width = 1.6, step = 0.05; // m, wide enough for the widest girders (e.g. deck bulb tees)
      for (int i = 0; i <= 20; i++)
      {
         Float64 t = 0.02 + 0.96 * i / 20.0;
         Eigen::Vector2d cl = girder.start + t * (girder.end - girder.start);
         std::optional<Float64> top, bottom;
         for (Float64 offset = -half_width; offset <= half_width; offset += step)
         {
            Eigen::Vector2d p = cl + offset * left;
            if (auto range = girder.mesh->z_range_at(p.x(), p.y()))
            {
               top = std::max(top.value_or(range->second), range->second);
               bottom = std::min(bottom.value_or(range->first), range->first);
            }
         }
         if (top)
            m_Points.push_back({ t * frame.length, *top, *bottom });
      }
   }

   std::optional<Float64> Top(Float64 x) const { return At(x, &Point::top); }
   std::optional<Float64> Bottom(Float64 x) const { return At(x, &Point::bottom); }

private:
   struct Point
   {
      Float64 x; // distance from the start of the girder in plan
      Float64 top, bottom; // elevations
   };
   std::vector<Point> m_Points;

   std::optional<Float64> At(Float64 x, Float64 Point::*value) const
   {
      if (m_Points.size() < 2)
         return std::nullopt;
      size_t i = 1;
      while (i + 1 < m_Points.size() && m_Points[i].x < x)
         i++;
      const auto& a = m_Points[i - 1];
      const auto& b = m_Points[i];
      return a.*value + (x - a.x) * (b.*value - a.*value) / (b.x - a.x);
   }
};

ReinforcementCenterlineReader::ReinforcementCenterlineReader(ifcopenshell::file& file)
{
   ifcopenshell::geom::settings settings;
   m_Mapping.reset(ifcopenshell::geom::impl::mapping_implementations().construct(&file, settings));
}

ReinforcementCenterlineReader::~ReinforcementCenterlineReader() = default;

std::vector<ReinforcementCenterline> ReinforcementCenterlineReader::Read(IfcSchema::IfcElement element) const
{
   std::vector<ReinforcementCenterline> centerlines;
   if (!element)
      return centerlines;
   for (const auto& centerline : get_centerlines(*m_Mapping, element))
   {
      if (2 <= centerline.points.size())
         centerlines.push_back({ centerline.points, centerline.diameter });
   }
   return centerlines;
}

GirderReinforcementReader::GirderReinforcementReader(ifcopenshell::file& file, const GirderLayout& layout)
{
   for (const auto& span : layout)
   {
      if (span.size() < 2)
         continue; // the girders of the span weren't imported
      for (const auto& girder : span)
      {
         m_Girders[girder.beam_id] = &girder;
         m_BeamIds.insert(girder.beam_id);
      }
   }

   ifcopenshell::geom::settings settings;
   m_Mapping.reset(ifcopenshell::geom::impl::mapping_implementations().construct(&file, settings));
}

GirderReinforcementReader::~GirderReinforcementReader() = default;

const GirderReinforcementReader::Profile& GirderReinforcementReader::GetProfile(int beam_id) const
{
   auto found = m_Profiles.find(beam_id);
   if (found == m_Profiles.end())
      found = m_Profiles.emplace(beam_id, std::make_unique<Profile>(*m_Girders.at(beam_id))).first;
   return *found->second;
}

const Eigen::AlignedBox3d& GirderReinforcementReader::GetBounds(int beam_id) const
{
   auto found = m_Bounds.find(beam_id);
   if (found == m_Bounds.end())
   {
      Eigen::AlignedBox3d bounds;
      for (const auto& v : m_Girders.at(beam_id)->mesh->verts)
         bounds.extend(v);
      found = m_Bounds.emplace(beam_id, bounds).first;
   }
   return found->second;
}

Float64 GirderReinforcementReader::GetLength(int beam_id) const
{
   return GirderFrame(*m_Girders.at(beam_id)).length;
}

std::optional<Float64> GirderReinforcementReader::GetDepth(int beam_id, Float64 x) const
{
   const auto& profile = GetProfile(beam_id);
   auto top = profile.Top(x);
   auto bottom = profile.Bottom(x);
   if (top && bottom)
      return *top - *bottom;
   return std::nullopt;
}

GirderReinforcementReader::Status GirderReinforcementReader::Read(IfcSchema::IfcElement element, std::vector<GirderReinforcement>& reinforcement)
{
   reinforcement.clear();
   auto centerlines = element ? get_centerlines(*m_Mapping, element) : std::vector<Centerline>();
   if (centerlines.empty())
      return Status::NoGeometry;

   auto assembly_girder = get_assembly_girder(element, m_BeamIds);
   auto status = Status::NoGirder;
   for (const auto& centerline : centerlines)
   {
      GirderReinforcement piece;
      auto piece_status = ReadPiece(centerline.points, assembly_girder, piece);
      if (piece_status == Status::Read)
      {
         piece.diameter = centerline.diameter;
         reinforcement.push_back(piece);
         status = Status::Read;
      }
      else if (status != Status::Read && piece_status == Status::Outside)
      {
         status = Status::Outside;
      }
   }
   return status;
}

GirderReinforcementReader::Status GirderReinforcementReader::ReadPiece(std::vector<Eigen::Vector3d> centerline, std::optional<int> beam_id, GirderReinforcement& reinforcement)
{
   if (!beam_id)
   {
      // the girder whose geometry contains the point halfway along the centerline
      Eigen::Vector3d mid = halfway(centerline);
      for (const auto& [id, girder] : m_Girders)
      {
         auto range = girder->mesh ? girder->mesh->z_range_at(mid.x(), mid.y()) : std::nullopt;
         if (range && range->first <= mid.z() && mid.z() <= range->second)
         {
            beam_id = id;
            break;
         }
      }
   }
   if (!beam_id)
   {
      // the girder that contains the most of it (e.g. a stirrup whose middle is above the girder, in the deck)
      auto points = densify(centerline, 0.1 /*m*/);
      IndexType most = 0;
      for (const auto& [id, girder] : m_Girders)
      {
         if (!girder->mesh)
            continue;
         const auto& bounds = GetBounds(id);
         IndexType n = 0;
         for (const auto& p : points)
         {
            if (!bounds.contains(p))
               continue;
            auto range = girder->mesh->z_range_at(p.x(), p.y());
            if (range && range->first <= p.z() && p.z() <= range->second)
               n++;
         }
         if (most < n)
         {
            most = n;
            beam_id = id;
         }
      }
   }
   if (!beam_id)
      return Status::NoGirder;

   // the part in the girder: strands and bars may project beyond the girder ends (e.g. bent up into the diaphragms)
   const auto& girder = *m_Girders.at(*beam_id);
   std::array<bool, 2> projects{ false, false }; // at the first and last point of the centerline
   if (girder.mesh)
      centerline = clip_to_girder(densify(centerline, 0.25 /*m*/), *girder.mesh, projects);
   if (centerline.size() < 2)
      return Status::Outside;

   // how far each end is from the girder's end face, along the centerline (0 where it projects)
   std::array<std::optional<Float64>, 2> end_distance;
   for (int end = 0; end < 2; end++)
   {
      if (projects[end])
      {
         end_distance[end] = 0.0;
         continue;
      }
      if (!girder.mesh)
         continue;
      const auto& p = (end == 0 ? centerline.front() : centerline.back());
      const auto& q = (end == 0 ? centerline[1] : centerline[centerline.size() - 2]);
      Eigen::Vector3d direction = (p - q).normalized();
      auto inside = [&girder](const Eigen::Vector3d& point) {return girder.mesh->z_range_at(point.x(), point.y()).has_value(); };
      Float64 in = 0, out = end_distance_limit;
      if (inside(p + out * direction))
         continue; // farther than the limit
      for (int i = 0; i < 30; i++)
      {
         Float64 mid = (in + out) / 2;
         (inside(p + mid * direction) ? in : out) = mid;
      }
      end_distance[end] = (in + out) / 2;
   }

   GirderFrame frame(girder);
   const auto& profile = GetProfile(*beam_id);
   reinforcement.beam_id = *beam_id;
   reinforcement.path.clear();
   for (const auto& p : centerline)
   {
      auto point = frame.ToGirder(p);
      if (auto top = profile.Top(point.x))
         point.y = *top - p.z(); // the depth below the girder's top there
      reinforcement.path.push_back(point);
   }
   if (reinforcement.path.back().x < reinforcement.path.front().x)
   {
      std::reverse(reinforcement.path.begin(), reinforcement.path.end());
      std::swap(projects[0], projects[1]);
      std::swap(end_distance[0], end_distance[1]);
   }
   reinforcement.path = simplify(reinforcement.path);
   reinforcement.projects = projects;
   reinforcement.end_distance = end_distance;
   return Status::Read;
}
