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
#include "Strands.h"
#include "IfcImporter.h"
#include "IfcImportUnits.h"
#include "IfcTargetReader.h"

#include <ifcgeom/abstract_mapping.h>

#include <psgLib/StrandData.h>
#include <psgLib/GirderLibraryEntry.h>
#include <IFace/BeamFactory.h>
#include <LRFD/StrandPool.h>

#include <iomanip>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   double inches(double length) { return WBFL::Units::ConvertFromSysUnits(length, WBFL::Units::Measure::Inch); }

   const Float64 position_tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 8.0, WBFL::Units::Measure::Inch);

   // a strand centerline in girder coordinates: x is the distance from the start of the girder in plan, z the offset
   // from the CL girder (positive to the right looking ahead), and y the depth below the top of the girder
   struct GirderPoint
   {
      Float64 x, z, y;
   };

   struct Strand
   {
      IfcSchema::IfcTendon tendon;
      std::vector<GirderPoint> path;
      pgsTypes::StrandType type = pgsTypes::Straight;
      std::array<Float64, 2> debond{ 0, 0 }; // debond length at the start and end (pgsTypes::MemberEndType)
   };

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

   // The 4x4 transformation of a placement, in metres
   Eigen::Matrix4d placement_matrix(ifcopenshell::geom::abstract_mapping& mapping, IfcSchema::IfcObjectPlacement placement)
   {
      if (!placement)
         return Eigen::Matrix4d::Identity();
      auto matrix = ifcopenshell::geom::taxonomy::dcast<ifcopenshell::geom::taxonomy::matrix4>(mapping.map(placement));
      return matrix ? matrix->ccomponents() : Eigen::Matrix4d::Identity();
   }

   // The centerline of a strand (global coordinates, system units), from an IfcSweptDiskSolid along a polyline in its body
   // representation, or empty if it doesn't have one
   std::vector<Eigen::Vector3d> get_centerline(ifcopenshell::geom::abstract_mapping& mapping, IfcSchema::IfcTendon tendon)
   {
      std::vector<Eigen::Vector3d> points;
      if (!tendon.Representation())
         return points;

      const auto& units = CIfcImporter::GetUnits();
      for (const auto& representation : tendon.Representation().Representations())
      {
         for (const auto& item : representation.Items())
         {
            auto solid = item.as<IfcSchema::IfcSweptDiskSolid>();
            auto polyline = solid ? solid.Directrix().as<IfcSchema::IfcPolyline>() : IfcSchema::IfcPolyline();
            if (!polyline)
               continue;

            Eigen::Matrix4d to_global = placement_matrix(mapping, tendon.ObjectPlacement());
            for (const auto& point : polyline.Points())
            {
               auto c = point.Coordinates();
               Eigen::Vector4d local(ConvertToSysUnits<IfcSchema::IfcLengthMeasure>(units, c[0]),
                                     ConvertToSysUnits<IfcSchema::IfcLengthMeasure>(units, 1 < c.size() ? c[1] : 0.0),
                                     ConvertToSysUnits<IfcSchema::IfcLengthMeasure>(units, 2 < c.size() ? c[2] : 0.0), 1.0);
               points.emplace_back((to_global * local).head<3>());
            }
            return points;
         }
      }
      return points;
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

   // the strand at x, by linear interpolation along its path and extrapolation of its first and last segments
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

   Float64 depth_at(const std::vector<GirderPoint>& path, Float64 x) { return point_at(path, x).y; }

   // true if the path bends: an interior point is off the line through the ends
   bool is_harped(const std::vector<GirderPoint>& path)
   {
      if (path.size() < 3)
         return false;
      const auto& a = path.front();
      const auto& b = path.back();
      for (size_t i = 1; i + 1 < path.size(); i++)
      {
         Float64 y = IsEqual(a.x, b.x) ? a.y : a.y + (path[i].x - a.x) * (b.y - a.y) / (b.x - a.x);
         if (position_tolerance < fabs(path[i].y - y))
            return true;
      }
      return false;
   }

   // A position in the girder section: X from the CL girder (positive to the right), Y up from the bottom of the girder
   struct SectionPoint
   {
      Float64 X, Y;
   };

   // A strand's or grid position's section points at the start, left harp point, right harp point, and end
   using SectionPath = std::array<SectionPoint, 4>;
   enum Location { Start = 0, LeftHP = 1, RightHP = 2, End = 3 };

   // The grid positions of a strand type of a library entry (X >= 0, a position with X > 0 is a pair of strands at +/-X).
   // Straight and temporary positions are at their start position at the left harp point and their end position at the right
   std::vector<SectionPath> get_grid(const GirderLibraryEntry* pEntry, pgsTypes::StrandType type, std::vector<bool>* pCanDebond = nullptr)
   {
      std::vector<SectionPath> grid;
      if (type == pgsTypes::Straight)
      {
         for (GridIndexType i = 0; i < pEntry->GetNumStraightStrandCoordinates(); i++)
         {
            Float64 Xs, Ys, Xe, Ye;
            bool canDebond;
            pEntry->GetStraightStrandCoordinates(i, &Xs, &Ys, &Xe, &Ye, &canDebond);
            grid.push_back({ SectionPoint{Xs, Ys}, SectionPoint{Xs, Ys}, SectionPoint{Xe, Ye}, SectionPoint{Xe, Ye} });
            if (pCanDebond)
               pCanDebond->push_back(canDebond);
         }
      }
      else if (type == pgsTypes::Harped)
      {
         for (GridIndexType i = 0; i < pEntry->GetNumHarpedStrandCoordinates(); i++)
         {
            Float64 Xs, Ys, Xhp, Yhp, Xe, Ye;
            pEntry->GetHarpedStrandCoordinates(i, &Xs, &Ys, &Xhp, &Yhp, &Xe, &Ye);
            grid.push_back({ SectionPoint{Xs, Ys}, SectionPoint{Xhp, Yhp}, SectionPoint{Xhp, Yhp}, SectionPoint{Xe, Ye} });
         }
      }
      else
      {
         for (GridIndexType i = 0; i < pEntry->GetNumTemporaryStrandCoordinates(); i++)
         {
            Float64 Xs, Ys, Xe, Ye;
            pEntry->GetTemporaryStrandCoordinates(i, &Xs, &Ys, &Xe, &Ye);
            grid.push_back({ SectionPoint{Xs, Ys}, SectionPoint{Xs, Ys}, SectionPoint{Xe, Ye}, SectionPoint{Xe, Ye} });
         }
      }
      return grid;
   }

   // A filled grid position: its index in the grid of its strand type and its strands (indices of one or two strands)
   struct GridFill
   {
      GridIndexType grid_index;
      std::vector<size_t> strands;
   };

   // a grid position with X > 0 at any location is a pair of strands at +/-X (PGSuper)
   bool is_pair(const SectionPath& p)
   {
      return std::any_of(p.begin(), p.end(), [](const auto& point) {return position_tolerance < point.X; });
   }

   // Fills grid positions with the strands at a location, with the grid shifted up by dy. A pair position holds two strands
   // at +/-X (both at the center where its X is 0, e.g. harped strands bundled at the harp points), or one strand at X = 0
   // (PGSuper puts a single strand of a pair at the center); another position holds one strand. Pairs are filled first,
   // then the other positions, then single strands in pair positions. nullopt unless every strand is at a grid position
   std::optional<std::vector<GridFill>> fill_grid_at(const std::vector<SectionPath>& strands, const std::vector<SectionPath>& grid, Location loc, Float64 dy)
   {
      std::vector<bool> used_strand(strands.size(), false), used_grid(grid.size(), false);
      std::vector<GridFill> fills;
      auto find = [&](Float64 X, Float64 Y, std::optional<size_t> other = std::nullopt) -> std::optional<size_t>
         {
            for (size_t s = 0; s < strands.size(); s++)
            {
               if (!used_strand[s] && s != other && IsEqual(strands[s][loc].X, X, position_tolerance) && IsEqual(strands[s][loc].Y, Y, position_tolerance))
                  return s;
            }
            return std::nullopt;
         };

      for (size_t g = 0; g < grid.size(); g++)
      {
         if (!is_pair(grid[g]))
            continue;
         const auto& p = grid[g][loc];
         auto left = find(-p.X, p.Y + dy);
         auto right = left ? find(p.X, p.Y + dy, left) : std::nullopt;
         if (left && right)
         {
            used_strand[*left] = used_strand[*right] = used_grid[g] = true;
            fills.push_back({ g, {*left, *right} });
         }
      }

      for (bool bPairPositions : {false, true})
      {
         for (size_t g = 0; g < grid.size(); g++)
         {
            const auto& p = grid[g][loc];
            if (used_grid[g] || is_pair(grid[g]) != bPairPositions)
               continue;
            if (auto s = find(0.0, p.Y + dy))
            {
               used_strand[*s] = used_grid[g] = true;
               fills.push_back({ g, {*s} });
            }
         }
      }

      if (std::find(used_strand.begin(), used_strand.end(), false) != used_strand.end())
         return std::nullopt;

      std::sort(fills.begin(), fills.end(), [](const auto& a, const auto& b) {return a.grid_index < b.grid_index; });
      return fills;
   }

   // The vertical shift of the grid at a location that puts the filled positions on their strands, or nullopt if there's none
   std::optional<Float64> grid_shift_at(const std::vector<GridFill>& fills, const std::vector<SectionPath>& strands, const std::vector<SectionPath>& grid, Location loc)
   {
      std::optional<Float64> dy;
      for (const auto& fill : fills)
      {
         const auto& p = grid[fill.grid_index][loc];
         Float64 X = (fill.strands.size() == 1 ? 0.0 : p.X);
         if (fill.strands.size() == 2 && position_tolerance < X && 0 < strands[fill.strands[0]][loc].X * strands[fill.strands[1]][loc].X)
            return std::nullopt; // a pair is on both sides of the CL girder (unless both are at the center)
         for (auto s : fill.strands)
         {
            const auto& q = strands[s][loc];
            if (!IsEqual(fabs(q.X), X, position_tolerance))
               return std::nullopt;
            if (!dy)
               dy = q.Y - p.Y;
            else if (!IsEqual(q.Y - p.Y, *dy, position_tolerance))
               return std::nullopt;
         }
      }
      return dy.value_or(0.0);
   }

   // The grid fill of a strand type: the filled positions and the grid's vertical shift at each location. Straight and
   // temporary strands aren't shifted. A harped grid is shifted at the harp points and at the ends: the shift at the left
   // harp point is found by putting the lowest strand there on each grid position, and the other shifts follow from the fill.
   // nullopt (with the reason) unless every strand is at a grid position
   struct TypeFill
   {
      std::vector<GridFill> fills;
      std::array<Float64, 4> shift{ 0, 0, 0, 0 };
   };
   std::optional<TypeFill> fill_type(const std::vector<SectionPath>& strands, const std::vector<SectionPath>& grid, pgsTypes::StrandType type, std::string& reason)
   {
      TypeFill result;
      if (strands.empty())
         return result;

      if (type != pgsTypes::Harped)
      {
         auto fills = fill_grid_at(strands, grid, Start, 0.0);
         auto shift = fills ? grid_shift_at(*fills, strands, grid, End) : std::nullopt;
         if (!fills || !shift || position_tolerance < fabs(*shift))
         {
            reason = std::string(type == pgsTypes::Straight ? "straight" : "temporary") + " strands aren't at the library entry's grid positions";
            return std::nullopt;
         }
         result.fills = *fills;
         return result;
      }

      // The fill is found at the start, where the strands of a pair are usually apart, else at the left harp point. The
      // shift there is found by putting the lowest strand on each grid position in turn
      for (auto fill_loc : { Start, LeftHP })
      {
         const auto& q = (*std::min_element(strands.begin(), strands.end(), [fill_loc](const auto& a, const auto& b) {return a[fill_loc].Y < b[fill_loc].Y; }))[fill_loc];
         for (const auto& p : grid)
         {
            if (!IsEqual(fabs(q.X), p[fill_loc].X, position_tolerance) && !IsZero(q.X, position_tolerance))
               continue;
            auto fills = fill_grid_at(strands, grid, fill_loc, q.Y - p[fill_loc].Y);
            if (!fills)
               continue;

            bool bFits = true;
            for (auto loc : { Start, LeftHP, RightHP, End })
            {
               auto shift = grid_shift_at(*fills, strands, grid, loc);
               if (!shift)
               {
                  bFits = false;
                  break;
               }
               result.shift[loc] = *shift;
            }
            if (bFits)
            {
               result.fills = *fills;
               return result;
            }
         }
      }
      reason = "harped strands aren't at the library entry's grid positions, with a common vertical offset at each end and harp point";
      return std::nullopt;
   }

   CDirectStrandFillCollection make_fill_collection(const std::vector<GridFill>& fills)
   {
      CDirectStrandFillCollection collection;
      for (const auto& fill : fills)
         collection.AddFill(CDirectStrandFillInfo(fill.grid_index, fill.strands.size()));
      return collection;
   }

   // Fills the library entry's strand grid with a girder's strands (pgsTypes::sdtDirectSelection). False, with the reason,
   // unless every strand is at a grid position (the harped strands with a common vertical offset at each end and harp
   // point that the entry allows), both strands of a pair are debonded alike, and debonded strands are at debondable
   // positions. The harped strand offsets are input as PGSuper measures them, from the strands: the top harped strand to
   // the top of the girder at the ends, and the bottom harped strand to the bottom of the girder at the harp points
   bool fill_strand_grid(const GirderLibraryEntry* pEntry, const std::array<std::vector<const Strand*>, 3>& by_type, const std::array<Float64, 4>& X, CStrandData& strand_data, std::string& reason)
   {
      auto factory = pEntry->GetBeamFactory();
      Float64 Hstart = factory->GetBeamHeight(pEntry->GetDimensions(), pgsTypes::metStart);
      Float64 Hend = factory->GetBeamHeight(pEntry->GetDimensions(), pgsTypes::metEnd);
      std::array<Float64, 4> Hg{ Hstart, Hstart, Hend, Hend };

      std::array<TypeFill, 3> type_fills;
      std::vector<bool> can_debond;
      for (int i = 0; i < 3; i++)
      {
         auto type = (pgsTypes::StrandType)i;
         auto grid = get_grid(pEntry, type, type == pgsTypes::Straight ? &can_debond : nullptr);
         std::vector<SectionPath> paths;
         for (const auto* strand : by_type[i])
         {
            SectionPath path;
            for (int loc = 0; loc < 4; loc++)
            {
               auto p = point_at(strand->path, X[loc]);
               path[loc] = { p.z, Hg[loc] - p.y };
            }
            paths.push_back(path);
         }
         auto fill = fill_type(paths, grid, type, reason);
         if (!fill)
            return false;
         type_fills[i] = *fill;
      }

      const auto& harped = by_type[pgsTypes::Harped];
      if (!harped.empty())
      {
         const auto& shift = type_fills[pgsTypes::Harped].shift;
         if (pEntry->GetAdjustableStrandType() == pgsTypes::asStraight)
         {
            reason = "the library entry's adjustable strands are straight, but the model's are harped";
            return false;
         }
         if (!pEntry->IsVerticalAdjustmentAllowedEnd() && (position_tolerance < fabs(shift[Start]) || position_tolerance < fabs(shift[End])))
         {
            reason = "the harped strands are offset at the girder ends, which the library entry doesn't allow";
            return false;
         }
         if (!pEntry->IsVerticalAdjustmentAllowedHP() && (position_tolerance < fabs(shift[LeftHP]) || position_tolerance < fabs(shift[RightHP])))
         {
            reason = "the harped strands are offset at the harp points, which the library entry doesn't allow";
            return false;
         }
      }

      std::array<std::vector<CDebondData>, 3> debonds;
      for (int i = 0; i < 3; i++)
      {
         for (const auto& fill : type_fills[i].fills)
         {
            const auto& debond = by_type[i][fill.strands.front()]->debond;
            for (auto s : fill.strands)
            {
               const auto& other = by_type[i][s]->debond;
               if (!IsEqual(debond[0], other[0], position_tolerance) || !IsEqual(debond[1], other[1], position_tolerance))
               {
                  reason = "the two strands of a grid position are debonded differently";
                  return false;
               }
            }
            if (IsZero(debond[0]) && IsZero(debond[1]))
               continue;
            if (i != pgsTypes::Straight || !can_debond[fill.grid_index])
            {
               reason = "a debonded strand is at a grid position that can't be debonded";
               return false;
            }
            CDebondData debond_data;
            debond_data.strandTypeGridIdx = fill.grid_index;
            debond_data.Length = debond;
            debonds[i].push_back(debond_data);
         }
      }

      strand_data.SetStrandDefinitionType(pgsTypes::sdtDirectSelection);
      strand_data.SetDirectStrandFillStraight(make_fill_collection(type_fills[pgsTypes::Straight].fills));
      strand_data.SetDirectStrandFillHarped(make_fill_collection(type_fills[pgsTypes::Harped].fills));
      strand_data.SetDirectStrandFillTemporary(make_fill_collection(type_fills[pgsTypes::Temporary].fills));
      for (int i = 0; i < 3; i++)
         strand_data.SetDebonding((pgsTypes::StrandType)i, debonds[i]);

      if (!harped.empty())
      {
         strand_data.SetAdjustableStrandType(pgsTypes::asHarped);
         strand_data.SetHarpStrandOffsetMeasurementAtEnd(hsoTOP2TOP);
         strand_data.SetHarpStrandOffsetMeasurementAtHarpPoint(hsoBOTTOM2BOTTOM);
         for (auto [end, loc] : { std::make_pair(pgsTypes::metStart, Start), std::make_pair(pgsTypes::metEnd, End) })
         {
            Float64 top = std::numeric_limits<Float64>::max();
            for (const auto* strand : harped)
               top = std::min(top, depth_at(strand->path, X[loc]));
            strand_data.SetHarpStrandOffsetAtEnd(end, top);
         }
         for (auto [end, loc] : { std::make_pair(pgsTypes::metStart, LeftHP), std::make_pair(pgsTypes::metEnd, RightHP) })
         {
            Float64 bottom = std::numeric_limits<Float64>::max();
            for (const auto* strand : harped)
               bottom = std::min(bottom, Hg[loc] - depth_at(strand->path, X[loc]));
            strand_data.SetHarpStrandOffsetAtHarpPoint(end, bottom);
         }
      }
      return true;
   }

   // Inputs a girder's strands as individual strands (pgsTypes::sdtDirectStrandInput) with common harp points
   void fill_individual_strands(const std::vector<Strand>& strands, const std::array<Float64, 4>& X, CStrandData& strand_data, const std::string& girder_name)
   {
      strand_data.SetStrandDefinitionType(pgsTypes::sdtDirectStrandInput);
      strand_data.SetHarpPoints(X[Start], X[LeftHP], X[RightHP], X[End]);

      CStrandRowCollection rows;
      IndexType nMoving = 0; // harped strands whose offset from the CL girder changes along the girder
      for (const auto& strand : strands)
      {
         CStrandRow row;
         row.m_StrandType = strand.type;
         row.m_nStrands = 1; // CStrandData counts the strands of the rows, also for individual strands
         for (int i = 0; i < 4; i++)
         {
            row.m_Y[i] = depth_at(strand.path, X[i]);
            row.m_Face[i] = pgsTypes::TopFace;
         }

         if (strand.type == pgsTypes::Harped)
         {
            // individual strand input has one offset for a strand: the offset at the harp points
            row.m_Z = (strand.path[1].z + strand.path[strand.path.size() - 2].z) / 2;
            if (position_tolerance < std::max(fabs(strand.path.front().z - row.m_Z), fabs(strand.path.back().z - row.m_Z)))
               nMoving++;
         }
         else
         {
            row.m_Z = (strand.path.front().z + strand.path.back().z) / 2;
         }

         for (auto end : { pgsTypes::metStart, pgsTypes::metEnd })
         {
            row.m_bIsDebonded[end] = IsGT(0.0, strand.debond[end]);
            row.m_DebondLength[end] = strand.debond[end];
         }
         rows.push_back(row);
      }
      strand_data.SetStrandRows(rows);

      if (0 < nMoving)
      {
         std::ostringstream os;
         os << "Strands of " << girder_name << ": " << nMoving << " harped strands change their offset from the CL girder along the girder. Individual strand input has one offset for each strand, so the offset at the harp points is used.";
         WBFL::System::Logger::Warning(os.str());
      }
   }

   // The strand material of the strands, from their diameter (IfcTendonType.NominalDiameter, else the swept disk) and the
   // tendon.fpu, tendon.grade, and tendon.corrosion_treatment targets. Low relaxation unless tendon.grade says stress
   // relieved. nullptr if no strand in the pool has the diameter
   const WBFL::Materials::PsStrand* get_strand_material(const std::vector<const Strand*>& strands, const std::string& girder_name)
   {
      using namespace WBFL::Materials;
      const auto& reader = CIfcImporter::GetTargetReader();
      const auto& units = CIfcImporter::GetUnits();
      auto tendon = strands.front()->tendon;

      std::optional<Float64> diameter;
      for (const auto& rel : tendon.IsTypedBy())
      {
         auto tendon_type = rel.RelatingType().as<IfcSchema::IfcTendonType>();
         if (tendon_type && tendon_type.NominalDiameter())
            diameter = ConvertToSysUnits<IfcSchema::IfcPositiveLengthMeasure>(units, *tendon_type.NominalDiameter());
      }
      if (!diameter && tendon.Representation())
      {
         for (const auto& representation : tendon.Representation().Representations())
            for (const auto& item : representation.Items())
               if (auto solid = item.as<IfcSchema::IfcSweptDiskSolid>(); solid && !diameter)
                  diameter = 2 * ConvertToSysUnits<IfcSchema::IfcPositiveLengthMeasure>(units, solid.Radius());
      }
      if (!diameter)
         return nullptr;

      auto grade = PsStrand::Grade::Gr1860;
      std::string grade_source = "the default";
      std::optional<std::string> grade_text;
      if (auto reading = reader.Read("tendon.grade", tendon))
         grade_text = std::get<std::string>(reading->value);
      if (auto reading = reader.Read("tendon.fpu", tendon))
      {
         Float64 fpu = WBFL::Units::ConvertFromSysUnits(std::get<Float64>(reading->value), WBFL::Units::Measure::KSI);
         grade = (fpu < 260 ? PsStrand::Grade::Gr1725 : fpu < 285 ? PsStrand::Grade::Gr1860 : PsStrand::Grade::Gr2070);
         grade_source = reading->Source();
      }
      else if (grade_text)
      {
         if (grade_text->find("250") != std::string::npos) grade = PsStrand::Grade::Gr1725;
         else if (grade_text->find("300") != std::string::npos) grade = PsStrand::Grade::Gr2070;
         grade_source = "tendon.grade \"" + *grade_text + "\"";
      }

      auto type = PsStrand::Type::LowRelaxation;
      if (grade_text)
      {
         std::string upper(*grade_text);
         std::transform(upper.begin(), upper.end(), upper.begin(), ::toupper);
         if (upper.find("STRESS RELIEVED") != std::string::npos)
            type = PsStrand::Type::StressRelieved;
      }

      auto coating = PsStrand::Coating::None;
      if (auto reading = reader.Read("tendon.corrosion_treatment", tendon))
      {
         if (std::get<std::string>(reading->value) == "EPOXYCOATED")
            coating = PsStrand::Coating::SmoothEpoxy; // IFC doesn't say smooth or grit impregnated
      }

      const PsStrand* best = nullptr;
      WBFL::LRFD::StrandIter iter(grade, type, coating);
      for (iter.Begin(); iter; iter.Next())
      {
         const auto* candidate = iter.GetCurrentStrand();
         if (!best || fabs(candidate->GetNominalDiameter() - *diameter) < fabs(best->GetNominalDiameter() - *diameter))
            best = candidate;
      }
      if (!best || WBFL::Units::ConvertToSysUnits(0.5, WBFL::Units::Measure::Millimeter) < fabs(best->GetNominalDiameter() - *diameter))
      {
         std::ostringstream os;
         os << "Strands of " << girder_name << ": no strand in the strand pool has a diameter of " << std::fixed << std::setprecision(3) << inches(*diameter) << " in.";
         WBFL::System::Logger::Warning(os.str());
         return nullptr;
      }

      USES_CONVERSION;
      std::ostringstream os;
      os << "Strands of " << girder_name << ": " << T2A(best->GetName().c_str()) << " (diameter " << std::fixed << std::setprecision(3) << inches(*diameter) << " in, grade from " << grade_source << ").";
      WBFL::System::Logger::Info(os.str());
      return best;
   }

   // The jacking force of the strands: the sum of IfcTendon.PreStress times the strand area, else of IfcTendon.TensionForce
   std::optional<Float64> get_jacking_force(const std::vector<const Strand*>& strands, const WBFL::Materials::PsStrand* pStrand)
   {
      const auto& units = CIfcImporter::GetUnits();
      Float64 force = 0;
      for (const auto* strand : strands)
      {
         if (auto prestress = strand->tendon.PreStress(); prestress && pStrand)
            force += ConvertToSysUnits<IfcSchema::IfcPressureMeasure>(units, *prestress) * pStrand->GetNominalArea();
         else if (auto tension = strand->tendon.TensionForce())
            force += ConvertToSysUnits<IfcSchema::IfcForceMeasure>(units, *tension);
         else
            return std::nullopt;
      }
      return force;
   }
}

void import_strands(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const GirderLayout& layout, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;
   const auto& reader = CIfcImporter::GetTargetReader();
   if (!reader.HasSelector(ElementKind::Tendon))
   {
      WBFL::System::Logger::Info("The mapping table has no tendon role, so strands aren't imported.");
      return;
   }

   auto tendons = reader.Select(ElementKind::Tendon, file);
   if (tendons.empty())
   {
      WBFL::System::Logger::Info("The model has no strands (" + reader.GetTable().GetSelector(ElementKind::Tendon)->Describe() + "). The girders keep the template's strands.");
      return;
   }

   std::map<int, const GirderCenterline*> girders; // by IfcBeam id
   std::set<int> beam_ids;
   for (const auto& span : layout)
   {
      if (span.size() < 2)
         continue; // the girders of the span weren't imported
      for (const auto& girder : span)
      {
         girders[girder.beam_id] = &girder;
         beam_ids.insert(girder.beam_id);
      }
   }

   ifcopenshell::geom::settings settings;
   std::unique_ptr<ifcopenshell::geom::abstract_mapping> mapping(ifcopenshell::geom::impl::mapping_implementations().construct(&file, settings));

   // the strands of each girder, in girder coordinates
   std::map<int, std::vector<Strand>> girder_strands;
   IndexType nWithoutGeometry = 0, nWithoutGirder = 0;
   for (const auto& object : tendons)
   {
      auto tendon = object.as<IfcSchema::IfcTendon>();
      auto centerline = tendon ? get_centerline(*mapping, tendon) : std::vector<Eigen::Vector3d>();
      if (centerline.size() < 2)
      {
         nWithoutGeometry++;
         continue;
      }

      auto beam_id = get_assembly_girder(tendon, beam_ids);
      if (!beam_id)
      {
         // the girder whose geometry contains the strand's mid-point
         Eigen::Vector3d mid = (centerline.front() + centerline.back()) / 2;
         for (const auto& [id, girder] : girders)
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
         nWithoutGirder++;
         continue;
      }

      GirderFrame frame(*girders[*beam_id]);
      Strand strand;
      strand.tendon = tendon;
      for (const auto& p : centerline)
         strand.path.push_back(frame.ToGirder(p));
      if (strand.path.back().x < strand.path.front().x)
         std::reverse(strand.path.begin(), strand.path.end());

      auto status = reader.Read("tendon.status", tendon);
      if (status && std::get<std::string>(status->value) == "TEMPORARY")
         strand.type = pgsTypes::Temporary;
      else if (is_harped(strand.path))
         strand.type = pgsTypes::Harped;

      for (auto [end, target] : { std::make_pair(pgsTypes::metStart, "tendon.debond_start"), std::make_pair(pgsTypes::metEnd, "tendon.debond_end") })
      {
         if (auto reading = reader.Read(target, tendon); reading && IsGT(0.0, std::get<Float64>(reading->value)))
            strand.debond[end] = std::get<Float64>(reading->value);
      }

      girder_strands[*beam_id].push_back(std::move(strand));
   }

   if (0 < nWithoutGeometry)
      WBFL::System::Logger::Warning(std::to_string(nWithoutGeometry) + " strands were skipped: their geometry isn't an IfcSweptDiskSolid along a polyline.");
   if (0 < nWithoutGirder)
      WBFL::System::Logger::Warning(std::to_string(nWithoutGirder) + " strands were skipped: they aren't part of an imported girder and aren't inside one.");

   for (const auto& [beam_id, strands] : girder_strands)
   {
      const auto& girder = *girders[beam_id];
      auto beam = file.instance_by_id(beam_id).as<IfcSchema::IfcBeam>();
      std::string girder_name = beam.Name().value_or("unnamed beam");
      Float64 Ls = GirderFrame(girder).length;

      std::array<std::vector<const Strand*>, 3> by_type;
      for (const auto& strand : strands)
         by_type[strand.type].push_back(&strand);

      // the harp points are common to the harped strands: the first and last bend of each strand, averaged
      Float64 X1 = 0, X2 = Ls;
      if (!by_type[pgsTypes::Harped].empty())
      {
         Float64 sum1 = 0, sum2 = 0, spread = 0;
         for (const auto* strand : by_type[pgsTypes::Harped])
         {
            sum1 += strand->path[1].x;
            sum2 += strand->path[strand->path.size() - 2].x;
         }
         auto n = (Float64)by_type[pgsTypes::Harped].size();
         X1 = sum1 / n;
         X2 = sum2 / n;
         for (const auto* strand : by_type[pgsTypes::Harped])
            spread = std::max({ spread, fabs(strand->path[1].x - X1), fabs(strand->path[strand->path.size() - 2].x - X2) });
         if (position_tolerance < spread)
         {
            std::ostringstream os;
            os << "Strands of " << girder_name << ": the harp points of the harped strands differ by up to " << std::fixed << std::setprecision(3) << inches(spread)
               << " in. PGSuper's harped strands share their harp points, so the average is used.";
            WBFL::System::Logger::Warning(os.str());
         }
      }

      // the template's strand data, for what the model doesn't say (e.g. the temporary strand usage)
      auto* pGirder = bridge_desc.GetGirderGroup(girder.girder_key.groupIndex)->GetGirder(girder.girder_key.girderIndex);
      auto* pSegment = pGirder->GetSegment(0);
      CStrandData strand_data(pSegment->Strands);

      // the library entry's strand grid when every strand is on it (G2, decision 2), else individual strands
      std::array<Float64, 4> X{ 0.0, X1, X2, Ls };
      const GirderLibraryEntry* pEntry = bridge_desc.UseSameGirderForEntireBridge() ? bridge_desc.GetGirderLibraryEntry() : pGirder->GetGirderLibraryEntry();
      std::string reason = "the girder has no library entry";
      bool bGrid = pEntry && fill_strand_grid(pEntry, by_type, X, strand_data, reason);
      if (!bGrid)
      {
         WBFL::System::Logger::Info("Strands of " + girder_name + ": " + reason + ", so they are input as individual strands.");
         fill_individual_strands(strands, X, strand_data, girder_name);
      }

      for (int i = 0; i < 3; i++)
      {
         auto strand_type = (pgsTypes::StrandType)i;
         const auto* pStrand = pSegment->Strands.GetStrandMaterial(strand_type); // the template's, unless the model says otherwise
         if (!by_type[i].empty())
         {
            if (auto* pModelStrand = get_strand_material(by_type[i], girder_name))
               pStrand = pModelStrand;
         }
         strand_data.SetStrandMaterial(strand_type, pStrand);

         if (auto force = by_type[i].empty() ? std::optional<Float64>(0.0) : get_jacking_force(by_type[i], pStrand))
         {
            strand_data.IsPjackCalculated(strand_type, false);
            strand_data.SetPjack(strand_type, *force);
         }
         else
         {
            strand_data.IsPjackCalculated(strand_type, true);
            WBFL::System::Logger::Warning("Strands of " + girder_name + ": the jacking force isn't in the model (IfcTendon.PreStress or TensionForce), so PGSuper computes it.");
         }
      }
      pSegment->Strands = strand_data;

      std::ostringstream os;
      os << "Strands of " << girder_name << ": " << by_type[pgsTypes::Straight].size() << " straight, " << by_type[pgsTypes::Harped].size() << " harped, "
         << by_type[pgsTypes::Temporary].size() << " temporary, ";
      if (bGrid)
         os << "filled into the strand grid of " << T2A(pEntry->GetName().c_str());
      else
         os << "input as individual strands";
      if (!by_type[pgsTypes::Harped].empty())
         os << ", harp points (from the model) at " << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(X1, WBFL::Units::Measure::Feet) << " and " << WBFL::Units::ConvertFromSysUnits(X2, WBFL::Units::Measure::Feet) << " ft";
      os << ".";
      WBFL::System::Logger::Info(os.str());
   }
}
