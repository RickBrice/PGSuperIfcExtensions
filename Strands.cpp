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

   // depth of the strand at x, by linear interpolation along its path and extrapolation of its first and last segments
   Float64 depth_at(const std::vector<GirderPoint>& path, Float64 x)
   {
      if (path.size() == 1)
         return path.front().y;
      size_t i = 1;
      while (i + 1 < path.size() && path[i].x < x)
         i++;
      const auto& a = path[i - 1];
      const auto& b = path[i];
      return IsEqual(a.x, b.x) ? b.y : a.y + (x - a.x) * (b.y - a.y) / (b.x - a.x);
   }

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
      strand_data.SetStrandDefinitionType(pgsTypes::sdtDirectStrandInput);
      strand_data.SetHarpPoints(0.0, X1, X2, Ls);

      CStrandRowCollection rows;
      IndexType nMoving = 0; // harped strands whose offset from the CL girder changes along the girder
      for (const auto& strand : strands)
      {
         CStrandRow row;
         row.m_StrandType = strand.type;
         row.m_nStrands = 1; // CStrandData counts the strands of the rows, also for individual strands
         std::array<Float64, 4> X{ 0.0, X1, X2, Ls };
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

         for (auto [end, target] : { std::make_pair(pgsTypes::metStart, "tendon.debond_start"), std::make_pair(pgsTypes::metEnd, "tendon.debond_end") })
         {
            if (auto reading = reader.Read(target, strand.tendon); reading && IsGT(0.0, std::get<Float64>(reading->value)))
            {
               row.m_bIsDebonded[end] = true;
               row.m_DebondLength[end] = std::get<Float64>(reading->value);
            }
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
         << by_type[pgsTypes::Temporary].size() << " temporary, input as individual strands";
      if (!by_type[pgsTypes::Harped].empty())
         os << ", harp points at " << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(X1, WBFL::Units::Measure::Feet) << " and " << WBFL::Units::ConvertFromSysUnits(X2, WBFL::Units::Measure::Feet) << " ft";
      os << ".";
      WBFL::System::Logger::Info(os.str());
   }
}
