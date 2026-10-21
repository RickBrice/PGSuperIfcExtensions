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
#include "LongitudinalRebar.h"
#include "GirderReinforcement.h"
#include "IfcImporter.h"
#include "IfcImportUnits.h"
#include "IfcTargetReader.h"

#include <psgLib/LongitudinalRebarData.h>
#include <LRFD/RebarPool.h>

#include <iomanip>
#include <regex>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   using WBFL::Materials::Rebar;

   const Float64 position_tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 8.0, WBFL::Units::Measure::Inch);
   const Float64 center_tolerance = WBFL::Units::ConvertToSysUnits(0.5, WBFL::Units::Measure::Inch); // bars this far off center are centered (the CL girder from the geometry isn't exact)
   const Float64 end_tolerance = WBFL::Units::ConvertToSysUnits(3.0, WBFL::Units::Measure::Inch); // a bar end this close to the girder end is at the end
   const Float64 min_length = WBFL::Units::ConvertToSysUnits(2.0, WBFL::Units::Measure::Feet);
   const Float64 max_slope = tan(WBFL::Units::ConvertToSysUnits(2.0, WBFL::Units::Measure::Degree));


   const std::array<Rebar::Size, 11> bar_sizes{ Rebar::Size::bs3, Rebar::Size::bs4, Rebar::Size::bs5, Rebar::Size::bs6, Rebar::Size::bs7, Rebar::Size::bs8,
      Rebar::Size::bs9, Rebar::Size::bs10, Rebar::Size::bs11, Rebar::Size::bs14, Rebar::Size::bs18 };

   // A longitudinal bar in its girder
   struct Bar
   {
      IfcSchema::IfcReinforcingBar element;
      Rebar::Size size = Rebar::Size::bsNone;
      Float64 x1, x2; // the ends, from the start of the girder in plan
      Float64 z, y;   // offset from the CL girder and depth below the top of the girder (of its longest straight piece, at mid-length)
      std::array<std::optional<Float64>, 2> end_distance; // to the girder's end faces
      std::array<bool, 2> projects{ false, false };
   };

   // The longest straight piece of a bar's path if it runs along the girder and is at least min_length long. Its end
   // distances and projections are the bar's where the piece is at the bar's ends
   std::optional<Bar> longitudinal_piece(const GirderReinforcement& reinforcement)
   {
      const auto& path = reinforcement.path;
      if (path.size() < 2)
         return std::nullopt;
      size_t longest = 1;
      auto length = [&path](size_t i) {return Eigen::Vector3d(path[i].x - path[i - 1].x, path[i].z - path[i - 1].z, path[i].y - path[i - 1].y).norm(); };
      for (size_t i = 2; i < path.size(); i++)
      {
         if (length(longest) < length(i))
            longest = i;
      }
      const auto& a = path[longest - 1];
      const auto& b = path[longest];
      Float64 dx = b.x - a.x;
      if (dx < min_length || max_slope * dx < sqrt(pow(b.z - a.z, 2) + pow(b.y - a.y, 2)))
         return std::nullopt;

      Bar bar;
      bar.x1 = a.x;
      bar.x2 = b.x;
      bar.z = (a.z + b.z) / 2;
      bar.y = (a.y + b.y) / 2;
      if (longest == 1)
      {
         bar.end_distance[0] = reinforcement.end_distance[0];
         bar.projects[0] = reinforcement.projects[0];
      }
      if (longest == path.size() - 1)
      {
         bar.end_distance[1] = reinforcement.end_distance[1];
         bar.projects[1] = reinforcement.projects[1];
      }
      return bar;
   }

   // the nominal diameter of a bar size (sizes are the same for every type and grade), 0 if the pool doesn't have it
   Float64 nominal_diameter(Rebar::Size size)
   {
      const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(Rebar::Type::A615, Rebar::Grade60, size);
      return pRebar ? pRebar->GetNominalDimension() : 0.0;
   }

   // The bar size: rebar.bar_size (a number, e.g. "#5", "5", "No. 5"), else the nominal diameter of the bar type or the
   // bar, else the diameter as drawn. bsNone if no size has the diameter (within 0.1 in)
   Rebar::Size get_bar_size(IfcSchema::IfcReinforcingBar bar, std::optional<Float64> drawn_diameter, std::string& source)
   {
      USES_CONVERSION;
      const auto& reader = CIfcImporter::GetTargetReader();
      const auto& units = CIfcImporter::GetUnits();
      if (auto reading = reader.Read("rebar.bar_size", bar))
      {
         std::smatch match;
         const auto& text = std::get<std::string>(reading->value);
         if (std::regex_search(text, match, std::regex("[0-9]+")))
         {
            auto size = WBFL::LRFD::RebarPool::GetBarSize(A2T(("#" + match.str()).c_str()));
            if (size != Rebar::Size::bsNone)
            {
               source = "rebar.bar_size";
               return size;
            }
         }
      }

      std::optional<Float64> diameter;
      for (const auto& rel : bar.IsTypedBy())
      {
         auto bar_type = rel.RelatingType().as<IfcSchema::IfcReinforcingBarType>();
         if (bar_type && bar_type.NominalDiameter())
            diameter = ConvertToSysUnits<IfcSchema::IfcPositiveLengthMeasure>(units, *bar_type.NominalDiameter());
      }
      source = "the bar type's nominal diameter";
      if (!diameter && bar.NominalDiameter())
      {
         diameter = ConvertToSysUnits<IfcSchema::IfcPositiveLengthMeasure>(units, *bar.NominalDiameter());
         source = "the bar's nominal diameter";
      }
      if (!diameter)
      {
         diameter = drawn_diameter;
         source = "the geometry";
      }
      if (!diameter)
         return Rebar::Size::bsNone;

      auto best = std::min_element(bar_sizes.begin(), bar_sizes.end(), [&diameter](auto a, auto b) {return fabs(nominal_diameter(a) - *diameter) < fabs(nominal_diameter(b) - *diameter); });
      return fabs(nominal_diameter(*best) - *diameter) <= WBFL::Units::ConvertToSysUnits(0.1, WBFL::Units::Measure::Inch) ? *best : Rebar::Size::bsNone;
   }

   // The bar material: the type from rebar.grade (A615 unless the text says A706 or A1035), the grade from rebar.fy, else a
   // grade number in rebar.grade. nullopt if the model says neither
   std::optional<std::pair<Rebar::Type, Rebar::Grade>> get_bar_material(IfcSchema::IfcReinforcingBar bar)
   {
      const auto& reader = CIfcImporter::GetTargetReader();
      std::optional<std::string> text;
      if (auto reading = reader.Read("rebar.grade", bar))
         text = std::get<std::string>(reading->value);
      std::optional<Float64> fy;
      if (auto reading = reader.Read("rebar.fy", bar))
         fy = WBFL::Units::ConvertFromSysUnits(std::get<Float64>(reading->value), WBFL::Units::Measure::KSI);
      if (!text && !fy)
         return std::nullopt;

      auto type = Rebar::Type::A615;
      if (text && text->find("A706") != std::string::npos)
         type = Rebar::Type::A706;
      else if (text && text->find("A1035") != std::string::npos)
         type = Rebar::Type::A1035;

      const std::array<std::pair<Float64, Rebar::Grade>, 6> grades{ { {40, Rebar::Grade40}, {60, Rebar::Grade60}, {75, Rebar::Grade75}, {80, Rebar::Grade80}, {100, Rebar::Grade100}, {120, Rebar::Grade120} } };
      if (!fy && text)
      {
         // the grade number in the text (e.g. "AASHTO M31 (A615) - Grade 60", or just "60")
         std::smatch match;
         std::string rest = std::regex_replace(*text, std::regex("A(615|706|1035)"), "");
         if (std::regex_search(rest, match, std::regex("[0-9]+")))
            fy = std::stod(match.str());
      }
      auto grade = Rebar::Grade60;
      if (fy)
         grade = std::min_element(grades.begin(), grades.end(), [&fy](const auto& a, const auto& b) {return fabs(a.first - *fy) < fabs(b.first - *fy); })->second;
      return std::make_pair(type, grade);
   }

   // Bars of the same size that lap (collinear within 1.5 bar diameters, overlapping along the girder) as one bar
   std::vector<Bar> merge_laps(std::vector<Bar> bars, IndexType& nLaps)
   {
      std::sort(bars.begin(), bars.end(), [](const auto& a, const auto& b) {return a.x1 < b.x1; });
      std::vector<Bar> merged;
      for (const auto& bar : bars)
      {
         Float64 tolerance = 1.5 * nominal_diameter(bar.size);
         auto lap = std::find_if(merged.begin(), merged.end(), [&](const auto& m)
            {
               return m.size == bar.size && fabs(m.z - bar.z) <= tolerance && fabs(m.y - bar.y) <= tolerance && bar.x1 <= m.x2 + position_tolerance;
            });
         if (lap == merged.end())
         {
            merged.push_back(bar);
         }
         else
         {
            nLaps++;
            if (lap->x2 < bar.x2)
            {
               lap->x2 = bar.x2;
               lap->end_distance[1] = bar.end_distance[1];
               lap->projects[1] = bar.projects[1];
            }
         }
      }
      return merged;
   }

   // A bar's row data except the number of bars and their spacing
   CLongitudinalRebarData::RebarRow row_of(const Bar& bar, Float64 Lg, std::optional<Float64> depth)
   {
      CLongitudinalRebarData::RebarRow row;
      row.BarSize = bar.size;
      Float64 db = nominal_diameter(bar.size);
      if (depth && *depth / 2 < bar.y)
      {
         row.Face = pgsTypes::BottomFace;
         row.Cover = *depth - bar.y - db / 2;
      }
      else
      {
         row.Face = pgsTypes::TopFace;
         row.Cover = bar.y - db / 2;
      }

      bool bAtStart = bar.end_distance[0] && *bar.end_distance[0] <= end_tolerance;
      bool bAtEnd = bar.end_distance[1] && *bar.end_distance[1] <= end_tolerance;
      if (bAtStart && bAtEnd)
      {
         row.BarLayout = pgsTypes::blFullLength;
      }
      else if (bAtStart)
      {
         row.BarLayout = pgsTypes::blFromLeft;
         row.BarLength = bar.x2;
      }
      else if (bAtEnd)
      {
         row.BarLayout = pgsTypes::blFromRight;
         row.BarLength = Lg - bar.x1;
      }
      else if (IsEqual(bar.x1, Lg - bar.x2, position_tolerance))
      {
         row.BarLayout = pgsTypes::blMidGirderLength;
         row.BarLength = bar.x2 - bar.x1;
      }
      else
      {
         bool bLeft = bar.x1 <= Lg - bar.x2;
         row.BarLayout = bLeft ? pgsTypes::blFromLeft : pgsTypes::blFromRight;
         row.DistFromEnd = bLeft ? bar.x1 : Lg - bar.x2;
         row.BarLength = bar.x2 - bar.x1;
      }

      // PGSuper extends bars that are full length or start at the end
      bool bFromLeftEnd = row.BarLayout == pgsTypes::blFullLength || (row.BarLayout == pgsTypes::blFromLeft && IsZero(row.DistFromEnd));
      bool bFromRightEnd = row.BarLayout == pgsTypes::blFullLength || (row.BarLayout == pgsTypes::blFromRight && IsZero(row.DistFromEnd));
      row.bExtendedLeft = bar.projects[0] && bFromLeftEnd;
      row.bExtendedRight = bar.projects[1] && bFromRightEnd;
      return row;
   }

   bool same_row(const CLongitudinalRebarData::RebarRow& a, const CLongitudinalRebarData::RebarRow& b)
   {
      return a.Face == b.Face && a.BarSize == b.BarSize && a.BarLayout == b.BarLayout && a.bExtendedLeft == b.bExtendedLeft && a.bExtendedRight == b.bExtendedRight &&
         IsEqual(a.Cover, b.Cover, position_tolerance) && IsEqual(a.DistFromEnd, b.DistFromEnd, position_tolerance) && IsEqual(a.BarLength, b.BarLength, position_tolerance);
   }

   // The rows of bars with the same row data at offsets zs from the CL girder: one row if they're equally spaced and centered
   // on the CL girder (within center_tolerance), else a row for each pair at +/-z and one for a bar at the center. nUnpaired
   // counts the bars that can't be input
   void add_rows(CLongitudinalRebarData::RebarRow row, std::vector<Float64> zs, std::vector<CLongitudinalRebarData::RebarRow>& rows, IndexType& nUnpaired)
   {
      std::sort(zs.begin(), zs.end());
      IndexType n = zs.size();
      Float64 spacing = (1 < n ? (zs.back() - zs.front()) / (n - 1) : 0.0);
      bool bOneRow = IsZero(zs.front() + zs.back(), 2 * center_tolerance);
      for (IndexType i = 0; bOneRow && i < n; i++)
         bOneRow = IsEqual(zs[i], zs.front() + i * spacing, position_tolerance);
      if (bOneRow && (n == 1 || position_tolerance < spacing))
      {
         row.NumberOfBars = n;
         row.BarSpacing = spacing;
         rows.push_back(row);
         return;
      }

      std::vector<bool> used(n, false);
      for (IndexType i = 0; i < n; i++)
      {
         if (used[i])
            continue;
         used[i] = true;
         if (IsZero(zs[i], center_tolerance))
         {
            row.NumberOfBars = 1;
            row.BarSpacing = 0;
            rows.push_back(row);
            continue;
         }
         IndexType mirror = i + 1;
         while (mirror < n && (used[mirror] || !IsEqual(zs[mirror], -zs[i], 2 * center_tolerance)))
            mirror++;
         if (mirror == n)
         {
            nUnpaired++;
            continue;
         }
         used[mirror] = true;
         row.NumberOfBars = 2;
         row.BarSpacing = fabs(zs[i] - zs[mirror]);
         rows.push_back(row);
      }
   }
}

void import_longitudinal_rebar(ifcopenshell::file& file, const GirderLayout& layout, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   // The girders' bars come from the template, which has the library entries' bars. Those are seed values, not the
   // girder's reinforcement, so a girder gets only the bars in the model.
   for (GroupIndexType grpIdx = 0; grpIdx < bridge_desc.GetGirderGroupCount(); grpIdx++)
   {
      auto* pGroup = bridge_desc.GetGirderGroup(grpIdx);
      for (GirderIndexType gdrIdx = 0; gdrIdx < pGroup->GetGirderCount(); gdrIdx++)
      {
         auto* pGirder = pGroup->GetGirder(gdrIdx);
         for (SegmentIndexType segIdx = 0; segIdx < pGirder->GetSegmentCount(); segIdx++)
            pGirder->GetSegment(segIdx)->LongitudinalRebarData.RebarRows.clear();
      }
   }

   const auto& reader = CIfcImporter::GetTargetReader();
   if (!reader.HasSelector(ElementKind::Rebar))
   {
      WBFL::System::Logger::Info("The mapping table has no rebar role, so longitudinal girder bars aren't imported. The girders have no longitudinal reinforcement.");
      return;
   }

   auto elements = reader.Select(ElementKind::Rebar, file);
   if (elements.empty())
   {
      WBFL::System::Logger::Info("The model has no reinforcing bars (" + reader.GetTable().GetSelector(ElementKind::Rebar)->Describe() + "). The girders have no longitudinal reinforcement.");
      return;
   }

   GirderReinforcementReader geometry(file, layout);

   // the longitudinal bars of each girder, and the girders that have bars
   std::map<int, std::vector<Bar>> girder_bars;
   std::set<int> girders_with_bars;
   std::map<int, std::map<std::pair<Rebar::Type, Rebar::Grade>, IndexType>> girder_materials; // number of bars of each material
   IndexType nWithoutGeometry = 0, nWithoutSize = 0;
   for (const auto& object : elements)
   {
      auto element = object.as<IfcSchema::IfcReinforcingBar>();
      std::vector<GirderReinforcement> centerlines; // more than one for a group of bars
      auto result = element ? geometry.Read(element, centerlines) : GirderReinforcementReader::Status::NoGeometry;
      if (result == GirderReinforcementReader::Status::NoGeometry)
         nWithoutGeometry++;
      if (result != GirderReinforcementReader::Status::Read)
         continue; // bars of other elements (deck, piers) aren't in a girder

      std::optional<std::pair<Rebar::Type, Rebar::Grade>> material;
      bool bMaterial = false;
      for (const auto& reinforcement : centerlines)
      {
         girders_with_bars.insert(reinforcement.beam_id);
         auto bar = longitudinal_piece(reinforcement);
         if (!bar)
            continue; // stirrups and other bars

         std::string size_source;
         bar->element = element;
         bar->size = get_bar_size(element, reinforcement.diameter, size_source);
         if (bar->size == Rebar::Size::bsNone)
         {
            nWithoutSize++;
            continue;
         }
         if (!bMaterial)
         {
            material = get_bar_material(element);
            bMaterial = true;
         }
         if (material)
            girder_materials[reinforcement.beam_id][*material]++;
         girder_bars[reinforcement.beam_id].push_back(*bar);
      }
   }

   if (0 < nWithoutGeometry)
      WBFL::System::Logger::Info(std::to_string(nWithoutGeometry) + " reinforcing bars were skipped: their geometry has no straight pieces that make a centerline (circular profiles extruded along a line, or swept disks along a polyline or an indexed poly curve).");
   if (0 < nWithoutSize)
      WBFL::System::Logger::Warning(std::to_string(nWithoutSize) + " longitudinal girder bars were skipped: their size isn't known (rebar.bar_size, the nominal diameter, or a drawn diameter of a bar size).");

   IndexType nWithoutBars = 0;
   for (const auto& [beam_id, girder] : geometry.GetGirders())
   {
      if (!girders_with_bars.count(beam_id))
      {
         nWithoutBars++;
         continue;
      }

      auto beam = file.instance_by_id(beam_id).as<IfcSchema::IfcBeam>();
      std::string girder_name = beam.Name().value_or("unnamed beam");
      auto* pSegment = bridge_desc.GetGirderGroup(girder->girder_key.groupIndex)->GetGirder(girder->girder_key.girderIndex)->GetSegment(0);
      CLongitudinalRebarData rebar_data(pSegment->LongitudinalRebarData); // the material if the model doesn't say
      rebar_data.RebarRows.clear();

      const auto& materials = girder_materials[beam_id];
      if (!materials.empty())
      {
         auto most = std::max_element(materials.begin(), materials.end(), [](const auto& a, const auto& b) {return a.second < b.second; });
         if (WBFL::LRFD::RebarPool::GetInstance()->GetRebar(most->first.first, most->first.second, Rebar::Size::bs5))
         {
            rebar_data.BarType = most->first.first;
            rebar_data.BarGrade = most->first.second;
         }
         if (1 < materials.size())
            WBFL::System::Logger::Warning("Longitudinal bars of " + girder_name + ": the bars have different materials. PGSuper has one for the girder, so the most common is used.");
      }

      IndexType nLaps = 0, nUnpaired = 0, nTop = 0, nBottom = 0, nBars = 0;
      std::vector<Bar> bars = merge_laps(girder_bars[beam_id], nLaps);
      Float64 Lg = geometry.GetLength(beam_id);

      // group the bars by their row data
      std::vector<std::pair<CLongitudinalRebarData::RebarRow, std::vector<Float64>>> groups;
      for (const auto& bar : bars)
      {
         auto row = row_of(bar, Lg, geometry.GetDepth(beam_id, (bar.x1 + bar.x2) / 2));
         auto group = std::find_if(groups.begin(), groups.end(), [&row](const auto& g) {return same_row(g.first, row); });
         if (group == groups.end())
            groups.emplace_back(row, std::vector<Float64>{ bar.z });
         else
            group->second.push_back(bar.z);
      }
      for (const auto& [row, zs] : groups)
         add_rows(row, zs, rebar_data.RebarRows, nUnpaired);
      for (const auto& row : rebar_data.RebarRows)
      {
         nBars += row.NumberOfBars;
         (row.Face == pgsTypes::TopFace ? nTop : nBottom) += row.NumberOfBars;
         if (row.Cover < 0)
            WBFL::System::Logger::Warning("Longitudinal bars of " + girder_name + ": a bar is outside the girder (negative cover).");
      }
      pSegment->LongitudinalRebarData = rebar_data;

      std::ostringstream os;
      os << "Longitudinal bars of " << girder_name << ": " << nBars << " bars (" << nTop << " referenced to the top of the girder, " << nBottom << " to the bottom) in "
         << rebar_data.RebarRows.size() << " rows, " << T2A(WBFL::LRFD::RebarPool::GetMaterialName(rebar_data.BarType, rebar_data.BarGrade).c_str()) << ".";
      if (0 < nLaps)
         os << " " << nLaps << " lapped bars are joined to the bars they lap (PGSuper doesn't model laps).";
      WBFL::System::Logger::Info(os.str());
      if (0 < nUnpaired)
         WBFL::System::Logger::Warning("Longitudinal bars of " + girder_name + ": " + std::to_string(nUnpaired) + " bars were skipped. PGSuper centers the bars of a row on the CL girder, and these bars have no bar opposite them.");
   }

   if (0 < nWithoutBars)
      WBFL::System::Logger::Info(std::to_string(nWithoutBars) + " girders have no reinforcing bars in the model, so they have no longitudinal reinforcement.");
}
