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
#include "IfcImporter.h"

#include <psgLib/LongitudinalRebarData.h>
#include <LRFD/RebarPool.h>

#include <iomanip>

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

void import_longitudinal_rebar(const GirderBars& girder_bars, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   // The girders' bars come from the template, which has the library entries' bars. Those are seed values, not the
   // girder's reinforcement, so a girder gets only the bars in the model.
   for_each_segment(bridge_desc, [](CPrecastSegmentData* pSegment) {pSegment->LongitudinalRebarData.RebarRows.clear(); });
   if (!girder_bars.HasBars())
   {
      WBFL::System::Logger::Info("Longitudinal girder bars aren't imported: " + girder_bars.GetReason() + ". The girders have no longitudinal reinforcement.");
      return;
   }

   const auto& geometry = girder_bars.GetGeometry();
   IndexType nWithoutBars = 0, nWithoutSize = 0;
   for (const auto& [beam_id, girder] : geometry.GetGirders())
   {
      const auto& all_bars = girder_bars.GetBars(beam_id);
      if (all_bars.empty())
      {
         nWithoutBars++;
         continue;
      }

      // the longitudinal bars (stirrups and other bars aren't)
      std::vector<Bar> longitudinal;
      std::map<BarMaterial, IndexType> materials; // number of bars of each material
      for (const auto& girder_bar : all_bars)
      {
         auto bar = longitudinal_piece(girder_bar.centerline);
         if (!bar)
            continue;
         if (girder_bar.size == Rebar::Size::bsNone)
         {
            nWithoutSize++;
            continue;
         }
         bar->element = girder_bar.element;
         bar->size = girder_bar.size;
         if (girder_bar.material)
            materials[*girder_bar.material]++;
         longitudinal.push_back(*bar);
      }

      std::string girder_name = girder_bars.GetGirderName(beam_id);
      auto* pSegment = GirderBars::GetSegment(bridge_desc, *girder);
      CLongitudinalRebarData rebar_data(pSegment->LongitudinalRebarData); // the material if the model doesn't say
      rebar_data.RebarRows.clear();

      bool bMixed;
      if (auto material = most_common_material(materials, bMixed))
      {
         rebar_data.BarType = material->first;
         rebar_data.BarGrade = material->second;
      }
      if (bMixed)
         WBFL::System::Logger::Warning("Longitudinal bars of " + girder_name + ": the bars have different materials. PGSuper has one for the girder, so the most common is used.");

      IndexType nLaps = 0, nUnpaired = 0, nTop = 0, nBottom = 0, nBars = 0;
      std::vector<Bar> bars = merge_laps(longitudinal, nLaps);
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

   if (0 < nWithoutSize)
      WBFL::System::Logger::Warning(std::to_string(nWithoutSize) + " longitudinal girder bars were skipped: their size isn't known (rebar.bar_size, the nominal diameter, or a drawn diameter of a bar size).");
   if (0 < nWithoutBars)
      WBFL::System::Logger::Info(std::to_string(nWithoutBars) + " girders have no reinforcing bars in the model, so they have no longitudinal reinforcement.");
}
