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
#include "DeckRebar.h"
#include "GirderBars.h"
#include "GirderReinforcement.h"
#include "IfcImporter.h"
#include "IfcTargetReader.h"

#include <psgLib/BridgeDescription2.h>
#include <psgLib/DeckDescription2.h>
#include <IFace/Alignment.h>
#include <LRFD/RebarPool.h>
#include <psgLib/GirderLabel.h>

#include <iomanip>
#include <numeric>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   using WBFL::Materials::Rebar;

   const Float64 inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
   const Float64 inside_fraction = 0.9;           // a deck bar's centerline is this much inside the deck slab
   const Float64 inside_tolerance = 0.25 * inch;  // a centerline point this close to the deck's surface is inside it
   const Float64 line_offset_tolerance = 2 * inch; // bars of a line (e.g. lapped end to end) are this close laterally
   const Float64 line_depth_tolerance = 1 * inch;  // and in depth
   const Float64 primary_coverage = 0.75;          // a primary line covers this much of the bridge length
   const Float64 mat_gap = 1 * inch;               // primary lines farther apart in depth than this are in different mats

   // A longitudinal deck bar: its size and material, where it starts and ends (stations), its offset from the alignment
   // (the mean of its ends), and the depth of its middle below the top of the deck and its height above the bottom
   struct DeckBar
   {
      Rebar::Size size;
      std::optional<BarMaterial> material;
      Float64 start, end, offset, depth, height;
      IndexType line = INVALID_INDEX;
      CDeckRebarData::RebarMat mat = CDeckRebarData::TopMat;
   };

   // A line of bars: bars at the same offset and depth
   struct Line
   {
      std::vector<IndexType> bars;
      Float64 offset, depth;
      bool bPrimary = false;
   };

   std::vector<Eigen::Vector3d> densify(const std::vector<Eigen::Vector3d>& points, Float64 step)
   {
      std::vector<Eigen::Vector3d> result{ points.front() };
      for (size_t i = 1; i < points.size(); i++)
      {
         auto n = std::max<size_t>(1, (size_t)std::ceil((points[i] - points[i - 1]).norm() / step));
         for (size_t j = 1; j <= n; j++)
            result.push_back(points[i - 1] + (points[i] - points[i - 1]) * (Float64)j / (Float64)n);
      }
      return result;
   }

   // the point halfway along a path
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

   Float64 median(std::vector<Float64> values)
   {
      std::sort(values.begin(), values.end());
      return values.empty() ? 0.0 : values[values.size() / 2];
   }

   // the most common size of the bars of a line
   Rebar::Size line_size(const Line& line, const std::vector<DeckBar>& bars)
   {
      std::map<Rebar::Size, IndexType> sizes;
      for (auto i : line.bars)
         sizes[bars[i].size]++;
      return std::max_element(sizes.begin(), sizes.end(), [](const auto& a, const auto& b) {return a.second < b.second; })->first;
   }

   // The area per width of lines: one bar of each over the width they cover (from the first to the last plus one spacing,
   // the median spacing between them, else spacing_if_one). Also returns the spacing
   std::pair<Float64, Float64> area_per_width(const std::vector<const Line*>& lines, const std::vector<DeckBar>& bars, Float64 spacing_if_one)
   {
      std::vector<Float64> offsets;
      Float64 area = 0;
      for (const auto* line : lines)
      {
         offsets.push_back(line->offset);
         area += nominal_area(line_size(*line, bars));
      }
      std::sort(offsets.begin(), offsets.end());
      std::vector<Float64> spacings;
      for (size_t i = 1; i < offsets.size(); i++)
         spacings.push_back(offsets[i] - offsets[i - 1]);
      Float64 spacing = spacings.empty() ? spacing_if_one : median(spacings);
      Float64 width = offsets.back() - offsets.front() + spacing;
      return { IsZero(width) ? 0.0 : area / width, spacing };
   }

   std::string in2_per_ft(Float64 value)
   {
      std::ostringstream os;
      os << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Inch2PerFoot) << " in^2/ft";
      return os.str();
   }

   std::string inches(Float64 value)
   {
      std::ostringstream os;
      os << std::fixed << std::setprecision(2) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Inch) << " in";
      return os.str();
   }

   std::string feet(Float64 value)
   {
      std::ostringstream os;
      os << std::fixed << std::setprecision(2) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Feet) << " ft";
      return os.str();
   }

   std::string size_names(const std::vector<const Line*>& lines, const std::vector<DeckBar>& bars)
   {
      USES_CONVERSION;
      std::set<Rebar::Size> sizes;
      for (const auto* line : lines)
         sizes.insert(line_size(*line, bars));
      std::string names;
      for (auto size : sizes)
         names += (names.empty() ? "" : ", ") + std::string(T2A(WBFL::LRFD::RebarPool::GetBarSize(size).c_str()));
      return names;
   }
}

void import_deck_rebar(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const Mesh& deck, const std::vector<Float64>& pier_stations, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;
   auto& rebar_data = bridge_desc.GetDeckDescription()->DeckRebarData;
   const auto& reader = CIfcImporter::GetTargetReader();
   if (!reader.HasSelector(ElementKind::Rebar) || pier_stations.size() < 2)
   {
      WBFL::System::Logger::Info("The mapping table has no rebar role, so the deck keeps the existing deck reinforcement.");
      return;
   }

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> point;
   point.CoCreateInstance(CLSID_Point2d);
   auto station_and_offset = [&](const Eigen::Vector3d& p)
      {
         point->Move(p.x(), p.y());
         Float64 station, offset;
         pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &station, &offset);
         return std::make_pair(station, offset);
      };
   auto in_deck = [&deck](const Eigen::Vector3d& p)
      {
         auto range = deck.z_range_at(p.x(), p.y());
         return range && range->first - inside_tolerance <= p.z() && p.z() <= range->second + inside_tolerance;
      };

   // the longitudinal bars in the deck
   ReinforcementCenterlineReader centerline_reader(file);
   std::vector<DeckBar> bars;
   IndexType nTransverse = 0, nWithoutSize = 0;
   for (const auto& object : reader.Select(ElementKind::Rebar, file))
   {
      auto element = object.as<IfcSchema::IfcReinforcingBar>();
      if (!element)
         continue;
      for (const auto& centerline : centerline_reader.Read(element))
      {
         auto points = densify(centerline.points, 0.25 /*m*/);
         auto nInside = std::count_if(points.begin(), points.end(), in_deck);
         if (nInside < inside_fraction * points.size())
            continue; // not a deck bar (e.g. a girder's or a barrier's)

         auto [s0, o0] = station_and_offset(centerline.points.front());
         auto [s1, o1] = station_and_offset(centerline.points.back());
         if (fabs(s1 - s0) <= fabs(o1 - o0))
         {
            nTransverse++;
            continue;
         }

         auto size = get_bar_size(element, centerline.diameter);
         if (size == Rebar::Size::bsNone)
         {
            nWithoutSize++;
            continue;
         }

         auto middle = halfway(centerline.points);
         auto range = deck.z_range_at(middle.x(), middle.y());
         if (!range)
            continue;
         bars.push_back({ size, get_bar_material(element), std::min(s0, s1), std::max(s0, s1), (o0 + o1) / 2, range->second - middle.z(), middle.z() - range->first });
      }
   }

   if (bars.empty())
   {
      std::ostringstream os;
      os << "The model has no longitudinal bars in the deck";
      if (0 < nTransverse)
         os << " (" << nTransverse << " transverse bars, which aren't PGSuper data)";
      os << ", so the deck keeps the existing deck reinforcement.";
      WBFL::System::Logger::Info(os.str());
      return;
   }

   // lines of bars at the same offset and depth
   std::vector<IndexType> order(bars.size());
   std::iota(order.begin(), order.end(), 0);
   std::sort(order.begin(), order.end(), [&bars](auto a, auto b) {return bars[a].offset < bars[b].offset; });
   std::vector<Line> lines;
   for (auto i : order)
   {
      auto& bar = bars[i];
      auto found = std::find_if(lines.rbegin(), lines.rend(), [&bar](const Line& line)
         {
            return fabs(line.offset - bar.offset) < line_offset_tolerance && fabs(line.depth - bar.depth) < line_depth_tolerance;
         });
      if (found == lines.rend())
      {
         lines.push_back({ {}, bar.offset, bar.depth });
         found = lines.rbegin();
      }
      found->bars.push_back(i);
      bar.line = std::distance(lines.begin(), found.base()) - 1;
   }

   // the primary lines cover most of the bridge
   Float64 first = pier_stations.front(), last = pier_stations.back();
   for (auto& line : lines)
   {
      std::vector<std::pair<Float64, Float64>> intervals;
      for (auto i : line.bars)
         intervals.emplace_back(std::max(bars[i].start, first), std::min(bars[i].end, last));
      std::sort(intervals.begin(), intervals.end());
      Float64 covered = 0, end = first;
      for (auto [a, b] : intervals)
      {
         a = std::max(a, end);
         if (a < b)
         {
            covered += b - a;
            end = b;
         }
      }
      line.bPrimary = primary_coverage * (last - first) <= covered;
   }

   std::vector<const Line*> primary;
   for (const auto& line : lines)
   {
      if (line.bPrimary)
         primary.push_back(&line);
   }
   if (primary.empty())
   {
      WBFL::System::Logger::Warning(std::to_string(bars.size()) + " longitudinal bars are in the deck, but none of them make lines that run the length of the bridge, so the deck keeps the existing deck reinforcement.");
      return;
   }

   // the mats, from the depths of the primary lines
   std::vector<Float64> depths;
   for (const auto* line : primary)
      depths.push_back(line->depth);
   std::sort(depths.begin(), depths.end());
   std::optional<Float64> split;
   Float64 gap = 0;
   for (size_t i = 1; i < depths.size(); i++)
   {
      if (gap < depths[i] - depths[i - 1])
      {
         gap = depths[i] - depths[i - 1];
         split = (depths[i] + depths[i - 1]) / 2;
      }
   }
   std::array<std::optional<Float64>, 2> mat_depth; // the mean depth of the primary lines of each mat
   if (split && mat_gap < gap)
   {
      std::array<std::vector<Float64>, 2> d;
      for (auto x : depths)
         d[x < *split ? CDeckRebarData::TopMat : CDeckRebarData::BottomMat].push_back(x);
      for (auto mat : { CDeckRebarData::TopMat, CDeckRebarData::BottomMat })
         mat_depth[mat] = std::accumulate(d[mat].begin(), d[mat].end(), 0.0) / d[mat].size();
   }
   else
   {
      // one mat: the top mat if it's in the upper half of the deck
      Float64 mean = std::accumulate(depths.begin(), depths.end(), 0.0) / depths.size();
      const auto& bar = bars[primary.front()->bars.front()];
      mat_depth[mean < (bar.depth + bar.height) / 2 ? CDeckRebarData::TopMat : CDeckRebarData::BottomMat] = mean;
   }
   for (auto& bar : bars)
   {
      const auto& line = lines[bar.line];
      if (!mat_depth[CDeckRebarData::BottomMat])
         bar.mat = CDeckRebarData::TopMat;
      else if (!mat_depth[CDeckRebarData::TopMat])
         bar.mat = CDeckRebarData::BottomMat;
      else
         bar.mat = fabs(line.depth - *mat_depth[CDeckRebarData::TopMat]) <= fabs(line.depth - *mat_depth[CDeckRebarData::BottomMat]) ? CDeckRebarData::TopMat : CDeckRebarData::BottomMat;
   }

   // the primary reinforcement of each mat
   std::array<Float64, 2> primary_spacing{ 0, 0 };
   std::ostringstream os;
   os << "Deck reinforcement (as lump sums, the area of bars per width of deck): " << bars.size() << " longitudinal bars in the deck in " << lines.size() << " lines";
   if (0 < nTransverse)
      os << " (and " << nTransverse << " transverse bars, which aren't PGSuper data)";
   os << ".";
   for (auto mat : { CDeckRebarData::TopMat, CDeckRebarData::BottomMat })
   {
      std::vector<const Line*> mat_lines;
      for (const auto* line : primary)
      {
         if (bars[line->bars.front()].mat == mat)
            mat_lines.push_back(line);
      }
      std::string mat_name(mat == CDeckRebarData::TopMat ? "top mat" : "bottom mat");
      Float64& lump_sum = (mat == CDeckRebarData::TopMat ? rebar_data.TopLumpSum : rebar_data.BottomLumpSum);
      if (mat == CDeckRebarData::TopMat)
         rebar_data.TopRebarSize = Rebar::Size::bsNone;
      else
         rebar_data.BottomRebarSize = Rebar::Size::bsNone;
      if (mat_lines.empty())
      {
         lump_sum = 0;
         os << " The " << mat_name << " has no bars that run the length of the bridge.";
         continue;
      }

      auto [area, spacing] = area_per_width(mat_lines, bars, 0.0);
      primary_spacing[mat] = spacing;
      lump_sum = area;

      // the cover: the depth (top mat) or height (bottom mat) of the bars' centroid
      Float64 sum_a = 0, sum_ay = 0;
      std::map<BarMaterial, IndexType> materials;
      for (const auto* line : mat_lines)
      {
         for (auto i : line->bars)
         {
            Float64 a = nominal_area(bars[i].size);
            sum_a += a;
            sum_ay += a * (mat == CDeckRebarData::TopMat ? bars[i].depth : bars[i].height);
            if (bars[i].material)
               materials[*bars[i].material]++;
         }
      }
      Float64 cover = sum_ay / sum_a;
      (mat == CDeckRebarData::TopMat ? rebar_data.TopCover : rebar_data.BottomCover) = cover;

      bool bMixed;
      if (auto material = most_common_material(materials, bMixed))
      {
         (mat == CDeckRebarData::TopMat ? rebar_data.TopRebarType : rebar_data.BottomRebarType) = material->first;
         (mat == CDeckRebarData::TopMat ? rebar_data.TopRebarGrade : rebar_data.BottomRebarGrade) = material->second;
      }

      os << " The " << mat_name << ": " << mat_lines.size() << " lines of " << size_names(mat_lines, bars) << " at " << inches(spacing) << ", " << in2_per_ft(area)
         << ", " << (mat == CDeckRebarData::TopMat ? "centroid " + inches(cover) + " below the top" : "centroid " + inches(cover) + " above the bottom") << " of the deck.";
   }

   // the supplemental reinforcement at the piers: the bars of the other lines, at the pier they cross nearest their middle
   std::map<std::pair<PierIndexType, CDeckRebarData::RebarMat>, std::vector<IndexType>> at_pier;
   IndexType nUnplaced = 0;
   for (IndexType i = 0; i < bars.size(); i++)
   {
      const auto& bar = bars[i];
      if (lines[bar.line].bPrimary)
         continue;
      Float64 middle = (bar.start + bar.end) / 2;
      std::optional<PierIndexType> pier;
      for (PierIndexType pierIdx = 0; pierIdx < pier_stations.size(); pierIdx++)
      {
         if (bar.start < pier_stations[pierIdx] && pier_stations[pierIdx] < bar.end && (!pier || fabs(pier_stations[pierIdx] - middle) < fabs(pier_stations[*pier] - middle)))
            pier = pierIdx;
      }
      if (pier)
         at_pier[{ *pier, bar.mat }].push_back(i);
      else
         nUnplaced++;
   }

   rebar_data.NegMomentRebar.clear();
   for (const auto& [key, indices] : at_pier)
   {
      auto [pierIdx, mat] = key;
      std::set<IndexType> line_indices;
      std::vector<Float64> left, right;
      std::map<BarMaterial, IndexType> materials;
      for (auto i : indices)
      {
         line_indices.insert(bars[i].line);
         left.push_back(pier_stations[pierIdx] - bars[i].start);
         right.push_back(bars[i].end - pier_stations[pierIdx]);
         if (bars[i].material)
            materials[*bars[i].material]++;
      }
      std::vector<const Line*> pier_lines;
      for (auto l : line_indices)
         pier_lines.push_back(&lines[l]);
      auto [area, spacing] = area_per_width(pier_lines, bars, primary_spacing[mat]);

      CDeckRebarData::NegMomentRebarData nm;
      nm.PierIdx = pierIdx;
      nm.Mat = mat;
      nm.LumpSum = area;
      nm.RebarSize = Rebar::Size::bsNone;
      nm.Spacing = 0;
      nm.LeftCutoff = median(left);
      nm.RightCutoff = median(right);
      nm.RebarType = (mat == CDeckRebarData::TopMat ? rebar_data.TopRebarType : rebar_data.BottomRebarType);
      nm.RebarGrade = (mat == CDeckRebarData::TopMat ? rebar_data.TopRebarGrade : rebar_data.BottomRebarGrade);
      bool bMixed;
      if (auto material = most_common_material(materials, bMixed))
      {
         nm.RebarType = material->first;
         nm.RebarGrade = material->second;
      }
      rebar_data.NegMomentRebar.push_back(nm);

      os << " At " << T2A(LABEL_PIER_EX(pierIdx == 0 || pierIdx == pier_stations.size() - 1, pierIdx)) << ", " << (mat == CDeckRebarData::TopMat ? "top" : "bottom") << " mat: "
         << indices.size() << " bars in " << pier_lines.size() << " lines of " << size_names(pier_lines, bars) << " at " << inches(spacing) << ", " << in2_per_ft(area)
         << ", from " << feet(nm.LeftCutoff) << " back to " << feet(nm.RightCutoff) << " ahead of the CL pier.";
   }
   WBFL::System::Logger::Info(os.str());

   if (0 < nUnplaced)
      WBFL::System::Logger::Warning(std::to_string(nUnplaced) + " longitudinal deck bars don't run the length of the bridge and don't cross a pier. PGSuper's supplemental deck reinforcement is at the piers, so they aren't imported.");
   if (0 < nWithoutSize)
      WBFL::System::Logger::Warning(std::to_string(nWithoutSize) + " longitudinal deck bars were skipped: their size isn't known (rebar.bar_size, the nominal diameter, or a drawn diameter of a bar size).");
}
