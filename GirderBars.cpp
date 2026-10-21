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
#include "GirderBars.h"
#include "IfcImporter.h"
#include "IfcImportUnits.h"
#include "IfcTargetReader.h"

#include <psgLib/PrecastSegmentData.h>
#include <LRFD/RebarPool.h>

#include <regex>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   using WBFL::Materials::Rebar;

   const std::array<Rebar::Size, 11> bar_sizes{ Rebar::Size::bs3, Rebar::Size::bs4, Rebar::Size::bs5, Rebar::Size::bs6, Rebar::Size::bs7, Rebar::Size::bs8,
      Rebar::Size::bs9, Rebar::Size::bs10, Rebar::Size::bs11, Rebar::Size::bs14, Rebar::Size::bs18 };

   // The bar size: rebar.bar_size (a number, e.g. "#5", "5", "No. 5"), else the nominal diameter of the bar type or the
   // bar, else the diameter as drawn. bsNone if no size has the diameter (within 0.1 in)
   Rebar::Size get_bar_size(IfcSchema::IfcReinforcingBar bar, std::optional<Float64> drawn_diameter)
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
               return size;
         }
      }

      std::optional<Float64> diameter;
      for (const auto& rel : bar.IsTypedBy())
      {
         auto bar_type = rel.RelatingType().as<IfcSchema::IfcReinforcingBarType>();
         if (bar_type && bar_type.NominalDiameter())
            diameter = ConvertToSysUnits<IfcSchema::IfcPositiveLengthMeasure>(units, *bar_type.NominalDiameter());
      }
      if (!diameter && bar.NominalDiameter())
         diameter = ConvertToSysUnits<IfcSchema::IfcPositiveLengthMeasure>(units, *bar.NominalDiameter());
      if (!diameter)
         diameter = drawn_diameter;
      if (!diameter)
         return Rebar::Size::bsNone;

      auto best = std::min_element(bar_sizes.begin(), bar_sizes.end(), [&diameter](auto a, auto b) {return fabs(nominal_diameter(a) - *diameter) < fabs(nominal_diameter(b) - *diameter); });
      return fabs(nominal_diameter(*best) - *diameter) <= WBFL::Units::ConvertToSysUnits(0.1, WBFL::Units::Measure::Inch) ? *best : Rebar::Size::bsNone;
   }

   // The bar material: the type from rebar.grade (A615 unless the text says A706 or A1035), the grade from rebar.fy, else a
   // grade number in rebar.grade. nullopt if the model says neither
   std::optional<BarMaterial> get_bar_material(IfcSchema::IfcReinforcingBar bar)
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
}

Float64 nominal_diameter(Rebar::Size size)
{
   const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(Rebar::Type::A615, Rebar::Grade60, size);
   return pRebar ? pRebar->GetNominalDimension() : 0.0;
}

std::optional<BarMaterial> most_common_material(const std::map<BarMaterial, IndexType>& materials, bool& bMixed)
{
   bMixed = 1 < materials.size();
   if (materials.empty())
      return std::nullopt;
   auto most = std::max_element(materials.begin(), materials.end(), [](const auto& a, const auto& b) {return a.second < b.second; });
   if (!WBFL::LRFD::RebarPool::GetInstance()->GetRebar(most->first.first, most->first.second, Rebar::Size::bs5))
      return std::nullopt; // not a material PGSuper has
   return most->first;
}

void for_each_segment(CBridgeDescription2& bridge_desc, const std::function<void(CPrecastSegmentData*)>& fn)
{
   for (GroupIndexType grpIdx = 0; grpIdx < bridge_desc.GetGirderGroupCount(); grpIdx++)
   {
      auto* pGroup = bridge_desc.GetGirderGroup(grpIdx);
      for (GirderIndexType gdrIdx = 0; gdrIdx < pGroup->GetGirderCount(); gdrIdx++)
      {
         auto* pGirder = pGroup->GetGirder(gdrIdx);
         for (SegmentIndexType segIdx = 0; segIdx < pGirder->GetSegmentCount(); segIdx++)
            fn(pGirder->GetSegment(segIdx));
      }
   }
}

GirderBars::GirderBars(ifcopenshell::file& file, const GirderLayout& layout) :
   m_File(file), m_Geometry(file, layout)
{
   const auto& reader = CIfcImporter::GetTargetReader();
   if (!reader.HasSelector(ElementKind::Rebar))
   {
      m_Reason = "the mapping table has no rebar role";
      return;
   }

   auto elements = reader.Select(ElementKind::Rebar, file);
   if (elements.empty())
   {
      m_Reason = "the model has no reinforcing bars (" + reader.GetTable().GetSelector(ElementKind::Rebar)->Describe() + ")";
      return;
   }
   m_bHasBars = true;

   IndexType nWithoutGeometry = 0;
   for (const auto& object : elements)
   {
      auto element = object.as<IfcSchema::IfcReinforcingBar>();
      std::vector<GirderReinforcement> centerlines; // more than one for a group of bars
      auto result = element ? m_Geometry.Read(element, centerlines) : GirderReinforcementReader::Status::NoGeometry;
      if (result == GirderReinforcementReader::Status::NoGeometry)
         nWithoutGeometry++;
      if (result != GirderReinforcementReader::Status::Read)
         continue; // bars of other elements (deck, piers) aren't in a girder

      // the size and material are the element's (the bars of a group are of one type)
      auto size = get_bar_size(element, centerlines.front().diameter);
      auto material = get_bar_material(element);
      for (const auto& centerline : centerlines)
         m_Bars[centerline.beam_id].push_back({ element, centerline, size, material });
   }

   if (0 < nWithoutGeometry)
      WBFL::System::Logger::Info(std::to_string(nWithoutGeometry) + " reinforcing bars were skipped: their geometry has no straight pieces that make a centerline (circular profiles extruded along a line, or swept disks along a polyline or an indexed poly curve).");
}

const std::vector<GirderBar>& GirderBars::GetBars(int beam_id) const
{
   static const std::vector<GirderBar> none;
   auto found = m_Bars.find(beam_id);
   return found == m_Bars.end() ? none : found->second;
}

std::string GirderBars::GetGirderName(int beam_id) const
{
   auto beam = m_File.instance_by_id(beam_id).as<IfcSchema::IfcBeam>();
   return beam ? beam.Name().value_or("unnamed beam") : "unnamed beam";
}

CPrecastSegmentData* GirderBars::GetSegment(CBridgeDescription2& bridge_desc, const GirderCenterline& girder)
{
   return bridge_desc.GetGirderGroup(girder.girder_key.groupIndex)->GetGirder(girder.girder_key.girderIndex)->GetSegment(0);
}
