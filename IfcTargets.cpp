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
#include "IfcTargets.h"
#include "IfcExporter.h"
#include "SteelSpecifications.h"

#include <IFace/Bridge.h>
#include <IFace/Project.h>
#include <IFace/Intervals.h>
#include <IFace/PointOfInterest.h>
#include <IFace/AnalysisResults.h>
#include <EAF/EAFDisplayUnits.h>
#include <PsgLib/BridgeDescription2.h>
#include <PsgLib/GirderLabel.h>

#include <algorithm>
#include <sstream>

namespace
{
   using DU = ExportUnit;

   Float64 pier_skew(std::shared_ptr<WBFL::EAF::Broker> pBroker, PierIndexType pierIdx)
   {
      GET_IFACE2(pBroker, IBridge, pBridge);
      CComPtr<IAngle> angle;
      pBridge->GetPierSkew(pierIdx, &angle);
      Float64 skew;
      angle->get_Value(&skew);
      return skew;
   }

   // f'c formatted in display units, as Pset_ConcreteElementGeneral.StrengthClass.
   // The trailing newline reproduces the exporter it replaces (see devdocs/MappingTablesDesign.md, export findings)
   std::string strength_class(std::shared_ptr<WBFL::EAF::Broker> pBroker, Float64 fc)
   {
      USES_CONVERSION;
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      std::ostringstream os;
      os << T2A(::FormatDimension(fc, pDisplayUnits->GetStressUnit())) << std::endl;
      return os.str();
   }

   // the released segment mid-span point of interest, where cambers are reported
   pgsPointOfInterest released_midspan(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey)
   {
      GET_IFACE2(pBroker, IPointOfInterest, pPoi);
      PoiList vPoi;
      pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
      return vPoi.front();
   }

   Float64 deck_fc(const ExportContext& c) { GET_IFACE2(c.broker, IMaterials, pMaterials); return pMaterials->GetDeckFc28(); }
   Float64 deck_max_aggregate_size(const ExportContext& c) { GET_IFACE2(c.broker, IMaterials, pMaterials); return pMaterials->GetDeckMaxAggrSize(); }
}

const std::vector<TargetDef>& GetTargetDefs()
{
   static const std::vector<TargetDef> targets{
      //
      // Bridge
      //
      { .name = "bridge.number_of_spans", .element = ElementKind::Bridge, .kind = ValueKind::Count, .description = "number of spans",
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue((Int64)pBridge->GetSpanCount()); } },
      { .name = "bridge.length", .element = ElementKind::Bridge, .kind = ValueKind::Length, .description = "bridge length", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetLength()); } },
      { .name = "bridge.roadway_width", .element = ElementKind::Bridge, .kind = ValueKind::Length, .description = "roadway width at the start of the bridge", .display_unit = DU::SpanLength, .display_round = 0.1,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetCurbToCurbWidth(0.0)); } },
      { .name = "bridge.start_station", .element = ElementKind::Bridge, .kind = ValueKind::Length, .description = "station of the start of the bridge", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetPierStation(0)); } },
      { .name = "bridge.end_station", .element = ElementKind::Bridge, .kind = ValueKind::Length, .description = "station of the end of the bridge", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetPierStation(pBridge->GetSpanCount())); } },
      { .name = "bridge.max_skew", .element = ElementKind::Bridge, .kind = ValueKind::Angle, .description = "largest pier skew", .display_unit = DU::Angle, .display_round = 0.0001,
        .get = [](const ExportContext& c)
         {
            // the last pier isn't included, as in the exporter this replaces
            GET_IFACE2(c.broker, IBridge, pBridge);
            Float64 skew = 0;
            for (PierIndexType i = 0; i < pBridge->GetSpanCount(); i++)
            {
               Float64 pier = pier_skew(c.broker, i);
               if (std::fabs(skew) < std::fabs(pier))
                  skew = pier;
            }
            return ExportValue(skew);
         } },
      { .name = "bridge.design_live_load", .element = ElementKind::Bridge, .kind = ValueKind::Text, .description = "design vehicular live load",
        .get = [](const ExportContext& c)
         {
            USES_CONVERSION;
            GET_IFACE2(c.broker, ILiveLoads, pLiveLoads);
            if (!pLiveLoads->IsLiveLoadDefined(pgsTypes::lltDesign))
               return ExportValue::Absent();
            return ExportValue(std::string(T2A(pLiveLoads->GetLiveLoadNames(pgsTypes::lltDesign).front().c_str())));
         } },

      //
      // Piers and abutments
      //
      { .name = "pier.station", .element = ElementKind::Pier, .kind = ValueKind::Length, .description = "station of the CL pier", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetPierStation(c.pier)); } },
      { .name = "pier.skew", .element = ElementKind::Pier, .kind = ValueKind::Angle, .description = "pier skew", .display_unit = DU::Angle, .display_round = 0.0001,
        .get = [](const ExportContext& c) { return ExportValue(pier_skew(c.broker, c.pier)); } },

      //
      // Girders (one IfcBeam per segment)
      //
      { .name = "girder.designation", .element = ElementKind::Girder, .kind = ValueKind::Text, .description = "girder designation (e.g. \"Span 1, Girder 2\")",
        .get = [](const ExportContext& c)
         {
            USES_CONVERSION;
            pgsAutoGirderLabel autoLabel;
            pgsGirderLabel::UseAlphaLabel(false);
            return ExportValue(std::string(T2A(SEGMENT_LABEL(c.segment))));
         } },
      { .name = "girder.type_names", .element = ElementKind::Girder, .kind = ValueKind::TextList, .description = "girder type name",
        .get = [](const ExportContext& c) { USES_CONVERSION; GET_IFACE2(c.broker, IBridgeDescription, pIBridgeDesc); return ExportValue(std::string(T2A(pIBridgeDesc->GetGirder(c.segment)->GetGirderName()))); } },
      { .name = "girder.family_name", .element = ElementKind::Girder, .kind = ValueKind::Text, .description = "girder family (e.g. I-Beam)",
        .get = [](const ExportContext& c) { USES_CONVERSION; GET_IFACE2(c.broker, IBridgeDescription, pIBridgeDesc); return ExportValue(std::string(T2A(pIBridgeDesc->GetBridgeDescription()->GetGirderFamilyName()))); } },
      { .name = "girder.fc", .element = ElementKind::Girder, .kind = ValueKind::Stress, .description = "girder concrete strength, f'c", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IMaterials, pMaterials); return ExportValue(pMaterials->GetSegmentFc28(c.segment)); } },
      { .name = "girder.fci", .element = ElementKind::Girder, .kind = ValueKind::Stress, .description = "girder concrete strength at release, f'ci", .display_unit = DU::Stress,
        .get = [](const ExportContext& c)
         {
            GET_IFACE2(c.broker, IIntervals, pIntervals);
            GET_IFACE2(c.broker, IMaterials, pMaterials);
            return ExportValue(pMaterials->GetSegmentFc(c.segment, pIntervals->GetPrestressReleaseInterval(c.segment)));
         } },
      { .name = "girder.strength_class", .element = ElementKind::Girder, .kind = ValueKind::Text, .description = "girder concrete strength class",
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IMaterials, pMaterials); return ExportValue(strength_class(c.broker, pMaterials->GetSegmentFc28(c.segment))); } },
      { .name = "girder.max_aggregate_size", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "girder concrete maximum aggregate size", .display_unit = DU::Deflection,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IMaterials, pMaterials); return ExportValue(pMaterials->GetSegmentMaxAggrSize(c.segment)); } },
      { .name = "girder.assembly_place", .element = ElementKind::Girder, .kind = ValueKind::Text, .description = "girder assembly place (e.g. FACTORY)" },
      { .name = "girder.casting_method", .element = ElementKind::Girder, .kind = ValueKind::Text, .description = "girder casting method (e.g. PRECAST)" },
      { .name = "girder.span", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "girder span length", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetSegmentSpanLength(c.segment)); } },
      { .name = "girder.slope", .element = ElementKind::Girder, .kind = ValueKind::Angle, .description = "girder slope angle", .display_unit = DU::Angle,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(atan(pBridge->GetSegmentSlope(c.segment))); } },
      { .name = "girder.roll", .element = ElementKind::Girder, .kind = ValueKind::Angle, .description = "girder roll angle", .display_unit = DU::Angle,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IGirder, pGirder); return ExportValue(atan(pGirder->GetOrientation(c.segment))); } },
      { .name = "girder.batter", .element = ElementKind::Girder, .kind = ValueKind::Angle, .description = "girder end batter", .display_unit = DU::Angle,
        .get = [](const ExportContext& c)
         {
            if (!c.options->batter_ends)
               return ExportValue(0.0);
            GET_IFACE2(c.broker, IBridge, pBridge);
            return ExportValue(atan(pBridge->GetSegmentSlope(c.segment)));
         } },
      { .name = "girder.bunk_point", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "girder hauling support location", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c)
         {
            GET_IFACE2(c.broker, IBridgeDescription, pBridgeDesc);
            const auto& hauling = pBridgeDesc->GetBridgeDescription()->GetGirderGroup(c.segment.groupIndex)->GetGirder(c.segment.girderIndex)->GetSegment(c.segment.segmentIndex)->HandlingData;
            return ExportValue(std::max(hauling.LeadingSupportPoint, hauling.TrailingSupportPoint));
         } },
      { .name = "girder.jacking_stress", .element = ElementKind::Girder, .kind = ValueKind::Stress, .description = "jacking stress of the permanent strands", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IStrandGeometry, pStrandGeom); return ExportValue(pStrandGeom->GetJackingStress(c.segment, pgsTypes::Permanent)); } },
      { .name = "girder.camber_ratio", .element = ElementKind::Girder, .kind = ValueKind::Ratio, .description = "camber at release divided by the girder length",
        .get = [](const ExportContext& c)
         {
            if (!c.options->include_camber)
               return ExportValue::NoValue();

            GET_IFACE2(c.broker, IGirder, pGirder);
            GET_IFACE2(c.broker, IIntervals, pIntervals);
            GET_IFACE2(c.broker, IProductForces, pProduct);
            GET_IFACE2(c.broker, IBridge, pBridge);
            auto poi = released_midspan(c.broker, c.segment);
            auto releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(c.segment);
            auto bat = pProduct->GetBridgeAnalysisType(pgsTypes::Minimize); // the greatest downward deflection
            Float64 ps = pProduct->GetDeflection(releaseIntervalIdx, pgsTypes::pftPretension, poi, bat, rtCumulative, false);
            Float64 girder = pProduct->GetDeflection(releaseIntervalIdx, pgsTypes::pftGirder, poi, bat, rtCumulative, false);
            return ExportValue((ps + girder + pGirder->GetPrecamber(c.segment)) / pBridge->GetSegmentPlanLength(c.segment));
         } },
      { .name = "girder.camber_at_release", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "camber at release", .display_unit = DU::Deflection,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, ICamber, pCamber); return ExportValue(pCamber->GetInitialCamber(released_midspan(c.broker, c.segment))); } },
      { .name = "girder.camber_after_losses", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "excess camber", .display_unit = DU::Deflection,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, ICamber, pCamber); return ExportValue(pCamber->GetExcessCamber(released_midspan(c.broker, c.segment), pgsTypes::CreepTime::Max)); } },
      { .name = "girder.screed_camber", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "screed camber", .display_unit = DU::Deflection,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, ICamber, pCamber); return ExportValue(pCamber->GetScreedCamber(released_midspan(c.broker, c.segment), pgsTypes::CreepTime::Max)); } },
      // Qto_BeamBaseQuantities. Simple sections are assumed (no end blocks or variable depth)
      { .name = "girder.length", .element = ElementKind::Girder, .kind = ValueKind::Length, .description = "girder plan length", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c) { GET_IFACE2(c.broker, IBridge, pBridge); return ExportValue(pBridge->GetSegmentPlanLength(c.segment)); } },
      { .name = "girder.cross_section_area", .element = ElementKind::Girder, .kind = ValueKind::Area, .description = "girder cross section area at mid-span", .display_unit = DU::SmallArea,
        .get = [](const ExportContext& c)
         {
            GET_IFACE2(c.broker, IIntervals, pIntervals);
            GET_IFACE2(c.broker, ISectionProperties, pSectProps);
            return ExportValue(pSectProps->GetAg(pIntervals->GetPrestressReleaseInterval(c.segment), released_midspan(c.broker, c.segment)));
         } },
      { .name = "girder.outer_surface_area", .element = ElementKind::Girder, .kind = ValueKind::Area, .description = "girder outer surface area (length times perimeter)", .display_unit = DU::BigArea,
        .get = [](const ExportContext& c)
         {
            GET_IFACE2(c.broker, IBridge, pBridge);
            GET_IFACE2(c.broker, ISectionProperties, pSectProps);
            return ExportValue(pBridge->GetSegmentPlanLength(c.segment) * pSectProps->GetPerimeter(released_midspan(c.broker, c.segment)));
         } },
      { .name = "girder.gross_weight", .element = ElementKind::Girder, .kind = ValueKind::Mass, .description = "girder mass (weight and mass are the same in IfcQuantityWeight)", .display_unit = DU::Mass,
        .get = [](const ExportContext& c)
         {
            GET_IFACE2(c.broker, ISectionProperties, pSectProps);
            return ExportValue(pSectProps->GetSegmentWeight(c.segment) / WBFL::Units::System::GetGravitationalAcceleration());
         } },

      //
      // Deck
      //
      { .name = "deck.gross_depth", .element = ElementKind::Deck, .kind = ValueKind::Length, .description = "deck gross depth" },
      { .name = "deck.fc", .element = ElementKind::Deck, .kind = ValueKind::Stress, .description = "deck concrete strength, f'c", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { return ExportValue(deck_fc(c)); } },
      { .name = "deck.strength_class", .element = ElementKind::Deck, .kind = ValueKind::Text, .description = "deck concrete strength class",
        .get = [](const ExportContext& c) { return ExportValue(strength_class(c.broker, deck_fc(c))); } },
      { .name = "deck.max_aggregate_size", .element = ElementKind::Deck, .kind = ValueKind::Length, .description = "deck concrete maximum aggregate size", .display_unit = DU::Deflection,
        .get = [](const ExportContext& c) { return ExportValue(deck_max_aggregate_size(c)); } },

      //
      // Barriers (PGSuper assumes they are the same concrete as the deck)
      //
      { .name = "barrier.fc", .element = ElementKind::Barrier, .kind = ValueKind::Stress, .description = "barrier concrete strength, f'c", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { return ExportValue(deck_fc(c)); } },
      { .name = "barrier.strength_class", .element = ElementKind::Barrier, .kind = ValueKind::Text, .description = "barrier concrete strength class",
        .get = [](const ExportContext& c) { return ExportValue(strength_class(c.broker, deck_fc(c))); } },
      { .name = "barrier.max_aggregate_size", .element = ElementKind::Barrier, .kind = ValueKind::Length, .description = "barrier concrete maximum aggregate size", .display_unit = DU::Deflection,
        .get = [](const ExportContext& c) { return ExportValue(deck_max_aggregate_size(c)); } },

      //
      // Strands (IfcTendon) and their material
      //
      { .name = "tendon.fy", .element = ElementKind::Tendon, .kind = ValueKind::Stress, .description = "strand yield strength", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { return ExportValue(c.strand->GetYieldStrength()); } },
      { .name = "tendon.fpu", .element = ElementKind::Tendon, .kind = ValueKind::Stress, .description = "strand ultimate strength", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { return ExportValue(c.strand->GetUltimateStrength()); } },
      { .name = "tendon.grade", .element = ElementKind::Tendon, .kind = ValueKind::Text, .description = "strand grade (e.g. ASTM A416 Grade 270)",
        .get = [](const ExportContext& c)
         {
            // ASTM A416 is for low relaxation Grade 250 and 270 strand. PGSuper also has stress relieved strand and Grade 300
            USES_CONVERSION;
            return ExportValue(std::string("ASTM A416 Grade ") + T2A(WBFL::Materials::PsStrand::GetGrade(c.strand->GetGrade(), true/*US units*/).c_str()));
         } },
      { .name = "tendon.corrosion_treatment", .element = ElementKind::Tendon, .kind = ValueKind::Text, .description = "strand corrosion treatment (NONE or EPOXYCOATED)",
        .get = [](const ExportContext& c) { return ExportValue(std::string(c.strand->GetCoating() == WBFL::Materials::PsStrand::Coating::None ? "NONE" : "EPOXYCOATED")); } },
      { .name = "tendon.status", .element = ElementKind::Tendon, .kind = ValueKind::Text, .description = "strand status (TEMPORARY for temporary strands)",
        .get = [](const ExportContext& c) { return c.strand_type == pgsTypes::Temporary ? ExportValue(std::string("TEMPORARY")) : ExportValue::Absent(); } },
      { .name = "tendon.debond_start", .element = ElementKind::Tendon, .kind = ValueKind::Length, .description = "debond length at the start of a debonded strand", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c)
         {
            // only with usBridge classification, as in the exporter this replaces
            if (!c.options->classify)
               return ExportValue::Absent();
            GET_IFACE2(c.broker, IStrandGeometry, pStrandGeom);
            Float64 start, end;
            return pStrandGeom->IsStrandDebonded(c.segment, c.strand_index, c.strand_type, nullptr, &start, &end) ? ExportValue(start) : ExportValue::Absent();
         } },
      { .name = "tendon.debond_end", .element = ElementKind::Tendon, .kind = ValueKind::Length, .description = "debond length at the end of a debonded strand", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c)
         {
            if (!c.options->classify)
               return ExportValue::Absent();
            GET_IFACE2(c.broker, IStrandGeometry, pStrandGeom);
            Float64 start, end;
            return pStrandGeom->IsStrandDebonded(c.segment, c.strand_index, c.strand_type, nullptr, &start, &end) ? ExportValue(end) : ExportValue::Absent();
         } },
      { .name = "tendon.debond_midspan_to_start", .element = ElementKind::Tendon, .kind = ValueKind::Length, .description = "debonded length from mid-span toward the start of a temporary strand", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c)
         {
            // assumes the WSDOT detail: temporary top strands are debonded except for their end 10 ft
            if (c.strand_type != pgsTypes::Temporary)
               return ExportValue::Absent();
            GET_IFACE2(c.broker, IBridge, pBridge);
            return ExportValue(pBridge->GetSegmentPlanLength(c.segment) / 2 - WBFL::Units::ConvertToSysUnits(10.0, WBFL::Units::Measure::Feet));
         } },
      { .name = "tendon.debond_midspan_to_end", .element = ElementKind::Tendon, .kind = ValueKind::Length, .description = "debonded length from mid-span toward the end of a temporary strand", .display_unit = DU::SpanLength,
        .get = [](const ExportContext& c)
         {
            if (c.strand_type != pgsTypes::Temporary)
               return ExportValue::Absent();
            GET_IFACE2(c.broker, IBridge, pBridge);
            return ExportValue(pBridge->GetSegmentPlanLength(c.segment) / 2 - WBFL::Units::ConvertToSysUnits(10.0, WBFL::Units::Measure::Feet));
         } },

      //
      // Reinforcing bar material
      //
      { .name = "rebar.fy", .element = ElementKind::Rebar, .kind = ValueKind::Stress, .description = "reinforcing bar yield strength", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { return ExportValue(c.rebar->GetYieldStrength()); } },
      { .name = "rebar.fu", .element = ElementKind::Rebar, .kind = ValueKind::Stress, .description = "reinforcing bar ultimate strength", .display_unit = DU::Stress,
        .get = [](const ExportContext& c) { return ExportValue(c.rebar->GetUltimateStrength()); } },
      { .name = "rebar.elongation", .element = ElementKind::Rebar, .kind = ValueKind::Ratio, .description = "reinforcing bar elongation at fracture",
        .get = [](const ExportContext& c) { return ExportValue(c.rebar->GetElongation()); } },
      { .name = "rebar.grade", .element = ElementKind::Rebar, .kind = ValueKind::Text, .description = "reinforcing bar material (e.g. AASHTO M31 (A615) - Grade 60)",
        .get = [](const ExportContext& c) { USES_CONVERSION; return ExportValue(std::string(T2A(WBFL::LRFD::RebarPool::GetMaterialName(c.rebar->GetType(), c.rebar->GetGrade()).c_str()))); } },
      { .name = "rebar.specification", .element = ElementKind::Rebar, .kind = ValueKind::Text, .description = "reinforcing bar specification (e.g. ASTM A615 (AASHTO M31))",
        .get = [](const ExportContext& c) { return ExportValue(GetRebarSpecification(c.rebar)); } },
      { .name = "rebar.specification_edition", .element = ElementKind::Rebar, .kind = ValueKind::Text, .description = "edition of the reinforcing bar specification",
        .get = [](const ExportContext& c) { return ExportValue(GetRebarSpecificationEdition(c.rebar)); } },

      //
      // Bearings
      //
      { .name = "bearing.fixed_x", .element = ElementKind::Bearing, .kind = ValueKind::Boolean, .description = "bearing fixed along the girder" },
      { .name = "bearing.fixed_y", .element = ElementKind::Bearing, .kind = ValueKind::Boolean, .description = "bearing fixed across the girder" },
   };
   return targets;
}

const TargetDef* FindTargetDef(std::string_view name)
{
   const auto& targets = GetTargetDefs();
   auto found = std::find_if(targets.begin(), targets.end(), [name](const auto& target) {return target.name == name; });
   return found == targets.end() ? nullptr : &(*found);
}

namespace
{
   const std::vector<std::pair<ElementKind, std::string_view>>& element_role_names()
   {
      static const std::vector<std::pair<ElementKind, std::string_view>> names{
         { ElementKind::Project, "project" },
         { ElementKind::Site, "site" },
         { ElementKind::Bridge, "bridge" },
         { ElementKind::Superstructure, "superstructure" },
         { ElementKind::Substructure, "substructure" },
         { ElementKind::DeckPart, "deck_part" },
         { ElementKind::Pier, "pier" },
         { ElementKind::Abutment, "abutment" },
         { ElementKind::Foundation, "foundation" },
         { ElementKind::Alignment, "alignment" },
         { ElementKind::Referent, "referent" },
         { ElementKind::Girder, "girder" },
         { ElementKind::GirderAssembly, "girder_assembly" },
         { ElementKind::ClosureJoint, "closure_joint" },
         { ElementKind::Tendon, "tendon" },
         { ElementKind::TendonBundle, "tendon_bundle" },
         { ElementKind::Rebar, "rebar" },
         { ElementKind::ReinforcementCage, "reinforcement_cage" },
         { ElementKind::Deck, "deck" },
         { ElementKind::Haunch, "haunch" },
         { ElementKind::Bearing, "bearing" },
         { ElementKind::Barrier, "barrier" },
      };
      return names;
   }
}

std::string_view GetElementRoleName(ElementKind kind)
{
   for (const auto& [k, name] : element_role_names())
   {
      if (k == kind)
         return name;
   }
   ASSERT(false); // every kind should have a name
   return "";
}

bool GetElementKind(std::string_view role_name, ElementKind& kind)
{
   for (const auto& [k, name] : element_role_names())
   {
      if (name == role_name)
      {
         kind = k;
         return true;
      }
   }
   return false;
}

ElementKind GetTargetElement(ElementKind role)
{
   return role == ElementKind::Abutment ? ElementKind::Pier : role;
}

bool HasUnit(ValueKind kind)
{
   return kind == ValueKind::Stress || kind == ValueKind::Length || kind == ValueKind::Angle || kind == ValueKind::Force || kind == ValueKind::Area || kind == ValueKind::Mass;
}

bool IsNumeric(ValueKind kind)
{
   return HasUnit(kind) || kind == ValueKind::Ratio || kind == ValueKind::Count;
}

std::string_view GetValueKindName(ValueKind kind)
{
   switch (kind)
   {
   case ValueKind::Stress: return "stress";
   case ValueKind::Length: return "length";
   case ValueKind::Angle: return "angle";
   case ValueKind::Force: return "force";
   case ValueKind::Area: return "area";
   case ValueKind::Mass: return "mass";
   case ValueKind::Ratio: return "ratio";
   case ValueKind::Count: return "count";
   case ValueKind::Boolean: return "boolean";
   case ValueKind::Text: return "text";
   case ValueKind::TextList: return "text list";
   }
   ASSERT(false);
   return "";
}

std::string FormatTargetValue(const TargetValue& value)
{
   std::ostringstream os;
   std::visit([&os](const auto& v)
      {
         using T = std::decay_t<decltype(v)>;
         if constexpr (std::is_same_v<T, bool>)
            os << (v ? "true" : "false");
         else if constexpr (std::is_same_v<T, std::string>)
            os << "\"" << v << "\"";
         else
            os << v;
      }, value);
   return os.str();
}
