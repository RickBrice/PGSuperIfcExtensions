///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2025  Washington State Department of Transportation
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
#pragma once

#include <IFace\Bridge.h>
#include <IFace\AnalysisResults.h>
#include <IFace\Intervals.h>

#include "PsetEnum.h"

template <typename Schema>
void Add_TPF_Classification(IfcHierarchyHelper<Schema>& file)
{
   // we are using the TPFBridge bSDD for classifications
   auto classification = new Schema::IfcClassification(
      std::string("TPF Bridge")/*Source*/,
      std::string("2") /*Edition*/,
      std::string("2024-08-13") /*EditionDate*/,
      std::string("TPF Bridge (USA)"),
      boost::none /*Description*/,
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2") /*Specification*/,
      boost::none /*ReferenceTokens*/);
   file.addEntity(classification);

   auto project = file.getSingle<typename Schema::IfcProject>();

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr projects(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   projects->push(project);

   auto rel_associates_classification = new Schema::IfcRelAssociatesClassification(
      IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, projects, classification);
   file.addEntity(rel_associates_classification);
}

template <typename Schema>
void Create_Pset_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   auto project = file.getSingle<typename Schema::IfcProject>();

   // 5.1.8.1 PEnum_ProjectType
   std::vector<std::string> enum_values{ "MODIFICAITON","NEWBUILD","OPERATIONMAINTENANCE","RENOVATION","REPAIR" };
   auto project_type_enum = createPropertyEnumeration<Schema>("PEnum_ProjectType", enum_values);
   auto project_type_property = createPropertyEnumeratedValue<Schema>("ProjectType", project_type_enum, "NEWBUILD");

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(project_type_property);

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ProjectCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_projects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_projects->push(project);

   auto project_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_projects, property_set);
   file.addEntity(project_properties);
}

template <typename Schema>
void Create_Pset_TPFBridge_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   auto project = file.getSingle<typename Schema::IfcProject>();

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("ContractNumber"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("DesignNumber"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("ProjectNumber"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("ProjectWebsite"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_ProjectCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_projects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_projects->push(project);

   auto project_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_projects, property_set);
   file.addEntity(project_properties);
}

template <typename Schema>
void Create_Pset_BridgeCommon(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
{
   // 5.4.8.4 PEnum_StructureIndicator
   std::vector<std::string> enum_values{ "COATED","COMPOSITE","HOMOGENEOUS" };
   auto penum = createPropertyEnumeration<Schema>("PEnum_StructureIndicator", enum_values);
   auto property = createPropertyEnumeratedValue<Schema>("StructureIndicator", penum, "COMPOSITE");

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(property);

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_BridgeCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_bridges(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_bridges->push(bridge);

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_Pset_TPFBridge_BridgeCommon(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcBridge* bridge)
{
   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nSpans = pBridge->GetSpanCount();
   auto nPiers = pBridge->GetPierCount();

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("tpfBridge_NumberOfSpans"), boost::none, new Schema::IfcInteger((int)nSpans), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("tpfBridge_NumberOfSupports"), boost::none, new Schema::IfcInteger((int)nPiers), nullptr));

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_BridgeCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_bridges(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_bridges->push(bridge);

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_Pset_TPFBridge_RailingCommon(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*> railings)
{
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("tpfBridge_MASHCompliantRailing"), boost::none, new Schema::IfcBoolean(true), nullptr));

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_RailingCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_railings(new aggregate_of<typename Schema::IfcObjectDefinition>());
   for (auto& railing : railings)
   {
      related_railings->push(railing);
   }

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_railings, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_Pset_TPFBridge_GirderCommon(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options, const CSegmentKey& segmentKey, typename Schema::IfcElement* segment)
{
   GET_IFACE2(pBroker, IMaterials, pMaterials);
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   auto releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* displacement_unit = nullptr;

   if (pDisplayUnits->GetUnitMode() == eafTypes::umUS)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);
   }

   auto fci = pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx);
   auto fc = pMaterials->GetSegmentFc28(segmentKey);
   auto fpj = pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent);
   if (pDisplayUnits->GetUnitMode() == eafTypes::umUS)
   {
      fci = WBFL::Units::ConvertFromSysUnits(fci, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fpj = WBFL::Units::ConvertFromSysUnits(fpj, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   list_of_properties->push(new Schema::IfcPropertySingleValue(
      std::string("tpfBridge_GirderStrengthatTimeOfPrestress"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_ConcreteStrengthatTimeOfPrestressing"),
      new Schema::IfcInteger((int)fci), stress_unit));

   list_of_properties->push(new Schema::IfcPropertySingleValue(
      std::string("tpfBridge_ConceteStrengthat28Days"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_ConceteStrengthat28Days"),
      new Schema::IfcInteger((int)fc), stress_unit));

   list_of_properties->push(new Schema::IfcPropertySingleValue(
      std::string("tpfBridge_JackingForce"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_JackingForce"),
      new Schema::IfcReal(fpj), stress_unit));

   // https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/TPFBridge_GirderCommon
   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_GirderCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_segments->push(segment);

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, property_set);
   file.addEntity(related_properties);

   if (options.include_camber)
   {
      GET_IFACE2(pBroker, IGirder, pGirder);

      typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());

      Float64 precamber = pGirder->GetPrecamber(segmentKey);
      if (!IsZero(precamber))
      {
         if (pDisplayUnits->GetUnitMode() == eafTypes::umUS)
         {
            precamber = WBFL::Units::ConvertFromSysUnits(precamber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
         }
         list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("tpfBridge_BuiltInCamber"), boost::none, new Schema::IfcReal(precamber), displacement_unit));
      }

      GET_IFACE2(pBroker, IPointOfInterest, pPoi);
      PoiList vPoi;
      pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
      CHECK(vPoi.size() == 1);
      const pgsPointOfInterest& poiMS = vPoi.front();

      GET_IFACE2(pBroker, IProductForces, pProduct);
      auto bat = pProduct->GetBridgeAnalysisType(pgsTypes::Minimize); // minimize because we want the greatest downward deflection

      Float64 ps = pProduct->GetDeflection(releaseIntervalIdx, pgsTypes::pftPretension, poiMS, bat, rtCumulative, false);
      Float64 girder = pProduct->GetDeflection(releaseIntervalIdx, pgsTypes::pftGirder, poiMS, bat, rtCumulative, false);
      Float64 camber = ps + girder;
      if (pDisplayUnits->GetUnitMode() == eafTypes::umUS)
      {
         camber = WBFL::Units::ConvertFromSysUnits(camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
      }
      list_of_properties->push(new Schema::IfcPropertySingleValue(
         std::string("tpfBridge_CamberatPrestressingRelease"),
         std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_CamberatPrestressingRelease"),
         new Schema::IfcReal(camber + precamber), displacement_unit));

      auto lastIntervalIdx = pIntervals->GetIntervalCount() - 1;
      auto lastCompositeIntervalIdx = pIntervals->GetLastCompositeDeckInterval();
      GET_IFACE2(pBroker, ICombinedForces, pCombined);
      Float64 dc_final = pCombined->GetDeflection(lastIntervalIdx, lcDC, poiMS, bat, rtCumulative);
      Float64 dc_composite = pCombined->GetDeflection(lastCompositeIntervalIdx, lcDC, poiMS, bat, rtCumulative);
      Float64 dw_final = pCombined->GetDeflection(lastIntervalIdx, lcDW, poiMS, bat, rtCumulative);
      Float64 dw_composite = pCombined->GetDeflection(lastCompositeIntervalIdx, lcDW, poiMS, bat, rtCumulative);
      Float64 d = (dc_final - dc_composite) + (dw_final - dw_composite);
      list_of_properties->push(new Schema::IfcPropertySingleValue(
         std::string("tpfBridge_DeflectionDuetoCompositLoads"),
         std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_DeflectionDuetoCompositLoads"),
         new Schema::IfcReal(d), displacement_unit));

      // https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/TPFBridge_MemberCamber
      auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_MemberCamber"), boost::none, list_of_properties);

      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
      related_segments->push(segment);

      auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, property_set);
      file.addEntity(related_properties);
   }
}

template <typename Schema>
void Create_Pset_TPFBridge_ReinforcementCommon(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* tendon, Float64 Pjack, bool bDebonded, Float64 ldb)
{
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());

   list_of_properties->push(new Schema::IfcPropertySingleValue(
      std::string("tpfBridge_TendonJackingForce"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonJackingForce"),
      new Schema::IfcReal(Pjack),
      nullptr));

   list_of_properties->push(new Schema::IfcPropertySingleValue(
      std::string("tpfBridge_TendonBonding"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonBonding"),
#pragma Reminder("bSDD - should the applicable values be the string value or the URI reference to the string value?")
      // not sure if this should be URI reference or "Debonded" "Bonded" both are strings
      new Schema::IfcURIReference(bDebonded ? "https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonBonding/value/TendonBondingDebonded" : "https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonBonding/value/TendonBondingBonded"),
      nullptr));

   if (bDebonded)
   {
      std::ostringstream os;
      os << ldb;
      list_of_properties->push(new Schema::IfcPropertySingleValue(
         std::string("tpfBridge_TendonDebondedLength"),
         std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonDebondedLength"),
#pragma Reminder("bSDD - why is debond length a string?")
         new Schema::IfcText(os.str()), // this could be IfcIdentifier, IfcLabel, or IfcText none actually represent a value
         nullptr));
   }

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_ReinforcementCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_tendons(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_tendons->push(tendon);

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_tendons, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Classify_TPFBridge(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
{
   auto classification = file.getSingle<typename Schema::IfcClassification>();

   auto classification_reference = new Schema::IfcClassificationReference(
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Bridge"),
      std::string("tpfBridge_Bridge") /*Identification*/,
      std::string("Bridge") /*Name*/,
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);
   file.addEntity(classification_reference);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_bridges(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   related_bridges->push(bridge);

   auto related_classes = new Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, classification_reference);

   file.addEntity(related_classes);
}

template <typename Schema>
void Classify_TPFBridgePart(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* part, const std::string& uri, const std::string& code, const std::string& name, const std::string& part_type)
{
   std::vector<typename Schema::IfcProduct*> parts{ part };
   Classify_TPFBridgeParts(file, parts, uri, code, name, part_type);
}

template <typename Schema>
void Classify_TPFBridgeParts(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& parts, const std::string& uri, const std::string& code, const std::string& name, const std::string& part_type)
{
   if (parts.size() == 0)
      return;

   auto classification = file.getSingle<typename Schema::IfcClassification>();

   auto classification_reference = new Schema::IfcClassificationReference(
      uri, /*Class identifier (uri) = Location*/
      code /*Class code = Identification*/,
      name,/*Class name = name*/
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);
   file.addEntity(classification_reference);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_parts(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   for (auto& part : parts)
   {
      related_parts->push(part);
   }

   auto related_classes = new Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_parts, classification_reference);
   file.addEntity(related_classes);
}

template <typename Schema>
void Classify_TPFSuperstructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* superstructure)
{
   //Classify_TPFBridgePart(file, superstructure, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_BridgeSuperstructure"), std::string("tpfBridge_BridgeSuperstructure"), std::string("Bridge Superstructure"), std::string("IfcBridgePart.SUPERSTRUCTURE"));
   Classify_TPFBridgePart(file, superstructure, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_BridgeSuperstructure"), std::string("tpfBridge_BridgeSuperstructure"), std::string("IfcBridgePartSUPERSTRUCTURE"), std::string("IfcBridgePart.SUPERSTRUCTURE"));
}

template <typename Schema>
void Classify_TPFSubstructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* substructure)
{
   //Classify_TPFBridgePart(file, substructure, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_BridgeSubstructure"), std::string("tpfBridge_BridgeSubstructure"), std::string("Bridge Substructure"), std::string("IfcBridgePart.SUBSTRUCTURE"));
   Classify_TPFBridgePart(file, substructure, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_BridgeSubstructure"), std::string("tpfBridge_BridgeSubstructure"), std::string("IfcBridgePartSUBSTRUCTURE"), std::string("IfcBridgePart.SUBSTRUCTURE"));
}

template <typename Schema>
void Classify_TPFAbutment(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* abutment)
{
   //Classify_TPFBridgePart(file, abutment, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_AbutmentSpatial"), std::string("tpfBridge_AbutmentSpatial"), std::string("Abutment (Spatial)"), std::string("IfcBridgePart.ABUTMENT"));
   Classify_TPFBridgePart(file, abutment, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_AbutmentSpatial"), std::string("tpfBridge_AbutmentSpatial"), std::string("IfcBridgePartABUTMENT"), std::string("IfcBridgePart.ABUTMENT"));
}

template <typename Schema>
void Classify_TPFPier(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* pier)
{
   //Classify_TPFBridgePart(file, pier, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_PierSpatial"), std::string("tpfBridge_PierSpatial"), std::string("Pier (Spatial)"), std::string("IfcBridgePart.PIER"));
   Classify_TPFBridgePart(file, pier, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_PierSpatial"), std::string("tpfBridge_PierSpatial"), std::string("IfcBridgePartPIER"), std::string("IfcBridgePart.PIER"));
}

template <typename Schema>
void Classify_TPFFoundation(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* foundation)
{
   //Classify_TPFBridgePart(file, foundation, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Foundation"), std::string("tpfBridge_Foundation"), std::string("Foundation"), std::string("IfcBridgePart.FOUNDATION"));
   Classify_TPFBridgePart(file, foundation, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Foundation"), std::string("tpfBridge_Foundation"), std::string("IfcBridgePartFOUNDATION"), std::string("IfcBridgePart.FOUNDATION"));
}

template <typename Schema>
void Classify_TPFRailings(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& railings)
{
   //Classify_TPFBridgeParts(file, railings, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Railing"), std::string("tpfBridge_Railing"), std::string("Railing"), std::string("IfcRailing.BALUSTRADE"));
   Classify_TPFBridgeParts(file, railings, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Railing"), std::string("tpfBridge_Railing"), std::string("IfcRailingBALUSTRADE"), std::string("IfcRailing.BALUSTRADE"));
}

template <typename Schema>
void Classify_TPFGirders(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& girders)
{
   Classify_TPFBridgeParts(file, girders, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Girder"), std::string("tpfBridge_Girder"), std::string("IfcElementAssemblyGIRDER"), std::string("IfcElementAssembly.GIRDER"));
}

template <typename Schema>
void Classify_TPFPrecastGirderElements(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& girders)
{
   //Classify_TPFBridgeParts(file, girders, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_GirderPrestressedConcrete"), std::string("tpfBridge_GirderPrestressedConcrete"), std::string("IfcElementAssemblyGIRDER"), std::string("IfcElementAssembly.GIRDER"));
   Classify_TPFBridgeParts(file, girders, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_GirderPrestressedConcrete"), std::string("tpfBridge_GirderPrestressedConcrete"), std::string("IfcBeamGIRDER_SEGMENT"), std::string("IfcBeam.GIRDER_SEGMENT"));
}

template <typename Schema>
void Classify_TPFPrestressing(IfcHierarchyHelper<Schema>& file, typename Schema::IfcTendonType* tendon)
{
   auto classification = file.getSingle<typename Schema::IfcClassification>();

   auto classification_reference = new Schema::IfcClassificationReference(
      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Prestressing"),
      std::string("tpfBridge_Prestressing") /*Identification*/,
      std::string("Prestressing") /*Name*/,
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);
   file.addEntity(classification_reference);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_tendons(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   related_tendons->push(tendon);

   auto related_classes = new Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_tendons, classification_reference);

   file.addEntity(related_classes);
}
