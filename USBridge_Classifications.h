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

#include <IFace/Tools.h>
#include <IFace/PointOfInterest.h>
#include <IFace\Bridge.h>
#include <IFace\AnalysisResults.h>
#include <IFace\Intervals.h>
#include "Units.h"
#include "PsetEnum.h"

template <typename Schema>
void Add_USBridge_Classification(IfcHierarchyHelper<Schema>& file)
{
   // we are using the TPFBridge bSDD for classifications
   auto classification = new typename Schema::IfcClassification(
      std::string("US Bridge")/*Source*/,
      std::string("1") /*Edition*/,
      std::string("2024-11-12") /*EditionDate*/,
      std::string("US Bridge"),
      boost::none /*Description*/,
      std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1") /*Specification*/,
      boost::none /*ReferenceTokens*/);
   file.addEntity(classification);

   auto project = file.getSingle<typename Schema::IfcProject>();

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr projects(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   projects->push(project);

   auto rel_associates_classification = new typename Schema::IfcRelAssociatesClassification(
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

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ProjectCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_projects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_projects->push(project);

   auto project_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_projects, property_set);
   file.addEntity(project_properties);
}

template <typename Schema>
void Create_Pset_usBridge_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   // Was in TPFBridge, but not in USBridge
    
   //auto project = file.getSingle<typename Schema::IfcProject>();

   //typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ContractNumber"), boost::none, new typename Schema::IfcLabel(std::string("Unknown")), nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignNumber"), boost::none, new typename Schema::IfcLabel(std::string("Unknown")), nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ProjectNumber"), boost::none, new typename Schema::IfcLabel(std::string("Unknown")), nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ProjectWebsite"), boost::none, new typename Schema::IfcLabel(std::string("Unknown")), nullptr));

   //auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_ProjectCommon"), boost::none, list_of_properties);

   //typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_projects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   //related_projects->push(project);

   //auto project_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_projects, property_set);
   //file.addEntity(project_properties);
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

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_BridgeCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_bridges(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_bridges->push(bridge);

   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_Pset_usBridge_BridgeCommon(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nSpans = pBridge->GetSpanCount();
   auto nPiers = pBridge->GetPierCount();

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("usBridge_NumberOfSpans"), boost::none, new typename Schema::IfcInteger((int)nSpans), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("usBridge_NumberOfSupports"), boost::none, new typename Schema::IfcInteger((int)nPiers), nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBridge_BridgeCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_bridges(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_bridges->push(bridge);

   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_Pset_usBridge_RailingCommon(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*> railings)
{
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("usBridge_MASHCompliantRailing"), boost::none, new typename Schema::IfcBoolean(true), nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBridge_RailingCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_railings(new aggregate_of<typename Schema::IfcObjectDefinition>());
   for (auto& railing : railings)
   {
      related_railings->push(railing);
   }

   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_railings, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_Pset_usBridge_GirderCommon(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey, typename Schema::IfcElement* segment)
{
   GET_IFACE2(pBroker, IMaterials, pMaterials);
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   auto releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* displacement_unit = nullptr;

   if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);
   }

   auto fci = pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx);
   auto fc = pMaterials->GetSegmentFc28(segmentKey);
   auto fpj = pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent);
   //if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   //{
   //   fci = WBFL::Units::ConvertFromSysUnits(fci, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //   fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //   fpj = WBFL::Units::ConvertFromSysUnits(fpj, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //}

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(
      std::string("usBridge_ConcreteStrengthatTimeOfPrestressing"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/usBridge_ConcreteStrengthatTimeOfPrestressing"),
      new typename Schema::IfcInteger((int)fci), stress_unit));

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(
      std::string("usBridge_ConceteStrengthat28Days"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/usBridge_ConceteStrengthat28Days"),
      new typename Schema::IfcInteger((int)fc), stress_unit));

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(
      std::string("usBridge_JackingForce"),
      std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/usBridge_JackingForce"),
      new typename Schema::IfcReal(fpj), stress_unit));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBridge_GirderCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_segments->push(segment);

   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, property_set);
   file.addEntity(related_properties);

   if (options.include_camber)
   {
      GET_IFACE2(pBroker, IGirder, pGirder);

      typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());

      Float64 precamber = pGirder->GetPrecamber(segmentKey);
      if (!IsZero(precamber))
      {
         if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
         {
            precamber = WBFL::Units::ConvertFromSysUnits(precamber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
         }
         list_of_properties->push(new typename Schema::IfcPropertySingleValue(
            std::string("usBridge_BuiltInCamber"), 
            std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/usBridge_BuiltInCamber"), 
            new typename Schema::IfcReal(precamber), displacement_unit));
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
      //if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      //{
      //   camber = WBFL::Units::ConvertFromSysUnits(camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
      //}
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(
         std::string("usBridge_CamberatPrestressingRelease"),
         std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/usBridge_CamberatPrestressingRelease"),
         new typename Schema::IfcReal(camber + precamber), displacement_unit));

      auto lastIntervalIdx = pIntervals->GetIntervalCount() - 1;
      auto lastCompositeIntervalIdx = pIntervals->GetLastCompositeDeckInterval();
      GET_IFACE2(pBroker, ICombinedForces, pCombined);
      Float64 dc_final = pCombined->GetDeflection(lastIntervalIdx, lcDC, poiMS, bat, rtCumulative);
      Float64 dc_composite = pCombined->GetDeflection(lastCompositeIntervalIdx, lcDC, poiMS, bat, rtCumulative);
      Float64 dw_final = pCombined->GetDeflection(lastIntervalIdx, lcDW, poiMS, bat, rtCumulative);
      Float64 dw_composite = pCombined->GetDeflection(lastCompositeIntervalIdx, lcDW, poiMS, bat, rtCumulative);
      Float64 d = (dc_final - dc_composite) + (dw_final - dw_composite);
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(
         std::string("usBridge_DeflectionDuetoCompositLoads"),
         std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/usBridge_DeflectionDuetoCompositLoads"),
         new typename Schema::IfcReal(d), displacement_unit));

      // https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/TPFBridge_MemberCamber
      auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBridge_MemberCamber"), boost::none, list_of_properties);

      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
      related_segments->push(segment);

      auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, property_set);
      file.addEntity(related_properties);
   }
}

template <typename Schema>
void Create_Pset_usBridge_ReinforcementCommon(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* tendon, Float64 Pjack, bool bDebonded, Float64 ldb)
{
   // Was in TPF Bridge, but not in USBridge
   
//   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
//
//   list_of_properties->push(new typename Schema::IfcPropertySingleValue(
//      std::string("tpfBridge_TendonJackingForce"),
//      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonJackingForce"),
//      new typename Schema::IfcReal(Pjack),
//      nullptr));
//
//   list_of_properties->push(new typename Schema::IfcPropertySingleValue(
//      std::string("tpfBridge_TendonBonding"),
//      std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonBonding"),
//#pragma Reminder("bSDD - should the applicable values be the string value or the URI reference to the string value?")
//      // not sure if this should be URI reference or "Debonded" "Bonded" both are strings
//      new typename Schema::IfcURIReference(bDebonded ? "https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonBonding/value/TendonBondingDebonded" : "https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonBonding/value/TendonBondingBonded"),
//      nullptr));
//
//   if (bDebonded)
//   {
//      std::ostringstream os;
//      os << ldb;
//      list_of_properties->push(new typename Schema::IfcPropertySingleValue(
//         std::string("tpfBridge_TendonDebondedLength"),
//         std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/prop/tpfBridge_TendonDebondedLength"),
//#pragma Reminder("bSDD - why is debond length a string?")
//         new typename Schema::IfcText(os.str()), // this could be IfcIdentifier, IfcLabel, or IfcText none actually represent a value
//         nullptr));
//   }
//
//   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_ReinforcementCommon"), boost::none, list_of_properties);
//
//   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_tendons(new aggregate_of<typename Schema::IfcObjectDefinition>());
//   related_tendons->push(tendon);
//
//   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_tendons, property_set);
//   file.addEntity(related_properties);
}

template <typename Schema>
void Classify_Bridge(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
{
   auto classification = file.getSingle<typename Schema::IfcClassification>();

   auto classification_reference = new typename Schema::IfcClassificationReference(
      std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Bridge"),
      std::string("usBridge_Bridge") /*Identification*/,
      std::string("Bridge") /*Name*/,
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);
   file.addEntity(classification_reference);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_bridges(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   related_bridges->push(bridge);

   auto related_classes = new typename Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, classification_reference);

   file.addEntity(related_classes);
}

template <typename Schema>
void Classify_BridgePart(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* part, const std::string& uri, const std::string& code, const std::string& name, const std::string& part_type)
{
   std::vector<typename Schema::IfcProduct*> parts{ part };
   Classify_BridgeParts(file, parts, uri, code, name, part_type);
}

template <typename Schema>
void Classify_BridgeParts(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& parts, const std::string& uri, const std::string& code, const std::string& name, const std::string& part_type)
{
   if (parts.size() == 0)
      return;

   auto classification = file.getSingle<typename Schema::IfcClassification>();

   auto classification_reference = new typename Schema::IfcClassificationReference(
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

   auto related_classes = new typename Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_parts, classification_reference);
   file.addEntity(related_classes);
}

template <typename Schema>
void Classify_Superstructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* superstructure)
{
   Classify_BridgePart(file, superstructure, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_BridgeSuperstructure"), std::string("usBridge_BridgeSuperstructure"), std::string("IfcBridgePartSUPERSTRUCTURE"), std::string("IfcBridgePart.SUPERSTRUCTURE"));
}

template <typename Schema>
void Classify_Substructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* substructure)
{
   Classify_BridgePart(file, substructure, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_BridgeSubstructure"), std::string("usBridge_BridgeSubstructure"), std::string("IfcBridgePartSUBSTRUCTURE"), std::string("IfcBridgePart.SUBSTRUCTURE"));
}

template <typename Schema>
void Classify_Abutment(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* abutment)
{
   Classify_BridgePart(file, abutment, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_AbutmentSpatial"), std::string("usBridge_AbutmentSpatial"), std::string("IfcBridgePartABUTMENT"), std::string("IfcBridgePart.ABUTMENT"));
}

template <typename Schema>
void Classify_Pier(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* pier)
{
   Classify_BridgePart(file, pier, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_PierSpatial"), std::string("usBridge_PierSpatial"), std::string("IfcBridgePartPIER"), std::string("IfcBridgePart.PIER"));
}

template <typename Schema>
void Classify_Foundation(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* foundation)
{
   Classify_BridgePart(file, foundation, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Foundation"), std::string("usBridge_Foundation"), std::string("IfcBridgePartFOUNDATION"), std::string("IfcBridgePart.FOUNDATION"));
}

template <typename Schema>
void Classify_Railings(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& railings)
{
   Classify_BridgeParts(file, railings, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Railing"), std::string("usBridge_Railing"), std::string("IfcRailingBALUSTRADE"), std::string("IfcRailing.BALUSTRADE"));
}

template <typename Schema>
void Classify_Girders(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& girders)
{
   Classify_BridgeParts(file, girders, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Girder"), std::string("usBridge_Girder"), std::string("IfcElementAssemblyGIRDER"), std::string("IfcElementAssembly.GIRDER"));
}

template <typename Schema>
void Classify_PrecastGirderElements(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*>& girders)
{
   Classify_BridgeParts(file, girders, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_GirderPrestressedConcrete"), std::string("usBridge_GirderPrestressedConcrete"), std::string("IfcBeamGIRDER_SEGMENT"), std::string("IfcBeam.GIRDER_SEGMENT"));
}

template <typename Schema>
void Classify_Prestressing(IfcHierarchyHelper<Schema>& file, typename Schema::IfcTendonType* tendon)
{
   // Was in TPF Bridge, but not in USBridge

   //auto classification = file.getSingle<typename Schema::IfcClassification>();

   //auto classification_reference = new typename Schema::IfcClassificationReference(
   //   std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/2/class/tpfBridge_Prestressing"),
   //   std::string("tpfBridge_Prestressing") /*Identification*/,
   //   std::string("Prestressing") /*Name*/,
   //   classification,
   //   boost::none /*Description*/, boost::none /*Sort*/);
   //file.addEntity(classification_reference);

   //typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_tendons(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   //related_tendons->push(tendon);

   //auto related_classes = new typename Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_tendons, classification_reference);

   //file.addEntity(related_classes);
}
