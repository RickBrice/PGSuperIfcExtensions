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
#pragma once

#include <IFace/Tools.h>
#include <IFace/PointOfInterest.h>
#include <IFace\Bridge.h>
#include <IFace\AnalysisResults.h>
#include <IFace\Intervals.h>
#include "Units.h"
#include "IfcExporter.h"
#include "Properties.h"

#include <PsgLib\BridgeDescription2.h>


template <typename Schema>
typename Schema::IfcPropertySet* Create_Pset_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new Schema::IfcProperty::list);

   // 5.1.8.1 PEnum_ProjectType
   std::vector<std::string> enum_values{ "MODIFICAITON","NEWBUILD","OPERATIONMAINTENANCE","RENOVATION","REPAIR" };
   auto project_type_enum = createPropertyEnumeration<Schema>("PEnum_ProjectType", enum_values);
   auto project_type_property = createPropertyEnumeratedValue<Schema>("ProjectType", project_type_enum, "NEWBUILD");

   list_of_properties->push(project_type_property);

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ProjectCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}


template <typename Schema>
void Create_Pset_MaterialSteel_Strand(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial* material, const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;

   // Pset_MaterialSteel
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;

   // define strand properties
   auto fy = pStrand->GetYieldStrength();
   auto fpu = pStrand->GetUltimateStrength();
   auto eu = 0.035; // from ASTM A416 for prestressing strands

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      fy = WBFL::Units::ConvertFromSysUnits(fy, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fpu = WBFL::Units::ConvertFromSysUnits(fpu, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   // Need to clean this up
   // ASTM A416 is for low relaxation strand... PGSuper does low relaxation and stress relieved
   // ASTM A416 is for Grade 250 and Grade 270... PGSuper does grade 300 as well, but there doesn't seem to be an ASTM
   // We are assuming same material for all strands, but that is not the case in the PGSuper data model
   // straight, harped, and temporary can be different - Grade 250, Grade 270, Grade 300
   // Strand size/diameter is a property on IfcTendon
   std::ostringstream os;
   os << "ASTM A416 Grade " << T2A(WBFL::Materials::PsStrand::GetGrade(pStrand->GetGrade(), true/*US units*/).c_str());
   auto grade = os.str();

   typename Schema::IfcProperty::list::ptr material_steel_properties(new typename Schema::IfcProperty::list);
   //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("YieldStress"), boost::none, new typename Schema::IfcPressureMeasure(fy), stress_unit));
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("UltimateStress"), boost::none, new typename Schema::IfcPressureMeasure(fpu), stress_unit));
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("UltimateStrain"), boost::none, new typename Schema::IfcPositiveRatioMeasure(eu), nullptr));
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StructuralGrade"), boost::none, new typename Schema::IfcLabel(grade.c_str()), nullptr));
   auto pset_material_steel = new typename Schema::IfcMaterialProperties(std::string("Pset_MaterialSteel"), boost::none/*description*/, material_steel_properties, material);
   file.addEntity(pset_material_steel);
}

template <typename Schema>
void Create_Pset_MaterialSteel_ReinforcingBar(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial* material, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   // Pset_MaterialSteel
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;

   // define rebar properties
   auto fy = pRebar->GetYieldStrength();
   auto fpu = pRebar->GetUltimateStrength();
   auto eu = pRebar->GetElongation();

   if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      fy = WBFL::Units::ConvertFromSysUnits(fy, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fpu = WBFL::Units::ConvertFromSysUnits(fpu, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   std::ostringstream os;
   os << T2A(WBFL::LRFD::RebarPool::GetMaterialName(pRebar->GetType(),pRebar->GetGrade()).c_str());
   auto grade = os.str();

   typename Schema::IfcProperty::list::ptr material_steel_properties(new typename Schema::IfcProperty::list);
   //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("YieldStress"), boost::none, new typename Schema::IfcPressureMeasure(fy), stress_unit));
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("UltimateStress"), boost::none, new typename Schema::IfcPressureMeasure(fpu), stress_unit));
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("UltimateStrain"), boost::none, new typename Schema::IfcPositiveRatioMeasure(eu), nullptr));
   material_steel_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StructuralGrade"), boost::none, new typename Schema::IfcLabel(grade.c_str()), nullptr));
   auto material_properties = new typename Schema::IfcMaterialProperties(std::string("Pset_MaterialSteel"), boost::none/*description*/, material_steel_properties, material);
   file.addEntity(material_properties);
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_Pset_ConcreteElementGeneral(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, std::optional<std::string> assemblyPlace,std::optional<std::string> castingMethod,std::optional<Float64> fc)
{
   USES_CONVERSION;

   // Pset_ConcreteElementGeneral
   typename Schema::IfcProperty::list::ptr concrete_element_general_properties(new typename Schema::IfcProperty::list);
   
   // PEnum_AssemblyPlace
   if (assemblyPlace.has_value())
   {
      std::vector<std::string> assembly_place_enum_values{ "FACTORY","OFFSITE","SITE","OTHER","UNKNOWN","UNSET" };
      auto assembly_place_property_enum_values = createPropertyEnumeration<Schema>("PEnum_AssemblyPlace", assembly_place_enum_values);
      auto assembly_place = createPropertyEnumeratedValue<Schema>("AssemblyPlace", assembly_place_property_enum_values, *assemblyPlace);
      concrete_element_general_properties->push(assembly_place);
   }
   
   // PEnum_ConcreteCastingMethod
   if (castingMethod.has_value())
   {
      std::vector<std::string> casting_method_enum_values{ "INSITU","MIXED","PRECAST","PRINTED","OTHER","UNKNOWN","UNSET" };
      auto casting_method_property_enum_values = createPropertyEnumeration<Schema>("PEnum_ConcreteCastingMethod", casting_method_enum_values);
      auto casting_method = createPropertyEnumeratedValue<Schema>("CastingMethod", casting_method_property_enum_values, *castingMethod);
      concrete_element_general_properties->push(casting_method);
   }

   if (fc.has_value())
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
      auto fc_value = *fc;

      std::ostringstream os;
      os << T2A(::FormatDimension(fc_value, pDisplayUnits->GetStressUnit())) << std::endl;
      concrete_element_general_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StrengthClass"), boost::none, new typename Schema::IfcLabel(os.str()),nullptr));
   }
   
   // create Pset_ConcreteElementGeneral
   auto pset_concrete_element_general = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ConcreteElementGeneral"), boost::none, concrete_element_general_properties);
   file.addEntity(pset_concrete_element_general);
   return pset_concrete_element_general;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_Pset_BeamCommon(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
{
   GET_IFACE2(pBroker, IBridge, pBridge);

   auto span_length = pBridge->GetSegmentSpanLength(segmentKey);
   auto slope = pBridge->GetSegmentSlope(segmentKey);

   auto slope_angle = atan(slope);

   GET_IFACE2(pBroker, IGirder, pGirder);
   auto roll = pGirder->GetOrientation(segmentKey);
   auto roll_angle = atan(roll);


   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);
   typename Schema::IfcConversionBasedUnit* length_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* angle_unit = nullptr;

   if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);
      angle_unit = GetAngleUnit<Schema>(file, pBroker);

      span_length = WBFL::Units::ConvertFromSysUnits(span_length, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      slope_angle = WBFL::Units::ConvertFromSysUnits(slope_angle, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
      roll_angle = WBFL::Units::ConvertFromSysUnits(roll_angle, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
   }


   typename Schema::IfcProperty::list::ptr list_of_properties(new Schema::IfcProperty::list);

   // Depreciated in IFC4.3. Use the Name attribute of the relating type
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Reference"), boost::none, nullptr, nullptr)); 

   std::vector<std::string> enum_values{ "DEMOLISH","EXISTING","NEW","TEMPORARY", "OTHER","UNKNOWN","UNSET"};
   auto property_enum_values = createPropertyEnumeration<Schema>("PEnum_ElementStatus", enum_values);
   auto status = createPropertyEnumeratedValue<Schema>("Status", property_enum_values, "NEW");
   list_of_properties->push(status);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Span"), boost::none, new typename Schema::IfcPositiveLengthMeasure(span_length), length_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Slope"), boost::none, new typename Schema::IfcPlaneAngleMeasure(slope_angle), angle_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Roll"), boost::none, new typename Schema::IfcPlaneAngleMeasure(roll_angle), angle_unit));

   // This properties are not applicable to bridges
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("IsExternal"), boost::none, nullptr, nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ThermalTransmittance"), boost::none, nullptr, nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LoadBearing"), boost::none, nullptr, nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FireRating"), boost::none, nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_BeamCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_Pset_PrecastConcreteElementGeneral(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IBridge, pBridge);
   GET_IFACE2(pBroker, IMaterials, pMaterials);
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   auto releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
   auto liftingIntervalIdx = pIntervals->GetLiftSegmentInterval(segmentKey);
   auto haulingIntervalIdx = pIntervals->GetHaulSegmentInterval(segmentKey);

   auto fci = pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx);
   auto fcl = pMaterials->GetSegmentFc(segmentKey, liftingIntervalIdx);
   auto fch = pMaterials->GetSegmentFc(segmentKey, haulingIntervalIdx);
   auto fc = pMaterials->GetSegmentFc28(segmentKey);
   auto fpj = pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent);


   Float64 batter = 0.0;
   if (options.batter_ends)
   {
      Float64 slope = pBridge->GetSegmentSlope(segmentKey);
      batter = atan(slope);
   }

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* displacement_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* angle_unit = nullptr;

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file,pBroker);
      displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);
      angle_unit = GetAngleUnit<Schema>(file, pBroker);

      fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fci = WBFL::Units::ConvertFromSysUnits(fci, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fcl = WBFL::Units::ConvertFromSysUnits(fcl, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fch = WBFL::Units::ConvertFromSysUnits(fch, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fpj = WBFL::Units::ConvertFromSysUnits(fpj, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      batter = WBFL::Units::ConvertFromSysUnits(batter, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
   }

   Float64 camber = 0.;
   if (options.include_camber)
   {
      GET_IFACE2(pBroker, IGirder, pGirder);
      Float64 precamber = pGirder->GetPrecamber(segmentKey);

      GET_IFACE2(pBroker, IPointOfInterest, pPoi);
      PoiList vPoi;
      pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
      CHECK(vPoi.size() == 1);
      const pgsPointOfInterest& poiMS = vPoi.front();

      GET_IFACE2(pBroker, IProductForces, pProduct);
      auto bat = pProduct->GetBridgeAnalysisType(pgsTypes::Minimize); // minimize because we want the greatest downward deflection

      Float64 ps = pProduct->GetDeflection(releaseIntervalIdx, pgsTypes::pftPretension, poiMS, bat, rtCumulative, false);
      Float64 girder = pProduct->GetDeflection(releaseIntervalIdx, pgsTypes::pftGirder, poiMS, bat, rtCumulative, false);
      camber = ps + girder + precamber;

      Float64 Lg = pBridge->GetSegmentPlanLength(segmentKey);
      camber /= Lg; // convert camber to a ratio as per the CamberAtMidSpan property definition
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new Schema::IfcProperty::list);

   GET_IFACE2(pBroker, IBridgeDescription, pBridgeDesc);
   auto family_name = pBridgeDesc->GetBridgeDescription()->GetGirderFamilyName();

   auto hauling_data = pBridgeDesc->GetBridgeDescription()->GetGirderGroup(segmentKey.groupIndex)->GetGirder(segmentKey.girderIndex)->GetSegment(segmentKey.segmentIndex)->HandlingData;
   auto bunk_point = std::max(hauling_data.LeadingSupportPoint, hauling_data.TrailingSupportPoint);

   typename Schema::IfcConversionBasedUnit* length_unit = nullptr;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);
      bunk_point = WBFL::Units::ConvertFromSysUnits(bunk_point, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
   }

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TypeDesignation"), boost::none, new typename Schema::IfcLabel(T2A(family_name)), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CornerChamfer"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ManufacturingToleranceClass"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FormStrippingStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LiftingStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ReleaseStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumAllowableSupportLength"), boost::none, new typename Schema::IfcPositiveLengthMeasure(bunk_point), length_unit));

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("InitialTension"), boost::none, new typename Schema::IfcPressureMeasure(fpj), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TendonRelaxation"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TransportationStrength"), boost::none, new typename Schema::IfcPressureMeasure(fc), stress_unit));
   
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SupportDuringTransportDescription"), 
      boost::none, new typename Schema::IfcText("Assumed to be truck transported with bunking locations per MinimumAllowableSupportLength property"), nullptr));
   
   // this is the proper way to define this property, excpt that IfcPropertyReferenceValue isn't part of the AbV.
   list_of_properties->push(new typename Schema::IfcPropertyReferenceValue(std::string("SupportDuringTransportDocReference"), 
      boost::none, boost::none, 
      new typename Schema::IfcDocumentReference(
         std::string("https://www.pci.org/ItemDetail?iProductCode=CB-02-26H&Category=TRANSPORT&WebsiteKey=5a7b2064-98c2-4c8e-9b4b-18c80973da1e"), // Location
         boost::none, // Identification
         std::string("Recommended Practice for Lateral Stability of Precast, Prestressed Concrete Bridge Girders, 2nd Edition (CB-02-26H)"), // Name
         boost::none, // Description
         nullptr))); // Referenced Document

   // this is how the usBridge DD says to define this property, but it isn't valid IFC and results in validation service errors
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SupportDuringTransportDocReference"),
   //   boost::none, new typename Schema::IfcText(
   //      std::string("Recommended Practice for Lateral Stability of Precast, Prestressed Concrete Bridge Girders, 2nd Edition (CB-02-26H)")),
   //      nullptr));

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("HollowCorePlugging"), boost::none, nullptr, nullptr));
   
   if(options.include_camber)
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CamberAtMidspan"), boost::none, new typename Schema::IfcRatioMeasure(camber), nullptr));
   else
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CamberAtMidspan"), boost::none, nullptr, nullptr));
   
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BatterAtStart"), boost::none, new typename Schema::IfcPlaneAngleMeasure(batter), angle_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BatterAtEnd"), boost::none, new typename Schema::IfcPlaneAngleMeasure(batter), angle_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Twisting"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Shortening"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PieceMark"), boost::none, nullptr, nullptr));

   {
#pragma Reminder("WORKING HERE - need to use the usBridge schema for this property")
      pgsAutoGirderLabel autoLabel;
      pgsGirderLabel::UseAlphaLabel(false);
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignLocationNumber"), boost::none, new typename Schema::IfcLabel(T2A(SEGMENT_LABEL(segmentKey))), nullptr));
   }

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_PrecastConcreteElementGeneral"), boost::none, list_of_properties);
   file.addEntity(property_set);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   // we don't have these properties, but will set up the property set with blank properties so that it shows up in the file and can be filled in by hand or by a future version of the exporter.
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ApprovalStatus"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ApprovalStatus"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FileNumber"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/FileNumber"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LettingDate"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/LettingDate"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ModelPreparationDate"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ModelPreparationDate"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ModelVersion"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ModelVersion"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ProjectDirectory"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ProjectDirectory"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ProjectIdentification"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ProjectIdentification"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ProjectNumber"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ProjectNumber"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ProjectURL"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ProjectURL"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PSEData"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PSEData"), nullptr, nullptr));
   
   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ProjectCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_ProjectLocation(IfcHierarchyHelper<Schema>& file)
{
   // we don't have these properties, but will set up the property set with blank properties so that it shows up in the file and can be filled in by hand or by a future version of the exporter.
   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
#pragma Reminder("WORKING HERE - need to update the URLs - bSDD is down right now")
   // County and State are required properties.
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("City"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/City"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("County"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/County"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("District"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/District"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("State"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/State"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Section"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Section"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Township"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Township"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Range"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Range"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ProjectLocation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_ProjectLocation"), list_of_properties);
   file.addEntity(property_set);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_Common(IfcHierarchyHelper<Schema>& file)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   std::vector<std::string> enum_values{ "NEW","EXISTING - REMAIN","EXISTING - REMOVE","TEMPORARY", "OTHER"};
   auto property_enum_values = createPropertyEnumeration<Schema>("usBrPEnum_ElementStatus", enum_values); // creates an IfcPropertyEnumeration
   auto status = createPropertyEnumeratedValue<Schema>("Status", property_enum_values, enum_values.front()); // creates an IfcPropertyEnumeratedValue
   status->setSpecification(std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Status"));
   list_of_properties->push(status);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("AssociatedStandard"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/AssociatedStandard"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Note"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Note"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SNBIElementNumber"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SNBIElementNumber"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("WorkingDrawingApproval"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/WorkingDrawingApproval"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_Common"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_BridgePartCommon(IfcHierarchyHelper<Schema>& file)
{
   // this property set is required by the IDS, but none of its properties are required.

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndSkew"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndSkew"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndStation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndStationOffset"), nullptr, nullptr));

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartSkew"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartSkew"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartStation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartStationOffset"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_BridgePartCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_SlabCommon(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatingNote"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatingNote"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CrossSectionalArea"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CrossSectionalArea"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CrossSlope"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CrossSlope"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndSkew"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndSkew"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndStation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndStationOffset"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartSkew"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartSkew"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartStation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartStationOffset"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_SlabCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_SubstructureCommon(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, PierIndexType pierIdx = INVALID_INDEX)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   if (pierIdx != INVALID_INDEX)
   {
      GET_IFACE2(pBroker, IBridge, pBridge);
      auto pier_station = pBridge->GetPierStation(pierIdx);

      CComPtr<IAngle> angle;
      pBridge->GetPierSkew(pierIdx,&angle);
      Float64 skew;
      angle->get_Value(&skew);

      GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);
      typename Schema::IfcConversionBasedUnit* station_unit = nullptr;
      typename Schema::IfcConversionBasedUnit* angle_unit = nullptr;
      if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         station_unit = GetSpanLengthUnit<Schema>(file, pBroker);
         pier_station = WBFL::Units::ConvertFromSysUnits(pier_station, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);

         angle_unit = GetAngleUnit<Schema>(file, pBroker);
         skew = WBFL::Units::ConvertFromSysUnits(skew, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
         skew = RoundOff(skew, 0.0001);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StationAheadBearing"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StationAheadBearing"), nullptr, nullptr));
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StationAtCenterline"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PierStation"), new typename Schema::IfcLengthMeasure(pier_station), station_unit));
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StationBackBearing"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StationBackBearing"), nullptr, nullptr));
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StationOffsetAheadBearing"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StationOffsetAheadBearing"), nullptr, nullptr));
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StationOffsetAtCenterline"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StationOffsetAtCenterline"), nullptr, nullptr));
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StationOffsetBackBearing"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StationOffsetBackBearing"), nullptr, nullptr));
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SubstructureSkewAngle"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SubstructureSkewAngle"), new typename Schema::IfcPlaneAngleMeasure(skew), angle_unit));
   }

   std::vector<std::string> enum_values{ "New","Other" }; // see SNBI
   auto property_enum_values = createPropertyEnumeration<Schema>("usBrPEnum_SubstructureType", enum_values); // creates an IfcPropertyEnumeration
   auto substruture_type = createPropertyEnumeratedValue<Schema>("SubstructureType", property_enum_values, enum_values.front()); // creates an IfcPropertyEnumeratedValue
   substruture_type->setSpecification(std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SubstructureType"));
   list_of_properties->push(substruture_type);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ReturnInterval"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ReturnInterval"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ScourElevation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ScourElevation"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_SubstructureCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_PayItemQuantities(IfcHierarchyHelper<Schema>& file)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PayItemReference"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PayItemReference"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PayQuantity"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PayQuantity"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("UnitOfMeasure"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/UnitOfMeasure"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_PayItemQuantities"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_BridgeGeometry(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nSpans = pBridge->GetSpanCount();

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   auto length = pBridge->GetLength();
   auto width = pBridge->GetCurbToCurbWidth(0.0); // getting width a start of bridge. the actual requirements are the most restrictive width. Report to nearest 0.1 ft
   auto start_station = pBridge->GetPierStation(0);
   auto end_station = pBridge->GetPierStation(nSpans);

   double skew = 0;
   for (PierIndexType i = 0; i < nSpans; i++)
   {
      CComPtr<IAngle> angle;
      pBridge->GetPierSkew(i,&angle);
      Float64 pier_skew;
      angle->get_Value(&pier_skew);
      if (std::fabs(skew) < std::fabs(pier_skew))
         skew = pier_skew;
   }

   typename Schema::IfcConversionBasedUnit* length_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* angle_unit = nullptr;
   if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      length = WBFL::Units::ConvertFromSysUnits(length, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      width = WBFL::Units::ConvertFromSysUnits(width, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      width = RoundOff(width, 0.1); // round to nearest 0.1 ft

      start_station = WBFL::Units::ConvertFromSysUnits(start_station, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      end_station = WBFL::Units::ConvertFromSysUnits(end_station, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);

      angle_unit = GetAngleUnit<Schema>(file, pBroker);
      skew = WBFL::Units::ConvertFromSysUnits(skew, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
      skew = RoundOff(skew, 0.0001); // round to nearest 0.0001 degree
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeEndStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeEndStation"), new typename Schema::IfcLengthMeasure(end_station), length_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeEndStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeEndStationOffset"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeLength"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeLength"), new typename Schema::IfcPositiveLengthMeasure(length), length_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeSkew"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeSkew"), new typename Schema::IfcPlaneAngleMeasure(skew), angle_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeStartStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeStartStation"), new typename Schema::IfcLengthMeasure(start_station), length_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeStartStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeStartStationOffset"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LowBeamElevation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/LowBeamElevation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("NumberOfSpans"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/NumberOfSpans"), new typename Schema::IfcInteger((int)nSpans), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("RoadwayWidth"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/RoadwayWidth"), new typename Schema::IfcPositiveLengthMeasure(width), length_unit));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_BridgeGeometry"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_BridgeIdentification(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SNBIBridgeNumber"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SNBIBridgeNumber"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SNBIBridgeName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SNBIBridgeName"), nullptr, nullptr));

   // other optional attributes to be provided later

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_BridgeIdentification"), boost::none, list_of_properties);
   file.addEntity(property_set);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_DesignLoading(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, ILiveLoads, pLiveLoads);

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignMethodology"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DesignMethodology"), new typename Schema::IfcLabel("LRFD"), nullptr));

   if (pLiveLoads->IsLiveLoadDefined(pgsTypes::lltDesign))
   {
      auto live_load_names = pLiveLoads->GetLiveLoadNames(pgsTypes::LiveLoadType::lltDesign);
      auto live_load_name = live_load_names[0];
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("VehicularLiveLoad"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/VehicularLiveLoad"), new typename Schema::IfcLabel(T2A(live_load_name.c_str())), nullptr));
   }

#pragma Reminder("WORKING HERE - finish filling out these properties")
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PedestrianLiveLoad"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PedestrianLiveLoad"), nullptr, nullptr));
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignFutureWearingLoad"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DesignFutureWearingLoad"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_DesignLoading"), boost::none, list_of_properties);
   file.addEntity(property_set);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_FeatureIdentification(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FeatureType"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/FeatureType"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FeatureLocation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/FeatureLocation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FeatureName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/FeatureName"), nullptr, nullptr));

   // other optional attributes to be provided later

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_FeatureIdentification"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
void Create_usBrPset_HydraulicData(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_NavigableWaterway(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_Railroad(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_Roadway(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("AnnualAverageDailyTraffic"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/AnnualAverageDailyTraffic"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("AnnualAverageDailyTruckTraffic"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/AnnualAverageDailyTruckTraffic"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignAnnualAverageDailyTraffic"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DesignAnnualAverageDailyTraffic"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignAnnualAverageDailyTruckTraffic"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DesignAnnualAverageDailyTruckTraffic"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignTrafficYear"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DesignTrafficYear"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumHorizontalClearanceLeft"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MinimumHorizontalClearanceLeft"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumHorizontalClearanceRight"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MinimumHorizontalClearanceRight"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("RoadwayMinimumVerticalClearance"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/RoadwayMinimumVerticalClearance"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("RoadwayName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/RoadwayName"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("RoadwayType"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/RoadwayType"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_Roadway"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_Roadway"), list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_RoadwaySlab(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FormType"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/FormType"), new typename Schema::IfcLabel("Wood"), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SacrificialWearingSurface"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SacrificialWearingThickness"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SurfaceArea"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SurfaceArea"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SurfaceFinish"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SurfaceFinish"), new typename Schema::IfcLabel("Raked"), nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_RoadwaySlab"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_RoadwaySlab"), list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_MASH(IfcHierarchyHelper<Schema>& file)
{
#pragma Reminder("WORKING HERE - need to get real barrier properties and names for PGSuper - dummy values used for now")
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MASHTestingLevel"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MASHTestingLevel"), new typename Schema::IfcLabel("Unknown"), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("GuardrailType"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/GuardrailType"), new typename Schema::IfcLabel("Unknown"), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarrierType"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarrierType"), new typename Schema::IfcLabel("Unknown"), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarrierHeight"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarrierHeight"), new typename Schema::IfcPositiveLengthMeasure(0.01), nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_MASH"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_MASH"), list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_TendonDebondingAtEnds(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, double debond_start,double debond_end)
{
   // Debonding at ends of beam (normal debonding)

   typename Schema::IfcConversionBasedUnit* length_unit = nullptr;

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      debond_start = WBFL::Units::ConvertFromSysUnits(debond_start, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      debond_end = WBFL::Units::ConvertFromSysUnits(debond_end, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new Schema::IfcProperty::list);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DebondLengthStart"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DebondLengthStart"), new typename Schema::IfcLengthMeasure(debond_start), length_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DebondLengthEnd"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DebondLengthEnd"), new typename Schema::IfcLengthMeasure(debond_end), length_unit));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_TendonDebonding"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_TendonDebonding"), list_of_properties);
   file.addEntity(property_set);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_TendonDebondingInCenter(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, double debond_start, double debond_end)
{
   // Debonding in the middle of the beam (typically for temporary top strands)
   typename Schema::IfcConversionBasedUnit* length_unit = nullptr;

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      debond_start = WBFL::Units::ConvertFromSysUnits(debond_start, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      debond_end = WBFL::Units::ConvertFromSysUnits(debond_end, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new Schema::IfcProperty::list);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DebondLengthMidspanToStart"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DebondLengthMidspanToStart"), new typename Schema::IfcLengthMeasure(debond_start), length_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DebondLengthMidspanToEnd"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DebondLengthMidspanToEnd"), new typename Schema::IfcLengthMeasure(debond_end), length_unit));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_TendonDebonding"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_TendonDebonding"), list_of_properties);
   file.addEntity(property_set);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_PrecastConcreteBeam(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   auto pGirder = pIBridgeDesc->GetGirder(segmentKey);
   auto shape_name = pGirder->GetGirderName();

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);

   const pgsPointOfInterest& poi = vPoi.front();
   
   GET_IFACE2(pBroker, ICamber, pCamber);
   auto initial_camber = pCamber->GetInitialCamber(poi);
   auto final_camber = pCamber->GetExcessCamber(poi, pgsTypes::CreepTime::Max);
   auto screed_camber = pCamber->GetScreedCamber(poi, pgsTypes::CreepTime::Max);

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);
   typename Schema::IfcConversionBasedUnit* deflection_unit = nullptr;
   if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      deflection_unit = GetDisplacementUnit<Schema>(file, pBroker);
      initial_camber = WBFL::Units::ConvertFromSysUnits(initial_camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
      final_camber = WBFL::Units::ConvertFromSysUnits(final_camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
      screed_camber = WBFL::Units::ConvertFromSysUnits(screed_camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DeflectionLongTerm"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DeflectionLongTerm"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DeflectionShortTerm"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DeflectionShortTerm"), new typename Schema::IfcLengthMeasure(screed_camber), deflection_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ElasticShortening"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ElasticShortening"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LiftingLoopLocation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/LiftingLoopLocation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MidSpanCamberAfterLosses"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MidSpanCamberAfterLosses"), new typename Schema::IfcLengthMeasure(final_camber), deflection_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MidSpanCamberAtRelease"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MidSpanCamberAtRelease"), new typename Schema::IfcLengthMeasure(initial_camber), deflection_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumTimeToDeckPlacement"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MinimumTimetoDeckPlacement"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ShapeName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ShapeName"), new typename Schema::IfcLabel(T2A(shape_name)), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TopSurfaceFinish"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/TopSurfaceFinish"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_PrecastConcreteBeam"), boost::none, list_of_properties);
   file.addEntity(property_set);

   return property_set;
}


inline std::string GetRebarSpecification(const WBFL::Materials::Rebar* pRebar)
{
   std::string spec;
   switch (pRebar->GetType())
   {
   case WBFL::Materials::Rebar::Type::A615: spec = "ASTM A615 (AASHTO M31)"; break;
   case WBFL::Materials::Rebar::Type::A706: spec = "ASTM A706"; break;
   case WBFL::Materials::Rebar::Type::A1035: spec = "ASTM A1035"; break;
   default: spec = "Unknown"; break;
   }
   return spec;
}

inline std::string GetRebarSpecificationEdition(const WBFL::Materials::Rebar* pRebar)
{
   std::string edition;
   switch (pRebar->GetType())
   {
   case WBFL::Materials::Rebar::Type::A615: edition = "2026"; break;
   case WBFL::Materials::Rebar::Type::A706: edition = "2026"; break;
   case WBFL::Materials::Rebar::Type::A1035: edition = "2024"; break;
   default: edition = "Unknown"; break;
   }
   return edition;
}

inline std::string GetStrandSpecification(const WBFL::Materials::PsStrand* pStrand)
{
   std::string spec("ASTM A416 (AASHTO M203)");
   return spec;
}

inline std::string GetStrandSpecificationEdition(const WBFL::Materials::PsStrand* pStrand)
{
   std::string edition("Unknown");
   return edition;
}

template <typename Schema>
void Create_usBrPset_ACITendonMaterial(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial* material, const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;

   auto spec = GetStrandSpecification(pStrand);
   auto spec_edition = GetStrandSpecificationEdition(pStrand);

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   auto fpu = pStrand->GetUltimateStrength();

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      fpu = WBFL::Units::ConvertFromSysUnits(fpu, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Specification"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Specification"), new typename Schema::IfcLabel(spec.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SpecificationVersion"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SpecificationVersion"), new typename Schema::IfcLabel(spec_edition.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TendonGrade"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/TendonGrade"), new typename Schema::IfcPressureMeasure(fpu), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatingSpecification"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatingSpecification"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatingSpecificationVersion"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatingSpecificationVersion"), nullptr, nullptr));

   auto material_properties = new typename Schema::IfcMaterialProperties(std::string("usBrPset_ACITendonMaterial"), boost::none/*description*/, list_of_properties, material);
   file.addEntity(material_properties);
}

template <typename Schema>
void Create_usBrPset_ACIReinforcingMaterial(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial* material, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   auto spec = GetRebarSpecification(pRebar);
   auto spec_edition = GetRebarSpecificationEdition(pRebar);

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   auto fy = pRebar->GetYieldStrength();

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      fy = WBFL::Units::ConvertFromSysUnits(fy, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Specification"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Specification"), new typename Schema::IfcLabel(spec.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SpecificationVersion"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SpecificationVersion"), new typename Schema::IfcLabel(spec_edition.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ReinforcingGrade"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ReinforcingGrade"), new typename Schema::IfcPressureMeasure(fy), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Subtype"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Subtype"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatingSpecification"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatingSpecification"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatingSpecificationVersion"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatingSpecificationVersion"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatingSubtype"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatingSubtype"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CoatedBeforeFabrication"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/CoatedBeforeFabrication"), nullptr, nullptr));

   auto material_properties = new typename Schema::IfcMaterialProperties(std::string("usBrPset_ACIReinforcingMaterial"), boost::none/*description*/, list_of_properties, material);
   file.addEntity(material_properties);
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_ACIReinforcingBarType(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::string mark, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   auto size = WBFL::LRFD::RebarPool::GetBarSize(pRebar->GetSize());
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarMark"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarMark"), new typename Schema::IfcLabel(mark.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarMass"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarMass"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarSize"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarSize"), new typename Schema::IfcLabel(T2A(size.c_str())), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndEndPrep"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndEndPrep"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StartEndPrep"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StartEndPrep"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ACIReinforcingBarType"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
void addBarDimension(std::string name,std::string uri,double dim,typename Schema::IfcConversionBasedUnit* unit, typename Schema::IfcProperty::list::ptr list_of_properties)
{
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(name, uri, new typename Schema::IfcPositiveLengthMeasure(dim), unit));
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_ACIBarShape(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::string bend_shape_name,
   double bend_radius, const std::unordered_map<std::string, double>& dimensions)
{
   std::string standard = "ACI 315-99";
   std::string standard_version = "1999";

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StandardName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StandardName"), new typename Schema::IfcLabel(standard.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("StandardVersion"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/StandardVersion"), new typename Schema::IfcLabel(standard_version.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BendShapeName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BendShapeName"), new typename Schema::IfcLabel(bend_shape_name.c_str()), nullptr));
   
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* dimension_unit = nullptr;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      dimension_unit = GetComponentDimUnit<Schema>(file, pBroker);
      bend_radius = WBFL::Units::ConvertFromSysUnits(bend_radius, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
   }
   addBarDimension<Schema>("DefaultInsideBendRadius", "https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DefaultInsideBendRadius", bend_radius, dimension_unit, list_of_properties);

   for (auto & [name, dim] : dimensions)
   {
      auto value = dim;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
      }
      addBarDimension<Schema>(name, "https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/" + name, value, dimension_unit, list_of_properties);
   }

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ACIBarShape"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_ACIReinforcingBar(IfcHierarchyHelper<Schema>& file, std::string element,std::string use,std::string position)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarElement"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarElement"), new typename Schema::IfcLabel(element.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarUse"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarUse"), new typename Schema::IfcLabel(use.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarPosition"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarPosition"), new typename Schema::IfcLabel(position.c_str()), nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ACIReinforcingBar"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_ACIReinforcingBar"), list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_Reinforcing(IfcHierarchyHelper<Schema>& file)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EmbedmentAtEnd"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EmbedmentAtEnd"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EmbedmentAtStart"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EmbedmentAtStart"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumSpliceLength"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MinimumSpliceLength"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PlacementMethod"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PlacementMethod"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_Reinforcing"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBrPset_Reinforcing"), list_of_properties);
   file.addEntity(property_set);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet* Create_usBrPset_ReinforcingCover(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::optional<double> top, std::optional<double> side, std::optional<double> bottom, std::optional<double> end)
{
   if (!top && !side && !end && !bottom)
      return nullptr;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* cover_unit = nullptr;
   auto length_unit = pDisplayUnits->GetComponentDimUnit();
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      cover_unit = GetComponentDimUnit<Schema>(file, pBroker);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   if (top)
   {
      double value = *top;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*top, length_unit.UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TopFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/TopFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   if (side)
   {
      double value = *side;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*side, length_unit.UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SideFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SideFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   if (end)
   {
      double value = *end;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*end, length_unit.UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   if (bottom)
   {
      double value = *bottom;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*bottom, length_unit.UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BottomFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BottomFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ACIReinforcingCover"), boost::none, list_of_properties);
   file.addEntity(property_set);
   return property_set;
}