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
void Create_Pset_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   auto project = file.getSingle<typename Schema::IfcProject>();

   typename Schema::IfcProperty::list::ptr list_of_properties(new Schema::IfcProperty::list);

   // 5.1.8.1 PEnum_ProjectType
   std::vector<std::string> enum_values{ "MODIFICAITON","NEWBUILD","OPERATIONMAINTENANCE","RENOVATION","REPAIR" };
   auto project_type_enum = createPropertyEnumeration<Schema>("PEnum_ProjectType", enum_values);
   auto project_type_property = createPropertyEnumeratedValue<Schema>("ProjectType", project_type_enum, "NEWBUILD");

   list_of_properties->push(project_type_property);

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ProjectCommon"), boost::none, list_of_properties);
   file.addEntity(property_set);

   AddPropertySet(file, project, property_set);
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
typename Schema::IfcPropertySet* Create_Pset_ConcreteElementGeneral(IfcHierarchyHelper<Schema>& file,std::string assemblyPlace,std::string castingMethod)
{
   // Pset_ConcreteElementGeneral
   typename Schema::IfcProperty::list::ptr concrete_element_general_properties(new typename Schema::IfcProperty::list);
   
   // PEnum_AssemblyPlace
   std::vector<std::string> assembly_place_enum_values{ "FACTORY","OFFSITE","SITE","OTHER","UNKNOWN","UNSET" };
   auto assembly_place_property_enum_values = createPropertyEnumeration<Schema>("PEnum_AssemblyPlace", assembly_place_enum_values);
   auto assembly_place = createPropertyEnumeratedValue<Schema>("AssemblyPlace", assembly_place_property_enum_values, "FACTORY");
   concrete_element_general_properties->push(assembly_place);
   
   // PEnum_ConcreteCastingMethod
   std::vector<std::string> casting_method_enum_values{ "INSITU","MIXED","PRECAST","PRINTED","OTHER","UNKNOWN","UNSET" };
   auto casting_method_property_enum_values = createPropertyEnumeration<Schema>("PEnum_ConcreteCastingMethod", casting_method_enum_values);
   auto casting_method = createPropertyEnumeratedValue<Schema>("CastingMethod", casting_method_property_enum_values, "PRECAST");
   concrete_element_general_properties->push(casting_method);
   
   // create Pset_ConcreteElementGeneral
   auto pset_concrete_element_general = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ConcreteElementGeneral"), boost::none, concrete_element_general_properties);
   file.addEntity(pset_concrete_element_general);
   return pset_concrete_element_general;
}

template <typename Schema>
void Create_Pset_BeamCommon(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey, typename Schema::IfcElement* segment)
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

   AddPropertySet(file, segment, property_set);
}

template <typename Schema>
void Create_Pset_PrecastConcreteElementGeneral(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey, typename Schema::IfcElement* segment)
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
      GET_IFACE2(pBroker, IBridge, pBridge);
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

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TypeDesignation"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CornerChamfer"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ManufacturingToleranceClass"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FormStrippingStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LiftingStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ReleaseStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumAllowableSupportLength"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("InitialTension"), boost::none, new typename Schema::IfcPressureMeasure(fpj), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TendonRelaxation"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TransportationStrength"), boost::none, new typename Schema::IfcPressureMeasure(fc), stress_unit));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SupportDuringTransportDescription"), boost::none, nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertyReferenceValue(std::string("SupportDuringTransportDocReference"), boost::none, boost::none, nullptr));
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

   AddPropertySet(file, segment, property_set);
}

template <typename Schema>
void Create_usBrPset_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   auto project = file.getSingle<typename Schema::IfcProject>();

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

   AddPropertySet(file, project, property_set);
}

// This Pset is not required by usBridge.
//template <typename Schema>
//void Create_Pset_BridgeCommon(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
//{
//   // 5.4.8.4 PEnum_StructureIndicator
//   std::vector<std::string> enum_values{ "COATED","COMPOSITE","HOMOGENEOUS" };
//   auto penum = createPropertyEnumeration<Schema>("PEnum_StructureIndicator", enum_values);
//   auto property = createPropertyEnumeratedValue<Schema>("StructureIndicator", penum, "COMPOSITE");
//
//   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
//   list_of_properties->push(property);
//
//   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_BridgeCommon"), boost::none, list_of_properties);
//
//   typename Schema::IfcObjectDefinition::list::ptr related_bridges(new typename Schema::IfcObjectDefinition::list);
//   related_bridges->push(bridge);
//
//   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
//   file.addEntity(related_properties);
//}

template <typename Schema>
void Create_usBrPset_Common(IfcHierarchyHelper<Schema>& file, typename Schema::IfcObject* object)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   std::vector<std::string> enum_values{ "New","Existing - Remain","Existing - Remove","Temporary", "Other"};
   auto property_enum_values = createPropertyEnumeration<Schema>("usBrPEnum_ElementStatus", enum_values); // creates an IfcPropertyEnumeration
   auto status = createPropertyEnumeratedValue<Schema>("Status", property_enum_values, enum_values.front()); // creates an IfcPropertyEnumeratedValue
   status->setSpecification(std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Status"));
   list_of_properties->push(status);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("AssociatedStandard"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/AssociatedStandard"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Note"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/Note"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SNBIElementNumber"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SNBIElementNumber"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("WorkingDrawingApproval"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/WorkingDrawingApproval"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_Common"), boost::none, list_of_properties);

   AddPropertySet(file, object, property_set);
}

template <typename Schema>
void Create_usBrPset_PayItemQuantities(IfcHierarchyHelper<Schema>& file, typename Schema::IfcObject* object)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PayItemReference"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PayItemReference"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("PayQuantity"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/PayQuantity"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("UnitOfMeasure"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/UnitOfMeasure"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_PayItemQuantities"), boost::none, list_of_properties);
   file.addEntity(property_set);

   AddPropertySet(file, object, property_set);
}

template <typename Schema>
void Create_usBrPset_BridgeGeometry(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nSpans = pBridge->GetSpanCount();

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeEndStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeEndStation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeEndStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeEndStationOffset"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeLength"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeLength"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeSkewAngle"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeSkewAngle"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeStartStation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeStartStation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BridgeStartStationOffset"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BridgeStartStationOffset"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LowBeamElevation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/LowBeamElevation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("NumberOfSpans"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/NumberOfSpans"), new typename Schema::IfcInteger((int)nSpans), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("RoadwayWidth"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/RoadwayWidth"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_BridgeGeometry"), boost::none, list_of_properties);
   file.addEntity(property_set);

   AddPropertySet(file, bridge, property_set);
}

template <typename Schema>
void Create_usBrPset_BridgeIdentification(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_DesignLoad(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_FeatureIdentification(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
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
void Create_usBrPset_PayItemQuantities(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_Railroad(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_Roadway(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge* bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_MASH(IfcHierarchyHelper<Schema>& file, std::vector<typename Schema::IfcProduct*> barriers)
{
   // placeholder - to be implemented later

   // commented out because we don't have the properties for the updated usBrPset_MASH property set.
   // MASHComplianBarrier is no longer a property in the usBridge DD
   //typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   //list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("usBridge_MASHCompliantBarrier"), boost::none, new typename Schema::IfcBoolean(true), nullptr));

   //auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_MASH"), boost::none, list_of_properties);

   //typename Schema::IfcObjectDefinition::list::ptr related_barriers(new typename Schema::IfcObjectDefinition::list);
   //for (auto& barrier : barriers)
   //{
   //   related_barriers->push(barrier);
   //}

   //auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_barriers, property_set);
   //file.addEntity(related_properties);
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
void Create_usBrPset_PrecastConcreteBeam(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey, typename Schema::IfcElement* segment)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   auto pGirder = pIBridgeDesc->GetGirder(segmentKey);
   auto shape_name = pGirder->GetGirderName();

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DeflectionLongTerm"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DeflectionLongTerm"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DeflectionShortTerm"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/DeflectionShortTerm"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ElasticShortening"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ElasticShortening"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LiftingLoopLocation"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/LiftingLoopLocation"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MidSpanCamberAfterLosses"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MidSpanCamberAfterLosses"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MidSpanCamberAfterRelease"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MidSpanCamberAtRelease"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("MinimumTimeToDeckPlacement"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/MinimumTimetoDeckPlacement"), nullptr, nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ShapeName"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/ShapeName"), new typename Schema::IfcLabel(T2A(shape_name)), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TopSurfaceFinish"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/TopSurfaceFinish"), nullptr, nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_PrecastConcreteBeam"), boost::none, list_of_properties);
   file.addEntity(property_set);

   AddPropertySet(file, segment, property_set);
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

template <typename Schema,typename Reinforcing>
void Create_usBrPset_ACIReinforcingBarType(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Reinforcing* rebar_type, std::string mark, const WBFL::Materials::Rebar* pRebar)
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

   if constexpr (std::is_same_v<Reinforcing, typename Schema::IfcReinforcingBarType>)
   {
      AddPropertySetToTypeObject(file, rebar_type, property_set);
   }
   else
   {
      AddPropertySet(file, rebar_type, property_set);
   }
}

template <typename Schema>
void Create_usBrPset_ACIBarShape(IfcHierarchyHelper<Schema>& file, typename Schema::IfcReinforcingBar* rebar)
{
}

template <typename Schema,typename Reinforcing>
void Create_usBrPset_ACIReinforcingBar(IfcHierarchyHelper<Schema>& file, typename Reinforcing* rebar_type,std::string element,std::string use,std::string position)
{
   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarElement"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarElement"), new typename Schema::IfcLabel(element.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarUse"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarUse"), new typename Schema::IfcLabel(use.c_str()), nullptr));
   list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BarPosition"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BarPosition"), new typename Schema::IfcLabel(position.c_str()), nullptr));

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ACIReinforcingBar"), boost::none, list_of_properties);
   file.addEntity(property_set);

   if constexpr (std::is_same_v<Reinforcing, typename Schema::IfcReinforcingBarType>)
   {
      AddPropertySetToTypeObject(file, rebar_type, property_set);
   }
   else
   {
      AddPropertySet(file, rebar_type, property_set);
   }
}

template <typename Schema>
void Create_usBrPset_Reinforcing(IfcHierarchyHelper<Schema>& file, typename Schema::IfcReinforcingBar* rebar)
{
}

template <typename Schema,typename Reinforcing>
void Create_usBrPset_ReinforcingCover(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Reinforcing* rebar_type, std::optional<double> top, std::optional<double> side, std::optional<double> bottom, std::optional<double> end)
{
   if (!top && !side && !end && !bottom)
      return;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit* cover_unit = nullptr;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      cover_unit = GetXSectionDimUnit<Schema>(file, pBroker);
   }

   typename Schema::IfcProperty::list::ptr list_of_properties(new typename Schema::IfcProperty::list);

   if (top)
   {
      double value = *top;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*top, pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TopFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/TopFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   if (side)
   {
      double value = *side;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*side, pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("SideFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/SideFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   if (end)
   {
      double value = *end;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*end, pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("EndFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/EndFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   if (bottom)
   {
      double value = *bottom;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*bottom, pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure);
      }
      list_of_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BottomFaceCover"), std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/prop/BottomFaceCover"), new typename Schema::IfcPositiveLengthMeasure(value), cover_unit));
   }

   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("usBrPset_ACIReinforcingCover"), boost::none, list_of_properties);
   file.addEntity(property_set);

   if constexpr (std::is_same_v<Reinforcing, typename Schema::IfcReinforcingBarType>)
   {
      AddPropertySetToTypeObject(file, rebar_type, property_set);
   }
   else
   {
      AddPropertySet(file, rebar_type, property_set);
   }

}