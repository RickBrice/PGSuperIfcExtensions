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
#include "bSDD.h"
#include <PsgLib\BridgeDescription2.h>


template <typename Schema>
typename Schema::IfcPropertySet Create_Pset_ProjectCommon(hierarchy_helper<Schema>& file)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;

   // 5.1.8.1 PEnum_ProjectType
   std::vector<std::string> enum_values{ "MODIFICAITON","NEWBUILD","OPERATIONMAINTENANCE","RENOVATION","REPAIR" };
   auto project_type_enum = createPropertyEnumeration<Schema>(file, "PEnum_ProjectType", enum_values);
   auto project_type_property = createPropertyEnumeratedValue<Schema>(file, "ProjectType", project_type_enum, "NEWBUILD");

   list_of_properties.push_back(project_type_property);

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_ProjectCommon"), std::nullopt, list_of_properties);
   return property_set;
}


template <typename Schema>
void Create_Pset_MaterialSteel_Strand(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial material, const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;

   // Pset_MaterialSteel
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit stress_unit;

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

   std::vector<typename Schema::IfcProperty> material_steel_properties;
   //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("YieldStress"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fy), stress_unit));
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("UltimateStress"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fpu), stress_unit));
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("UltimateStrain"), std::nullopt, file.create<typename Schema::IfcPositiveRatioMeasure>().initialize(eu), typename Schema::IfcUnit{}));
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StructuralGrade"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(grade.c_str()), typename Schema::IfcUnit{}));
   auto pset_material_steel = file.create<typename Schema::IfcMaterialProperties>().initialize(std::string("Pset_MaterialSteel"), std::nullopt/*description*/, material_steel_properties, material);
}

template <typename Schema>
void Create_Pset_MaterialSteel_ReinforcingBar(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial material, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   // Pset_MaterialSteel
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit stress_unit;

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

   std::vector<typename Schema::IfcProperty> material_steel_properties;
   //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("YieldStress"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fy), stress_unit));
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("UltimateStress"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fpu), stress_unit));
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("UltimateStrain"), std::nullopt, file.create<typename Schema::IfcPositiveRatioMeasure>().initialize(eu), typename Schema::IfcUnit{}));
   material_steel_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StructuralGrade"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(grade.c_str()), typename Schema::IfcUnit{}));
   auto material_properties = file.create<typename Schema::IfcMaterialProperties>().initialize(std::string("Pset_MaterialSteel"), std::nullopt/*description*/, material_steel_properties, material);
}

template <typename Schema>
typename Schema::IfcPropertySet Create_Pset_ConcreteElementGeneral(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, std::optional<std::string> assemblyPlace,std::optional<std::string> castingMethod,std::optional<Float64> fc)
{
   USES_CONVERSION;

   // Pset_ConcreteElementGeneral
   std::vector<typename Schema::IfcProperty> concrete_element_general_properties;

   // PEnum_AssemblyPlace
   if (assemblyPlace.has_value())
   {
      std::vector<std::string> assembly_place_enum_values{ "FACTORY","OFFSITE","SITE","OTHER","UNKNOWN","UNSET" };
      auto assembly_place_property_enum_values = createPropertyEnumeration<Schema>(file, "PEnum_AssemblyPlace", assembly_place_enum_values);
      auto assembly_place = createPropertyEnumeratedValue<Schema>(file, "AssemblyPlace", assembly_place_property_enum_values, *assemblyPlace);
      concrete_element_general_properties.push_back(assembly_place);
   }

   // PEnum_ConcreteCastingMethod
   if (castingMethod.has_value())
   {
      std::vector<std::string> casting_method_enum_values{ "INSITU","MIXED","PRECAST","PRINTED","OTHER","UNKNOWN","UNSET" };
      auto casting_method_property_enum_values = createPropertyEnumeration<Schema>(file, "PEnum_ConcreteCastingMethod", casting_method_enum_values);
      auto casting_method = createPropertyEnumeratedValue<Schema>(file, "CastingMethod", casting_method_property_enum_values, *castingMethod);
      concrete_element_general_properties.push_back(casting_method);
   }

   if (fc.has_value())
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto fc_value = *fc;

      std::ostringstream os;
      os << T2A(::FormatDimension(fc_value, pDisplayUnits->GetStressUnit())) << std::endl;
      concrete_element_general_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StrengthClass"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(os.str()), typename Schema::IfcUnit{}));
   }

   // create Pset_ConcreteElementGeneral
   auto pset_concrete_element_general = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_ConcreteElementGeneral"), std::nullopt, concrete_element_general_properties);
   return pset_concrete_element_general;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_Pset_BeamCommon(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
{
   GET_IFACE2(pBroker, IBridge, pBridge);

   auto span_length = pBridge->GetSegmentSpanLength(segmentKey);
   auto slope = pBridge->GetSegmentSlope(segmentKey);

   auto slope_angle = atan(slope);

   GET_IFACE2(pBroker, IGirder, pGirder);
   auto roll = pGirder->GetOrientation(segmentKey);
   auto roll_angle = atan(roll);


   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);
   typename Schema::IfcConversionBasedUnit length_unit;
   typename Schema::IfcConversionBasedUnit angle_unit;

   if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);
      angle_unit = GetAngleUnit<Schema>(file, pBroker);

      span_length = WBFL::Units::ConvertFromSysUnits(span_length, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      slope_angle = WBFL::Units::ConvertFromSysUnits(slope_angle, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
      roll_angle = WBFL::Units::ConvertFromSysUnits(roll_angle, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
   }


   std::vector<typename Schema::IfcProperty> list_of_properties;

   // Depreciated in IFC4.3. Use the Name attribute of the relating type
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Reference"), std::nullopt, {}, {}));

   std::vector<std::string> enum_values{ "DEMOLISH","EXISTING","NEW","TEMPORARY", "OTHER","UNKNOWN","UNSET"};
   auto property_enum_values = createPropertyEnumeration<Schema>(file, "PEnum_ElementStatus", enum_values);
   auto status = createPropertyEnumeratedValue<Schema>(file, "Status", property_enum_values, "NEW");
   list_of_properties.push_back(status);

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Span"), std::nullopt, file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(span_length), length_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Slope"), std::nullopt, file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(slope_angle), angle_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Roll"), std::nullopt, file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(roll_angle), angle_unit));

   // This properties are not applicable to bridges
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("IsExternal"), std::nullopt, {}, {}));
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ThermalTransmittance"), std::nullopt, {}, {}));
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LoadBearing"), std::nullopt, {}, {}));
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FireRating"), std::nullopt, {}, {}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_BeamCommon"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_Pset_PrecastConcreteElementGeneral(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
{
   USES_CONVERSION;

   GET_IFACE2_NOCHECK(pBroker, IBridge, pBridge);
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

   typename Schema::IfcConversionBasedUnit stress_unit;
   typename Schema::IfcConversionBasedUnit displacement_unit;
   typename Schema::IfcConversionBasedUnit angle_unit;

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

   GET_IFACE2(pBroker, IBridgeDescription, pBridgeDesc);
   auto family_name = pBridgeDesc->GetBridgeDescription()->GetGirderFamilyName();

   auto hauling_data = pBridgeDesc->GetBridgeDescription()->GetGirderGroup(segmentKey.groupIndex)->GetGirder(segmentKey.girderIndex)->GetSegment(segmentKey.segmentIndex)->HandlingData;
   auto bunk_point = std::max(hauling_data.LeadingSupportPoint, hauling_data.TrailingSupportPoint);

   typename Schema::IfcConversionBasedUnit length_unit;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);
      bunk_point = WBFL::Units::ConvertFromSysUnits(bunk_point, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TypeDesignation"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(T2A(family_name)), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CornerChamfer"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ManufacturingToleranceClass"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FormStrippingStrength"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fci), stress_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LiftingStrength"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fci), stress_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ReleaseStrength"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fci), stress_unit));

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MinimumAllowableSupportLength"), std::nullopt, file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(bunk_point), length_unit));

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("InitialTension"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fpj), stress_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TendonRelaxation"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TransportationStrength"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fc), stress_unit));

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SupportDuringTransportDescription"),
      std::nullopt, file.create<typename Schema::IfcText>().initialize("Assumed to be truck transported with bunking locations per MinimumAllowableSupportLength property"), typename Schema::IfcUnit{}));

   // not in bSDD
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertyReferenceValue>().initialize(std::string("SupportDuringTransportDocReference"),
   //   std::nullopt, std::nullopt,
   //   file.create<typename Schema::IfcDocumentReference>().initialize(
   //      std::string("https://www.pci.org/ItemDetail?iProductCode=CB-02-26H&Category=TRANSPORT&WebsiteKey=5a7b2064-98c2-4c8e-9b4b-18c80973da1e"), // Location
   //      std::nullopt, // Identification
   //      std::string("Recommended Practice for Lateral Stabiilty of Precast, Prestressed Concrete Bridge Girders, 2nd Edition (CB-02-26H)"), // Name
   //      std::nullopt, // Description
   //      {}))); // Referenced Document

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("HollowCorePlugging"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   if(options.include_camber)
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CamberAtMidspan"), std::nullopt, file.create<typename Schema::IfcRatioMeasure>().initialize(camber), typename Schema::IfcUnit{}));
   else
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CamberAtMidspan"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BatterAtStart"), std::nullopt, file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(batter), angle_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BatterAtEnd"), std::nullopt, file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(batter), angle_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Twisting"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Shortening"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("PieceMark"), std::nullopt, typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   {
#pragma Reminder("WORKING HERE - need to use the usBridge schema for this property")
      pgsAutoGirderLabel autoLabel;
      pgsGirderLabel::UseAlphaLabel(false);
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DesignLocationNumber"), std::nullopt, file.create<typename Schema::IfcLabel>().initialize(T2A(SEGMENT_LABEL(segmentKey))), typename Schema::IfcUnit{}));
   }

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("Pset_PrecastConcreteElementGeneral"), std::nullopt, list_of_properties);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ProjectCommon(hierarchy_helper<Schema>& file)
{
   // we don't have these properties, but will set up the property set with blank properties so that it shows up in the file and can be filled in by hand or by a future version of the exporter.
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ApprovalStatus"), BSDD_PROPERTY("ApprovalStatus"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FileNumber"), BSDD_PROPERTY("FileNumber"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LettingDate"), BSDD_PROPERTY("LettingDate"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ModelPreparationDate"), BSDD_PROPERTY("ModelPreparationDate"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ModelVersion"), BSDD_PROPERTY("ModelVersion"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ProjectDirectory"), BSDD_PROPERTY("ProjectDirectory"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ProjectIdentification"), BSDD_PROPERTY("ProjectIdentification"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ProjectNumber"), BSDD_PROPERTY("ProjectNumber"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ProjectURL"), BSDD_PROPERTY("ProjectURL"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("PSEData"), BSDD_PROPERTY("PSEData"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ProjectCommon"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ProjectLocation(hierarchy_helper<Schema>& file)
{
   // we don't have these properties, but will set up the property set with blank properties so that it shows up in the file and can be filled in by hand or by a future version of the exporter.
   std::vector<typename Schema::IfcProperty> list_of_properties;
#pragma Reminder("WORKING HERE - need to update the URLs - bSDD is down right now")
   // County and State are required properties.
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("City"), BSDD_PROPERTY("City"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("County"), BSDD_PROPERTY("County"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("District"), BSDD_PROPERTY("District"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("State"), BSDD_PROPERTY("State"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Section"), BSDD_PROPERTY("Section"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Township"), BSDD_PROPERTY("Township"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Range"), BSDD_PROPERTY("Range"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ProjectLocation"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_ProjectLocation"), list_of_properties);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_Common(hierarchy_helper<Schema>& file)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;

   std::vector<std::string> enum_values{ "NEW","EXISTING - REMAIN","EXISTING - REMOVE","TEMPORARY", "OTHER"};
   auto property_enum_values = createPropertyEnumeration<Schema>(file, "usBrPEnum_ElementStatus", enum_values); // creates an IfcPropertyEnumeration
   auto status = createPropertyEnumeratedValue<Schema>(file, "Status", property_enum_values, enum_values.front()); // creates an IfcPropertyEnumeratedValue
   status.setSpecification(BSDD_PROPERTY("Status"));
   list_of_properties.push_back(status);

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("AssociatedStandard"), BSDD_PROPERTY("AssociatedStandard"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Note"), BSDD_PROPERTY("Note"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SNBIElementNumber"), BSDD_PROPERTY("SNBIElementNumber"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("WorkingDrawingApproval"), BSDD_PROPERTY("WorkingDrawingApproval"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_Common"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_BridgePartCommon(hierarchy_helper<Schema>& file)
{
   // this property set is required by the IDS, but none of its properties are required.

   std::vector<typename Schema::IfcProperty> list_of_properties;

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndSkew"), BSDD_PROPERTY("EndSkew"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndStation"), BSDD_PROPERTY("EndStation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndStationOffset"), BSDD_PROPERTY("EndStationOffset"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartSkew"), BSDD_PROPERTY("StartSkew"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartStation"), BSDD_PROPERTY("StartStation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartStationOffset"), BSDD_PROPERTY("StartStationOffset"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_BridgePartCommon"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_SlabCommon(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatingNote"), BSDD_PROPERTY("CoatingNote"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CrossSectionalArea"), BSDD_PROPERTY("CrossSectionalArea"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CrossSlope"), BSDD_PROPERTY("CrossSlope"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndSkew"), BSDD_PROPERTY("EndSkew"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndStation"), BSDD_PROPERTY("EndStation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndStationOffset"), BSDD_PROPERTY("EndStationOffset"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartSkew"), BSDD_PROPERTY("StartSkew"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartStation"), BSDD_PROPERTY("StartStation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartStationOffset"), BSDD_PROPERTY("StartStationOffset"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_SlabCommon"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_SubstructureCommon(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, PierIndexType pierIdx = INVALID_INDEX)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;

   if (pierIdx != INVALID_INDEX)
   {
      GET_IFACE2(pBroker, IBridge, pBridge);
      auto pier_station = pBridge->GetPierStation(pierIdx);

      CComPtr<IAngle> angle;
      pBridge->GetPierSkew(pierIdx,&angle);
      Float64 skew;
      angle->get_Value(&skew);

      GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);
      typename Schema::IfcConversionBasedUnit station_unit;
      typename Schema::IfcConversionBasedUnit angle_unit;
      if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         station_unit = GetSpanLengthUnit<Schema>(file, pBroker);
         pier_station = WBFL::Units::ConvertFromSysUnits(pier_station, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);

         angle_unit = GetAngleUnit<Schema>(file, pBroker);
         skew = WBFL::Units::ConvertFromSysUnits(skew, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
         skew = RoundOff(skew, 0.0001);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StationAheadBearing"), BSDD_PROPERTY("StationAheadBearing"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StationAtCenterline"), BSDD_PROPERTY("PierStation"), file.create<typename Schema::IfcLengthMeasure>().initialize(pier_station), station_unit));
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StationBackBearing"), BSDD_PROPERTY("StationBackBearing"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StationOffsetAheadBearing"), BSDD_PROPERTY("StationOffsetAheadBearing"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StationOffsetAtCenterline"), BSDD_PROPERTY("StationOffsetAtCenterline"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StationOffsetBackBearing"), BSDD_PROPERTY("StationOffsetBackBearing"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SubstructureSkewAngle"), BSDD_PROPERTY("SubstructureSkewAngle"), file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(skew), angle_unit));
   }

   std::vector<std::string> enum_values{ "New","Other" }; // see SNBI
   auto property_enum_values = createPropertyEnumeration<Schema>(file, "usBrPEnum_SubstructureType", enum_values); // creates an IfcPropertyEnumeration
   auto substruture_type = createPropertyEnumeratedValue<Schema>(file, "SubstructureType", property_enum_values, enum_values.front()); // creates an IfcPropertyEnumeratedValue
   substruture_type.setSpecification(BSDD_PROPERTY("SubstructureType"));
   list_of_properties.push_back(substruture_type);

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ReturnInterval"), BSDD_PROPERTY("ReturnInterval"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ScourElevation"), BSDD_PROPERTY("ScourElevation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_SubstructureCommon"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_PayItemQuantities(hierarchy_helper<Schema>& file)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("PayItemReference"), BSDD_PROPERTY("PayItemReference"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("PayQuantity"), BSDD_PROPERTY("PayQuantity"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("UnitOfMeasure"), BSDD_PROPERTY("UnitOfMeasure"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_PayItemQuantities"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_BridgeGeometry(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
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

   typename Schema::IfcConversionBasedUnit length_unit;
   typename Schema::IfcConversionBasedUnit angle_unit;
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

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BridgeEndStation"), BSDD_PROPERTY("BridgeEndStation"), file.create<typename Schema::IfcLengthMeasure>().initialize(end_station), length_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BridgeEndStationOffset"), BSDD_PROPERTY("BridgeEndStationOffset"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BridgeLength"), BSDD_PROPERTY("BridgeLength"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(length), length_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BridgeSkew"), BSDD_PROPERTY("BridgeSkew"), file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(skew), angle_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BridgeStartStation"), BSDD_PROPERTY("BridgeStartStation"), file.create<typename Schema::IfcLengthMeasure>().initialize(start_station), length_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BridgeStartStationOffset"), BSDD_PROPERTY("BridgeStartStationOffset"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LowBeamElevation"), BSDD_PROPERTY("LowBeamElevation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("NumberOfSpans"), BSDD_PROPERTY("NumberOfSpans"), file.create<typename Schema::IfcInteger>().initialize((int64_t)nSpans), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("RoadwayWidth"), BSDD_PROPERTY("RoadwayWidth"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(width), length_unit));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_BridgeGeometry"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_BridgeIdentification(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SNBIBridgeNumber"), BSDD_PROPERTY("SNBIBridgeNumber"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SNBIBridgeName"), BSDD_PROPERTY("SNBIBridgeName"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   // other optional attributes to be provided later

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_BridgeIdentification"), std::nullopt, list_of_properties);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_DesignLoading(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, ILiveLoads, pLiveLoads);

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DesignMethodology"), BSDD_PROPERTY("DesignMethodology"), file.create<typename Schema::IfcLabel>().initialize("LRFD"), typename Schema::IfcUnit{}));

   if (pLiveLoads->IsLiveLoadDefined(pgsTypes::lltDesign))
   {
      auto live_load_names = pLiveLoads->GetLiveLoadNames(pgsTypes::LiveLoadType::lltDesign);
      auto live_load_name = live_load_names[0];
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("VehicularLiveLoad"), BSDD_PROPERTY("VehicularLiveLoad"), file.create<typename Schema::IfcLabel>().initialize(T2A(live_load_name.c_str())), typename Schema::IfcUnit{}));
   }

#pragma Reminder("WORKING HERE - finish filling out these properties")
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("PedestrianLiveLoad"), BSDD_PROPERTY("PedestrianLiveLoad"), {}, {}));
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DesignFutureWearingLoad"), BSDD_PROPERTY("DesignFutureWearingLoad"), {}, {}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_DesignLoading"), std::nullopt, list_of_properties);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_FeatureIdentification(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FeatureType"), BSDD_PROPERTY("FeatureType"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FeatureLocation"), BSDD_PROPERTY("FeatureLocation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FeatureName"), BSDD_PROPERTY("FeatureName"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   // other optional attributes to be provided later

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_FeatureIdentification"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
void Create_usBrPset_HydraulicData(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_NavigableWaterway(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
void Create_usBrPset_Railroad(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge bridge)
{
   // placeholder - to be implemented later
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_Roadway(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridge bridge)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("AnnualAverageDailyTraffic"), BSDD_PROPERTY("AnnualAverageDailyTraffic"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("AnnualAverageDailyTruckTraffic"), BSDD_PROPERTY("AnnualAverageDailyTruckTraffic"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DesignAnnualAverageDailyTraffic"), BSDD_PROPERTY("DesignAnnualAverageDailyTraffic"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DesignAnnualAverageDailyTruckTraffic"), BSDD_PROPERTY("DesignAnnualAverageDailyTruckTraffic"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DesignTrafficYear"), BSDD_PROPERTY("DesignTrafficYear"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MinimumHorizontalClearanceLeft"), BSDD_PROPERTY("MinimumHorizontalClearanceLeft"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MinimumHorizontalClearanceRight"), BSDD_PROPERTY("MinimumHorizontalClearanceRight"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("RoadwayMinimumVerticalClearance"), BSDD_PROPERTY("RoadwayMinimumVerticalClearance"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("RoadwayName"), BSDD_PROPERTY("RoadwayName"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("RoadwayType"), BSDD_PROPERTY("RoadwayType"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_Roadway"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_Roadway"), list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_RoadwaySlab(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("FormType"), BSDD_PROPERTY("FormType"), file.create<typename Schema::IfcLabel>().initialize("Wood"), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SacrificialWearingSurface"), BSDD_PROPERTY("SacrificialWearingThickness"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SurfaceArea"), BSDD_PROPERTY("SurfaceArea"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SurfaceFinish"), BSDD_PROPERTY("SurfaceFinish"), file.create<typename Schema::IfcLabel>().initialize("Raked"), typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_RoadwaySlab"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_RoadwaySlab"), list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_MASH(hierarchy_helper<Schema>& file)
{
#pragma Reminder("WORKING HERE - need to get real barrier properties and names for PGSuper - dummy values used for now")
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MASHTestingLevel"), BSDD_PROPERTY("MASHTestingLevel"), file.create<typename Schema::IfcLabel>().initialize("Unknown"), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("GuardrailType"), BSDD_PROPERTY("GuardrailType"), file.create<typename Schema::IfcLabel>().initialize("Unknown"), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarrierType"), BSDD_PROPERTY("BarrierType"), file.create<typename Schema::IfcLabel>().initialize("Unknown"), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarrierHeight"), BSDD_PROPERTY("BarrierHeight"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(0.01), typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_MASH"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_MASH"), list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_TendonDebondingAtEnds(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, double debond_start,double debond_end)
{
   // Debonding at ends of beam (normal debonding)

   typename Schema::IfcConversionBasedUnit length_unit;

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      debond_start = WBFL::Units::ConvertFromSysUnits(debond_start, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      debond_end = WBFL::Units::ConvertFromSysUnits(debond_end, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DebondLengthStart"), BSDD_PROPERTY("DebondLengthStart"), file.create<typename Schema::IfcLengthMeasure>().initialize(debond_start), length_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DebondLengthEnd"), BSDD_PROPERTY("DebondLengthEnd"), file.create<typename Schema::IfcLengthMeasure>().initialize(debond_end), length_unit));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_TendonDebonding"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_TendonDebonding"), list_of_properties);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_TendonDebondingInCenter(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, double debond_start, double debond_end)
{
   // Debonding in the middle of the beam (typically for temporary top strands)
   typename Schema::IfcConversionBasedUnit length_unit;

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      debond_start = WBFL::Units::ConvertFromSysUnits(debond_start, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      debond_end = WBFL::Units::ConvertFromSysUnits(debond_end, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;

   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DebondLengthMidspanToStart"), BSDD_PROPERTY("DebondLengthMidspanToStart"), file.create<typename Schema::IfcLengthMeasure>().initialize(debond_start), length_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DebondLengthMidspanToEnd"), BSDD_PROPERTY("DebondLengthMidspanToEnd"), file.create<typename Schema::IfcLengthMeasure>().initialize(debond_end), length_unit));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_TendonDebonding"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_TendonDebonding"), list_of_properties);

   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_PrecastConcreteBeam(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
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
   typename Schema::IfcConversionBasedUnit deflection_unit;
   if(options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      deflection_unit = GetDisplacementUnit<Schema>(file, pBroker);
      initial_camber = WBFL::Units::ConvertFromSysUnits(initial_camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
      final_camber = WBFL::Units::ConvertFromSysUnits(final_camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
      screed_camber = WBFL::Units::ConvertFromSysUnits(screed_camber, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DeflectionLongTerm"), BSDD_PROPERTY("DeflectionLongTerm"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DeflectionShortTerm"), BSDD_PROPERTY("DeflectionShortTerm"), file.create<typename Schema::IfcLengthMeasure>().initialize(screed_camber), deflection_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ElasticShortening"), BSDD_PROPERTY("ElasticShortening"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("LiftingLoopLocation"), BSDD_PROPERTY("LiftingLoopLocation"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MidSpanCamberAfterLosses"), BSDD_PROPERTY("MidSpanCamberAfterLosses"), file.create<typename Schema::IfcLengthMeasure>().initialize(final_camber), deflection_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MidSpanCamberAtRelease"), BSDD_PROPERTY("MidSpanCamberAtRelease"), file.create<typename Schema::IfcLengthMeasure>().initialize(initial_camber), deflection_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MinimumTimeToDeckPlacement"), BSDD_PROPERTY("MinimumTimetoDeckPlacement"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ShapeName"), BSDD_PROPERTY("ShapeName"), file.create<typename Schema::IfcLabel>().initialize(T2A(shape_name)), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TopSurfaceFinish"), BSDD_PROPERTY("TopSurfaceFinish"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_PrecastConcreteBeam"), std::nullopt, list_of_properties);

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
void Create_usBrPset_ACI_TendonMaterial(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial material, const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;

   auto spec = GetStrandSpecification(pStrand);
   auto spec_edition = GetStrandSpecificationEdition(pStrand);

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit stress_unit;
   auto fpu = pStrand->GetUltimateStrength();

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      fpu = WBFL::Units::ConvertFromSysUnits(fpu, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Specification"), BSDD_PROPERTY("Specification"), file.create<typename Schema::IfcLabel>().initialize(spec.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SpecificationVersion"), BSDD_PROPERTY("SpecificationVersion"), file.create<typename Schema::IfcLabel>().initialize(spec_edition.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TendonGrade"), BSDD_PROPERTY("TendonGrade"), file.create<typename Schema::IfcPressureMeasure>().initialize(fpu), stress_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatingSpecification"), BSDD_PROPERTY("CoatingSpecification"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatingSpecificationVersion"), BSDD_PROPERTY("CoatingSpecificationVersion"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto material_properties = file.create<typename Schema::IfcMaterialProperties>().initialize(std::string("usBrPset_ACI_TendonMaterial"), std::nullopt/*description*/, list_of_properties, material);
}

template <typename Schema>
void Create_usBrPset_ACI_ReinforcingMaterial(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, typename Schema::IfcMaterial material, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   auto spec = GetRebarSpecification(pRebar);
   auto spec_edition = GetRebarSpecificationEdition(pRebar);

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit stress_unit;
   auto fy = pRebar->GetYieldStrength();

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      fy = WBFL::Units::ConvertFromSysUnits(fy, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Specification"), BSDD_PROPERTY("Specification"), file.create<typename Schema::IfcLabel>().initialize(spec.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SpecificationVersion"), BSDD_PROPERTY("SpecificationVersion"), file.create<typename Schema::IfcLabel>().initialize(spec_edition.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("ReinforcingGrade"), BSDD_PROPERTY("ReinforcingGrade"), file.create<typename Schema::IfcPressureMeasure>().initialize(fy), stress_unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("Subtype"), BSDD_PROPERTY("Subtype"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatingSpecification"), BSDD_PROPERTY("CoatingSpecification"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatingSpecificationVersion"), BSDD_PROPERTY("CoatingSpecificationVersion"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatingSubtype"), BSDD_PROPERTY("CoatingSubtype"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CoatedBeforeFabrication"), BSDD_PROPERTY("CoatedBeforeFabrication"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto material_properties = file.create<typename Schema::IfcMaterialProperties>().initialize(std::string("usBrPset_ACI_ReinforcingMaterial"), std::nullopt/*description*/, list_of_properties, material);
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ACI_ReinforcingBarType(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::string mark, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   auto size = WBFL::LRFD::RebarPool::GetBarSize(pRebar->GetSize());
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarMark"), BSDD_PROPERTY("BarMark"), file.create<typename Schema::IfcLabel>().initialize(mark.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarMass"), BSDD_PROPERTY("BarMass"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarSize"), BSDD_PROPERTY("BarSize"), file.create<typename Schema::IfcLabel>().initialize(T2A(size.c_str())), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndEndPrep"), BSDD_PROPERTY("EndEndPrep"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartEndPrep"), BSDD_PROPERTY("StartEndPrep"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_ReinforcingBarType"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
void addBarDimension(std::string name,std::string uri,double dim,typename Schema::IfcConversionBasedUnit unit, hierarchy_helper<Schema>& file, std::vector<typename Schema::IfcProperty>& list_of_properties)
{
   // ACI 131 says to use IfcLengthMeasure for distances and IfcPlaneAngleMeasure for angles. IfcReal is for nondimensional real values
   // but usBridge has distances as IfcReal.
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(name, uri, file.create<typename Schema::IfcLengthMeasure>().initialize(dim), unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(name, uri, file.create<typename Schema::IfcReal>().initialize(dim), unit));
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ACI_BarShape(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::string bend_shape_name,
   double bend_radius, const std::unordered_map<std::string, double>& dimensions)
{
   std::string standard = "ACI 315-99";
   std::string standard_version = "1999";

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StandardName"), BSDD_PROPERTY("StandardName"), file.create<typename Schema::IfcLabel>().initialize(standard.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StandardVersion"), BSDD_PROPERTY("StandardVersion"), file.create<typename Schema::IfcLabel>().initialize(standard_version.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BendShapeName"), BSDD_PROPERTY("BendShapeName"), file.create<typename Schema::IfcLabel>().initialize(bend_shape_name.c_str()), typename Schema::IfcUnit{}));

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit dimension_unit;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      dimension_unit = GetComponentDimUnit<Schema>(file, pBroker);
      bend_radius = WBFL::Units::ConvertFromSysUnits(bend_radius, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
   }
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DefaultInsideBendRadius"), BSDD_PROPERTY("DefaultInsideBendRadius"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(bend_radius), dimension_unit));

   for (auto & [name, dim] : dimensions)
   {
      auto value = dim;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
      }
      addBarDimension<Schema>(name, "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/prop/" + name, value, dimension_unit, file, list_of_properties);
   }

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_BarShape"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ACI_ReinforcingBar(hierarchy_helper<Schema>& file, std::string element,std::string use,std::string position)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarElement"), BSDD_PROPERTY("BarElement"), file.create<typename Schema::IfcLabel>().initialize(element.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarUse"), BSDD_PROPERTY("BarUse"), file.create<typename Schema::IfcLabel>().initialize(use.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarPosition"), BSDD_PROPERTY("BarPosition"), file.create<typename Schema::IfcLabel>().initialize(position.c_str()), typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_ReinforcingBar"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_ACI_ReinforcingBar"), list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_Reinforcing(hierarchy_helper<Schema>& file)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EmbedmentAtEnd"), BSDD_PROPERTY("EmbedmentAtEnd"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EmbedmentAtStart"), BSDD_PROPERTY("EmbedmentAtStart"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MinimumSpliceLength"), BSDD_PROPERTY("MinimumSpliceLength"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("PlacementMethod"), BSDD_PROPERTY("PlacementMethod"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_Reinforcing"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_Reinforcing"), list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ReinforcingCover(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::optional<double> top, std::optional<double> side, std::optional<double> bottom, std::optional<double> end)
{
   if (!top && !side && !end && !bottom)
      return {};

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit cover_unit;
   auto length_unit = pDisplayUnits->GetComponentDimUnit();
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      cover_unit = GetComponentDimUnit<Schema>(file, pBroker);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;

   if (top)
   {
      double value = *top;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*top, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TopFaceCover"), BSDD_PROPERTY("TopFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   if (side)
   {
      double value = *side;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*side, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SideFaceCover"), BSDD_PROPERTY("SideFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   if (end)
   {
      double value = *end;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*end, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndFaceCover"), BSDD_PROPERTY("EndFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   if (bottom)
   {
      double value = *bottom;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*bottom, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BottomFaceCover"), BSDD_PROPERTY("BottomFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_ReinforcingCover"), std::nullopt, list_of_properties);
   return property_set;
}
