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

#include "PGSuperColors.h"

std::string GetStrandMaterialName(const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;
   return T2A(pStrand->GetName().c_str());
}

template <typename Schema>
typename Schema::IfcStyledRepresentation* CreateMaterialRepresentation(IfcHierarchyHelper<Schema>& file, std::string name, COLORREF clr)
{
   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   double r = (double)GetRValue(clr) / 255.;
   double g = (double)GetGValue(clr) / 255.;
   double b = (double)GetBValue(clr) / 255.;

   auto color = new Schema::IfcColourRgb(name, r, g, b);
   file.addEntity(color);

   auto ssr = new Schema::IfcSurfaceStyleRendering(color, boost::none, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, Schema::IfcReflectanceMethodEnum::IfcReflectanceMethod_NOTDEFINED);
   file.addEntity(ssr);

   typename aggregate_of<typename Schema::IfcSurfaceStyleElementSelect>::ptr list_of_surface_styles(new aggregate_of<typename Schema::IfcSurfaceStyleElementSelect>());
   list_of_surface_styles->push(ssr);
   
   auto ss = new Schema::IfcSurfaceStyle(name, Schema::IfcSurfaceSide::IfcSurfaceSide_BOTH, list_of_surface_styles);
   file.addEntity(ss);

   typename aggregate_of<typename Schema::IfcPresentationStyle>::ptr list_of_presentation_styles(new aggregate_of<typename Schema::IfcPresentationStyle>());
   list_of_presentation_styles->push(ss);
   
   auto styled_item = new Schema::IfcStyledItem(nullptr, list_of_presentation_styles, boost::none);
   file.addEntity(styled_item);

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr styled_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   styled_items->push(styled_item);
   auto styled_representation = new Schema::IfcStyledRepresentation(geometric_representation_context, boost::none, boost::none, styled_items);
   file.addEntity(styled_representation);

   return styled_representation;
}

template <typename Schema>
void AssignPset_MaterialSteel(IfcHierarchyHelper<Schema>& file, typename Schema::IfcMaterial* material,Float64 fy, Float64 fpu, Float64 eu, const std::string& grade)
{
   // Pset_MaterialSteel
   typename aggregate_of<typename Schema::IfcProperty>::ptr material_steel_properties(new aggregate_of<typename Schema::IfcProperty>());
   //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
   material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("YieldStress"), boost::none, new Schema::IfcPressureMeasure(fy), nullptr));
   material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStress"), boost::none, new Schema::IfcPressureMeasure(fpu), nullptr));
   material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStrain"), boost::none, new Schema::IfcPositiveRatioMeasure(eu), nullptr));
   material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("StructuralGrade"), boost::none, new Schema::IfcLabel(grade.c_str()), nullptr));
   auto pset_material_steel = new Schema::IfcMaterialProperties(std::string("Pset_MaterialSteel"), boost::none/*description*/, material_steel_properties, material);
   file.addEntity(pset_material_steel);
}

template <typename Schema>
typename Schema::IfcMaterial* GetStrandMaterial(IfcHierarchyHelper<Schema>& file,const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;

   auto name = GetStrandMaterialName(pStrand);

   // search to see if an IfcMaterial for this kind of strand has already been created
   auto materials = file.instances_by_type<typename Schema::IfcMaterial>();
   for (auto material : *materials)
   {
      if (material->Name() == name)
         return material;
   }

   // if we got this far, the material was not previously created
   // create it now
   auto strand_material = new Schema::IfcMaterial(name, boost::none/*description*/, std::string("steel")/*category*/);
   file.addEntity(strand_material);

   // define strand properties
   auto fy = pStrand->GetYieldStrength();
   auto fpu = pStrand->GetUltimateStrength();
   auto eu = 0.035; // from ASTM A416 spec

   // Need to clean this up
   // ASTM A416 is for low relaxation strand... PGSuper does low relaxation and stress relieved
   // ASTM A416 is for Grade 250 and Grade 270... PGSuper does grade 300 as well, but there doesn't seem to be an ASTM
   // We are assuming same material for all strands, but that is not the case in the PGSuper data model
   // straight, harped, and temporary can be different - Grade 250, Grade 270, Grade 300
   // Strand size/diameter is a property on IfcTendon
   std::ostringstream os;
   os << "ASTM A416 Grade " << T2A(WBFL::Materials::PsStrand::GetGrade(pStrand->GetGrade(), true/*US units*/).c_str());
   auto grade = os.str();

   // Pset_MaterialSteel
   AssignPset_MaterialSteel(file, strand_material, fy, fpu, eu, grade);

   // create the representation style
   auto material_representation = CreateMaterialRepresentation<Schema>(file, "Strand", STRAND_BORDER_COLOR);

   // assigns the presentation styles to the material
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
   list_of_representations->push(material_representation);
   auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, strand_material);
   file.addEntity(material_defintion_representation);

   return strand_material;
}

std::string GetRebarMaterialName(const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;
   return T2A(WBFL::LRFD::RebarPool::GetMaterialName(pRebar->GetType(),pRebar->GetGrade()).c_str());
}

template <typename Schema>
typename Schema::IfcMaterial* GetRebarMaterial(IfcHierarchyHelper<Schema>& file, const WBFL::Materials::Rebar* pRebar,const std::string& styleName,COLORREF color)
{
   USES_CONVERSION;

   auto name = GetRebarMaterialName(pRebar);

   // search to see if an IfcMaterial for this kind of strand has already been created
   auto materials = file.instances_by_type<typename Schema::IfcMaterial>();
   for (auto material : *materials)
   {
      if (material->Name() == name)
         return material;
   }


   // if we got this far, the material was not previously created
   // create it now
   auto rebar_material = new Schema::IfcMaterial(name, boost::none/*description*/, std::string("steel")/*category*/);
   file.addEntity(rebar_material);

   // define rebar properties
   auto fy = pRebar->GetYieldStrength();
   auto fpu = pRebar->GetUltimateStrength();
   auto eu = pRebar->GetElongation(); // depends on bar size and we are using a dummy #3 bar

   std::ostringstream os;
   os << T2A(pRebar->GetName().c_str());
   auto grade = os.str();

   // Pset_MaterialSteel
   AssignPset_MaterialSteel(file, rebar_material, fy, fpu, eu, grade);

   // create the representation style
   auto material_representation = CreateMaterialRepresentation<Schema>(file, styleName, color);

   // assigns the presentation styles to the material
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
   list_of_representations->push(material_representation);
   auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, rebar_material);
   file.addEntity(material_defintion_representation);

   return rebar_material;
}

template <typename Schema>
typename Schema::IfcMaterial* GetConcreteMaterial(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, Float64 fc, Float64 max_agg_size, const std::string& styleName, COLORREF color)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   std::ostringstream os;
   os << "Precast Concrete, f'c = " << T2A((LPCTSTR)(::FormatDimension(fc, pDisplayUnits->GetStressUnit())));

   auto name = os.str();

   // search to see if an IfcMaterial for this kind of strand has already been created
   auto materials = file.instances_by_type<typename Schema::IfcMaterial>();
   for (auto material : *materials)
   {
      if (material->Name() == name)
         return material;
   }


   // if we got this far, the material was not previously created
   // create it now
   auto concrete_material = new Schema::IfcMaterial(name, boost::none/*description*/, std::string("concrete")/*category*/);
   file.addEntity(concrete_material);

   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* displacement_unit = nullptr;
   if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);

      fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      max_agg_size = WBFL::Units::ConvertFromSysUnits(max_agg_size, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
   }

   // Pset_MaterialConcrete
   typename aggregate_of<typename Schema::IfcProperty>::ptr material_concrete_properties(new aggregate_of<typename Schema::IfcProperty>());
   material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("CompressiveStrength"), boost::none, new Schema::IfcPressureMeasure(fc), stress_unit));
   material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("MaxAggregateSize"), boost::none, new Schema::IfcPositiveLengthMeasure(max_agg_size), displacement_unit));
   auto pset_material_concrete = new Schema::IfcMaterialProperties(std::string("Pset_MaterialConcrete"), boost::none/*description*/, material_concrete_properties, concrete_material);
   file.addEntity(pset_material_concrete);

   // create the representation style
   auto material_representation = CreateMaterialRepresentation<Schema>(file, styleName, color);

   // assigns the presentation styles to the material
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
   list_of_representations->push(material_representation);
   auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, concrete_material);
   file.addEntity(material_defintion_representation);

   return concrete_material;
}
