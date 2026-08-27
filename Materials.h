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

#include "PGSuperColors.h"

std::string GetStrandMaterialName(const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;
   return T2A(pStrand->GetName().c_str());
}

template <typename Schema>
typename Schema::IfcStyledRepresentation CreateMaterialRepresentation(hierarchy_helper<Schema>& file, std::string name, COLORREF clr)
{
   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   double r = (double)GetRValue(clr) / 255.;
   double g = (double)GetGValue(clr) / 255.;
   double b = (double)GetBValue(clr) / 255.;

   auto color = file.create<typename Schema::IfcColourRgb>().initialize(name, r, g, b);


   auto ssr = file.create<typename Schema::IfcSurfaceStyleRendering>().initialize(color, std::nullopt, {}, {}, {}, {}, {}, {}, Schema::IfcReflectanceMethodEnum::IfcReflectanceMethod_NOTDEFINED);


   std::vector<typename Schema::IfcSurfaceStyleElementSelect> list_of_surface_styles;
   list_of_surface_styles.push_back(ssr);
   
   auto ss = file.create<typename Schema::IfcSurfaceStyle>().initialize(name, Schema::IfcSurfaceSide::IfcSurfaceSide_BOTH, list_of_surface_styles);


   std::vector<typename Schema::IfcPresentationStyle> list_of_presentation_styles;
   list_of_presentation_styles.push_back(ss);
   
   auto styled_item = file.create<typename Schema::IfcStyledItem>().initialize({}, list_of_presentation_styles, std::nullopt);


   std::vector<typename Schema::IfcRepresentationItem> styled_items;
   styled_items.push_back(styled_item);
   auto styled_representation = file.create<typename Schema::IfcStyledRepresentation>().initialize(geometric_representation_context, std::nullopt, std::nullopt, styled_items);


   return styled_representation;
}

template <typename Schema>
typename Schema::IfcMaterial GetStrandMaterial(hierarchy_helper<Schema>& file,std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;

   auto name = GetStrandMaterialName(pStrand);

   // search to see if an IfcMaterial for this kind of strand has already been created
   auto materials = file.instances_by_type<typename Schema::IfcMaterial>();
   for (auto& material : materials)
   {
      if (material.Name() == name)
         return material;
   }

   // if we got this far, the material was not previously created
   // create it now
   auto strand_material = file.create<typename Schema::IfcMaterial>().initialize(name, std::nullopt/*description*/, std::string("steel")/*category*/);


   // Pset_MaterialSteel
   Create_Pset_MaterialSteel_Strand(file, pBroker, options, strand_material, pStrand);

   if(options.classify)
   {
      Create_usBrPset_ACI_TendonMaterial(file, pBroker, options, strand_material, pStrand);
   }

   if (!options.classify)
   {
      // Material style representations are not part of AbV so don't include 
      // if we are exporting with usBridge classification system

      // create the representation style
      auto material_representation = CreateMaterialRepresentation<Schema>(file, "Strand", STRAND_BORDER_COLOR);

      // assigns the presentation styles to the material
      std::vector<typename Schema::IfcRepresentation> list_of_representations;
      list_of_representations.push_back(material_representation);
      auto material_defintion_representation = file.create<typename Schema::IfcMaterialDefinitionRepresentation>().initialize(std::nullopt, std::nullopt, list_of_representations, strand_material);

   }

   return strand_material;
}

std::string GetRebarMaterialName(const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;
   return T2A(WBFL::LRFD::RebarPool::GetMaterialName(pRebar->GetType(),pRebar->GetGrade()).c_str());
}

template <typename Schema>
typename Schema::IfcMaterial GetRebarMaterial(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const WBFL::Materials::Rebar* pRebar,const std::string& styleName,COLORREF color)
{
   USES_CONVERSION;

   auto name = GetRebarMaterialName(pRebar);

   // search to see if an IfcMaterial for this kind of strand has already been created
   auto materials = file.instances_by_type<typename Schema::IfcMaterial>();
   for (auto& material : materials)
   {
      if (material.Name() == name)
         return material;
   }


   // if we got this far, the material was not previously created
   // create it now
   auto rebar_material = file.create<typename Schema::IfcMaterial>().initialize(name, std::nullopt/*description*/, std::string("steel")/*category*/);



   // Pset_MaterialSteel
   Create_Pset_MaterialSteel_ReinforcingBar(file, pBroker, options, rebar_material, pRebar);
   Create_usBrPset_ACI_ReinforcingMaterial(file, pBroker, options, rebar_material, pRebar);


   if (!options.classify)
   {
      // Material style representations are not part of AbV so don't include 
      // if we are exporting with usBridge classification system
       
      // create the representation style
      auto material_representation = CreateMaterialRepresentation<Schema>(file, styleName, color);

      // assigns the presentation styles to the material
      std::vector<typename Schema::IfcRepresentation> list_of_representations;
      list_of_representations.push_back(material_representation);
      auto material_defintion_representation = file.create<typename Schema::IfcMaterialDefinitionRepresentation>().initialize(std::nullopt, std::nullopt, list_of_representations, rebar_material);

   }

   return rebar_material;
}

template <typename Schema>
typename Schema::IfcMaterial GetConcreteMaterial(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, Float64 fc, Float64 max_agg_size, const std::string& styleName, COLORREF color)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   std::ostringstream os;
   os << "Precast Concrete, f'c = " << T2A((LPCTSTR)(::FormatDimension(fc, pDisplayUnits->GetStressUnit())));

   auto name = os.str();

   // search to see if an IfcMaterial for this kind of strand has already been created
   auto materials = file.instances_by_type<typename Schema::IfcMaterial>();
   for (auto& material : materials)
   {
      if (material.Name() == name)
         return material;
   }


   // if we got this far, the material was not previously created
   // create it now
   auto concrete_material = file.create<typename Schema::IfcMaterial>().initialize(name, std::nullopt/*description*/, std::string("concrete")/*category*/);


   typename Schema::IfcConversionBasedUnit stress_unit;
   typename Schema::IfcConversionBasedUnit displacement_unit;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      stress_unit = GetStressUnit<Schema>(file, pBroker);
      displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);

      fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      max_agg_size = WBFL::Units::ConvertFromSysUnits(max_agg_size, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
   }

   // Pset_MaterialConcrete
   std::vector<typename Schema::IfcProperty> material_concrete_properties;
   material_concrete_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("CompressiveStrength"), std::nullopt, file.create<typename Schema::IfcPressureMeasure>().initialize(fc), stress_unit));
   material_concrete_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("MaxAggregateSize"), std::nullopt, file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(max_agg_size), displacement_unit));
   auto pset_material_concrete = file.create<typename Schema::IfcMaterialProperties>().initialize(std::string("Pset_MaterialConcrete"), std::nullopt/*description*/, material_concrete_properties, concrete_material);


   if (!options.classify)
   {
      // Material style representations are not part of AbV so don't include 
      // if we are exporting with usBridge classification system

      // create the representation style
      auto material_representation = CreateMaterialRepresentation<Schema>(file, styleName, color);

      // assigns the presentation styles to the material
      std::vector<typename Schema::IfcRepresentation> list_of_representations;
      list_of_representations.push_back(material_representation);
      auto material_defintion_representation = file.create<typename Schema::IfcMaterialDefinitionRepresentation>().initialize(std::nullopt, std::nullopt, list_of_representations, concrete_material);

   }

   return concrete_material;
}

template <typename Schema> 
void AssociateMaterial(hierarchy_helper<Schema>& file, typename Schema::IfcMaterial material, typename Schema::IfcProduct product)
{
   auto associations = product.HasAssociations();
   for (auto& association : associations)
   {
      auto rel_material = association.template as<typename Schema::IfcRelAssociatesMaterial>();
      if (rel_material)
      {
         auto related_objects = rel_material.RelatedObjects();
         // this product already has a material association, so we will just update it to point to the new material
         rel_material.setRelatingMaterial(material);
         return;
      }
   }

   // if we get this far, there was not already a material association for this product, so we will create a new one
   std::vector<typename Schema::IfcDefinitionSelect> related_objects;
   related_objects.push_back(product);
   auto rel_material = file.create<typename Schema::IfcRelAssociatesMaterial>().initialize(
      ifcopenshell::global_id(),
      {},
      std::nullopt, // Name
      std::nullopt, // Description
      related_objects, // RelatedObjects
      material // RelatingMaterial
   );

}
