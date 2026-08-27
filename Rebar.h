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

inline double getMinBendRadius(const WBFL::Materials::Rebar* pRebar,bool bStirrup)
{
   auto radius = 0.0;
   auto db = pRebar->GetNominalDimension();
   if (pRebar->GetSize() <= WBFL::Materials::Rebar::Size::bs5 && bStirrup)
      radius = 4 * db;
   else if (pRebar->GetSize() <= WBFL::Materials::Rebar::Size::bs8)
      radius = 6 * db;
   else if (pRebar->GetSize() <= WBFL::Materials::Rebar::Size::bs11)
      radius = 8 * db;
   else
      radius = 10 * db;

   return radius;
}

template <typename Schema>
typename Schema::IfcReinforcingBarType GetReinforcingBarType(hierarchy_helper<Schema>& file, const std::string& name, bool bStirrup, const WBFL::Materials::Rebar* pRebar)
{
   // search to see if an IfcReinforcingBarType has already been created
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   for (auto& rel_declares : rel_declares_instances)
   {
      if (rel_declares.RelatingContext().template as<typename Schema::IfcProject>())
      {
         auto related_definitions = rel_declares.RelatedDefinitions();
         for (auto& reldef : related_definitions)
         {
            auto rebar_type = reldef.template as<typename Schema::IfcReinforcingBarType>();
            if (rebar_type && rebar_type.Name() == name)
            {
               return rebar_type;
            }
         }
      }
   }

   return {};
}

template <typename Schema>
typename Schema::IfcReinforcingBarType CreateReinforcingBarType(hierarchy_helper<Schema>& file, const CIfcExportOptions& options, const std::string& name, const WBFL::Materials::Rebar* pRebar, typename Schema::IfcReinforcingBarTypeEnum::Value type, typename Schema::IfcShapeRepresentation shape_representation)
{
   auto placement = file.addPlacement3d();
   auto representation_map = file.create<typename Schema::IfcRepresentationMap>().initialize(placement, shape_representation);
   std::vector<typename Schema::IfcRepresentationMap> representation_maps;
   representation_maps.push_back(representation_map);

   // if we get this far, we need a new IfcReinforcingBarType
   auto rebar_type = file.create<typename Schema::IfcReinforcingBarType>().initialize(
      ifcopenshell::global_id(),
      {},
      name, /*Name*/
      std::nullopt, /*Description*/
      std::nullopt, /*ApplicableOccurrence*/
      std::nullopt, /*HasPropertySets*/
      representation_maps, /*RepresentationMaps*/
      std::nullopt, /*Tag*/
      std::nullopt, /*ElementType*/
      type, /*PredefinedType*/
      pRebar->GetNominalDimension(), /*NominalDiameter*/
      pRebar->GetNominalArea(), /*CrossSectionArea*/
      std::nullopt, /*BarLength*/
      std::nullopt, /*BarSurface*/
      std::nullopt, /*BendingShapeCode*/
      std::nullopt /*BendingParameters*/
   );



   if (options.classify)
   {
      Classify_usBridge_ReinforcingBarType(file, rebar_type);
   }

   // add the new definition to the project
   auto project = file.getSingle<typename Schema::IfcProject>();
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   if (rel_declares_instances.size() == 0)
   {
      std::vector<typename Schema::IfcDefinitionSelect> related_definitions;
      related_definitions.push_back(rebar_type);

      auto rel_declares = file.create<typename Schema::IfcRelDeclares>().initialize(
         ifcopenshell::global_id(),
         {},
         std::nullopt,
         std::nullopt,
         project,
         related_definitions);


   }
   else
   {
      for (auto& rel_declares : rel_declares_instances)
      {
         if (rel_declares.RelatingContext().template as<typename Schema::IfcProject>())
         {
            auto related_definitions = rel_declares.RelatedDefinitions();
            related_definitions.push_back(rebar_type);
            rel_declares.setRelatedDefinitions(related_definitions);
            break;
         }
      }
   }

   return rebar_type;
}
