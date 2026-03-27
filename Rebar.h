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


template <typename Schema>
typename Schema::IfcReinforcingBarType* GetReinforcingBarType(IfcHierarchyHelper<Schema>& file, const std::string& name, bool bStirrup, const WBFL::Materials::Rebar* pRebar)
{
   // search to see if an IfcReinforcingBarType has already been created
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   for (auto& rel_declares : *rel_declares_instances)
   {
      if (rel_declares->RelatingContext()->as<typename Schema::IfcProject>())
      {
         auto related_definitions = rel_declares->RelatedDefinitions();
         for (auto& reldef : *related_definitions)
         {
            auto rebar_type = reldef->as<typename Schema::IfcReinforcingBarType>();
            if (rebar_type && rebar_type->Name() == name)
            {
               return rebar_type;
            }
         }
      }
   }

   return nullptr;
}

template <typename Schema>
typename Schema::IfcReinforcingBarType* CreateReinforcingBarType(IfcHierarchyHelper<Schema>& file, const std::string& name, const WBFL::Materials::Rebar* pRebar, typename Schema::IfcReinforcingBarTypeEnum type, typename Schema::IfcShapeRepresentation* shape_representation)
{
   auto placement = file.addPlacement3d();
   auto representation_map = new typename Schema::IfcRepresentationMap(placement, shape_representation);
   typename Schema::IfcRepresentationMap::list::ptr representation_maps(new typename Schema::IfcRepresentationMap::list);
   representation_maps->push(representation_map);

   // if we get this far, we need a new IfcReinforcingBarType
   auto rebar_type = new typename Schema::IfcReinforcingBarType(
      IfcParse::IfcGlobalId(),
      nullptr,
      name, /*Name*/
      boost::none, /*Description*/
      boost::none, /*ApplicableOccurrence*/
      boost::none, /*HasPropertySets*/
      representation_maps, /*RepresentationMaps*/
      boost::none, /*Tag*/
      boost::none, /*ElementType*/
      type, /*PredefinedType*/
      pRebar->GetNominalDimension(), /*NominalDiameter*/
      pRebar->GetNominalArea(), /*CrossSectionArea*/
      boost::none, /*BarLength*/
      boost::none, /*BarSurface*/
      boost::none, /*BendingShapeCode*/
      boost::none /*BendingParameters*/
   );

   file.addEntity(rebar_type);

   // add the new definition to the project
   auto project = file.getSingle<typename Schema::IfcProject>();
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   if (rel_declares_instances->size() == 0)
   {
      typename Schema::IfcDefinitionSelect::list::ptr related_definitions(new typename Schema::IfcDefinitionSelect::list);
      related_definitions->push(rebar_type);

      auto rel_declares = new typename Schema::IfcRelDeclares(
         IfcParse::IfcGlobalId(),
         nullptr,
         boost::none,
         boost::none,
         project,
         related_definitions);

      file.addEntity(rel_declares);
   }
   else
   {
      for (auto& rel_declares : *rel_declares_instances)
      {
         if (rel_declares->RelatingContext()->as<typename Schema::IfcProject>())
         {
            auto related_definitions = rel_declares->RelatedDefinitions();
            related_definitions->push(rebar_type);
            rel_declares->setRelatedDefinitions(related_definitions);
            break;
         }
      }
   }

   return rebar_type;
}
