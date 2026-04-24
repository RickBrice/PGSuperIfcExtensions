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

template <typename Schema>
void Add_usBridge_Classification(IfcHierarchyHelper<Schema>& file)
{
   // we are using the usBridge bSDD for classifications
   auto classification = new typename Schema::IfcClassification(
      std::string("usBridge")/*Source*/,
      std::string("1") /*Edition*/,
      std::string("2026-04-20") /*EditionDate*/,
      std::string("usBridge"),
      boost::none /*Description*/,
      std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1") /*Specification*/,
      boost::none /*ReferenceTokens*/);
   file.addEntity(classification);

   auto project = file.getSingle<typename Schema::IfcProject>();

   typename Schema::IfcDefinitionSelect::list::ptr projects(new Schema::IfcDefinitionSelect::list);
   projects->push(project);

   auto rel_associates_classification = new typename Schema::IfcRelAssociatesClassification(
      IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, projects, classification);
   file.addEntity(rel_associates_classification);
}

/// @brief Checks if an object has the specified classification
/// @tparam Schema 
/// @param object 
/// @param identifier 
/// @return 
template <typename Schema>
bool HasClassification(typename Schema::IfcObjectDefinition* object, std::string identifier)
{
   auto associations = object->HasAssociations();
   if (associations)
   {
      for (auto rel : *associations)
      {
         auto rel_associates_classification = rel->as<typename Schema::IfcRelAssociatesClassification>();
         if (rel_associates_classification)
         {
            auto classification_reference = rel_associates_classification->RelatingClassification()->as<typename Schema::IfcClassificationReference>();
            if (classification_reference && classification_reference->Identification().value_or("") == identifier)
               return true;
         }
      }
   }

   return false;
}


template <typename Schema>
void AssociateClassification(IfcHierarchyHelper<Schema>& file, typename Schema::IfcClassificationSelect* classification, typename Schema::IfcDefinitionSelect* related_object)
{
   auto associations = related_object->as<typename Schema::IfcObjectDefinition>()->HasAssociations();
   if (associations)
   {
      for (auto association : *associations)
      {
         auto rel_classification = association->as<typename Schema::IfcRelAssociatesClassification>();
         if (rel_classification && rel_classification->RelatingClassification() == classification)
         {
            auto related_objects = rel_classification->RelatedObjects();
            related_objects->push(related_object);
            rel_classification->setRelatedObjects(related_objects);
            return;
         }
      }
   }

   // if we get this far, there was not already a classicification association for this product, so we will create a new one
   typename Schema::IfcDefinitionSelect::list::ptr related_objects(new typename Schema::IfcDefinitionSelect::list);
   related_objects->push(related_object);
   auto rel_classification = new typename Schema::IfcRelAssociatesClassification(
      IfcParse::IfcGlobalId(),
      nullptr,
      boost::none, // Name
      boost::none, // Description
      related_objects, // RelatedObjects
      classification // RelatingClassification
   );
   file.addEntity(rel_classification);
}

template <typename Schema>
void Classify_ObjectDefinition(IfcHierarchyHelper<Schema>& file, typename Schema::IfcObjectDefinition* object, const std::string& uri, const std::string& code, const std::string& name)
{
   auto classification = file.getSingle<typename Schema::IfcClassification>();

   auto classification_reference = new typename Schema::IfcClassificationReference(
      uri, /*Class identifier (uri) = Location*/
      code, /*Class code = Identification*/
      name,/*Class name = name*/
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);
   file.addEntity(classification_reference);

   AssociateClassification<Schema>(file, classification_reference, object);
}

template <typename Schema>
void Classify_usBridge_BridgeProject(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProject* project)
{
   Classify_ObjectDefinition(file, project, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_BridgeProject"), std::string("usBridge_BridgeProject"), std::string("BridgeProject"));
}

template <typename Schema>
void Classify_usBridge_BridgeSite(IfcHierarchyHelper<Schema>& file, typename Schema::IfcSite* site)
{
   Classify_ObjectDefinition<Schema>(file, site, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_BridgeSite"), std::string("usBridge_BridgeSite"), std::string("BridgeSite"));
}

template <typename Schema>
void Classify_usBridge_GirderBridge(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
{
   Classify_ObjectDefinition<Schema>(file, bridge, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_GirderBridge"), std::string("usBridge_GirderBridge"), std::string("GirderBridge"));
}

template <typename Schema>
void Classify_usBridge_Superstructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* superstructure)
{
   Classify_ObjectDefinition(file, superstructure, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_BridgeSuperstructure"), std::string("usBridge_BridgeSuperstructure"), std::string("BridgeSuperstructure"));
}

template <typename Schema>
void Classify_usBridge_Substructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* substructure)
{
   Classify_ObjectDefinition(file, substructure, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_BridgeSubstructure"), std::string("usBridge_BridgeSubstructure"), std::string("BridgeSubstructure"));
}

template <typename Schema>
void Classify_usBridge_Deck(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* deck)
{
   Classify_ObjectDefinition(file, deck, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Deck"), std::string("usBridge_Deck"), std::string("Deck"));
}

template <typename Schema>
void Classify_usBridge_Abutment(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* abutment)
{
   Classify_ObjectDefinition(file, abutment, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Abutment"), std::string("usBridge_Abutment"), std::string("Abutment"));
}

template <typename Schema>
void Classify_usBridge_Pier(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* pier)
{
   Classify_ObjectDefinition(file, pier, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Pier"), std::string("usBridge_Pier"), std::string("Pier"));
}

template <typename Schema>
void Classify_usBridge_Foundation(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* foundation)
{
   Classify_ObjectDefinition(file, foundation, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Foundation"), std::string("usBridge_Foundation"), std::string("Foundation"));
}

template <typename Schema>
void Classify_usBridge_Slab(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* slab)
{
   Classify_ObjectDefinition(file, slab, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_DeckSlab"), std::string("usBridge_DeckSlab"), std::string("DeckSlab"));
}

template <typename Schema>
void Classify_usBridge_Barrier(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* barrier)
{
   Classify_ObjectDefinition(file, barrier, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Barrier"), std::string("usBridge_Barrier"), std::string("Barrier"));
}

template <typename Schema>
void Classify_usBridge_Girder(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* girder)
{
   Classify_ObjectDefinition(file, girder, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Girder"), std::string("usBridge_Girder"), std::string("Girder"));
}

template <typename Schema>
void Classify_usBridge_PrecastGirderElement(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* girder)
{
   Classify_ObjectDefinition(file, girder, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_GirderPrecastConcrete"), std::string("usBridge_GirderPrecastConcrete"), std::string("GirderPrecastConcrete"));
}

template <typename Schema>
void Classify_usBridge_Tendon(IfcHierarchyHelper<Schema>& file, typename Schema::IfcTendon* tendon)
{
   Classify_ObjectDefinition(file, tendon, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_Tendon"), std::string("usBridge_Tendon"), std::string("Tendon"));
}

template <typename Schema>
void Classify_usBridge_ReinforcementCage(IfcHierarchyHelper<Schema>& file, typename Schema::IfcProduct* rebar_assembly)
{
   Classify_ObjectDefinition(file, rebar_assembly, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_ReinforcementCage"), std::string("usBridge_ReinforcementCage"), std::string("ReinforcementCage"));
}

template <typename Schema>
void Classify_usBridge_ReinforcingBar(IfcHierarchyHelper<Schema>& file, typename Schema::IfcReinforcingBar* rebar)
{
   Classify_ObjectDefinition(file, rebar, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_ReinforcingBar"), std::string("usBridge_ReinforcingBar"), std::string("ReinforcingBar"));
}

template <typename Schema>
void Classify_usBridge_ReinforcingBarType(IfcHierarchyHelper<Schema>& file, typename Schema::IfcReinforcingBarType* rebar_type)
{
   // usBridge does not define a classification for IfcReinforcingBarType.
   // Classification it is unclear if classification association is inherited from the type to the instance.
   //Classify_ObjectDefinition(file, rebar_type, std::string("https://identifier.buildingsmart.org/uri/aashto/usBridge/1/class/usBridge_ReinforcingBar"), std::string("usBridge_ReinforcingBar"), std::string("ReinforcingBar"));
}