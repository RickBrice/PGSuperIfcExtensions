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

template <typename Schema>
void Add_usBridge_Classification(hierarchy_helper<Schema>& file)
{
   // we are using the usBridge bSDD for classifications
   auto classification = file.create<typename Schema::IfcClassification>().initialize(
      std::string("usBridge")/*Source*/,
      std::string("1") /*Edition*/,
      std::string("2026-04-20") /*EditionDate*/,
      std::string("usBridge"),
      std::nullopt /*Description*/,
      BSDD_URI /*Specification*/,
      std::nullopt /*ReferenceTokens*/);


   auto project = file.getSingle<typename Schema::IfcProject>();

   std::vector<typename Schema::IfcDefinitionSelect> projects;
   projects.push_back(project);

   auto rel_associates_classification = file.create<typename Schema::IfcRelAssociatesClassification>();
   rel_associates_classification.setGlobalId(ifcopenshell::global_id());
   rel_associates_classification.setRelatedObjects(projects);
   rel_associates_classification.setRelatingClassification(classification);
}

/// @brief Checks if an object has the specified classification
/// @tparam Schema 
/// @param object 
/// @param identifier 
/// @return 
template <typename Schema>
bool HasClassification(typename Schema::IfcObjectDefinition object, std::string identifier)
{
   auto associations = object.HasAssociations();
   for (auto& rel : associations)
   {
      auto rel_associates_classification = rel.template as<typename Schema::IfcRelAssociatesClassification>();
      if (rel_associates_classification)
      {
         auto classification_reference = rel_associates_classification.RelatingClassification().template as<typename Schema::IfcClassificationReference>();
         if (classification_reference && classification_reference.Identification().value_or("") == identifier)
            return true;
      }
   }

   return false;
}


template <typename Schema>
void AssociateClassification(hierarchy_helper<Schema>& file, typename Schema::IfcClassificationSelect classification, typename Schema::IfcDefinitionSelect related_object)
{
   auto associations = related_object.template as<typename Schema::IfcObjectDefinition>().HasAssociations();
   for (auto& association : associations)
   {
      auto rel_classification = association.template as<typename Schema::IfcRelAssociatesClassification>();
      if (rel_classification && rel_classification.RelatingClassification() == classification)
      {
         auto related_objects = rel_classification.RelatedObjects();
         related_objects.push_back(related_object);
         rel_classification.setRelatedObjects(related_objects);
         return;
      }
   }

   // if we get this far, there was not already a classicification association for this product, so we will create a new one
   std::vector<typename Schema::IfcDefinitionSelect> related_objects;
   related_objects.push_back(related_object);
   auto rel_classification = file.create<typename Schema::IfcRelAssociatesClassification>();
   rel_classification.setGlobalId(ifcopenshell::global_id());
   rel_classification.setRelatedObjects(related_objects);
   rel_classification.setRelatingClassification(classification);
}

template <typename Schema>
void Classify_ObjectDefinition(hierarchy_helper<Schema>& file, typename Schema::IfcObjectDefinition object, const std::string& name)
{
   auto classification = file.getSingle<typename Schema::IfcClassification>();

   std::string code("usBridge_");
   code += name;
   auto classification_reference = file.create<typename Schema::IfcClassificationReference>();
   classification_reference.setLocation(BSDD_URI + std::string("class/") + code);
   classification_reference.setIdentification(code);
   classification_reference.setName(name);
   classification_reference.setReferencedSource(classification);

   AssociateClassification<Schema>(file, classification_reference, object);
}

template <typename Schema>
void Classify_usBridge_BridgeProject(hierarchy_helper<Schema>& file, typename Schema::IfcProject project)
{
   Classify_ObjectDefinition(file, project, std::string("BridgeProject"));
}

template <typename Schema>
void Classify_usBridge_BridgeSite(hierarchy_helper<Schema>& file, typename Schema::IfcSite site)
{
   Classify_ObjectDefinition<Schema>(file, site, std::string("BridgeSite"));
}

template <typename Schema>
void Classify_usBridge_GirderBridge(hierarchy_helper<Schema>& file, typename Schema::IfcBridge bridge)
{
   Classify_ObjectDefinition<Schema>(file, bridge, std::string("GirderBridge"));
}

template <typename Schema>
void Classify_usBridge_Superstructure(hierarchy_helper<Schema>& file, typename Schema::IfcProduct superstructure)
{
   Classify_ObjectDefinition(file, superstructure, std::string("BridgeSuperstructure"));
}

template <typename Schema>
void Classify_usBridge_Substructure(hierarchy_helper<Schema>& file, typename Schema::IfcProduct substructure)
{
   Classify_ObjectDefinition(file, substructure, std::string("BridgeSubstructure"));
}

template <typename Schema>
void Classify_usBridge_Deck(hierarchy_helper<Schema>& file, typename Schema::IfcProduct deck)
{
   Classify_ObjectDefinition(file, deck, std::string("Deck"));
}

template <typename Schema>
void Classify_usBridge_Abutment(hierarchy_helper<Schema>& file, typename Schema::IfcProduct abutment)
{
   Classify_ObjectDefinition(file, abutment, std::string("Abutment"));
}

template <typename Schema>
void Classify_usBridge_Pier(hierarchy_helper<Schema>& file, typename Schema::IfcProduct pier)
{
   Classify_ObjectDefinition(file, pier, std::string("Pier"));
}

template <typename Schema>
void Classify_usBridge_Foundation(hierarchy_helper<Schema>& file, typename Schema::IfcProduct foundation)
{
   Classify_ObjectDefinition(file, foundation, std::string("Foundation"));
}

template <typename Schema>
void Classify_usBridge_Slab(hierarchy_helper<Schema>& file, typename Schema::IfcProduct slab)
{
   Classify_ObjectDefinition(file, slab, std::string("DeckSlab"));
}

template <typename Schema>
void Classify_usBridge_Barrier(hierarchy_helper<Schema>& file, typename Schema::IfcProduct barrier)
{
   Classify_ObjectDefinition(file, barrier, std::string("Barrier"));
}

template <typename Schema>
void Classify_usBridge_Girder(hierarchy_helper<Schema>& file, typename Schema::IfcProduct girder)
{
   Classify_ObjectDefinition(file, girder, std::string("Girder"));
}

template <typename Schema>
void Classify_usBridge_PrecastGirderElement(hierarchy_helper<Schema>& file, typename Schema::IfcProduct girder)
{
   Classify_ObjectDefinition(file, girder, std::string("GirderPrecastConcrete"));
}

template <typename Schema>
void Classify_usBridge_Tendon(hierarchy_helper<Schema>& file, typename Schema::IfcTendon tendon)
{
   Classify_ObjectDefinition(file, tendon, std::string("Tendon"));
}

template <typename Schema>
void Classify_usBridge_TendonBundle(hierarchy_helper<Schema>& file, typename Schema::IfcElementAssembly tendon_bundle)
{
   Classify_ObjectDefinition(file, tendon_bundle, std::string("TendonBundle"));
}

template <typename Schema>
void Classify_usBridge_ReinforcementCage(hierarchy_helper<Schema>& file, typename Schema::IfcProduct rebar_assembly)
{
   Classify_ObjectDefinition(file, rebar_assembly, std::string("ReinforcementCage"));
}

template <typename Schema>
void Classify_usBridge_ReinforcingBar(hierarchy_helper<Schema>& file, typename Schema::IfcReinforcingBar rebar)
{
   Classify_ObjectDefinition(file, rebar, std::string("ReinforcingBar"));
}

template <typename Schema>
void Classify_usBridge_ReinforcingBarType(hierarchy_helper<Schema>& file, typename Schema::IfcReinforcingBarType rebar_type)
{
   // usBridge does not define a classification for IfcReinforcingBarType.
   // Classification it is unclear if classification association is inherited from the type to the instance.
   //Classify_ObjectDefinition(file, rebar_type, std::string("usBridge_ReinforcingBarType"), std::string("ReinforcingBarType"));
}