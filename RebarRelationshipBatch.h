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

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "bSDD.h"

// RebarRelationshipBatch buffers reinforcement relationship membership while the
// bar/strand occurrences are authored and emits each relationship exactly once in
// Flush().
//
// It replaces the per-occurrence calls to hierarchy_helper::addRelatedObject(),
// AssociateMaterial(), Classify_usBridge_*() and AddPropertySet()/AddQto(). In
// IfcOpenShell 0.9 every one of those calls does a full instances_by_type<>() file
// scan and then copies + rewrites the target relationship's entire RelatedObjects
// list, so calling them once per bar is O(n^2) in bar count (a real PGSuper model
// has ~7000 bars). Accumulating here and writing each RelatedObjects list a single
// time makes it O(n).
template <typename Schema>
class RebarRelationshipBatch
{
public:
   // occurrence -> the IfcElementAssembly it belongs to (IfcRelAggregates)
   void Aggregate(typename Schema::IfcObjectDefinition assembly, typename Schema::IfcObjectDefinition occurrence)
   {
      auto& entry = m_Aggregates[assembly.id()];
      entry.first = assembly;
      entry.second.push_back(occurrence);
   }

   // occurrence -> its type, e.g. IfcReinforcingBarType / IfcTendonType (IfcRelDefinesByType)
   void Type(typename Schema::IfcObjectDefinition type, typename Schema::IfcObjectDefinition occurrence)
   {
      auto& entry = m_Types[type.id()];
      entry.first = type;
      entry.second.push_back(occurrence);
   }

   // occurrence -> its material (IfcRelAssociatesMaterial)
   void Material(typename Schema::IfcMaterial material, typename Schema::IfcObjectDefinition occurrence)
   {
      auto& entry = m_Materials[material.id()];
      entry.first = material;
      entry.second.push_back(occurrence);
   }

   // occurrence -> a shared property set or quantity set (IfcRelDefinesByProperties).
   // The caller creates the definition once (outside its bar loop) and registers
   // every occurrence against that one instance.
   void Properties(typename Schema::IfcPropertySetDefinition definition, typename Schema::IfcObjectDefinition occurrence)
   {
      if (!definition)
         return;

      auto& entry = m_Properties[definition.id()];
      entry.first = definition;
      entry.second.push_back(occurrence);
   }

   // occurrence -> a usBridge bSDD classification (IfcRelAssociatesClassification).
   // The shared IfcClassificationReference is built on first use of each name and
   // mirrors Classify_ObjectDefinition() in USBridge_Classifications.h.
   void Classify(hierarchy_helper<Schema>& file, const std::string& name, typename Schema::IfcObjectDefinition occurrence)
   {
      auto found = m_Classifications.find(name);
      if (found == m_Classifications.end())
      {
         std::string code("usBridge_");
         code += name;

         auto classification_reference = file.template create<typename Schema::IfcClassificationReference>();
         classification_reference.setLocation(BSDD_URI + std::string("class/") + code);
         classification_reference.setIdentification(code);
         classification_reference.setName(name);
         classification_reference.setReferencedSource(file.template getSingle<typename Schema::IfcClassification>());

         found = m_Classifications.emplace(name,
            std::make_pair(classification_reference, std::vector<typename Schema::IfcObjectDefinition>{})).first;
      }

      found->second.second.push_back(occurrence);
   }

   // Emit one relationship instance per relating object with its full member list.
   void Flush(hierarchy_helper<Schema>& file)
   {
      typename Schema::IfcOwnerHistory no_owner_history{};

      for (auto& [id, entry] : m_Aggregates)
      {
         file.template create<typename Schema::IfcRelAggregates>().initialize(
            ifcopenshell::global_id(), no_owner_history, std::nullopt, std::nullopt,
            entry.first, entry.second);
      }

      for (auto& [id, entry] : m_Types)
      {
         std::vector<typename Schema::IfcObject> related_objects;
         related_objects.reserve(entry.second.size());
         for (auto& occurrence : entry.second)
            related_objects.push_back(occurrence.template as<typename Schema::IfcObject>());

         file.template create<typename Schema::IfcRelDefinesByType>().initialize(
            ifcopenshell::global_id(), no_owner_history, std::nullopt, std::nullopt,
            related_objects, entry.first.template as<typename Schema::IfcTypeObject>());
      }

      for (auto& [id, entry] : m_Materials)
      {
         file.template create<typename Schema::IfcRelAssociatesMaterial>().initialize(
            ifcopenshell::global_id(), no_owner_history, std::nullopt, std::nullopt,
            ToDefinitionSelect(entry.second), entry.first);
      }

      for (auto& [name, entry] : m_Classifications)
      {
         file.template create<typename Schema::IfcRelAssociatesClassification>().initialize(
            ifcopenshell::global_id(), no_owner_history, std::nullopt, std::nullopt,
            ToDefinitionSelect(entry.second), entry.first);
      }

      for (auto& [id, entry] : m_Properties)
      {
         file.template create<typename Schema::IfcRelDefinesByProperties>().initialize(
            ifcopenshell::global_id(), no_owner_history, std::nullopt, std::nullopt,
            entry.second, entry.first);
      }

      m_Aggregates.clear();
      m_Types.clear();
      m_Materials.clear();
      m_Classifications.clear();
      m_Properties.clear();
   }

private:
   static std::vector<typename Schema::IfcDefinitionSelect> ToDefinitionSelect(const std::vector<typename Schema::IfcObjectDefinition>& objects)
   {
      std::vector<typename Schema::IfcDefinitionSelect> selects;
      selects.reserve(objects.size());
      for (auto& object : objects)
         selects.push_back(object);
      return selects;
   }

   // keyed by the relating instance's id (std::string name for classifications)
   std::map<uint32_t, std::pair<typename Schema::IfcObjectDefinition, std::vector<typename Schema::IfcObjectDefinition>>> m_Aggregates;
   std::map<uint32_t, std::pair<typename Schema::IfcObjectDefinition, std::vector<typename Schema::IfcObjectDefinition>>> m_Types;
   std::map<uint32_t, std::pair<typename Schema::IfcMaterial, std::vector<typename Schema::IfcObjectDefinition>>> m_Materials;
   std::map<std::string, std::pair<typename Schema::IfcClassificationReference, std::vector<typename Schema::IfcObjectDefinition>>> m_Classifications;
   std::map<uint32_t, std::pair<typename Schema::IfcPropertySetDefinition, std::vector<typename Schema::IfcObjectDefinition>>> m_Properties;
};
