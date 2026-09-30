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

/*****************************************************************************
   Table-driven export of property sets (devdocs/MappingTablesDesign.md, M3)

   The exporter decides which elements play which element role and when their
   property sets are written. These functions write the property sets the
   mapping table declares for the role, with values from the targets' getters.

   CIfcExportSession loads the mapping table for one export and makes it
   available to the functions below, as CIfcImporter::GetTargetReader() does
   for the import.
*****************************************************************************/

#include "IfcMappingTable.h"
#include "IfcExporter.h"
#include "Properties.h"
#include "Units.h"
#include "USBridge_Classifications.h"

#include <EAF/EAFDisplayUnits.h>

class CIfcExportSession
{
public:
   // Loads the mapping table for an export: options.mapping_file, or the installed standard table.
   // Throws CIfcMappingTableException if the table can't be used
   CIfcExportSession(const CIfcExportOptions& options);
   ~CIfcExportSession();

   // The session of the export in progress
   static CIfcExportSession& Current();

   const CIfcMappingTable& GetTable() const { return *m_pTable; }

   // IfcPropertyEnumeration instances already written, by name and values, so each is written once
   std::map<std::string, uint32_t> enumerations;

   // IfcClassification instances written by WriteClassificationSystems, by name
   std::map<std::string, uint32_t> classification_systems;

   // IfcClassificationReference instances already written, by system and identification. Elements with the same
   // classification share a reference and its IfcRelAssociatesClassification
   std::map<std::string, uint32_t> classification_references;

private:
   std::unique_ptr<CIfcMappingTable> m_pTable;
   static CIfcExportSession* ms_pCurrent;
};

// True if the export options include the property set
bool IncludePropertySet(const PropertySetDeclaration& pset, const CIfcExportOptions& options);

// Value in the display unit, rounded to the target's increment
Float64 ConvertToDisplayUnits(std::shared_ptr<WBFL::EAF::Broker> pBroker, const TargetDef& target, Float64 value);

// True if properties are written in display units with their own unit (see CIfcExportOptions::display_units_for_properties)
bool UseDisplayUnits(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options);

namespace ifc_property_writer
{
   inline Float64 as_number(const TargetValue& value)
   {
      if (std::holds_alternative<Float64>(value)) return std::get<Float64>(value);
      if (std::holds_alternative<Int64>(value)) return (Float64)std::get<Int64>(value);
      if (std::holds_alternative<bool>(value)) return std::get<bool>(value) ? 1.0 : 0.0;
      return 0.0;
   }

   inline std::string as_text(const TargetValue& value)
   {
      if (std::holds_alternative<std::string>(value)) return std::get<std::string>(value);
      return FormatTargetValue(value);
   }

   template <typename Schema>
   typename Schema::IfcUnit display_unit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, ExportUnit unit)
   {
      switch (unit)
      {
      case ExportUnit::SpanLength: return GetSpanLengthUnit<Schema>(file, pBroker);
      case ExportUnit::Deflection: return GetDisplacementUnit<Schema>(file, pBroker);
      case ExportUnit::Stress: return GetStressUnit<Schema>(file, pBroker);
      case ExportUnit::Angle: return GetAngleUnit<Schema>(file, pBroker);
      case ExportUnit::SmallArea: return GetSmallAreaUnit<Schema>(file, pBroker);
      case ExportUnit::BigArea: return GetBigAreaUnit<Schema>(file, pBroker);
      case ExportUnit::Mass: return GetMassUnit<Schema>(file, pBroker);
      default: return {};
      }
   }

   // an IFC value of a type the table names (see IsExportValueType)
   template <typename Schema>
   typename Schema::IfcValue create_value(hierarchy_helper<Schema>& file, const std::string& type, const TargetValue& value)
   {
      Float64 v = as_number(value);
      if (type == "IFCLABEL") return file.template create<typename Schema::IfcLabel>().initialize(as_text(value));
      if (type == "IFCTEXT") return file.template create<typename Schema::IfcText>().initialize(as_text(value));
      if (type == "IFCIDENTIFIER") return file.template create<typename Schema::IfcIdentifier>().initialize(as_text(value));
      if (type == "IFCBOOLEAN") return file.template create<typename Schema::IfcBoolean>().initialize(v != 0.0);
      if (type == "IFCINTEGER") return file.template create<typename Schema::IfcInteger>().initialize((int64_t)std::llround(v));
      if (type == "IFCCOUNTMEASURE") return file.template create<typename Schema::IfcCountMeasure>().initialize((int64_t)std::llround(v));
      if (type == "IFCREAL") return file.template create<typename Schema::IfcReal>().initialize(v);
      if (type == "IFCPRESSUREMEASURE") return file.template create<typename Schema::IfcPressureMeasure>().initialize(v);
      if (type == "IFCLENGTHMEASURE") return file.template create<typename Schema::IfcLengthMeasure>().initialize(v);
      if (type == "IFCPOSITIVELENGTHMEASURE") return file.template create<typename Schema::IfcPositiveLengthMeasure>().initialize(v);
      if (type == "IFCNONNEGATIVELENGTHMEASURE") return file.template create<typename Schema::IfcNonNegativeLengthMeasure>().initialize(v);
      if (type == "IFCPLANEANGLEMEASURE") return file.template create<typename Schema::IfcPlaneAngleMeasure>().initialize(v);
      if (type == "IFCPOSITIVEPLANEANGLEMEASURE") return file.template create<typename Schema::IfcPositivePlaneAngleMeasure>().initialize(v);
      if (type == "IFCRATIOMEASURE") return file.template create<typename Schema::IfcRatioMeasure>().initialize(v);
      if (type == "IFCPOSITIVERATIOMEASURE") return file.template create<typename Schema::IfcPositiveRatioMeasure>().initialize(v);
      if (type == "IFCAREAMEASURE") return file.template create<typename Schema::IfcAreaMeasure>().initialize(v);
      if (type == "IFCVOLUMEMEASURE") return file.template create<typename Schema::IfcVolumeMeasure>().initialize(v);
      if (type == "IFCMASSMEASURE") return file.template create<typename Schema::IfcMassMeasure>().initialize(v);
      if (type == "IFCFORCEMEASURE") return file.template create<typename Schema::IfcForceMeasure>().initialize(v);
      ASSERT(false); // the table loader only accepts the types above
      return {};
   }

   template <typename Schema>
   typename Schema::IfcPropertyEnumeration enumeration(hierarchy_helper<Schema>& file, const PropertyDeclaration& property)
   {
      auto& enumerations = CIfcExportSession::Current().enumerations;
      std::string key = property.enumeration_name;
      for (const auto& v : property.enumeration_values)
         key += "|" + v;

      auto found = enumerations.find(key);
      if (found != enumerations.end())
         return file.instance_by_id(found->second).template as<typename Schema::IfcPropertyEnumeration>();

      auto values = property.enumeration_values;
      auto e = createPropertyEnumeration<Schema>(file, property.enumeration_name, values);
      enumerations.emplace(key, e.id());
      return e;
   }

   // The value of a property or quantity, from its constant or its target, in display or system units, and the unit of a
   // display unit value. Returns false if the target's value is absent
   template <typename Schema>
   bool get_value(hierarchy_helper<Schema>& file, const PropertyDeclaration& property, const ExportContext& context, std::optional<TargetValue>& value, typename Schema::IfcUnit& unit)
   {
      value = property.value;
      if (!value && property.target && property.target->get)
      {
         auto result = property.target->get(context);
         if (result.state == ExportValue::State::Absent)
            return false;

         if (result.state == ExportValue::State::Value)
         {
            value = result.value;
            if (std::holds_alternative<Float64>(*value) && property.target->display_unit != ExportUnit::None && UseDisplayUnits(context.broker, *context.options))
            {
               value = ConvertToDisplayUnits(context.broker, *property.target, std::get<Float64>(*value));
               unit = display_unit<Schema>(file, context.broker, property.target->display_unit);
            }
         }
      }
      return true;
   }

   // The property, or an empty property if the target's value is absent
   template <typename Schema>
   typename Schema::IfcProperty create_property(hierarchy_helper<Schema>& file, const PropertyDeclaration& property, const ExportContext& context)
   {
      std::optional<TargetValue> value;
      typename Schema::IfcUnit unit;
      if (!get_value<Schema>(file, property, context, value, unit))
         return {};

      std::optional<std::string> specification;
      if (!property.uri.empty())
         specification = property.uri;

      if (!property.enumeration_name.empty())
      {
         std::optional<std::vector<typename Schema::IfcValue>> selected;
         if (value)
            selected = std::vector<typename Schema::IfcValue>{ file.template create<typename Schema::IfcLabel>().initialize(as_text(*value)) };
         return file.template create<typename Schema::IfcPropertyEnumeratedValue>().initialize(property.name, specification, selected, enumeration<Schema>(file, property));
      }

      typename Schema::IfcValue ifc_value;
      if (value)
         ifc_value = create_value<Schema>(file, property.type, *value);

      return file.template create<typename Schema::IfcPropertySingleValue>().initialize(property.name, specification, ifc_value, unit);
   }

   template <typename Schema>
   std::vector<typename Schema::IfcProperty> create_properties(hierarchy_helper<Schema>& file, const PropertySetDeclaration& pset, const ExportContext& context)
   {
      std::vector<typename Schema::IfcProperty> properties;
      for (const auto& property : pset.properties)
      {
         if (auto p = create_property<Schema>(file, property, context))
            properties.push_back(p);
      }
      return properties;
   }

   // The quantity, or an empty quantity if the target's value is absent. A quantity without a value is written as 0,
   // because IfcPhysicalSimpleQuantity values aren't optional
   template <typename Schema>
   typename Schema::IfcPhysicalQuantity create_quantity(hierarchy_helper<Schema>& file, const PropertyDeclaration& quantity, const ExportContext& context)
   {
      std::optional<TargetValue> value;
      typename Schema::IfcUnit unit;
      if (!get_value<Schema>(file, quantity, context, value, unit))
         return {};

      Float64 v = value ? as_number(*value) : 0.0;
      auto named_unit = unit.template as<typename Schema::IfcNamedUnit>();
      const auto& type = quantity.type;
      if (type == "IFCQUANTITYLENGTH") return file.template create<typename Schema::IfcQuantityLength>().initialize(quantity.name, std::nullopt, named_unit, v, std::nullopt);
      if (type == "IFCQUANTITYAREA") return file.template create<typename Schema::IfcQuantityArea>().initialize(quantity.name, std::nullopt, named_unit, v, std::nullopt);
      if (type == "IFCQUANTITYVOLUME") return file.template create<typename Schema::IfcQuantityVolume>().initialize(quantity.name, std::nullopt, named_unit, v, std::nullopt);
      if (type == "IFCQUANTITYWEIGHT") return file.template create<typename Schema::IfcQuantityWeight>().initialize(quantity.name, std::nullopt, named_unit, v, std::nullopt);
      if (type == "IFCQUANTITYCOUNT") return file.template create<typename Schema::IfcQuantityCount>().initialize(quantity.name, std::nullopt, named_unit, (int64_t)std::llround(v), std::nullopt);
      ASSERT(false); // the table loader only accepts the types above
      return {};
   }

   template <typename Schema>
   typename Schema::IfcElementQuantity create_quantity_set(hierarchy_helper<Schema>& file, const PropertySetDeclaration& qto, const ExportContext& context)
   {
      std::vector<typename Schema::IfcPhysicalQuantity> quantities;
      for (const auto& quantity : qto.properties)
      {
         if (auto q = create_quantity<Schema>(file, quantity, context))
            quantities.push_back(q);
      }
      if (quantities.empty())
         return {};

      std::optional<std::string> description;
      if (!qto.uri.empty())
         description = qto.uri;

      std::optional<std::string> method;
      if (!qto.method.empty())
         method = qto.method;

      return file.template create<typename Schema::IfcElementQuantity>().initialize(ifcopenshell::global_id(), {}, qto.name, description, method, quantities);
   }

   template <typename Schema>
   typename Schema::IfcPropertySet create_property_set(hierarchy_helper<Schema>& file, const PropertySetDeclaration& pset, const ExportContext& context)
   {
      auto properties = create_properties<Schema>(file, pset, context);
      if (properties.empty())
         return {}; // IfcPropertySet.HasProperties needs at least one property

      std::optional<std::string> description;
      if (!pset.uri.empty())
         description = pset.uri;

      return file.template create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, pset.name, description, properties);
   }
}

// Writes the property sets the mapping table declares for the element role that attach to occurrences.
// Each property set is created once and related to all of the objects. With a condition, only the property sets
// with that condition are written (for elements whose property sets are written in more than one place)
template <typename Schema>
void WritePropertySets(hierarchy_helper<Schema>& file, ElementKind role, const std::vector<typename Schema::IfcObjectDefinition>& objects, const ExportContext& context, std::optional<PropertySetDeclaration::Condition> only = std::nullopt)
{
   for (const auto* pset : CIfcExportSession::Current().GetTable().GetPropertySets(role, PropertyOwner::Occurrence))
   {
      if (!IncludePropertySet(*pset, *context.options) || (only && pset->condition != *only) || pset->shared)
         continue;

      if (pset->quantities)
      {
         if (auto qto = ifc_property_writer::create_quantity_set<Schema>(file, *pset, context))
            file.template create<typename Schema::IfcRelDefinesByProperties>().initialize(ifcopenshell::global_id(), {}, std::nullopt, std::nullopt, objects, qto);
      }
      else if (auto property_set = ifc_property_writer::create_property_set<Schema>(file, *pset, context))
      {
         AddPropertySet(file, objects, property_set);
      }
   }
}

template <typename Schema>
void WritePropertySets(hierarchy_helper<Schema>& file, ElementKind role, typename Schema::IfcObjectDefinition object, const ExportContext& context, std::optional<PropertySetDeclaration::Condition> only = std::nullopt)
{
   std::vector<typename Schema::IfcObjectDefinition> objects{ object };
   WritePropertySets<Schema>(file, role, objects, context, only);
}

// The shared property sets and quantity sets the mapping table declares for the element role ("shared": true). They have no
// element-specific values, so they are created once and related to many elements, e.g. every bar through RebarRelationshipBatch.
// WritePropertySets leaves them out
template <typename Schema>
std::vector<typename Schema::IfcPropertySetDefinition> CreateSharedPropertySets(hierarchy_helper<Schema>& file, ElementKind role, const ExportContext& context)
{
   std::vector<typename Schema::IfcPropertySetDefinition> property_sets;
   for (const auto* pset : CIfcExportSession::Current().GetTable().GetPropertySets(role, PropertyOwner::Occurrence))
   {
      if (!pset->shared || !IncludePropertySet(*pset, *context.options))
         continue;

      if (pset->quantities)
      {
         if (auto qto = ifc_property_writer::create_quantity_set<Schema>(file, *pset, context))
            property_sets.push_back(qto);
      }
      else if (auto property_set = ifc_property_writer::create_property_set<Schema>(file, *pset, context))
      {
         property_sets.push_back(property_set);
      }
   }
   return property_sets;
}

// The property sets the mapping table declares for the element role that attach to type objects (IfcTypeObject.HasPropertySets)
template <typename Schema>
std::vector<typename Schema::IfcPropertySetDefinition> CreateTypePropertySets(hierarchy_helper<Schema>& file, ElementKind role, const ExportContext& context)
{
   std::vector<typename Schema::IfcPropertySetDefinition> property_sets;
   for (const auto* pset : CIfcExportSession::Current().GetTable().GetPropertySets(role, PropertyOwner::Type))
   {
      if (!IncludePropertySet(*pset, *context.options))
         continue;

      if (auto property_set = ifc_property_writer::create_property_set<Schema>(file, *pset, context))
         property_sets.push_back(property_set);
   }
   return property_sets;
}

// Writes the material properties the mapping table declares for the element role, for a material the role's element uses
template <typename Schema>
void WriteMaterialProperties(hierarchy_helper<Schema>& file, ElementKind role, typename Schema::IfcMaterial material, const ExportContext& context)
{
   for (const auto* pset : CIfcExportSession::Current().GetTable().GetPropertySets(role, PropertyOwner::Material))
   {
      if (!IncludePropertySet(*pset, *context.options))
         continue;

      auto properties = ifc_property_writer::create_properties<Schema>(file, *pset, context);
      if (properties.empty())
         continue;

      std::optional<std::string> description;
      if (!pset->uri.empty())
         description = pset->uri;

      file.template create<typename Schema::IfcMaterialProperties>().initialize(pset->name, description, properties, material);
   }
}

// Writes the classification systems (IfcClassification) the mapping table declares, associated with the project
template <typename Schema>
void WriteClassificationSystems(hierarchy_helper<Schema>& file)
{
   auto& session = CIfcExportSession::Current();
   auto project = file.template getSingle<typename Schema::IfcProject>();
   for (const auto* system : session.GetTable().GetClassificationSystems())
   {
      auto optional = [](const std::string& text) {return text.empty() ? std::optional<std::string>() : std::optional<std::string>(text); };
      auto classification = file.template create<typename Schema::IfcClassification>().initialize(
         optional(system->source), optional(system->edition), optional(system->edition_date), system->name,
         std::nullopt /*Description*/, optional(system->specification), std::nullopt /*ReferenceTokens*/);
      session.classification_systems[system->name] = classification.id();

      std::vector<typename Schema::IfcDefinitionSelect> projects{ project };
      auto rel = file.template create<typename Schema::IfcRelAssociatesClassification>();
      rel.setGlobalId(ifcopenshell::global_id());
      rel.setRelatedObjects(projects);
      rel.setRelatingClassification(classification);
   }
}

// The classification references the mapping table declares for an element role, created the first time they're needed.
// For elements whose classification is batched (RebarRelationshipBatch::Classify); other elements use Classify
template <typename Schema>
std::vector<typename Schema::IfcClassificationReference> GetClassificationReferences(hierarchy_helper<Schema>& file, ElementKind role, const CIfcExportOptions& options)
{
   std::vector<typename Schema::IfcClassificationReference> references;
   auto& session = CIfcExportSession::Current();
   for (const auto* classification : session.GetTable().GetClassifications(role))
   {
      PropertySetDeclaration condition;
      condition.condition = classification->condition;
      if (!IncludePropertySet(condition, options))
         continue;

      std::string key = classification->system + "|" + classification->identification;
      typename Schema::IfcClassificationReference reference;
      auto found = session.classification_references.find(key);
      if (found != session.classification_references.end())
      {
         reference = file.instance_by_id(found->second).template as<typename Schema::IfcClassificationReference>();
      }
      else
      {
         reference = file.template create<typename Schema::IfcClassificationReference>();
         if (!classification->location.empty())
            reference.setLocation(classification->location);
         reference.setIdentification(classification->identification);
         if (!classification->name.empty())
            reference.setName(classification->name);

         // a system that wasn't written in this export (e.g. a girder only export) leaves the source empty, as before
         auto system = session.classification_systems.find(classification->system);
         if (system != session.classification_systems.end())
            reference.setReferencedSource(file.instance_by_id(system->second).template as<typename Schema::IfcClassification>());

         session.classification_references.emplace(key, reference.id());
      }

      references.push_back(reference);
   }
   return references;
}

// Classifies an object with the classification references the mapping table declares for its element role
template <typename Schema>
void Classify(hierarchy_helper<Schema>& file, ElementKind role, typename Schema::IfcObjectDefinition object, const CIfcExportOptions& options)
{
   for (auto& reference : GetClassificationReferences<Schema>(file, role, options))
      AssociateClassification<Schema>(file, reference, object);
}
