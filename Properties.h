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

//#include "Utilities.h"

template <typename Schema>
void AddPropertySet(hierarchy_helper<Schema>& file, typename Schema::IfcObjectDefinition object, typename Schema::IfcPropertySet pset)
{
   if (!pset)
      return;

   std::vector<typename Schema::IfcObjectDefinition> related_objects;
   related_objects.push_back(object);

   AddPropertySet(file, related_objects, pset);
}

template <typename Schema>
void AddPropertySet(hierarchy_helper<Schema>& file, std::vector<typename Schema::IfcObjectDefinition> related_objects, typename Schema::IfcPropertySet pset)
{
   if (!pset)
      return;

   file.create<typename Schema::IfcRelDefinesByProperties>().initialize(ifcopenshell::global_id(), {}, std::nullopt, std::nullopt, related_objects, pset);
}

template <typename Schema>
void AddPropertySet(hierarchy_helper<Schema>& file, typename Schema::IfcTypeObject type, typename Schema::IfcPropertySet pset)
{
   if (!pset)
      return;

   auto has_property_sets = type.HasPropertySets().value_or(std::vector<typename Schema::IfcPropertySetDefinition>{});
   has_property_sets.push_back(pset.as<typename Schema::IfcPropertySetDefinition>());
   type.setHasPropertySets(has_property_sets);
}

//#pragma Reminder("TODO - generalize the property enum methods and move to IfcHierarchyHelper")
// Need to cache the IfcPropertyEnumeration for lookup - it can be used multiple times by reference
// Need to have a getPropertyEnumeration method
// Need to generalize the enumValues from strings to IfcValue
// createPropertyEnumeratedValue needs two forms, a single value and a vector of values
template <typename Schema>
typename Schema::IfcPropertyEnumeration createPropertyEnumeration(hierarchy_helper<Schema>& file, const std::string& name, std::vector<std::string>& enumValues, typename Schema::IfcUnit unit = {})
{
   std::vector<typename Schema::IfcValue> enum_values;
   for (const auto& value : enumValues)
   {
      enum_values.push_back(file.create<typename Schema::IfcLabel>().initialize(value));
   }

   auto property_enum = file.create<typename Schema::IfcPropertyEnumeration>().initialize(name, enum_values, unit);
   return property_enum;
}

template <typename Schema>
typename Schema::IfcPropertyEnumeratedValue createPropertyEnumeratedValue(hierarchy_helper<Schema>& file, const std::string& property_name, typename Schema::IfcPropertyEnumeration enumeration, const std::string& value)
{
   std::vector<typename Schema::IfcValue> list_of_selected_enum_values;
   list_of_selected_enum_values.push_back(file.create<typename Schema::IfcLabel>().initialize(value));
   auto property_enum_value = file.create<typename Schema::IfcPropertyEnumeratedValue>().initialize(property_name, std::nullopt, list_of_selected_enum_values, enumeration);
   return property_enum_value;
}


template <typename Schema>
typename Schema::IfcLabel getPropertyEnumeratedValue(typename Schema::IfcPropertyEnumeratedValue enum_value)
{
   if (enum_value)
   {
      auto list_of_selected_enum_values = enum_value.EnumerationValues();
      if (list_of_selected_enum_values)
      {
         auto& values = *list_of_selected_enum_values;
         ASSERT(values.size() == 1); // only expecting one, but there could be more. This is a limitation of this function
         auto value = values.front();
         return value.as<typename Schema::IfcLabel>();
      }
   }
   return {};
}

template <typename Schema>
typename Schema::IfcPropertySet GetPropertySet(typename Schema::IfcObject object, std::string name)
{
   // First check for property sets on the object itself since they override
   // properties defined on the object type (if used).
   // See 5.1.3.6 IfcObject
   auto rel_defines_by_properties = object.IsDefinedBy();
   for (auto& rel : rel_defines_by_properties)
   {
      auto prop_set = rel.RelatingPropertyDefinition().as<typename Schema::IfcPropertySet>();
      if (prop_set && prop_set.Name() == name)
      {
         return prop_set;
      }
   }

   // Now check the object types (if used)
   auto rel_defines_by_type = object.IsTypedBy();
   for (auto& rel : rel_defines_by_type)
   {
      auto relating_type = rel.RelatingType();
      auto property_set_definitions = relating_type.HasPropertySets();

      if (property_set_definitions)
      {
         for (auto& prop_set_definition : *property_set_definitions)
         {
            auto prop_set = prop_set_definition.as<typename Schema::IfcPropertySet>();
            if (prop_set && prop_set.Name() == name)
               return prop_set;
         }
      }
   }

   return {};
}

template <typename Schema, typename T>
std::optional<T> GetProperty(typename Schema::IfcObject object, std::string pset_name, std::string property_name)
{
   auto pset = GetPropertySet<Schema>(object, pset_name);
   if (pset)
   {
      auto properties = pset.HasProperties();
      for (auto& property : properties)
      {
         if (property.Name() == property_name)
         {
            auto p = property.as<typename Schema::IfcPropertySingleValue>();
            if (p)
               return p.NominalValue().as<T>();

            //TRACE(GetEntityType(property).c_str());
         }
      }
   }
   return std::nullopt;
}

template <typename Schema,typename T>
std::vector<T> GetPropertyList(typename Schema::IfcObject object, std::string pset_name, std::string property_name)
{
   std::vector<T> result;
   auto pset = GetPropertySet<Schema>(object, pset_name);
   if (pset)
   {
      auto properties = pset.HasProperties();
      for (auto& property : properties)
      {
         if (property.Name() == property_name)
         {
            auto list_value = property.as<typename Schema::IfcPropertyListValue>();
            if (list_value)
            {
               auto list = list_value.ListValues();
               if (list)
               {
                  for (auto& value : *list)
                  {
                     result.push_back(value.as<T>());
                  }
               }
            }
         }
      }
   }
   return result;
}

template <typename Schema, typename T>
std::optional<T> GetPropertyEnum(typename Schema::IfcObject object, std::string pset_name, std::string property_name)
{
   auto pset = GetPropertySet<Schema>(object, pset_name);
   if (pset)
   {
      auto properties = pset.HasProperties();
      for (auto& property : properties)
      {
         if (property.Name() == property_name)
         {
            auto enum_value = property.as<typename Schema::IfcPropertyEnumeratedValue>();
            if (enum_value)
               return getPropertyEnumeratedValue<Schema>(enum_value);

            TRACE(GetEntityType(property).c_str());
         }
      }
   }
   return std::nullopt;
}


template <typename Schema>
typename Schema::IfcMaterialProperties GetMaterialPropertySet(typename Schema::IfcMaterialDefinition matdef, std::string name)
{
   auto has_properties = matdef.HasProperties();
   for (auto& material_properties : has_properties)
   {
      if (material_properties.Name() == name)
      {
         return material_properties;
      }
   }

   return {};
}

template <typename Schema, typename T>
std::optional<T> GetMaterialProperty(typename Schema::IfcMaterialDefinition matdef, std::string pset_name, std::string property_name)
{
   auto pset = GetMaterialPropertySet<Schema>(matdef, pset_name);
   if (pset)
   {
      auto properties = pset.Properties();
      for (auto& property : properties)
      {
         if (property.Name() == property_name)
         {
            auto p = property.as<typename Schema::IfcPropertySingleValue>();
            if (p)
               return p.NominalValue().as<T>();
         }
      }
   }
   return std::nullopt;
}
