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

#include "Utilities.h"

//#pragma Reminder("TODO - generalize the property enum methods and move to IfcHierarchyHelper")
// Need to cache the IfcPropertyEnumeration for lookup - it can be used multiple times by reference
// Need to have a getPropertyEnumeration method
// Need to generalize the enumValues from strings to IfcValue
// createPropertyEnumeratedValue needs two forms, a single value and a vector of values
template <typename Schema>
typename Schema::IfcPropertyEnumeration* createPropertyEnumeration(const std::string& name, std::vector<std::string>& enumValues, typename Schema::IfcUnit* unit = nullptr)
{
   typename Schema::IfcValue::list::ptr enum_values(new typename Schema::IfcValue::list);
   for (const auto& value : enumValues)
   {
      enum_values->push(new typename Schema::IfcLabel(value));
   }

   auto property_enum = new typename Schema::IfcPropertyEnumeration(name, enum_values, unit);
   return property_enum;
}

template <typename Schema>
typename Schema::IfcPropertyEnumeratedValue* createPropertyEnumeratedValue(const std::string& property_name, typename Schema::IfcPropertyEnumeration* enumeration, const std::string& value)
{
   typename Schema::IfcValue::list::ptr list_of_selected_enum_values(new typename Schema::IfcValue::list);
   list_of_selected_enum_values->push(new typename Schema::IfcLabel(value));
   auto property_enum_value = new typename Schema::IfcPropertyEnumeratedValue(property_name, boost::none, list_of_selected_enum_values, enumeration);
   return property_enum_value;
}


template <typename Schema>
typename Schema::IfcLabel* getPropertyEnumeratedValue(typename Schema::IfcPropertyEnumeratedValue* enum_value)
{
   if (enum_value)
   {
      auto list_of_selected_enum_values = enum_value->EnumerationValues();
      if (list_of_selected_enum_values)
      {
         auto ptr = *list_of_selected_enum_values;
         ASSERT(ptr->size() == 1); // only expecting one, but there could be more. This is a limitation of this function
         auto value = *(ptr->begin());
         return value->as<typename Schema::IfcLabel>();
      }
   }
   return nullptr;
}

template <typename Schema>
typename Schema::IfcPropertySet* GetPropertySet(typename Schema::IfcObject* object, std::string name)
{
   // First check for property sets on the object itself since they override
   // properties defined on the object type (if used).
   // See 5.1.3.6 IfcObject
   auto rel_defines_by_properties = object->IsDefinedBy();
   for (auto rel : *rel_defines_by_properties)
   {
      auto prop_set = rel->RelatingPropertyDefinition()->as<typename Schema::IfcPropertySet>();
      if(prop_set && prop_set->Name() == name)
      {
         return prop_set;
      }
   }

   // Now check the object types (if used)
   auto rel_defines_by_type = object->IsTypedBy();
   for (auto rel : *rel_defines_by_type)
   {
      auto relating_type = rel->RelatingType();
      auto property_set_definitions = relating_type->HasPropertySets().value_or(nullptr);

      if (property_set_definitions)
      {
         for (auto prop_set_definition : *property_set_definitions)
         {
            auto prop_set = prop_set_definition->as<typename Schema::IfcPropertySet>();
            if (prop_set && prop_set->Name() == name)
               return prop_set;
         }
      }
   }

   return nullptr;
}

template <typename Schema, typename T>
typename T* GetProperty(typename Schema::IfcObject* object, std::string pset_name, std::string property_name)
{
   auto pset = GetPropertySet<Schema>(object, pset_name);
   if (pset)
   {
      auto properties = pset->HasProperties();
      for (auto property : *properties)
      {
         if (property->Name() == property_name)
         {
            auto p = property->as<typename Schema::IfcPropertySingleValue>();
            if (p)
               return p->NominalValue()->as<T>();
            
            TRACE(GetEntityType(property).c_str());
         }
      }
   }
   return nullptr;
}

template <typename Schema,typename T>
std::vector<T*> GetPropertyList(typename Schema::IfcObject* object, std::string pset_name, std::string property_name)
{
   std::vector<T*> result;
   auto pset = GetPropertySet<Schema>(object, pset_name);
   if (pset)
   {
      auto properties = pset->HasProperties();
      for (auto property : *properties)
      {
         if (property->Name() == property_name)
         {
            auto list = *(property->as<typename Schema::IfcPropertyListValue>()->ListValues());
            for (auto value : *list)
            {
               result.push_back(value->as<T>());
            }
         }
      }
   }
   return result;
}

template <typename Schema, typename T>
typename T* GetPropertyEnum(typename Schema::IfcObject* object, std::string pset_name, std::string property_name)
{
   auto pset = GetPropertySet<Schema>(object, pset_name);
   if (pset)
   {
      auto properties = pset->HasProperties();
      for (auto property : *properties)
      {
         if (property->Name() == property_name)
         {
            auto enum_value = property->as<typename Schema::IfcPropertyEnumeratedValue>();
            if (enum_value)
               return getPropertyEnumeratedValue<Schema>(enum_value);

            TRACE(GetEntityType(property).c_str());
         }
      }
   }
   return nullptr;
}


template <typename Schema>
typename Schema::IfcMaterialProperties* GetMaterialPropertySet(typename Schema::IfcMaterialDefinition* matdef, std::string name)
{
   auto has_properties = matdef->HasProperties();
   for (auto material_properties : *has_properties)
   {
      if (material_properties->Name() == name)
      {
         return material_properties;
      }
   }

   return nullptr;
}

template <typename Schema, typename T>
typename T* GetMaterialProperty(typename Schema::IfcMaterialDefinition* matdef, std::string pset_name, std::string property_name)
{
   auto pset = GetMaterialPropertySet<Schema>(matdef, pset_name);
   if (pset)
   {
      auto properties = pset->Properties();
      for (auto property : *properties)
      {
         if (property->Name() == property_name)
         {
            auto p = property->as<typename Schema::IfcPropertySingleValue>();
            return p->NominalValue()->as<T>();
         }
      }
   }
   return nullptr;
}