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


template <typename Schema>
typename Schema::IfcPropertySet* GetPropertySet(typename Schema::IfcObject* object, std::string name)
{
   auto rels = object->IsDefinedBy();
   for (auto rel : *rels)
   {
      auto prop_set = rel->RelatingPropertyDefinition()->as<typename Schema::IfcPropertySet>();
      if (prop_set->Name() == name)
      {
         return prop_set;
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
            return p->NominalValue()->as<T>();
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