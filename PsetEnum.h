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

//#pragma Reminder("TODO - generalize the property enum methods and move to IfcHierarchyHelper")
// Need to cache the IfcPropertyEnumeration for lookup - it can be used multiple times by reference
// Need to have a getPropertyEnumeration method
// Need to generalize the enumValues from strings to IfcValue
// createPropertyEnumeratedValue needs two forms, a single value and a vector of values
template <typename Schema>
typename Schema::IfcPropertyEnumeration* createPropertyEnumeration(const std::string& name, std::vector<std::string>& enumValues, typename Schema::IfcUnit* unit = nullptr)
{
   typename aggregate_of<typename Schema::IfcValue>::ptr enum_values(new aggregate_of<typename Schema::IfcValue>());
   for (const auto& value : enumValues)
   {
      enum_values->push(new Schema::IfcLabel(value));
   }

   auto property_enum = new Schema::IfcPropertyEnumeration(name, enum_values, unit);
   return property_enum;
}

template <typename Schema>
typename Schema::IfcPropertyEnumeratedValue* createPropertyEnumeratedValue(const std::string& property_name, typename Schema::IfcPropertyEnumeration* enumeration, const std::string& value)
{
   typename aggregate_of<typename Schema::IfcValue>::ptr list_of_selected_enum_values(new aggregate_of<typename Schema::IfcValue>());
   list_of_selected_enum_values->push(new Schema::IfcLabel(value));
   auto property_enum_value = new Schema::IfcPropertyEnumeratedValue(property_name, boost::none, list_of_selected_enum_values, enumeration);
   return property_enum_value;
}
