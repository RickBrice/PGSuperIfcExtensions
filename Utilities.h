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

#include <string>
#include <numbers>

#include "IfcImporterException.h"
#include "Properties.h"

static CGirderKey girder_key_from_string(const std::string& s)
{
   GroupIndexType grpIdx = INVALID_INDEX;
   GirderIndexType gdrIdx = INVALID_INDEX;
   std::istringstream iss(s);
   std::string word;

   while (iss >> word)
   {
      if (word == "Span")
         iss >> grpIdx;
      else if (word == "Girder")
         iss >> gdrIdx;
      else if (word == ",")
      { // do nothing
      }
      else
      {
         return CGirderKey(); // values are INVALID_INDEX
         //USES_CONVERSION;
         //std::_tostringstream os;
         //os << _T("Unexpected DesignLocationNumber property in Pset_PrecastConcreteElementGeneral property set (") << A2T(s.c_str()) << _T(")");
         //IFC_THROW(os.str().c_str());
      }
   }

   return { grpIdx - 1,gdrIdx - 1 };
}

static CGirderKey get_girder_key(const IfcSchema::IfcBeam* beam)
{
   CGirderKey girder_key;
   auto design_location_number = GetProperty<IfcSchema, IfcSchema::IfcLabel>(beam, "Pset_PrecastConcreteElementGeneral", "DesignLocationNumber");
   if (design_location_number)
      girder_key = girder_key_from_string(*design_location_number);
   else if (beam->Name())
      girder_key = girder_key_from_string(*(beam->Name()));

   return girder_key;
}

template <typename T>
constexpr T deg2rad(T deg) {
   return deg * std::numbers::pi_v<T> / 180.;
}


inline std::string GetEntityType(IfcUtil::IfcBaseInterface* entity)
{
   return entity->declaration().name();
}


template <typename E>
E* GetType(IfcSchema::IfcObject* object)
{
   auto types = object->IsTypedBy();
   for (auto type : *types)
   {
      return type->RelatingType()->as<E>();
   }

   return nullptr;
}

/// @brief Returns the predefined type of an object.
/// The predefined type is taken from any associated IfcTypeObject.
/// If the object is not associated with IfcTypeObject, then its type is from its PredefinedType attribute
/// @tparam O object class
/// @tparam T type object class
/// @tparam E predefined type enumeration
/// @param object 
/// @return 
template <typename O, typename T, typename E>
boost::optional<typename E> GetPredefinedType(IfcSchema::IfcObject* object)
{
   // first check if the object is typed
   auto types = object->IsTypedBy();
   for (auto type : *types)
   {
      return type->RelatingType()->as<T>()->PredefinedType();
   }

   return object->as<O>()->PredefinedType();
}
