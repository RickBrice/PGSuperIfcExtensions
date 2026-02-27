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
inline std::string GetEntityType(IfcUtil::IfcBaseInterface* entity)
{
   return entity->declaration().name();
}


template <typename E>
E* GetType(Ifc4x3_add2::IfcObject* object)
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
boost::optional<typename E> GetPredefinedType(Ifc4x3_add2::IfcObject* object)
{
   // first check if the object is typed
   auto types = object->IsTypedBy();
   for (auto type : *types)
   {
      return type->RelatingType()->as<T>()->PredefinedType();
   }

   return object->as<O>()->PredefinedType();
}
