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

static PierIndexType get_pier_count(IfcParse::IfcFile& file)
{
   auto parts = file.instances_by_type<IfcSchema::IfcBridgePart>();
   PierIndexType nPiers = 0;
   for (auto& part : *parts)
   {
      if (part->PredefinedType().has_value() && (part->PredefinedType().get() == IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_ABUTMENT || part->PredefinedType().get() == IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_PIER))
         nPiers++;
   }
   return nPiers;
}
#include <string>
#include <regex>
#include <utility>
#include <cctype>
#include <stdexcept>

// Convert letter sequence (A, B, ..., Z, AA, AB, ...) to zero-based index
static int letterToIndex(const std::string& s) {
   int value = 0;
   for (char c : s) {
      if (!std::isalpha(c))
         throw std::runtime_error("Invalid letter sequence");

      value = value * 26 + (std::toupper(c) - 'A' + 1);
   }
   return value - 1; // zero-based
}

static CGirderKey girder_key_from_string(const std::string& input) {
   // input string is expected to be in one of these formats:
   // Span (number) Beam|Girder (number or alpha)
   // Beam|Girder (number or alpha) Span (number)
   // Examples
   // Span 1 Girder 2
   // Span 3 Beam C
   // Beam 2 Span 4
   // Girder AA Span 3

   // Case-insensitive regex
   std::regex pattern(
      R"((?:Span\s+(\d+)\s+(?:Beam|Girder)\s+([A-Za-z0-9]+))|(?:Beam|Girder)\s+([A-Za-z0-9]+)\s+Span\s+(\d+))",
      std::regex_constants::icase
   );


   std::smatch match;
   if (!std::regex_match(input, match, pattern)) {
      std::ostringstream os;
      os << "Beam designation: " << input << " is not valid";
      WBFL::System::Logger::Debug(os.str().c_str());
      return CGirderKey();
   }

   int span;
   std::string member;

   if (match[1].matched) {
      span = std::stoi(match[1].str());
      member = match[2].str();
   }
   else {
      member = match[3].str();
      span = std::stoi(match[4].str());
   }

   // Convert span to zero-based index
   SpanIndexType spanIndex = span - 1;

   // Convert girder/beam designation to zero-based index
   GirderIndexType girderIndex;
   if (std::isalpha(member[0])) {
      try {
         girderIndex = letterToIndex(member);
      }
      catch (...)
      {
         WBFL::System::Logger::Debug("Invalid beam designation");
         spanIndex = INVALID_INDEX;
         girderIndex = INVALID_INDEX;
      }
   }
   else {
      girderIndex = std::stoi(member) - 1;
   }

   return CGirderKey(spanIndex, girderIndex );
}

static CGirderKey get_girder_key(const IfcSchema::IfcBeam* beam)
{
   CGirderKey girder_key;
   auto design_location_number = GetProperty<IfcSchema, IfcSchema::IfcLabel>(beam, "Pset_PrecastConcreteElementGeneral", "DesignLocationNumber");
   if (design_location_number)
   {
      girder_key = girder_key_from_string(*design_location_number);
      if (girder_key == CGirderKey())
      {
         WBFL::System::Logger::Debug("DesignLocationNumber property not found in Pset_PrecastConcreteElementGeneral, or the property was not formatted as expected");
      }
   }
   else
   {
      WBFL::System::Logger::Debug("Pset_PrecastConcreteElementGeneral not found.");
   }

   if (girder_key == CGirderKey())
   {
      if (beam->Name())
      {
         girder_key = girder_key_from_string(*(beam->Name()));
         if (girder_key == CGirderKey())
         {
            WBFL::System::Logger::Debug("IfcBeam::Name was not formatted as expected");
         }
      }
      else
      {
         WBFL::System::Logger::Debug("IfcBeam::Name attribute not found");
      }
   }

   if (girder_key == CGirderKey())
   {
      WBFL::System::Logger::Debug("Could not determine girder key from IfcBeam");
   }

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


static IfcSchema::IfcBridgePart* GetBridgePart(IfcParse::IfcFile& file, IfcSchema::IfcBridgePartTypeEnum part_type)
{
   auto parts = file.instances_by_type<IfcSchema::IfcBridgePart>();
   for (auto part : *parts)
   {
      if (part->PredefinedType().has_value() && part->PredefinedType().get() == part_type)
      {
         return part;
      }
   }
   return nullptr;
}

static std::vector<IfcSchema::IfcBridgePart*> GetBridgeParts(IfcParse::IfcFile& file, IfcSchema::IfcBridgePartTypeEnum part_type)
{
   std::vector<IfcSchema::IfcBridgePart*> parts_found;
   auto parts = file.instances_by_type<IfcSchema::IfcBridgePart>();
   for (auto part : *parts)
   {
      if (part->PredefinedType().has_value() && part->PredefinedType().get() == part_type)
      {
         parts_found.push_back(part);
      }
   }
   return parts_found;
}

template <typename Schema>
typename Schema::IfcMaterial* GetMaterial(typename Schema::IfcObjectDefinition* objectdef)
{
   auto associations = objectdef->HasAssociations();
   if (associations)
   {
      for (auto rel : *associations)
      {
         auto rel_associates_material = rel->as<IfcSchema::IfcRelAssociatesMaterial>();
         if (rel_associates_material)
         {
            auto material = rel_associates_material->RelatingMaterial();
            return material->as<typename Schema::IfcMaterial>();
         }
      }
   }

   return nullptr;
}

static int GetBeamTypeCount(IfcParse::IfcFile& file)
{
   int count = 0;
   auto beam_types = file.instances_by_type<IfcSchema::IfcBeamType>();
   for (auto beam_type : *beam_types)
   {
      if (beam_type->PredefinedType() == IfcSchema::IfcBeamTypeEnum::IfcBeamType_BEAM)
         count++;
   }

   return count;
}
