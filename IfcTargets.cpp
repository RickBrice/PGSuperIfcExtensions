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
#include "stdafx.h"
#include "IfcTargets.h"

#include <algorithm>
#include <sstream>

const std::vector<TargetDef>& GetTargetDefs()
{
   // Targets used by the importer (M1). Export targets are added with the export engine (M3)
   static const std::vector<TargetDef> targets{
      { "bridge.number_of_spans", ElementKind::Bridge, ValueKind::Count, "number of spans" },
      { "girder.designation", ElementKind::Girder, ValueKind::Text, "girder designation (e.g. \"Span 1, Girder 2\")" },
      { "girder.type_names", ElementKind::Girder, ValueKind::TextList, "girder type name" },
      { "girder.fc", ElementKind::Girder, ValueKind::Stress, "girder concrete strength, f'c" },
      { "girder.fci", ElementKind::Girder, ValueKind::Stress, "girder concrete strength at release, f'ci" },
      { "girder.assembly_place", ElementKind::Girder, ValueKind::Text, "girder assembly place (e.g. FACTORY)" },
      { "girder.casting_method", ElementKind::Girder, ValueKind::Text, "girder casting method (e.g. PRECAST)" },
      { "bearing.fixed_x", ElementKind::Bearing, ValueKind::Boolean, "bearing fixed along the girder" },
      { "bearing.fixed_y", ElementKind::Bearing, ValueKind::Boolean, "bearing fixed across the girder" },
   };
   return targets;
}

const TargetDef* FindTargetDef(std::string_view name)
{
   const auto& targets = GetTargetDefs();
   auto found = std::find_if(targets.begin(), targets.end(), [name](const auto& target) {return target.name == name; });
   return found == targets.end() ? nullptr : &(*found);
}

namespace
{
   const std::vector<std::pair<ElementKind, std::string_view>>& element_role_names()
   {
      static const std::vector<std::pair<ElementKind, std::string_view>> names{
         { ElementKind::Project, "project" },
         { ElementKind::Site, "site" },
         { ElementKind::Bridge, "bridge" },
         { ElementKind::BridgePart, "bridge_part" },
         { ElementKind::Pier, "pier" },
         { ElementKind::Foundation, "foundation" },
         { ElementKind::Alignment, "alignment" },
         { ElementKind::Referent, "referent" },
         { ElementKind::Girder, "girder" },
         { ElementKind::Deck, "deck" },
         { ElementKind::Haunch, "haunch" },
         { ElementKind::Bearing, "bearing" },
         { ElementKind::Barrier, "barrier" },
      };
      return names;
   }
}

std::string_view GetElementRoleName(ElementKind kind)
{
   for (const auto& [k, name] : element_role_names())
   {
      if (k == kind)
         return name;
   }
   ASSERT(false); // every kind should have a name
   return "";
}

bool GetElementKind(std::string_view role_name, ElementKind& kind)
{
   for (const auto& [k, name] : element_role_names())
   {
      if (name == role_name)
      {
         kind = k;
         return true;
      }
   }
   return false;
}

bool HasUnit(ValueKind kind)
{
   return kind == ValueKind::Stress || kind == ValueKind::Length || kind == ValueKind::Angle || kind == ValueKind::Force;
}

bool IsNumeric(ValueKind kind)
{
   return HasUnit(kind) || kind == ValueKind::Ratio || kind == ValueKind::Count;
}

std::string_view GetValueKindName(ValueKind kind)
{
   switch (kind)
   {
   case ValueKind::Stress: return "stress";
   case ValueKind::Length: return "length";
   case ValueKind::Angle: return "angle";
   case ValueKind::Force: return "force";
   case ValueKind::Ratio: return "ratio";
   case ValueKind::Count: return "count";
   case ValueKind::Boolean: return "boolean";
   case ValueKind::Text: return "text";
   case ValueKind::TextList: return "text list";
   }
   ASSERT(false);
   return "";
}

std::string FormatTargetValue(const TargetValue& value)
{
   std::ostringstream os;
   std::visit([&os](const auto& v)
      {
         using T = std::decay_t<decltype(v)>;
         if constexpr (std::is_same_v<T, bool>)
            os << (v ? "true" : "false");
         else if constexpr (std::is_same_v<T, std::string>)
            os << "\"" << v << "\"";
         else
            os << v;
      }, value);
   return os.str();
}
