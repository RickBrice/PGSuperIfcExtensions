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

// Targets are the PGSuper data items that mapping tables bind to IFC locations.
// The list is fixed in code. Tables refer to targets by name. See devdocs/MappingTablesDesign.md

#include <string>
#include <string_view>
#include <variant>
#include <vector>

// The kind of PGSuper item a target belongs to. Mapping table element roles (e.g. "girder") use these names
enum class ElementKind
{
   Project,
   Site,
   Bridge,
   BridgePart,
   Pier,
   Foundation,
   Alignment,
   Referent,
   Girder,
   Deck,
   Haunch,
   Bearing,
   Barrier
};

// The kind of value a target holds. It sets the unit type and the IFC value types that can hold it
enum class ValueKind
{
   Stress,
   Length,
   Angle,
   Force,
   Ratio,
   Count,
   Boolean,
   Text,
   TextList // every location's value is collected (e.g. candidate girder type names)
};

// A target value. Numbers are in PGSuper system units (SI)
using TargetValue = std::variant<Float64, Int64, bool, std::string>;

struct TargetDef
{
   std::string_view name;  // e.g. "girder.fci"
   ElementKind element;    // the element the target belongs to
   ValueKind kind;
   std::string_view description; // for messages
};

// All targets
const std::vector<TargetDef>& GetTargetDefs();

// The target with the given name, or nullptr
const TargetDef* FindTargetDef(std::string_view name);

// Mapping table name of an element kind (e.g. "girder"), and the reverse. Returns false if the name isn't an element role
std::string_view GetElementRoleName(ElementKind kind);
bool GetElementKind(std::string_view role_name, ElementKind& kind);

// True if the value kind is a number with a unit (stress, length, angle, force)
bool HasUnit(ValueKind kind);

// True if the value kind is a number (with or without a unit)
bool IsNumeric(ValueKind kind);

// Name of a value kind, for messages
std::string_view GetValueKindName(ValueKind kind);

// Value as text, for messages. Numbers are in system units
std::string FormatTargetValue(const TargetValue& value);
