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

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <PsgLib\Keys.h>

namespace WBFL { namespace EAF { class Broker; }; };
class CIfcExportOptions;

// The kind of PGSuper item a target belongs to. Mapping table element roles (e.g. "girder") use these names
enum class ElementKind
{
   Project,
   Site,
   Bridge,
   Superstructure, // IfcBridgePart.SUPERSTRUCTURE
   Substructure,   // IfcBridgePart.SUBSTRUCTURE
   DeckPart,       // IfcBridgePart.DECK
   Pier,
   Abutment,       // an abutment is a pier in PGSuper: pier targets apply to it
   Foundation,
   Alignment,
   Referent,
   Girder,         // IfcBeam, one per segment
   GirderAssembly, // IfcElementAssembly.GIRDER of the segments of a spliced girder
   ClosureJoint,
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
   Area,
   Mass,
   Ratio,
   Count,
   Boolean,
   Text,
   TextList // every location's value is collected (e.g. candidate girder type names)
};

// A target value. Numbers are in PGSuper system units (SI)
using TargetValue = std::variant<Float64, Int64, bool, std::string>;

// The PGSuper display unit of an exported number, when properties are exported in display units
enum class ExportUnit
{
   None,       // written as it is, without a unit
   SpanLength, // e.g. ft
   Deflection, // e.g. in
   Stress,     // e.g. ksi
   Angle,      // degrees
   SmallArea,  // e.g. in^2
   BigArea,    // ft^2
   Mass        // lb
};

// The element whose target values are exported. The exporter fills in the keys the element's targets need
struct ExportContext
{
   std::shared_ptr<WBFL::EAF::Broker> broker;
   const CIfcExportOptions* options = nullptr;
   CSegmentKey segment;                 // Girder, ClosureJoint
   PierIndexType pier = INVALID_INDEX;  // Pier, Foundation
};

// The value a target's getter gives the exporter
struct ExportValue
{
   enum class State
   {
      Value,   // write the value
      NoValue, // write the property without a value
      Absent   // leave the property out
   };

   State state = State::NoValue;
   TargetValue value;

   ExportValue() = default;
   ExportValue(const TargetValue& v) : state(State::Value), value(v) {}
   ExportValue(Float64 v) : state(State::Value), value(v) {}
   ExportValue(Int64 v) : state(State::Value), value(v) {}
   ExportValue(const std::string& v) : state(State::Value), value(v) {}
   static ExportValue NoValue() { return ExportValue(); }
   static ExportValue Absent() { ExportValue v; v.state = State::Absent; return v; }
};

using ExportGetter = std::function<ExportValue(const ExportContext&)>;

struct TargetDef
{
   std::string_view name;  // e.g. "girder.fci"
   ElementKind element;    // the element the target belongs to
   ValueKind kind;
   std::string_view description; // for messages
   ExportUnit display_unit = ExportUnit::None; // export: the display unit of a number
   Float64 display_round = 0; // export: rounding increment in display units (0: none)
   ExportGetter get; // export: the value for an element. Empty for import-only targets
};

// All targets
const std::vector<TargetDef>& GetTargetDefs();

// The target with the given name, or nullptr
const TargetDef* FindTargetDef(std::string_view name);

// Mapping table name of an element kind (e.g. "girder"), and the reverse. Returns false if the name isn't an element role
std::string_view GetElementRoleName(ElementKind kind);
bool GetElementKind(std::string_view role_name, ElementKind& kind);

// The element whose targets an element role uses (e.g. an abutment uses the pier targets)
ElementKind GetTargetElement(ElementKind role);

// True if the value kind is a number with a unit (stress, length, angle, force, area, mass)
bool HasUnit(ValueKind kind);

// True if the value kind is a number (with or without a unit)
bool IsNumeric(ValueKind kind);

// Name of a value kind, for messages
std::string_view GetValueKindName(ValueKind kind);

// Value as text, for messages. Numbers are in system units
std::string FormatTargetValue(const TargetValue& value);
