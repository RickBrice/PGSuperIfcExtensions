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

// Types shared by the IFC extension and the plug-ins that extend its export (see IfcModelBuilder.h and
// devdocs/ExtensionsDesign.md). Targets are the PGSuper data items that mapping tables bind to IFC locations.

#include <any>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <PsgLib\Keys.h>

// Changes when a type or interface of the extension SDK (Include\IfcExport) changes. The IFC extension refuses an
// export extension built with another version
#define IFC_EXPORT_EXTENSION_API_VERSION 2 // 2: ValueKind::Moment

namespace WBFL { namespace EAF { class Broker; }; };
namespace WBFL { namespace Materials { class PsStrand; class Rebar; }; };
class CIfcExportOptions;

// The kind of PGSuper item a target belongs to. Mapping table element roles (e.g. "girder") use these names.
// Export extensions register their own element roles, numbered from FirstExtension
enum class ElementKind : int
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
   Tendon,            // IfcTendon, a prestressing strand
   TendonBundle,      // IfcElementAssembly of the strands of a segment, and its type
   Rebar,             // IfcReinforcingBar
   ReinforcementCage, // IfcElementAssembly of the reinforcing bars of a segment
   Deck,
   Haunch,
   Bearing,
   Barrier,

   FirstExtension = 1000 // element roles registered by export extensions
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
   TextList, // every location's value is collected (e.g. candidate girder type names)
   Time,        // s
   LinearForce, // N/m
   MassDensity, // kg/m^3
   Real,        // a number without a unit that isn't a ratio (e.g. a coefficient)
   Moment       // N m (IfcTorqueMeasure)
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
   CSegmentKey segment;                 // Girder, ClosureJoint, Tendon
   PierIndexType pier = INVALID_INDEX;  // Pier, Foundation
   const WBFL::Materials::PsStrand* strand = nullptr; // Tendon (and its material)
   pgsTypes::StrandType strand_type = pgsTypes::Straight; // Tendon
   StrandIndexType strand_index = INVALID_INDEX; // Tendon
   const WBFL::Materials::Rebar* rebar = nullptr; // Rebar material
   pgsTypes::TrafficBarrierOrientation barrier_side = pgsTypes::tboLeft; // Barrier
   SpanIndexType span = INVALID_INDEX; // extension roles that belong to a span
   IndexType index = INVALID_INDEX;    // extension roles: the element's position in its list (e.g. a stirrup zone)
   IndexType row = INVALID_INDEX;      // Rebar: the longitudinal reinforcement row of the bar
   std::any extension;                 // data an export extension passes to its own targets
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
   std::string_view provider; // the export extension that registered the target. Empty for the IFC extension's own targets
};
