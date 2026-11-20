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
//
// The SDK is kept small, so that the IFC extension can change without changing it: extensions read an element's keys
// through IfcExport::Element (accessors, not the IFC extension's export context) and pass their own elements' keys in
// IfcExport::ElementKeys. The IFC extension's own export context and options aren't part of it

#include <any>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <PsgLib\Keys.h>

// Changes when a type or interface of the extension SDK (Include\IfcExport) changes, from the first release that
// includes the SDK on. The IFC extension refuses an export extension built with another version. Version 1 is
// provisional until it's released: changes before then don't change it
#define IFC_EXPORT_EXTENSION_API_VERSION 1

namespace WBFL { namespace EAF { class Broker; }; };
namespace WBFL { namespace Materials { class PsStrand; class Rebar; }; };
class CBearingData2;
struct BearingElevationDetails;

// An intermediate diaphragm the IFC export writes: the cast-in-place diaphragm of one bay of a span, made of the halves
// PGSuper gives the girders on either side (IBridge::GetCastInPlaceDiaphragms), from web face to web face
struct IntermediateDiaphragmDetails
{
   SpanIndexType span = INVALID_INDEX;
   GirderIndexType left_girder = INVALID_INDEX; // the bay is between this girder and the next
   IndexType index = INVALID_INDEX;             // the diaphragm's position in the girders' lists
   bool bGeometry = false;  // false: PGSuper has only the diaphragm's weight (the library rule's weight is input)
   Float64 height = 0;      // when bGeometry
   Float64 thickness = 0;   // along the girders, when bGeometry
   Float64 length = 0;      // web face to web face
   Float64 volume = 0;      // when bGeometry: the solid's volume (its ends follow the girders)
   Float64 mass = 0;        // the mass of PGSuper's diaphragm load for the bay
};

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
   ErectionTower,  // IfcBridgePart of a PGSplice erection tower (a temporary support)
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
   Diaphragm,      // IfcBeam.DIAPHRAGM, a pier diaphragm or an intermediate diaphragm
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
   Moment,      // N m (IfcTorqueMeasure)
   Volume       // m^3
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
   Mass,       // lb
   Volume,     // ft^3
   ComponentDim // e.g. in (PGSuper's component dimension unit)
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
   std::vector<TargetValue> list; // a list value (written as IfcPropertyListValue) when not empty

   ExportValue() = default;
   ExportValue(const TargetValue& v) : state(State::Value), value(v) {}
   ExportValue(Float64 v) : state(State::Value), value(v) {}
   ExportValue(Int64 v) : state(State::Value), value(v) {}
   ExportValue(const std::string& v) : state(State::Value), value(v) {}
   static ExportValue List(const std::vector<TargetValue>& values) { ExportValue v; v.state = State::Value; v.list = values; return v; }
   static ExportValue NoValue() { return ExportValue(); }
   static ExportValue Absent() { ExportValue v; v.state = State::Absent; return v; }
};

namespace IfcExport
{
   // The element whose target values are exported, as an export extension's targets see it. The IFC extension gives
   // the keys its element roles' targets need (e.g. Segment() for a girder, Pier() for a pier); an extension's own
   // elements have the ElementKeys its hooks gave. A key the element doesn't have is INVALID_INDEX or nullptr.
   // New keys are added at the end
   class Element
   {
   public:
      virtual std::shared_ptr<WBFL::EAF::Broker> Broker() const = 0;
      virtual const CSegmentKey& Segment() const = 0;            // girder, closure joint, tendon
      virtual PierIndexType Pier() const = 0;                    // pier, abutment, foundation
      virtual SpanIndexType Span() const = 0;                    // an extension's elements that belong to a span; an intermediate diaphragm
      virtual IndexType Index() const = 0;                       // the element's position in its list (e.g. a stirrup zone, an intermediate diaphragm)
      virtual const WBFL::Materials::PsStrand* Strand() const = 0; // tendon (and its material)
      virtual pgsTypes::StrandType StrandType() const = 0;       // tendon
      virtual StrandIndexType StrandIndex() const = 0;           // tendon
      virtual const WBFL::Materials::Rebar* Rebar() const = 0;   // rebar material
      virtual pgsTypes::TrafficBarrierOrientation BarrierSide() const = 0; // barrier
      virtual const CBearingData2* Bearing() const = 0;          // bearing
      virtual const BearingElevationDetails* BearingElevation() const = 0; // bearing: its location and elevations
      virtual const IntermediateDiaphragmDetails* IntermediateDiaphragm() const = 0; // diaphragm: an intermediate diaphragm (Pier() is INVALID_INDEX)
      virtual const std::any& Extension() const = 0;             // data an export extension passes to its own targets (ElementKeys::extension)

   protected:
      ~Element() = default;
   };

   // An export extension target's value for an element
   using Getter = std::function<ExportValue(const Element&)>;

   // The keys of an export extension's element, for the targets of its property sets (IServices::WritePropertySets)
   struct ElementKeys
   {
      CSegmentKey segment;                // the girder it belongs to
      PierIndexType pier = INVALID_INDEX; // the pier it belongs to
      SpanIndexType span = INVALID_INDEX; // the span it belongs to
      IndexType index = INVALID_INDEX;    // its position in its list (e.g. a stirrup zone)
      std::any extension;                 // anything else its targets need
   };
}
