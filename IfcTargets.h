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

#include "Include\IfcExport\IfcExtensionTypes.h"

class CIfcExportOptions;
struct BarDefinition;

// The element whose target values are exported. The exporter fills in the keys the element's targets need. Export
// extensions see it through IfcExport::Element (IfcExtensionRegistry.cpp), so it can change without changing the SDK
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
   const BarDefinition* bar = nullptr; // Rebar: the bar of a bar type (see BarShape.h)
   pgsTypes::TrafficBarrierOrientation barrier_side = pgsTypes::tboLeft; // Barrier
   const CBearingData2* bearing = nullptr; // Bearing
   const BearingElevationDetails* bearing_elevation = nullptr; // Bearing: its location and elevations
   const IntermediateDiaphragmDetails* intermediate_diaphragm = nullptr; // Diaphragm: an intermediate diaphragm (pier is INVALID_INDEX)
   SpanIndexType span = INVALID_INDEX; // extension roles that belong to a span
   IndexType index = INVALID_INDEX;    // extension roles: the element's position in its list (e.g. a stirrup zone)
   std::any extension;                 // data an export extension passes to its own targets
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
   bool needs_analysis = false; // export: the value comes from analysis results. Left out (no value) when the export
                                // options say so, e.g. for PGSplice, whose time-step analysis is slow
   ExportGetter get; // export: the value for an element. Empty for import-only targets
   std::string_view provider; // the export extension that registered the target. Empty for the IFC extension's own targets
};

// All of the IFC extension's own targets (not those of export extensions)
const std::vector<TargetDef>& GetTargetDefs();

// The target with the given name, or nullptr. During an export with an export extension, the extension's targets are
// found too (see IfcExtensionRegistry.h)
const TargetDef* FindTargetDef(std::string_view name);

// Mapping table name of an element kind (e.g. "girder"), and the reverse. Returns false if the name isn't an element role
std::string_view GetElementRoleName(ElementKind kind);
bool GetElementKind(std::string_view role_name, ElementKind& kind);

// The element whose targets an element role uses (e.g. an abutment uses the pier targets)
ElementKind GetTargetElement(ElementKind role);

// True if the value kind is a number with a unit (stress, length, angle, force, area, mass, time, linear force, mass density, moment)
bool HasUnit(ValueKind kind);

// True if the value kind is a number (with or without a unit)
bool IsNumeric(ValueKind kind);

// Name of a value kind, for messages
std::string_view GetValueKindName(ValueKind kind);

// Value as text, for messages. Numbers are in system units
std::string FormatTargetValue(const TargetValue& value);
