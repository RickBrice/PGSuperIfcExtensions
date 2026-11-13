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
