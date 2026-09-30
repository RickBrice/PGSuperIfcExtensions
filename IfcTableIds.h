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

/*****************************************************************************
   Mapping tables and general-purpose IDS (devdocs/MappingTablesDesign.md, M5)

   WriteTableAsIds writes what a mapping table exports as a general IDS: one
   specification per element role, applicable to the role's elements, with a
   requirement for each exported property, quantity, and classification. It
   has no design values: properties with a target are required, constants
   are required with their value, and placeholders are optional.

   GenerateTableFromIds does the reverse: a mapping table from a general IDS
   (e.g. an agency IDS). A facet's target comes from the PGSuper instructions
   the writer puts on it, else the binding file, else the standard table's
   property with the same property set and name. Facets without a target are
   reported, with a binding entry to fill in. Manual assignments live in the
   binding file, so regenerating the table keeps them. The agency IDS itself
   is never modified (A4).

   The design-value IDS (IdsExporter.h) is not an input: it pins values to
   single elements of one model.
*****************************************************************************/

#include "IfcMappingTable.h"

#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

// Writes the table as a general IDS. notes gets what the IDS can't express (e.g. material property sets).
// Throws std::runtime_error if the IDS can't be written
void WriteTableAsIds(const CIfcMappingTable& table, std::ostream& os, std::vector<std::string>& notes);

struct IdsToTableResult
{
   std::string table_json;          // the generated mapping table, in the style of Standard.json (FormatMappingTable)
   std::vector<std::string> report; // specifications and facets without an element role or target (with binding entries), unsupported constructs
   size_t bound = 0;                // property facets with a target
   size_t unbound = 0;              // property facets without a target
};

// Generates a mapping table from a general IDS. binding may be empty.
// bExtendStandard: the table extends the standard table and has only what differs from it: element roles with another
// selector, new property sets and quantity sets, standard ones with more properties (merged with the standard ones),
// and new classifications. Otherwise the table stands alone with everything the IDS has (e.g. to check a round trip).
// Throws std::runtime_error with a message if the IDS or the binding file can't be read
IdsToTableResult GenerateTableFromIds(const std::filesystem::path& ids, const std::filesystem::path& binding, const CIfcMappingTable& standard, const std::string& table_name, bool bExtendStandard = true);
