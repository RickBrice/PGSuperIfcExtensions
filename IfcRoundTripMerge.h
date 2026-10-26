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

#include "IfcRoundTrip.h"

#include <filesystem>
#include <map>
#include <string>
#include <vector>

class CIfcMappingTable;

struct RoundTripMergeResult
{
   std::vector<std::string> log;                          // what the merge changed and added, and what it left alone
   std::map<std::string, std::string> global_ids;         // GlobalId in PGSuper's model -> GlobalId in the result, of the elements the source has
   IndexType nUpdated = 0;                                // source elements whose PGSuper values changed
   IndexType nPropertiesChanged = 0, nPropertiesAdded = 0;
   IndexType nAdded = 0;                                  // elements added from PGSuper's model (with their parts)
   bool bModelTendons = false, bModelBars = false;        // some girders keep the model's strands, bars
};

// Writes a copy of the source model with PGSuper's design in it (devdocs/RoundTripDesign.md, decided 2026-10-26):
// - The source's elements that PGSuper models (the registry's) keep their GlobalIds, geometry, names, and everything
//   related to them. Their declared properties are updated: a property with a target gets PGSuper's value, in the source
//   property's unit; a constant or a placeholder keeps the source's value and is added only where it's missing; the
//   table's classifications are added where missing. A shared property set, type, or material is copied for an element
//   whose values differ from the others'.
// - What PGSuper designs and the source doesn't have is added from PGSuper's model: the strands and bars of girders that
//   have none in the source (girders with strands or bars keep the source's), bearings, and pier diaphragms. They're placed
//   in the source's coordinates and units and related to the source's elements.
// - Everything else is kept as it is. Elements of PGSuper's model that the source doesn't have (e.g. a girder added in
//   PGSuper) aren't added yet; the log lists them.
// exported is PGSuper's model as the exporter built it (it's changed). Throws std::runtime_error if the source model can't
// be read or has another schema
RoundTripMergeResult MergeIntoSourceModel(ifcopenshell::file& exported, const IfcExportIdentities& identities, const CRoundTripRegistry& registry, const CIfcMappingTable& table,
   const std::filesystem::path& source_path, const std::filesystem::path& result_path);
