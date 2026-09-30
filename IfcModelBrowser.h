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
   An IFC model opened by the mapping table editor (M7 stage 4)

   The editor browses a model to pick the properties a table reads ("pick
   from a model"), and reads every target from the model with a table to show
   what an import would find ("try on a model"). It uses the import engine
   (CIfcTargetReader, CIfcTargetHints), so what it shows is what an import
   reads.
*****************************************************************************/

#include "IfcMappingTable.h"
#include "IfcImportUnits.h"

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class CIfcTargetReader;

// An element of the model that plays an element role
struct ModelElement
{
   int id = 0;
   std::string label; // e.g. "IfcBeam #123 'Span 1, Girder A'"
};

// A property of an element, as a location can refer to it
struct ModelProperty
{
   PropertyOwner owner = PropertyOwner::Occurrence; // where the property set is: the element, its type, or its material
   std::string pset;
   std::string name;
   std::string value;    // as text, e.g. "5.0" or "BTB45". A list or enumerated value shows all its values
   std::string ifc_type; // e.g. "IFCREAL", or the kind of property ("list", "enumeration") and the type of its values
};

class CIfcModel
{
public:
   // Opens an IFC 4.3 (IFC4X3_ADD2) model. Throws std::runtime_error with a message if it can't be read
   static std::unique_ptr<CIfcModel> Open(const std::filesystem::path& path);
   ~CIfcModel();

   const std::filesystem::path& GetPath() const { return m_Path; }
   ifcopenshell::file& GetFile() { return *m_pFile; }
   const CIfcImportUnits& GetUnits() const { return m_Units; }

   // The elements that play the role with the table, in file order. Empty if the table doesn't define the role
   std::vector<ModelElement> GetElements(const CIfcMappingTable& table, ElementKind role);

   // The properties of an element: its own property sets, its type's, and its material's
   std::vector<ModelProperty> GetProperties(int id);

   // The element with the id, or a null object
   IfcSchema::IfcObject GetObject(int id);

   // What the table reads for a target from an element, e.g. "5.000 ksi (34.47 MPa): "5.0" from IaDOT_PPCB... [ksi]",
   // or "not found (tried: ...)", with the importer's hints
   std::string DescribeReading(const CIfcMappingTable& table, std::string_view target, int id);

   // Reads every target from the elements of every element role with the table, and reports what an import would find
   std::string TryTable(const CIfcMappingTable& table);

private:
   CIfcModel() = default;
   std::filesystem::path m_Path;
   std::unique_ptr<ifcopenshell::file> m_pFile;
   CIfcImportUnits m_Units;
};

// A target value in SI and in a US display unit, e.g. "34473786 Pa (5.000 ksi)"
std::string FormatTargetValue(const TargetDef& target, const TargetValue& value);
