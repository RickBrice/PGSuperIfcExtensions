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

#include <memory>
#include <string>

namespace WBFL { namespace EAF { class Broker; }; };
class CIfcMappingTable;

// The mapping tables a project remembers, saved with the project (CIfcExtensionAgent).
// A table is a file path, or CIfcMappingTable::STANDARD_TABLE for the installed standard table. Empty: none.
// An export of the project uses its last export table, else its import table, else the configured table,
// unless another table is chosen for the export.
// {6C2B2D31-0AF5-418D-B038-D865A1EE660F}
DEFINE_GUID(IID_IIfcProjectMappingTable,
   0x6c2b2d31, 0x0af5, 0x418d, 0xb0, 0x38, 0xd8, 0x65, 0xa1, 0xee, 0x66, 0x0f);
class IIfcProjectMappingTable
{
public:
   // The table the IFC import that made this project used
   virtual std::wstring GetImportMappingTable() const = 0;
   virtual void SetImportMappingTable(const std::wstring& table) = 0;

   // The table of the project's last IFC export
   virtual std::wstring GetExportMappingTable() const = 0;
   virtual void SetExportMappingTable(const std::wstring& table) = 0;
};

// The project's table as IIfcProjectMappingTable stores it: the path of the table that was selected (not the tables it
// extends), or CIfcMappingTable::STANDARD_TABLE
std::wstring TableChoiceOf(const CIfcMappingTable& table);

// The table the project's exports default to: its last export table, else its import table (*pbFromImport is then
// true). Empty: the configured table
std::wstring GetProjectExportTable(std::shared_ptr<WBFL::EAF::Broker> pBroker, bool* pbFromImport = nullptr);

// Remembers the table an IFC import or export used, and adds it to the recent tables. An export that changes the
// remembered table marks the project modified, so it's saved
void RememberImportTable(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcMappingTable& table);
void RememberExportTable(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcMappingTable& table);
