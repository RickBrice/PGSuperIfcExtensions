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
#include "stdafx.h"
#include "IfcProjectMappingTable.h"
#include "IfcMappingTable.h"

#include <EAF\EAFUIIntegration.h>

std::wstring TableChoiceOf(const CIfcMappingTable& table)
{
   const auto& files = table.GetFiles();
   if (files.empty() || files.front().source == MappingTableSource::InstalledStandard)
      return CIfcMappingTable::STANDARD_TABLE;
   return files.front().path.wstring();
}

std::wstring GetProjectExportTable(std::shared_ptr<WBFL::EAF::Broker> pBroker, bool* pbFromImport)
{
   GET_IFACE2_NOCHECK(pBroker, IIfcProjectMappingTable, pProjectTable);
   std::wstring table = pProjectTable ? pProjectTable->GetExportMappingTable() : std::wstring();
   bool bFromImport = false;
   if (table.empty() && pProjectTable)
   {
      table = pProjectTable->GetImportMappingTable();
      bFromImport = !table.empty();
   }
   if (pbFromImport)
      *pbFromImport = bFromImport;
   return table;
}

void RememberImportTable(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcMappingTable& table)
{
   GET_IFACE2_NOCHECK(pBroker, IIfcProjectMappingTable, pProjectTable);
   if (!pProjectTable)
      return;

   auto choice = TableChoiceOf(table);
   pProjectTable->SetImportMappingTable(choice);
   CIfcMappingTable::AddRecentTable(choice);
}

void RememberExportTable(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcMappingTable& table)
{
   GET_IFACE2_NOCHECK(pBroker, IIfcProjectMappingTable, pProjectTable);
   if (!pProjectTable)
      return;

   auto choice = TableChoiceOf(table);
   CIfcMappingTable::AddRecentTable(choice);
   if (pProjectTable->GetExportMappingTable() != choice)
   {
      pProjectTable->SetExportMappingTable(choice);

      // the project must be saved for it to remember the table
      GET_IFACE2_NOCHECK(pBroker, IEAFDocument, pDoc);
      if (pDoc)
         pDoc->SetModified();
   }
}
