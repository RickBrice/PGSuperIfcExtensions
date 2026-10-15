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

#include "IfcMappingTable.h"
#include <filesystem>
#include <vector>

// The name of a table file for display ("PennDOT"), from the table's "name". The standard table (empty path or
// CIfcMappingTable::STANDARD_TABLE) is the installed standard table. If the file can't be read, its file name
CString GetMappingTableDisplayName(const std::filesystem::path& table);

// Drives the controls that choose a mapping table, used by the Configure BridgeLink page, the IFC Import window,
// and the IFC export options:
// - IDC_MAPPING_STANDARD: radio button for the default table (its text is set by Init)
// - IDC_MAPPING_AGENCY: radio button for another table
// - IDC_MAPPING_FILE: drop-down list of the recently used tables (CIfcMappingTable::GetRecentTables), by name, after
//   the standard table (CIfcMappingTable::STANDARD_TABLE) unless it's the default
// - IDC_MAPPING_BROWSE: chooses a table file, which is added to the list
// - IDC_MAPPING_STATUS: read-only multi-line edit that shows the chosen table and the tables it extends, or why it can't be used
// The dialog routes the controls' notifications to OnTableTypeChanged, OnRecentChanged, and OnBrowse.
class CMappingTableChooser
{
public:
   // default_label: text of the default radio button, e.g. "The configured table: Standard"
   // default_table and default_source: the default table as LoadActive takes it (an empty path is the configured table)
   // other_source: how a table chosen with the other radio button is described in logs and messages
   void Init(CWnd* pParent, const CString& default_label, const std::filesystem::path& default_table, MappingTableSource default_source, MappingTableSource other_source);

   // Selects the default table, or another table
   void SelectDefault();
   void SelectOther(const std::filesystem::path& table);

   bool IsDefault() const;
   std::filesystem::path GetOtherTable() const; // empty if none is selected

   // The choice as CIfcMappingTable::LoadActive takes it
   std::filesystem::path GetTable() const;
   MappingTableSource GetSource() const;

   // Loads the chosen table and shows what it is, or why it can't be used. Returns true if it can be used
   bool Check();

   void OnTableTypeChanged();
   void OnRecentChanged();
   void OnBrowse();

private:
   CWnd* m_pParent = nullptr;
   std::filesystem::path m_DefaultTable;
   MappingTableSource m_DefaultSource = MappingTableSource::ConfigurationSetting;
   MappingTableSource m_OtherSource = MappingTableSource::ChosenForImport;
   std::vector<std::filesystem::path> m_Tables; // the list's tables, in list order

   void FillList(const std::filesystem::path& select);
   void UpdateControls();
};
