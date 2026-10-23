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

// The windows that choose the mapping table IFC imports and exports use (devdocs/MappingTablesDesign.md, G5):
// - the Configure BridgeLink wizard page for the configured table
// - the IFC Import window, which can choose another table for one import
// - the Mapping Table tab of the IFC export options, which can choose another table for one export
// A table is loaded and validated when it's chosen, and a table that can't be used can't be chosen (no silent fallback).

#include "MappingTableChooser.h"
#include <BridgeLinkConfiguration.h>

// The Configure BridgeLink wizard page
class CMappingTableConfigPage : public CPropertyPage
{
public:
   CMappingTableConfigPage();
   ~CMappingTableConfigPage();

   std::filesystem::path m_Table; // in: the configured table; out: the chosen table. Empty: the installed standard table

protected:
   BOOL OnInitDialog() override;
   LRESULT OnWizardNext() override;
   LRESULT OnWizardBack() override;
   BOOL OnWizardFinish() override;

   afx_msg void OnTableTypeChanged() { m_Chooser.OnTableTypeChanged(); }
   afx_msg void OnRecentChanged() { m_Chooser.OnRecentChanged(); }
   afx_msg void OnBrowse() { m_Chooser.OnBrowse(); }
   afx_msg BOOL OnHelpInfo(HELPINFO* pHelpInfo);
   DECLARE_MESSAGE_MAP()

private:
   CMappingTableChooser m_Chooser;
   HICON m_hLogo = nullptr; // the bSI logo (IDI_BSI, at the size of its control)
   bool TakeChoice(bool bReport); // validates the choice and puts it in m_Table
};

// Adds the page to the Configure BridgeLink wizard and saves the configured table when the wizard finishes
class CMappingTableConfigCallback : public IBridgeLinkConfigurationCallback
{
public:
   CPropertyPage* CreatePropertyPage() override;
   void OnOK(CPropertyPage* pPage) override;
};

// The IFC Import window, after the IFC file is chosen
class CIfcImportTableDlg : public CDialog
{
public:
   CIfcImportTableDlg(const CString& strModelFile, CWnd* pParent = nullptr);

   // the choice as CIfcImportOptions::mapping_file and mapping_source take it
   CString m_strTable;
   MappingTableSource m_Source = MappingTableSource::ConfigurationSetting;

protected:
   BOOL OnInitDialog() override;
   void OnOK() override;

   afx_msg void OnTableTypeChanged() { m_Chooser.OnTableTypeChanged(); }
   afx_msg void OnRecentChanged() { m_Chooser.OnRecentChanged(); }
   afx_msg void OnBrowse() { m_Chooser.OnBrowse(); }
   afx_msg BOOL OnHelpInfo(HELPINFO* pHelpInfo);
   DECLARE_MESSAGE_MAP()

private:
   CString m_strModelFile;
   CMappingTableChooser m_Chooser;
};

// The Mapping Table tab of the IFC export options
class CExportTablePage : public CPropertyPage
{
public:
   CExportTablePage();

   // The default table (the project's, else the configured table), set before the options are shown
   void SetDefault(const CString& label, const std::filesystem::path& table, MappingTableSource source);

   // the choice as CIfcExportOptions::mapping_file and mapping_source take it (the default if the tab isn't visited)
   CString m_strTable;
   MappingTableSource m_Source = MappingTableSource::ConfigurationSetting;

protected:
   BOOL OnInitDialog() override;
   BOOL OnKillActive() override;
   void OnOK() override;

   afx_msg void OnTableTypeChanged() { m_Chooser.OnTableTypeChanged(); }
   afx_msg void OnRecentChanged() { m_Chooser.OnRecentChanged(); }
   afx_msg void OnBrowse() { m_Chooser.OnBrowse(); }
   afx_msg BOOL OnHelpInfo(HELPINFO* pHelpInfo);
   DECLARE_MESSAGE_MAP()

private:
   CString m_DefaultLabel;
   std::filesystem::path m_DefaultTable;
   MappingTableSource m_DefaultSource = MappingTableSource::ConfigurationSetting;
   CMappingTableChooser m_Chooser;
   bool TakeChoice();
};
