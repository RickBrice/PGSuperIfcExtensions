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
#include "resource.h"
#include "MappingTableWindows.h"
#include "MappingEditorViews.h" // ShowMappingHelp
#include "Documentation\MappingTableEditor\IfcMappingTableEditor.hh"
#include "MenuLogo.h"

#include <CommCtrl.h>
#pragma comment(lib, "comctl32.lib") // LoadIconWithScaleDown

namespace
{
   const TCHAR* const CANT_BE_USED = _T("This mapping table can't be used. See the message in the window, or choose another table.");

   CString configured_table_label()
   {
      return _T("The configured table: ") + GetMappingTableDisplayName(CIfcMappingTable::GetTableSetting());
   }
}

/////////////////////////////////////////////////////////////////////////////
// CMappingTableConfigPage

BEGIN_MESSAGE_MAP(CMappingTableConfigPage, CPropertyPage)
   ON_BN_CLICKED(IDC_MAPPING_STANDARD, &CMappingTableConfigPage::OnTableTypeChanged)
   ON_BN_CLICKED(IDC_MAPPING_AGENCY, &CMappingTableConfigPage::OnTableTypeChanged)
   ON_CBN_SELCHANGE(IDC_MAPPING_FILE, &CMappingTableConfigPage::OnRecentChanged)
   ON_BN_CLICKED(IDC_MAPPING_BROWSE, &CMappingTableConfigPage::OnBrowse)
   ON_WM_HELPINFO()
END_MESSAGE_MAP()

CMappingTableConfigPage::CMappingTableConfigPage() : CPropertyPage(IDD_MAPPING_CONFIG_PAGE)
{
   SetPageLogo(*this);
   m_psp.dwFlags |= PSP_USETITLE | PSP_USEHEADERTITLE | PSP_USEHEADERSUBTITLE;
   m_psp.dwFlags &= ~PSP_HASHELP;
   m_psp.pszTitle = _T("Configure IFC Extension");
   m_psp.pszHeaderTitle = _T("IFC Extension");
   m_psp.pszHeaderSubTitle = _T("Choose the mapping table that IFC imports and exports use");
}

CMappingTableConfigPage::~CMappingTableConfigPage()
{
   if (m_hLogo)
      ::DestroyIcon(m_hLogo);
}

BOOL CMappingTableConfigPage::OnInitDialog()
{
   CPropertyPage::OnInitDialog();

   // the bSI logo, from the multi-size icon at the size of its control (which follows the dialog font and DPI), so it
   // isn't stretched. The page runs in the application's module state, so the icon comes from the DLL's resources explicitly
   if (!m_hLogo)
   {
      CRect rect;
      GetDlgItem(IDC_MAPPING_LOGO)->GetClientRect(&rect);
      int size = std::min(rect.Width(), rect.Height());
      if (FAILED(::LoadIconWithScaleDown(AfxGetStaticModuleState()->m_hCurrentResourceHandle, MAKEINTRESOURCE(IDI_BSI), size, size, &m_hLogo)))
         m_hLogo = nullptr;
   }
   if (m_hLogo)
      SendDlgItemMessage(IDC_MAPPING_LOGO, STM_SETICON, (WPARAM)m_hLogo);

   m_Chooser.Init(this, _T("Standard table (installed with the IFC extension)"), CIfcMappingTable::STANDARD_TABLE,
      MappingTableSource::InstalledStandard, MappingTableSource::ConfigurationSetting);
   if (!m_Table.empty())
      m_Chooser.SelectOther(m_Table);

   return TRUE;
}

bool CMappingTableConfigPage::TakeChoice(bool bReport)
{
   if (!GetSafeHwnd())
      return true; // the page wasn't shown, m_Table is unchanged

   if (!m_Chooser.Check())
   {
      // no silent fallback: a table that can't be used isn't saved (G5)
      if (bReport)
         AfxMessageBox(CANT_BE_USED, MB_OK | MB_ICONEXCLAMATION);
      return false;
   }
   m_Table = m_Chooser.IsDefault() ? std::filesystem::path() : m_Chooser.GetOtherTable();
   return true;
}

// The Configure BridgeLink wizard tracks the page order, so its pages tell it when they change (as PGSuper's page does)
LRESULT CMappingTableConfigPage::OnWizardNext()
{
   if (!TakeChoice(true))
      return -1;

   AFX_MANAGE_STATE(AfxGetAppModuleState());
   dynamic_cast<IBridgeLinkConfigurationParent*>(GetParent())->GetNextPage();
   return CPropertyPage::OnWizardNext();
}

LRESULT CMappingTableConfigPage::OnWizardBack()
{
   TakeChoice(false); // a table that can't be used isn't kept, but doesn't stop going back

   AFX_MANAGE_STATE(AfxGetAppModuleState());
   dynamic_cast<IBridgeLinkConfigurationParent*>(GetParent())->GetBackPage();
   return CPropertyPage::OnWizardBack();
}

BOOL CMappingTableConfigPage::OnWizardFinish()
{
   return TakeChoice(true) ? CPropertyPage::OnWizardFinish() : FALSE;
}

BOOL CMappingTableConfigPage::OnHelpInfo(HELPINFO* pHelpInfo)
{
   ShowMappingHelp(IDH_MAPPING_CHOOSING_TABLE);
   return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CMappingTableConfigCallback

CPropertyPage* CMappingTableConfigCallback::CreatePropertyPage()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   auto* pPage = new CMappingTableConfigPage; // deleted by the wizard
   pPage->m_Table = CIfcMappingTable::GetTableSetting();
   return pPage;
}

void CMappingTableConfigCallback::OnOK(CPropertyPage* pPage)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CIfcMappingTable::SetTableSetting(static_cast<CMappingTableConfigPage*>(pPage)->m_Table);
}

/////////////////////////////////////////////////////////////////////////////
// CIfcImportTableDlg

BEGIN_MESSAGE_MAP(CIfcImportTableDlg, CDialog)
   ON_BN_CLICKED(IDC_MAPPING_STANDARD, &CIfcImportTableDlg::OnTableTypeChanged)
   ON_BN_CLICKED(IDC_MAPPING_AGENCY, &CIfcImportTableDlg::OnTableTypeChanged)
   ON_CBN_SELCHANGE(IDC_MAPPING_FILE, &CIfcImportTableDlg::OnRecentChanged)
   ON_BN_CLICKED(IDC_MAPPING_BROWSE, &CIfcImportTableDlg::OnBrowse)
   ON_WM_HELPINFO()
END_MESSAGE_MAP()

CIfcImportTableDlg::CIfcImportTableDlg(const CString& strModelFile, CWnd* pParent)
   : CDialog(IDD_IFC_IMPORT_TABLE, pParent), m_strModelFile(strModelFile)
{
}

BOOL CIfcImportTableDlg::OnInitDialog()
{
   CDialog::OnInitDialog();
   SetWindowLogo(this);

   SetDlgItemText(IDC_IMPORT_MODEL_FILE, m_strModelFile);
   // the default is the configured table (an empty path: LoadActive uses the setting, else the standard table)
   m_Chooser.Init(this, configured_table_label(), std::filesystem::path(), MappingTableSource::ConfigurationSetting, MappingTableSource::ChosenForImport);

   return TRUE;
}

void CIfcImportTableDlg::OnOK()
{
   if (!m_Chooser.Check())
   {
      AfxMessageBox(CANT_BE_USED, MB_OK | MB_ICONEXCLAMATION);
      return;
   }
   m_strTable = m_Chooser.GetTable().c_str();
   m_Source = m_Chooser.GetSource();
   CDialog::OnOK();
}

BOOL CIfcImportTableDlg::OnHelpInfo(HELPINFO* pHelpInfo)
{
   ShowMappingHelp(IDH_MAPPING_CHOOSING_TABLE);
   return TRUE;
}

/////////////////////////////////////////////////////////////////////////////
// CExportTablePage

BEGIN_MESSAGE_MAP(CExportTablePage, CPropertyPage)
   ON_BN_CLICKED(IDC_MAPPING_STANDARD, &CExportTablePage::OnTableTypeChanged)
   ON_BN_CLICKED(IDC_MAPPING_AGENCY, &CExportTablePage::OnTableTypeChanged)
   ON_CBN_SELCHANGE(IDC_MAPPING_FILE, &CExportTablePage::OnRecentChanged)
   ON_BN_CLICKED(IDC_MAPPING_BROWSE, &CExportTablePage::OnBrowse)
   ON_WM_HELPINFO()
END_MESSAGE_MAP()

CExportTablePage::CExportTablePage() : CPropertyPage(IDD_EXPORT_TABLE_PAGE)
{
   SetPageLogo(*this);
   m_DefaultLabel = configured_table_label();
}

void CExportTablePage::SetDefault(const CString& label, const std::filesystem::path& table, MappingTableSource source)
{
   m_DefaultLabel = label;
   m_DefaultTable = table;
   m_DefaultSource = source;
   m_strTable = table.c_str();
   m_Source = source;
}

BOOL CExportTablePage::OnInitDialog()
{
   CPropertyPage::OnInitDialog();
   m_Chooser.Init(this, m_DefaultLabel, m_DefaultTable, m_DefaultSource, MappingTableSource::ChosenForExport);
   return TRUE;
}

bool CExportTablePage::TakeChoice()
{
   if (!m_Chooser.Check())
   {
      AfxMessageBox(CANT_BE_USED, MB_OK | MB_ICONEXCLAMATION);
      return false;
   }
   m_strTable = m_Chooser.GetTable().c_str();
   m_Source = m_Chooser.GetSource();
   return true;
}

BOOL CExportTablePage::OnKillActive()
{
   return TakeChoice() ? CPropertyPage::OnKillActive() : FALSE;
}

void CExportTablePage::OnOK()
{
   TakeChoice(); // OnKillActive has checked it
   CPropertyPage::OnOK();
}

BOOL CExportTablePage::OnHelpInfo(HELPINFO* pHelpInfo)
{
   ShowMappingHelp(IDH_MAPPING_CHOOSING_TABLE);
   return TRUE;
}
