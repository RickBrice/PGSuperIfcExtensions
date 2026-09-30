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

// MappingTableDlg.cpp : implementation file
//

#include "stdafx.h"
#include "resource.h"
#include "MappingTableDlg.h"
#include "IfcMappingTable.h"

#include <sstream>

IMPLEMENT_DYNAMIC(CMappingTableDlg, CDialog)

CMappingTableDlg::CMappingTableDlg(CWnd* pParent /*=nullptr*/)
   : CDialog(IDD_MAPPING_TABLE, pParent), m_TableType(0)
{
}

void CMappingTableDlg::DoDataExchange(CDataExchange* pDX)
{
   CDialog::DoDataExchange(pDX);
   DDX_Radio(pDX, IDC_MAPPING_STANDARD, m_TableType);
   DDX_Text(pDX, IDC_MAPPING_FILE, m_strTable);
}

BEGIN_MESSAGE_MAP(CMappingTableDlg, CDialog)
   ON_BN_CLICKED(IDC_MAPPING_STANDARD, &CMappingTableDlg::OnTableTypeChanged)
   ON_BN_CLICKED(IDC_MAPPING_AGENCY, &CMappingTableDlg::OnTableTypeChanged)
   ON_BN_CLICKED(IDC_MAPPING_BROWSE, &CMappingTableDlg::OnBrowse)
   ON_EN_KILLFOCUS(IDC_MAPPING_FILE, &CMappingTableDlg::OnFileChanged)
END_MESSAGE_MAP()

BOOL CMappingTableDlg::OnInitDialog()
{
   m_TableType = m_strTable.IsEmpty() ? 0 : 1;

   CDialog::OnInitDialog();

   UpdateControls();
   CheckTable();

   return TRUE;
}

void CMappingTableDlg::UpdateControls()
{
   BOOL bAgency = IsDlgButtonChecked(IDC_MAPPING_AGENCY) == BST_CHECKED;
   GetDlgItem(IDC_MAPPING_FILE)->EnableWindow(bAgency);
   GetDlgItem(IDC_MAPPING_BROWSE)->EnableWindow(bAgency);
}

void CMappingTableDlg::OnTableTypeChanged()
{
   UpdateControls();
   CheckTable();
}

void CMappingTableDlg::OnFileChanged()
{
   CheckTable();
}

void CMappingTableDlg::OnBrowse()
{
   CString strFile;
   GetDlgItemText(IDC_MAPPING_FILE, strFile);

   CFileDialog dlg(TRUE, _T("json"), strFile, OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,
      _T("IFC Mapping Tables (*.json)|*.json|All Files (*.*)|*.*||"), this);
   if (dlg.DoModal() == IDOK)
   {
      SetDlgItemText(IDC_MAPPING_FILE, dlg.GetPathName());
      CheckTable();
   }
}

bool CMappingTableDlg::CheckTable()
{
   CString strFile;
   GetDlgItemText(IDC_MAPPING_FILE, strFile);
   strFile.Trim();
   bool bAgency = IsDlgButtonChecked(IDC_MAPPING_AGENCY) == BST_CHECKED;

   std::ostringstream os;
   bool bOK = false;
   if (bAgency && strFile.IsEmpty())
   {
      os << "Choose a mapping table file.";
   }
   else
   {
      try
      {
         auto pTable = CIfcMappingTable::Load(bAgency ? std::filesystem::path(strFile.GetString()) : std::filesystem::path(),
            bAgency ? MappingTableSource::ConfigurationSetting : MappingTableSource::InstalledStandard);

         // the table and the tables it extends
         for (const auto& file : pTable->GetFiles())
         {
            os << (&file == &pTable->GetFiles().front() ? "" : "Extends: ") << file.name << " (version " << file.version << ")" << std::endl;
            os << "   " << PathToString(file.path) << std::endl;
         }
         os << std::endl << "The table can be used.";
         bOK = true;
      }
      catch (const std::exception& e)
      {
         os << e.what();
      }
   }

   CString strText(os.str().c_str());
   strText.Replace(_T("\r\n"), _T("\n"));
   strText.Replace(_T("\n"), _T("\r\n")); // multi-line edit controls need CR LF
   SetDlgItemText(IDC_MAPPING_STATUS, strText);
   return bOK;
}

void CMappingTableDlg::OnOK()
{
   if (!UpdateData(TRUE))
      return;

   m_strTable.Trim();
   if (!CheckTable())
   {
      // no silent fallback: a table that can't be used isn't saved (G5)
      AfxMessageBox(_T("This mapping table can't be used. See the message in the dialog, or choose the standard table."), MB_OK | MB_ICONEXCLAMATION);
      return;
   }

   if (m_TableType == 0)
      m_strTable.Empty();

   CDialog::OnOK();
}
