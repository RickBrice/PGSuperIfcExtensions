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
#include "MappingTableChooser.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <sstream>

CString GetMappingTableDisplayName(const std::filesystem::path& table)
{
   bool bStandard = table.empty() || table == CIfcMappingTable::STANDARD_TABLE;
   std::filesystem::path path = bStandard ? CIfcMappingTable::GetStandardTablePath() : table;
   try
   {
      std::ifstream file(path);
      auto json = nlohmann::ordered_json::parse(file);
      if (json.is_object() && json.contains("name") && json["name"].is_string())
         return CString(CA2T(json["name"].get<std::string>().c_str(), CP_UTF8));
   }
   catch (...)
   {
   }
   return bStandard ? CString(_T("Standard")) : CString(path.filename().c_str());
}

void CMappingTableChooser::Init(CWnd* pParent, const CString& default_label, const std::filesystem::path& default_table, MappingTableSource default_source, MappingTableSource other_source)
{
   m_pParent = pParent;
   m_DefaultTable = default_table;
   m_DefaultSource = default_source;
   m_OtherSource = other_source;
   m_pParent->SetDlgItemText(IDC_MAPPING_STANDARD, default_label);
   FillList(std::filesystem::path());
   SelectDefault();
}

void CMappingTableChooser::FillList(const std::filesystem::path& select)
{
   m_Tables = CIfcMappingTable::GetRecentTables();
   if (!select.empty() && std::find(m_Tables.begin(), m_Tables.end(), select) == m_Tables.end())
      m_Tables.insert(m_Tables.begin(), select);

   // the standard table can be chosen for one import or export, unless it's the default
   if (m_DefaultTable != CIfcMappingTable::STANDARD_TABLE)
      m_Tables.insert(m_Tables.begin(), std::filesystem::path(CIfcMappingTable::STANDARD_TABLE));

   CComboBox* pList = (CComboBox*)m_pParent->GetDlgItem(IDC_MAPPING_FILE);
   pList->ResetContent();
   int sel = -1;
   for (int i = 0; i < (int)m_Tables.size(); i++)
   {
      CString text;
      bool bStandard = (m_Tables[i] == CIfcMappingTable::STANDARD_TABLE);
      text.Format(_T("%s  (%s)"), GetMappingTableDisplayName(m_Tables[i]).GetString(), bStandard ? _T("installed with the IFC extension") : m_Tables[i].c_str());
      pList->AddString(text);
      if (m_Tables[i] == select)
         sel = i;
   }
   pList->SetCurSel(sel < 0 && !m_Tables.empty() ? 0 : sel);
}

void CMappingTableChooser::SelectDefault()
{
   m_pParent->CheckRadioButton(IDC_MAPPING_STANDARD, IDC_MAPPING_AGENCY, IDC_MAPPING_STANDARD);
   UpdateControls();
   Check();
}

void CMappingTableChooser::SelectOther(const std::filesystem::path& table)
{
   FillList(table);
   m_pParent->CheckRadioButton(IDC_MAPPING_STANDARD, IDC_MAPPING_AGENCY, IDC_MAPPING_AGENCY);
   UpdateControls();
   Check();
}

bool CMappingTableChooser::IsDefault() const
{
   return m_pParent->IsDlgButtonChecked(IDC_MAPPING_STANDARD) == BST_CHECKED;
}

std::filesystem::path CMappingTableChooser::GetOtherTable() const
{
   int sel = ((CComboBox*)m_pParent->GetDlgItem(IDC_MAPPING_FILE))->GetCurSel();
   return (0 <= sel && sel < (int)m_Tables.size()) ? m_Tables[sel] : std::filesystem::path();
}

std::filesystem::path CMappingTableChooser::GetTable() const
{
   return IsDefault() ? m_DefaultTable : GetOtherTable();
}

MappingTableSource CMappingTableChooser::GetSource() const
{
   return IsDefault() ? m_DefaultSource : m_OtherSource;
}

void CMappingTableChooser::UpdateControls()
{
   BOOL bOther = !IsDefault();
   m_pParent->GetDlgItem(IDC_MAPPING_FILE)->EnableWindow(bOther);
   m_pParent->GetDlgItem(IDC_MAPPING_BROWSE)->EnableWindow(bOther);
}

bool CMappingTableChooser::Check()
{
   std::ostringstream os;
   bool bOK = false;
   if (!IsDefault() && GetOtherTable().empty())
   {
      os << "Choose a mapping table file with the ... button.";
   }
   else
   {
      try
      {
         auto pTable = CIfcMappingTable::LoadActive(GetTable(), GetSource()); // STANDARD_TABLE: the installed standard table
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

   CString strText(CA2T(os.str().c_str(), CP_UTF8));
   strText.Replace(_T("\r\n"), _T("\n"));
   strText.Replace(_T("\n"), _T("\r\n")); // multi-line edit controls need CR LF
   m_pParent->SetDlgItemText(IDC_MAPPING_STATUS, strText);
   return bOK;
}

void CMappingTableChooser::OnTableTypeChanged()
{
   UpdateControls();
   Check();
}

void CMappingTableChooser::OnRecentChanged()
{
   Check();
}

void CMappingTableChooser::OnBrowse()
{
   auto current = GetOtherTable();
   CFileDialog dlg(TRUE, _T("ifcmap"), current.empty() ? nullptr : current.c_str(), OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,
      _T("IFC Mapping Tables (*.ifcmap)|*.ifcmap|JSON Files (*.json)|*.json|All Files (*.*)|*.*||"), m_pParent);
   if (dlg.DoModal() == IDOK)
   {
      FillList(std::filesystem::path(dlg.GetPathName().GetString()));
      Check();
   }
}
