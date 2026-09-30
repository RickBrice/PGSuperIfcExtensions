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

// MappingEditorModelPane.cpp : the mapping table editor's Model item (stage 4)
//

#include "stdafx.h"
#include "resource.h"
#include "MappingEditorModelPane.h"
#include "MappingEditorExportPanes.h"
#include "MappingEditorDoc.h"
#include "MappingEditorUtil.h"

using namespace mapping_editor;

BEGIN_MESSAGE_MAP(CModelPane, CMappingPane)
   ON_BN_CLICKED(IDC_MODEL_OPEN, &CModelPane::OnOpenModel)
   ON_CBN_SELCHANGE(IDC_MODEL_ROLE, &CModelPane::OnRoleChanged)
   ON_CBN_SELCHANGE(IDC_MODEL_ELEMENT, &CModelPane::OnElementChanged)
   ON_CBN_SELCHANGE(IDC_MODEL_TARGET, &CModelPane::OnTargetChanged)
   ON_BN_CLICKED(IDC_MODEL_ADD, &CModelPane::OnAdd)
   ON_NOTIFY(LVN_ITEMCHANGED, IDC_MODEL_PROPERTIES, &CModelPane::OnPropertiesChanged)
   ON_NOTIFY(NM_DBLCLK, IDC_MODEL_PROPERTIES, &CModelPane::OnPropertiesDblClk)
END_MESSAGE_MAP()

const CIfcMappingTable* CModelPane::GetTable(std::string* pError)
{
   std::string error;
   const auto* pTable = m_pDoc->GetCurrentTable(error);
   if (pError)
      *pError = error;
   return pTable;
}

BOOL CModelPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();

   m_Properties.SubclassDlgItem(IDC_MODEL_PROPERTIES, this);
   m_Properties.SetExtendedStyle(LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES);
   CRect rect;
   m_Properties.GetClientRect(&rect);
   int w = rect.Width();
   m_Properties.InsertColumn(0, _T("On"), LVCFMT_LEFT, w * 12 / 100);
   m_Properties.InsertColumn(1, _T("Property set"), LVCFMT_LEFT, w * 26 / 100);
   m_Properties.InsertColumn(2, _T("Property"), LVCFMT_LEFT, w * 24 / 100);
   m_Properties.InsertColumn(3, _T("Value"), LVCFMT_LEFT, w * 20 / 100);
   m_Properties.InsertColumn(4, _T("Type"), LVCFMT_LEFT, w * 18 / 100);

   CIfcModel* pModel = m_pDoc->GetModel();
   std::string error;
   const auto* pTable = GetTable(&error);
   if (!pModel)
   {
      SetDlgItemText(IDC_MODEL_INFO, _T("No model is open. Open an IFC 4.3 model to browse its elements and their properties, and to see what the table reads from them."));
   }
   else if (!pTable)
   {
      SetDlgItemText(IDC_MODEL_INFO, Utf8ToCString("Model: " + PathToString(pModel->GetPath())));
      SetDlgItemText(IDC_MODEL_READING, ToWindowsLines(_T("The table must be valid to browse the model with it.\n\n") + Utf8ToCString(error)));
   }
   else
   {
      SetDlgItemText(IDC_MODEL_INFO, Utf8ToCString("Model: " + PathToString(pModel->GetPath())));

      // the element roles the table defines, with how many elements play them
      CWaitCursor wait;
      CComboBox* pRole = (CComboBox*)GetDlgItem(IDC_MODEL_ROLE);
      for (int kind = (int)ElementKind::Project; kind <= (int)ElementKind::Barrier; kind++)
      {
         ElementKind role = (ElementKind)kind;
         if (!pTable->GetSelector(role))
            continue;
         auto count = pModel->GetElements(*pTable, role).size();
         CString label;
         label.Format(_T("%s (%d)"), Utf8ToCString(std::string(GetElementRoleName(role))).GetString(), (int)count);
         pRole->AddString(label);
         m_Roles.push_back(role);
      }

      // start with the girders, the role agency tables change most
      auto girder = std::find(m_Roles.begin(), m_Roles.end(), ElementKind::Girder);
      pRole->SetCurSel(girder == m_Roles.end() ? 0 : (int)(girder - m_Roles.begin()));
      OnRoleChanged();
   }

   for (int nID : { IDC_MODEL_ROLE, IDC_MODEL_ELEMENT, IDC_MODEL_PROPERTIES, IDC_MODEL_TARGET })
      GetDlgItem(nID)->EnableWindow(pModel && pTable);
   UpdateButtons();
   return TRUE;
}

void CModelPane::OnOpenModel()
{
   // the document replaces this pane when the model is open
   m_pDoc->OpenModel();
}

void CModelPane::OnRoleChanged()
{
   CIfcModel* pModel = m_pDoc->GetModel();
   const auto* pTable = GetTable();
   int sel = ((CComboBox*)GetDlgItem(IDC_MODEL_ROLE))->GetCurSel();
   if (!pModel || !pTable || sel < 0 || (int)m_Roles.size() <= sel)
      return;
   ElementKind role = m_Roles[sel];

   CWaitCursor wait;
   m_Elements = pModel->GetElements(*pTable, role);
   CComboBox* pElement = (CComboBox*)GetDlgItem(IDC_MODEL_ELEMENT);
   pElement->ResetContent();
   for (const auto& element : m_Elements)
      pElement->AddString(Utf8ToCString(element.label));
   pElement->SetCurSel(m_Elements.empty() ? -1 : 0);

   // the targets of the role's elements
   CString strTarget;
   GetDlgItemText(IDC_MODEL_TARGET, strTarget);
   CComboBox* pTarget = (CComboBox*)GetDlgItem(IDC_MODEL_TARGET);
   pTarget->ResetContent();
   m_Targets.clear();
   for (const auto& target : GetTargetDefs())
   {
      if (target.element == GetTargetElement(role))
      {
         m_Targets.emplace_back(target.name);
         pTarget->AddString(Utf8ToCString(std::string(target.name)));
      }
   }
   if (pTarget->SelectString(-1, strTarget) == CB_ERR)
      pTarget->SetCurSel(m_Targets.empty() ? -1 : 0);

   OnElementChanged();
}

void CModelPane::OnElementChanged()
{
   m_Properties.DeleteAllItems();
   m_ElementProperties.clear();
   CIfcModel* pModel = m_pDoc->GetModel();
   int id = GetElementId();
   if (pModel && id != 0)
   {
      m_ElementProperties = pModel->GetProperties(id);
      int i = 0;
      for (const auto& p : m_ElementProperties)
      {
         m_Properties.InsertItem(i, Utf8ToCString(OwnerName(p.owner)));
         m_Properties.SetItemText(i, 1, Utf8ToCString(p.pset));
         m_Properties.SetItemText(i, 2, Utf8ToCString(p.name));
         m_Properties.SetItemText(i, 3, Utf8ToCString(p.value));
         m_Properties.SetItemText(i, 4, Utf8ToCString(p.ifc_type));
         i++;
      }
   }
   UpdateReading();
   UpdateButtons();
}

void CModelPane::OnTargetChanged()
{
   UpdateReading();
   UpdateButtons();
}

int CModelPane::GetElementId()
{
   int sel = ((CComboBox*)GetDlgItem(IDC_MODEL_ELEMENT))->GetCurSel();
   return (0 <= sel && sel < (int)m_Elements.size()) ? m_Elements[sel].id : 0;
}

std::string CModelPane::GetTargetName()
{
   int sel = ((CComboBox*)GetDlgItem(IDC_MODEL_TARGET))->GetCurSel();
   return (0 <= sel && sel < (int)m_Targets.size()) ? m_Targets[sel] : std::string();
}

void CModelPane::UpdateReading()
{
   CIfcModel* pModel = m_pDoc->GetModel();
   const auto* pTable = GetTable();
   int id = GetElementId();
   std::string target = GetTargetName();
   if (!pModel || !pTable || id == 0 || target.empty())
   {
      if (pModel && pTable)
         SetDlgItemText(IDC_MODEL_READING, id == 0 ? _T("No element of this role in the model.") : _T(""));
      return;
   }
   SetDlgItemText(IDC_MODEL_READING, ToWindowsLines(Utf8ToCString(pModel->DescribeReading(*pTable, target, id))));
}

void CModelPane::UpdateButtons()
{
   int sel = m_Properties.GetSafeHwnd() ? m_Properties.GetNextItem(-1, LVNI_SELECTED) : -1;
   GetDlgItem(IDC_MODEL_ADD)->EnableWindow(0 <= sel && !GetTargetName().empty());
}

void CModelPane::OnPropertiesChanged(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   UpdateButtons();
}

void CModelPane::OnPropertiesDblClk(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   if (GetDlgItem(IDC_MODEL_ADD)->IsWindowEnabled())
      OnAdd();
}

void CModelPane::OnAdd()
{
   int sel = m_Properties.GetNextItem(-1, LVNI_SELECTED);
   std::string target_name = GetTargetName();
   const TargetDef* target = FindTargetDef(target_name);
   if (sel < 0 || (int)m_ElementProperties.size() <= sel || !target)
      return;

   // the property as a location, to complete in the location dialog (e.g. the unit of a plain number)
   const auto& p = m_ElementProperties[sel];
   ordered_json location = ordered_json::object();
   location["property"] = { { "pset", p.pset }, { "name", p.name } };
   if (p.owner != PropertyOwner::Occurrence)
      location["on"] = OwnerName(p.owner);

   CLocationDlg dlg(*target, location, this);
   if (dlg.DoModal() != IDOK)
      return;

   // appended: tried after this table's other locations for the target
   auto& table = m_pDoc->GetTable();
   if (!table.contains("targets"))
      table["targets"] = ordered_json::object();
   auto& item = table["targets"][target_name];
   if (item.is_object() && item.contains("locations") && item["locations"].is_array())
      item["locations"].push_back(dlg.m_Location);
   else if (item.is_array())
      item.push_back(dlg.m_Location);
   else
      item = ordered_json::array({ dlg.m_Location });

   m_pDoc->TableChanged(true);
   UpdateReading();
}
