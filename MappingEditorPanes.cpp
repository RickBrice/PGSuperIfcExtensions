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

// MappingEditorPanes.cpp : the mapping table editor's panes and dialogs
//

#include "stdafx.h"
#include "resource.h"
#include "MappingEditorPanes.h"
#include "MappingEditorDoc.h"
#include "IfcMappingTable.h"

#include <sstream>

using ordered_json = nlohmann::ordered_json;

namespace
{
   CString GetText(CWnd* pWnd, int nID)
   {
      CString text;
      pWnd->GetDlgItemText(nID, text);
      return text;
   }

   // Multi-line edit controls need CR LF
   CString ToWindowsLines(CString text)
   {
      text.Replace(_T("\r\n"), _T("\n"));
      text.Replace(_T("\n"), _T("\r\n"));
      return text;
   }

   std::string string_value(const ordered_json& j, const char* key)
   {
      return (j.is_object() && j.contains(key) && j[key].is_string()) ? j[key].get<std::string>() : std::string();
   }

   // Sets a text member, or removes it when the text is empty. Returns true if the object changed
   bool set_or_erase(ordered_json& j, const char* key, const std::string& value)
   {
      if (value.empty())
      {
         if (!j.contains(key))
            return false;
         j.erase(key);
         return true;
      }
      if (j.contains(key) && j[key] == value)
         return false;
      j[key] = value;
      return true;
   }

   const char* value_kind_name(ValueKind kind)
   {
      switch (kind)
      {
      case ValueKind::Stress: return "stress";
      case ValueKind::Length: return "length";
      case ValueKind::Angle: return "angle";
      case ValueKind::Force: return "force";
      case ValueKind::Area: return "area";
      case ValueKind::Mass: return "mass";
      case ValueKind::Ratio: return "ratio";
      case ValueKind::Count: return "count";
      case ValueKind::Boolean: return "true or false";
      case ValueKind::Text: return "text";
      case ValueKind::TextList: return "text (every location's value is collected)";
      }
      return "";
   }

   std::vector<std::string> unit_names(ValueKind kind)
   {
      std::vector<std::string> names;
      for (const auto& unit : GetTableUnits())
      {
         if (unit.kind == kind)
            names.emplace_back(unit.name);
      }
      return names;
   }

   bool has_map(ValueKind kind)
   {
      return kind == ValueKind::Boolean || kind == ValueKind::Text || kind == ValueKind::TextList;
   }

   ordered_json selector_json(const ElementSelector& selector)
   {
      ordered_json j = ordered_json::object();
      if (!selector.any_of.empty())
      {
         j["any_of"] = ordered_json::array();
         for (const auto& alternative : selector.any_of)
            j["any_of"].push_back(selector_json(alternative));
         return j;
      }
      j["entity"] = selector.entity;
      if (!selector.predefined_type.empty())
         j["predefined_type"] = selector.predefined_type;
      if (!selector.attributes.empty())
      {
         j["attributes"] = ordered_json::object();
         for (const auto& [name, value] : selector.attributes)
            j["attributes"][name] = value;
      }
      if (!selector.classification.empty())
         j["classification"] = selector.classification;
      return j;
   }

   // A one-line JSON value, as the table file shows it
   std::string inline_json(const ordered_json& j)
   {
      return j.dump(-1, ' ', false);
   }
}

CString DescribeLocation(const ordered_json& location)
{
   std::string text;
   if (location.contains("property") && location["property"].is_object())
   {
      text = string_value(location["property"], "pset") + "." + string_value(location["property"], "name");
      auto on = string_value(location, "on");
      if (!on.empty() && on != "occurrence")
         text += " (" + on + ")";
   }
   else if (location.contains("attribute"))
   {
      text = "attribute " + string_value(location, "attribute");
   }
   else if (location.contains("type_attribute"))
   {
      text = "type attribute " + string_value(location, "type_attribute");
   }
   else if (location.contains("classification"))
   {
      const auto& c = location["classification"];
      auto system = string_value(c, "system");
      auto identification = string_value(c, "identification");
      auto field = string_value(location, "field");
      text = "classification " + (system.empty() ? std::string("(any system)") : system) + (identification.empty() ? "" : " " + identification) +
         ", " + (field.empty() ? std::string("Name") : field);
   }
   else
   {
      text = inline_json(location);
   }
   return Utf8ToCString(text);
}

CString DescribeLocationReading(const ordered_json& location)
{
   std::vector<std::string> parts;
   if (location.contains("unit"))
      parts.push_back("unit " + string_value(location, "unit"));
   if (location.contains("parse"))
   {
      const auto& parse = location["parse"];
      parts.push_back(parse.is_string() ? "parse " + parse.get<std::string>() : "parse regex " + string_value(parse, "regex"));
   }
   if (location.contains("list_index") && location["list_index"].is_number())
      parts.push_back("list item " + std::to_string(location["list_index"].get<long long>()));
   if (location.contains("value_types"))
      parts.push_back("value types " + inline_json(location["value_types"]));
   if (location.contains("map") && location["map"].is_object())
      parts.push_back("map " + inline_json(location["map"]));
   if (location.contains("comment"))
      parts.push_back("comment: " + string_value(location, "comment"));

   std::string text;
   for (const auto& part : parts)
      text += (text.empty() ? "" : "; ") + part;
   return Utf8ToCString(text);
}

/////////////////////////////////////////////////////////////////////////////
// CreateMappingPane

namespace
{
   CString section_text(const ordered_json& table, const char* key, const CString& intro)
   {
      CString text(intro);
      text += _T("\n\n");
      if (!table.contains(key) || !table[key].is_array() || table[key].empty())
      {
         text += _T("This table has none.");
         return text;
      }

      std::ostringstream os;
      for (const auto& item : table[key])
      {
         if (item.is_object() && (item.contains("properties") || item.contains("quantities")))
         {
            const char* items_key = item.contains("properties") ? "properties" : "quantities";
            os << string_value(item, "name") << "   applies to " << inline_json(item.value("applies_to", ordered_json())) <<
               (item.contains("attach") ? ", attach " + string_value(item, "attach") : "") <<
               (item.contains("condition") ? ", condition " + string_value(item, "condition") : "") <<
               (item.value("remove", false) ? ", removed" : "") << std::endl;
            for (const auto& p : item[items_key])
            {
               os << "      " << string_value(p, "name") << "   " << string_value(p, "type");
               if (p.contains("target")) os << "   target " << string_value(p, "target");
               if (p.contains("value")) os << "   value " << inline_json(p["value"]);
               if (p.contains("import") && p["import"] == false) os << "   export only";
               os << std::endl;
            }
         }
         else
         {
            os << inline_json(item) << std::endl;
         }
      }
      text += Utf8ToCString(os.str());
      return text;
   }
}

std::unique_ptr<CMappingPane> CreateMappingPane(const MappingNode& node, CMappingEditorDoc* pDoc)
{
   const auto& table = pDoc->GetTable();
   const CString later(_T("Editing this section comes in a later stage of the editor; this is what the table has. A table that extends another table has only its own additions and changes here."));

   switch (node.kind)
   {
   case MappingNode::Kind::Table:
      return std::make_unique<CTablePane>(pDoc);

   case MappingNode::Kind::Role:
      return std::make_unique<CSelectorPane>(pDoc, node.key);

   case MappingNode::Kind::Target:
      if (const auto* target = FindTargetDef(node.key))
         return std::make_unique<CTargetPane>(pDoc, *target);
      return nullptr;

   case MappingNode::Kind::Roles:
      return std::make_unique<CInfoPane>(pDoc, _T("Element roles say which IFC elements play each role (girder, deck, bearing, ...). Select a role to see or change its selector.\n\nRoles in bold have a selector in this table; the others use the base table's selector."));

   case MappingNode::Kind::Targets:
      return std::make_unique<CInfoPane>(pDoc, _T("Targets are the PGSuper values an IFC import reads. Each has an ordered list of locations; the first location that has a value wins. Select a target to see or change its locations.\n\nTargets in bold have locations in this table. Their locations are tried before the base table's, unless they replace them.\n\nA target also reads the properties declared with it in the property sets (unless \"import\" is false)."));

   case MappingNode::Kind::TargetGroup:
      return std::make_unique<CInfoPane>(pDoc, Utf8ToCString("The targets of the " + node.key + " elements. Select a target to see or change its locations."));

   case MappingNode::Kind::PropertySets:
      return std::make_unique<CInfoPane>(pDoc, section_text(table, "property_sets", _T("Property sets exported for each element role, and the targets bound to their properties. ") + later));

   case MappingNode::Kind::QuantitySets:
      return std::make_unique<CInfoPane>(pDoc, section_text(table, "quantity_sets", _T("Quantity sets exported for each element role. ") + later));

   case MappingNode::Kind::ClassificationSystems:
      return std::make_unique<CInfoPane>(pDoc, section_text(table, "classification_systems", _T("Classification systems written to the project. ") + later));

   case MappingNode::Kind::Classifications:
      return std::make_unique<CInfoPane>(pDoc, section_text(table, "classifications", _T("Classification references of each element role. ") + later));
   }
   return nullptr;
}

/////////////////////////////////////////////////////////////////////////////
// CTablePane

BEGIN_MESSAGE_MAP(CTablePane, CMappingPane)
   ON_EN_KILLFOCUS(IDC_TABLE_NAME, &CTablePane::OnNameChanged)
   ON_EN_KILLFOCUS(IDC_TABLE_COMMENT, &CTablePane::OnCommentChanged)
   ON_CBN_SELCHANGE(IDC_TABLE_EXTENDS, &CTablePane::OnExtendsSelected)
   ON_CBN_KILLFOCUS(IDC_TABLE_EXTENDS, &CTablePane::OnExtendsChanged)
   ON_BN_CLICKED(IDC_TABLE_EXTENDS_BROWSE, &CTablePane::OnBrowseExtends)
END_MESSAGE_MAP()

static const TCHAR* const NO_BASE_TABLE = _T("(none)");

BOOL CTablePane::OnInitDialog()
{
   CMappingPane::OnInitDialog();

   const auto& table = m_pDoc->GetTable();
   SetDlgItemText(IDC_TABLE_NAME, Utf8ToCString(string_value(table, "name")));
   SetDlgItemText(IDC_TABLE_COMMENT, ToWindowsLines(Utf8ToCString(string_value(table, "comment"))));

   CComboBox* pExtends = (CComboBox*)GetDlgItem(IDC_TABLE_EXTENDS);
   pExtends->AddString(_T("standard"));
   pExtends->AddString(NO_BASE_TABLE);
   auto extends = string_value(table, "extends");
   pExtends->SetWindowText(extends.empty() ? CString(NO_BASE_TABLE) : Utf8ToCString(extends));

   UpdateInfo();
   return TRUE;
}

void CTablePane::UpdateInfo()
{
   std::ostringstream os;
   os << "File: " << (m_pDoc->GetPathName().IsEmpty() ? std::string("not saved yet") : PathToString(m_pDoc->GetTablePath())) << std::endl;
   if (const auto* pBase = m_pDoc->GetBaseTable())
   {
      os << "Base table: ";
      bool bFirst = true;
      for (const auto& file : pBase->GetFiles())
      {
         os << (bFirst ? "" : ", which extends ") << "\"" << file.name << "\"";
         bFirst = false;
      }
      os << std::endl;
   }
   else if (!m_pDoc->GetBaseTableError().empty())
   {
      os << "The base table can't be loaded:" << std::endl << m_pDoc->GetBaseTableError();
   }
   else
   {
      os << "The table stands alone." << std::endl;
   }
   SetDlgItemText(IDC_TABLE_INFO, ToWindowsLines(Utf8ToCString(os.str()))); // a read-only edit, so long paths wrap and can be copied
}

void CTablePane::OnNameChanged()
{
   auto& table = m_pDoc->GetTable();
   auto name = CStringToUtf8(GetText(this, IDC_TABLE_NAME).Trim());
   if (string_value(table, "name") != name)
   {
      table["name"] = name;
      m_pDoc->TableChanged(false);
   }
}

void CTablePane::OnCommentChanged()
{
   CString text = GetText(this, IDC_TABLE_COMMENT);
   text.Replace(_T("\r\n"), _T("\n"));
   text.Trim();
   if (set_or_erase(m_pDoc->GetTable(), "comment", CStringToUtf8(text)))
      m_pDoc->TableChanged(false);
}

void CTablePane::OnExtendsSelected()
{
   CComboBox* pExtends = (CComboBox*)GetDlgItem(IDC_TABLE_EXTENDS);
   int sel = pExtends->GetCurSel();
   if (sel != CB_ERR)
   {
      CString text;
      pExtends->GetLBText(sel, text);
      SetExtends(text);
   }
}

void CTablePane::OnExtendsChanged()
{
   SetExtends(GetText(this, IDC_TABLE_EXTENDS));
}

void CTablePane::OnBrowseExtends()
{
   CFileDialog dlg(TRUE, _T("json"), nullptr, OFN_HIDEREADONLY | OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,
      _T("IFC Mapping Tables (*.json)|*.json|All Files (*.*)|*.*||"), this);
   if (dlg.DoModal() != IDOK)
      return;

   // relative to this table's folder, so the tables can be moved together
   std::error_code ec;
   std::filesystem::path file(dlg.GetPathName().GetString());
   auto relative = std::filesystem::relative(file, m_pDoc->GetTablePath().parent_path(), ec);
   CString strExtends = (ec || relative.empty()) ? CString(file.c_str()) : CString(relative.c_str());
   strExtends.Replace(_T('\\'), _T('/'));
   SetDlgItemText(IDC_TABLE_EXTENDS, strExtends);
   SetExtends(strExtends);
}

void CTablePane::SetExtends(CString strExtends)
{
   strExtends.Trim();
   if (strExtends == NO_BASE_TABLE)
      strExtends.Empty();

   if (set_or_erase(m_pDoc->GetTable(), "extends", CStringToUtf8(strExtends)))
   {
      CWaitCursor wait;
      m_pDoc->LoadBaseTable();
      m_pDoc->TableChanged(true);
      UpdateInfo();
   }
}

/////////////////////////////////////////////////////////////////////////////
// CSelectorPane

BEGIN_MESSAGE_MAP(CSelectorPane, CMappingPane)
   ON_BN_CLICKED(IDC_SELECTOR_INHERITED, &CSelectorPane::OnInherited)
   ON_BN_CLICKED(IDC_SELECTOR_DEFINED, &CSelectorPane::OnDefined)
   ON_EN_KILLFOCUS(IDC_SELECTOR_ENTITY, &CSelectorPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SELECTOR_PREDEFINED, &CSelectorPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SELECTOR_CLASSIFICATION, &CSelectorPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SELECTOR_ATTRIBUTES, &CSelectorPane::OnFieldChanged)
END_MESSAGE_MAP()

BOOL CSelectorPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();

   ElementKind role;
   GetElementKind(m_Role, role);
   std::ostringstream os;
   os << "Element role: " << m_Role;
   ElementKind target_element = GetTargetElement(role);
   if (target_element != role)
      os << " (its elements use the " << GetElementRoleName(target_element) << " targets)";
   SetDlgItemText(IDC_SELECTOR_HEADER, Utf8ToCString(os.str()));

   Fill();
   return TRUE;
}

void CSelectorPane::Fill()
{
   const auto& table = m_pDoc->GetTable();
   bool bDefined = table.contains("elements") && table["elements"].is_object() && table["elements"].contains(m_Role);

   ElementKind role;
   GetElementKind(m_Role, role);
   const ElementSelector* pBase = m_pDoc->GetBaseTable() ? m_pDoc->GetBaseTable()->GetSelector(role) : nullptr;
   CString strInherited = pBase ? _T("Use the base table's selector: ") + Utf8ToCString(pBase->Describe()) :
      (m_pDoc->GetBaseTable() ? CString(_T("Use the base table's selector (it has none: the role isn't used)")) : CString(_T("No selector (the role isn't used)")));
   SetDlgItemText(IDC_SELECTOR_INHERITED, strInherited);
   CheckRadioButton(IDC_SELECTOR_INHERITED, IDC_SELECTOR_DEFINED, bDefined ? IDC_SELECTOR_DEFINED : IDC_SELECTOR_INHERITED);

   ordered_json selector = bDefined ? table["elements"][m_Role] : (pBase ? selector_json(*pBase) : ordered_json::object());
   bool bAlternatives = selector.contains("any_of");

   SetDlgItemText(IDC_SELECTOR_ENTITY, Utf8ToCString(string_value(selector, "entity")));
   SetDlgItemText(IDC_SELECTOR_PREDEFINED, Utf8ToCString(string_value(selector, "predefined_type")));
   SetDlgItemText(IDC_SELECTOR_CLASSIFICATION, Utf8ToCString(string_value(selector, "classification")));
   CString strAttributes;
   if (selector.contains("attributes") && selector["attributes"].is_object())
   {
      for (const auto& [name, value] : selector["attributes"].items())
         strAttributes += Utf8ToCString(name + "=" + (value.is_string() ? value.get<std::string>() : value.dump())) + _T("\r\n");
   }
   SetDlgItemText(IDC_SELECTOR_ATTRIBUTES, strAttributes);

   BOOL bEdit = bDefined && !bAlternatives;
   for (int nID : { IDC_SELECTOR_ENTITY, IDC_SELECTOR_PREDEFINED, IDC_SELECTOR_CLASSIFICATION, IDC_SELECTOR_ATTRIBUTES })
      ((CEdit*)GetDlgItem(nID))->SetReadOnly(!bEdit);

   SetDlgItemText(IDC_SELECTOR_NOTE, bAlternatives ?
      Utf8ToCString("This selector has alternatives (any_of), which this version of the editor shows but doesn't edit:\n" + inline_json(selector)) :
      CString(_T("The loader checks entity and attribute names against the IFC schema (Table > Validate).")));
}

void CSelectorPane::OnInherited()
{
   auto& table = m_pDoc->GetTable();
   if (!table.contains("elements") || !table["elements"].contains(m_Role))
      return;

   table["elements"].erase(m_Role);
   if (table["elements"].empty())
      table.erase("elements");
   m_pDoc->TableChanged(true);
   Fill();
}

void CSelectorPane::OnDefined()
{
   auto& table = m_pDoc->GetTable();
   if (table.contains("elements") && table["elements"].contains(m_Role))
      return;

   // start from the base table's selector
   ElementKind role;
   GetElementKind(m_Role, role);
   const ElementSelector* pBase = m_pDoc->GetBaseTable() ? m_pDoc->GetBaseTable()->GetSelector(role) : nullptr;
   if (!table.contains("elements"))
      table["elements"] = ordered_json::object();
   table["elements"][m_Role] = pBase ? selector_json(*pBase) : ordered_json{ { "entity", "" } };
   m_pDoc->TableChanged(true);
   Fill();
}

void CSelectorPane::OnFieldChanged()
{
   auto& table = m_pDoc->GetTable();
   if (!table.contains("elements") || !table["elements"].contains(m_Role) || table["elements"][m_Role].contains("any_of"))
      return;

   auto& selector = table["elements"][m_Role];
   ordered_json before = selector;

   selector["entity"] = CStringToUtf8(GetText(this, IDC_SELECTOR_ENTITY).Trim());
   CString strPredefined = GetText(this, IDC_SELECTOR_PREDEFINED).Trim();
   strPredefined.MakeUpper();
   set_or_erase(selector, "predefined_type", CStringToUtf8(strPredefined));
   set_or_erase(selector, "classification", CStringToUtf8(GetText(this, IDC_SELECTOR_CLASSIFICATION).Trim()));

   // Name=Value lines
   ordered_json attributes = ordered_json::object();
   CString strAttributes = GetText(this, IDC_SELECTOR_ATTRIBUTES);
   int pos = 0;
   CString line = strAttributes.Tokenize(_T("\r\n"), pos);
   while (!line.IsEmpty() || pos != -1)
   {
      line.Trim();
      int eq = line.Find(_T('='));
      if (0 < eq)
         attributes[CStringToUtf8(line.Left(eq).Trim())] = CStringToUtf8(line.Mid(eq + 1).Trim());
      if (pos == -1)
         break;
      line = strAttributes.Tokenize(_T("\r\n"), pos);
   }
   if (attributes.empty())
      selector.erase("attributes");
   else
      selector["attributes"] = attributes;

   if (selector != before)
      m_pDoc->TableChanged(false);
}

/////////////////////////////////////////////////////////////////////////////
// CTargetPane

BEGIN_MESSAGE_MAP(CTargetPane, CMappingPane)
   ON_BN_CLICKED(IDC_LOCATION_ADD, &CTargetPane::OnAdd)
   ON_BN_CLICKED(IDC_LOCATION_EDIT, &CTargetPane::OnEdit)
   ON_BN_CLICKED(IDC_LOCATION_REMOVE, &CTargetPane::OnRemove)
   ON_BN_CLICKED(IDC_LOCATION_UP, &CTargetPane::OnUp)
   ON_BN_CLICKED(IDC_LOCATION_DOWN, &CTargetPane::OnDown)
   ON_BN_CLICKED(IDC_REPLACE_BASE, &CTargetPane::OnReplace)
   ON_NOTIFY(NM_DBLCLK, IDC_LOCATIONS, &CTargetPane::OnLocationsDblClk)
   ON_NOTIFY(LVN_ITEMCHANGED, IDC_LOCATIONS, &CTargetPane::OnLocationsChanged)
END_MESSAGE_MAP()

BOOL CTargetPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();

   m_Locations.SubclassDlgItem(IDC_LOCATIONS, this);
   m_BaseLocations.SubclassDlgItem(IDC_BASE_LOCATIONS, this);
   m_Locations.SetExtendedStyle(LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES);
   m_BaseLocations.SetExtendedStyle(LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES);

   CRect rect;
   m_Locations.GetClientRect(&rect);
   m_Locations.InsertColumn(0, _T("Location"), LVCFMT_LEFT, rect.Width() / 2);
   m_Locations.InsertColumn(1, _T("Reading"), LVCFMT_LEFT, rect.Width() / 2);
   m_BaseLocations.GetClientRect(&rect);
   m_BaseLocations.InsertColumn(0, _T("Location"), LVCFMT_LEFT, rect.Width() * 2 / 3);
   m_BaseLocations.InsertColumn(1, _T("Table"), LVCFMT_LEFT, rect.Width() / 3);

   std::ostringstream os;
   os << m_Target.name << ": " << m_Target.description << std::endl;
   os << "Value: " << value_kind_name(m_Target.kind);
   auto units = unit_names(m_Target.kind);
   if (!units.empty())
   {
      os << ". A plain number needs a unit:";
      for (const auto& unit : units)
         os << " " << unit;
   }
   os << std::endl << "Element role: " << GetElementRoleName(m_Target.element);
   SetDlgItemText(IDC_TARGET_HEADER, Utf8ToCString(os.str()));

   // the base table's locations, in the order the importer tries them
   if (const auto* pBase = m_pDoc->GetBaseTable())
   {
      int i = 0;
      for (const auto& location : pBase->GetLocations(m_Target))
      {
         m_BaseLocations.InsertItem(i, Utf8ToCString(location.Describe()));
         m_BaseLocations.SetItemText(i, 1, Utf8ToCString(location.table));
         i++;
      }
   }
   else
   {
      m_BaseLocations.InsertItem(0, m_pDoc->GetBaseTableError().empty() ? _T("(no base table)") : _T("(the base table can't be loaded; see the Table item)"));
   }

   CheckDlgButton(IDC_REPLACE_BASE, IsReplace() ? BST_CHECKED : BST_UNCHECKED);
   GetDlgItem(IDC_REPLACE_BASE)->EnableWindow(m_pDoc->GetBaseTable() != nullptr || IsReplace());

   FillLocations(-1);
   return TRUE;
}

ordered_json CTargetPane::GetLocations() const
{
   const auto& table = m_pDoc->GetTable();
   std::string name(m_Target.name);
   if (!table.contains("targets") || !table["targets"].contains(name))
      return ordered_json::array();

   const auto& item = table["targets"][name];
   if (item.is_array())
      return item;
   if (item.is_object() && item.contains("locations") && item["locations"].is_array())
      return item["locations"];
   return ordered_json::array();
}

bool CTargetPane::IsReplace() const
{
   const auto& table = m_pDoc->GetTable();
   std::string name(m_Target.name);
   return table.contains("targets") && table["targets"].contains(name) && table["targets"][name].is_object() &&
      string_value(table["targets"][name], "mode") == "replace";
}

void CTargetPane::SetLocations(const ordered_json& locations, bool bReplace)
{
   auto& table = m_pDoc->GetTable();
   std::string name(m_Target.name);

   if (locations.empty() && !bReplace)
   {
      if (table.contains("targets"))
      {
         table["targets"].erase(name);
         if (table["targets"].empty())
            table.erase("targets");
      }
   }
   else
   {
      if (!table.contains("targets"))
         table["targets"] = ordered_json::object();

      if (bReplace)
      {
         // keep the object's other members (e.g. a comment)
         ordered_json item = table["targets"].contains(name) && table["targets"][name].is_object() ? table["targets"][name] : ordered_json::object();
         item["mode"] = "replace";
         item["locations"] = locations;
         table["targets"][name] = item;
      }
      else
      {
         table["targets"][name] = locations;
      }
   }
   m_pDoc->TableChanged(true);
}

void CTargetPane::FillLocations(int select)
{
   m_Locations.DeleteAllItems();
   int i = 0;
   for (const auto& location : GetLocations())
   {
      m_Locations.InsertItem(i, DescribeLocation(location));
      m_Locations.SetItemText(i, 1, DescribeLocationReading(location));
      i++;
   }
   if (0 <= select && select < m_Locations.GetItemCount())
   {
      m_Locations.SetItemState(select, LVIS_SELECTED | LVIS_FOCUSED, LVIS_SELECTED | LVIS_FOCUSED);
      m_Locations.EnsureVisible(select, FALSE);
   }
   UpdateButtons();
}

void CTargetPane::UpdateButtons()
{
   int sel = m_Locations.GetNextItem(-1, LVNI_SELECTED);
   int count = m_Locations.GetItemCount();
   GetDlgItem(IDC_LOCATION_EDIT)->EnableWindow(0 <= sel);
   GetDlgItem(IDC_LOCATION_REMOVE)->EnableWindow(0 <= sel);
   GetDlgItem(IDC_LOCATION_UP)->EnableWindow(0 < sel);
   GetDlgItem(IDC_LOCATION_DOWN)->EnableWindow(0 <= sel && sel < count - 1);
}

void CTargetPane::OnLocationsChanged(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   UpdateButtons();
}

void CTargetPane::OnLocationsDblClk(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   OnEdit();
}

void CTargetPane::OnAdd()
{
   CLocationDlg dlg(m_Target, ordered_json::object(), this);
   if (dlg.DoModal() == IDOK)
   {
      auto locations = GetLocations();
      locations.push_back(dlg.m_Location);
      SetLocations(locations, IsReplace());
      FillLocations((int)locations.size() - 1);
   }
}

void CTargetPane::OnEdit()
{
   int sel = m_Locations.GetNextItem(-1, LVNI_SELECTED);
   auto locations = GetLocations();
   if (sel < 0 || (int)locations.size() <= sel)
      return;

   CLocationDlg dlg(m_Target, locations[sel], this);
   if (dlg.DoModal() == IDOK && dlg.m_Location != locations[sel])
   {
      locations[sel] = dlg.m_Location;
      SetLocations(locations, IsReplace());
      FillLocations(sel);
   }
}

void CTargetPane::OnRemove()
{
   int sel = m_Locations.GetNextItem(-1, LVNI_SELECTED);
   auto locations = GetLocations();
   if (sel < 0 || (int)locations.size() <= sel)
      return;

   locations.erase(locations.begin() + sel);
   SetLocations(locations, IsReplace());
   FillLocations(std::min(sel, (int)locations.size() - 1));
}

void CTargetPane::OnUp()
{
   int sel = m_Locations.GetNextItem(-1, LVNI_SELECTED);
   auto locations = GetLocations();
   if (sel <= 0 || (int)locations.size() <= sel)
      return;

   std::swap(locations[sel - 1], locations[sel]);
   SetLocations(locations, IsReplace());
   FillLocations(sel - 1);
}

void CTargetPane::OnDown()
{
   int sel = m_Locations.GetNextItem(-1, LVNI_SELECTED);
   auto locations = GetLocations();
   if (sel < 0 || (int)locations.size() - 1 <= sel)
      return;

   std::swap(locations[sel], locations[sel + 1]);
   SetLocations(locations, IsReplace());
   FillLocations(sel + 1);
}

void CTargetPane::OnReplace()
{
   bool bReplace = IsDlgButtonChecked(IDC_REPLACE_BASE) == BST_CHECKED;
   if (bReplace != IsReplace())
      SetLocations(GetLocations(), bReplace);
}

/////////////////////////////////////////////////////////////////////////////
// CInfoPane

BEGIN_MESSAGE_MAP(CInfoPane, CMappingPane)
   ON_WM_SIZE()
END_MESSAGE_MAP()

BOOL CInfoPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();
   SetDlgItemText(IDC_INFO_TEXT, ToWindowsLines(m_Text));
   return TRUE;
}

void CInfoPane::OnSize(UINT nType, int cx, int cy)
{
   CMappingPane::OnSize(nType, cx, cy);
   if (CWnd* pText = GetDlgItem(IDC_INFO_TEXT))
   {
      CRect margin(0, 0, 7, 7);
      MapDialogRect(&margin);
      pText->SetWindowPos(nullptr, margin.right, margin.bottom, std::max(0, (int)(cx - 2 * margin.right)), std::max(0, (int)(cy - 2 * margin.bottom)), SWP_NOZORDER | SWP_NOACTIVATE);
   }
}

/////////////////////////////////////////////////////////////////////////////
// CLocationDlg

BEGIN_MESSAGE_MAP(CLocationDlg, CDialog)
   ON_BN_CLICKED(IDC_LOC_PROPERTY, &CLocationDlg::OnKindChanged)
   ON_BN_CLICKED(IDC_LOC_ATTRIBUTE, &CLocationDlg::OnKindChanged)
   ON_BN_CLICKED(IDC_LOC_TYPE_ATTRIBUTE, &CLocationDlg::OnKindChanged)
   ON_BN_CLICKED(IDC_LOC_CLASSIFICATION, &CLocationDlg::OnKindChanged)
   ON_CBN_SELCHANGE(IDC_LOC_PARSE, &CLocationDlg::OnParseChanged)
END_MESSAGE_MAP()

namespace
{
   const TCHAR* const NO_UNIT = _T("(none)");
   const TCHAR* const ON_VALUES[] = { _T("occurrence"), _T("type"), _T("material") };
   const TCHAR* const ON_LABELS[] = { _T("the element, then its type"), _T("the element's type only"), _T("the element's material") };
   const TCHAR* const PARSE_VALUES[] = { _T(""), _T("number"), _T("feet_inches"), _T("regex") };
   const TCHAR* const PARSE_LABELS[] = { _T("(no parsing)"), _T("number"), _T("feet and inches"), _T("regex") };
}

CLocationDlg::CLocationDlg(const TargetDef& target, const ordered_json& location, CWnd* pParent)
   : CDialog(IDD_MAPPING_LOCATION, pParent), m_Location(location), m_Target(target)
{
}

int CLocationDlg::GetKind() const
{
   return GetCheckedRadioButton(IDC_LOC_PROPERTY, IDC_LOC_CLASSIFICATION);
}

BOOL CLocationDlg::OnInitDialog()
{
   CDialog::OnInitDialog();

   CString strCaption;
   strCaption.Format(_T("Location of %s"), Utf8ToCString(std::string(m_Target.name)).GetString());
   SetWindowText(strCaption);

   const auto& loc = m_Location;
   int kind = IDC_LOC_PROPERTY;
   if (loc.contains("attribute")) kind = IDC_LOC_ATTRIBUTE;
   else if (loc.contains("type_attribute")) kind = IDC_LOC_TYPE_ATTRIBUTE;
   else if (loc.contains("classification")) kind = IDC_LOC_CLASSIFICATION;
   CheckRadioButton(IDC_LOC_PROPERTY, IDC_LOC_CLASSIFICATION, kind);

   if (loc.contains("property"))
   {
      SetDlgItemText(IDC_LOC_PSET, Utf8ToCString(string_value(loc["property"], "pset")));
      SetDlgItemText(IDC_LOC_NAME, Utf8ToCString(string_value(loc["property"], "name")));
   }
   else if (loc.contains("attribute") || loc.contains("type_attribute"))
   {
      SetDlgItemText(IDC_LOC_NAME, Utf8ToCString(string_value(loc, loc.contains("attribute") ? "attribute" : "type_attribute")));
   }

   CComboBox* pOn = (CComboBox*)GetDlgItem(IDC_LOC_ON);
   for (auto label : ON_LABELS)
      pOn->AddString(label);
   auto on = Utf8ToCString(string_value(loc, "on"));
   pOn->SetCurSel(0);
   for (int i = 0; i < _countof(ON_VALUES); i++)
   {
      if (on == ON_VALUES[i])
         pOn->SetCurSel(i);
   }

   if (loc.contains("classification"))
   {
      SetDlgItemText(IDC_LOC_SYSTEM, Utf8ToCString(string_value(loc["classification"], "system")));
      SetDlgItemText(IDC_LOC_IDENTIFICATION, Utf8ToCString(string_value(loc["classification"], "identification")));
   }
   CComboBox* pField = (CComboBox*)GetDlgItem(IDC_LOC_FIELD);
   pField->AddString(_T("Name"));
   pField->AddString(_T("Identification"));
   pField->SetCurSel(string_value(loc, "field") == "Identification" ? 1 : 0);

   CComboBox* pUnit = (CComboBox*)GetDlgItem(IDC_LOC_UNIT);
   pUnit->AddString(NO_UNIT);
   for (const auto& unit : unit_names(m_Target.kind))
      pUnit->AddString(Utf8ToCString(unit));
   auto unit = Utf8ToCString(string_value(loc, "unit"));
   if (pUnit->SelectString(-1, unit.IsEmpty() ? CString(NO_UNIT) : unit) == CB_ERR)
   {
      // a unit that doesn't fit the target: keep it, so validation reports it
      pUnit->SetCurSel(pUnit->AddString(unit));
   }

   CComboBox* pParse = (CComboBox*)GetDlgItem(IDC_LOC_PARSE);
   for (auto label : PARSE_LABELS)
      pParse->AddString(label);
   int parse = 0;
   if (loc.contains("parse"))
   {
      const auto& jparse = loc["parse"];
      if (jparse.is_string())
      {
         auto text = Utf8ToCString(jparse.get<std::string>());
         for (int i = 1; i < 3; i++)
         {
            if (text == PARSE_VALUES[i])
               parse = i;
         }
      }
      else if (jparse.is_object())
      {
         parse = 3;
         SetDlgItemText(IDC_LOC_REGEX, Utf8ToCString(string_value(jparse, "regex")));
         if (jparse.contains("group") && jparse["group"].is_number())
            SetDlgItemInt(IDC_LOC_GROUP, (UINT)jparse["group"].get<long long>());
      }
   }
   pParse->SetCurSel(parse);
   if (parse != 3)
      SetDlgItemInt(IDC_LOC_GROUP, 1);

   if (loc.contains("list_index") && loc["list_index"].is_number())
      SetDlgItemInt(IDC_LOC_LIST_INDEX, (UINT)loc["list_index"].get<long long>());

   if (loc.contains("value_types") && loc["value_types"].is_array())
   {
      CString types;
      for (const auto& type : loc["value_types"])
         types += (types.IsEmpty() ? _T("") : _T(", ")) + Utf8ToCString(type.is_string() ? type.get<std::string>() : type.dump());
      SetDlgItemText(IDC_LOC_VALUE_TYPES, types);
   }

   if (loc.contains("map") && loc["map"].is_object())
   {
      CString map;
      for (const auto& [key, value] : loc["map"].items())
         map += Utf8ToCString(key + " = " + (value.is_string() ? value.get<std::string>() : value.dump())) + _T("\r\n");
      SetDlgItemText(IDC_LOC_MAP, map);
   }

   SetDlgItemText(IDC_LOC_COMMENT, Utf8ToCString(string_value(loc, "comment")));

   UpdateControls();
   return TRUE;
}

void CLocationDlg::OnKindChanged()
{
   UpdateControls();
}

void CLocationDlg::OnParseChanged()
{
   UpdateControls();
}

void CLocationDlg::UpdateControls()
{
   int kind = GetKind();
   GetDlgItem(IDC_LOC_PSET)->EnableWindow(kind == IDC_LOC_PROPERTY);
   GetDlgItem(IDC_LOC_ON)->EnableWindow(kind == IDC_LOC_PROPERTY);
   GetDlgItem(IDC_LOC_NAME)->EnableWindow(kind != IDC_LOC_CLASSIFICATION);
   GetDlgItem(IDC_LOC_SYSTEM)->EnableWindow(kind == IDC_LOC_CLASSIFICATION);
   GetDlgItem(IDC_LOC_IDENTIFICATION)->EnableWindow(kind == IDC_LOC_CLASSIFICATION);
   GetDlgItem(IDC_LOC_FIELD)->EnableWindow(kind == IDC_LOC_CLASSIFICATION);

   CComboBox* pUnit = (CComboBox*)GetDlgItem(IDC_LOC_UNIT);
   pUnit->EnableWindow(1 < pUnit->GetCount());

   bool bRegex = ((CComboBox*)GetDlgItem(IDC_LOC_PARSE))->GetCurSel() == 3;
   GetDlgItem(IDC_LOC_REGEX)->EnableWindow(bRegex);
   GetDlgItem(IDC_LOC_GROUP)->EnableWindow(bRegex);

   GetDlgItem(IDC_LOC_MAP)->EnableWindow(has_map(m_Target.kind));
}

void CLocationDlg::OnOK()
{
   // start from the location as it was, so its keys keep their order
   ordered_json loc = m_Location.is_object() ? m_Location : ordered_json::object();
   for (const char* key : { "property", "attribute", "type_attribute", "classification", "field", "on" })
      loc.erase(key);

   int kind = GetKind();
   CString strName = GetText(this, IDC_LOC_NAME).Trim();
   if (kind == IDC_LOC_PROPERTY)
   {
      CString strPset = GetText(this, IDC_LOC_PSET).Trim();
      if (strPset.IsEmpty() || strName.IsEmpty())
      {
         AfxMessageBox(_T("Enter the property set and the property."), MB_OK | MB_ICONEXCLAMATION);
         return;
      }
      loc["property"] = { { "pset", CStringToUtf8(strPset) }, { "name", CStringToUtf8(strName) } };
      int on = ((CComboBox*)GetDlgItem(IDC_LOC_ON))->GetCurSel();
      if (0 < on)
         loc["on"] = CStringToUtf8(ON_VALUES[on]);
   }
   else if (kind == IDC_LOC_ATTRIBUTE || kind == IDC_LOC_TYPE_ATTRIBUTE)
   {
      if (strName.IsEmpty())
      {
         AfxMessageBox(_T("Enter the attribute."), MB_OK | MB_ICONEXCLAMATION);
         return;
      }
      loc[kind == IDC_LOC_ATTRIBUTE ? "attribute" : "type_attribute"] = CStringToUtf8(strName);
   }
   else
   {
      ordered_json classification = ordered_json::object();
      set_or_erase(classification, "system", CStringToUtf8(GetText(this, IDC_LOC_SYSTEM).Trim()));
      set_or_erase(classification, "identification", CStringToUtf8(GetText(this, IDC_LOC_IDENTIFICATION).Trim()));
      loc["classification"] = classification;
      if (((CComboBox*)GetDlgItem(IDC_LOC_FIELD))->GetCurSel() == 1)
         loc["field"] = "Identification";
   }

   CString strUnit;
   CComboBox* pUnit = (CComboBox*)GetDlgItem(IDC_LOC_UNIT);
   if (pUnit->GetCurSel() != CB_ERR)
      pUnit->GetLBText(pUnit->GetCurSel(), strUnit);
   set_or_erase(loc, "unit", strUnit == NO_UNIT ? std::string() : CStringToUtf8(strUnit));

   int parse = ((CComboBox*)GetDlgItem(IDC_LOC_PARSE))->GetCurSel();
   if (parse == 1 || parse == 2)
   {
      loc["parse"] = CStringToUtf8(PARSE_VALUES[parse]);
   }
   else if (parse == 3)
   {
      CString strRegex = GetText(this, IDC_LOC_REGEX);
      if (strRegex.IsEmpty())
      {
         AfxMessageBox(_T("Enter the regular expression."), MB_OK | MB_ICONEXCLAMATION);
         return;
      }
      BOOL bGroup;
      UINT group = GetDlgItemInt(IDC_LOC_GROUP, &bGroup, FALSE);
      loc["parse"] = { { "regex", CStringToUtf8(strRegex) }, { "group", bGroup ? group : 1 } };
   }
   else
   {
      loc.erase("parse");
   }

   CString strListIndex = GetText(this, IDC_LOC_LIST_INDEX).Trim();
   if (strListIndex.IsEmpty())
      loc.erase("list_index");
   else
      loc["list_index"] = (unsigned)_ttoi(strListIndex);

   CString strTypes = GetText(this, IDC_LOC_VALUE_TYPES);
   ordered_json types = ordered_json::array();
   int pos = 0;
   for (CString type = strTypes.Tokenize(_T(", "), pos); !type.IsEmpty(); type = strTypes.Tokenize(_T(", "), pos))
      types.push_back(CStringToUtf8(type.MakeUpper()));
   if (types.empty())
      loc.erase("value_types");
   else
      loc["value_types"] = types;

   if (has_map(m_Target.kind))
   {
      ordered_json map = ordered_json::object();
      CString strMap = GetText(this, IDC_LOC_MAP);
      pos = 0;
      for (CString line = strMap.Tokenize(_T("\r\n"), pos); !line.IsEmpty(); line = strMap.Tokenize(_T("\r\n"), pos))
      {
         int eq = line.ReverseFind(_T('='));
         if (eq <= 0)
         {
            AfxMessageBox(_T("Each map line is: model value = target value\n\n") + line, MB_OK | MB_ICONEXCLAMATION);
            return;
         }
         std::string key = CStringToUtf8(line.Left(eq).Trim());
         CString value = line.Mid(eq + 1).Trim();
         if (m_Target.kind == ValueKind::Boolean)
         {
            if (value.CompareNoCase(_T("true")) == 0)
               map[key] = true;
            else if (value.CompareNoCase(_T("false")) == 0)
               map[key] = false;
            else
            {
               AfxMessageBox(_T("A map value of this target is true or false:\n\n") + line, MB_OK | MB_ICONEXCLAMATION);
               return;
            }
         }
         else
         {
            map[key] = CStringToUtf8(value);
         }
      }
      if (map.empty())
         loc.erase("map");
      else
         loc["map"] = map;
   }

   set_or_erase(loc, "comment", CStringToUtf8(GetText(this, IDC_LOC_COMMENT).Trim()));

   m_Location = loc;
   CDialog::OnOK();
}

/////////////////////////////////////////////////////////////////////////////
// CMappingMessagesDlg

BEGIN_MESSAGE_MAP(CMappingMessagesDlg, CDialog)
END_MESSAGE_MAP()

CMappingMessagesDlg::CMappingMessagesDlg(const CString& strCaption, const CString& strText, CWnd* pParent)
   : CDialog(IDD_IMPORT_RESULTS, pParent), m_strCaption(strCaption), m_strText(strText)
{
}

BOOL CMappingMessagesDlg::OnInitDialog()
{
   CDialog::OnInitDialog();
   SetWindowText(m_strCaption);
   SetDlgItemText(IDC_EDIT, ToWindowsLines(m_strText));
   return TRUE;
}
