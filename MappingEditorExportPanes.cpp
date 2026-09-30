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

// MappingEditorExportPanes.cpp : the mapping table editor's export sections (stage 2)
//

#include "stdafx.h"
#include "resource.h"
#include "MappingEditorExportPanes.h"
#include "MappingEditorDoc.h"
#include "MappingEditorUtil.h"
#include "bSDD.h"

#include <sstream>

using namespace mapping_editor;

namespace
{
   const TCHAR* const VALUE_TYPES[] = { _T("IfcLabel"), _T("IfcText"), _T("IfcIdentifier"), _T("IfcBoolean"), _T("IfcInteger"), _T("IfcReal"),
      _T("IfcCountMeasure"), _T("IfcPressureMeasure"), _T("IfcLengthMeasure"), _T("IfcPositiveLengthMeasure"), _T("IfcNonNegativeLengthMeasure"),
      _T("IfcPlaneAngleMeasure"), _T("IfcPositivePlaneAngleMeasure"), _T("IfcRatioMeasure"), _T("IfcPositiveRatioMeasure"), _T("IfcAreaMeasure"),
      _T("IfcVolumeMeasure"), _T("IfcMassMeasure"), _T("IfcForceMeasure"), _T("IfcDate") };
   const TCHAR* const QUANTITY_TYPES[] = { _T("IfcQuantityLength"), _T("IfcQuantityArea"), _T("IfcQuantityVolume"), _T("IfcQuantityWeight"), _T("IfcQuantityCount") };

   // "IFCLABEL" (as the loader keeps it) -> "IfcLabel" (as tables are written)
   std::string type_spelling(const std::string& type)
   {
      return IfcSpelling(type);
   }

   const char* section_key(MappingNode::Kind kind)
   {
      switch (kind)
      {
      case MappingNode::Kind::PropertySets: case MappingNode::Kind::PropertySetGroup: case MappingNode::Kind::PropertySet: return "property_sets";
      case MappingNode::Kind::QuantitySets: case MappingNode::Kind::QuantitySetGroup: case MappingNode::Kind::QuantitySet: return "quantity_sets";
      case MappingNode::Kind::ClassificationSystems: case MappingNode::Kind::ClassificationSystem: return "classification_systems";
      default: return "classifications";
      }
   }

   ElementKind role_kind(const std::string& role)
   {
      ElementKind kind = ElementKind::Bridge;
      GetElementKind(role, kind);
      return kind;
   }

   bool has_role(const ordered_json& item, const std::string& role)
   {
      auto roles = roles_of(item);
      return std::find(roles.begin(), roles.end(), role) != roles.end();
   }

   const char* condition_name(PropertySetDeclaration::Condition condition)
   {
      switch (condition)
      {
      case PropertySetDeclaration::Condition::Quantities: return "quantities";
      case PropertySetDeclaration::Condition::Always: return "always";
      default: return "classify";
      }
   }

   // bSDD URIs as a table writes them: "bsdd" when the URI is built from the name, "bsdd:Other" for another bSDD name
   std::string short_uri(const std::string& uri, const std::string& name, const char* kind)
   {
      std::string prefix = BSDD_URI + kind + "/";
      if (uri == prefix + name)
         return "bsdd";
      if (uri.compare(0, prefix.size(), prefix) == 0)
         return "bsdd:" + uri.substr(prefix.size());
      return uri;
   }

   ordered_json value_json(const TargetValue& value)
   {
      return std::visit([](const auto& v) { return ordered_json(v); }, value);
   }

   // A base table's property set as an entry of this table, for the given role
   ordered_json set_json(const PropertySetDeclaration& pset)
   {
      ordered_json j = ordered_json::object();
      j["name"] = pset.name;
      j["applies_to"] = std::string(GetElementRoleName(pset.applies_to));
      if (!pset.quantities && pset.attach != PropertyOwner::Occurrence)
         j["attach"] = OwnerName(pset.attach);
      if (pset.condition != PropertySetDeclaration::Condition::Classify)
         j["condition"] = condition_name(pset.condition);
      if (pset.quantities && !pset.method.empty())
         j["method"] = pset.method;
      if (!pset.uri.empty())
         j["uri"] = short_uri(pset.uri, pset.name, "class");
      if (pset.shared)
         j["shared"] = true;

      ordered_json items = ordered_json::array();
      for (const auto& p : pset.properties)
      {
         ordered_json jp = ordered_json::object();
         jp["name"] = p.name;
         jp["type"] = type_spelling(p.type);
         if (p.target)
            jp["target"] = std::string(p.target->name);
         if (p.value)
            jp["value"] = value_json(*p.value);
         if (p.target && p.import == pset.quantities) // the default is import for properties, not for quantities
            jp["import"] = p.import;
         if (!p.uri.empty())
            jp["uri"] = short_uri(p.uri, p.name, "prop");
         if (!p.enumeration_name.empty())
            jp["enumeration"] = { { "name", p.enumeration_name }, { "values", p.enumeration_values } };
         items.push_back(jp);
      }
      j[pset.quantities ? "quantities" : "properties"] = items;
      return j;
   }

   ordered_json system_json(const ClassificationSystemDeclaration& system)
   {
      ordered_json j = ordered_json::object();
      j["name"] = system.name;
      if (!system.source.empty()) j["source"] = system.source;
      if (!system.edition.empty()) j["edition"] = system.edition;
      if (!system.edition_date.empty()) j["edition_date"] = system.edition_date;
      if (!system.specification.empty()) j["specification"] = (system.specification == BSDD_URI ? std::string("bsdd") : system.specification);
      return j;
   }

   ordered_json classification_json(const ClassificationDeclaration& c)
   {
      ordered_json j = ordered_json::object();
      j["applies_to"] = std::string(GetElementRoleName(c.applies_to));
      if (c.condition != PropertySetDeclaration::Condition::Classify)
         j["condition"] = condition_name(c.condition);
      j["system"] = c.system;
      j["identification"] = c.identification;
      if (!c.name.empty()) j["name"] = c.name;
      if (!c.location.empty())
         j["location"] = (c.location == BSDD_URI + "class/" + c.identification ? std::string("bsdd") : c.location);
      return j;
   }

   const PropertySetDeclaration* find_base_set(CMappingEditorDoc* pDoc, bool bQuantities, const std::string& key)
   {
      for (int kind = (int)ElementKind::Project; kind <= (int)ElementKind::Barrier; kind++)
      {
         for (const auto* pset : GetBaseSets(pDoc, bQuantities, (ElementKind)kind))
         {
            if (BaseSetKey(*pset) == key)
               return pset;
         }
      }
      return nullptr;
   }

   const ClassificationSystemDeclaration* find_base_system(CMappingEditorDoc* pDoc, const std::string& name)
   {
      for (const auto* system : GetBaseSystems(pDoc))
      {
         if (system->name == name)
            return system;
      }
      return nullptr;
   }

   const ClassificationDeclaration* find_base_classification(CMappingEditorDoc* pDoc, const std::string& key)
   {
      for (int kind = (int)ElementKind::Project; kind <= (int)ElementKind::Barrier; kind++)
      {
         for (const auto* c : GetBaseClassifications(pDoc, (ElementKind)kind))
         {
            if (BaseClassificationKey(*c) == key)
               return c;
         }
      }
      return nullptr;
   }

   std::string describe_properties(const ordered_json& items)
   {
      std::ostringstream os;
      for (const auto& p : items)
      {
         os << "      " << string_value(p, "name") << "   " << string_value(p, "type");
         if (p.contains("target")) os << "   target " << string_value(p, "target");
         if (p.contains("value")) os << "   value " << inline_json(p["value"]);
         if (p.contains("import")) os << (p["import"] == true ? "   imported" : "   export only");
         os << std::endl;
      }
      return os.str();
   }

   // Appends an entry to a section of the table and returns its index
   int append_entry(ordered_json& table, const char* key, const ordered_json& entry)
   {
      if (!table.contains(key) || !table[key].is_array())
         table[key] = ordered_json::array();
      table[key].push_back(entry);
      return (int)table[key].size() - 1;
   }

   // Removes an entry from a section, and the section when it's empty
   void erase_entry(ordered_json& table, const char* key, int index)
   {
      if (!table.contains(key) || !table[key].is_array() || index < 0 || (int)table[key].size() <= index)
         return;
      table[key].erase(table[key].begin() + index);
      if (table[key].empty())
         table.erase(key);
   }

   // The known classification system names: this table's and the base table's
   std::vector<std::string> system_names(CMappingEditorDoc* pDoc)
   {
      std::vector<std::string> names;
      const auto& table = pDoc->GetTable();
      if (table.contains("classification_systems") && table["classification_systems"].is_array())
      {
         for (const auto& system : table["classification_systems"])
            names.push_back(string_value(system, "name"));
      }
      if (const auto* pBase = pDoc->GetBaseTable())
      {
         for (const auto* system : pBase->GetClassificationSystems())
         {
            if (std::find(names.begin(), names.end(), system->name) == names.end())
               names.push_back(system->name);
         }
      }
      return names;
   }

   void fill_condition(CComboBox* pCombo, const std::string& condition)
   {
      for (auto name : { _T("classify"), _T("quantities"), _T("always") })
         pCombo->AddString(name);
      pCombo->SelectString(-1, condition.empty() ? CString(_T("classify")) : Utf8ToCString(condition));
   }

   // The text of a combo box (the selected item's, or what was typed)
   CString combo_text(CWnd* pWnd, int nID)
   {
      CString text;
      pWnd->GetDlgItem(nID)->GetWindowText(text);
      return text.Trim();
   }

   void resize_to_width(CDialog* pPane, int nID, int cx)
   {
      if (CWnd* pControl = pPane->GetDlgItem(nID))
      {
         CRect rect;
         pControl->GetWindowRect(&rect);
         pPane->ScreenToClient(&rect);
         CRect margin(0, 0, 7, 7);
         pPane->MapDialogRect(&margin);
         pControl->SetWindowPos(nullptr, 0, 0, std::max(50, (int)(cx - rect.left - margin.right)), rect.Height(), SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
      }
   }
}

const char* OwnerName(PropertyOwner owner)
{
   switch (owner)
   {
   case PropertyOwner::Type: return "type";
   case PropertyOwner::Material: return "material";
   default: return "occurrence";
   }
}

std::string BaseSetKey(const PropertySetDeclaration& pset)
{
   return pset.name + "|" + std::string(GetElementRoleName(pset.applies_to)) + "|" + OwnerName(pset.attach);
}

std::string BaseClassificationKey(const ClassificationDeclaration& c)
{
   return std::string(GetElementRoleName(c.applies_to)) + "|" + c.identification;
}

std::vector<const PropertySetDeclaration*> GetBaseSets(CMappingEditorDoc* pDoc, bool bQuantities, ElementKind role)
{
   std::vector<const PropertySetDeclaration*> sets;
   const auto* pBase = pDoc->GetBaseTable();
   if (!pBase)
      return sets;

   const auto& table = pDoc->GetTable();
   const char* key = bQuantities ? "quantity_sets" : "property_sets";
   std::string role_name(GetElementRoleName(role));
   for (auto attach : { PropertyOwner::Occurrence, PropertyOwner::Type, PropertyOwner::Material })
   {
      for (const auto* pset : pBase->GetPropertySets(role, attach))
      {
         if (pset->quantities != bQuantities)
            continue;

         // overridden (changed or left out) by an entry of this table
         bool bOverridden = false;
         if (table.contains(key) && table[key].is_array())
         {
            for (const auto& entry : table[key])
            {
               std::string entry_attach = string_value(entry, "attach");
               bOverridden |= (string_value(entry, "name") == pset->name && has_role(entry, role_name) &&
                  (entry_attach.empty() ? std::string("occurrence") : entry_attach) == OwnerName(attach));
            }
         }
         if (!bOverridden)
            sets.push_back(pset);
      }
   }
   return sets;
}

std::vector<const ClassificationSystemDeclaration*> GetBaseSystems(CMappingEditorDoc* pDoc)
{
   std::vector<const ClassificationSystemDeclaration*> systems;
   const auto* pBase = pDoc->GetBaseTable();
   if (!pBase)
      return systems;

   const auto& table = pDoc->GetTable();
   for (const auto* system : pBase->GetClassificationSystems())
   {
      bool bOverridden = false;
      if (table.contains("classification_systems") && table["classification_systems"].is_array())
      {
         for (const auto& entry : table["classification_systems"])
            bOverridden |= (string_value(entry, "name") == system->name);
      }
      if (!bOverridden)
         systems.push_back(system);
   }
   return systems;
}

std::vector<const ClassificationDeclaration*> GetBaseClassifications(CMappingEditorDoc* pDoc, ElementKind role)
{
   std::vector<const ClassificationDeclaration*> classifications;
   const auto* pBase = pDoc->GetBaseTable();
   if (!pBase)
      return classifications;

   const auto& table = pDoc->GetTable();
   std::string role_name(GetElementRoleName(role));
   for (const auto* c : pBase->GetClassifications(role))
   {
      bool bOverridden = false;
      if (table.contains("classifications") && table["classifications"].is_array())
      {
         for (const auto& entry : table["classifications"])
            bOverridden |= (string_value(entry, "identification") == c->identification && has_role(entry, role_name));
      }
      if (!bOverridden)
         classifications.push_back(c);
   }
   return classifications;
}

std::unique_ptr<CMappingPane> CreateExportPane(const MappingNode& node, CMappingEditorDoc* pDoc)
{
   const CString bold(_T("\n\nEntries in bold are in this table. The others come from the base table: select one to change it in this table, or to leave it out of the export."));
   switch (node.kind)
   {
   case MappingNode::Kind::PropertySets:
      return std::make_unique<CSectionPane>(pDoc, node, _T("Property sets exported for each element role, and the targets bound to their properties. A property bound to a target is also where an import reads it, unless \"import\" is off.") + bold);
   case MappingNode::Kind::QuantitySets:
      return std::make_unique<CSectionPane>(pDoc, node, _T("Quantity sets exported for each element role (when the export includes quantities).") + bold);
   case MappingNode::Kind::PropertySetGroup:
   case MappingNode::Kind::QuantitySetGroup:
      return std::make_unique<CSectionPane>(pDoc, node, Utf8ToCString("The " + std::string(node.kind == MappingNode::Kind::PropertySetGroup ? "property" : "quantity") + " sets of the " + node.key + " elements.") + bold);
   case MappingNode::Kind::ClassificationSystems:
      return std::make_unique<CSectionPane>(pDoc, node, _T("Classification systems (IfcClassification) written to the project. A classification refers to its system by name.") + bold);
   case MappingNode::Kind::Classifications:
      return std::make_unique<CSectionPane>(pDoc, node, _T("Classification references of each element role.") + bold);
   case MappingNode::Kind::ClassificationGroup:
      return std::make_unique<CSectionPane>(pDoc, node, Utf8ToCString("The classifications of the " + node.key + " elements.") + bold);

   case MappingNode::Kind::PropertySet:
   case MappingNode::Kind::QuantitySet:
      if (0 <= node.index)
         return std::make_unique<CSetPane>(pDoc, node.kind == MappingNode::Kind::QuantitySet, node.index);
      return std::make_unique<CBaseItemPane>(pDoc, node);

   case MappingNode::Kind::ClassificationSystem:
      if (0 <= node.index)
         return std::make_unique<CSystemPane>(pDoc, node.index);
      return std::make_unique<CBaseItemPane>(pDoc, node);

   case MappingNode::Kind::Classification:
      if (0 <= node.index)
         return std::make_unique<CClassificationPane>(pDoc, node.index);
      return std::make_unique<CBaseItemPane>(pDoc, node);

   default:
      return nullptr;
   }
}

/////////////////////////////////////////////////////////////////////////////
// CSectionPane

BEGIN_MESSAGE_MAP(CSectionPane, CMappingPane)
   ON_BN_CLICKED(IDC_SECTION_ADD, &CSectionPane::OnAdd)
   ON_WM_SIZE()
END_MESSAGE_MAP()

BOOL CSectionPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();
   SetDlgItemText(IDC_SECTION_TEXT, ToWindowsLines(m_Text));

   CString strAdd;
   switch (m_Node.kind)
   {
   case MappingNode::Kind::PropertySets: case MappingNode::Kind::PropertySetGroup: strAdd = _T("Add Property Set"); break;
   case MappingNode::Kind::QuantitySets: case MappingNode::Kind::QuantitySetGroup: strAdd = _T("Add Quantity Set"); break;
   case MappingNode::Kind::ClassificationSystems: strAdd = _T("Add System"); break;
   default: strAdd = _T("Add Classification"); break;
   }
   SetDlgItemText(IDC_SECTION_ADD, strAdd);
   return TRUE;
}

void CSectionPane::OnSize(UINT nType, int cx, int cy)
{
   CMappingPane::OnSize(nType, cx, cy);
   resize_to_width(this, IDC_SECTION_TEXT, cx);
}

void CSectionPane::OnAdd()
{
   auto& table = m_pDoc->GetTable();
   std::string role = m_Node.key.empty() ? std::string("girder") : m_Node.key;
   const char* key = section_key(m_Node.kind);

   ordered_json entry;
   MappingNode::Kind item_kind;
   switch (m_Node.kind)
   {
   case MappingNode::Kind::PropertySets: case MappingNode::Kind::PropertySetGroup:
      entry = { { "name", "NewPropertySet" }, { "applies_to", role }, { "properties", ordered_json::array({ { { "name", "NewProperty" }, { "type", "IfcLabel" } } }) } };
      item_kind = MappingNode::Kind::PropertySet;
      break;
   case MappingNode::Kind::QuantitySets: case MappingNode::Kind::QuantitySetGroup:
      entry = { { "name", "NewQuantitySet" }, { "applies_to", role }, { "quantities", ordered_json::array({ { { "name", "NewQuantity" }, { "type", "IfcQuantityLength" } } }) } };
      item_kind = MappingNode::Kind::QuantitySet;
      break;
   case MappingNode::Kind::ClassificationSystems:
      entry = { { "name", "NewSystem" } };
      item_kind = MappingNode::Kind::ClassificationSystem;
      break;
   default:
   {
      auto systems = system_names(m_pDoc);
      entry = { { "applies_to", role }, { "system", systems.empty() ? std::string() : systems.front() }, { "identification", "NewClassification" } };
      item_kind = MappingNode::Kind::Classification;
      break;
   }
   }

   int index = append_entry(table, key, entry);
   m_pDoc->TableChanged(true);
   m_pDoc->SelectNode(MappingNode{ item_kind, "", index });
}

/////////////////////////////////////////////////////////////////////////////
// CBaseItemPane

BEGIN_MESSAGE_MAP(CBaseItemPane, CMappingPane)
   ON_BN_CLICKED(IDC_BASE_CHANGE, &CBaseItemPane::OnChange)
   ON_BN_CLICKED(IDC_BASE_LEAVE_OUT, &CBaseItemPane::OnLeaveOut)
   ON_WM_SIZE()
END_MESSAGE_MAP()

BOOL CBaseItemPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();

   std::ostringstream os;
   std::string table_name = m_pDoc->GetBaseTable() ? m_pDoc->GetBaseTable()->GetFiles().front().name : std::string();
   bool bLeaveOut = true;
   switch (m_Node.kind)
   {
   case MappingNode::Kind::PropertySet:
   case MappingNode::Kind::QuantitySet:
      if (const auto* pset = find_base_set(m_pDoc, m_Node.kind == MappingNode::Kind::QuantitySet, m_Node.key))
      {
         auto j = set_json(*pset);
         os << (pset->quantities ? "Quantity set " : "Property set ") << pset->name << ", from the base table" << std::endl << std::endl;
         os << "Applies to: " << GetElementRoleName(pset->applies_to) << std::endl;
         if (!pset->quantities) os << "Attach to: " << OwnerName(pset->attach) << std::endl;
         os << "Condition: " << condition_name(pset->condition) << std::endl;
         if (!pset->method.empty()) os << "Method: " << pset->method << std::endl;
         if (!pset->uri.empty()) os << "URI: " << pset->uri << std::endl;
         if (pset->shared) os << "Shared by all the elements of the role" << std::endl;
         os << (pset->quantities ? "Quantities:" : "Properties:") << std::endl << describe_properties(j[pset->quantities ? "quantities" : "properties"]);
      }
      break;
   case MappingNode::Kind::ClassificationSystem:
      bLeaveOut = false; // a system that no classification uses isn't written anyway
      if (const auto* system = find_base_system(m_pDoc, m_Node.key))
      {
         os << "Classification system " << system->name << ", from the base table" << std::endl << std::endl;
         if (!system->source.empty()) os << "Source: " << system->source << std::endl;
         if (!system->edition.empty()) os << "Edition: " << system->edition << std::endl;
         if (!system->edition_date.empty()) os << "Edition date: " << system->edition_date << std::endl;
         if (!system->specification.empty()) os << "Specification: " << system->specification << std::endl;
      }
      break;
   case MappingNode::Kind::Classification:
      if (const auto* c = find_base_classification(m_pDoc, m_Node.key))
      {
         os << "Classification " << c->identification << ", from the base table" << std::endl << std::endl;
         os << "Applies to: " << GetElementRoleName(c->applies_to) << std::endl;
         os << "Condition: " << condition_name(c->condition) << std::endl;
         os << "System: " << c->system << std::endl;
         if (!c->name.empty()) os << "Name: " << c->name << std::endl;
         if (!c->location.empty()) os << "Location: " << c->location << std::endl;
      }
      break;
   default:
      break;
   }
   if (!table_name.empty())
      os << std::endl << "Base table: " << table_name;

   SetDlgItemText(IDC_BASE_TEXT, ToWindowsLines(Utf8ToCString(os.str())));
   GetDlgItem(IDC_BASE_LEAVE_OUT)->ShowWindow(bLeaveOut ? SW_SHOW : SW_HIDE);
   return TRUE;
}

void CBaseItemPane::OnSize(UINT nType, int cx, int cy)
{
   CMappingPane::OnSize(nType, cx, cy);
   resize_to_width(this, IDC_BASE_TEXT, cx);
}

void CBaseItemPane::OnChange()
{
   auto& table = m_pDoc->GetTable();
   const char* key = section_key(m_Node.kind);
   ordered_json entry;
   switch (m_Node.kind)
   {
   case MappingNode::Kind::PropertySet:
   case MappingNode::Kind::QuantitySet:
      if (const auto* pset = find_base_set(m_pDoc, m_Node.kind == MappingNode::Kind::QuantitySet, m_Node.key))
         entry = set_json(*pset);
      break;
   case MappingNode::Kind::ClassificationSystem:
      if (const auto* system = find_base_system(m_pDoc, m_Node.key))
         entry = system_json(*system);
      break;
   case MappingNode::Kind::Classification:
      if (const auto* c = find_base_classification(m_pDoc, m_Node.key))
         entry = classification_json(*c);
      break;
   default:
      break;
   }
   if (entry.is_null())
      return;

   std::string role = entry.contains("applies_to") ? entry["applies_to"].get<std::string>() : std::string();
   int index = append_entry(table, key, entry);
   m_pDoc->TableChanged(true);
   m_pDoc->SelectNode(MappingNode{ m_Node.kind, role, index });
}

void CBaseItemPane::OnLeaveOut()
{
   auto& table = m_pDoc->GetTable();
   const char* key = section_key(m_Node.kind);
   ordered_json entry;
   std::string role;
   if (m_Node.kind == MappingNode::Kind::PropertySet || m_Node.kind == MappingNode::Kind::QuantitySet)
   {
      if (const auto* pset = find_base_set(m_pDoc, m_Node.kind == MappingNode::Kind::QuantitySet, m_Node.key))
      {
         role = GetElementRoleName(pset->applies_to);
         entry = { { "name", pset->name }, { "applies_to", role } };
         if (!pset->quantities && pset->attach != PropertyOwner::Occurrence)
            entry["attach"] = OwnerName(pset->attach);
         entry["remove"] = true;
      }
   }
   else if (m_Node.kind == MappingNode::Kind::Classification)
   {
      if (const auto* c = find_base_classification(m_pDoc, m_Node.key))
      {
         role = GetElementRoleName(c->applies_to);
         entry = { { "applies_to", role }, { "identification", c->identification }, { "remove", true } };
      }
   }
   if (entry.is_null())
      return;

   int index = append_entry(table, key, entry);
   m_pDoc->TableChanged(true);
   m_pDoc->SelectNode(MappingNode{ m_Node.kind, role, index });
}

/////////////////////////////////////////////////////////////////////////////
// CSetPane

BEGIN_MESSAGE_MAP(CSetPane, CMappingPane)
   ON_EN_KILLFOCUS(IDC_SET_NAME, &CSetPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SET_METHOD, &CSetPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SET_URI, &CSetPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SET_COMMENT, &CSetPane::OnFieldChanged)
   ON_CBN_SELCHANGE(IDC_SET_ATTACH, &CSetPane::OnFieldChanged)
   ON_CBN_SELCHANGE(IDC_SET_CONDITION, &CSetPane::OnFieldChanged)
   ON_BN_CLICKED(IDC_SET_SHARED, &CSetPane::OnFieldChanged)
   ON_LBN_SELCHANGE(IDC_SET_ROLES, &CSetPane::OnRolesChanged)
   ON_BN_CLICKED(IDC_PROP_ADD, &CSetPane::OnAdd)
   ON_BN_CLICKED(IDC_PROP_EDIT, &CSetPane::OnEdit)
   ON_BN_CLICKED(IDC_PROP_REMOVE, &CSetPane::OnRemove)
   ON_BN_CLICKED(IDC_PROP_UP, &CSetPane::OnUp)
   ON_BN_CLICKED(IDC_PROP_DOWN, &CSetPane::OnDown)
   ON_BN_CLICKED(IDC_SET_DELETE, &CSetPane::OnDelete)
   ON_NOTIFY(NM_DBLCLK, IDC_SET_PROPERTIES, &CSetPane::OnPropertiesDblClk)
   ON_NOTIFY(LVN_ITEMCHANGED, IDC_SET_PROPERTIES, &CSetPane::OnPropertiesChanged)
END_MESSAGE_MAP()

ordered_json& CSetPane::GetSet()
{
   return m_pDoc->GetTable()[m_bQuantities ? "quantity_sets" : "property_sets"][m_Index];
}

BOOL CSetPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();

   const auto& set = GetSet();
   bool bRemove = set.value("remove", false);

   SetDlgItemText(IDC_SET_NAME, Utf8ToCString(string_value(set, "name")));
   FillRoleList((CListBox*)GetDlgItem(IDC_SET_ROLES), roles_of(set));

   if (m_bQuantities)
   {
      SetDlgItemText(IDC_SET_ATTACH_LABEL, _T("Method:"));
      GetDlgItem(IDC_SET_ATTACH)->ShowWindow(SW_HIDE);
      GetDlgItem(IDC_SET_SHARED)->ShowWindow(SW_HIDE);
      SetDlgItemText(IDC_SET_METHOD, Utf8ToCString(string_value(set, "method")));
      SetDlgItemText(IDC_SET_ITEMS_LABEL, _T("Quantities:"));
   }
   else
   {
      GetDlgItem(IDC_SET_METHOD)->ShowWindow(SW_HIDE);
      CComboBox* pAttach = (CComboBox*)GetDlgItem(IDC_SET_ATTACH);
      for (auto name : { _T("occurrence"), _T("type"), _T("material") })
         pAttach->AddString(name);
      auto attach = string_value(set, "attach");
      pAttach->SelectString(-1, attach.empty() ? CString(_T("occurrence")) : Utf8ToCString(attach));
      CheckDlgButton(IDC_SET_SHARED, set.value("shared", false) ? BST_CHECKED : BST_UNCHECKED);
   }
   fill_condition((CComboBox*)GetDlgItem(IDC_SET_CONDITION), string_value(set, "condition"));
   SetDlgItemText(IDC_SET_URI, Utf8ToCString(string_value(set, "uri")));
   SetDlgItemText(IDC_SET_COMMENT, Utf8ToCString(string_value(set, "comment")));

   m_Properties.SubclassDlgItem(IDC_SET_PROPERTIES, this);
   m_Properties.SetExtendedStyle(LVS_EX_FULLROWSELECT | LVS_EX_GRIDLINES);
   CRect rect;
   m_Properties.GetClientRect(&rect);
   m_Properties.InsertColumn(0, _T("Name"), LVCFMT_LEFT, rect.Width() * 3 / 10);
   m_Properties.InsertColumn(1, _T("Type"), LVCFMT_LEFT, rect.Width() * 3 / 10);
   m_Properties.InsertColumn(2, _T("Target or value"), LVCFMT_LEFT, rect.Width() * 4 / 10);

   if (bRemove)
   {
      for (int nID : { IDC_SET_NAME, IDC_SET_ROLES, IDC_SET_ATTACH, IDC_SET_METHOD, IDC_SET_CONDITION, IDC_SET_URI, IDC_SET_SHARED, IDC_SET_COMMENT, IDC_SET_PROPERTIES, IDC_PROP_ADD })
         GetDlgItem(nID)->EnableWindow(FALSE);
      SetDlgItemText(IDC_SET_NOTE, _T("This entry leaves the base table's set with this name, element role, and attach out of the export. Delete it to export the base table's set again."));
   }
   else
   {
      SetDlgItemText(IDC_SET_NOTE, _T("URI: \"bsdd\" (built from the name), \"bsdd:Name\", or a URI. A shared set is one instance for all the elements of the role; its properties can't have targets."));
   }

   FillProperties(-1);
   return TRUE;
}

void CSetPane::FillProperties(int select)
{
   m_Properties.DeleteAllItems();
   const auto& set = GetSet();
   if (set.contains(ItemsKey()) && set[ItemsKey()].is_array())
   {
      int i = 0;
      for (const auto& p : set[ItemsKey()])
      {
         m_Properties.InsertItem(i, Utf8ToCString(string_value(p, "name")));
         m_Properties.SetItemText(i, 1, Utf8ToCString(string_value(p, "type")));
         std::string binding;
         if (p.contains("target"))
            binding = string_value(p, "target") + (p.contains("import") ? (p["import"] == true ? " (imported)" : " (export only)") : "");
         else if (p.contains("value"))
            binding = "= " + inline_json(p["value"]);
         if (p.contains("enumeration"))
            binding += " (enumeration)";
         m_Properties.SetItemText(i, 2, Utf8ToCString(binding));
         i++;
      }
   }
   if (0 <= select && select < m_Properties.GetItemCount())
   {
      m_Properties.SetItemState(select, LVIS_SELECTED | LVIS_FOCUSED, LVIS_SELECTED | LVIS_FOCUSED);
      m_Properties.EnsureVisible(select, FALSE);
   }
   UpdateButtons();
}

void CSetPane::UpdateButtons()
{
   bool bRemove = GetSet().value("remove", false);
   int sel = m_Properties.GetNextItem(-1, LVNI_SELECTED);
   int count = m_Properties.GetItemCount();
   GetDlgItem(IDC_PROP_EDIT)->EnableWindow(!bRemove && 0 <= sel);
   GetDlgItem(IDC_PROP_REMOVE)->EnableWindow(!bRemove && 0 <= sel);
   GetDlgItem(IDC_PROP_UP)->EnableWindow(!bRemove && 0 < sel);
   GetDlgItem(IDC_PROP_DOWN)->EnableWindow(!bRemove && 0 <= sel && sel < count - 1);
}

void CSetPane::OnPropertiesChanged(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   UpdateButtons();
}

void CSetPane::OnPropertiesDblClk(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   OnEdit();
}

void CSetPane::OnFieldChanged()
{
   auto& set = GetSet();
   if (set.value("remove", false))
      return;

   ordered_json before = set;
   set["name"] = CStringToUtf8(GetText(this, IDC_SET_NAME).Trim());
   if (m_bQuantities)
   {
      set_or_erase(set, "method", CStringToUtf8(GetText(this, IDC_SET_METHOD).Trim()));
   }
   else
   {
      CString attach = combo_text(this, IDC_SET_ATTACH);
      set_or_erase(set, "attach", attach == _T("occurrence") ? std::string() : CStringToUtf8(attach));
      if (IsDlgButtonChecked(IDC_SET_SHARED) == BST_CHECKED)
         set["shared"] = true;
      else
         set.erase("shared");
   }
   CString condition = combo_text(this, IDC_SET_CONDITION);
   set_or_erase(set, "condition", condition == _T("classify") ? std::string() : CStringToUtf8(condition));
   set_or_erase(set, "uri", CStringToUtf8(GetText(this, IDC_SET_URI).Trim()));
   set_or_erase(set, "comment", CStringToUtf8(GetText(this, IDC_SET_COMMENT).Trim()));

   if (set != before)
      m_pDoc->TableChanged(true); // the name and attach show in the tree
}

void CSetPane::OnRolesChanged()
{
   auto& set = GetSet();
   auto roles = GetRoleList((CListBox*)GetDlgItem(IDC_SET_ROLES));
   if (set["applies_to"] != roles)
   {
      set["applies_to"] = roles;
      m_pDoc->TableChanged(true);
   }
}

void CSetPane::SetProperties(const ordered_json& properties, int select)
{
   GetSet()[ItemsKey()] = properties;
   m_pDoc->TableChanged(false);
   FillProperties(select);
}

void CSetPane::OnAdd()
{
   CPropertyDlg dlg(m_bQuantities, roles_of(GetSet()), ordered_json::object(), this);
   if (dlg.DoModal() == IDOK)
   {
      auto items = GetSet().contains(ItemsKey()) ? GetSet()[ItemsKey()] : ordered_json::array();
      items.push_back(dlg.m_Property);
      SetProperties(items, (int)items.size() - 1);
   }
}

void CSetPane::OnEdit()
{
   int sel = m_Properties.GetNextItem(-1, LVNI_SELECTED);
   auto items = GetSet()[ItemsKey()];
   if (sel < 0 || (int)items.size() <= sel)
      return;

   CPropertyDlg dlg(m_bQuantities, roles_of(GetSet()), items[sel], this);
   if (dlg.DoModal() == IDOK && dlg.m_Property != items[sel])
   {
      items[sel] = dlg.m_Property;
      SetProperties(items, sel);
   }
}

void CSetPane::OnRemove()
{
   int sel = m_Properties.GetNextItem(-1, LVNI_SELECTED);
   auto items = GetSet()[ItemsKey()];
   if (sel < 0 || (int)items.size() <= sel)
      return;
   items.erase(items.begin() + sel);
   SetProperties(items, std::min(sel, (int)items.size() - 1));
}

void CSetPane::OnUp()
{
   int sel = m_Properties.GetNextItem(-1, LVNI_SELECTED);
   auto items = GetSet()[ItemsKey()];
   if (sel <= 0 || (int)items.size() <= sel)
      return;
   std::swap(items[sel - 1], items[sel]);
   SetProperties(items, sel - 1);
}

void CSetPane::OnDown()
{
   int sel = m_Properties.GetNextItem(-1, LVNI_SELECTED);
   auto items = GetSet()[ItemsKey()];
   if (sel < 0 || (int)items.size() - 1 <= sel)
      return;
   std::swap(items[sel], items[sel + 1]);
   SetProperties(items, sel + 1);
}

void CSetPane::OnDelete()
{
   erase_entry(m_pDoc->GetTable(), m_bQuantities ? "quantity_sets" : "property_sets", m_Index);
   m_pDoc->TableChanged(true);
   m_pDoc->SelectNode(MappingNode{ m_bQuantities ? MappingNode::Kind::QuantitySets : MappingNode::Kind::PropertySets, "" });
}

/////////////////////////////////////////////////////////////////////////////
// CPropertyDlg

BEGIN_MESSAGE_MAP(CPropertyDlg, CDialog)
   ON_BN_CLICKED(IDC_PROP_NONE, &CPropertyDlg::OnBindingChanged)
   ON_BN_CLICKED(IDC_PROP_BIND_TARGET, &CPropertyDlg::OnBindingChanged)
   ON_BN_CLICKED(IDC_PROP_BIND_VALUE, &CPropertyDlg::OnBindingChanged)
END_MESSAGE_MAP()

CPropertyDlg::CPropertyDlg(bool bQuantity, const std::vector<std::string>& roles, const ordered_json& property, CWnd* pParent)
   : CDialog(IDD_MAPPING_PROPERTY, pParent), m_Property(property), m_bQuantity(bQuantity), m_Roles(roles)
{
}

BOOL CPropertyDlg::OnInitDialog()
{
   CDialog::OnInitDialog();
   SetWindowText(m_bQuantity ? _T("Quantity") : _T("Property"));

   const auto& p = m_Property;
   SetDlgItemText(IDC_PROP_NAME, Utf8ToCString(string_value(p, "name")));

   CComboBox* pType = (CComboBox*)GetDlgItem(IDC_PROP_TYPE);
   if (m_bQuantity)
      for (auto type : QUANTITY_TYPES) pType->AddString(type);
   else
      for (auto type : VALUE_TYPES) pType->AddString(type);
   CString type = Utf8ToCString(string_value(p, "type"));
   pType->SetWindowText(type.IsEmpty() ? (m_bQuantity ? CString(_T("IfcQuantityLength")) : CString(_T("IfcLabel"))) : type);

   // the targets of the elements the set applies to
   CComboBox* pTarget = (CComboBox*)GetDlgItem(IDC_PROP_TARGET);
   std::string current = string_value(p, "target");
   for (const auto& target : GetTargetDefs())
   {
      bool bFits = false;
      for (const auto& role : m_Roles)
         bFits |= (target.element == GetTargetElement(role_kind(role)));
      if (bFits || target.name == current)
         pTarget->AddString(Utf8ToCString(std::string(target.name)));
   }
   if (!current.empty())
      pTarget->SelectString(-1, Utf8ToCString(current));

   int binding = p.contains("target") ? IDC_PROP_BIND_TARGET : (p.contains("value") ? IDC_PROP_BIND_VALUE : IDC_PROP_NONE);
   CheckRadioButton(IDC_PROP_NONE, IDC_PROP_BIND_VALUE, binding);
   if (p.contains("value"))
      SetDlgItemText(IDC_PROP_VALUE, Utf8ToCString(p["value"].is_string() ? p["value"].get<std::string>() : p["value"].dump()));
   CheckDlgButton(IDC_PROP_IMPORT, p.value("import", !m_bQuantity) ? BST_CHECKED : BST_UNCHECKED);

   SetDlgItemText(IDC_PROP_URI, Utf8ToCString(string_value(p, "uri")));
   if (p.contains("enumeration") && p["enumeration"].is_object())
   {
      SetDlgItemText(IDC_PROP_ENUM_NAME, Utf8ToCString(string_value(p["enumeration"], "name")));
      CString values;
      if (p["enumeration"].contains("values") && p["enumeration"]["values"].is_array())
      {
         for (const auto& v : p["enumeration"]["values"])
            values += Utf8ToCString(v.is_string() ? v.get<std::string>() : v.dump()) + _T("\r\n");
      }
      SetDlgItemText(IDC_PROP_ENUM_VALUES, values);
   }
   SetDlgItemText(IDC_PROP_COMMENT, Utf8ToCString(string_value(p, "comment")));

   for (int nID : { IDC_PROP_URI, IDC_PROP_ENUM_NAME, IDC_PROP_ENUM_VALUES })
      GetDlgItem(nID)->EnableWindow(!m_bQuantity);

   UpdateControls();
   return TRUE;
}

void CPropertyDlg::OnBindingChanged()
{
   UpdateControls();
}

void CPropertyDlg::UpdateControls()
{
   int binding = GetCheckedRadioButton(IDC_PROP_NONE, IDC_PROP_BIND_VALUE);
   GetDlgItem(IDC_PROP_TARGET)->EnableWindow(binding == IDC_PROP_BIND_TARGET);
   GetDlgItem(IDC_PROP_IMPORT)->EnableWindow(binding == IDC_PROP_BIND_TARGET);
   GetDlgItem(IDC_PROP_VALUE)->EnableWindow(binding == IDC_PROP_BIND_VALUE);
}

void CPropertyDlg::OnOK()
{
   ordered_json p = m_Property.is_object() ? m_Property : ordered_json::object();
   for (const char* key : { "target", "value", "import" })
      p.erase(key);

   CString strName = GetText(this, IDC_PROP_NAME).Trim();
   CString strType = combo_text(this, IDC_PROP_TYPE);
   if (strName.IsEmpty() || strType.IsEmpty())
   {
      AfxMessageBox(_T("Enter the name and the type."), MB_OK | MB_ICONEXCLAMATION);
      return;
   }
   p["name"] = CStringToUtf8(strName);
   p["type"] = CStringToUtf8(strType);

   int binding = GetCheckedRadioButton(IDC_PROP_NONE, IDC_PROP_BIND_VALUE);
   if (binding == IDC_PROP_BIND_TARGET)
   {
      CString strTarget = combo_text(this, IDC_PROP_TARGET);
      if (strTarget.IsEmpty())
      {
         AfxMessageBox(_T("Choose the target."), MB_OK | MB_ICONEXCLAMATION);
         return;
      }
      p["target"] = CStringToUtf8(strTarget);
      bool bImport = IsDlgButtonChecked(IDC_PROP_IMPORT) == BST_CHECKED;
      if (bImport == m_bQuantity) // the default is import for properties, not for quantities
         p["import"] = bImport;
   }
   else if (binding == IDC_PROP_BIND_VALUE)
   {
      // a constant of the property's type
      CString strValue = GetText(this, IDC_PROP_VALUE).Trim();
      CString upper(strType);
      upper.MakeUpper();
      if (m_bQuantity || (upper != _T("IFCLABEL") && upper != _T("IFCTEXT") && upper != _T("IFCIDENTIFIER") && upper != _T("IFCDATE") && upper != _T("IFCBOOLEAN")))
      {
         TCHAR* end = nullptr;
         double value = _tcstod(strValue, &end);
         if (strValue.IsEmpty() || (end && *end != 0))
         {
            AfxMessageBox(_T("The value of a ") + strType + _T(" is a number."), MB_OK | MB_ICONEXCLAMATION);
            return;
         }
         if (upper == _T("IFCINTEGER") || upper == _T("IFCCOUNTMEASURE"))
            p["value"] = (long long)value;
         else
            p["value"] = value;
      }
      else if (upper == _T("IFCBOOLEAN"))
      {
         if (strValue.CompareNoCase(_T("true")) == 0)
            p["value"] = true;
         else if (strValue.CompareNoCase(_T("false")) == 0)
            p["value"] = false;
         else
         {
            AfxMessageBox(_T("The value of an IfcBoolean is true or false."), MB_OK | MB_ICONEXCLAMATION);
            return;
         }
      }
      else
      {
         p["value"] = CStringToUtf8(strValue);
      }
   }

   if (!m_bQuantity)
   {
      set_or_erase(p, "uri", CStringToUtf8(GetText(this, IDC_PROP_URI).Trim()));
      CString strEnum = GetText(this, IDC_PROP_ENUM_NAME).Trim();
      auto lines = GetLines(this, IDC_PROP_ENUM_VALUES);
      if (strEnum.IsEmpty() != lines.empty())
      {
         AfxMessageBox(_T("An enumeration has a name and values (one per line)."), MB_OK | MB_ICONEXCLAMATION);
         return;
      }
      if (strEnum.IsEmpty())
      {
         p.erase("enumeration");
      }
      else
      {
         ordered_json values = ordered_json::array();
         for (const auto& line : lines)
            values.push_back(CStringToUtf8(line));
         p["enumeration"] = { { "name", CStringToUtf8(strEnum) }, { "values", values } };
      }
   }
   set_or_erase(p, "comment", CStringToUtf8(GetText(this, IDC_PROP_COMMENT).Trim()));

   m_Property = p;
   CDialog::OnOK();
}

/////////////////////////////////////////////////////////////////////////////
// CSystemPane

BEGIN_MESSAGE_MAP(CSystemPane, CMappingPane)
   ON_EN_KILLFOCUS(IDC_SYS_NAME, &CSystemPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SYS_SOURCE, &CSystemPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SYS_EDITION, &CSystemPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SYS_DATE, &CSystemPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SYS_SPEC, &CSystemPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_SYS_COMMENT, &CSystemPane::OnFieldChanged)
   ON_BN_CLICKED(IDC_SYS_DELETE, &CSystemPane::OnDelete)
END_MESSAGE_MAP()

BOOL CSystemPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();
   const auto& system = m_pDoc->GetTable()["classification_systems"][m_Index];
   SetDlgItemText(IDC_SYS_NAME, Utf8ToCString(string_value(system, "name")));
   SetDlgItemText(IDC_SYS_SOURCE, Utf8ToCString(string_value(system, "source")));
   SetDlgItemText(IDC_SYS_EDITION, Utf8ToCString(string_value(system, "edition")));
   SetDlgItemText(IDC_SYS_DATE, Utf8ToCString(string_value(system, "edition_date")));
   SetDlgItemText(IDC_SYS_SPEC, Utf8ToCString(string_value(system, "specification")));
   SetDlgItemText(IDC_SYS_COMMENT, Utf8ToCString(string_value(system, "comment")));
   return TRUE;
}

void CSystemPane::OnFieldChanged()
{
   auto& system = m_pDoc->GetTable()["classification_systems"][m_Index];
   ordered_json before = system;
   system["name"] = CStringToUtf8(GetText(this, IDC_SYS_NAME).Trim());
   set_or_erase(system, "source", CStringToUtf8(GetText(this, IDC_SYS_SOURCE).Trim()));
   set_or_erase(system, "edition", CStringToUtf8(GetText(this, IDC_SYS_EDITION).Trim()));
   set_or_erase(system, "edition_date", CStringToUtf8(GetText(this, IDC_SYS_DATE).Trim()));
   set_or_erase(system, "specification", CStringToUtf8(GetText(this, IDC_SYS_SPEC).Trim()));
   set_or_erase(system, "comment", CStringToUtf8(GetText(this, IDC_SYS_COMMENT).Trim()));
   if (system != before)
      m_pDoc->TableChanged(true);
}

void CSystemPane::OnDelete()
{
   erase_entry(m_pDoc->GetTable(), "classification_systems", m_Index);
   m_pDoc->TableChanged(true);
   m_pDoc->SelectNode(MappingNode{ MappingNode::Kind::ClassificationSystems, "" });
}

/////////////////////////////////////////////////////////////////////////////
// CClassificationPane

BEGIN_MESSAGE_MAP(CClassificationPane, CMappingPane)
   ON_LBN_SELCHANGE(IDC_CLS_ROLES, &CClassificationPane::OnRolesChanged)
   ON_CBN_SELCHANGE(IDC_CLS_CONDITION, &CClassificationPane::OnFieldChanged)
   ON_CBN_SELCHANGE(IDC_CLS_SYSTEM, &CClassificationPane::OnFieldChanged)
   ON_CBN_KILLFOCUS(IDC_CLS_SYSTEM, &CClassificationPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_CLS_IDENTIFICATION, &CClassificationPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_CLS_NAME, &CClassificationPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_CLS_LOCATION, &CClassificationPane::OnFieldChanged)
   ON_EN_KILLFOCUS(IDC_CLS_COMMENT, &CClassificationPane::OnFieldChanged)
   ON_BN_CLICKED(IDC_CLS_DELETE, &CClassificationPane::OnDelete)
END_MESSAGE_MAP()

BOOL CClassificationPane::OnInitDialog()
{
   CMappingPane::OnInitDialog();
   const auto& c = m_pDoc->GetTable()["classifications"][m_Index];
   bool bRemove = c.value("remove", false);

   FillRoleList((CListBox*)GetDlgItem(IDC_CLS_ROLES), roles_of(c));
   fill_condition((CComboBox*)GetDlgItem(IDC_CLS_CONDITION), string_value(c, "condition"));
   CComboBox* pSystem = (CComboBox*)GetDlgItem(IDC_CLS_SYSTEM);
   for (const auto& name : system_names(m_pDoc))
      pSystem->AddString(Utf8ToCString(name));
   pSystem->SetWindowText(Utf8ToCString(string_value(c, "system")));
   SetDlgItemText(IDC_CLS_IDENTIFICATION, Utf8ToCString(string_value(c, "identification")));
   SetDlgItemText(IDC_CLS_NAME, Utf8ToCString(string_value(c, "name")));
   SetDlgItemText(IDC_CLS_LOCATION, Utf8ToCString(string_value(c, "location")));
   SetDlgItemText(IDC_CLS_COMMENT, Utf8ToCString(string_value(c, "comment")));

   if (bRemove)
   {
      for (int nID : { IDC_CLS_ROLES, IDC_CLS_CONDITION, IDC_CLS_SYSTEM, IDC_CLS_IDENTIFICATION, IDC_CLS_NAME, IDC_CLS_LOCATION, IDC_CLS_COMMENT })
         GetDlgItem(nID)->EnableWindow(FALSE);
      SetDlgItemText(IDC_CLS_NOTE, _T("This entry leaves the base table's classification with this identification out of the export for these element roles. Delete it to export the base table's classification again."));
   }
   else
   {
      SetDlgItemText(IDC_CLS_NOTE, _T("Location: \"bsdd\" (built from the identification) or a URI. An entry with the identification of a base table's classification replaces it for these element roles."));
   }
   return TRUE;
}

void CClassificationPane::OnFieldChanged()
{
   auto& c = m_pDoc->GetTable()["classifications"][m_Index];
   if (c.value("remove", false))
      return;

   ordered_json before = c;
   CString condition = combo_text(this, IDC_CLS_CONDITION);
   set_or_erase(c, "condition", condition == _T("classify") ? std::string() : CStringToUtf8(condition));
   c["system"] = CStringToUtf8(combo_text(this, IDC_CLS_SYSTEM));
   c["identification"] = CStringToUtf8(GetText(this, IDC_CLS_IDENTIFICATION).Trim());
   set_or_erase(c, "name", CStringToUtf8(GetText(this, IDC_CLS_NAME).Trim()));
   set_or_erase(c, "location", CStringToUtf8(GetText(this, IDC_CLS_LOCATION).Trim()));
   set_or_erase(c, "comment", CStringToUtf8(GetText(this, IDC_CLS_COMMENT).Trim()));
   if (c != before)
      m_pDoc->TableChanged(true); // the identification shows in the tree
}

void CClassificationPane::OnRolesChanged()
{
   auto& c = m_pDoc->GetTable()["classifications"][m_Index];
   auto roles = GetRoleList((CListBox*)GetDlgItem(IDC_CLS_ROLES));
   if (c["applies_to"] != roles)
   {
      c["applies_to"] = roles;
      m_pDoc->TableChanged(true);
   }
}

void CClassificationPane::OnDelete()
{
   erase_entry(m_pDoc->GetTable(), "classifications", m_Index);
   m_pDoc->TableChanged(true);
   m_pDoc->SelectNode(MappingNode{ MappingNode::Kind::Classifications, "" });
}
