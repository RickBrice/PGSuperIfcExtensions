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

// Mapping table editor, stage 2: the export sections (property sets, quantity sets, classification systems, classifications).
//
// The tree lists this table's entries (bold) and the base table's entries that this table doesn't override. A base entry
// can be copied into this table to change it ("Change in this table"), or left out of the export with a "remove" entry.
// An extending table's entry replaces the base table's entry with the same identity:
//    property and quantity sets: name, element role, and attach
//    classification systems: name
//    classifications: element role and identification

#include "MappingEditorPanes.h"
#include "IfcMappingTable.h"

// The base table's entries that this table doesn't override
std::vector<const PropertySetDeclaration*> GetBaseSets(CMappingEditorDoc* pDoc, bool bQuantities, ElementKind role);
std::vector<const ClassificationSystemDeclaration*> GetBaseSystems(CMappingEditorDoc* pDoc);
std::vector<const ClassificationDeclaration*> GetBaseClassifications(CMappingEditorDoc* pDoc, ElementKind role);

// Keys of base entries in MappingNode::key
std::string BaseSetKey(const PropertySetDeclaration& pset);         // name|role|attach
std::string BaseClassificationKey(const ClassificationDeclaration& c); // role|identification

// "occurrence", "type", or "material"
const char* OwnerName(PropertyOwner owner);

// A section header (e.g. Property sets, or the property sets of one element role) with an Add button
class CSectionPane : public CMappingPane
{
public:
   CSectionPane(CMappingEditorDoc* pDoc, const MappingNode& node, const CString& text) : CMappingPane(IDD_MAPPING_SECTION_PANE, pDoc), m_Node(node), m_Text(text) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnAdd();
   afx_msg void OnSize(UINT nType, int cx, int cy);
   DECLARE_MESSAGE_MAP()

private:
   MappingNode m_Node;
   CString m_Text;
};

// A base table entry: what it is, and buttons to change it in this table or leave it out of the export
class CBaseItemPane : public CMappingPane
{
public:
   CBaseItemPane(CMappingEditorDoc* pDoc, const MappingNode& node) : CMappingPane(IDD_MAPPING_BASE_ITEM_PANE, pDoc), m_Node(node) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnChange();
   afx_msg void OnLeaveOut();
   afx_msg void OnSize(UINT nType, int cx, int cy);
   DECLARE_MESSAGE_MAP()

private:
   MappingNode m_Node;
};

// A property set or quantity set of this table
class CSetPane : public CMappingPane
{
public:
   CSetPane(CMappingEditorDoc* pDoc, bool bQuantities, int index) : CMappingPane(IDD_MAPPING_SET_PANE, pDoc), m_bQuantities(bQuantities), m_Index(index) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnFieldChanged();
   afx_msg void OnRolesChanged();
   afx_msg void OnAdd();
   afx_msg void OnEdit();
   afx_msg void OnRemove();
   afx_msg void OnUp();
   afx_msg void OnDown();
   afx_msg void OnDelete();
   afx_msg void OnPropertiesDblClk(NMHDR* pNMHDR, LRESULT* pResult);
   afx_msg void OnPropertiesChanged(NMHDR* pNMHDR, LRESULT* pResult);
   DECLARE_MESSAGE_MAP()

private:
   bool m_bQuantities;
   int m_Index;
   CListCtrl m_Properties;

   nlohmann::ordered_json& GetSet();
   const char* ItemsKey() const { return m_bQuantities ? "quantities" : "properties"; }
   void FillProperties(int select);
   void UpdateButtons();
   void SetProperties(const nlohmann::ordered_json& properties, int select);
};

// Edits one property (or quantity) of a set
class CPropertyDlg : public CDialog
{
public:
   CPropertyDlg(bool bQuantity, const std::vector<std::string>& roles, const nlohmann::ordered_json& property, CWnd* pParent = nullptr);

   nlohmann::ordered_json m_Property;

protected:
   BOOL OnInitDialog() override;
   void OnOK() override;
   afx_msg void OnBindingChanged();
   DECLARE_MESSAGE_MAP()

private:
   bool m_bQuantity;
   std::vector<std::string> m_Roles;
   void UpdateControls();
};

// A classification system of this table
class CSystemPane : public CMappingPane
{
public:
   CSystemPane(CMappingEditorDoc* pDoc, int index) : CMappingPane(IDD_MAPPING_SYSTEM_PANE, pDoc), m_Index(index) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnFieldChanged();
   afx_msg void OnDelete();
   DECLARE_MESSAGE_MAP()

private:
   int m_Index;
};

// A classification of this table
class CClassificationPane : public CMappingPane
{
public:
   CClassificationPane(CMappingEditorDoc* pDoc, int index) : CMappingPane(IDD_MAPPING_CLASSIFICATION_PANE, pDoc), m_Index(index) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnFieldChanged();
   afx_msg void OnRolesChanged();
   afx_msg void OnDelete();
   DECLARE_MESSAGE_MAP()

private:
   int m_Index;
};

// Stage 2 panes for a tree node, or nullptr if the node isn't an export section
std::unique_ptr<CMappingPane> CreateExportPane(const MappingNode& node, CMappingEditorDoc* pDoc);
