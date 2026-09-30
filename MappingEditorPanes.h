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

// Mapping table editor panes (the right side of the frame) and dialogs.
// A pane edits the document's JSON directly: a field is written to the table when it loses the focus.

#include <afxcmn.h>
#include "MappingEditorViews.h"
#include "IfcTargets.h"
#include <nlohmann/json.hpp>

class CMappingEditorDoc;

class CMappingPane : public CDialog
{
public:
   CMappingPane(UINT nIDTemplate, CMappingEditorDoc* pDoc) : CDialog(nIDTemplate), m_nIDTemplate(nIDTemplate), m_pDoc(pDoc) {}

   BOOL CreatePane(CWnd* pParent) { return Create(m_nIDTemplate, pParent); }

protected:
   UINT m_nIDTemplate;
   CMappingEditorDoc* m_pDoc;

   // Enter and Esc don't close a pane
   void OnOK() override {}
   void OnCancel() override {}
};

// The pane for a tree node
std::unique_ptr<CMappingPane> CreateMappingPane(const MappingNode& node, CMappingEditorDoc* pDoc);

// Name, extends, and comment of the table
class CTablePane : public CMappingPane
{
public:
   CTablePane(CMappingEditorDoc* pDoc) : CMappingPane(IDD_MAPPING_TABLE_PANE, pDoc) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnNameChanged();
   afx_msg void OnCommentChanged();
   afx_msg void OnExtendsSelected();
   afx_msg void OnExtendsChanged();
   afx_msg void OnBrowseExtends();
   DECLARE_MESSAGE_MAP()

private:
   void SetExtends(CString strExtends);
   void UpdateInfo();
};

// The selector of an element role
class CSelectorPane : public CMappingPane
{
public:
   CSelectorPane(CMappingEditorDoc* pDoc, const std::string& role) : CMappingPane(IDD_MAPPING_SELECTOR_PANE, pDoc), m_Role(role) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnInherited();
   afx_msg void OnDefined();
   afx_msg void OnFieldChanged();
   DECLARE_MESSAGE_MAP()

private:
   std::string m_Role;
   void Fill();
};

// The locations of a target
class CTargetPane : public CMappingPane
{
public:
   CTargetPane(CMappingEditorDoc* pDoc, const TargetDef& target) : CMappingPane(IDD_MAPPING_TARGET_PANE, pDoc), m_Target(target) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnAdd();
   afx_msg void OnEdit();
   afx_msg void OnRemove();
   afx_msg void OnUp();
   afx_msg void OnDown();
   afx_msg void OnReplace();
   afx_msg void OnLocationsDblClk(NMHDR* pNMHDR, LRESULT* pResult);
   afx_msg void OnLocationsChanged(NMHDR* pNMHDR, LRESULT* pResult);
   DECLARE_MESSAGE_MAP()

private:
   const TargetDef& m_Target;
   CListCtrl m_Locations;
   CListCtrl m_BaseLocations;

   nlohmann::ordered_json GetLocations() const;
   bool IsReplace() const;
   void SetLocations(const nlohmann::ordered_json& locations, bool bReplace);
   void FillLocations(int select);
   void UpdateButtons();
};

// A read-only view of a section (property sets, quantity sets, and classifications are edited in a later stage)
class CInfoPane : public CMappingPane
{
public:
   CInfoPane(CMappingEditorDoc* pDoc, const CString& text) : CMappingPane(IDD_MAPPING_INFO_PANE, pDoc), m_Text(text) {}

protected:
   BOOL OnInitDialog() override;
   afx_msg void OnSize(UINT nType, int cx, int cy);
   DECLARE_MESSAGE_MAP()

private:
   CString m_Text;
};

// Edits one location of a target
class CLocationDlg : public CDialog
{
public:
   CLocationDlg(const TargetDef& target, const nlohmann::ordered_json& location, CWnd* pParent = nullptr);

   nlohmann::ordered_json m_Location;

protected:
   BOOL OnInitDialog() override;
   void OnOK() override;
   afx_msg void OnKindChanged();
   afx_msg void OnParseChanged();
   DECLARE_MESSAGE_MAP()

private:
   const TargetDef& m_Target;
   int GetKind() const;
   void UpdateControls();
};

// Shows a message (e.g. the validation result) in a resizable window
class CMappingMessagesDlg : public CDialog
{
public:
   CMappingMessagesDlg(const CString& strCaption, const CString& strText, CWnd* pParent = nullptr);

protected:
   BOOL OnInitDialog() override;
   DECLARE_MESSAGE_MAP()

private:
   CString m_strCaption;
   CString m_strText;
};

// A location of the table's JSON as text, e.g. "IaDOT_PPCB.5_Final Concrete Strength, Fc" and "unit ksi"
CString DescribeLocation(const nlohmann::ordered_json& location);
CString DescribeLocationReading(const nlohmann::ordered_json& location);
