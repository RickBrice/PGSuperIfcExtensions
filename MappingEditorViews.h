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

// Mapping table editor frame: the table's sections in a tree (left) and the selected item's pane (right)

#include <afxcview.h>
#include <memory>
#include <optional>
#include <string>
#include <vector>

class CMappingEditorDoc;
class CMappingPane;

// An item of the tree
struct MappingNode
{
   enum class Kind
   {
      Table, Roles, Role, Targets, TargetGroup, Target,
      PropertySets, PropertySetGroup, PropertySet,
      QuantitySets, QuantitySetGroup, QuantitySet,
      ClassificationSystems, ClassificationSystem,
      Classifications, ClassificationGroup, Classification,
      Model
   };
   Kind kind = Kind::Table;
   // Role: role name; TargetGroup, PropertySetGroup, QuantitySetGroup, ClassificationGroup: element role name; Target: target name.
   // This table's entry (index >= 0): the role of the group it's listed in (an entry can apply to several roles).
   // A base table's entry (index < 0): its identity (BaseSetKey, BaseClassificationKey, or the system name)
   std::string key;
   int index = -1; // this table's entry: its index in the section's JSON list

   bool operator==(const MappingNode& other) const { return kind == other.kind && key == other.key && index == other.index; }

   // A node to select: the key is ignored when it's empty (e.g. a new entry, selected wherever it's listed first)
   bool Matches(const MappingNode& target) const { return kind == target.kind && index == target.index && (target.key.empty() || key == target.key); }
};

// UpdateAllViews hint object for CMappingEditorDoc::HINT_SELECT
class CMappingSelectHint : public CObject
{
public:
   MappingNode node;
};

class CMappingTreeView : public CTreeView
{
protected:
   CMappingTreeView() = default;
   DECLARE_DYNCREATE(CMappingTreeView)

public:
   CMappingEditorDoc* GetDocument() const;

protected:
   void OnInitialUpdate() override;
   void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) override;
   BOOL PreCreateWindow(CREATESTRUCT& cs) override;

   afx_msg void OnSelChanged(NMHDR* pNMHDR, LRESULT* pResult);
   afx_msg LRESULT OnRebuild(WPARAM wParam, LPARAM lParam);
   DECLARE_MESSAGE_MAP()

private:
   std::vector<MappingNode> m_Nodes; // item data is the index
   bool m_bBuilding = false;
   bool m_bInSelChange = false; // a pane may commit an edit (and change the table) while the selection changes
   std::optional<MappingNode> m_PendingSelect; // selected after the next (posted) rebuild

   // Adds the export sections (property sets, quantity sets, classification systems, classifications)
   void AddExportSections();

   // (Re)builds the tree and selects the node that was selected
   void Build();
   HTREEITEM Add(HTREEITEM hParent, const CString& label, MappingNode node, bool bBold);
};

// Hosts the pane of the selected node. It scrolls when the pane is bigger than the view
class CMappingDetailView : public CScrollView
{
protected:
   CMappingDetailView() = default;
   DECLARE_DYNCREATE(CMappingDetailView)

public:
   CMappingEditorDoc* GetDocument() const;

   // Shows the pane for a tree node
   void ShowNode(MappingNode node);

protected:
   void OnDraw(CDC* pDC) override;
   void OnInitialUpdate() override;
   void OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint) override;
   BOOL PreTranslateMessage(MSG* pMsg) override;

   afx_msg void OnSize(UINT nType, int cx, int cy);
   afx_msg BOOL OnEraseBkgnd(CDC* pDC);
   DECLARE_MESSAGE_MAP()

private:
   std::unique_ptr<CMappingPane> m_pPane;
   CSize m_PaneSize{ 0, 0 }; // the pane's size as designed (its dialog template)
   void SizePane();
};

class CMappingEditorFrame : public CMDIChildWnd
{
   DECLARE_DYNCREATE(CMappingEditorFrame)

public:
   CMappingDetailView* GetDetailView();
   void ActivateFrame(int nCmdShow = -1) override;

protected:
   BOOL OnCreateClient(LPCREATESTRUCT lpcs, CCreateContext* pContext) override;

   DECLARE_MESSAGE_MAP()

private:
   CSplitterWnd m_Splitter;
};
