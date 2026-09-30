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

// MappingEditorViews.cpp : the mapping table editor's frame and views
//

#include "stdafx.h"
#include "resource.h"
#include "MappingEditorViews.h"
#include "MappingEditorDoc.h"
#include "MappingEditorPanes.h"
#include "IfcTargets.h"

/////////////////////////////////////////////////////////////////////////////
// CMappingEditorFrame

IMPLEMENT_DYNCREATE(CMappingEditorFrame, CMDIChildWnd)

BEGIN_MESSAGE_MAP(CMappingEditorFrame, CMDIChildWnd)
END_MESSAGE_MAP()

BOOL CMappingEditorFrame::OnCreateClient(LPCREATESTRUCT lpcs, CCreateContext* pContext)
{
   if (!m_Splitter.CreateStatic(this, 1, 2))
      return FALSE;

   if (!m_Splitter.CreateView(0, 0, RUNTIME_CLASS(CMappingTreeView), CSize(280, 0), pContext) ||
      !m_Splitter.CreateView(0, 1, RUNTIME_CLASS(CMappingDetailView), CSize(0, 0), pContext))
      return FALSE;

   SetActiveView((CView*)m_Splitter.GetPane(0, 0));
   return TRUE;
}

CMappingDetailView* CMappingEditorFrame::GetDetailView()
{
   return (CMappingDetailView*)m_Splitter.GetPane(0, 1);
}

/////////////////////////////////////////////////////////////////////////////
// CMappingTreeView

IMPLEMENT_DYNCREATE(CMappingTreeView, CTreeView)

static const UINT WM_MAPPING_REBUILD = WM_APP + 1;

BEGIN_MESSAGE_MAP(CMappingTreeView, CTreeView)
   ON_NOTIFY_REFLECT(TVN_SELCHANGED, &CMappingTreeView::OnSelChanged)
   ON_MESSAGE(WM_MAPPING_REBUILD, &CMappingTreeView::OnRebuild)
END_MESSAGE_MAP()

CMappingEditorDoc* CMappingTreeView::GetDocument() const
{
   return (CMappingEditorDoc*)m_pDocument;
}

BOOL CMappingTreeView::PreCreateWindow(CREATESTRUCT& cs)
{
   cs.style |= TVS_HASLINES | TVS_LINESATROOT | TVS_HASBUTTONS | TVS_SHOWSELALWAYS;
   return CTreeView::PreCreateWindow(cs);
}

void CMappingTreeView::OnInitialUpdate()
{
   CTreeView::OnInitialUpdate();
   Build();
}

void CMappingTreeView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint)
{
   if (lHint == CMappingEditorDoc::HINT_STRUCTURE)
   {
      if (m_bInSelChange)
         PostMessage(WM_MAPPING_REBUILD); // don't delete the items while the tree is changing the selection
      else
         Build();
   }
}

LRESULT CMappingTreeView::OnRebuild(WPARAM wParam, LPARAM lParam)
{
   Build();
   return 0;
}

HTREEITEM CMappingTreeView::Add(HTREEITEM hParent, const CString& label, MappingNode node, bool bBold)
{
   CTreeCtrl& tree = GetTreeCtrl();
   HTREEITEM hItem = tree.InsertItem(label, hParent);
   tree.SetItemData(hItem, (DWORD_PTR)m_Nodes.size());
   if (bBold)
      tree.SetItemState(hItem, TVIS_BOLD, TVIS_BOLD);
   m_Nodes.push_back(std::move(node));
   return hItem;
}

void CMappingTreeView::Build()
{
   CTreeCtrl& tree = GetTreeCtrl();
   const auto& table = GetDocument()->GetTable();

   // keep the selection and the expanded sections
   std::optional<MappingNode> selected;
   std::vector<MappingNode> expanded;
   if (HTREEITEM hSelected = tree.GetSelectedItem())
      selected = m_Nodes[tree.GetItemData(hSelected)];
   for (HTREEITEM hItem = tree.GetRootItem(); hItem; hItem = tree.GetNextItem(hItem, TVGN_NEXTVISIBLE))
   {
      if (tree.GetItemState(hItem, TVIS_EXPANDED) & TVIS_EXPANDED)
         expanded.push_back(m_Nodes[tree.GetItemData(hItem)]);
   }

   m_bBuilding = true;
   tree.SetRedraw(FALSE);
   tree.DeleteAllItems();
   m_Nodes.clear();

   const auto elements = table.contains("elements") ? table["elements"] : nlohmann::ordered_json::object();
   const auto targets = table.contains("targets") ? table["targets"] : nlohmann::ordered_json::object();

   Add(TVI_ROOT, _T("Table"), { MappingNode::Kind::Table, "" }, false);

   // element roles, bold if this table defines the selector
   HTREEITEM hRoles = Add(TVI_ROOT, _T("Element roles"), { MappingNode::Kind::Roles, "" }, elements.is_object() && !elements.empty());
   for (int kind = (int)ElementKind::Project; kind <= (int)ElementKind::Barrier; kind++)
   {
      std::string role(GetElementRoleName((ElementKind)kind));
      Add(hRoles, Utf8ToCString(role), { MappingNode::Kind::Role, role }, elements.is_object() && elements.contains(role));
   }

   // targets by element, bold if this table lists locations for them
   HTREEITEM hTargets = Add(TVI_ROOT, _T("Targets"), { MappingNode::Kind::Targets, "" }, targets.is_object() && !targets.empty());
   std::vector<std::pair<std::string, HTREEITEM>> groups;
   for (const auto& target : GetTargetDefs())
   {
      std::string role(GetElementRoleName(target.element));
      auto group = std::find_if(groups.begin(), groups.end(), [&role](const auto& g) {return g.first == role; });
      if (group == groups.end())
      {
         bool bGroupBold = false;
         for (const auto& t : GetTargetDefs())
            bGroupBold |= (t.element == target.element && targets.is_object() && targets.contains(std::string(t.name)));
         groups.emplace_back(role, Add(hTargets, Utf8ToCString(role), { MappingNode::Kind::TargetGroup, role }, bGroupBold));
         group = groups.end() - 1;
      }
      std::string name(target.name);
      Add(group->second, Utf8ToCString(name), { MappingNode::Kind::Target, name }, targets.is_object() && targets.contains(name));
   }

   auto count_label = [&table](const char* key, LPCTSTR label)
   {
      CString text(label);
      if (table.contains(key) && table[key].is_array())
         text.AppendFormat(_T(" (%d)"), (int)table[key].size());
      return text;
   };
   auto has = [&table](const char* key) {return table.contains(key) && table[key].is_array() && !table[key].empty(); };
   Add(TVI_ROOT, count_label("property_sets", _T("Property sets")), { MappingNode::Kind::PropertySets, "" }, has("property_sets"));
   Add(TVI_ROOT, count_label("quantity_sets", _T("Quantity sets")), { MappingNode::Kind::QuantitySets, "" }, has("quantity_sets"));
   Add(TVI_ROOT, count_label("classification_systems", _T("Classification systems")), { MappingNode::Kind::ClassificationSystems, "" }, has("classification_systems"));
   Add(TVI_ROOT, count_label("classifications", _T("Classifications")), { MappingNode::Kind::Classifications, "" }, has("classifications"));

   // restore
   HTREEITEM hSelect = nullptr;
   std::vector<HTREEITEM> items;
   for (HTREEITEM hItem = tree.GetRootItem(); hItem; )
   {
      items.push_back(hItem);
      HTREEITEM hNext = tree.GetChildItem(hItem);
      if (!hNext)
      {
         while (hItem && !(hNext = tree.GetNextSiblingItem(hItem)))
            hItem = tree.GetParentItem(hItem);
      }
      hItem = hNext;
   }
   for (HTREEITEM hItem : items)
   {
      const auto& node = m_Nodes[tree.GetItemData(hItem)];
      if (std::find(expanded.begin(), expanded.end(), node) != expanded.end())
         tree.Expand(hItem, TVE_EXPAND);
      if (selected && node == *selected)
         hSelect = hItem;
   }

   tree.SetRedraw(TRUE);
   m_bBuilding = false;

   if (hSelect)
   {
      // the pane shown for it stays (it made the change)
      m_bBuilding = true;
      tree.SelectItem(hSelect);
      m_bBuilding = false;
   }
   else
   {
      tree.SelectItem(tree.GetRootItem()); // shows the table pane
   }
}

void CMappingTreeView::OnSelChanged(NMHDR* pNMHDR, LRESULT* pResult)
{
   *pResult = 0;
   if (m_bBuilding)
      return;

   NMTREEVIEW* pNMTreeView = reinterpret_cast<NMTREEVIEW*>(pNMHDR);
   if (!pNMTreeView->itemNew.hItem)
      return;

   MappingNode node = m_Nodes[GetTreeCtrl().GetItemData(pNMTreeView->itemNew.hItem)];
   auto pFrame = (CMappingEditorFrame*)GetParentFrame();
   m_bInSelChange = true;
   pFrame->GetDetailView()->ShowNode(node);
   m_bInSelChange = false;
}

/////////////////////////////////////////////////////////////////////////////
// CMappingDetailView

IMPLEMENT_DYNCREATE(CMappingDetailView, CView)

BEGIN_MESSAGE_MAP(CMappingDetailView, CView)
   ON_WM_SIZE()
   ON_WM_ERASEBKGND()
END_MESSAGE_MAP()

CMappingEditorDoc* CMappingDetailView::GetDocument() const
{
   return (CMappingEditorDoc*)m_pDocument;
}

void CMappingDetailView::OnDraw(CDC* pDC)
{
}

BOOL CMappingDetailView::OnEraseBkgnd(CDC* pDC)
{
   CRect rect;
   GetClientRect(&rect);
   pDC->FillSolidRect(&rect, ::GetSysColor(COLOR_BTNFACE));
   return TRUE;
}

void CMappingDetailView::OnUpdate(CView* pSender, LPARAM lHint, CObject* pHint)
{
   // the pane that changed the table shows the change itself
}

void CMappingDetailView::ShowNode(MappingNode node)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   if (m_pPane)
   {
      // a field being edited is written to the table when it loses the focus
      SetFocus();
      m_pPane->DestroyWindow();
      m_pPane.reset();
   }

   m_pPane = CreateMappingPane(node, GetDocument());
   if (m_pPane && m_pPane->CreatePane(this))
   {
      SizePane();
      m_pPane->ShowWindow(SW_SHOW);
   }
   else
   {
      m_pPane.reset();
   }
}

void CMappingDetailView::OnSize(UINT nType, int cx, int cy)
{
   CView::OnSize(nType, cx, cy);
   SizePane();
}

void CMappingDetailView::SizePane()
{
   if (m_pPane && m_pPane->GetSafeHwnd())
   {
      CRect rect;
      GetClientRect(&rect);
      m_pPane->SetWindowPos(nullptr, 0, 0, rect.Width(), rect.Height(), SWP_NOZORDER | SWP_NOACTIVATE);
   }
}

BOOL CMappingDetailView::PreTranslateMessage(MSG* pMsg)
{
   // tab between the pane's controls
   if (m_pPane && m_pPane->GetSafeHwnd() && m_pPane->IsDialogMessage(pMsg))
      return TRUE;
   return CView::PreTranslateMessage(pMsg);
}
