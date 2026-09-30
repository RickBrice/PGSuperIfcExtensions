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
#include "MappingEditorExportPanes.h"
#include "MappingEditorModelPane.h"
#include "MappingEditorUtil.h"
#include "IfcTargets.h"

using namespace mapping_editor;

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
   else if (lHint == CMappingEditorDoc::HINT_SELECT && pHint)
   {
      // the pane that asked is replaced, so not while its handler runs
      m_PendingSelect = ((CMappingSelectHint*)pHint)->node;
      PostMessage(WM_MAPPING_REBUILD);
   }
}

LRESULT CMappingTreeView::OnRebuild(WPARAM wParam, LPARAM lParam)
{
   Build();

   if (m_PendingSelect)
   {
      MappingNode target = *m_PendingSelect;
      m_PendingSelect.reset();

      CTreeCtrl& tree = GetTreeCtrl();
      std::vector<HTREEITEM> stack{ tree.GetRootItem() };
      while (!stack.empty())
      {
         HTREEITEM hItem = stack.back();
         stack.pop_back();
         for (; hItem; hItem = tree.GetNextSiblingItem(hItem))
         {
            if (m_Nodes[tree.GetItemData(hItem)].Matches(target))
            {
               tree.EnsureVisible(hItem);
               if (tree.GetSelectedItem() == hItem)
               {
                  // already selected (e.g. the whole table was replaced): show its pane again
                  auto pFrame = (CMappingEditorFrame*)GetParentFrame();
                  pFrame->GetDetailView()->ShowNode(m_Nodes[tree.GetItemData(hItem)]);
                  return 0;
               }
               tree.SelectItem(hItem); // shows its pane
               if (tree.GetSelectedItem() == hItem)
                  GetParentFrame()->SetFocus();
               return 0;
            }
            if (HTREEITEM hChild = tree.GetChildItem(hItem))
               stack.push_back(hChild);
         }
      }
   }
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

   AddExportSections();

   // the IFC model, for picking properties and trying the table (stage 4)
   CString strModel(_T("Model"));
   if (const auto* pModel = GetDocument()->GetModel())
      strModel += _T(": ") + CString(pModel->GetPath().filename().c_str());
   Add(TVI_ROOT, strModel, { MappingNode::Kind::Model, "" }, GetDocument()->GetModel() != nullptr);

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
   else if (!selected)
   {
      tree.SelectItem(tree.GetRootItem()); // a new tree: show the table pane
   }
   else if (!m_PendingSelect)
   {
      // the selected item is gone (e.g. its entry no longer applies to that role). Its pane may be the one that changed
      // the table, so it's replaced after its handler returns
      m_PendingSelect = MappingNode{ MappingNode::Kind::Table, "" };
      PostMessage(WM_MAPPING_REBUILD);
   }
}

void CMappingTreeView::AddExportSections()
{
   const auto& table = GetDocument()->GetTable();
   CMappingEditorDoc* pDoc = GetDocument();
   auto section = [&table](const char* key) { return (table.contains(key) && table[key].is_array()) ? table[key] : ordered_json::array(); };
   auto kind_of = [](const char* name) { ElementKind kind = ElementKind::Bridge; GetElementKind(name, kind); return kind; };

   // property sets and quantity sets, by element role
   for (bool bQuantities : { false, true })
   {
      const char* key = bQuantities ? "quantity_sets" : "property_sets";
      auto sets = section(key);
      auto header_kind = bQuantities ? MappingNode::Kind::QuantitySets : MappingNode::Kind::PropertySets;
      auto group_kind = bQuantities ? MappingNode::Kind::QuantitySetGroup : MappingNode::Kind::PropertySetGroup;
      auto item_kind = bQuantities ? MappingNode::Kind::QuantitySet : MappingNode::Kind::PropertySet;

      CString header = bQuantities ? _T("Quantity sets") : _T("Property sets");
      if (!sets.empty())
         header.AppendFormat(_T(" (%d in this table)"), (int)sets.size());
      HTREEITEM hHeader = Add(TVI_ROOT, header, { header_kind, "" }, !sets.empty());

      // entries without a valid element role, so they can still be selected and corrected
      for (int i = 0; i < (int)sets.size(); i++)
      {
         bool bRole = false;
         for (const auto& role : roles_of(sets[i]))
         {
            ElementKind kind;
            bRole |= GetElementKind(role, kind);
         }
         if (!bRole)
            Add(hHeader, Utf8ToCString(string_value(sets[i], "name") + " (no element role)"), { item_kind, "", i }, true);
      }

      for (const auto& role : all_roles())
      {
         std::vector<int> own;
         for (int i = 0; i < (int)sets.size(); i++)
         {
            auto roles = roles_of(sets[i]);
            if (std::find(roles.begin(), roles.end(), role) != roles.end())
               own.push_back(i);
         }
         auto base = GetBaseSets(pDoc, bQuantities, kind_of(role.c_str()));
         if (own.empty() && base.empty())
            continue;

         HTREEITEM hGroup = Add(hHeader, Utf8ToCString(role), { group_kind, role }, !own.empty());
         for (int i : own)
         {
            const auto& set = sets[i];
            std::string label = string_value(set, "name");
            auto attach = string_value(set, "attach");
            if (!attach.empty() && attach != "occurrence")
               label += " [" + attach + "]";
            if (set.value("remove", false))
               label += " (left out)";
            Add(hGroup, Utf8ToCString(label), { item_kind, role, i }, true);
         }
         for (const auto* pset : base)
         {
            std::string label = pset->name;
            if (pset->attach != PropertyOwner::Occurrence)
               label += std::string(" [") + OwnerName(pset->attach) + "]";
            Add(hGroup, Utf8ToCString(label), { item_kind, BaseSetKey(*pset) }, false);
         }
      }
   }

   // classification systems
   {
      auto systems = section("classification_systems");
      CString header(_T("Classification systems"));
      if (!systems.empty())
         header.AppendFormat(_T(" (%d in this table)"), (int)systems.size());
      HTREEITEM hHeader = Add(TVI_ROOT, header, { MappingNode::Kind::ClassificationSystems, "" }, !systems.empty());
      for (int i = 0; i < (int)systems.size(); i++)
         Add(hHeader, Utf8ToCString(string_value(systems[i], "name")), { MappingNode::Kind::ClassificationSystem, "", i }, true);
      for (const auto* system : GetBaseSystems(pDoc))
         Add(hHeader, Utf8ToCString(system->name), { MappingNode::Kind::ClassificationSystem, system->name }, false);
   }

   // classifications, by element role
   {
      auto classifications = section("classifications");
      CString header(_T("Classifications"));
      if (!classifications.empty())
         header.AppendFormat(_T(" (%d in this table)"), (int)classifications.size());
      HTREEITEM hHeader = Add(TVI_ROOT, header, { MappingNode::Kind::Classifications, "" }, !classifications.empty());
      for (int i = 0; i < (int)classifications.size(); i++)
      {
         bool bRole = false;
         for (const auto& role : roles_of(classifications[i]))
         {
            ElementKind kind;
            bRole |= GetElementKind(role, kind);
         }
         if (!bRole)
            Add(hHeader, Utf8ToCString(string_value(classifications[i], "identification") + " (no element role)"), { MappingNode::Kind::Classification, "", i }, true);
      }
      for (const auto& role : all_roles())
      {
         std::vector<int> own;
         for (int i = 0; i < (int)classifications.size(); i++)
         {
            auto roles = roles_of(classifications[i]);
            if (std::find(roles.begin(), roles.end(), role) != roles.end())
               own.push_back(i);
         }
         auto base = GetBaseClassifications(pDoc, kind_of(role.c_str()));
         if (own.empty() && base.empty())
            continue;

         HTREEITEM hGroup = Add(hHeader, Utf8ToCString(role), { MappingNode::Kind::ClassificationGroup, role }, !own.empty());
         for (int i : own)
         {
            std::string label = string_value(classifications[i], "identification");
            if (classifications[i].value("remove", false))
               label += " (left out)";
            Add(hGroup, Utf8ToCString(label), { MappingNode::Kind::Classification, role, i }, true);
         }
         for (const auto* c : base)
            Add(hGroup, Utf8ToCString(c->identification), { MappingNode::Kind::Classification, BaseClassificationKey(*c) }, false);
      }
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

IMPLEMENT_DYNCREATE(CMappingDetailView, CScrollView)

BEGIN_MESSAGE_MAP(CMappingDetailView, CScrollView)
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

void CMappingDetailView::OnInitialUpdate()
{
   SetScrollSizes(MM_TEXT, CSize(1, 1));
   CScrollView::OnInitialUpdate();
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

   ScrollToPosition(CPoint(0, 0));
   m_pPane = CreateMappingPane(node, GetDocument());
   if (m_pPane && m_pPane->CreatePane(this))
   {
      CRect rect;
      m_pPane->GetWindowRect(&rect); // the dialog template's size
      m_PaneSize = rect.Size();
      SetScrollSizes(MM_TEXT, m_PaneSize);
      SizePane();
      m_pPane->ShowWindow(SW_SHOW);
   }
   else
   {
      m_pPane.reset();
      m_PaneSize = CSize(0, 0);
      SetScrollSizes(MM_TEXT, CSize(1, 1));
   }
}

void CMappingDetailView::OnSize(UINT nType, int cx, int cy)
{
   CScrollView::OnSize(nType, cx, cy);
   SizePane();
}

void CMappingDetailView::SizePane()
{
   if (m_pPane && m_pPane->GetSafeHwnd())
   {
      // fills the view, and is at least its designed size (the view scrolls)
      CRect rect;
      GetClientRect(&rect);
      CPoint scroll = GetScrollPosition();
      m_pPane->SetWindowPos(nullptr, -scroll.x, -scroll.y, std::max(rect.Width(), (int)m_PaneSize.cx), std::max(rect.Height(), (int)m_PaneSize.cy), SWP_NOZORDER | SWP_NOACTIVATE);
   }
}

BOOL CMappingDetailView::PreTranslateMessage(MSG* pMsg)
{
   // tab between the pane's controls
   if (m_pPane && m_pPane->GetSafeHwnd() && m_pPane->IsDialogMessage(pMsg))
      return TRUE;
   return CScrollView::PreTranslateMessage(pMsg);
}
