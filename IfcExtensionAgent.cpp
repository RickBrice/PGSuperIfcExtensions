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

// IfcExtensionAgent.cpp : Implementation of CIfcExtensionAgent

#include "stdafx.h"
#include "IfcExtensions.h"
#include "IfcExtensionAgent.h"

#include <IFace\Tools.h>
#include <IFace\EditByUI.h>
#include <EAF\Transaction.h>
#include "EditGeoreferencing.h"

BEGIN_MESSAGE_MAP(CIfcExtensionAgent,CCmdTarget)
   ON_COMMAND(ID_EDIT_GEOREFERENCING,&CIfcExtensionAgent::OnEditGeoreferencing)
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////
// IAgentEx

bool CIfcExtensionAgent::RegisterInterfaces()
{
   EAF_AGENT_REGISTER_INTERFACES;
   REGISTER_INTERFACE(IGeoreferencing);

   return true;
}

bool CIfcExtensionAgent::Init()
{
   EAF_AGENT_INIT;
   CREATE_LOGFILE(_T("IfcExtensionAgent"));

   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_bmpMenu.LoadBitmap(IDB_BSI));

   return true;
}

bool CIfcExtensionAgent::Reset()
{
   EAF_AGENT_RESET;
   m_bmpMenu.DeleteObject();
   return true;
}

bool CIfcExtensionAgent::ShutDown()
{
   EAF_AGENT_SHUTDOWN;
   return true;
}

CLSID CIfcExtensionAgent::GetCLSID() const
{
   return CLSID_PGSuperIfcExtensionAgent;
}

////////////////////////////////////////////////////////////////////
// IAgentPersist

WBFL::EAF::Broker::LoadResult CIfcExtensionAgent::Load(WBFL::System::IStructuredLoad* pStrLoad)
{
   if ( !pStrLoad->BeginUnit(_T("IfcExtensionAgent")) )
      return WBFL::EAF::Broker::LoadResult::Error;

   m_GeoreferencingData.Load(pStrLoad);

   if ( !pStrLoad->EndUnit() )
      return WBFL::EAF::Broker::LoadResult::Error;

   return WBFL::EAF::Broker::LoadResult::Success;
}

bool CIfcExtensionAgent::Save(WBFL::System::IStructuredSave* pStrSave)
{
   pStrSave->BeginUnit(_T("IfcExtensionAgent"),1.0);
   m_GeoreferencingData.Save(pStrSave);
   pStrSave->EndUnit();
   return true;
}

////////////////////////////////////////////////////////////////////
// IAgentUIIntegration

bool CIfcExtensionAgent::IntegrateWithUI(bool bIntegrate)
{
   if ( bIntegrate )
   {
      RegisterUIExtensions();
      CreateMenus();
   }
   else
   {
      RemoveMenus();
      UnregisterUIExtensions();
   }

   return true;
}

void CIfcExtensionAgent::CreateMenus()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   GET_IFACE(IEAFMainMenu,pMainMenu);
   auto pMenu = pMainMenu->GetMainMenu();

   UINT editPos = pMenu->FindMenuItem(_T("Edit"));
   m_pEditMenu = pMenu->GetSubMenu(editPos);

   UINT alignmentPos = m_pEditMenu->FindMenuItem(_T("Alignment..."));

   auto callback = std::dynamic_pointer_cast<WBFL::EAF::ICommandCallback>(shared_from_this());
   m_pEditMenu->InsertMenu(alignmentPos, ID_EDIT_GEOREFERENCING, _T("&Georeferencing..."), callback);
   m_pEditMenu->SetMenuItemBitmaps(ID_EDIT_GEOREFERENCING, MF_BYCOMMAND, &m_bmpMenu, nullptr, callback);
}

void CIfcExtensionAgent::RemoveMenus()
{
   if ( m_pEditMenu )
   {
      auto callback = std::dynamic_pointer_cast<WBFL::EAF::ICommandCallback>(shared_from_this());
      m_pEditMenu->RemoveMenu(ID_EDIT_GEOREFERENCING, MF_BYCOMMAND, callback);
      m_pEditMenu.reset();
   }
}

////////////////////////////////////////////////////////////////////
// ICommandCallback

BOOL CIfcExtensionAgent::OnCommandMessage(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo)
{
   return OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
}

BOOL CIfcExtensionAgent::GetStatusBarMessageString(UINT nID, CString& rMessage) const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   if (rMessage.LoadString(nID))
   {
      // first newline terminates actual string
      rMessage.Replace('\n', '\0');
   }
   else
   {
      TRACE1("Warning (CIfcExtensionAgent): no message line prompt for ID 0x%04X.\n", nID);
   }

   return TRUE;
}

BOOL CIfcExtensionAgent::GetToolTipMessageString(UINT nID, CString& rMessage) const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CString string;
   if (string.LoadString(nID))
   {
      // tip is after first newline
      int pos = string.Find('\n');
      if (0 < pos)
         rMessage = string.Mid(pos + 1);
   }
   else
   {
      TRACE1("Warning (CIfcExtensionAgent): no tool tip for ID 0x%04X.\n", nID);
   }

   return TRUE;
}

void CIfcExtensionAgent::OnEditGeoreferencing()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   GET_IFACE(IEditByUI,pEditByUI);
   pEditByUI->EditAlignmentDescription(_T("Georeferencing"));
}

void CIfcExtensionAgent::RegisterUIExtensions()
{
   // RegisterEditAlignmentCallback is declared on the shared IExtendUI base, but only
   // IExtendPGSuperUI/IExtendPGSpliceUI are ever registered on the broker (see
   // CPGSuperDocProxyAgent::RegisterInterfaces) - IID_IExtendUI itself is never
   // queryable directly. This plugin is PGSuper-specific, so IExtendPGSuperUI is the
   // right one to go through (same pattern ExampleExtensionAgent uses).
   GET_IFACE(IExtendPGSuperUI,pExtendPGSuperUI);
   m_EditAlignmentCallbackID = pExtendPGSuperUI->RegisterEditAlignmentCallback(this);
}

void CIfcExtensionAgent::UnregisterUIExtensions()
{
   GET_IFACE(IExtendPGSuperUI,pExtendPGSuperUI);
   pExtendPGSuperUI->UnregisterEditAlignmentCallback(m_EditAlignmentCallbackID);
}

////////////////////////////////////////////////////////////////////
// IGeoreferencing
void CIfcExtensionAgent::SetGeoreferencingData(const GeoreferencingData& data)
{
   m_GeoreferencingData = data;
}

const GeoreferencingData& CIfcExtensionAgent::GetGeoreferencingData() const
{
   return m_GeoreferencingData;
}


////////////////////////////////////////////////////////////////////
// IEditAlignmentCallback

CPropertyPage* CIfcExtensionAgent::CreatePropertyPage(IEditAlignmentData* pAlignmentData)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CGeoReferencingPage* pPage = new CGeoReferencingPage();
   pPage->m_GeoRefData = m_GeoreferencingData;
   return pPage;
}

std::unique_ptr<WBFL::EAF::Transaction> CIfcExtensionAgent::OnOK(CPropertyPage* pPage,IEditAlignmentData* pAlignmentData)
{
   CGeoReferencingPage* pMyPage = (CGeoReferencingPage*)pPage;
   auto pTxn = std::make_unique<txnEditGeoreferencing>(m_GeoreferencingData, pMyPage->m_GeoRefData);
   m_GeoreferencingData = pMyPage->m_GeoRefData;
   return pTxn;
}
