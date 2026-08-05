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
#include <EAF\Transaction.h>

BEGIN_MESSAGE_MAP(CIfcExtensionAgent,CCmdTarget)
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
   return true;
}

bool CIfcExtensionAgent::Reset()
{
   EAF_AGENT_RESET;
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

   //std::_tstring epsgCode;
   //if ( !pStrLoad->Property(_T("EPSGCode"),&epsgCode) )
   //   return WBFL::EAF::Broker::LoadResult::Error;

   //m_EPSGCode = epsgCode.c_str();

   if ( !pStrLoad->EndUnit() )
      return WBFL::EAF::Broker::LoadResult::Error;

   return WBFL::EAF::Broker::LoadResult::Success;
}

bool CIfcExtensionAgent::Save(WBFL::System::IStructuredSave* pStrSave)
{
   pStrSave->BeginUnit(_T("IfcExtensionAgent"),1.0);
   //pStrSave->Property(_T("EPSGCode"),m_EPSGCode);
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
   }
   else
   {
      UnregisterUIExtensions();
   }

   return true;
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
   m_GeoreferencingData = pMyPage->m_GeoRefData;
#pragma Reminder("TODO: Implement a transaction to save the georeferencing data to the alignment data")
   return nullptr;
}
