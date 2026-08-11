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

// IfcExtensionAgent.h : Declaration of the CIfcExtensionAgent
//
// A general PGSuperIfcExtensions agent - not a single-purpose "GeoReferencing
// agent". Capturing the alignment's EPSG code is its first responsibility, but this
// is the place any future IFC-related project settings for this plugin should live,
// under the same persisted unit (see Save()/Load()), rather than adding another
// agent/CLSID per setting.
//
// Modeled on PGSuper\ExtensionAgentExample\ExampleExtensionAgent.h/.cpp, trimmed to
// just what this plugin needs: persistence (IAgentPersist) and one UI extension
// (IEditAlignmentCallback). No menus, toolbar, views, reports, or graphs.

#pragma once

#include <EAF\Agent.h>
#include <EAF\EAFUIIntegration.h>
#include <IFace\ExtendUI.h>
#include "GeoReferencingPage.h"
#include "GeoReferencing.h"

class CIfcExtensionAgent : public CCmdTarget, // must be first parent for inheritance, see Warning C4407
   public WBFL::EAF::Agent,
   public WBFL::EAF::IAgentPersist,
   public WBFL::EAF::IAgentUIIntegration,
   public IGeoreferencing,
   public IEditAlignmentCallback
{
public:
   CIfcExtensionAgent()
   {
   }

// IAgentEx
public:
   std::_tstring GetName() const override { return _T("IfcExtensionAgent"); }
   bool RegisterInterfaces() override;
   bool Init() override;
   bool Reset() override;
   bool ShutDown() override;
   CLSID GetCLSID() const override;

// IAgentPersist
public:
   WBFL::EAF::Broker::LoadResult Load(WBFL::System::IStructuredLoad* pStrLoad) override;
   bool Save(WBFL::System::IStructuredSave* pStrSave) override;

// IAgentUIIntegration
public:
   bool IntegrateWithUI(bool bIntegrate) override;

// IGeoreferencing
public:
   void SetGeoreferencingData(const GeoreferencingData& data) override;
   const GeoreferencingData& GetGeoreferencingData() const override;

// IEditAlignmentCallback
public:
   CPropertyPage* CreatePropertyPage(IEditAlignmentData* pAlignmentData) override;
   std::unique_ptr<WBFL::EAF::Transaction> OnOK(CPropertyPage* pPage,IEditAlignmentData* pAlignmentData) override;
   // This is what actually makes the GeoReferencing page the first tab - the dialog
   // itself has no "first" special-casing, it just asks every registered
   // IEditAlignmentCallback where its page belongs (see CExtensionPageManager).
   ExtensionPagePosition GetPropertyPagePosition() override { return ExtensionPagePosition::AtStart(); }

   // Names this page "Georeferencing" instead of the framework's auto-generated "ExtensionN",
   // so a future menu command could look it up by a fixed name instead of computing it.
   std::_tstring GetPropertyPageName() override { return _T("Georeferencing"); }

   DECLARE_MESSAGE_MAP()

private:
   EAF_DECLARE_AGENT_DATA;

   void RegisterUIExtensions();
   void UnregisterUIExtensions();
   IDType m_EditAlignmentCallbackID;

   GeoreferencingData m_GeoreferencingData;
};
