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


// Note: Proxy/Stub Information
//      To build a separate proxy/stub DLL, 
//      run nmake -f IEPluginExampleps.mk in the project directory.

#include "stdafx.h"
#include "resource.h"
#include <initguid.h>
#include <Plugins\PGSuperIEPlugin.h>
#include "IfcExtensions.h"
#include "PGSuperCatCom.h"
#include <WBFLCogo_i.c>
#include <WBFLGeometry_i.c>

// IID_IEAFStatusCenter (used by EAF_AGENT_INIT in IfcExtensionAgent.cpp, the first agent
// this project has ever had) is declared via DEFINE_GUID in EAF\EAFStatusCenter.h, which
// only actually emits storage under INITGUID. The usual fix - #include <initguid.h> before
// the header, in some .cpp of the project - doesn't work here: stdafx.h's <PGSuperAll.h>
// already pulls EAFStatusCenter.h in ahead of the <initguid.h> above, and its #pragma once
// locks that non-storage form in for every translation unit in this project. So give the
// symbol real storage directly instead, matching EAFStatusCenter.h's DEFINE_GUID exactly:
// {77977E9B-B074-401f-8994-73A418FC4FFF}
const GUID IID_IEAFStatusCenter = { 0x77977e9b, 0xb074, 0x401f, { 0x89, 0x94, 0x73, 0xa4, 0x18, 0xfc, 0x4f, 0xff } };

#include "PGSuperDataImporter.h"
#include "PGSuperDataExporter.h"
#include "PGSuperProjectImporter.h"
#include "IfcExtensionAgent.h"

#include <IFace\Project.h>
#include <IFace\VersionInfo.h>
#include <IFace\Alignment.h>
#include <IFace\DocumentType.h>
#include <IFace/PointOfInterest.h>
#include <IFace\Bridge.h>
#include <IFace\Intervals.h>
#include <IFace\AnalysisResults.h>

#include <EAF\EAFDisplayUnits.h>
#include <EAF/EAFUIIntegration.h>
#include <EAF/EAFProgress.h>
#include <Plugins\BeamFamilyCLSID.h>


#include <EAF\ComponentModule.h>
WBFL::EAF::ComponentModule _Module;
EAF_BEGIN_OBJECT_MAP(ObjectMap)
   EAF_OBJECT_ENTRY(CLSID_PGSuperIfcImporter, CPGSuperDataImporter)
   EAF_OBJECT_ENTRY(CLSID_PGSuperIfcExporter, CPGSuperDataExporter)
   EAF_OBJECT_ENTRY(CLSID_PGSuperIfcProjectImporter, CPGSuperProjectImporter)
   EAF_OBJECT_ENTRY(CLSID_PGSuperIfcExtensionAgent, CIfcExtensionAgent)
EAF_END_OBJECT_MAP()

class CIFCExtensionsApp : public CWinApp
{
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CIFCExtensionsApp)
	public:
    virtual BOOL InitInstance() override;
    virtual int ExitInstance() override;
	//}}AFX_VIRTUAL

	//{{AFX_MSG(CIFCExtensionsApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

BEGIN_MESSAGE_MAP(CIFCExtensionsApp, CWinApp)
	//{{AFX_MSG_MAP(CIFCExtensionsApp)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

CIFCExtensionsApp theApp;

BOOL CIFCExtensionsApp::InitInstance()
{
	_Module.Init(ObjectMap);
    return CWinApp::InitInstance();
}

int CIFCExtensionsApp::ExitInstance()
{
    _Module.Term();
    return CWinApp::ExitInstance();
}
