///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright � 1999-2023  Washington State Department of Transportation
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
#include <WBFLCore_i.c>
#include <WBFLCogo_i.c>
#include <WBFLGeometry_i.c>

#include "PGSuperDataImporter.h"
#include "PGSuperDataExporter.h"

#include <IFace\Project.h>
#include <IFace\VersionInfo.h>
#include <IFace\Alignment.h>
#include <IFace\DocumentType.h>

#include <IFace\Bridge.h>
#include <IFace\Intervals.h>

#include <EAF\EAFDisplayUnits.h>

// Build environment setup
// Define the environment variable IFCOPENSHELL_DIR with the root location of IfcOpenShell (e.g. F:\IfcOpenShell)

//#if defined _DEBUG
//#pragma comment(lib,"Debug/IfcParse.lib")
//#else
//#pragma comment(lib,"Release/IfcParse.lib")
//#endif

CComModule _Module;

BEGIN_OBJECT_MAP(ObjectMap)
   OBJECT_ENTRY(CLSID_PGSuperDataImporter,    CPGSuperDataImporter)
   OBJECT_ENTRY(CLSID_PGSuperDataExporter,    CPGSuperDataExporter)
END_OBJECT_MAP()

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
    _Module.Init(ObjectMap, m_hInstance, &LIBID_IFCExtensions);
    return CWinApp::InitInstance();
}

int CIFCExtensionsApp::ExitInstance()
{
    _Module.Term();
    return CWinApp::ExitInstance();
}

/////////////////////////////////////////////////////////////////////////////
// Used to determine whether the DLL can be unloaded by OLE

STDAPI DllCanUnloadNow(void)
{
    AFX_MANAGE_STATE(AfxGetStaticModuleState());
    return (AfxDllCanUnloadNow()==S_OK && _Module.GetLockCount()==0) ? S_OK : S_FALSE;
}

/////////////////////////////////////////////////////////////////////////////
// Returns a class factory to create an object of the requested type

STDAPI DllGetClassObject(REFCLSID rclsid, REFIID riid, LPVOID* ppv)
{
    return _Module.GetClassObject(rclsid, riid, ppv);
}

void RegisterPlugins(bool bRegister)
{
   // Importer/Exporter Plugins

   // PGSuper
   //sysComCatMgr::RegWithCategory(CLSID_PGSuperProjectImporter, CATID_PGSuperProjectImporter, bRegister);
   sysComCatMgr::RegWithCategory(CLSID_PGSuperDataImporter,    CATID_PGSuperDataImporter,    bRegister);
   sysComCatMgr::RegWithCategory(CLSID_PGSuperDataExporter,    CATID_PGSuperDataExporter,    bRegister);

   // PGSplice
   //sysComCatMgr::RegWithCategory(CLSID_PGSpliceProjectImporter, CATID_PGSpliceProjectImporter, bRegister);
   sysComCatMgr::RegWithCategory(CLSID_PGSuperDataImporter, CATID_PGSpliceDataImporter, bRegister);
   sysComCatMgr::RegWithCategory(CLSID_PGSuperDataExporter, CATID_PGSpliceDataExporter, bRegister);
}

/////////////////////////////////////////////////////////////////////////////
// DllRegisterServer - Adds entries to the system registry

STDAPI DllRegisterServer(void)
{
	// registers object, typelib and all interfaces in typelib
	HRESULT hr = _Module.RegisterServer(FALSE);
   if ( FAILED(hr) )
      return hr;

   RegisterPlugins(true);

   return S_OK;
}

/////////////////////////////////////////////////////////////////////////////
// DllUnregisterServer - Removes entries from the system registry

STDAPI DllUnregisterServer(void)
{
   RegisterPlugins(false);
   
   _Module.UnregisterServer();
	return S_OK;
}


