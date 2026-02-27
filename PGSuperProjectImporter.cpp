///////////////////////////////////////////////////////////////////////
// IEPluginExample
// Copyright © 1999-2025  Washington State Department of Transportation
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

// PGSuperProjectImporter.cpp : Implementation of CPGSuperProjectImporter
#include "stdafx.h"
#include "IfcExtensions.h"
#include "PGSuperProjectImporter.h"
#include "IfcImporter.h"
#include <EAF\EAFApp.h>

CPGSuperProjectImporter::CPGSuperProjectImporter()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_Bitmap.LoadBitmap(IDB_BSI));
}

CString CPGSuperProjectImporter::GetItemText() const
{
   return CString("IFC Model Importer");
}

CLSID CPGSuperProjectImporter::GetCLSID() const
{
   return CLSID_PGSuperIfcProjectImporter;
}

HICON CPGSuperProjectImporter::GetIcon() const
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   return AfxGetApp()->LoadIcon(IDI_BSI);
}

CString CPGSuperProjectImporter::GetTemplateFilePath() const
{
   CEAFApp* pApp = EAFGetApp();
   CString strFileName = pApp->GetAppLocation();
   strFileName += CString(_T("IfcImportTemplate.pgt"));
   return strFileName;
}

#include <EAF\EAFUtilities.h>
HRESULT CPGSuperProjectImporter::Import(std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   HRESULT hr = E_FAIL;
   CFileDialog dlg(TRUE, _T("ifc"), NULL, OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST, _T("IFC Files (*.ifc)|*.ifc||"));
   if (dlg.DoModal() == IDOK)
   {
      CString fileName = dlg.GetPathName();

      CIfcImportOptions options;
      options.model_elements = CIfcImportOptions::ModelElements::AlignmentAndBridge;

      hr = CIfcImporter(pBroker).ImportFromIFC(fileName, options);
   }

   return hr;
}
