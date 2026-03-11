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

// PGSuperDataImporter.cpp : Implementation of CPGSuperDataImporter
#include "stdafx.h"
#include "IfcExtensions.h"
#include "IfcImporter.h"
#include "PGSuperDataImporter.h"
#include <EAF/AutoProgress.h>
#include <IFace\Project.h>
#include "ImportOptions.h"

/////////////////////////////////////////////////////////////////////////////
// CPGSuperDataImporter
CPGSuperDataImporter::CPGSuperDataImporter()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_Bitmap.LoadBitmap(IDB_BSI));
}

STDMETHODIMP CPGSuperDataImporter::Init(UINT nCmdID)
{
   return S_OK;
}

CString CPGSuperDataImporter::GetMenuText() const
{
   return CString("Alignment from IFC");
}

HBITMAP CPGSuperDataImporter::GetBitmapHandle() const
{
   return m_Bitmap;
}

CString CPGSuperDataImporter::GetCommandHintText() const
{
   return CString("Status line hint text\nTool tip text");
}

HRESULT CPGSuperDataImporter::Import(std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CFileDialog dlg(TRUE, _T("ifc"),NULL,OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,_T("IFC Files (*.ifc)|*.ifc||"));
   if (dlg.DoModal() == IDOK)
   {
      CString fileName = dlg.GetPathName();

      CIfcImportOptions options;
      options.model_elements = CIfcImportOptions::ModelElements::AlignmentOnly;

      HRESULT hr = CIfcImporter(pBroker).ImportFromIFC(fileName, options);
   }
   return S_OK;
}

