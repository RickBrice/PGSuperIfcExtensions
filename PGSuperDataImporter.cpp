///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2023  Washington State Department of Transportation
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
#include "PGSuperDataImporter.h"
#include <EAF\EAFAutoProgress.h>
#include <IFace\Project.h>
#include "IfcAlignmentConverter.h"

/////////////////////////////////////////////////////////////////////////////
// CPGSuperDataImporter
HRESULT CPGSuperDataImporter::FinalConstruct()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_Bitmap.LoadBitmap(IDB_BSI));
   return S_OK;
}
STDMETHODIMP CPGSuperDataImporter::Init(UINT nCmdID)
{
   return S_OK;
}

STDMETHODIMP CPGSuperDataImporter::GetMenuText(BSTR*  bstrText) const
{
   *bstrText = CComBSTR("Alignment from IFC File");
   return S_OK;
}

STDMETHODIMP CPGSuperDataImporter::GetBitmapHandle(HBITMAP* phBmp) const
{
   *phBmp = m_Bitmap;
   return S_OK;
}

STDMETHODIMP CPGSuperDataImporter::GetCommandHintText(BSTR*  bstrText) const
{
   *bstrText = CComBSTR("Status line hint text\nTool tip text");
   return S_OK;   
}

STDMETHODIMP CPGSuperDataImporter::Import(IBroker* pBroker)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CFileDialog dlg(TRUE, _T("ifc"),NULL,OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,_T("IFC Files (*.ifc)|*.ifc||"));
   if (dlg.DoModal() == IDOK)
   {
      CString fileName = dlg.GetPathName();

      HRESULT hr = m_IfcConverter.ConvertToPGSuper(pBroker, fileName);
   }
   return S_OK;
}

