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

// PGSuperDataExporter.cpp : Implementation of CPGSuperDataExporter
#include "stdafx.h"
#include "IfcExtensions.h"
#include "PGSuperDataExporter.h"
#include "IfcModelBuilder.h"
#include <MFCTools\Prompts.h>

HRESULT CPGSuperDataExporter::FinalConstruct()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   VERIFY(m_Bitmap.LoadBitmap(IDB_BSI));
   return S_OK;
}

/////////////////////////////////////////////////////////////////////////////
// CPGSuperDataExporter

STDMETHODIMP CPGSuperDataExporter::Init(UINT nCmdID)
{
   return S_OK;
}

STDMETHODIMP CPGSuperDataExporter::GetMenuText(BSTR*  bstrText) const
{
   *bstrText = CComBSTR("Bridge Model to IFC");
   return S_OK;
}

STDMETHODIMP CPGSuperDataExporter::GetBitmapHandle(HBITMAP* phBmp) const
{
   *phBmp = m_Bitmap;
   return S_OK;
}

STDMETHODIMP CPGSuperDataExporter::GetCommandHintText(BSTR*  bstrText) const
{
   *bstrText = CComBSTR("Status line hint text\nTool tip text");
   return S_OK;   
}

STDMETHODIMP CPGSuperDataExporter::Export(IBroker* pBroker)
{
   // write some bridge data to a text file
   int alignment_type = AfxRBChoose(_T("Alignment Geometric Representation"), _T("Select a geometric representation of the alignment"), _T("IfcPolyline (3D Wire)\nIfcGradientCurve"),1);
   CIfcModelBuilder::SchemaType schemaType = (CIfcModelBuilder::SchemaType)AfxRBChoose(_T("Format"), _T("Select IFC Format"), _T("IFC 4x3 RC3\nIFC 4x3 RC4"),1);
	CFileDialog  dlg(FALSE,_T("ifc"),_T("PGSuperExport.ifc"),OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, _T("IFC File(*.ifc)|*.ifc||"));
	if (dlg.DoModal() == IDOK)
	{
		CString file_path = dlg.GetPathName();

       CIfcModelBuilder builder;
       builder.BuildModel(pBroker, file_path, schemaType, alignment_type == 0 ? true : false);
   }

   return S_OK;
}
