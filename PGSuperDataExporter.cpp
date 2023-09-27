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
#include "ExportOptions.h"

#include <EAF/EAFUIIntegration.h>

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
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   // write some bridge data to a text file
   CExportOptions options_dlg;
   if (options_dlg.DoModal() == IDCLOSE)
      return S_OK;

   GET_IFACE2(pBroker, IEAFDocument, pDoc);
   auto file_title = pDoc->GetFileTitle();
   auto file_root = pDoc->GetFileRoot();
   CString default_file_name;
   default_file_name.Format(_T("%s%s.ifc"), file_root, file_title);

	CFileDialog  dlg(FALSE,_T("ifc"),default_file_name,OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, _T("IFC File (*.ifc)|*.ifc||"));
	if (dlg.DoModal() == IDOK)
	{
		CString file_path = dlg.GetPathName();

       CIfcModelBuilder builder;
       bool bResult = builder.BuildModel(pBroker, file_path, options_dlg.options);
       CString strMsg;
       strMsg.Format(_T("Model export %s for %s"), (bResult ? _T("successful") : _T("failed")), file_path);
       AfxMessageBox(strMsg, MB_OK | (bResult ? MB_ICONEXCLAMATION : MB_ICONSTOP));
   }

   return S_OK;
}
