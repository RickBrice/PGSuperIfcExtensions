
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

