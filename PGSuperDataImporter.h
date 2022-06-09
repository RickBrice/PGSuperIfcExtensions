
// PGSuperDataImporter.h : Declaration of the CPGSuperDataImporter

#pragma once

#include <Plugins\PGSuperIEPlugin.h>
#include "resource.h"       // main symbols
#include "IfcAlignmentConverter.h"


/////////////////////////////////////////////////////////////////////////////
// CPGSuperDataImporter
class ATL_NO_VTABLE CPGSuperDataImporter : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPGSuperDataImporter, &CLSID_PGSuperDataImporter>,
   public IPGSDataImporter
{
public:
	CPGSuperDataImporter()
	{
	}

   HRESULT FinalConstruct();

   CBitmap m_Bitmap;

DECLARE_REGISTRY_RESOURCEID(IDR_PGSUPERDATAIMPORTER)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CPGSuperDataImporter)
	COM_INTERFACE_ENTRY(IPGSDataImporter)
END_COM_MAP()

// IPGSDataImporter
public:
   STDMETHOD(Init)(UINT nCmdID) override;
   STDMETHOD(GetMenuText)(/*[out,retval]*/BSTR*  bstrText) const override;
   STDMETHOD(GetBitmapHandle)(/*[out]*/HBITMAP* phBmp) const override;
   STDMETHOD(GetCommandHintText)(BSTR*  bstrText) const override;
   STDMETHOD(Import)(/*[in]*/IBroker* pBroker) override;

private:
   CIfcAlignmentConverter m_IfcConverter;
};
