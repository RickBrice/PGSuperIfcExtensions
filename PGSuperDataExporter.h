
// PGSuperExporter.h : Declaration of the CPGSuperExporter

#ifndef __PGSUPEREXPORTER_H_
#define __PGSUPEREXPORTER_H_

#include <Plugins\PGSuperIEPlugin.h>
#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CPGSuperDataExporter
class ATL_NO_VTABLE CPGSuperDataExporter : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPGSuperDataExporter, &CLSID_PGSuperDataExporter>,
   public IPGSDataExporter
{
public:
	CPGSuperDataExporter()
	{
	}

   HRESULT FinalConstruct();

   CBitmap m_Bitmap;

DECLARE_REGISTRY_RESOURCEID(IDR_PGSUPERDATAEXPORTER)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CPGSuperDataExporter)
   COM_INTERFACE_ENTRY(IPGSDataExporter)
END_COM_MAP()

// IPGSDataExporter
public:
   STDMETHOD(Init)(UINT nCmdID) override;
   STDMETHOD(GetMenuText)(/*[out,retval]*/BSTR*  bstrText) const override;
   STDMETHOD(GetBitmapHandle)(/*[out]*/HBITMAP* phBmp) const override;
   STDMETHOD(GetCommandHintText)(BSTR*  bstrText) const override;
   STDMETHOD(Export)(/*[in]*/IBroker* pBroker) override;
};

#endif //__PGSUPEREXPORTER_H_
