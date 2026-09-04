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

#pragma once
#include "afxdialogex.h"
#include "IdsExporter.h"

// CIdsExportOptionsDlg - the "IDS" tab of CExportOptionsSheet (ExportOptionsSheet.h).
// IDS export only ever runs alongside an IFC build, so options.enabled (this page's
// first control) is the sole on/off switch - there is no standalone IDS export.
class CIdsExportOptionsDlg : public CPropertyPage
{
   DECLARE_DYNAMIC(CIdsExportOptionsDlg)

public:
   CIdsExportOptionsDlg();
   virtual ~CIdsExportOptionsDlg();

   CIdsExportOptions options;

#ifdef AFX_DESIGN_TIME
   enum { IDD = IDD_IDS_EXPORT_OPTIONS };
#endif

protected:
   virtual BOOL OnInitDialog() override;
   virtual void DoDataExchange(CDataExchange* pDX) override;

   afx_msg void OnEnabledClicked();
   void UpdateControlStates(); // enable/disable everything else based on options.enabled

   DECLARE_MESSAGE_MAP()
};
