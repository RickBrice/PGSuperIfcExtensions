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

// IdsExportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "resource.h"
#include "IdsExportOptions.h"

#include <MfcTools/CustomDDX.h>

IMPLEMENT_DYNAMIC(CIdsExportOptionsDlg, CPropertyPage)

CIdsExportOptionsDlg::CIdsExportOptionsDlg()
   : CPropertyPage(IDD_IDS_EXPORT_OPTIONS)
{
}

CIdsExportOptionsDlg::~CIdsExportOptionsDlg()
{
}

namespace
{
   // DDX_Check needs a BOOL lvalue; bridge it to the bool option field in both directions.
   void DDX_CheckBool(CDataExchange* pDX, int nIDC, bool& value)
   {
      BOOL b = value ? TRUE : FALSE;
      DDX_Check(pDX, nIDC, b);
      value = (b != FALSE);
   }

   // Controls that only make sense once IDS export is enabled.
   const int g_GatedControls[] = {
      IDC_IDS_ID_NAME, IDC_IDS_ID_GLOBALID,
      IDC_IDS_INC_BEAM_PROPS, IDC_IDS_INC_QUANTITIES, IDC_IDS_INC_CONCRETE_MATERIAL,
      IDC_IDS_INC_STRANDS, IDC_IDS_INC_REBAR
   };
}

void CIdsExportOptionsDlg::DoDataExchange(CDataExchange* pDX)
{
   CPropertyPage::DoDataExchange(pDX);
   DDX_CheckBool(pDX, IDC_IDS_ENABLED, options.enabled);

   DDX_RadioEnum<CIdsExportOptions::BeamId>(pDX, IDC_IDS_ID_NAME, options.beam_id);

   DDX_CheckBool(pDX, IDC_IDS_INC_BEAM_PROPS, options.include_beam_properties);
   DDX_CheckBool(pDX, IDC_IDS_INC_QUANTITIES, options.include_quantities);
   DDX_CheckBool(pDX, IDC_IDS_INC_CONCRETE_MATERIAL, options.include_concrete_material);
   DDX_CheckBool(pDX, IDC_IDS_INC_STRANDS, options.include_strands);
   DDX_CheckBool(pDX, IDC_IDS_INC_REBAR, options.include_rebar);
}

BEGIN_MESSAGE_MAP(CIdsExportOptionsDlg, CPropertyPage)
   ON_BN_CLICKED(IDC_IDS_ENABLED, OnEnabledClicked)
END_MESSAGE_MAP()

BOOL CIdsExportOptionsDlg::OnInitDialog()
{
   CPropertyPage::OnInitDialog();
   UpdateControlStates();
   return TRUE;
}

void CIdsExportOptionsDlg::OnEnabledClicked()
{
   UpdateData(TRUE); // capture the click into options.enabled before using it
   UpdateControlStates();
}

void CIdsExportOptionsDlg::UpdateControlStates()
{
   for (int nIDC : g_GatedControls)
   {
      if (CWnd* pWnd = GetDlgItem(nIDC))
         pWnd->EnableWindow(options.enabled);
   }
}
