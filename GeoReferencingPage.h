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

#include "resource.h"

// CGeoReferencingPage dialog
//
// The property page CIfcExtensionAgent injects as the first tab of PGSuper's Edit
// Alignment Description dialog (CAlignmentDescriptionDlg). Holds a single EPSG code
// value. Follows the same seed/harvest pattern as PGSuper\ExtensionAgentExample's
// CEditPierPage: the agent constructs this page, assigns m_EPSGCode directly
// (CIfcExtensionAgent::CreatePropertyPage), then reads it back after the user clicks
// OK (CIfcExtensionAgent::OnOK). There's no constructor parameter for this, because -
// unlike CEditPierPage, which is handed a live IEditPierData* to query per-instance
// pier data (girder counts, etc.) - IEditAlignmentData is a dummy interface with
// nothing to query; the EPSG code is the agent's own persisted state, not something
// that comes from the alignment dialog itself.
//
// The tab also carries the bSI logo (IDI_BSI), the same icon used on this plugin's
// "Export IFC Model" command, via PSP_USEICONID/m_psp.pszIcon - see the constructor.

class CGeoReferencingPage : public CPropertyPage
{
	DECLARE_DYNAMIC(CGeoReferencingPage)

public:
	CGeoReferencingPage();
	virtual ~CGeoReferencingPage();

   CString m_EPSGCode; // bound to IDC_EPSG_CODE via DDX_Text

// Dialog Data
	enum { IDD = IDD_GEOREFERENCING };

protected:
	virtual void DoDataExchange(CDataExchange* pDX) override;    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
