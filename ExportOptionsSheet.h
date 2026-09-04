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
#include "ExportOptions.h"
#include "IdsExportOptions.h"

// CExportOptionsSheet - the single "Bridge Model to IFC" export options dialog. Hosts
// the IFC options and IDS options as two tabs of one property sheet, replacing what
// used to be two dialogs shown one after the other (IFC options, then - only if the
// user opted in - a separate IDS options dialog). IDS export has no standalone entry
// point any more: m_IdsPage.options.enabled (the first control on the IDS tab) is the
// only switch for whether PGSuperDataExporter.cpp writes a .ids file at all.
class CExportOptionsSheet : public CPropertySheet
{
public:
   CExportOptionsSheet(LPCTSTR pszCaption, CWnd* pParentWnd = nullptr);

   CExportOptions       m_IfcPage;
   CIdsExportOptionsDlg m_IdsPage;
};
