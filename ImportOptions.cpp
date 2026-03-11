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

// CImportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "ImportOptions.h"
#include <MfcTools/CustomDDX.h>

// CImportOptions dialog

#pragma Reminder("WORKING HERE - This import options dialog is currently unused. Keeping it as a placeholder for now.")
// When the data importer was changed from a project importer to a true data importer the needs for options
// went away. Now the project importer creates a new project entirely from IFC and the data importer imports data into
// an existing model (and alignment data is the only option for this right now).

IMPLEMENT_DYNAMIC(CImportOptions, CDialog)

CImportOptions::CImportOptions(CWnd* pParent /*=nullptr*/)
	: CDialog(IDD_IMPORT_OPTIONS, pParent)
{

}

CImportOptions::~CImportOptions()
{
}

void CImportOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_RadioEnum<CIfcImportOptions::ModelElements>(pDX, IDC_ALIGNMENT_ONLY, options.model_elements);
}


BEGIN_MESSAGE_MAP(CImportOptions, CDialog)
END_MESSAGE_MAP()


// CImportOptions message handlers
