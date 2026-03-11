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
#include "ImportResults.h"

// CImportResults dialog

IMPLEMENT_DYNAMIC(CImportResults, CDialog)

CImportResults::CImportResults(std::ostringstream& log,CWnd* pParent /*=nullptr*/)
	: m_Log(log), CDialog(IDD_IMPORT_RESULTS, pParent)
{

}

CImportResults::~CImportResults()
{
}

void CImportResults::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	CString txt(m_Log.str().c_str());
	if (!pDX->m_bSaveAndValidate)
	{
	   txt.Replace(_T("\n"), _T("\r\n"));
	}
	DDX_Text(pDX, IDC_EDIT, txt );
}


BEGIN_MESSAGE_MAP(CImportResults, CDialog)
END_MESSAGE_MAP()


// CImportResults message handlers
