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

#include <filesystem>
#include <fstream>

// CImportResults dialog

IMPLEMENT_DYNAMIC(CImportResults, CDialog)

CImportResults::CImportResults(const std::string& log, const CString& model, CWnd* pParent /*=nullptr*/)
	: CDialog(IDD_IMPORT_RESULTS, pParent), m_Log(log), m_Model(model)
{
}

CImportResults::~CImportResults()
{
}

void CImportResults::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	CString txt(m_Log.c_str());
	if (!pDX->m_bSaveAndValidate)
	{
	   txt.Replace(_T("\n"), _T("\r\n"));
	}
	DDX_Text(pDX, IDC_EDIT, txt );
}


BEGIN_MESSAGE_MAP(CImportResults, CDialog)
   ON_BN_CLICKED(IDC_SAVE_LOG, OnSave)
END_MESSAGE_MAP()


// CImportResults message handlers

void CImportResults::OnSave()
{
   // the default is beside the model, named after it
   CString default_name;
   if (!m_Model.IsEmpty())
   {
      std::filesystem::path model(m_Model.GetString());
      default_name = (model.parent_path() / (model.stem().wstring() + L".import.log")).c_str();
   }

   CFileDialog dlg(FALSE, _T("log"), default_name, OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, _T("Log Files (*.log)|*.log|Text Files (*.txt)|*.txt|All Files (*.*)|*.*||"), this);
   if (dlg.DoModal() != IDOK)
      return;

   std::ofstream file(std::filesystem::path(dlg.GetPathName().GetString()), std::ios::binary);
   if (file)
      file << m_Log;
   if (!file)
      AfxMessageBox(_T("The log can't be written to ") + dlg.GetPathName(), MB_OK | MB_ICONSTOP);
}
