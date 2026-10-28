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
#include "IfcImporter.h"

// CImportResults dialog

// The log of the IFC import that made the project (or the project's last IFC import), kept for the session so View > IFC
// Import Log can show it again. It isn't saved with the project
// {4E7A2C91-6B3D-4F58-9A1E-8C2D5B7F0A63}
DEFINE_GUID(IID_IIfcImportLog,
   0x4e7a2c91, 0x6b3d, 0x4f58, 0x9a, 0x1e, 0x8c, 0x2d, 0x5b, 0x7f, 0x0a, 0x63);
class IIfcImportLog
{
public:
   virtual void SetImportLog(const std::string& log, const std::wstring& model) = 0;
   virtual const std::string& GetImportLog() const = 0; // empty if the project wasn't imported in this session
   virtual const std::wstring& GetImportLogModel() const = 0; // the IFC model the import read
};

// Shows an import log, which can be saved to a file
class CImportResults : public CDialog
{
	DECLARE_DYNAMIC(CImportResults)

public:
   // log: the import log. model: the IFC model the import read; the log is saved beside it by default
   CImportResults(const std::string& log, const CString& model, CWnd* pParent = nullptr);
	virtual ~CImportResults();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_IMPORT_RESULTS };
#endif

protected:
   std::string m_Log;
   CString m_Model;

	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
   BOOL OnInitDialog() override;
   afx_msg void OnSave();
   afx_msg void OnDestroy();

   void SaveWindowPosition();
   void RestoreWindowPosition();

	DECLARE_MESSAGE_MAP()
};
