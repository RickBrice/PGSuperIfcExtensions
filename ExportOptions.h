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
#include "IfcExporter.h"

// CExportOptions - the "IFC" tab of CExportOptionsSheet (ExportOptionsSheet.h).

class CExportOptions : public CPropertyPage
{
	DECLARE_DYNAMIC(CExportOptions)

public:
	CExportOptions();
	virtual ~CExportOptions();

	CIfcExportOptions options;

	// Round trips: the model the project was imported from (empty if none), why it can't be used as it is (e.g. it changed
	// since the import; empty if it can), and whether the export updates a copy of it. Set before the page is shown; the
	// user can locate the model
	CString m_strSourceModel;
	CString m_strSourceModelProblem;
	bool m_bUpdateSourceModel = false;

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_EXPORT_OPTIONS };
#endif

protected:
	virtual BOOL OnInitDialog() override;
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	void OnSpanChanged();
	afx_msg void OnIncludeRebarClicked();
	afx_msg void OnBrowseSourceModel();

	void UpdateSourceModelControls();

	// Rebar Model radios only mean something when rebar is being exported at all.
	void UpdateRebarRepresentationEnable();

	DECLARE_MESSAGE_MAP()
};
