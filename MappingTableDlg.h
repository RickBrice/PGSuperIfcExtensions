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

// Options > IFC Mapping Table: chooses the mapping table IFC imports and exports use (devdocs/MappingTablesDesign.md, M6).
// The table is loaded and validated when it's chosen, and the setting is only saved for a table that can be used.
class CMappingTableDlg : public CDialog
{
   DECLARE_DYNAMIC(CMappingTableDlg)

public:
   CMappingTableDlg(CWnd* pParent = nullptr);

   CString m_strTable; // the chosen table file, empty for the installed standard table

#ifdef AFX_DESIGN_TIME
   enum { IDD = IDD_MAPPING_TABLE };
#endif

protected:
   virtual BOOL OnInitDialog() override;
   virtual void DoDataExchange(CDataExchange* pDX) override;
   virtual void OnOK() override;

   afx_msg void OnTableTypeChanged();
   afx_msg void OnBrowse();
   afx_msg void OnFileChanged();

   int m_TableType; // 0 = standard, 1 = agency table file

   // Loads the chosen table and shows what it is, or why it can't be used. Returns true if it can be used
   bool CheckTable();
   void UpdateControls();

   DECLARE_MESSAGE_MAP()
};
