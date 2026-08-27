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
#include "EpsgCatalog.h"

// CEpsgPickerDlg dialog
//
// A searchable picker over an EPSG coordinate reference system list (horizontal or vertical -
// the caller passes in whichever cached list from EpsgCatalog.h is relevant). The search edit
// box filters the list control in-memory as the user types; OK is only enabled once a row is
// selected. Reused for both the "EPSG Code" and "Vertical Datum" browse buttons on
// CGeoReferencingPage.
//
// If pReferenceCrs is supplied (e.g. the horizontal CRS already picked, when this dialog is
// browsing vertical datums), entries whose area of use overlaps its area of use are sorted to
// the top of the list - the full list is still shown below, this is a sort order hint only.
// Entries that do NOT overlap are also custom-drawn in gray (via NM_CUSTOMDRAW on the list
// control) so the relevant ones visually stand out, with IDC_EPSG_LEGEND explaining the color
// coding - both only shown/active when pReferenceCrs is non-null.
//
// If strCurrentCode is supplied (the EPSG code already stored for this field, if any), that row
// is pre-selected and scrolled into view when the dialog opens, so re-opening the picker on an
// already-set field doesn't discard the current selection.

class CEpsgPickerDlg : public CDialog
{
	DECLARE_DYNAMIC(CEpsgPickerDlg)

public:
	CEpsgPickerDlg(const CString& strCaption, const std::vector<EpsgCrsInfo>& crsList, const CString& strCurrentCode = _T(""), const EpsgCrsInfo* pReferenceCrs = nullptr, CWnd* pParent = nullptr);
	virtual ~CEpsgPickerDlg();

	const EpsgCrsInfo& GetSelectedCrs() const { return m_SelectedCrs; }

#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_EPSG_PICKER };
#endif

protected:
	virtual BOOL OnInitDialog() override;
	virtual void DoDataExchange(CDataExchange* pDX) override;
	virtual void OnOK() override;

	afx_msg void OnSearchChanged();
	afx_msg void OnListItemChanged(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg void OnListDblClick(NMHDR* pNMHDR, LRESULT* pResult);
	afx_msg void OnListCustomDraw(NMHDR* pNMHDR, LRESULT* pResult);
	DECLARE_MESSAGE_MAP()

private:
	void RefreshList(bool selectCurrentCode = false);

	CString m_strCaption;
	const std::vector<EpsgCrsInfo>& m_CrsList;
	CString m_strCurrentCode;
	const EpsgCrsInfo* m_pReferenceCrs;
	int m_SelectedIndex = -1; // index into m_CrsList, or -1 if nothing is selected
	EpsgCrsInfo m_SelectedCrs;

	CEdit m_ctrlSearch;
	CListCtrl m_ctrlList;
};
