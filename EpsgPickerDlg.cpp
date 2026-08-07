///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright � 1999-2026  Washington State Department of Transportation
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

#include "stdafx.h"
#include "afxdialogex.h"
#include "EpsgPickerDlg.h"

#include <algorithm>

IMPLEMENT_DYNAMIC(CEpsgPickerDlg, CDialog)

CEpsgPickerDlg::CEpsgPickerDlg(const CString& strCaption, const std::vector<EpsgCrsInfo>& crsList, const CString& strCurrentCode, const EpsgCrsInfo* pReferenceCrs, CWnd* pParent)
	: CDialog(IDD_EPSG_PICKER, pParent), m_strCaption(strCaption), m_CrsList(crsList), m_strCurrentCode(strCurrentCode), m_pReferenceCrs(pReferenceCrs)
{
}

CEpsgPickerDlg::~CEpsgPickerDlg()
{
}

void CEpsgPickerDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_EPSG_SEARCH, m_ctrlSearch);
	DDX_Control(pDX, IDC_EPSG_LIST, m_ctrlList);
}

BEGIN_MESSAGE_MAP(CEpsgPickerDlg, CDialog)
	ON_EN_CHANGE(IDC_EPSG_SEARCH, OnSearchChanged)
	ON_NOTIFY(LVN_ITEMCHANGED, IDC_EPSG_LIST, OnListItemChanged)
	ON_NOTIFY(NM_DBLCLK, IDC_EPSG_LIST, OnListDblClick)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_EPSG_LIST, OnListCustomDraw)
END_MESSAGE_MAP()

BOOL CEpsgPickerDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	if (!m_strCaption.IsEmpty())
		SetWindowText(m_strCaption);

	m_ctrlList.SetExtendedStyle(m_ctrlList.GetExtendedStyle() | LVS_EX_FULLROWSELECT);
	m_ctrlList.InsertColumn(0, _T("EPSG Code"), LVCFMT_LEFT, 70);
	m_ctrlList.InsertColumn(1, _T("Name"), LVCFMT_LEFT, 260);
	m_ctrlList.InsertColumn(2, _T("Area of Use"), LVCFMT_LEFT, 260);

	// The gray/relevant color coding (and its legend) only means anything once we have a
	// reference CRS to compare area-of-use against.
	GetDlgItem(IDC_EPSG_LEGEND)->ShowWindow(m_pReferenceCrs ? SW_SHOW : SW_HIDE);

	RefreshList(true);

	if (m_SelectedIndex >= 0)
	{
		// Give the list control focus so the pre-selected current setting renders with the
		// active (blue) highlight immediately, instead of the easy-to-miss unfocused-selection
		// gray it would get if focus stayed on the search box by default.
		m_ctrlList.SetFocus();
		return FALSE; // we set focus ourselves
	}

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void CEpsgPickerDlg::RefreshList(bool selectCurrentCode)
{
	CString strFilter;
	m_ctrlSearch.GetWindowText(strFilter);
	strFilter.MakeLower();

	std::vector<int> matches;
	for (int i = 0; i < (int)m_CrsList.size(); i++)
	{
		const auto& crs = m_CrsList[i];
		if (!strFilter.IsEmpty())
		{
			CString code(crs.Code), name(crs.Name), area(crs.AreaName);
			code.MakeLower();
			name.MakeLower();
			area.MakeLower();
			if (code.Find(strFilter) < 0 && name.Find(strFilter) < 0 && area.Find(strFilter) < 0)
				continue;
		}
		matches.push_back(i);
	}

	if (m_pReferenceCrs)
	{
		// stable_sort: entries whose area of use overlaps the reference CRS move to the top,
		// preserving the existing (alphabetical-by-code) order within each group.
		std::stable_sort(matches.begin(), matches.end(), [this](int a, int b)
			{
				bool overlapsA = EpsgAreaOfUseOverlaps(m_CrsList[a], *m_pReferenceCrs);
				bool overlapsB = EpsgAreaOfUseOverlaps(m_CrsList[b], *m_pReferenceCrs);
				return overlapsA && !overlapsB;
			});
	}

	m_ctrlList.DeleteAllItems();
	int rowToSelect = -1;
	for (int i : matches)
	{
		const auto& crs = m_CrsList[i];
		int row = m_ctrlList.InsertItem(m_ctrlList.GetItemCount(), crs.Code);
		m_ctrlList.SetItemText(row, 1, crs.Name);
		m_ctrlList.SetItemText(row, 2, crs.AreaName);
		m_ctrlList.SetItemData(row, (DWORD_PTR)i);

		if (selectCurrentCode && rowToSelect < 0 && !m_strCurrentCode.IsEmpty() && crs.Code.CompareNoCase(m_strCurrentCode) == 0)
			rowToSelect = row;
	}

	m_SelectedIndex = -1;
	GetDlgItem(IDOK)->EnableWindow(FALSE);

	if (rowToSelect >= 0)
	{
		// Triggers LVN_ITEMCHANGED, which sets m_SelectedIndex and re-enables OK via OnListItemChanged.
		m_ctrlList.SetItemState(rowToSelect, LVIS_SELECTED | LVIS_FOCUSED, LVIS_SELECTED | LVIS_FOCUSED);
		m_ctrlList.EnsureVisible(rowToSelect, FALSE);
	}
}

void CEpsgPickerDlg::OnSearchChanged()
{
	RefreshList();
}

void CEpsgPickerDlg::OnListItemChanged(NMHDR* pNMHDR, LRESULT* pResult)
{
	NM_LISTVIEW* pNMListView = (NM_LISTVIEW*)pNMHDR;
	if ((pNMListView->uChanged & LVIF_STATE) != 0 && (pNMListView->uNewState & LVIS_SELECTED) != 0)
	{
		m_SelectedIndex = (int)m_ctrlList.GetItemData(pNMListView->iItem);
		GetDlgItem(IDOK)->EnableWindow(TRUE);
	}
	*pResult = 0;
}

void CEpsgPickerDlg::OnListDblClick(NMHDR* pNMHDR, LRESULT* pResult)
{
	if (m_SelectedIndex >= 0)
		OnOK();
	*pResult = 0;
}

void CEpsgPickerDlg::OnListCustomDraw(NMHDR* pNMHDR, LRESULT* pResult)
{
	NMLVCUSTOMDRAW* pLVCD = reinterpret_cast<NMLVCUSTOMDRAW*>(pNMHDR);

	switch (pLVCD->nmcd.dwDrawStage)
	{
	case CDDS_PREPAINT:
		*pResult = CDRF_NOTIFYITEMDRAW; // ask for a per-item notification so we can color rows
		return;

	case CDDS_ITEMPREPAINT:
		if (m_pReferenceCrs)
		{
			int index = (int)pLVCD->nmcd.lItemlParam; // set via SetItemData() in RefreshList
			if (!EpsgAreaOfUseOverlaps(m_CrsList[index], *m_pReferenceCrs))
				pLVCD->clrText = RGB(128, 128, 128); // gray = outside the reference CRS's area of use
		}
		*pResult = CDRF_DODEFAULT;
		return;

	default:
		*pResult = CDRF_DODEFAULT;
		return;
	}
}

void CEpsgPickerDlg::OnOK()
{
	if (m_SelectedIndex < 0)
		return; // shouldn't happen since IDOK is disabled without a selection

	m_SelectedCrs = m_CrsList[m_SelectedIndex];
	CDialog::OnOK();
}
