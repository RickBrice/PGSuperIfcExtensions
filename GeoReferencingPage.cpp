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

// GeoReferencingPage.cpp : implementation file

#include "stdafx.h"
#include "IfcExtensions.h"
#include "GeoReferencingPage.h"
#include "EpsgCatalog.h"
#include "EpsgPickerDlg.h"

#include <MfcTools/MfcTools.h>

#include <algorithm>

// CGeoReferencingPage dialog

IMPLEMENT_DYNAMIC(CGeoReferencingPage, CPropertyPage)

CGeoReferencingPage::CGeoReferencingPage()
	: CPropertyPage(CGeoReferencingPage::IDD)
{
   m_psp.dwFlags |= PSP_HASHELP | PSP_USEICONID;
   m_psp.pszIcon = MAKEINTRESOURCE(IDI_BSI); // same bSI logo used on the Export IFC Model command
}

CGeoReferencingPage::~CGeoReferencingPage()
{
}

void CGeoReferencingPage::DoDataExchange(CDataExchange* pDX)
{
   USES_CONVERSION;

	CPropertyPage::DoDataExchange(pDX);

   DDX_Text(pDX, IDC_EPSG_CODE, m_GeoRefData.Name);
   DDX_Text(pDX, IDC_DESCRIPTION, m_GeoRefData.Description);
   DDX_Text(pDX, IDC_GEODETIC_DATUM, m_GeoRefData.GeodeticDatum);
   DDX_Text(pDX, IDC_VERTICAL_DATUM, m_GeoRefData.VerticalDatum);
   DDX_Text(pDX, IDC_MAP_PROJECTION, m_GeoRefData.MapProjection);

   DDX_Text(pDX, IDC_EASTINGS, m_GeoRefData.Eastings);
   DDX_Text(pDX, IDC_NORTHINGS, m_GeoRefData.Northings);
   DDX_Text(pDX, IDC_ORTHOGONAL_HEIGHT, m_GeoRefData.OrthogonalHeight);
   DDX_Text(pDX, IDC_XAXIS_ABSCISSA, m_GeoRefData.XAxisAbscissa);
   DDX_Text(pDX, IDC_XAXIS_ORDINATE, m_GeoRefData.XAxisOrdinate);
   DDX_Text(pDX, IDC_SCALE, m_GeoRefData.Scale);
}

BOOL CGeoReferencingPage::OnInitDialog()
{
   CPropertyPage::OnInitDialog();

   if (m_GeoRefData.IsCRSValid)
   {
      static const int crsControlIds[] = { IDC_EPSG_CODE, IDC_DESCRIPTION, IDC_GEODETIC_DATUM, IDC_VERTICAL_DATUM, IDC_MAP_PROJECTION, IDC_BROWSE_EPSG_CODE, IDC_BROWSE_VERTICAL_DATUM };
      for (int id : crsControlIds)
      {
         GetDlgItem(id)->EnableWindow(FALSE);
      }
   }

   static const int mapConversionControlIds[] = { IDC_EASTINGS, IDC_NORTHINGS, IDC_ORTHOGONAL_HEIGHT, IDC_XAXIS_ABSCISSA, IDC_XAXIS_ORDINATE, IDC_SCALE };
   static const int mapConversionLabelIds[] = { IDC_EASTINGS_LABEL, IDC_NORTHINGS_LABEL, IDC_ORTHOGONAL_HEIGHT_LABEL, IDC_XAXIS_ABSCISSA_LABEL, IDC_XAXIS_ORDINATE_LABEL, IDC_SCALE_LABEL };
   if (m_GeoRefData.IsMapConversionValid)
   {
      for (int id : mapConversionControlIds)
      {
         GetDlgItem(id)->EnableWindow(FALSE);
      }
   }
   else
   {
      GetDlgItem(IDC_MAP_CONVERSION_GROUPBOX)->ShowWindow(SW_HIDE);
      for (int id : mapConversionControlIds)
      {
         GetDlgItem(id)->ShowWindow(SW_HIDE);
      }
      for (int id : mapConversionLabelIds)
      {
         GetDlgItem(id)->ShowWindow(SW_HIDE);
      }
   }

   return TRUE;
}

BOOL CGeoReferencingPage::OnSetActive()
{
   BOOL bResult = CPropertyPage::OnSetActive();

   // The dialog manager auto-selects all the text in the first tab-stop edit control whenever
   // this page becomes active - too easy to accidentally overtype/delete the EPSG code, and no
   // other dialog in this UI behaves that way. Post (not send) the deselect so it runs after
   // that default focus/select-all handling rather than being clobbered by it.
   GetDlgItem(IDC_EPSG_CODE)->PostMessage(EM_SETSEL, 0, 0);

   return bResult;
}

void CGeoReferencingPage::OnBrowseEpsgCode()
{
   UpdateData(TRUE); // commit any pending edits before overwriting m_GeoRefData's CRS fields

   CEpsgPickerDlg dlg(_T("Select Coordinate Reference System"), GetHorizontalCrsList(), m_GeoRefData.Name, nullptr, this);
   if (dlg.DoModal() == IDOK)
   {
      const auto& crs = dlg.GetSelectedCrs();
      m_GeoRefData.Name = crs.Code;
      m_GeoRefData.Description = crs.Name;
      m_GeoRefData.MapProjection = crs.ProjectionMethod;

      CString datumName = GetGeodeticDatumName(crs.Code);
      if (!datumName.IsEmpty())
         m_GeoRefData.GeodeticDatum = datumName;

      UpdateData(FALSE);
   }
}

void CGeoReferencingPage::OnBrowseVerticalDatum()
{
   UpdateData(TRUE);

   // If a horizontal CRS has already been picked, sort vertical datums whose area of use
   // overlaps it to the top of the list - purely a sort hint, the full list is still shown.
   const auto& horizontalList = GetHorizontalCrsList();
   auto it = std::find_if(horizontalList.begin(), horizontalList.end(), [this](const EpsgCrsInfo& crs)
      {
         return crs.Code == m_GeoRefData.Name;
      });
   const EpsgCrsInfo* pHorizontalCrs = (it != horizontalList.end()) ? &(*it) : nullptr;

   CEpsgPickerDlg dlg(_T("Select Vertical Datum"), GetVerticalCrsList(), m_GeoRefData.VerticalDatum, pHorizontalCrs, this);
   if (dlg.DoModal() == IDOK)
   {
      m_GeoRefData.VerticalDatum = dlg.GetSelectedCrs().Code;
      UpdateData(FALSE);
   }
}

BEGIN_MESSAGE_MAP(CGeoReferencingPage, CPropertyPage)
   ON_BN_CLICKED(IDC_BROWSE_EPSG_CODE, OnBrowseEpsgCode)
   ON_BN_CLICKED(IDC_BROWSE_VERTICAL_DATUM, OnBrowseVerticalDatum)
END_MESSAGE_MAP()
