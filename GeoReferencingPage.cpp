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

#include <MfcTools/MfcTools.h>

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

BEGIN_MESSAGE_MAP(CGeoReferencingPage, CPropertyPage)
END_MESSAGE_MAP()
