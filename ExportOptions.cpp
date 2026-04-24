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

// ExportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "ExportOptions.h"
#include <MfcTools/CustomDDX.h>
#include <PsgLib\BridgeDescription2.h>
#include <psgLib\GirderLabel.h>

#include <IFace/Tools.h>
#include <IFace\Project.h>

#include <EAF\EAFUtilities.h>

// CExportOptions dialog

IMPLEMENT_DYNAMIC(CExportOptions, CDialog)

CExportOptions::CExportOptions(CWnd* pParent /*=nullptr*/)
	: CDialog(IDD_EXPORT_OPTIONS, pParent)
{

}

CExportOptions::~CExportOptions()
{
}

void CExportOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_RadioEnum<CIfcExportOptions::Schema>(pDX, IDC_4X3_ADD2, options.schema);
	DDX_CBIndex(pDX, IDC_SPAN, (int&)options.girderKey.groupIndex);
   DDX_CBIndex(pDX, IDC_GIRDER, (int&)options.girderKey.girderIndex);
	DDX_Check_Bool(pDX, IDC_BSDD, options.classify);
	DDX_RadioEnum<CIfcExportOptions::ModelElements>(pDX, IDC_ALIGNMENT_ONLY, options.model_elements);
	DDX_RadioEnum<CIfcExportOptions::AlignmentModel>(pDX, IDC_POLYLINE, options.alignment_model);
	DDX_RadioEnum<CIfcExportOptions::Tangents>(pDX, IDC_POLYLINE_TANGENT, options.tangents);
	DDX_RadioEnum<CIfcExportOptions::Representations>(pDX, IDC_REPRESENTATION_CURVE, options.representations);
	DDX_RadioEnum<CIfcExportOptions::SweepProfile>(pDX, IDC_SWEEP_POLYLINE, options.sweep_profile);

	DDX_Check_Bool(pDX, IDC_INCLUDE_REBAR, options.include_rebar);
	DDX_Check_Bool(pDX, IDC_INCLUDE_CAMBER, options.include_camber);
	DDX_Check_Bool(pDX, IDC_QUANTITIES, options.include_quantities);
   DDX_Check_Bool(pDX, IDC_PROPERTY_UNITS, options.display_units_for_properties);

	DDX_RadioEnum<CIfcExportOptions::BeamPlacement>(pDX, IDC_LINEAR_PLACEMENT, options.beam_placement);
	DDX_RadioEnum<CIfcExportOptions::BeamModel>(pDX, IDC_MODEL_SSH, options.beam_model);
}


BEGIN_MESSAGE_MAP(CExportOptions, CDialog)
	ON_CBN_SELCHANGE(IDC_SPAN, OnSpanChanged)
END_MESSAGE_MAP()


// CExportOptions message handlers

BOOL CExportOptions::OnInitDialog()
{
	CDialog::OnInitDialog();

   auto broker = EAFGetBroker();

	GET_IFACE2(broker,IBridgeDescription, pIBridgeDesc);
	const CBridgeDescription2* pBridgeDesc = pIBridgeDesc->GetBridgeDescription();

	CComboBox* pGroupBox = (CComboBox*)GetDlgItem(IDC_SPAN);
	CComboBox* pGirderBox = (CComboBox*)GetDlgItem(IDC_GIRDER);

	GroupIndexType nGroups = pBridgeDesc->GetGirderGroupCount();
	for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
	{
		CString strLabel;
		strLabel.Format(_T("Span %s"), LABEL_SPAN(grpIdx));
		pGroupBox->AddString(strLabel);
	}

	OnSpanChanged();

	/* Intialize each combo selections */
	if (pGroupBox->SetCurSel((int)options.girderKey.groupIndex) == CB_ERR)
		pGroupBox->SetCurSel(0);

	if (pGirderBox->SetCurSel((int)options.girderKey.girderIndex) == CB_ERR)
		pGirderBox->SetCurSel(0);

	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void CExportOptions::OnSpanChanged()
{
	auto broker = EAFGetBroker();

	GET_IFACE2(broker, IBridgeDescription, pIBridgeDesc);
	const CBridgeDescription2* pBridgeDesc = pIBridgeDesc->GetBridgeDescription();

	CComboBox* pGroupBox = (CComboBox*)GetDlgItem(IDC_SPAN);
	CComboBox* pGirderBox = (CComboBox*)GetDlgItem(IDC_GIRDER);

	int group = pGroupBox->GetCurSel();
	if (group == CB_ERR)
		group = 0;

	int girder = pGirderBox->GetCurSel();
	if (girder == CB_ERR)
		girder = 0;

	pGirderBox->ResetContent();

	/* Get count of girders (same number of girders in all spans) */
	GirderIndexType nGirders = pBridgeDesc->GetGirderGroup(group)->GetGirderCount();
	for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
	{
		/* Add current girder string to girder list */
		CString strLabel;
		strLabel.Format(_T("Girder %s"), LABEL_GIRDER(gdrIdx));
		pGirderBox->AddString(strLabel);
	}

	girder = (int)Min(GirderIndexType(girder), nGirders - 1); // don't allow out of bounds if ng decreases between spans

	pGirderBox->SetCurSel(girder);
}
