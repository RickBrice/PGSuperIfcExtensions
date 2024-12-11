// ExportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "ExportOptions.h"
#include <MfcTools/CustomDDX.h>

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
	DDX_RadioEnum<CIfcModelBuilderOptions::Schema>(pDX, IDC_4X3_ADD2, options.schema);
	DDX_Check_Bool(pDX, IDC_BSDD, options.classify);
	DDX_RadioEnum<CIfcModelBuilderOptions::ModelElements>(pDX, IDC_ALIGNMENT_ONLY, options.model_elements);
	DDX_RadioEnum<CIfcModelBuilderOptions::AlignmentModel>(pDX, IDC_POLYLINE, options.alignment_model);
	DDX_RadioEnum<CIfcModelBuilderOptions::Tangents>(pDX, IDC_POLYLINE_TANGENT, options.tangents);
	DDX_RadioEnum<CIfcModelBuilderOptions::Representations>(pDX, IDC_REPRESENTATION_CURVE, options.representations);
	DDX_RadioEnum<CIfcModelBuilderOptions::SweepProfile>(pDX, IDC_SWEEP_POLYLINE, options.sweep_profile);
	DDX_RadioEnum<CIfcModelBuilderOptions::Railings>(pDX, IDC_WALL, options.railings);
	DDX_Check_Bool(pDX, IDC_WORK_PLAN, options.include_work_plan);
	DDX_Check_Bool(pDX, IDC_INCLUDE_REBAR, options.include_rebar);
	DDX_Check_Bool(pDX, IDC_INCLUDE_CAMBER, options.include_camber);
	DDX_Check_Bool(pDX, IDC_QUANTITIES, options.include_quantities);
}


BEGIN_MESSAGE_MAP(CExportOptions, CDialog)
END_MESSAGE_MAP()


// CExportOptions message handlers
