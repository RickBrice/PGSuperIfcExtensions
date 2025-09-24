// CImportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "ImportOptions.h"
#include <MfcTools/CustomDDX.h>

// CImportOptions dialog

IMPLEMENT_DYNAMIC(CImportOptions, CDialog)

CImportOptions::CImportOptions(CWnd* pParent /*=nullptr*/)
	: CDialog(IDD_IMPORT_OPTIONS, pParent)
{

}

CImportOptions::~CImportOptions()
{
}

void CImportOptions::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_RadioEnum<CIfcImportOptions::ModelElements>(pDX, IDC_ALIGNMENT_ONLY, options.model_elements);
}


BEGIN_MESSAGE_MAP(CImportOptions, CDialog)
END_MESSAGE_MAP()


// CImportOptions message handlers
