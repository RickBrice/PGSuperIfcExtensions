// CImportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "ImportOptions.h"
#include <MfcTools/CustomDDX.h>

// CImportOptions dialog

#pragma Reminder("WORKING HERE - This import options dialog is currently unused. Keeping it as a placeholder for now.")
// When the data importer was changed from a project importer to a true data importer the needs for options
// went away. Now the project importer creates a new project entirely from IFC and the data importer imports data into
// an existing model (and alignment data is the only option for this right now).

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
