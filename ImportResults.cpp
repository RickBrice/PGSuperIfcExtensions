// CImportOptions.cpp : implementation file
//

#include "stdafx.h"
#include "afxdialogex.h"
#include "ImportResults.h"

// CImportResults dialog

IMPLEMENT_DYNAMIC(CImportResults, CDialog)

CImportResults::CImportResults(std::ostringstream& log,CWnd* pParent /*=nullptr*/)
	: m_Log(log), CDialog(IDD_IMPORT_RESULTS, pParent)
{

}

CImportResults::~CImportResults()
{
}

void CImportResults::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	CString txt(m_Log.str().c_str());
	if (!pDX->m_bSaveAndValidate)
	{
	   txt.Replace(_T("\n"), _T("\r\n"));
	}
	DDX_Text(pDX, IDC_EDIT, txt );
}


BEGIN_MESSAGE_MAP(CImportResults, CDialog)
END_MESSAGE_MAP()


// CImportResults message handlers
