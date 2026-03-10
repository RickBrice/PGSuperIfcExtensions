#pragma once
#include "afxdialogex.h"
#include "IfcImporter.h"

// CImportResults dialog

class CImportResults : public CDialog
{
	DECLARE_DYNAMIC(CImportResults)

public:
   CImportResults(std::ostringstream& log,CWnd* pParent = nullptr);   // standard constructor
	virtual ~CImportResults();

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_IMPORT_RESULTS };
#endif

protected:
   std::ostringstream& m_Log;

	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
