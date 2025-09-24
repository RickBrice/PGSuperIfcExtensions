#pragma once
#include "afxdialogex.h"
#include "IfcImporter.h"

// CImportOptions dialog

class CImportOptions : public CDialog
{
	DECLARE_DYNAMIC(CImportOptions)

public:
	CImportOptions(CWnd* pParent = nullptr);   // standard constructor
	virtual ~CImportOptions();

	CIfcImportOptions options;

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_IMPORT_OPTIONS };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
