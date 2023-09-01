#pragma once
#include "afxdialogex.h"
#include "IfcModelBuilder.h"

// CExportOptions dialog

class CExportOptions : public CDialog
{
	DECLARE_DYNAMIC(CExportOptions)

public:
	CExportOptions(CWnd* pParent = nullptr);   // standard constructor
	virtual ~CExportOptions();

	CIfcModelBuilderOptions options;

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_EXPORT_OPTIONS };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

	DECLARE_MESSAGE_MAP()
};
