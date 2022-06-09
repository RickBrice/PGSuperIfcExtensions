#include "StdAfx.h"
#include "IfcAlignmentConverterException.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CIfcAlignmentConverterException::CIfcAlignmentConverterException(LPCTSTR strWhat) :
   m_strWhat(strWhat)
{
}

CIfcAlignmentConverterException::~CIfcAlignmentConverterException(void)
{
}

LPCTSTR CIfcAlignmentConverterException::What() const
{
   return m_strWhat.c_str();
}
