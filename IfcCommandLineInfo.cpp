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
#include "stdafx.h"
#include "IfcCommandLineInfo.h"

namespace
{
   // if strParam is "<key>=<value>" (key is case insensitive), puts value in strValue and returns true
   bool GetKeyValue(const CString& strParam, LPCTSTR lpszKey, CString& strValue)
   {
      CString strKey(lpszKey);
      strKey += _T("=");
      if (strParam.Left(strKey.GetLength()).CompareNoCase(strKey) != 0)
         return false;

      strValue = strParam.Mid(strKey.GetLength());
      strValue.Trim(_T(" \""));
      return true;
   }

   CString ReplaceExtension(const CString& strFile, LPCTSTR lpszExt)
   {
      int dot = strFile.ReverseFind(_T('.'));
      int slash = std::max(strFile.ReverseFind(_T('\\')), strFile.ReverseFind(_T('/')));
      CString strBase = (slash < dot ? strFile.Left(dot) : strFile);
      return strBase + lpszExt;
   }
}

CIfcCommandLineInfo::CIfcCommandLineInfo() :
   CEAFCommandLineInfo(),
   m_bIfcImport(false),
   m_bIfcExport(false),
   m_bDisplayUnitsForProperties(true)
{
}

CIfcCommandLineInfo::~CIfcCommandLineInfo()
{
}

void CIfcCommandLineInfo::ParseParam(LPCTSTR lpszParam, BOOL bFlag, BOOL bLast)
{
   bool bMyParameter = false;

   if (bFlag)
   {
      CString strParam(lpszParam);
      CString strValue;
      if (GetKeyValue(strParam, _T("IfcImport"), strValue))
      {
         m_strIfcFile = strValue;
         m_bIfcImport = true;
         m_bCommandLineMode = TRUE; // batch run - the application shuts down when we are done
         bMyParameter = true;
      }
      else if (GetKeyValue(strParam, _T("IfcExport"), strValue))
      {
         m_strIfcFile = strValue;
         m_bIfcExport = true;
         m_bCommandLineMode = TRUE; // batch run - the application shuts down when we are done
         bMyParameter = true;
      }
      else if (GetKeyValue(strParam, _T("IfcPropertyUnits"), strValue))
      {
         if (strValue.CompareNoCase(_T("Display")) == 0)
            m_bDisplayUnitsForProperties = true;
         else if (strValue.CompareNoCase(_T("System")) == 0)
            m_bDisplayUnitsForProperties = false;
         else
            m_bError = TRUE;
         bMyParameter = true;
      }
      else if (GetKeyValue(strParam, _T("IfcOut"), strValue))
      {
         m_strOutFile = strValue;
         bMyParameter = true;
      }
      else if (GetKeyValue(strParam, _T("IfcLog"), strValue))
      {
         m_strLogFile = strValue;
         bMyParameter = true;
      }
   }

   if (!bMyParameter)
   {
      CEAFCommandLineInfo::ParseParam(lpszParam, bFlag, bLast);
   }

   if (bLast && m_bIfcImport && m_bIfcExport)
   {
      m_bError = TRUE; // one or the other
      return;
   }

   if (bLast && m_bIfcExport)
   {
      if (m_strIfcFile.IsEmpty() || m_strIfcFile.Right(4).CompareNoCase(_T(".ifc")) != 0 || m_strFileName.IsEmpty())
      {
         m_bError = TRUE;
         return;
      }

      if (m_strLogFile.IsEmpty())
         m_strLogFile = ReplaceExtension(m_strIfcFile, _T(".export.log"));
   }

   if (bLast && m_bIfcImport)
   {
      if (m_strIfcFile.IsEmpty() || (!m_strOutFile.IsEmpty() && m_strOutFile.Right(4).CompareNoCase(_T(".pgs")) != 0))
      {
         m_bError = TRUE;
         return;
      }

      if (m_strOutFile.IsEmpty())
         m_strOutFile = ReplaceExtension(m_strIfcFile, _T(".pgs"));

      if (m_strLogFile.IsEmpty())
         m_strLogFile = ReplaceExtension(m_strOutFile, _T(".log"));
   }
}

CString CIfcCommandLineInfo::GetUsageMessage()
{
   CString strMsg;
   strMsg.Format(_T("Usage: BridgeLink.exe /IfcImport=<model.ifc> <template.pgt> [/IfcOut=<project.pgs>] [/IfcLog=<import.log>]\n\n")
                 _T("<template.pgt> - PGSuper project template used to create the new project (e.g. IfcImportTemplate.pgt)\n")
                 _T("/IfcImport - IFC model to import\n")
                 _T("/IfcOut - PGSuper project file to create. Defaults to the IFC file name with a .pgs extension\n")
                 _T("/IfcLog - Import log file. Defaults to the project file name with a .log extension"));
   return strMsg;
}

CString CIfcCommandLineInfo::GetErrorMessage()
{
   return GetUsageMessage();
}
