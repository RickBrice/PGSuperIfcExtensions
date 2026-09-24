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
#pragma once

#include <EAF\EAFCommandLineInfo.h>

/*****************************************************************************
CLASS
   CIfcCommandLineInfo

   Command line parser for batch IFC import and export.

   Import an IFC model into a new PGSuper project:
   BridgeLink.exe /IfcImport=<model.ifc> <template.pgt> [/IfcOut=<project.pgs>] [/IfcLog=<import.log>]

   The template file is opened by the application as a new project. The IFC model
   is then imported into that project and the project is saved.
   /IfcOut defaults to the IFC file name with a .pgs extension.
   /IfcLog defaults to the output file name with a .log extension.

   Export a PGSuper project to an IFC model:
   BridgeLink.exe /IfcExport=<model.ifc> <project.pgs> [/IfcPropertyUnits=Display|System] [/IfcLog=<export.log>]

   The project is exported with the default export options.
   /IfcPropertyUnits=Display (default) exports property values in display units with their units,
   /IfcPropertyUnits=System exports property values in system units using the project units.
   /IfcLog defaults to the IFC file name with a .export.log extension.

   Put /IfcImport or /IfcExport first. PGSuper builds before the Test and TxDOT agent command line
   fix claim (and reject) any command line whose first parameter is not a flag.
*****************************************************************************/

class CIfcCommandLineInfo : public CEAFCommandLineInfo
{
public:
   CIfcCommandLineInfo();
   virtual ~CIfcCommandLineInfo();

   virtual void ParseParam(LPCTSTR lpszParam, BOOL bFlag, BOOL bLast) override;

   virtual CString GetUsageMessage() override;
   virtual CString GetErrorMessage() override;

   bool m_bIfcImport; // true if /IfcImport was given
   bool m_bIfcExport; // true if /IfcExport was given
   CString m_strIfcFile; // IFC file to import or export
   CString m_strOutFile; // PGSuper project file created by an import
   CString m_strLogFile;
   bool m_bDisplayUnitsForProperties; // export option, see CIfcExportOptions::display_units_for_properties

private:
   // Prevent accidental copying and assignment
   CIfcCommandLineInfo(const CIfcCommandLineInfo&) = delete;
   CIfcCommandLineInfo& operator=(const CIfcCommandLineInfo&) = delete;
};
