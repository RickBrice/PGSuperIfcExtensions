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

#include <map>
#include <string>
#include <PsgLib\Keys.h>

#include "Include\IfcExport\IfcExportOptions.h"

namespace WBFL { namespace EAF { class Broker; }; };
class IIfcExportExtension;

///////////////////////////////////////////////////////////////////////////
// CIfcExporter
class CIfcExporter
{
public:
    CIfcExporter(void);
    ~CIfcExporter(void);

    // Throws CIfcMappingTableException if the mapping table can't be used.
    // extension: an export extension (see Include\IfcExport\IfcModelBuilder.h), or nullptr for the standard export.
    // With an extension, its targets and element roles are available to the mapping table and its hooks add to the model.
    // Throws CIfcExtensionException if the extension can't be used or one of its hooks fails
    bool BuildModel(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CString& strFilePath, IIfcExportExtension* extension = nullptr);

    // The mapping table files of the last BuildModel, the selected table first, e.g. for the export log
    const std::vector<std::string>& GetMappingTableFiles() const { return m_MappingTableFiles; }

    // The lines the export extension's hooks logged in the last BuildModel
    const std::vector<std::string>& GetExtensionLog() const { return m_ExtensionLog; }

private:
    std::vector<std::string> m_MappingTableFiles;
    std::vector<std::string> m_ExtensionLog;

    template <typename Schema>
    bool BuildModel(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CString& strFilePath, IIfcExportExtension* extension);
};


template <typename Schema>
typename Schema::IfcCartesianPoint ConvertPoint(hierarchy_helper<Schema>& file, IPoint2d* pPoint, bool bMirror = false)
{
   Float64 x, y;
   pPoint->Location(&x, &y);
   x = IsZero(x) ? 0.0 : x;
   y = IsZero(y) ? 0.0 : y;
   return file.create<typename Schema::IfcCartesianPoint>().initialize(std::vector<double>{bMirror ? -x : x, y});
}
