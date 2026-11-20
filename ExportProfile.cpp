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
#include "ExportProfile.h"
#include "IfcExportOptions.h"
#include "IfcExtensionRegistry.h"
#include "IfcMappingTable.h"
#include "Include\IfcExport\IfcModelBuilder.h"

CExportProfile CExportProfile::Standard()
{
   return CExportProfile(nullptr);
}

CExportProfile CExportProfile::Extension(IIfcExportExtension& extension)
{
   return CExportProfile(&extension);
}

std::string CExportProfile::Describe() const
{
   return m_pExtension ? "the export extension \"" + m_pExtension->GetId() + "\"" : std::string("the standard export");
}

std::unique_ptr<IfcExport::IHooks> CExportProfile::CreateHooks() const
{
   std::unique_ptr<IfcExport::IHooks> hooks;
   if (m_pExtension)
      hooks = m_pExtension->CreateHooks();

   if (!hooks)
      hooks = std::make_unique<IfcExport::IHooks>(); // every hook does nothing

   return hooks;
}

CIfcExportOptions CExportProfile::Apply(const CIfcExportOptions& options) const
{
   CIfcExportOptions result(options);
   if (m_pExtension)
   {
      if (!result.source_model.IsEmpty())
         throw CIfcExtensionException("An export with an export extension can't be a copy of the model the project was imported from (a round trip).");

      result.mapping_source = MappingTableSource::ExportExtension; // the extension's table, for messages
   }
   return result;
}
