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

/*****************************************************************************
   Export profiles (devdocs/ExportProfilesPlan.md)

   What makes one kind of IFC export: its mapping table, its hooks, and its
   option overrides. Every profile's model is built the same way
   (CIfcExporter::BuildModel) through the same session and hook code:

   - the standard export (usBridge, or any table that extends the standard
     table): no hooks (a null object), the table of the export options
   - an export extension's export (IIfcModelBuilder::BuildModel): the
     extension's targets, element roles, and hooks, and its table
*****************************************************************************/

#include <memory>
#include <string>

class IIfcExportExtension;
class CIfcExportOptions;
namespace IfcExport { class IHooks; }

class CExportProfile
{
public:
   // The standard export: no extension and no hooks. The mapping table is the export options'
   static CExportProfile Standard();

   // An export extension's export: its targets, element roles, and hooks, and the table of the options (the extension's).
   // The extension must outlive the profile
   static CExportProfile Extension(IIfcExportExtension& extension);

   // The export extension, or nullptr for the standard export
   IIfcExportExtension* GetExtension() const { return m_pExtension; }

   // For messages, e.g. "the standard export" or "the export extension "MyAgency.MyExchanger""
   std::string Describe() const;

   // The hooks for one export: the extension's, or hooks that do nothing. Never nullptr
   std::unique_ptr<IfcExport::IHooks> CreateHooks() const;

   // The export options with the profile's overrides. Throws CIfcExtensionException if the options can't be used with the
   // profile (an extension's export can't be a round trip)
   CIfcExportOptions Apply(const CIfcExportOptions& options) const;

private:
   explicit CExportProfile(IIfcExportExtension* pExtension) : m_pExtension(pExtension) {}

   IIfcExportExtension* m_pExtension = nullptr;
};
