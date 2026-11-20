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
   A test-only export extension (devdocs/ExportProfilesPlan.md, step 2)

   /IfcTestExtension[=<mode>] on the /IfcExport and /IfcTableToIds command lines
   uses it, so the extension API is tested by this repository's validation
   (Tests/ImportValidation/check_extension.py), not only by BridgeSight's builds.
   It goes through the same SDK as an extension in another DLL:

   - it registers an element role (ifctest.marker) and targets on its own role
     and on the IFC extension's (girder, pier)
   - its hooks record every call, in order, and write a marker element
     (IfcBuildingElementProxy, ObjectType TEST_MARKER) whose property set has
     the calls (Tests/ImportValidation/mappings/TestExtension.ifcmap, which
     "requires" the extension)
   - its modes test the refusals: a hook that throws (throw), another API
     version (api), and another build (fingerprint)
*****************************************************************************/

#include "Include\IfcExport\IfcModelBuilder.h"

class CTestExportExtension : public IIfcExportExtension
{
public:
   enum class Mode
   {
      Hooks,      // every hook records its call; the marker element has them
      Throw,      // OnPier throws: the export stops without writing a file
      ApiVersion, // reports another API version: the extension is refused
      Fingerprint // reports another build fingerprint: the extension is refused
   };

   explicit CTestExportExtension(Mode mode) : m_Mode(mode) {}
   virtual ~CTestExportExtension() = default;

   // The mode of /IfcTestExtension=<mode>: empty or "hooks", "throw", "api", or "fingerprint". Returns false for another name
   static bool GetMode(const CString& name, Mode& mode);

   std::string GetId() const override;
   std::string GetVersion() const override;
   std::string GetPrefix() const override;
   int GetApiVersion() const override;
   std::string GetBuildFingerprint() const override;
   void RegisterContent(IfcExport::IRegistry& registry) override;
   std::unique_ptr<IfcExport::IHooks> CreateHooks() override;

private:
   Mode m_Mode;
};
