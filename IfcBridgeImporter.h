///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2025  Washington State Department of Transportation
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

#include "IfcExporter.h"
#include <IFace\Project.h>

class CIfcImporter;

class CIfcBridgeImporter
{
public:
   CIfcBridgeImporter(CIfcImporter& importer);
   bool Import(IfcParse::IfcFile& file);

private:
   CIfcImporter& m_Importer;

   Ifc4x3_add2::IfcBridge* GetBridge(IfcParse::IfcFile& file);
   bool IsValidBridge(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridge* bridge);

   void SetGirderProperties(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc);
   void ImportSlab(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc);
};