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
CLASS
   CIfcTargetHints

   Hints for the author of a mapping table: when a target isn't found, places in
   the model that look like they might hold it. Hints only go to the log. They
   are never used as values (design decision G3, devdocs/MappingTablesDesign.md).
*****************************************************************************/

#include "IfcTargets.h"

class CIfcTargetHints
{
public:
   // Places on the element that may hold the target, with their values, e.g. "IaDOT_PPCB.2_Type = \"BTB45\"".
   // Empty if there are none, or the target has no hint rule
   static std::vector<std::string> Find(const TargetDef& target, IfcSchema::IfcObject object);

   // Logs elements that look like haunches (ObjectType or Name contains "haunch") when the mapping table
   // has no haunch element role or it selects none. excluded_ids are elements that aren't hinted (e.g. the deck)
   static void LogHaunchHints(ifcopenshell::file& file, const std::set<int>& excluded_ids);
};
