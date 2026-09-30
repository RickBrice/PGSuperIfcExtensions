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

#include <IFace/Tools.h>
#include <IFace/PointOfInterest.h>
#include <IFace\Bridge.h>
#include <IFace\AnalysisResults.h>
#include <IFace\Intervals.h>
#include "Units.h"
#include "IfcExporter.h"
#include "Properties.h"

#include <PsgLib\BridgeDescription2.h>

template <typename Schema>
void AddQto(hierarchy_helper<Schema>& file, typename Schema::IfcObjectDefinition object,typename Schema::IfcElementQuantity qto)
{
   if (qto == nullptr)
      return;

   std::vector<typename Schema::IfcObjectDefinition> related_objects;
   related_objects.push_back(object);

   auto rel_defines_by_properties = file.create<typename Schema::IfcRelDefinesByProperties>().initialize(ifcopenshell::global_id(), {}, std::nullopt, std::nullopt, related_objects, qto);

}


// For CIfcExportOptions::RebarRepresentation::Mapped, one IfcReinforcingBar represents
// several physical bars, so it needs a real Count quantity. This is per bar-group (the count
// differs between groups), so it must NOT be registered as a shared/batched property set the way
// the mapping table's "shared" sets are (see CreateSharedPropertySets). The other quantity sets
// come from the mapping table (devdocs/MappingTablesDesign.md, M3); this one stays in code
// because its value comes from the rebar detailing.
template <typename Schema>
typename Schema::IfcElementQuantity Create_Qto_ReinforcingElementGroupQuantities(hierarchy_helper<Schema>& file, IndexType count)
{
   std::vector<typename Schema::IfcPhysicalQuantity> quantities;
   quantities.push_back(file.create<typename Schema::IfcQuantityCount>().initialize(std::string("Count"), std::nullopt, std::nullopt, (int64_t)count, std::nullopt));

   auto qto = file.create<typename Schema::IfcElementQuantity>().initialize(ifcopenshell::global_id(), {}, std::string("Qto_ReinforcingElementBaseQuantities"), std::nullopt, std::string("BaseQuantities"), quantities);

   return qto;
}

