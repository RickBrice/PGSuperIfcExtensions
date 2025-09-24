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

#include "Properties.h"
#include <IFace\Bridge.h>
#include <PsgLib\GirderLabel.h>

#include <WBFLCogo\CogoHelpers.h>

template <typename Schema>
Schema::IfcRelNests* GetReferentNest(IfcHierarchyHelper<Schema>& file, typename Schema::IfcAlignment* alignment)
{
   auto nests = alignment->IsNestedBy();
   for (auto nest : *nests)
   {
      auto related_objects = nest->RelatedObjects();
      for (auto related_object : *related_objects)
      {
         if (auto referent = related_object->as<Ifc4x3_add2::IfcReferent>())
         {
            return nest;
         }
      }
   }

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr referents(new aggregate_of<typename Schema::IfcObjectDefinition>());
   auto rel_nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Nests referents"), alignment, referents);
   file.addEntity(rel_nests);
   return rel_nests;
}

template <typename Schema>
void AddReferent(IfcHierarchyHelper<Schema>& file, typename Schema::IfcAlignment* alignment, typename Schema::IfcReferent* referent)
{
   auto nest = GetReferentNest<Schema>(file, alignment);
   auto related_objects = nest->RelatedObjects();
   related_objects->push(referent);
   //std::sort(related_objects->begin(), related_objects->end(),
   //   [](typename Schema::IfcObjectDefinition* obj1, typename Schema::IfcObjectDefinition* obj2)
   //   {
   //      auto ref1 = obj1->as<typename Schema::IfcReferent>();
   //      auto ref2 = obj2->as<typename Schema::IfcReferent>();
   //      if (ref1 && ref2)
   //      {
   //         auto value1 = GetProperty<Schema, Schema::IfcReal>(ref1, "Pset_Stationing", "Station");
   //         auto value2 = GetProperty<Schema, Schema::IfcReal>(ref2, "Pset_Stationing", "Station");
   //         if (value1 && value2)
   //         {
   //            return (Float64)(*value1) < (Float64)(*value2);
   //         }
   //      }
   //      return false;
   //   });
   nest->setRelatedObjects(related_objects);
}

template <typename Schema>
void CreateAlignmentStartStationReferent(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   USES_CONVERSION;

   auto directrix = GetAlignmentDirectrix(file, options);

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   // get stationing information
   Float64 startStation, startElevation, startGrade;
   auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   // Referents must be in order so start with the start of alignment referent

   //
   // Referent at start of alignment
   //

   // Referent position
   auto point_on_alignment = new Schema::IfcPointByDistanceExpression(
      new Schema::IfcLengthMeasure(0.0),
      boost::none, boost::none, boost::none,
      directrix);
   auto relative_placement = new Schema::IfcAxis2PlacementLinear(point_on_alignment, nullptr, nullptr);
   auto referent_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, nullptr);

   // Create referent
   auto start_station_referent = new Schema::IfcReferent(IfcParse::IfcGlobalId(), nullptr, std::string("Start of alignment station"), boost::none, boost::none, referent_placement, nullptr, Schema::IfcReferentTypeEnum::IfcReferentType_STATION);

   // Define properties for Pset_Stationing
   typename aggregate_of<typename Schema::IfcProperty>::ptr pset_station_properties(new aggregate_of<typename Schema::IfcProperty>());
   pset_station_properties->push(new Schema::IfcPropertySingleValue(std::string("Station"), boost::none, new Schema::IfcLengthMeasure(startStation), nullptr));

   // Create Pset and assign properties
   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_Stationing"), boost::none, pset_station_properties);
   file.addEntity(property_set);

   // Assign the property set to the referent
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr referents(new aggregate_of<typename Schema::IfcObjectDefinition>());
   referents->push(start_station_referent);

   auto rel_defines_by_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Relates start station properties to referent"), boost::none, referents, property_set);
   file.addEntity(rel_defines_by_properties);

   //
   // Nest the referent to alignment
   //
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();
   AddReferent(file, alignment, start_station_referent);
}
