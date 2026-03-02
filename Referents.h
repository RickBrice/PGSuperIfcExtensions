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
typename Schema::IfcRelNests* GetReferentNest(IfcHierarchyHelper<Schema>& file, typename Schema::IfcAlignment* alignment)
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

   typename Schema::IfcObjectDefinition::list::ptr referents(new typename Schema::IfcObjectDefinition::list);
   auto rel_nests = new typename Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Nests referents"), alignment, referents);
   file.addEntity(rel_nests);
   return rel_nests;
}

template <typename Schema>
void AddReferent(IfcHierarchyHelper<Schema>& file, typename Schema::IfcAlignment* alignment, typename Schema::IfcReferent* referent)
{
   typename Schema::IfcRelNests* nest = GetReferentNest<Schema>(file, alignment);
   auto related_objects = nest->RelatedObjects();
   related_objects->push(referent);
   //std::sort(related_objects->begin(), related_objects->end(),
   //   [](typename Schema::IfcObjectDefinition* obj1, typename Schema::IfcObjectDefinition* obj2)
   //   {
   //      typename Schema::IfcReferent* ref1 = obj1->as<typename Schema::IfcReferent>();
   //      typename Schema::IfcReferent* ref2 = obj2->as<typename Schema::IfcReferent>();
   //      if (ref1 && ref2)
   //      {
   //         typename Schema::IfcReal* value1 = GetProperty<Schema, Schema::IfcReal>(ref1, "Pset_Stationing", "Station");
   //         typename Schema::IfcReal* value2 = GetProperty<Schema, Schema::IfcReal>(ref2, "Pset_Stationing", "Station");
   //         if (value1 && value2)
   //         {
   //            return (double)(*value1) < (double)(*value2);
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
   auto point_on_alignment = new typename Schema::IfcPointByDistanceExpression(
      new typename Schema::IfcLengthMeasure(0.0),
      boost::none, boost::none, boost::none,
      directrix);
   auto relative_placement = new typename Schema::IfcAxis2PlacementLinear(point_on_alignment, nullptr, nullptr);
   auto referent_placement = new typename Schema::IfcLinearPlacement(nullptr, relative_placement, nullptr);

   // Create referent
   auto start_station_referent = new typename Schema::IfcReferent(IfcParse::IfcGlobalId(), nullptr, std::string("Start of alignment station"), boost::none, boost::none, referent_placement, nullptr, Schema::IfcReferentTypeEnum::IfcReferentType_STATION);

   // Define properties for Pset_Stationing
   typename Schema::IfcProperty::list::ptr pset_station_properties(new typename Schema::IfcProperty::list);
   pset_station_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Station"), boost::none, new typename Schema::IfcLengthMeasure(startStation), nullptr));

   // Create Pset and assign properties
   auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_Stationing"), boost::none, pset_station_properties);
   file.addEntity(property_set);

   // Assign the property set to the referent
   typename Schema::IfcObjectDefinition::list::ptr referents(new Schema::IfcObjectDefinition::list);
   referents->push(start_station_referent);

   auto rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Relates start station properties to referent"), boost::none, referents, property_set);
   file.addEntity(rel_defines_by_properties);

   //
   // Nest the referent to alignment
   //
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();
   AddReferent(file, alignment, start_station_referent);
}
