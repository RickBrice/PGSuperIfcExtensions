///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2024  Washington State Department of Transportation
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

#include <IFace\Bridge.h>
#include <PgsExt\GirderLabel.h>

#include <WBFLCogo\CogoHelpers.h>

template <typename Schema>
void CreateReferents(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr alignment_referents(new aggregate_of<typename Schema::IfcObjectDefinition>());

   auto directrix = GetAlignmentDirectrix(file, options);

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   // get stationing information
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

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
   file.addEntity(start_station_referent);
   alignment_referents->push(start_station_referent); // add to list of all alignment referents

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

   // now do referents for each pier

   if (options.model_elements == CIfcModelBuilderOptions::ModelElements::AlignmentAndBridge)
   {
      //
   // Referents for pier locations
   //
      GET_IFACE2(pBroker, IBridge, pBridge);
      auto nPiers = pBridge->GetPierCount();
      for (auto pierIdx = 0; pierIdx < nPiers; pierIdx++)
      {
         // referent position
         auto pierStation = pBridge->GetPierStation(pierIdx);

         auto point_on_alignment = new Schema::IfcPointByDistanceExpression(
            new Schema::IfcLengthMeasure(pierStation - startStation),
            boost::none, boost::none, boost::none,
            directrix);
         auto relative_placement = new Schema::IfcAxis2PlacementLinear(point_on_alignment, nullptr, nullptr);
         auto referent_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, nullptr);

         // create referent
         std::ostringstream os;
         os << "Station " << T2A(WBFL::COGO::Station(pierStation).AsString(station_format).c_str()) << " " << T2A(LABEL_PIER_EX(pBridge->IsAbutment(pierIdx), pierIdx));
         auto referent = new Schema::IfcReferent(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, referent_placement, nullptr, Schema::IfcReferentTypeEnum::IfcReferentType_POSITION);
         file.addEntity(referent);
         alignment_referents->push(referent);

         // create and assign Pset_Stationing
         typename aggregate_of<typename Schema::IfcProperty>::ptr pset_station_properties(new aggregate_of<typename Schema::IfcProperty>());
         pset_station_properties->push(new Schema::IfcPropertySingleValue(std::string("Station"), boost::none, new Schema::IfcLengthMeasure(pierStation), nullptr));

         auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_Stationing"), boost::none, pset_station_properties);
         file.addEntity(property_set);

         typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr referents(new aggregate_of<typename Schema::IfcObjectDefinition>());
         referents->push(referent);

         auto rel_defines_by_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Relates pier station properties to referent"), boost::none, referents, property_set);
         file.addEntity(rel_defines_by_properties);
      }
   }

   //
   // Nest referents to alignment
   //
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();
   auto nests_stationing = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, std::string("Nests Referents with Alignment"), boost::none, alignment, alignment_referents);
   file.addEntity(nests_stationing);
}
