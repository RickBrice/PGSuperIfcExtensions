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
void AddQto(IfcHierarchyHelper<Schema>& file, typename Schema::IfcObjectDefinition* object,typename Schema::IfcElementQuantity* qto)
{
   typename Schema::IfcObjectDefinition::list::ptr related_objects(new typename Schema::IfcObjectDefinition::list);
   related_objects->push(object);

   auto rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_objects, qto);
   file.addEntity(rel_defines_by_properties);
}


template <typename Schema>
void Create_Qto_BeamBaseQuantities(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey, typename Schema::IfcElement* segment)
{
#pragma Reminder("NOTE: These are a little bit dummy quantities - updated in the future")
   // assuming simple sections (no change in cross section or depth like end blocks are variable depth hammerhead segments)
   // need to update PGSuper so we can get the different surface areas directly instead of having to compute them here
   GET_IFACE2(pBroker, IBridge, pBridge);
   GET_IFACE2(pBroker, ISectionProperties, pSectProps);
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   auto releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
   CHECK(vPoi.size() == 1);
   const pgsPointOfInterest& poiMS = vPoi.front();

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   auto L = pBridge->GetSegmentPlanLength(segmentKey);
   auto A = pSectProps->GetAg(releaseIntervalIdx, poiMS);
   auto P = pSectProps->GetPerimeter(poiMS);
   auto OSA = L * P;
   auto GSA = OSA + 2 * A;
   auto GV = L * A;
   auto W = pSectProps->GetSegmentWeight(segmentKey); // this is a unit of mass - see https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/IfcQuantityWeight.htm , WR21 - weight and mass are the same thing in IFC
   auto g = WBFL::Units::System::GetGravitationalAcceleration();
   W /= g;

   typename Schema::IfcConversionBasedUnit* big_area_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* small_area_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* volume_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* mass_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* length_unit = nullptr;

   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      big_area_unit = GetBigAreaUnit<Schema>(file, pBroker);
      small_area_unit = GetSmallAreaUnit<Schema>(file, pBroker);
      volume_unit = GetVolumeUnit<Schema>(file, pBroker);
      mass_unit = GetMassUnit<Schema>(file, pBroker);
      length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      GSA = WBFL::Units::ConvertFromSysUnits(GSA, WBFL::Units::Measure::Feet2);
      GV = WBFL::Units::ConvertFromSysUnits(GV, WBFL::Units::Measure::Feet3);
      L = WBFL::Units::ConvertFromSysUnits(L, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      A = WBFL::Units::ConvertFromSysUnits(A, pDisplayUnits->GetAreaUnit().UnitOfMeasure);
      W = WBFL::Units::ConvertFromSysUnits(W, WBFL::Units::Measure::PoundMass);
   }


   typename Schema::IfcPhysicalQuantity::list::ptr quantities(new typename Schema::IfcPhysicalQuantity::list);

   quantities->push(new typename Schema::IfcQuantityLength(std::string("Length"), boost::none, length_unit, L, boost::none));
   quantities->push(new typename Schema::IfcQuantityArea(std::string("CrossSectionArea"), boost::none, small_area_unit, A, boost::none));
   quantities->push(new typename Schema::IfcQuantityArea(std::string("OuterSurfaceArea"), boost::none, big_area_unit, OSA, boost::none));
   //quantities->push(new typename Schema::IfcQuantityArea(std::string("GrossSurfaceArea"), boost::none, big_area_unit, GSA, boost::none));
   //quantities->push(new typename Schema::IfcQuantityArea(std::string("NetSurfaceArea"), boost::none, big_area_unit, NSA, boost::none));
   //quantities->push(new typename Schema::IfcQuantityArea(std::string("GrossVolume"), boost::none, volume_unit, GV, boost::none));
   //quantities->push(new typename Schema::IfcQuantityArea(std::string("NetVolume"), boost::none, volume_unit, NV, boost::none)); // optional per usBridge
   quantities->push(new typename Schema::IfcQuantityWeight(std::string("GrossWeight"), boost::none, mass_unit, W, boost::none));
   //quantities->push(new typename Schema::IfcQuantityWeight(std::string("NetWeight"), boost::none, mass_unit, NetMass, boost::none)); // optional per usBridge

   auto qto = new typename Schema::IfcElementQuantity(IfcParse::IfcGlobalId(), nullptr, std::string("Qto_BeamBaseQuantities"), boost::none, std::string("BaseQuantities"), quantities);
   file.addEntity(qto);

   AddQto(file,segment, qto);
}


template <typename Schema>
void Create_Qto_ReinforcingElementBaseQuantities(IfcHierarchyHelper<Schema>& file, typename Schema::IfcReinforcingBar* rebar)
{
   //typename Schema::IfcPhysicalQuantity::list::ptr quantities(new typename Schema::IfcPhysicalQuantity::list);

   //quantities->push(new typename Schema::IfcQuantityCount(std::string("Count"), boost::none, nullptr, boost::none));
   //quantities->push(new typename Schema::IfcQuantityLength(std::string("Length"), boost::none, length_unit, boost::none));
   //quantities->push(new typename Schema::IfcQuantityWeight(std::string("Weight"), boost::none, weight_unit, boost::none));

   //auto qto = new typename Schema::IfcElementQuantity(IfcParse::IfcGlobalId(), nullptr, std::string("Qto_ReinforcingElementBaseQuantities"), boost::none, std::string("BaseQuantities"), quantities);
   //file.addEntity(qto);

   //AddQto(file, segment, qto);
}
