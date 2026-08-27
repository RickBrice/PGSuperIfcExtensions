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


template <typename Schema>
typename Schema::IfcElementQuantity Create_Qto_BeamBaseQuantities(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, const CSegmentKey& segmentKey)
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

   typename Schema::IfcConversionBasedUnit big_area_unit;
   typename Schema::IfcConversionBasedUnit small_area_unit;
   typename Schema::IfcConversionBasedUnit volume_unit;
   typename Schema::IfcConversionBasedUnit mass_unit;
   typename Schema::IfcConversionBasedUnit length_unit;

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


   std::vector<typename Schema::IfcPhysicalQuantity> quantities;

   quantities.push_back(file.create<typename Schema::IfcQuantityLength>().initialize(std::string("Length"), std::nullopt, length_unit, L, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("CrossSectionArea"), std::nullopt, small_area_unit, A, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("OuterSurfaceArea"), std::nullopt, big_area_unit, OSA, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("GrossSurfaceArea"), std::nullopt, big_area_unit, GSA, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("NetSurfaceArea"), std::nullopt, big_area_unit, NSA, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("GrossVolume"), std::nullopt, volume_unit, GV, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("NetVolume"), std::nullopt, volume_unit, NV, std::nullopt)); // optional per usBridge
   quantities.push_back(file.create<typename Schema::IfcQuantityWeight>().initialize(std::string("GrossWeight"), std::nullopt, mass_unit, W, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityWeight>().initialize(std::string("NetWeight"), std::nullopt, mass_unit, NetMass, std::nullopt)); // optional per usBridge

   auto qto = file.create<typename Schema::IfcElementQuantity>().initialize(ifcopenshell::global_id(), {}, std::string("Qto_BeamBaseQuantities"), std::nullopt, std::string("BaseQuantities"), quantities);

   return qto;
}


template <typename Schema>
typename Schema::IfcElementQuantity Create_Qto_ReinforcingElementBaseQuantities(hierarchy_helper<Schema>& file)
{
   //std::vector<typename Schema::IfcPhysicalQuantity> quantities;

   //quantities.push_back(file.create<typename Schema::IfcQuantityCount>().initialize(std::string("Count"), std::nullopt, {}, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityLength>().initialize(std::string("Length"), std::nullopt, length_unit, std::nullopt));
   //quantities.push_back(file.create<typename Schema::IfcQuantityWeight>().initialize(std::string("Weight"), std::nullopt, weight_unit, std::nullopt));

   //auto qto = file.create<typename Schema::IfcElementQuantity>().initialize(ifcopenshell::global_id(), {}, std::string("Qto_ReinforcingElementBaseQuantities"), std::nullopt, std::string("BaseQuantities"), quantities);
   //

   //AddQto(file, segment, qto);
   return {};
}

template <typename Schema>
typename Schema::IfcElementQuantity Create_Qto_SlabBaseQuantatities(hierarchy_helper<Schema>& file)
{
   std::vector<typename Schema::IfcPhysicalQuantity> quantities;

   quantities.push_back(file.create<typename Schema::IfcQuantityLength>().initialize(std::string("Width"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityLength>().initialize(std::string("Length"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityLength>().initialize(std::string("Depth"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityLength>().initialize(std::string("Perimeter"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("GrossArea"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityArea>().initialize(std::string("NetArea"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityVolume>().initialize(std::string("GrossVolume"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityVolume>().initialize(std::string("NetVolume"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityWeight>().initialize(std::string("GrossWeight"), std::nullopt, {}, 0.0, std::nullopt));
   quantities.push_back(file.create<typename Schema::IfcQuantityWeight>().initialize(std::string("NetWeight"), std::nullopt, {}, 0.0, std::nullopt));

   auto qto = file.create<typename Schema::IfcElementQuantity>().initialize(ifcopenshell::global_id(), {}, std::string("Qto_SlabBaseQuantities"), std::nullopt, std::string("BaseQuantities"), quantities);


   return qto;
}