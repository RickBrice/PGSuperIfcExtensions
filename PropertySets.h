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
#include "bSDD.h"
#include <PsgLib\BridgeDescription2.h>

// Reinforcing bar property sets whose values come from the rebar detailing code (bar marks, bar
// element/use/position, covers, bar shapes). All other property sets are declared in the mapping
// table and written by IfcPropertyWriter.h (devdocs/MappingTablesDesign.md, M3 stage 3).


template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ACI_ReinforcingBarType(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::string mark, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;

   auto size = WBFL::LRFD::RebarPool::GetBarSize(pRebar->GetSize());
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarMark"), BSDD_PROPERTY("BarMark"), file.create<typename Schema::IfcLabel>().initialize(mark.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarMass"), BSDD_PROPERTY("BarMass"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarSize"), BSDD_PROPERTY("BarSize"), file.create<typename Schema::IfcLabel>().initialize(T2A(size.c_str())), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndEndPrep"), BSDD_PROPERTY("EndEndPrep"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StartEndPrep"), BSDD_PROPERTY("StartEndPrep"), typename Schema::IfcValue{}, typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_ReinforcingBarType"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
void addBarDimension(std::string name,std::string uri,double dim,typename Schema::IfcConversionBasedUnit unit, hierarchy_helper<Schema>& file, std::vector<typename Schema::IfcProperty>& list_of_properties)
{
   // ACI 131 says to use IfcLengthMeasure for distances and IfcPlaneAngleMeasure for angles. IfcReal is for nondimensional real values
   // but usBridge has distances as IfcReal.
   //list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(name, uri, file.create<typename Schema::IfcLengthMeasure>().initialize(dim), unit));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(name, uri, file.create<typename Schema::IfcReal>().initialize(dim), unit));
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ACI_BarShape(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::string bend_shape_name,
   double bend_radius, const std::unordered_map<std::string, double>& dimensions)
{
   std::string standard = "ACI 315-99";
   std::string standard_version = "1999";

   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StandardName"), BSDD_PROPERTY("StandardName"), file.create<typename Schema::IfcLabel>().initialize(standard.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("StandardVersion"), BSDD_PROPERTY("StandardVersion"), file.create<typename Schema::IfcLabel>().initialize(standard_version.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BendShapeName"), BSDD_PROPERTY("BendShapeName"), file.create<typename Schema::IfcLabel>().initialize(bend_shape_name.c_str()), typename Schema::IfcUnit{}));

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit dimension_unit;
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      dimension_unit = GetComponentDimUnit<Schema>(file, pBroker);
      bend_radius = WBFL::Units::ConvertFromSysUnits(bend_radius, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
   }
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("DefaultInsideBendRadius"), BSDD_PROPERTY("DefaultInsideBendRadius"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(bend_radius), dimension_unit));

   for (auto & [name, dim] : dimensions)
   {
      auto value = dim;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
      }
      addBarDimension<Schema>(name, "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/prop/" + name, value, dimension_unit, file, list_of_properties);
   }

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_BarShape"), std::nullopt, list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ACI_ReinforcingBar(hierarchy_helper<Schema>& file, std::string element,std::string use,std::string position)
{
   std::vector<typename Schema::IfcProperty> list_of_properties;
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarElement"), BSDD_PROPERTY("BarElement"), file.create<typename Schema::IfcLabel>().initialize(element.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarUse"), BSDD_PROPERTY("BarUse"), file.create<typename Schema::IfcLabel>().initialize(use.c_str()), typename Schema::IfcUnit{}));
   list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BarPosition"), BSDD_PROPERTY("BarPosition"), file.create<typename Schema::IfcLabel>().initialize(position.c_str()), typename Schema::IfcUnit{}));

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_ReinforcingBar"), std::string("https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/class/usBrPset_ACI_ReinforcingBar"), list_of_properties);
   return property_set;
}

template <typename Schema>
typename Schema::IfcPropertySet Create_usBrPset_ReinforcingCover(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options, std::optional<double> top, std::optional<double> side, std::optional<double> bottom, std::optional<double> end)
{
   if (!top && !side && !end && !bottom)
      return {};

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);

   typename Schema::IfcConversionBasedUnit cover_unit;
   auto length_unit = pDisplayUnits->GetComponentDimUnit();
   if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   {
      cover_unit = GetComponentDimUnit<Schema>(file, pBroker);
   }

   std::vector<typename Schema::IfcProperty> list_of_properties;

   if (top)
   {
      double value = *top;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*top, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("TopFaceCover"), BSDD_PROPERTY("TopFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   if (side)
   {
      double value = *side;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*side, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("SideFaceCover"), BSDD_PROPERTY("SideFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   if (end)
   {
      double value = *end;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*end, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("EndFaceCover"), BSDD_PROPERTY("EndFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   if (bottom)
   {
      double value = *bottom;
      if (options.display_units_for_properties && pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      {
         value = WBFL::Units::ConvertFromSysUnits(*bottom, length_unit.UnitOfMeasure);
      }
      list_of_properties.push_back(file.create<typename Schema::IfcPropertySingleValue>().initialize(std::string("BottomFaceCover"), BSDD_PROPERTY("BottomFaceCover"), file.create<typename Schema::IfcPositiveLengthMeasure>().initialize(value), cover_unit));
   }

   auto property_set = file.create<typename Schema::IfcPropertySet>().initialize(ifcopenshell::global_id(), {}, std::string("usBrPset_ACI_ReinforcingCover"), std::nullopt, list_of_properties);
   return property_set;
}
