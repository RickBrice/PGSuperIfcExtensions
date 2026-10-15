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

// The ACI bar shape property set of a bar type, from the bend shape and its dimensions (see BarDefinition in
// BarShape.h). Its properties depend on the shape, so it's written here. The other property sets of bars are
// declared in the mapping table (devdocs/MappingTablesDesign.md).

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
