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
#include <ifcparse/Ifc4x3_add2.h>

#include <EAF\EAFDisplayUnits.h>

template <typename Schema>
typename Schema::IfcConversionBasedUnit* FindUnitByName(IfcHierarchyHelper<Schema>& file, const std::string& name)
{
   auto units = file.instances_by_type<typename Schema::IfcConversionBasedUnit>();
   for (auto u : *units)
   {
      if (u->Name() == name)
         return u;
   }

   return nullptr;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetStressUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("ksi");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetStressUnit().UnitOfMeasure.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(-1/*length*/, 1/*mass*/, -2/*time*/, 0, 0, 0, 0), // pressure = force/area = (force = mass*length*time^-2) / (area = length^2) = mass*length^-1*time^-2
         Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcPressureMeasure(cf), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_PASCAL))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetDisplacementUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("inch");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetDeflectionUnit().UnitOfMeasure.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(1/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcLengthMeasure(cf), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetSpanLengthUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("foot");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(1/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcLengthMeasure(cf), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetBigAreaUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("square foot");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      auto cf = WBFL::Units::Measure::Feet2.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(2/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_AREAUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcAreaMeasure(cf), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetSmallAreaUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("square inch");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetAreaUnit().UnitOfMeasure.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(2/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_AREAUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcAreaMeasure(cf), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetVolumeUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("cubic foot");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      auto cf = WBFL::Units::Measure::Feet3.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(3/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcVolumeMeasure(cf), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetMassUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("pound");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      auto cf = WBFL::Units::Measure::PoundMass.GetConvFactor(); // converts to base mass units which is KG
      auto cf2 = WBFL::Units::Measure::Gram.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(0/*length*/, 1/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_MASSUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcMassMeasure(cf * cf2), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_MASSUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_GRAM))
      );
   }
   return unit;
}
