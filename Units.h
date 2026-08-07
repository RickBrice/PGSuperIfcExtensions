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
#include <EAF\EAFDisplayUnits.h>
#include <CoordGeom\Angle.h>

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
double GetConversionFactor(typename Schema::IfcConversionBasedUnit* conversion_based_unit)
{
   double conversion_factor = 1.0;
   auto measure_with_unit = conversion_based_unit->ConversionFactor();

   // this way used to work until rocksdb support was added.
   // get_attribute_value has a lot more parameters that I don't know how to use
   try
   {
      auto value_component = measure_with_unit->ValueComponent();
      // here we know we're using in-memory so 'nullptr, nullptr, 0' is safe
      conversion_factor = (Float64)(value_component->data().get_attribute_value(nullptr,nullptr,0,0));
      //CHECK(value_component); // not dealing with anything but simple conversion factors
      //auto real = value_component->as<typename Schema::IfcReal>();
      //auto ratio = value_component->as<typename Schema::IfcRatioMeasure>();
      //auto length = value_component->as<typename Schema::IfcLengthMeasure>();
      //auto area = value_component->as<typename Schema::IfcAreaMeasure>();
      //auto volume = value_component->as<typename Schema::IfcVolumeMeasure>();
      //if (real)
      //   conversion_factor = *real;
      //else if (ratio)
      //   conversion_factor = *ratio;
      //else if (length)
      //   conversion_factor = *length;
      //else if (area)
      //   conversion_factor = *area;
      //else if (volume)
      //   conversion_factor = *volume;
      //else
      //   ASSERT(false);
   }
   catch (IfcParse::IfcInvalidTokenException& e)
   {
      // Was expecting something like 
      // #15 = IFCMEASUREWITHUNIT(IFCLENGTHMEASURE(3.28083333333333), #16);
      // where the expected token is IFCLENGTHMEASURE, but instead found something like
      // #15=IFCMEASUREWITHUNIT(3.28083333333333,#16);
      // we'll just get the value and keep going
      TRACE(e.what());
      auto pArgument = measure_with_unit->get("ValueComponent");
      CHECK(pArgument.type() == IfcUtil::Argument_DOUBLE);
      conversion_factor = double(pArgument);
   }
   return conversion_factor;
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
typename Schema::IfcConversionBasedUnit* GetXSectionDimUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("inch");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure.GetConvFactor();
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
typename Schema::IfcConversionBasedUnit* GetComponentDimUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("inch");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetComponentDimUnit().UnitOfMeasure.GetConvFactor();
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

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetForceUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("kip");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      auto cf = WBFL::Units::Measure::Kip.GetConvFactor(); // converts to base force units which is N
      auto cf2 = WBFL::Units::Measure::Newton.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(1/*length*/, 1/*mass*/, -2/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_FORCEUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcForceMeasure(cf * cf2), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_NEWTON))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit* GetAngleUnit(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("degree");
   typename Schema::IfcConversionBasedUnit* unit = FindUnitByName<Schema>(file, name);
   if (unit == nullptr)
   {
      auto cf = WBFL::Units::Measure::Degree.GetConvFactor(); // converts to base angle units which is rad
      auto cf2 = WBFL::Units::Measure::Radian.GetConvFactor();
      unit = new typename Schema::IfcConversionBasedUnit(
         new typename Schema::IfcDimensionalExponents(0/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT,
         name,
         new typename Schema::IfcMeasureWithUnit(new typename Schema::IfcPlaneAngleMeasure(cf * cf2), new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN))
      );
   }
   return unit;
}

template <typename Schema>
std::vector<int> GetCompoundPlaneAngleMeasure(Float64 angle_deg)
{
   WBFL::COGO::Angle angle(WBFL::Units::Convert(angle_deg,WBFL::Units::Measure::Degree,WBFL::Units::Measure::Radian));
   auto [d, m, s] = angle.GetDMS();
   return { d, m, static_cast<int>(s) };
}
