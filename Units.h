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
typename Schema::IfcConversionBasedUnit FindUnitByName(hierarchy_helper<Schema>& file, const std::string& name)
{
   auto units = file.instances_by_type<typename Schema::IfcConversionBasedUnit>();
   for (auto& u : units)
   {
      if (u.Name() == name)
         return u;
   }

   return {};
}

template <typename Schema>
double GetConversionFactor(typename Schema::IfcConversionBasedUnit conversion_based_unit)
{
   auto measure_with_unit = conversion_based_unit.ConversionFactor();
   auto value_component = measure_with_unit.ValueComponent();
   // ValueComponent is an IfcValue SELECT - regardless of which measure type it
   // actually is, the wrapped numeric value is always attribute index 0.
   return (double)value_component.get_attribute_value(0);
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetStressUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("ksi");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetStressUnit().UnitOfMeasure.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(-1/*length*/, 1/*mass*/, -2/*time*/, 0, 0, 0, 0), // pressure = force/area = (force = mass*length*time^-2) / (area = length^2) = mass*length^-1*time^-2
         Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcPressureMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_PASCAL))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetDisplacementUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("inch");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetDeflectionUnit().UnitOfMeasure.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(1/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcLengthMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetXSectionDimUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("inch");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(1/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcLengthMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetComponentDimUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("inch");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetComponentDimUnit().UnitOfMeasure.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(1/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcLengthMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetSpanLengthUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("foot");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(1/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcLengthMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetBigAreaUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("square foot");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      auto cf = WBFL::Units::Measure::Feet2.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(2/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_AREAUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcAreaMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetSmallAreaUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("square inch");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
      auto cf = pDisplayUnits->GetAreaUnit().UnitOfMeasure.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(2/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_AREAUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcAreaMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_AREAUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetVolumeUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("cubic foot");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      auto cf = WBFL::Units::Measure::Feet3.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(3/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcVolumeMeasure>().initialize(cf), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetMassUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("pound");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      auto cf = WBFL::Units::Measure::PoundMass.GetConvFactor(); // converts to base mass units which is KG
      auto cf2 = WBFL::Units::Measure::Gram.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(0/*length*/, 1/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_MASSUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcMassMeasure>().initialize(cf * cf2), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_MASSUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_GRAM))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetForceUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("kip");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      auto cf = WBFL::Units::Measure::Kip.GetConvFactor(); // converts to base force units which is N
      auto cf2 = WBFL::Units::Measure::Newton.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(1/*length*/, 1/*mass*/, -2/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_FORCEUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcForceMeasure>().initialize(cf * cf2), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_NEWTON))
      );
   }
   return unit;
}

template <typename Schema>
typename Schema::IfcConversionBasedUnit GetAngleUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   std::string name("degree");
   typename Schema::IfcConversionBasedUnit unit = FindUnitByName<Schema>(file, name);
   if (!unit)
   {
      auto cf = WBFL::Units::Measure::Degree.GetConvFactor(); // converts to base angle units which is rad
      auto cf2 = WBFL::Units::Measure::Radian.GetConvFactor();
      unit = file.create<typename Schema::IfcConversionBasedUnit>().initialize(
         file.create<typename Schema::IfcDimensionalExponents>().initialize(0/*length*/, 0/*mass*/, 0/*time*/, 0, 0, 0, 0),
         Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT,
         name,
         file.create<typename Schema::IfcMeasureWithUnit>().initialize(file.create<typename Schema::IfcPlaneAngleMeasure>().initialize(cf * cf2), file.create<typename Schema::IfcSIUnit>().initialize(Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN))
      );
   }
   return unit;
}

template <typename Schema>
std::vector<int64_t> GetCompoundPlaneAngleMeasure(Float64 angle_deg)
{
   WBFL::COGO::Angle angle(WBFL::Units::Convert(angle_deg,WBFL::Units::Measure::Degree,WBFL::Units::Measure::Radian));
   auto [d, m, s] = angle.GetDMS();
   return { d, m, static_cast<int64_t>(s) };
}
