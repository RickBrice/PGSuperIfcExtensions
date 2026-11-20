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

namespace units_detail
{
   // How a PGSuper unit (its WBFL tag) is written in IFC. A US unit is a conversion-based unit with the name IFC gives it;
   // a metric unit is an SI unit with a prefix. Checkers convert a property value by its unit's name (ifctester:
   // ifcopenshell.util.unit.convert), not by its conversion factor, so the name must be the unit's own. Another unit is a
   // conversion-based unit named with its tag
   template <typename Schema>
   struct IfcUnitName
   {
      std::string conversion_name; // a conversion-based unit, when not empty
      std::optional<typename Schema::IfcSIPrefix::Value> prefix; // an SI unit otherwise
      typename Schema::IfcSIUnitName::Value si_name = Schema::IfcSIUnitName::IfcSIUnitName_METRE;
   };

   template <typename Schema>
   IfcUnitName<Schema> GetIfcUnitName(const std::string& tag)
   {
      using P = typename Schema::IfcSIPrefix;
      using N = typename Schema::IfcSIUnitName;
      static const std::map<std::string, std::string> us_names{
         {"KSI", "ksi"}, {"in", "inch"}, {"ft", "foot"}, {"in^2", "square inch"}, {"ft^2", "square foot"}, {"ft^3", "cubic foot"},
         {"lbm", "pound"}, {"kip", "kip"}, {"deg", "degree"} };
      static const std::map<std::string, std::pair<std::optional<typename P::Value>, typename N::Value>> si_names{
         {"Pa", {std::nullopt, N::IfcSIUnitName_PASCAL}}, {"kPa", {P::IfcSIPrefix_KILO, N::IfcSIUnitName_PASCAL}}, {"MPa", {P::IfcSIPrefix_MEGA, N::IfcSIUnitName_PASCAL}},
         {"m", {std::nullopt, N::IfcSIUnitName_METRE}}, {"cm", {P::IfcSIPrefix_CENTI, N::IfcSIUnitName_METRE}}, {"mm", {P::IfcSIPrefix_MILLI, N::IfcSIUnitName_METRE}},
         {"m^2", {std::nullopt, N::IfcSIUnitName_SQUARE_METRE}}, {"cm^2", {P::IfcSIPrefix_CENTI, N::IfcSIUnitName_SQUARE_METRE}}, {"mm^2", {P::IfcSIPrefix_MILLI, N::IfcSIUnitName_SQUARE_METRE}},
         {"m^3", {std::nullopt, N::IfcSIUnitName_CUBIC_METRE}},
         {"kg", {P::IfcSIPrefix_KILO, N::IfcSIUnitName_GRAM}},
         {"N", {std::nullopt, N::IfcSIUnitName_NEWTON}}, {"kN", {P::IfcSIPrefix_KILO, N::IfcSIUnitName_NEWTON}},
         {"rad", {std::nullopt, N::IfcSIUnitName_RADIAN}} };

      IfcUnitName<Schema> name;
      if (auto found = si_names.find(tag); found != si_names.end())
      {
         name.prefix = found->second.first;
         name.si_name = found->second.second;
      }
      else
      {
         auto us = us_names.find(tag);
         name.conversion_name = (us == us_names.end() ? tag : us->second);
      }
      return name;
   }

   // A conversion-based unit's factor is to the IFC SI unit of its unit type, which isn't always PGSuper's system unit
   // (mass: IFC's is the gram, PGSuper's the kg), so it comes from WBFL's conversion, not from GetConvFactor()
   inline Float64 ToIfcSIUnit(const WBFL::Units::Length& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Meter); }
   inline Float64 ToIfcSIUnit(const WBFL::Units::Length2& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Meter2); }
   inline Float64 ToIfcSIUnit(const WBFL::Units::Length3& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Meter3); }
   inline Float64 ToIfcSIUnit(const WBFL::Units::Mass& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Gram); }
   inline Float64 ToIfcSIUnit(const WBFL::Units::Force& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Newton); }
   inline Float64 ToIfcSIUnit(const WBFL::Units::Pressure& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Pa); }
   inline Float64 ToIfcSIUnit(const WBFL::Units::Angle& unit) { return WBFL::Units::Convert(1.0, unit, WBFL::Units::Measure::Radian); }
}

// The IFC unit of values in a PGSuper unit (e.g. a display unit), created the first time it's needed. The values must be
// converted to the same unit (see ExportDisplayValue)
template <typename Schema, typename UnitOfMeasure>
typename Schema::IfcNamedUnit GetPGSuperUnit(hierarchy_helper<Schema>& file, typename Schema::IfcUnitEnum::Value unit_type, const UnitOfMeasure& unit)
{
   const auto& unit_tag = unit.UnitTag();
   std::string tag(unit_tag.begin(), unit_tag.end()); // WBFL tags are ASCII
   auto name = units_detail::GetIfcUnitName<Schema>(tag);

   if (name.conversion_name.empty())
   {
      for (auto& u : file.template instances_by_type<typename Schema::IfcSIUnit>())
      {
         if (u.UnitType() == unit_type && u.Name() == name.si_name && u.Prefix() == name.prefix)
            return u;
      }
      return file.template create<typename Schema::IfcSIUnit>().initialize(unit_type, name.prefix, name.si_name);
   }

   for (auto& u : file.template instances_by_type<typename Schema::IfcConversionBasedUnit>())
   {
      if (u.UnitType() == unit_type && u.Name() == name.conversion_name)
         return u;
   }

   using E = typename Schema::IfcUnitEnum;
   using N = typename Schema::IfcSIUnitName;
   auto cf = units_detail::ToIfcSIUnit(unit);
   std::array<int, 3> exponents{ 0, 0, 0 }; // length, mass, time
   typename N::Value si_name = N::IfcSIUnitName_METRE;
   typename Schema::IfcValue factor;
   switch (unit_type)
   {
   case E::IfcUnit_LENGTHUNIT:     exponents = { 1, 0, 0 };   si_name = N::IfcSIUnitName_METRE;        factor = file.template create<typename Schema::IfcLengthMeasure>().initialize(cf); break;
   case E::IfcUnit_AREAUNIT:       exponents = { 2, 0, 0 };   si_name = N::IfcSIUnitName_SQUARE_METRE; factor = file.template create<typename Schema::IfcAreaMeasure>().initialize(cf); break;
   case E::IfcUnit_VOLUMEUNIT:     exponents = { 3, 0, 0 };   si_name = N::IfcSIUnitName_CUBIC_METRE;  factor = file.template create<typename Schema::IfcVolumeMeasure>().initialize(cf); break;
   case E::IfcUnit_MASSUNIT:       exponents = { 0, 1, 0 };   si_name = N::IfcSIUnitName_GRAM;         factor = file.template create<typename Schema::IfcMassMeasure>().initialize(cf); break;
   case E::IfcUnit_FORCEUNIT:      exponents = { 1, 1, -2 };  si_name = N::IfcSIUnitName_NEWTON;       factor = file.template create<typename Schema::IfcForceMeasure>().initialize(cf); break;
   case E::IfcUnit_PRESSUREUNIT:   exponents = { -1, 1, -2 }; si_name = N::IfcSIUnitName_PASCAL;       factor = file.template create<typename Schema::IfcPressureMeasure>().initialize(cf); break;
   case E::IfcUnit_PLANEANGLEUNIT: exponents = { 0, 0, 0 };   si_name = N::IfcSIUnitName_RADIAN;       factor = file.template create<typename Schema::IfcPlaneAngleMeasure>().initialize(cf); break;
   default: ASSERT(false); factor = file.template create<typename Schema::IfcReal>().initialize(cf); break; // a unit type that isn't handled here - add it
   }

   return file.template create<typename Schema::IfcConversionBasedUnit>().initialize(
      file.template create<typename Schema::IfcDimensionalExponents>().initialize(exponents[0], exponents[1], exponents[2], 0, 0, 0, 0),
      unit_type,
      name.conversion_name,
      file.template create<typename Schema::IfcMeasureWithUnit>().initialize(factor, file.template create<typename Schema::IfcSIUnit>().initialize(unit_type, std::nullopt, si_name)));
}

// The units of property values in display units (CIfcExportOptions::display_units_for_properties)
template <typename Schema>
typename Schema::IfcNamedUnit GetStressUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT, pDisplayUnits->GetStressUnit().UnitOfMeasure);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetDisplacementUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetXSectionDimUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, pDisplayUnits->GetXSectionDimUnit().UnitOfMeasure);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetComponentDimUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, pDisplayUnits->GetComponentDimUnit().UnitOfMeasure);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetSpanLengthUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetSmallAreaUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_AREAUNIT, pDisplayUnits->GetAreaUnit().UnitOfMeasure);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetAngleUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, pDisplayUnits->GetAngleUnit().UnitOfMeasure);
}

// PGSuper has no display units for these, so they're always US units (ExportDisplayValue converts to the same ones)
template <typename Schema>
typename Schema::IfcNamedUnit GetBigAreaUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_AREAUNIT, WBFL::Units::Measure::Feet2);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetVolumeUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, WBFL::Units::Measure::Feet3);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetMassUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_MASSUNIT, WBFL::Units::Measure::PoundMass);
}

template <typename Schema>
typename Schema::IfcNamedUnit GetForceUnit(hierarchy_helper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   return GetPGSuperUnit<Schema>(file, Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, WBFL::Units::Measure::Kip);
}

template <typename Schema>
std::vector<int64_t> GetCompoundPlaneAngleMeasure(Float64 angle_deg)
{
   WBFL::COGO::Angle angle(WBFL::Units::Convert(angle_deg,WBFL::Units::Measure::Degree,WBFL::Units::Measure::Radian));
   auto [d, m, s] = angle.GetDMS();
   return { d, m, static_cast<int64_t>(s) };
}

// Values without their own unit use the project unit for their unit type. The project unit assignment
// starts with the fundamental units only. This adds the SI project unit (system units) for every other
// unit type that is used by a property or quantity without its own unit, and nothing more, so the model
// doesn't declare units it doesn't use. Call after the model is built.
template <typename Schema>
void AddUsedProjectUnits(hierarchy_helper<Schema>& file)
{
   using UnitType = typename Schema::IfcUnitEnum::Value;
   using DerivedUnitType = typename Schema::IfcDerivedUnitEnum::Value;
   std::set<UnitType> used;
   std::set<DerivedUnitType> used_derived;

   auto use_measure = [&used, &used_derived](const std::string& measure)
      {
         static const std::map<std::string, UnitType> unit_types{
            {"IfcLengthMeasure", Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
            {"IfcPositiveLengthMeasure", Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
            {"IfcNonNegativeLengthMeasure", Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
            {"IfcPlaneAngleMeasure", Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT},
            {"IfcPositivePlaneAngleMeasure", Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT},
            {"IfcAreaMeasure", Schema::IfcUnitEnum::IfcUnit_AREAUNIT},
            {"IfcVolumeMeasure", Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT},
            {"IfcMassMeasure", Schema::IfcUnitEnum::IfcUnit_MASSUNIT},
            {"IfcTimeMeasure", Schema::IfcUnitEnum::IfcUnit_TIMEUNIT},
            {"IfcForceMeasure", Schema::IfcUnitEnum::IfcUnit_FORCEUNIT},
            {"IfcPressureMeasure", Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT},
            {"IfcQuantityLength", Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT},
            {"IfcQuantityArea", Schema::IfcUnitEnum::IfcUnit_AREAUNIT},
            {"IfcQuantityVolume", Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT},
            {"IfcQuantityWeight", Schema::IfcUnitEnum::IfcUnit_MASSUNIT},
            {"IfcQuantityTime", Schema::IfcUnitEnum::IfcUnit_TIMEUNIT},
         };
         // measures of derived units (e.g. values of export extension targets). derived_units, below, says what they're made of
         static const std::map<std::string, DerivedUnitType> derived_unit_types{
            {"IfcLinearForceMeasure", Schema::IfcDerivedUnitEnum::IfcDerivedUnit_LINEARFORCEUNIT},
            {"IfcMassDensityMeasure", Schema::IfcDerivedUnitEnum::IfcDerivedUnit_MASSDENSITYUNIT},
            {"IfcTorqueMeasure", Schema::IfcDerivedUnitEnum::IfcDerivedUnit_TORQUEUNIT},
            {"IfcModulusOfElasticityMeasure", Schema::IfcDerivedUnitEnum::IfcDerivedUnit_MODULUSOFELASTICITYUNIT},
         };
         // measures without units don't need a project unit
         static const std::set<std::string> unitless{
            "IfcRatioMeasure",
            "IfcPositiveRatioMeasure",
            "IfcNormalisedRatioMeasure",
            "IfcCountMeasure",
            "IfcNumericMeasure",
            "IfcContextDependentMeasure",
            "IfcDescriptiveMeasure",
         };
         auto found = unit_types.find(measure);
         auto found_derived = derived_unit_types.find(measure);
         if (found != unit_types.end())
            used.insert(found->second);
         else if (found_derived != derived_unit_types.end())
            used_derived.insert(found_derived->second);
         else
            ASSERT(!measure.ends_with("Measure") || unitless.contains(measure)); // a measure type with units that isn't handled here - add it to unit_types
      };

   for (auto& property : file.template instances_by_type<typename Schema::IfcPropertySingleValue>())
   {
      if (property.NominalValue() && !property.Unit())
         use_measure(property.NominalValue().declaration().name());
   }

   for (auto& property : file.template instances_by_type<typename Schema::IfcPropertyListValue>())
   {
      auto values = property.ListValues();
      if (values && !property.Unit())
      {
         for (auto& value : *values)
            use_measure(value.declaration().name());
      }
   }

   for (auto& quantity : file.template instances_by_type<typename Schema::IfcPhysicalSimpleQuantity>())
   {
      if (!quantity.Unit())
         use_measure(quantity.declaration().name());
   }

   auto unit_assignment = file.template instances_by_type<typename Schema::IfcUnitAssignment>().front();
   auto units = unit_assignment.Units();
   for (auto& unit : units)
   {
      if (auto named_unit = unit.template as<typename Schema::IfcNamedUnit>())
         used.erase(named_unit.UnitType()); // already defined
      else if (auto derived_unit = unit.template as<typename Schema::IfcDerivedUnit>())
         used_derived.erase(derived_unit.UnitType());
   }

   if (used.empty() && used_derived.empty())
      return;

   // PGSuper system units are SI (kg, m, Pa, etc)
   static const std::map<UnitType, std::pair<std::optional<typename Schema::IfcSIPrefix::Value>, typename Schema::IfcSIUnitName::Value>> si_units{
      {Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_METRE}},
      {Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN}},
      {Schema::IfcUnitEnum::IfcUnit_AREAUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE}},
      {Schema::IfcUnitEnum::IfcUnit_VOLUMEUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE}},
      {Schema::IfcUnitEnum::IfcUnit_MASSUNIT, {Schema::IfcSIPrefix::IfcSIPrefix_KILO, Schema::IfcSIUnitName::IfcSIUnitName_GRAM}},
      {Schema::IfcUnitEnum::IfcUnit_TIMEUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_SECOND}},
      {Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_NEWTON}},
      {Schema::IfcUnitEnum::IfcUnit_PRESSUREUNIT, {std::nullopt, Schema::IfcSIUnitName::IfcSIUnitName_PASCAL}},
   };

   auto si_unit = [&file](UnitType unit_type)
      {
         const auto& [prefix, name] = si_units.at(unit_type);
         return file.template create<typename Schema::IfcSIUnit>().initialize(unit_type, prefix, name);
      };

   for (auto unit_type : used)
      units.push_back(si_unit(unit_type));

   // derived units are products of SI units
   static const std::map<DerivedUnitType, std::vector<std::pair<UnitType, int64_t>>> derived_units{
      {Schema::IfcDerivedUnitEnum::IfcDerivedUnit_LINEARFORCEUNIT, {{Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, 1}, {Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, -1}}}, // N/m
      {Schema::IfcDerivedUnitEnum::IfcDerivedUnit_MASSDENSITYUNIT, {{Schema::IfcUnitEnum::IfcUnit_MASSUNIT, 1}, {Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, -3}}}, // kg/m^3
      {Schema::IfcDerivedUnitEnum::IfcDerivedUnit_TORQUEUNIT, {{Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, 1}, {Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, 1}}},        // N m
      {Schema::IfcDerivedUnitEnum::IfcDerivedUnit_MODULUSOFELASTICITYUNIT, {{Schema::IfcUnitEnum::IfcUnit_FORCEUNIT, 1}, {Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, -2}}}, // N/m^2
   };
   for (auto unit_type : used_derived)
   {
      std::vector<typename Schema::IfcDerivedUnitElement> elements;
      for (const auto& [element_type, exponent] : derived_units.at(unit_type))
         elements.push_back(file.template create<typename Schema::IfcDerivedUnitElement>().initialize(si_unit(element_type), exponent));

      units.push_back(file.template create<typename Schema::IfcDerivedUnit>().initialize(elements, unit_type, std::nullopt, std::nullopt));
   }
   unit_assignment.setUnits(units);
}
