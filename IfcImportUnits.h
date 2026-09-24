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

#include "Properties.h"

/*****************************************************************************
CLASS
   CIfcImportUnits

   Converts IFC values to PGSuper system units.

   Every IFC unit is reduced to a conversion factor to the fundamental SI units
   (kg, m, s, C, rad), which is how WBFL::Units defines a unit of measure. The
   value is then converted to system units with WBFL::Units::ConvertToSysUnits.

   A value uses its own unit if it has one (e.g. IfcPropertySingleValue.Unit),
   otherwise the project unit for its unit type (IfcProject.UnitsInContext).
*****************************************************************************/
class CIfcImportUnits
{
public:
   using UnitType = IfcSchema::IfcUnitEnum::Value;

   // Reads the project units. Must be called before converting values
   void Init(ifcopenshell::file& file);

   // Returns the conversion factor to fundamental SI units for unit. If unit is empty,
   // the project unit for unit_type is used. If the project doesn't define a unit for
   // unit_type, SI is assumed and the assumption is logged.
   Float64 GetConversionFactor(IfcSchema::IfcUnit unit, UnitType unit_type) const;

   // Returns the conversion factor to fundamental SI units for a specific unit,
   // or nullopt if the unit is not supported
   static std::optional<Float64> GetSIConversionFactor(IfcSchema::IfcUnit unit);

   // Returns the name of the project unit for unit_type (for logging)
   std::string GetProjectUnitName(UnitType unit_type) const;

private:
   std::map<UnitType, std::pair<Float64, std::string>> m_ProjectUnits; // unit type -> conversion factor, unit name
   mutable std::set<UnitType> m_AssumedSI; // unit types without a project unit that have been logged
};

/*****************************************************************************
   IfcMeasureTraits

   Maps an IFC measure type to its unit type and the WBFL::Units unit of measure type
*****************************************************************************/
template <typename Measure> struct IfcMeasureTraits;

#define IFC_MEASURE_TRAITS(_measure_, _unit_type_, _unit_) \
template <> struct IfcMeasureTraits<IfcSchema::_measure_> \
{ \
   static constexpr CIfcImportUnits::UnitType unit_type = IfcSchema::IfcUnitEnum::_unit_type_; \
   using Unit = WBFL::Units::_unit_; \
};

IFC_MEASURE_TRAITS(IfcLengthMeasure, IfcUnit_LENGTHUNIT, Length)
IFC_MEASURE_TRAITS(IfcPositiveLengthMeasure, IfcUnit_LENGTHUNIT, Length)
IFC_MEASURE_TRAITS(IfcNonNegativeLengthMeasure, IfcUnit_LENGTHUNIT, Length)
IFC_MEASURE_TRAITS(IfcPlaneAngleMeasure, IfcUnit_PLANEANGLEUNIT, Angle)
IFC_MEASURE_TRAITS(IfcPositivePlaneAngleMeasure, IfcUnit_PLANEANGLEUNIT, Angle)
IFC_MEASURE_TRAITS(IfcAreaMeasure, IfcUnit_AREAUNIT, Area)
IFC_MEASURE_TRAITS(IfcVolumeMeasure, IfcUnit_VOLUMEUNIT, Volume)
IFC_MEASURE_TRAITS(IfcMassMeasure, IfcUnit_MASSUNIT, Mass)
IFC_MEASURE_TRAITS(IfcForceMeasure, IfcUnit_FORCEUNIT, Force)
IFC_MEASURE_TRAITS(IfcPressureMeasure, IfcUnit_PRESSUREUNIT, Pressure)

#undef IFC_MEASURE_TRAITS

// Converts value, expressed in unit (or the project unit if unit is empty), to system units
template <typename Measure>
Float64 ConvertToSysUnits(const CIfcImportUnits& units, Float64 value, IfcSchema::IfcUnit unit = {})
{
   using Traits = IfcMeasureTraits<Measure>;
   typename Traits::Unit unit_of_measure(units.GetConversionFactor(unit, Traits::unit_type), _T("IFC"));
   return WBFL::Units::ConvertToSysUnits(value, unit_of_measure);
}

// Converts a property value to system units. Returns nullopt if the property doesn't have a value of type Measure
template <typename Measure>
std::optional<Float64> ConvertToSysUnits(const CIfcImportUnits& units, IfcSchema::IfcPropertySingleValue property)
{
   if (!property || !property.NominalValue())
      return std::nullopt;

   auto measure = property.NominalValue().as<Measure>();
   if (!measure)
      return std::nullopt;

   return ConvertToSysUnits<Measure>(units, (Float64)measure, property.Unit());
}

// Returns the named property, or an empty property if not found
IfcSchema::IfcPropertySingleValue FindPropertySingleValue(IfcSchema::IfcObject object, const std::string& pset_name, const std::string& property_name);
IfcSchema::IfcPropertySingleValue FindMaterialPropertySingleValue(IfcSchema::IfcMaterialDefinition material, const std::string& pset_name, const std::string& property_name);

// Returns the value of a measure property in system units, or nullopt if the property isn't found
template <typename Measure>
std::optional<Float64> GetMeasureProperty(const CIfcImportUnits& units, IfcSchema::IfcObject object, const std::string& pset_name, const std::string& property_name)
{
   return ConvertToSysUnits<Measure>(units, FindPropertySingleValue(object, pset_name, property_name));
}

// Returns the value of a material measure property in system units, or nullopt if the property isn't found
template <typename Measure>
std::optional<Float64> GetMaterialMeasureProperty(const CIfcImportUnits& units, IfcSchema::IfcMaterialDefinition material, const std::string& pset_name, const std::string& property_name)
{
   return ConvertToSysUnits<Measure>(units, FindMaterialPropertySingleValue(material, pset_name, property_name));
}
