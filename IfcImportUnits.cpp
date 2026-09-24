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
#include "stdafx.h"
#include "IfcImportUnits.h"
#include "Units.h"

namespace
{
   Float64 GetPrefixFactor(std::optional<IfcSchema::IfcSIPrefix::Value> prefix)
   {
      if (!prefix)
         return 1.0;

      switch (*prefix)
      {
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_EXA:   return 1.0e18;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_PETA:  return 1.0e15;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_TERA:  return 1.0e12;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_GIGA:  return 1.0e9;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_MEGA:  return 1.0e6;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_KILO:  return 1.0e3;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_HECTO: return 1.0e2;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_DECA:  return 1.0e1;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_DECI:  return 1.0e-1;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_CENTI: return 1.0e-2;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_MILLI: return 1.0e-3;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_MICRO: return 1.0e-6;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_NANO:  return 1.0e-9;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_PICO:  return 1.0e-12;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_FEMTO: return 1.0e-15;
      case IfcSchema::IfcSIPrefix::IfcSIPrefix_ATTO:  return 1.0e-18;
      default: CHECK(false); return 1.0;
      }
   }

   std::string GetUnitName(IfcSchema::IfcNamedUnit unit)
   {
      if (auto conversion_based_unit = unit.as<IfcSchema::IfcConversionBasedUnit>())
         return conversion_based_unit.Name();

      if (auto si_unit = unit.as<IfcSchema::IfcSIUnit>())
      {
         std::string name = IfcSchema::IfcSIUnitName::ToString(si_unit.Name());
         return si_unit.Prefix() ? std::string(IfcSchema::IfcSIPrefix::ToString(*si_unit.Prefix())) + name : name;
      }

      return unit.declaration().name();
   }
}

std::optional<Float64> CIfcImportUnits::GetSIConversionFactor(IfcSchema::IfcUnit unit)
{
   if (!unit)
      return std::nullopt;

   if (auto si_unit = unit.as<IfcSchema::IfcSIUnit>())
   {
      Float64 prefix = GetPrefixFactor(si_unit.Prefix());
      switch (si_unit.Name())
      {
      case IfcSchema::IfcSIUnitName::IfcSIUnitName_SQUARE_METRE: return prefix * prefix;
      case IfcSchema::IfcSIUnitName::IfcSIUnitName_CUBIC_METRE:  return prefix * prefix * prefix;
      case IfcSchema::IfcSIUnitName::IfcSIUnitName_GRAM:         return prefix * 1.0e-3; // fundamental mass unit is kg
      default:                                                   return prefix; // all other SI units are fundamental or derived from fundamental units
      }
   }

   if (auto conversion_based_unit = unit.as<IfcSchema::IfcConversionBasedUnit>())
   {
      // conversion based units are defined as a multiple of another unit
      auto factor = GetSIConversionFactor(conversion_based_unit.ConversionFactor().UnitComponent());
      if (!factor)
         return std::nullopt;

      return ::GetConversionFactor<IfcSchema>(conversion_based_unit) * (*factor);
   }

   if (auto derived_unit = unit.as<IfcSchema::IfcDerivedUnit>())
   {
      // derived units are the product of other units raised to a power
      Float64 factor = 1.0;
      for (auto& element : derived_unit.Elements())
      {
         auto element_factor = GetSIConversionFactor(element.Unit());
         if (!element_factor)
            return std::nullopt;

         factor *= pow(*element_factor, (Float64)element.Exponent());
      }
      return factor;
   }

   return std::nullopt; // monetary units, etc
}

void CIfcImportUnits::Init(ifcopenshell::file& file)
{
   m_ProjectUnits.clear();
   m_AssumedSI.clear();

   IfcSchema::IfcUnitAssignment unit_assignment;
   auto projects = file.instances_by_type<IfcSchema::IfcProject>();
   if (!projects.empty())
      unit_assignment = projects.front().UnitsInContext();

   if (!unit_assignment)
   {
      WBFL::System::Logger::Info(_T("IfcProject does not have a unit assignment. SI units are assumed."));
      return;
   }

   for (auto& unit : unit_assignment.Units())
   {
      auto named_unit = unit.as<IfcSchema::IfcNamedUnit>();
      if (!named_unit)
         continue; // derived and monetary project units are not needed

      auto factor = GetSIConversionFactor(unit);
      if (factor)
      {
         m_ProjectUnits[named_unit.UnitType()] = std::make_pair(*factor, GetUnitName(named_unit));
      }
      else
      {
         std::ostringstream os;
         os << "Project unit " << GetUnitName(named_unit) << " is not supported. SI units are assumed for " << IfcSchema::IfcUnitEnum::ToString(named_unit.UnitType());
         WBFL::System::Logger::Info(os.str().c_str());
      }
   }
}

Float64 CIfcImportUnits::GetConversionFactor(IfcSchema::IfcUnit unit, UnitType unit_type) const
{
   if (unit)
   {
      auto factor = GetSIConversionFactor(unit);
      if (factor)
         return *factor;

      std::ostringstream os;
      os << "Unit " << unit.declaration().name() << " is not supported. The project unit is used.";
      WBFL::System::Logger::Info(os.str().c_str());
   }

   auto found = m_ProjectUnits.find(unit_type);
   if (found != m_ProjectUnits.end())
      return found->second.first;

   if (m_AssumedSI.insert(unit_type).second)
   {
      std::ostringstream os;
      os << "The IFC model does not define a project unit for " << IfcSchema::IfcUnitEnum::ToString(unit_type) << ". SI units are assumed.";
      WBFL::System::Logger::Info(os.str().c_str());
   }
   return 1.0;
}

std::string CIfcImportUnits::GetProjectUnitName(UnitType unit_type) const
{
   auto found = m_ProjectUnits.find(unit_type);
   return found == m_ProjectUnits.end() ? std::string("SI") : found->second.second;
}

IfcSchema::IfcPropertySingleValue FindPropertySingleValue(IfcSchema::IfcObject object, const std::string& pset_name, const std::string& property_name)
{
   auto pset = GetPropertySet<IfcSchema>(object, pset_name);
   if (pset)
   {
      for (auto& property : pset.HasProperties())
      {
         if (property.Name() == property_name)
            return property.as<IfcSchema::IfcPropertySingleValue>();
      }
   }
   return {};
}

IfcSchema::IfcPropertySingleValue FindMaterialPropertySingleValue(IfcSchema::IfcMaterialDefinition material, const std::string& pset_name, const std::string& property_name)
{
   auto pset = GetMaterialPropertySet<IfcSchema>(material, pset_name);
   if (pset)
   {
      for (auto& property : pset.Properties())
      {
         if (property.Name() == property_name)
            return property.as<IfcSchema::IfcPropertySingleValue>();
      }
   }
   return {};
}
