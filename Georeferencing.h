///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright � 1999-2026  Washington State Department of Transportation
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

#include "Units.h" // for GetConversionFactor<Schema>, used by ImportMapUnit below

namespace WBFL
{
   namespace System
   {
      class IStructuredSave;
      class IStructuredLoad;
   };
};

struct GeoreferencingData
{
   // For IfcProjectedCRS
   CString Name = _T("2927"); // horizontal datum EPSG code, WITHOUT the "EPSG:" prefix, e.g. "4326"
   CString Description= _T("Washington South (ftUS)");
   CString GeodeticDatum = _T("NAD83(HARN)");
   CString VerticalDatum = _T("5703"); // EPSG code, WITHOUT the "EPSG:" prefix
   CString MapProjection = _T("Lambert Conformal Conic 2SP");
   CString MapZone = _T(""); // IfcProjectedCRS.MapZone - optional, empty means not set

   // Pragmatic representation of IfcProjectedCRS.MapUnit. IFC allows any IfcNamedUnit here, but
   // this exporter/importer pair only ever produces/consumes a plain SI unit (metre) or a single
   // conversion-factor IfcConversionBasedUnit (e.g. "US survey foot") - that's all we model.
   bool IsMapUnitSI = false;                    // false => IfcConversionBasedUnit, true => IfcSIUnit (metre)
   CString MapUnitName = _T("US survey foot");  // human-readable unit name
   Float64 MapUnitToMeters = 1200. / 3937.;      // conversion factor to metres (1.0 when IsMapUnitSI)

   // true only when this data was read from an imported IFC file's IfcProjectedCRS/IfcMapConversion;
   // otherwise it is authored fresh (defaults or user-entered) and IfcMapConversion is computed at export time.
   bool IsCRSValid = false;
   bool IsMapConversionValid = false;
   Float64 Eastings = 1041929.;
   Float64 Northings = 630714.;
   Float64 OrthogonalHeight = 0.;
   Float64 XAxisAbscissa = 1.;
   Float64 XAxisOrdinate = 0.;
   Float64 Scale = 3937. / 1200.; // us survey foot to meter conversion factor

   void Save(WBFL::System::IStructuredSave* pSave) const;
   void Load(WBFL::System::IStructuredLoad* pLoad);
};

// Reads the IfcProjectedCRS.MapUnit attribute (an IfcNamedUnit*) into the pragmatic fields above.
// Leaves georefdata's MapUnit fields untouched if map_unit is null or an unrecognized unit type.
template <typename Schema>
void ImportMapUnit(typename Schema::IfcNamedUnit* map_unit, GeoreferencingData& georefdata)
{
   if (!map_unit)
      return;

   if (auto si_unit = map_unit->template as<typename Schema::IfcSIUnit>())
   {
      if (si_unit->UnitType() == Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT)
      {
         georefdata.IsMapUnitSI = true;
         georefdata.MapUnitName = _T("metre");
         Float64 factor = 1.0;
         auto prefix = si_unit->Prefix();
         if (prefix)
         {
            switch (*prefix)
            {
            case Schema::IfcSIPrefix::IfcSIPrefix_MILLI: factor = 0.001; break;
            case Schema::IfcSIPrefix::IfcSIPrefix_CENTI: factor = 0.01; break;
            case Schema::IfcSIPrefix::IfcSIPrefix_KILO:  factor = 1000.; break;
            default: factor = 1.0; break; // full SI-prefix fidelity is not required for this workflow
            }
         }
         georefdata.MapUnitToMeters = factor;
         return;
      }
   }

   if (auto conversion_based_unit = map_unit->template as<typename Schema::IfcConversionBasedUnit>())
   {
      georefdata.IsMapUnitSI = false;
      georefdata.MapUnitName = CString(conversion_based_unit->Name().c_str());
      georefdata.MapUnitToMeters = GetConversionFactor<Schema>(conversion_based_unit);
      return;
   }

   WBFL::System::Logger::Info(_T("IfcProjectedCRS.MapUnit is not a recognized unit type. Using default map unit."));
}

// Builds an IfcNamedUnit for IfcProjectedCRS.MapUnit from the pragmatic fields above. Mirrors, in
// reverse, ImportMapUnit's simplifications (SI branch is always plain metre, no prefix reconstruction).
template <typename Schema>
typename Schema::IfcNamedUnit* CreateMapUnit(const GeoreferencingData& georefdata)
{
   USES_CONVERSION;

   if (georefdata.IsMapUnitSI)
   {
      return new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
   }

   auto dims = new typename Schema::IfcDimensionalExponents(1, 0, 0, 0, 0, 0, 0); // length dimension
   auto si_meter = new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
   auto factor_value = new typename Schema::IfcLengthMeasure(georefdata.MapUnitToMeters);
   auto measure_with_unit = new typename Schema::IfcMeasureWithUnit(factor_value, si_meter);
   return new typename Schema::IfcConversionBasedUnit(dims, Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, std::string(T2A(georefdata.MapUnitName)), measure_with_unit);
}

// {03917488-9929-41F1-AB85-8FED05E73009}
DEFINE_GUID(IID_IGeoreferencing,
   0x3917488, 0x9929, 0x41f1, 0xab, 0x85, 0x8f, 0xed, 0x5, 0xe7, 0x30, 0x9);
class IGeoreferencing
{
public:
   virtual void SetGeoreferencingData(const GeoreferencingData& data) = 0;
   virtual const GeoreferencingData& GetGeoreferencingData() const = 0;
};
