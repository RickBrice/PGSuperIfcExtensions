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
#include "IfcPropertyWriter.h"

CIfcExportSession* CIfcExportSession::ms_pCurrent = nullptr;

CIfcExportSession::CIfcExportSession(const CIfcExportOptions& options)
{
   ASSERT(ms_pCurrent == nullptr); // one export at a time
   std::filesystem::path path(options.mapping_file.GetString());
   m_pTable = CIfcMappingTable::Load(path, path.empty() ? MappingTableSource::InstalledStandard : MappingTableSource::CommandLine);
   ms_pCurrent = this;
}

CIfcExportSession::~CIfcExportSession()
{
   ms_pCurrent = nullptr;
}

CIfcExportSession& CIfcExportSession::Current()
{
   ASSERT(ms_pCurrent); // only available during an export
   return *ms_pCurrent;
}

bool IncludePropertySet(const PropertySetDeclaration& pset, const CIfcExportOptions& options)
{
   switch (pset.condition)
   {
   case PropertySetDeclaration::Condition::Classify: return options.classify;
   case PropertySetDeclaration::Condition::Quantities: return options.include_quantities;
   case PropertySetDeclaration::Condition::Always: return true;
   }
   ASSERT(false);
   return true;
}

bool UseDisplayUnits(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   if (!options.display_units_for_properties)
      return false;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   return pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US;
}

Float64 ConvertToDisplayUnits(std::shared_ptr<WBFL::EAF::Broker> pBroker, const TargetDef& target, Float64 value)
{
   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   switch (target.display_unit)
   {
   case ExportUnit::SpanLength: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure); break;
   case ExportUnit::Deflection: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure); break;
   case ExportUnit::Stress: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetStressUnit().UnitOfMeasure); break;
   case ExportUnit::Angle: value = WBFL::Units::ConvertFromSysUnits(value, pDisplayUnits->GetAngleUnit().UnitOfMeasure); break;
   default: break;
   }

   if (0 < target.display_round)
      value = RoundOff(value, target.display_round);

   return value;
}
