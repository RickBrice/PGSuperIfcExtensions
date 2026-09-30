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

/*****************************************************************************
CLASS
   CIfcTargetReader

   Reads target values from IFC elements using the locations of a mapping table,
   and selects the IFC elements that play an element role.

   A value is converted to PGSuper system units (SI) as soon as it's read: with
   its own IFC unit if it has one, else the IFC project unit for its measure,
   else the unit the table states for the location (plain numbers and text).
*****************************************************************************/

#include "IfcMappingTable.h"
#include "IfcImportUnits.h"

struct TargetReading
{
   TargetValue value;
   const MappingLocation* location = nullptr;
   std::string raw; // the value as found in the model

   // e.g. "\"6.8\" from IaDOT_PPCB.6_Concrete Release Strength, Fci (ksi) in mapping table \"Iowa DOT\""
   std::string Source() const;
};

class CIfcTargetReader
{
public:
   CIfcTargetReader(const CIfcMappingTable& table, const CIfcImportUnits& units);

   const CIfcMappingTable& GetTable() const { return m_Table; }

   // The first value found in the target's locations, in location order. Values that are found
   // but can't be used (e.g. text that isn't a number) are logged, and the next location is tried
   std::optional<TargetReading> Read(std::string_view target, IfcSchema::IfcObject object) const;

   // Every value found in the target's locations, in location order, without duplicates
   std::vector<TargetReading> ReadAll(std::string_view target, IfcSchema::IfcObject object) const;

   // Logs that a target wasn't found for an element, the locations that were tried, and hints: properties of the
   // element that look like they might hold the value (see CIfcTargetHints). Hints are never used as values.
   // Only the first element is logged for each target. The others are counted, see LogNotFoundSummary
   void ReportNotFound(std::string_view target, IfcSchema::IfcObject object, const std::string& element_name) const;

   // Logs how many more elements each target wasn't found for
   void LogNotFoundSummary() const;

   // True if a table defines the element role
   bool HasSelector(ElementKind role) const;

   // True if the element plays the role. False if no table defines the role
   bool Matches(ElementKind role, IfcSchema::IfcObject object) const;

   // The elements that play the role, in the order they are in the file
   std::vector<IfcSchema::IfcObject> Select(ElementKind role, ifcopenshell::file& file) const;

private:
   const CIfcMappingTable& m_Table;
   const CIfcImportUnits& m_Units;
   mutable std::map<std::string, size_t, std::less<>> m_NotFound; // target -> number of elements it wasn't found for
   mutable std::set<std::string> m_Problems; // logged problems, so each is logged once

   const TargetDef& GetTarget(std::string_view target) const;
   void Read(const TargetDef& target, IfcSchema::IfcObject object, bool bAll, std::vector<TargetReading>& readings) const;
   bool Matches(const ElementSelector& selector, IfcSchema::IfcObject object) const;
   void LogProblem(const std::string& problem) const;
};
