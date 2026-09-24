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
#include <IFace\Project.h>
#include "IfcImportUnits.h"

class CIfcImportOptions
{
public:
   enum class ModelElements
   {
      AlignmentOnly,
      AlignmentAndBridge
   };

   ModelElements model_elements = ModelElements::AlignmentAndBridge;

   // When false, the import runs without user interaction (e.g. from the command line).
   // No dialogs are shown, choices default to the first option, and the import
   // log is written to log_file instead of being displayed.
   bool interactive = true;
   CString log_file;
};


///////////////////////////////////////////////////////////////////////////
// CIfcImporter
//
// Converts data between IFC and PGSuper data structures
class CIfcImporter
{
public:
   enum class ImportResult
   {
      Success,
      Fail,
      NotFound
   };

   CIfcImporter(std::shared_ptr<WBFL::EAF::Broker> pBroker);
   ~CIfcImporter(void);

   // Converts Ifc data to PGSuper data
   HRESULT ImportFromIFC(CString& strFilePath, CIfcImportOptions options);

   static Float64 GetPrecision() { return m_Precision; }

   std::shared_ptr<WBFL::EAF::Broker> GetBroker() { return m_pBroker; }

   bool IsInteractive() const { return m_Options.interactive; }

   // Project length and angle units, for converting geometric values to system units
   const WBFL::Units::Length& GetLengthUnit() const { return m_LengthUnit; }
   const WBFL::Units::Angle& GetAngleUnit() const { return m_AngleUnit; }

   // Converts property values to system units. See CIfcImportUnits
   static const CIfcImportUnits& GetUnits() { return m_Units; }


private:
   std::shared_ptr<WBFL::EAF::Broker> m_pBroker;
   CIfcImportOptions m_Options;
   static Float64 m_Precision;
   static CIfcImportUnits m_Units;
   WBFL::Units::Length m_LengthUnit{ 1.0, _T("m") };
   WBFL::Units::Angle m_AngleUnit{ 1.0, _T("rad") };
   std::ostringstream m_LogStream;
   std::ostream* m_pOldLogStream = nullptr;

   ImportResult ImportGeoreferencing(ifcopenshell::file& file);
   ImportResult ImportAlignment(ifcopenshell::file& file);
   ImportResult ImportBridge(ifcopenshell::file& file, bool bDeriveAlignmentFromDeck);

   void InitUnits(ifcopenshell::file& file);
};

