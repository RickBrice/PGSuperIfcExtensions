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

#include <IFace/Tools.h>
#include <IFace\Project.h>

class CIfcImportOptions
{
public:
   enum class ModelElements
   {
      AlignmentOnly,
      AlignmentAndBridge
   };

   ModelElements model_elements = ModelElements::AlignmentAndBridge;
};


///////////////////////////////////////////////////////////////////////////
// CIfcImporter
//
// Converts data between IFC and PGSuper data structures
class CIfcImporter
{
public:
   CIfcImporter(std::shared_ptr<WBFL::EAF::Broker> pBroker);
   ~CIfcImporter(void);

   // Converts Ifc data to PGSuper data
   HRESULT ImportFromIFC(CString& strFilePath, CIfcImportOptions options);

   void AddNote(const std::_tstring& str) { m_Notes.push_back(str); }
   void AddNote(LPCTSTR str) { m_Notes.push_back(str); }

   // Returns a list of notes that were generated during the IFC to PGSuper conversion process
   std::vector<std::_tstring> GetNotes();

   static Float64 GetPrecision() { return m_Precision; }

   std::shared_ptr<WBFL::EAF::Broker> GetBroker() { return m_pBroker; }

   const WBFL::Units::Length& GetLengthUnit() const { return *m_pLengthUnit; }
   const WBFL::Units::Angle& GetAngleUnit() const { return *m_pAngleUnit; }


private:
   std::shared_ptr<WBFL::EAF::Broker> m_pBroker;
   static Float64 m_Precision;
   const WBFL::Units::Length* m_pLengthUnit;
   const WBFL::Units::Angle* m_pAngleUnit;
   std::vector<std::_tstring> m_Notes;

   bool ImportAlignment(IfcParse::IfcFile& file);
   bool ImportBridge(IfcParse::IfcFile& file);

   void InitUnits(IfcParse::IfcFile& file);

   void SetGirderProperties(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc);
   void ImportSlab(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc);
};

