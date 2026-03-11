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
#include "IfcImporter.h"
#include "IfcImporterException.h"
#include "IfcAlignmentImporter.h"
#include "IfcBridgeImporter.h"
#include "Units.h"
#include "ImportResults.h"

#include <EAF/AutoProgress.h>

std::string& to_lower(std::string& s)
{
   std::transform(s.begin(), s.end(), s.begin(), [](auto c) {return std::tolower(c); });
   return s;
}

Float64 CIfcImporter::m_Precision = 0.001;

CIfcImporter::CIfcImporter(std::shared_ptr<WBFL::EAF::Broker> pBroker) :
   m_pBroker(pBroker)
{
   m_pLengthUnit = nullptr;
   m_pAngleUnit = nullptr;
}

CIfcImporter::~CIfcImporter(void)
{
}

void CIfcImporter::InitUnits(IfcParse::IfcFile& file)
{
   auto geometric_representation_contexts = file.instances_by_type<IfcSchema::IfcGeometricRepresentationContext>();
   auto geometric_representation_context = (0 < geometric_representation_contexts->size()) ? *(geometric_representation_contexts->begin()) : nullptr;
#pragma Reminder("WORKING HERE - There could be multiple geometric representation contexts, how do we know if we have the right one?")
   if (geometric_representation_context && geometric_representation_context->Precision() != boost::none)
   {
      m_Precision = *(geometric_representation_context->Precision());
   }

#pragma Reminder("WORKING HERE - UNITS - THERE ARE MANY CASES THIS DOESN'T DEAL WITH")
   auto unit_assignment_instances = file.instances_by_type<IfcSchema::IfcUnitAssignment>();
   ATLASSERT(unit_assignment_instances->size() == 1);
   auto unit_assignment = *(unit_assignment_instances->begin());
   auto units = unit_assignment->Units();
   for (auto unit : *units)
   {
      auto derived_unit = unit->as<IfcSchema::IfcDerivedUnit>();
      auto monitary_unit = unit->as<IfcSchema::IfcMonetaryUnit>();
      auto si_unit = unit->as<IfcSchema::IfcSIUnit>();
      auto conversion_based_unit = unit->as<IfcSchema::IfcConversionBasedUnit>();
      auto conversion_based_unit_with_offset = unit->as<IfcSchema::IfcConversionBasedUnitWithOffset>();

      if (si_unit)
      {
         if (si_unit->Name() == IfcSchema::IfcSIUnitName::IfcSIUnitName_METRE)
         {
            if (si_unit->Prefix() != boost::none)
            {
               switch (*(si_unit->Prefix()))
               {
               case IfcSchema::IfcSIPrefix::IfcSIPrefix_KILO:
                  m_pLengthUnit = &WBFL::Units::Measure::Kilometer;
                  break;

               case IfcSchema::IfcSIPrefix::IfcSIPrefix_CENTI:
                  m_pLengthUnit = &WBFL::Units::Measure::Centimeter;
                  break;

               case IfcSchema::IfcSIPrefix::IfcSIPrefix_MILLI:
                  m_pLengthUnit = &WBFL::Units::Measure::Millimeter;
                  break;

               default:
                  ATLASSERT(false); // unit prefix isn't supported
               }
            }
            else
            {
               m_pLengthUnit = &WBFL::Units::Measure::Meter;
            }
            continue;
         }

         if (si_unit->Name() == IfcSchema::IfcSIUnitName::IfcSIUnitName_RADIAN)
         {
            ATLASSERT(si_unit->Prefix() == boost::none); // not expecting anything like Kilo-radians
            m_pAngleUnit = &WBFL::Units::Measure::Radian;
            continue;
         }
      }

      if (conversion_based_unit)
      {
         auto conversion_factor = GetConversionFactor<IfcSchema>(conversion_based_unit);

         auto measure_with_unit = conversion_based_unit->ConversionFactor();
         auto unit_component = measure_with_unit->UnitComponent()->as<IfcSchema::IfcSIUnit>();

         if (conversion_based_unit->UnitType() == IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT)
         {
            ATLASSERT(unit_component->Name() == IfcSchema::IfcSIUnitName::IfcSIUnitName_RADIAN);

            if (IsEqual(conversion_factor, WBFL::Units::Measure::Degree.GetConvFactor()))
            {
               m_pAngleUnit = &WBFL::Units::Measure::Degree;
            }
         }
         else if (conversion_based_unit->UnitType() == IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT)
         {
            if (unit_component->Prefix() == IfcSchema::IfcSIPrefix::IfcSIPrefix_MILLI)
            {
               // lengths are in millimeter, so divide the conversion factor by 1000.
               // so it is in meter so we can match the WBFL::Measure::Length conversion factors, which convert to/from meter
               conversion_factor /= 1000.0;
            }

            ATLASSERT(unit_component->Name() == IfcSchema::IfcSIUnitName::IfcSIUnitName_METRE);

            std::string name = conversion_based_unit->Name();
            to_lower(name);
            if (IsEqual(conversion_factor, WBFL::Units::Measure::Feet.GetConvFactor()) || name == std::string("foot"))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Feet;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::USSurveyFoot.GetConvFactor()) || name == std::string("us survey foot"))
            {
               m_pLengthUnit = &WBFL::Units::Measure::USSurveyFoot;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::Inch.GetConvFactor()) || name == std::string("inch"))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Inch;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::Mile.GetConvFactor()) || name == std::string("mile"))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Mile;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::Yard.GetConvFactor()) || name == std::string("yard"))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Yard;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::USSurveyYard.GetConvFactor()))
            {
               m_pLengthUnit = &WBFL::Units::Measure::USSurveyYard;
            }
            else
            {
               ATLASSERT(false); // we don't have a unit of measure for this
            }
         }
         continue;
      }
   }

   // if the file doesn't have unit information,
   // assume meter and radian
   if (m_pLengthUnit == nullptr)
      m_pLengthUnit = &WBFL::Units::Measure::Meter;

   if(m_pAngleUnit == nullptr)
      m_pAngleUnit = &WBFL::Units::Measure::Radian;

}

// Here is a good reference to redirecting cout
// https://stackoverflow.com/questions/4810516/c-redirecting-stdout
#include <mutex>
class ProgressStringBuf : public std::stringbuf
{
public:
    ProgressStringBuf() : _accum(""), _lineNum(0), m_pProgress(nullptr) {};
    void SetProgress(std::shared_ptr<IEAFProgress> pProgress) { m_pProgress = pProgress; }
protected:
    virtual std::streamsize xsputn(const char* s, std::streamsize num)
    {
        std::mutex m;
        std::lock_guard<std::mutex> lg(m);

        //// Prepend with the line number
        std::string str(s, (const uint32_t)num);
        //str = std::to_string(_lineNum) + ": " + str + "\r\n";

        //// Accumulate the latest text to the front
        //_accum = str + _accum;
        _accum = str;

        // Write to the Win32 dialog edit control.
        m_pProgress->UpdateMessage((LPCTSTR)(std::_tstring(_accum.begin(), _accum.end())).c_str());

        _lineNum++;
        return(num);
    }

private:
    std::string _accum;
    uint32_t _lineNum;
    std::shared_ptr<IEAFProgress> m_pProgress;
};

class ProgressStream : public std::ostream
{
public:
    ProgressStream() : std::ostream(&_progress) {};
    void SetProgress(std::shared_ptr<IEAFProgress> pProgress) { _progress.SetProgress(pProgress); }
private:
    ProgressStringBuf _progress;
};

HRESULT CIfcImporter::ImportFromIFC(CString& strFilePath, CIfcImportOptions options)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   USES_CONVERSION;

   HRESULT hr = S_OK;
   try
   {
      m_pOldLogStream = WBFL::System::Logger::SetOutput(&m_LogStream);

      std::unique_ptr<IfcParse::IfcFile> pFile = nullptr;

      GET_IFACE(IEAFProgress, pProgress);
      WBFL::EAF::AutoProgress ap(pProgress);

      auto del = [&](std::streambuf* p) {std::cout.rdbuf(p); };
      std::unique_ptr<std::streambuf, decltype(del)> origBuffer(std::cout.rdbuf(), del);
      ProgressStream p;
      p.SetProgress(pProgress);

      p.copyfmt(std::cout);
      std::cout.rdbuf(p.rdbuf());

      Logger::SetOutput(&std::cout, &std::cout);

      pFile = std::make_unique<IfcParse::IfcFile>(T2A(strFilePath.GetBuffer()));

      if (!pFile->good())
      {
         IFC_THROW(_T("Unable to parse .ifc file"));
      }

      auto strSchemaName = pFile->schema()->name();
      if (strSchemaName == std::string("IFC4X3_ADD2"))
      {
         GET_IFACE(IEvents, pEvents);
         pEvents->HoldEvents();

         InitUnits(*pFile);

         auto import_alignment_result = ImportAlignment(*pFile);
         if (import_alignment_result == ImportResult::Fail)
            hr = E_FAIL;

         pEvents->FirePendingEvents(); // update internal data for correct alignment

         if (options.model_elements == CIfcImportOptions::ModelElements::AlignmentAndBridge)
         {
            auto import_result = ImportBridge(*pFile,import_alignment_result == ImportResult::NotFound);
            if (import_result == ImportResult::Fail || import_result == ImportResult::NotFound)
               hr = E_FAIL;
         }

         pEvents->FirePendingEvents();
      }
      else
      {
         IFC_THROW(_T("Schema not supported"));
      }
   }
    catch (CIfcImporterException& e)
    {
       std::_tostringstream os;
       os << _T("IFC import failed:\n") << e.What();
       WBFL::System::Logger::Info(os.str().c_str());
       hr = E_FAIL;
    }

    CImportResults dlg(m_LogStream);
    dlg.DoModal();

    WBFL::System::Logger::SetOutput(m_pOldLogStream);

   return hr;
}

CIfcImporter::ImportResult CIfcImporter::ImportAlignment(IfcParse::IfcFile& file)
{
   return CIfcAlignmentImporter(*this).Import(file);
}

CIfcImporter::ImportResult CIfcImporter::ImportBridge(IfcParse::IfcFile& file,bool bDeriveAlignmentFromDeck)
{
   return CIfcBridgeImporter(*this).Import(file, bDeriveAlignmentFromDeck);
}
