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
#include "IfcGeoreferencingImporter.h"
#include "IfcAlignmentImporter.h"
#include "IfcBridgeImporter.h"
#include "Units.h"
#include "ImportResults.h"

#include <EAF/AutoProgress.h>

Float64 CIfcImporter::m_Precision = 0.001;
CIfcImportUnits CIfcImporter::m_Units;

CIfcImporter::CIfcImporter(std::shared_ptr<WBFL::EAF::Broker> pBroker) :
   m_pBroker(pBroker)
{
}

CIfcImporter::~CIfcImporter(void)
{
}

void CIfcImporter::InitUnits(ifcopenshell::file& file)
{
   WBFL::System::Logger::Info(_T("Initializing units and precision from IFC file."));
   auto geometric_representation_contexts = file.instances_by_type<IfcSchema::IfcGeometricRepresentationContext>();
   auto geometric_representation_context = (0 < geometric_representation_contexts.size()) ? geometric_representation_contexts.front() : IfcSchema::IfcGeometricRepresentationContext{};
#pragma Reminder("WORKING HERE - There could be multiple geometric representation contexts, how do we know if we have the right one?")
   if (geometric_representation_context && geometric_representation_context.Precision() != std::nullopt)
   {
      m_Precision = *(geometric_representation_context.Precision());
   }

   m_Units.Init(file);

   // Geometry is always in project units
   USES_CONVERSION;
   m_LengthUnit = WBFL::Units::Length(m_Units.GetConversionFactor({}, IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT), A2T(m_Units.GetProjectUnitName(IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT).c_str()));
   m_AngleUnit = WBFL::Units::Angle(m_Units.GetConversionFactor({}, IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT), A2T(m_Units.GetProjectUnitName(IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT).c_str()));
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


   m_Options = options;

   HRESULT hr = S_OK;
   try
   {
      m_pOldLogStream = WBFL::System::Logger::SetOutput(&m_LogStream);

      WBFL::System::Logger::Info(_T("Starting IFC import from file"));

      std::unique_ptr<ifcopenshell::file> pFile = nullptr;

      GET_IFACE(IEAFProgress, pProgress);
      WBFL::EAF::AutoProgress ap(pProgress);

      auto del = [&](std::streambuf* p) {std::cout.rdbuf(p); };
      std::unique_ptr<std::streambuf, decltype(del)> origBuffer(std::cout.rdbuf(), del);
      ProgressStream p;
      p.SetProgress(pProgress);

      p.copyfmt(std::cout);
      std::cout.rdbuf(p.rdbuf());

      ifcopenshell::logger::root().set_output(&std::cout, &std::cout);

      pFile = std::make_unique<ifcopenshell::file>(T2A(strFilePath.GetBuffer()));

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

         auto import_georef_result = ImportGeoreferencing(*pFile);
         if (import_georef_result == ImportResult::Fail)
            hr = E_FAIL;

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
    catch (const std::exception& e)
    {
       // IfcOpenShell reports errors by throwing std::exception (e.g. the geometry iterator)
       // Don't let them escape - that terminates the application
       std::ostringstream os;
       os << "IFC import failed:\n" << e.what();
       WBFL::System::Logger::Info(os.str().c_str());
       hr = E_FAIL;
    }
    catch (...)
    {
       WBFL::System::Logger::Info(_T("IFC import failed:\nUnknown exception"));
       hr = E_FAIL;
    }

    WBFL::System::Logger::Info(_T("Done IFC import from file"));

    if (m_Options.interactive)
    {
       CImportResults dlg(m_LogStream);
       dlg.DoModal();
    }
    else if (!m_Options.log_file.IsEmpty())
    {
       std::ofstream log(m_Options.log_file.GetString());
       log << m_LogStream.str();
    }

    WBFL::System::Logger::SetOutput(m_pOldLogStream);

   return hr;
}

CIfcImporter::ImportResult CIfcImporter::ImportGeoreferencing(ifcopenshell::file& file)
{
   return CIfcGeoreferencingImporter(*this).Import(file);
}

CIfcImporter::ImportResult CIfcImporter::ImportAlignment(ifcopenshell::file& file)
{
   return CIfcAlignmentImporter(*this).Import(file);
}

CIfcImporter::ImportResult CIfcImporter::ImportBridge(ifcopenshell::file& file,bool bDeriveAlignmentFromDeck)
{
   return CIfcBridgeImporter(*this).Import(file, bDeriveAlignmentFromDeck);
}
