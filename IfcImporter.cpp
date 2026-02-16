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
#include "stdafx.h"
#include "IfcImporter.h"
#include "IfcImporterException.h"
#include "IfcAlignmentImporter.h"
#include "Properties.h"

#include <MFCTools\Prompts.h>

#include <boost/range/combine.hpp>

#include <IFace/Tools.h>
#include <EAF/AutoProgress.h>

#include <Units\Units.h>

#include <psgLib/BridgeDescription2.h>
#include <psgLib/GirderLabel.h>
#include <psgLib/GirderLibraryEntry.h>

// Constants for tracking the state of converting the profile data
#define PROFILE_NOT_STARTED -1
#define PROFILE_ESTABLISHED  1

Float64 CIfcImporter::m_Precision = 0.001;


std::pair<GroupIndexType, GirderIndexType> ExtractSpanAndGirder(const std::string& s)
{
   GroupIndexType grpIdx = INVALID_INDEX;
   GirderIndexType gdrIdx = INVALID_INDEX;
   std::istringstream iss(s);
   std::string word;

   while (iss >> word)
   {
      if (word == "Span")
         iss >> grpIdx;
      else if (word == "Girder")
         iss >> gdrIdx;
      else if (word == ",")
      { // do nothing
      }
      else
      {
         USES_CONVERSION;
         std::_tostringstream os;
         os << _T("Unexpected DesignLocationNumber property in Pset_PrecastConcreteElementGeneral property set (") << A2T(s.c_str()) << _T(")");
         IFC_THROW(os.str().c_str());
      }
   }

   return { grpIdx - 1,gdrIdx - 1 };
}

Ifc4x3_add2::IfcBridgePart* GetBridgePart(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridgePartTypeEnum part_type)
{
   auto parts = file.instances_by_type<Ifc4x3_add2::IfcBridgePart>();
   for (auto part : *parts)
   {
      if (part->PredefinedType().has_value() && part->PredefinedType().get() == part_type)
      {
         return part;
      }
   }
   return nullptr;
}

std::vector<Ifc4x3_add2::IfcBridgePart*> GetBridgeParts(IfcParse::IfcFile& file, Ifc4x3_add2::IfcBridgePartTypeEnum part_type)
{
   std::vector<Ifc4x3_add2::IfcBridgePart*> parts_found;
   auto parts = file.instances_by_type<Ifc4x3_add2::IfcBridgePart>();
   for (auto part : *parts)
   {
      if (part->PredefinedType().has_value() && part->PredefinedType().get() == part_type)
      {
         parts_found.push_back(part);
      }
   }
   return parts_found;
}

bool HasClassification(Ifc4x3_add2::IfcObjectDefinition* object, std::string identifier)
{
   auto associations = object->HasAssociations();
   if (associations)
   {
      for (auto rel : *associations)
      {
         auto rel_associates_classification = rel->as < Ifc4x3_add2::IfcRelAssociatesClassification>();
         auto classification_reference = rel_associates_classification->RelatingClassification()->as<Ifc4x3_add2::IfcClassificationReference>();
         if (classification_reference && classification_reference->Identification().value_or("") == identifier)
            return true;
      }
   }

   return false;
}

template <typename Schema>
typename Schema::IfcMaterial* GetMaterial(typename Schema::IfcObjectDefinition* objectdef)
{
   auto associations = objectdef->HasAssociations();
   if (associations)
   {
      for (auto rel : *associations)
      {
         auto rel_associates_material = rel->as<Ifc4x3_add2::IfcRelAssociatesMaterial>();
         if (rel_associates_material)
         {
            auto material = rel_associates_material->RelatingMaterial();
            return material->as<typename Schema::IfcMaterial>();
         }
      }
   }

   return nullptr;
}

int GetBeamTypeCount(IfcParse::IfcFile& file)
{
   int count = 0;
   auto beam_types = file.instances_by_type<Ifc4x3_add2::IfcBeamType>();
   for (auto beam_type : *beam_types)
   {
      if (beam_type->PredefinedType() == Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_BEAM)
         count++;
   }

   return count;
}

template <typename E>
E* GetType(Ifc4x3_add2::IfcObject* object)
{
   auto types = object->IsTypedBy();
   for (auto type : *types)
   {
      return type->RelatingType()->as<E>();
   }

   return nullptr;
}

template <typename T, typename E>
E GetPredefinedType(Ifc4x3_add2::IfcObject* object)
{
   auto types = object->IsTypedBy();
   for (auto type : *types)
   {
      return type->RelatingType()->as<T>()->PredefinedType();
   }

   return object->as<T>()->PredefinedType();
}

///////////////////////////////////////////////////////
///////////////////////////////////////////////////////
///////////////////////////////////////////////////////

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
   auto geometric_representation_contexts = file.instances_by_type<Ifc4x3_add2::IfcGeometricRepresentationContext>();
   auto geometric_representation_context = (0 < geometric_representation_contexts->size()) ? *(geometric_representation_contexts->begin()) : nullptr;
#pragma Reminder("WORKING HERE - There could be multiple geometric representation contexts, how do we know if we have the right one?")
   if (geometric_representation_context && geometric_representation_context->Precision() != boost::none)
   {
      m_Precision = *(geometric_representation_context->Precision());
   }

#pragma Reminder("WORKING HERE - UNITS - THERE ARE MANY CASES THIS DOESN'T DEAL WITH")
   auto unit_assignment_instances = file.instances_by_type<Ifc4x3_add2::IfcUnitAssignment>();
   ATLASSERT(unit_assignment_instances->size() == 1);
   auto unit_assignment = *(unit_assignment_instances->begin());
   auto units = unit_assignment->Units();
   for (auto unit : *units)
   {
      auto derived_unit = unit->as<Ifc4x3_add2::IfcDerivedUnit>();
      auto monitary_unit = unit->as<Ifc4x3_add2::IfcMonetaryUnit>();
      auto si_unit = unit->as<Ifc4x3_add2::IfcSIUnit>();
      auto conversion_based_unit = unit->as<Ifc4x3_add2::IfcConversionBasedUnit>();
      auto conversion_based_unit_with_offset = unit->as<Ifc4x3_add2::IfcConversionBasedUnitWithOffset>();

      if (si_unit)
      {
         if (si_unit->Name() == Ifc4x3_add2::IfcSIUnitName::IfcSIUnitName_METRE)
         {
            if (si_unit->Prefix() != boost::none)
            {
               switch (*(si_unit->Prefix()))
               {
               case Ifc4x3_add2::IfcSIPrefix::IfcSIPrefix_KILO:
                  m_pLengthUnit = &WBFL::Units::Measure::Kilometer;
                  break;

               case Ifc4x3_add2::IfcSIPrefix::IfcSIPrefix_CENTI:
                  m_pLengthUnit = &WBFL::Units::Measure::Centimeter;
                  break;

               case Ifc4x3_add2::IfcSIPrefix::IfcSIPrefix_MILLI:
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

         if (si_unit->Name() == Ifc4x3_add2::IfcSIUnitName::IfcSIUnitName_RADIAN)
         {
            ATLASSERT(si_unit->Prefix() == boost::none); // not expecting anything like Kilo-radians
            m_pAngleUnit = &WBFL::Units::Measure::Radian;
            continue;
         }
      }

      if (conversion_based_unit)
      {
         auto measure_with_unit = conversion_based_unit->ConversionFactor();
         auto unit_component = measure_with_unit->UnitComponent()->as<Ifc4x3_add2::IfcSIUnit>();
         Float64 conversion_factor;
         try
         {
            auto value_component = measure_with_unit->ValueComponent();
            ATLASSERT(value_component); // not dealing with anything but simple conversion factors
            conversion_factor = (Float64)(*value_component->as<Ifc4x3_add2::IfcReal>());
            //conversion_factor = (Float64)(value_component->data().get_attribute_value(0));
         }
         catch (IfcParse::IfcInvalidTokenException& e)
         {
            // Was expecting something like 
            // #15 = IFCMEASUREWITHUNIT(IFCLENGTHMEASURE(3.28083333333333), #16);
            // where the expected token is IFCLENGTHMEASURE, but instead found something like
            // #15=IFCMEASUREWITHUNIT(3.28083333333333,#16);
            // we'll just get the value and keep going
            TRACE(e.what());
            auto pArgument = measure_with_unit->get("ValueComponent");
            ATLASSERT(pArgument.type() == IfcUtil::Argument_DOUBLE);
            conversion_factor = double(pArgument);
         }

         if (unit_component->Prefix() == Ifc4x3_add2::IfcSIPrefix::IfcSIPrefix_MILLI)
         {
            // lengths are in millimeter, so divide the conversion factor by 1000.
            // so it is in meter so we can match the WBFL::Measure::Length conversion factors, which convert to/from meter
            conversion_factor /= 1000.0;
         }

         if (conversion_based_unit->UnitType() == Ifc4x3_add2::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT)
         {
            ATLASSERT(unit_component->Name() == Ifc4x3_add2::IfcSIUnitName::IfcSIUnitName_RADIAN);

            if (IsEqual(conversion_factor, WBFL::Units::Measure::Degree.GetConvFactor()))
            {
               m_pAngleUnit = &WBFL::Units::Measure::Degree;
            }
         }
         else if (conversion_based_unit->UnitType() == Ifc4x3_add2::IfcUnitEnum::IfcUnit_LENGTHUNIT)
         {
            ATLASSERT(unit_component->Name() == Ifc4x3_add2::IfcSIUnitName::IfcSIUnitName_METRE);

            if (IsEqual(conversion_factor, WBFL::Units::Measure::Feet.GetConvFactor()))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Feet;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::USSurveyFoot.GetConvFactor()))
            {
               m_pLengthUnit = &WBFL::Units::Measure::USSurveyFoot;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::Inch.GetConvFactor()))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Inch;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::Mile.GetConvFactor()))
            {
               m_pLengthUnit = &WBFL::Units::Measure::Mile;
            }
            else if (IsEqual(conversion_factor, WBFL::Units::Measure::Yard.GetConvFactor()))
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
    USES_CONVERSION;

    try
    {
       std::unique_ptr<IfcParse::IfcFile> pFile = nullptr;

       { // scope the progress window so it closes automatically when we are done with it
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
       }

       m_Notes.clear();

       auto strSchemaName = pFile->schema()->name();
       if (strSchemaName == std::string("IFC4X3_ADD2"))
       {
          GET_IFACE(IEvents, pEvents);
          pEvents->HoldEvents();

          InitUnits(*pFile);

          if (options.model_elements == CIfcImportOptions::ModelElements::AlignmentOnly)
          {
             ImportAlignment(*pFile);
          }
          else
          {
             if (ImportAlignment(*pFile))
                ImportBridge(*pFile);
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
       AfxMessageBox(os.str().c_str());
       return E_FAIL;
    }

   auto notes = GetNotes();
   std::_tstring strNotes;
   for (auto note : notes)
   {
      strNotes += note + _T("\n\n");
   }
   if (0 < strNotes.size())
   {
      AfxMessageBox(strNotes.c_str(), MB_OK);
   }

   return S_OK;
}

bool CIfcImporter::ImportAlignment(IfcParse::IfcFile& file)
{
   return CIfcAlignmentImporter(*this).Import(file);
}

bool CIfcImporter::ImportBridge(IfcParse::IfcFile& file)
{
   auto bridges = file.instances_by_type<Ifc4x3_add2::IfcBridge>();
   auto bridge = (*bridges->begin());

   SpanIndexType nSpans = INVALID_INDEX;
   auto value = GetProperty<Ifc4x3_add2,Ifc4x3_add2::IfcInteger>(bridge, "usBridge_BridgeCommon", "usBridge_NumberOfSpans");
   if (value)
   {
      nSpans = (SpanIndexType)(*value);
   }
   else
   {
      IFC_THROW(_T("usBridge_NumberOfSpans property not found in usBridge_BridgeCommon property set"));
   }

   std::vector<GirderIndexType> nGirders;
   nGirders.assign(nSpans, 0);
   auto superstructure = GetBridgePart(file, Ifc4x3_add2::IfcBridgePartTypeEnum::Value::IfcBridgePartType_SUPERSTRUCTURE);
   auto rel_contained_elements = superstructure->ContainsElements();
   for (auto contained_element : *rel_contained_elements)
   {
      auto related_elements = contained_element->RelatedElements();
      for (auto related_element : *related_elements)
      {
         auto beam = related_element->as<Ifc4x3_add2::IfcBeam>();
         if (beam)
         {
            auto value = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLabel>(beam, "Pset_PrecastConcreteElementGeneral", "DesignLocationNumber");
            if (!value)
            {
               IFC_THROW(_T("DesignLocationNumber property not found in Pset_PrecastConcreteElementGeneral"));
            }
            auto [spanIdx, gdrIdx] = ExtractSpanAndGirder(*value);
            nGirders[spanIdx]++;
         }
      }
   }

   bool bSameNumGirdersInAllSpans = std::adjacent_find(nGirders.begin(), nGirders.end(), std::not_equal_to<>()) == nGirders.end() ? true : false;

   // Get the existing bridge description. We are going to modify the bridge description with
   // information extracted from the IFC file.
   GET_IFACE(IBridgeDescription, pIBridgeDesc);
   auto bridge_desc = *(pIBridgeDesc->GetBridgeDescription());

   SpanIndexType nSpansToAdd = nSpans - pIBridgeDesc->GetSpanCount();
   if (0 < nSpansToAdd)
   {
      for (SpanIndexType i = 0; i < nSpansToAdd; i++)
      {
         bridge_desc.AppendSpan(nullptr, nullptr, true, 0);
      }

      bridge_desc.UseSameGirderForEntireBridge(true);

      if (bSameNumGirdersInAllSpans)
      {
         bridge_desc.UseSameNumberOfGirdersInAllGroups(true);
         bridge_desc.SetGirderCount(nGirders.front());
      }
   }

   if (!bSameNumGirdersInAllSpans)
   {
      bridge_desc.UseSameNumberOfGirdersInAllGroups(false);
      for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
      {
         bridge_desc.GetGirderGroup(spanIdx)->SetGirderCount(nGirders[spanIdx]);
      }
   }

   //
   // Position the abutments and piers
   //

   // Per TPF modeling guide, piers and abutments are different types.
   // Get the abutments and piers and put into a single vector because we need to treat them the same in PGSuper
   std::vector<Ifc4x3_add2::IfcBridgePart*> abutments = GetBridgeParts(file, Ifc4x3_add2::IfcBridgePartTypeEnum::Value::IfcBridgePartType_ABUTMENT);
   std::vector<Ifc4x3_add2::IfcBridgePart*> piers = GetBridgeParts(file, Ifc4x3_add2::IfcBridgePartTypeEnum::Value::IfcBridgePartType_PIER);
   piers.insert(piers.begin(), abutments.front());
   piers.insert(piers.end(), abutments.back());

   PierIndexType nPiers = bridge_desc.GetPierCount();
   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      // get the pier station from the positioning element and set it on the PGSuper pier
      auto pPier = bridge_desc.GetPier(pierIdx);
      auto pier = piers[pierIdx];
      auto rel_positions = pier->PositionedRelativeTo();
      if (rel_positions->size() == 0)
      {
         std::_tostringstream os;
         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": must be positioned with an IfcReferent");
         IFC_THROW(os.str().c_str());
      }

      auto positioning_element = (*rel_positions->begin())->RelatingPositioningElement();
      auto ref = positioning_element->as<Ifc4x3_add2::IfcReferent>();
      if (!ref)
      {
         std::_tostringstream os;
         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": must be positioned with an IfcReferent");
         IFC_THROW(os.str().c_str());
      }

      auto station = GetProperty<Ifc4x3_add2,Ifc4x3_add2::IfcLengthMeasure>(ref, "Pset_Stationing", "Station");
      if (!station)
      {
         std::_tostringstream os;
         os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Station property in Pset_Stationing not found");
         IFC_THROW(os.str().c_str());
      }

      pPier->SetStation(*station);
   }

   // NOTE: This is not the cleanest way to do this, but it gets the job done for now.
   // I want to keep girder spacing from a custom property set separate from the pier stationing.
   bridge_desc.SetGirderSpacingType(pgsTypes::SupportedBeamSpacing::sbsGeneral);
   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      auto pier = piers[pierIdx];
      auto pPier = bridge_desc.GetPier(pierIdx);
      if (0 < pierIdx)
      {
         auto spacing = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(pier, "pgsSpacing", "Back_Spacing");
         if (spacing.empty())
         {
            std::_tostringstream os;
            os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Back_Spacing property in pgsSpacing property set not found");
            IFC_THROW(os.str().c_str());
         }

         auto girder_spacing = pPier->GetGirderSpacing(pgsTypes::Back);
         girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem); // this is how the spacing is defined in the exporter
         girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtCenterlineBearing);
         girder_spacing->ExpandAll();
         IndexType idx = 0;
         for (auto s : spacing)
         {
            girder_spacing->SetGirderSpacing(idx++, *s);
         }
      }

      if (pierIdx < nPiers - 1)
      {
         auto spacing = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(pier, "pgsSpacing", "Ahead_Spacing");
         if (spacing.empty())
         {
            std::_tostringstream os;
            os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": Ahead_Spacing property in pgsSpacing property set not found");
            IFC_THROW(os.str().c_str());
         }

         auto girder_spacing = pPier->GetGirderSpacing(pgsTypes::Ahead);
         girder_spacing->SetMeasurementType(pgsTypes::MeasurementType::NormalToItem); // this is how the spacing is defined in the exporter
         girder_spacing->SetMeasurementLocation(pgsTypes::MeasurementLocation::AtCenterlineBearing);

         girder_spacing->ExpandAll();
         IndexType idx = 0;
         for (auto s : spacing)
         {
            girder_spacing->SetGirderSpacing(idx++, *s);
         }
      }
   }

   SetGirderProperties(file, bridge_desc);

   ImportSlab(file, bridge_desc);

   pIBridgeDesc->SetBridgeDescription(bridge_desc);

   return true;
}

std::vector<std::_tstring> CIfcImporter::GetNotes()
{
   return m_Notes;
}


void CIfcImporter::SetGirderProperties(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   auto beams = file.instances_by_type<Ifc4x3_add2::IfcBeam>();
   aggregate_of<Ifc4x3_add2::IfcBeam>::ptr prestressed_beams(new aggregate_of<Ifc4x3_add2::IfcBeam>());
   for (auto beam : *beams)
   {
      auto predefined_type = GetPredefinedType<Ifc4x3_add2::IfcBeamType, Ifc4x3_add2::IfcBeamTypeEnum>(beam);
      if (predefined_type == Ifc4x3_add2::IfcBeamTypeEnum::IfcBeamType_BEAM && HasClassification(beam,"usBridge_GirderPrestressedConcrete"))
         prestressed_beams->push(beam);
   }

   int beam_type_count = GetBeamTypeCount(file);

   bridge_desc.UseSameGirderForEntireBridge(beam_type_count == 1 ? true : false);
   GET_IFACE(ILibrary, pLibrary);
   if (bridge_desc.UseSameGirderForEntireBridge())
   {
      auto type = GetType<Ifc4x3_add2::IfcBeamType>(*beams->begin());
      auto girder_name = type->Name().get_value_or(std::string("Unknown"));
      auto girder_library_entry = pLibrary->GetGirderEntry(A2T(girder_name.c_str()));
      if (!girder_library_entry)
      {
         std::ostringstream os;
         os << "Girder type \"" << girder_name << "\" not found in the library";
         IFC_THROW(A2T(os.str().c_str()));
      }
      bridge_desc.SetGirderLibraryEntry(girder_library_entry);
      bridge_desc.SetGirderFamilyName(girder_library_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
      auto factory = girder_library_entry->GetBeamFactory();
      auto orientations = factory->GetSupportedGirderOrientation();
      bridge_desc.SetGirderOrientation(orientations.front());
   }

   for (auto beam : *prestressed_beams)
   {
      auto value = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLabel>(beam, "Pset_PrecastConcreteElementGeneral", "DesignLocationNumber");
      if (!value)
      {
         IFC_THROW(_T("DesignLocationNumber in Pset_PrecastConcreteElement property set not found"));
      }

      auto [spanIdx, gdrIdx] = ExtractSpanAndGirder(*value);
      auto fci = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcPressureMeasure>(beam, "Pset_PrecastConcreteElementGeneral", "ReleaseStrength");
      if (!fci)
      {
         IFC_THROW(_T("ReleaseStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
      }
      bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->GetSegment(0)->Material.Concrete.Fci = *fci;

      auto material = GetMaterial<Ifc4x3_add2>(beam);
      auto fc = GetMaterialProperty<Ifc4x3_add2, Ifc4x3_add2::IfcPressureMeasure>(material, "Pset_MaterialConcrete", "CompressiveStrength");
      if (!fc)
      {
         IFC_THROW(_T("CompressiveStrength property in Pset_PrecastConcreteElementGeneral property set not found"));
      }
      bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->GetSegment(0)->Material.Concrete.Fc = *fc;

      if (!bridge_desc.UseSameGirderForEntireBridge())
      {
         auto type = GetType<Ifc4x3_add2::IfcBeamType>(beam);
         auto girder_name = type->Name().get_value_or(std::string("Unknown"));
         auto girder_library_entry = pLibrary->GetGirderEntry(A2T(girder_name.c_str()));
         if (!girder_library_entry)
         {
            std::ostringstream os;
            os << "Girder type \"" << girder_name << "\" not found in the library";
            IFC_THROW(A2T(os.str().c_str()));
         }
         bridge_desc.GetGirderGroup(spanIdx)->GetGirder(gdrIdx)->SetGirderLibraryEntry(girder_library_entry);
         bridge_desc.SetGirderFamilyName(girder_library_entry->GetGirderFamilyName().c_str());

#pragma Reminder("WORKING HERE - This is assuming the first supported orientation. The IFC file doesn't have this information.")
#pragma Reminder("WORKING HERE - This is duplicate code from above. Simplify")
         auto factory = girder_library_entry->GetBeamFactory();
         auto orientations = factory->GetSupportedGirderOrientation();
         bridge_desc.SetGirderOrientation(orientations.front());
      }
   }
}

void CIfcImporter::ImportSlab(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc)
{
   auto slabs = file.instances_by_type<Ifc4x3_add2::IfcSlab>();
   auto it = std::find_if(slabs->begin(), slabs->end(), [](const auto& slab) {return slab->PredefinedType() == Ifc4x3_add2::IfcSlabTypeEnum::IfcSlabType_FLOOR; });
   if (it == slabs->end())
      return; // no slabs

   auto slab = *it;

   auto stations = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "Stations");
   if (stations.empty())
   {
      IFC_THROW(_T("Stations property in pgsDeck property set not found"));
   }

   auto left_edges = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "LeftEdges");
   if (left_edges.empty())
   {
      IFC_THROW(_T("LeftEdges property in pgsDeck property set not found"));
   }

   auto right_edges = GetPropertyList<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "RightEdges");
   if (right_edges.empty())
   {
      IFC_THROW(_T("RightEdges property in pgsDeck property set not found"));
   }

   if (stations.size() != left_edges.size() || stations.size() != right_edges.size())
   {
      IFC_THROW(_T("Stations, LeftEdges, and RightEdges properties in pgsDeck property set must have the same number of values"));
   }

   auto* pDeck = bridge_desc.GetDeckDescription();

   auto* gross_depth = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "GrossDepth");
   if (gross_depth)
   {
      pDeck->GrossDepth = *gross_depth;
   }
   else
   {
      IFC_THROW(_T("GrossDepth property in pgsDeck property set not found"));
   }

   auto* left_edge_depth = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "LeftEdgeDepth");
   if(left_edge_depth)
   {
      pDeck->OverhangEdgeDepth[pgsTypes::stLeft] = *left_edge_depth;
   }
   else
   {
      IFC_THROW(_T("LeftEdgeDepth property in pgsDeck property set not found"));
   }

   auto* right_edge_depth = GetProperty<Ifc4x3_add2, Ifc4x3_add2::IfcLengthMeasure>(slab, "pgsDeck", "RightEdgeDepth");
   if (right_edge_depth)
   {
      pDeck->OverhangEdgeDepth[pgsTypes::stRight] = *right_edge_depth;
   }
   else
   {
      IFC_THROW(_T("RightEdgeDepth property in pgsDeck property set not found"));
   }

   pDeck->DeckEdgePoints.clear();
   for (auto&& [station, left_edge, right_edge] : boost::combine(stations, left_edges, right_edges))
   {
      CDeckPoint deck_point;
      
      // dummy, default values
      deck_point.MeasurementType = pgsTypes::OffsetMeasurementType::omtBridge;
      deck_point.LeftTransitionType = stations.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;
      deck_point.RightTransitionType = stations.size() == 1 ? pgsTypes::DeckPointTransitionType::dptParallel : pgsTypes::DeckPointTransitionType::dptLinear;

      deck_point.Station = *station;
      deck_point.LeftEdge = *left_edge;
      deck_point.RightEdge = *right_edge;
      pDeck->DeckEdgePoints.emplace_back(deck_point);
   }
}
