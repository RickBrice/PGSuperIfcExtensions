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
#include "Properties.h"

#include <MFCTools\Prompts.h>

#include <IFace/Tools.h>
#include <EAF/AutoProgress.h>

#include <Units\Units.h>

#include <psgLib/BridgeDescription2.h>

// Constants for tracking the state of converting the profile data
#define PROFILE_NOT_STARTED -1
#define PROFILE_ESTABLISHED  1

// Use this throw macro when the data conversion cannot continue
// The catcher, or other, is responsible for deleting it
#define IFC_THROW(_s_) throw new CIfcImporterException(_s_);

Float64 CIfcImporter::m_Precision = 0.001;


HRESULT SameLocation(IPoint2d* pnt1, IPoint2d* pnt2,Float64 tolerance)
{
   Float64 x1, y1, x2, y2;
   pnt1->Location(&x1, &y1);
   pnt2->Location(&x2, &y2);

   return IsEqual(x1, x2, tolerance) && IsEqual(y1, y2, tolerance) ? S_OK : S_FALSE;
}


Ifc4x3_add2::IfcAlignmentHorizontalSegment* GetHorizontalAlignmentSegment(Ifc4x3_add2::IfcAlignmentSegment* alignment_segment)
{
   return alignment_segment->DesignParameters()->as<Ifc4x3_add2::IfcAlignmentHorizontalSegment>();
}

Ifc4x3_add2::IfcAlignmentVerticalSegment* GetVerticalAlignmentSegment(Ifc4x3_add2::IfcAlignmentSegment* alignment_segment)
{
   return alignment_segment->DesignParameters()->as<Ifc4x3_add2::IfcAlignmentVerticalSegment>();
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

CIfcImporter::CIfcImporter(void)
{
   m_pLengthUnit = nullptr;
   m_pAngleUnit = nullptr;

   m_bAlignmentStarted = false;
   m_ProfileState = PROFILE_NOT_STARTED;

   m_CogoEngine.CoCreateInstance(CLSID_CogoEngine);
   m_GeomUtil.CoCreateInstance(CLSID_GeomUtil);

   m_LastAlignmentType = Unknown;
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

bool IsTransitionCurve(Ifc4x3_add2::IfcAlignmentHorizontalSegment* horizontal_segment)
{
   static std::set<Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::Value> transition_curve_types
   {
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE,
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID,
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE,
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC,
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE,
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE,
      Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND
   };

   auto found = transition_curve_types.find(horizontal_segment->PredefinedType());
   return (found == transition_curve_types.end() ? false : true);
}



bool IsCircularCurve(Ifc4x3_add2::IfcAlignmentHorizontalSegment* horizontal_segment)
{
    return horizontal_segment->PredefinedType() == Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC ? true : false;
}

Ifc4x3_add2::IfcAlignmentHorizontal* GetAlignmentHorizontal(Ifc4x3_add2::IfcAlignment* pAlignment)
{
   Ifc4x3_add2::IfcAlignmentHorizontal* horizontal_alignment = nullptr;
   auto nested = pAlignment->IsNestedBy(); // these are the things that are nested by the alignment
   for (auto rel_nests : *nested)
   {
      ATLASSERT(rel_nests->RelatingObject() == pAlignment);
      auto related_objects = rel_nests->RelatedObjects();
      for (auto related_object : *related_objects)
      {
         horizontal_alignment = related_object->as<Ifc4x3_add2::IfcAlignmentHorizontal>();
         if (horizontal_alignment) break;
      }
      if (horizontal_alignment) break;
   }

   return horizontal_alignment;
}

Ifc4x3_add2::IfcAlignmentVertical* GetAlignmentVertical(Ifc4x3_add2::IfcAlignment* pAlignment)
{
   Ifc4x3_add2::IfcAlignmentVertical* vertical_alignment = nullptr;
   auto nested = pAlignment->IsNestedBy(); // these are the things that are nested by the alignment
   for (auto rel_nests : *nested)
   {
      ATLASSERT(rel_nests->RelatingObject() == pAlignment);
      auto related_objects = rel_nests->RelatedObjects();
      for (auto related_object : *related_objects)
      {
         vertical_alignment = related_object->as<Ifc4x3_add2::IfcAlignmentVertical>();
         if (vertical_alignment) break;
      }
      if (vertical_alignment) break;
   }

   return vertical_alignment;
}

bool CIfcImporter::IsValidAlignment(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment)
{
    auto horizontal_alignment = GetAlignmentHorizontal(pAlignment);
    ATLASSERT(horizontal_alignment); // should have found one
    if (!horizontal_alignment)
        return false;

    auto nested = horizontal_alignment->IsNestedBy();
    auto related_objects = (*nested->begin())->RelatedObjects();
    auto nSegments = related_objects->size();

    if (nSegments == 0)
    {
        // alignment doesn't have any segments
        return false;
    }
    else if (nSegments == 1)
    {
        // our model doesn't support isolated transition segments
        // transition curves must be adjacent to circular curves
        auto alignment_segment = (*(related_objects->begin()))->as<Ifc4x3_add2::IfcAlignmentSegment>();
        auto horizontal_segment = GetHorizontalAlignmentSegment(alignment_segment);
        return !IsTransitionCurve(horizontal_segment);
    }
    else if (nSegments == 2)
    {
        // our model doesn't support isolated transition segments
        // transition curves must be adjacent to circular curves
        // can't have two transitions adjacent to each other either
        auto alignment_segment1 = (*(related_objects->begin()))->as<Ifc4x3_add2::IfcAlignmentSegment>();
        auto horizontal_segment1 = GetHorizontalAlignmentSegment(alignment_segment1);
        bool bIsTransitionCurve1 = IsTransitionCurve(horizontal_segment1);


        auto alignment_segment2 = (*(related_objects->begin() + 1))->as<Ifc4x3_add2::IfcAlignmentSegment>();
        auto horizontal_segment2 = GetHorizontalAlignmentSegment(alignment_segment2);
        bool bIsTransitionCurve2 = IsTransitionCurve(horizontal_segment1);

        return (bIsTransitionCurve1 && bIsTransitionCurve2 ? false : true);
    }
    else
    {
        // walk the segments - if we run into a transition curve we have to check the following
        // * transition is adjacent to a circular curve
        // * common radius with circular curve and transition curve are equal
        // * radius of transition curve away from circular curve is infinite
        auto begin = related_objects->begin();
        auto iter = begin;
        auto end = related_objects->end();
        for (; iter != end; iter++)
        {
            auto alignment_segment = (*iter)->as<Ifc4x3_add2::IfcAlignmentSegment>();
            auto horizontal_segment = GetHorizontalAlignmentSegment(alignment_segment);
            if (IsTransitionCurve(horizontal_segment))
            {
                Float64 start_radius = horizontal_segment->StartRadiusOfCurvature();
                Float64 end_radius = horizontal_segment->EndRadiusOfCurvature();
                if (!IsZero(start_radius) && !IsZero(end_radius))
                    return false; // one radius must be zero

                if ((iter == begin && !IsZero(start_radius)) || (iter == end - 1 && !IsZero(end_radius)))
                    return false; // if starting with a transition curve, start radius must be zero or if ending with a transition curve, end radius must be zero

                if (iter != begin && !IsZero(start_radius))
                {
                    // transition starts with a radius so a circular curve must precede this transition curve
                    auto prev_alignment_segment = (*(iter - 1))->as<Ifc4x3_add2::IfcAlignmentSegment>();
                    auto prev_horizontal_segment = GetHorizontalAlignmentSegment(prev_alignment_segment);
                    if (!IsCircularCurve(prev_horizontal_segment)) return false; // previous is not a circular curve

                    Float64 circular_curve_radius = prev_horizontal_segment->StartRadiusOfCurvature();
                    ATLASSERT(IsEqual(circular_curve_radius, prev_horizontal_segment->EndRadiusOfCurvature()));
                    if (!IsEqual(circular_curve_radius, start_radius)) return false; // common radii must be equal

                    if (::BinarySign(start_radius) != ::BinarySign(circular_curve_radius)) return false; // curves must be same direction
                }

                if (iter != end - 1 && !IsZero(end_radius))
                {
                    // transition ends with a radius so a circular curve most come after this transition curve
                    auto next_alignment_segment = (*(iter + 1))->as<Ifc4x3_add2::IfcAlignmentSegment>();
                    auto next_horizontal_segment = GetHorizontalAlignmentSegment(next_alignment_segment);
                    if (!IsCircularCurve(next_horizontal_segment)) return false; // next is not a circular curve

                    Float64 circular_curve_radius = next_horizontal_segment->StartRadiusOfCurvature();
                    ATLASSERT(IsEqual(circular_curve_radius, next_horizontal_segment->EndRadiusOfCurvature()));
                    if (!IsEqual(circular_curve_radius, end_radius)) return false; // common radii must be equal

                    if (::BinarySign(end_radius) != ::BinarySign(circular_curve_radius)) return false; // curves must be same direction
                }
            }
        }
    }

   return true;
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

HRESULT CIfcImporter::ImportFromIFC(std::shared_ptr<WBFL::EAF::Broker> pBroker, CString& strFilePath, CIfcImportOptions options)
{
    USES_CONVERSION;

    std::unique_ptr<IfcParse::IfcFile> pFile = nullptr;

    { // scope the progress window so it closes automatically when we are done with it
        GET_IFACE2(pBroker, IEAFProgress, pProgress);
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
            AfxMessageBox(_T("Unable to parse .ifc file"));
            return S_OK;
        }
    }

   m_Notes.clear();

   auto strSchemaName = pFile->schema()->name();
   if (strSchemaName == std::string("IFC4X3_ADD2"))
   {
      GET_IFACE2(pBroker, IEvents, pEvents);
      pEvents->HoldEvents();

      InitUnits(*pFile);

      if (options.model_elements == CIfcImportOptions::ModelElements::AlignmentOnly)
      {
         ImportAlignment(pBroker, *pFile);
      }
      else
      {
         if (ImportAlignment(pBroker, *pFile))
            ImportBridge(pBroker, *pFile);
      }

      pEvents->FirePendingEvents();
   }
   else
   {
      AfxMessageBox(_T("Schema not supported"));
      ATLASSERT(false); // is there a new schema?
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

bool CIfcImporter::ImportAlignment(std::shared_ptr<WBFL::EAF::Broker> pBroker, IfcParse::IfcFile& file)
{
   AlignmentData2 alignment_data;
   ProfileData2 profile_data;
   RoadwaySectionData section_data;

   if (GetAlignmentParameters(file, &alignment_data, &profile_data, &section_data))
   {
      GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
      pRoadwayData->SetAlignmentData2(alignment_data);
      pRoadwayData->SetProfileData2(profile_data);
      pRoadwayData->SetRoadwaySectionData(section_data);
      return true;
   }
   else
   {
      return false;
   }
}

bool CIfcImporter::ImportBridge(std::shared_ptr<WBFL::EAF::Broker> pBroker, IfcParse::IfcFile& file)
{
   auto bridges = file.instances_by_type<Ifc4x3_add2::IfcBridge>();
   auto bridge = (*bridges->begin());

   SpanIndexType nSpans = INVALID_INDEX;
   GirderIndexType nGirders = INVALID_INDEX;
   auto value = GetProperty<Ifc4x3_add2,Ifc4x3_add2::IfcInteger>(bridge, "usBridge_BridgeCommon", "usBridge_NumberOfSpans");
   if (value)
   {
      nSpans = (SpanIndexType)(*value);
      m_Notes.push_back(std::_tstring(_T("Bridge Number of Spans: ")) + std::to_tstring(nSpans));
   }

   GirderIndexType nTotalGirders = 0;
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
            nTotalGirders++;
         }
      }
   }

   nGirders = nTotalGirders / nSpans;
   m_Notes.push_back(std::_tstring(_T("Bridge Number of Girders: ")) + std::to_tstring(nGirders));

   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   auto bridge_desc = *(pIBridgeDesc->GetBridgeDescription());

   SpanIndexType nSpansToAdd = nSpans - pIBridgeDesc->GetSpanCount();
   if (0 < nSpansToAdd)
   {
      for (SpanIndexType i = 0; i < nSpansToAdd; i++)
      {
         bridge_desc.AppendSpan(nullptr, nullptr, true, 0);
      }

      bridge_desc.UseSameNumberOfGirdersInAllGroups(true);
      bridge_desc.UseSameGirderForEntireBridge(true);
      bridge_desc.SetGirderCount(nGirders);
   }

   //
   // Position the abutments and piers
   //

   // get the abutments and piers and put into a single vector
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
      auto positioning_element = (*rel_positions->begin())->RelatingPositioningElement();
      auto ref = positioning_element->as<Ifc4x3_add2::IfcReferent>();
      auto station = GetProperty<Ifc4x3_add2,Ifc4x3_add2::IfcLengthMeasure>(ref, "Pset_Stationing", "Station");
      pPier->SetStation(*station);
   }


   pIBridgeDesc->SetBridgeDescription(bridge_desc);

   return true;
}

bool CIfcImporter::GetAlignmentParameters(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData)
{
   auto alignment = GetAlignment(file);
   if (alignment == nullptr)
      return false;

   Float64 alignment_adjustment = LoadAlignment(file,alignment);
   *pAlignmentData = m_AlignmentData;

   LoadProfile(file,alignment, alignment_adjustment);
   *pProfileData = m_ProfileData;

#pragma Reminder("WORKING HERE - Roadway Section Data")
   // this is dummy data
   RoadwaySectionTemplate roadway_template;
   roadway_template.LeftSlope = -0.02;
   roadway_template.RightSlope = -0.02;
   roadway_template.Station = 0;
   m_RoadwaySectionData.slopeMeasure = RoadwaySectionData::RelativeToAlignmentPoint;
   m_RoadwaySectionData.NumberOfSegmentsPerSection = 2;
   m_RoadwaySectionData.AlignmentPointIdx = 1;
   m_RoadwaySectionData.ProfileGradePointIdx = 1;
   m_RoadwaySectionData.RoadwaySectionTemplates.push_back(roadway_template);
   *pRoadwaySectionData = m_RoadwaySectionData;

   return true;
}

std::vector<std::_tstring> CIfcImporter::GetNotes()
{
   return m_Notes;
}

Ifc4x3_add2::IfcAlignment* CIfcImporter::GetAlignment(IfcParse::IfcFile& file)
{
   USES_CONVERSION;

   auto alignments = file.instances_by_type<Ifc4x3_add2::IfcAlignment>();
   std::vector<Ifc4x3_add2::IfcAlignment*> valid_alignments;

   for (auto alignment : *alignments)
   {
      if (IsValidAlignment(file,alignment))
         valid_alignments.push_back(alignment);
   }

   if (valid_alignments.size() == 0)
   {
      AfxMessageBox(_T("File does not contain alignments that are compatible with this software."), MB_OK);
   }
   else
   {
      std::ostringstream os;
      for (auto alignment : valid_alignments)
      {
         auto strLabel = (alignment->Name() ? *(alignment->Name()) : alignment->Description() ? *(alignment->Description()) : "Unnamed");
         os << strLabel << std::endl;
      }
      
      int result = 0;
      if (1 < valid_alignments.size()) // prompt to select if more than one alignment
         result = AfxChoose(_T("Select Alignment"), _T("Select alignment to import"), A2T(os.str().c_str()), 0, TRUE);

      if (result < 0)
         return nullptr; // dialog was canceled
      else
         return valid_alignments[result];
   }


   return nullptr;
}

void CIfcImporter::GetStations(Ifc4x3_add2::IfcAlignment* pAlignment, std::vector<std::pair<Float64, Float64>>& vStations, std::vector<std::tuple<Float64, Float64, Float64>>& vStationEquations)
{
   auto nested = pAlignment->IsNestedBy();
   for (auto rel_nests : *nested)
   {
      ATLASSERT(rel_nests->RelatingObject() == pAlignment);
      auto related_objects = rel_nests->RelatedObjects();
      for (auto related_object : *related_objects)
      {
         auto referent = related_object->as<Ifc4x3_add2::IfcReferent>();
         if (referent && referent->PredefinedType() && *(referent->PredefinedType()) == Ifc4x3_add2::IfcReferentTypeEnum::IfcReferentType_STATION)
         {
            Float64 distance_along = 0;
            if (referent->ObjectPlacement())
            {
               auto object_placement = referent->ObjectPlacement();
               auto linear_placement = object_placement->as<Ifc4x3_add2::IfcLinearPlacement>();
               if (linear_placement)
               {
                  // get the distance along the curve for the placement of the referent
                  auto axis2placementlinear = linear_placement->RelativePlacement();
                  auto location = axis2placementlinear->Location();
                  auto point_by_distance_expression = location->as<Ifc4x3_add2::IfcPointByDistanceExpression>();
                  if (point_by_distance_expression)
                  {
                     distance_along = *(point_by_distance_expression->DistanceAlong()->as<Ifc4x3_add2::IfcLengthMeasure>());
                  }
               }
            }

            auto rel_defines_by_properties = referent->IsDefinedBy();
            for (auto rel_defines_property : *rel_defines_by_properties)
            {
               auto property_set = rel_defines_property->RelatingPropertyDefinition()->as<Ifc4x3_add2::IfcPropertySet>();
               if (property_set->Name() == std::string("Pset_Stationing"))
               {
                  bool bHasStation = false;
                  bool bHasIncomingStation = false;
                  Float64 station, incoming_station;
                  auto properties = property_set->HasProperties();
                  for (auto prop : *properties)
                  {
                     if (prop->Name() == "Station")
                     {
                        auto single_value_property = prop->as<Ifc4x3_add2::IfcPropertySingleValue>();
                        if (single_value_property->NominalValue())
                        {
                           bHasStation = true;
                           station = *(single_value_property->NominalValue()->as<Ifc4x3_add2::IfcLengthMeasure>());
                        }
                     }
                     else if (prop->Name() == "IncomingStation")
                     {
                        auto single_value_property = prop->as<Ifc4x3_add2::IfcPropertySingleValue>();
                        if (single_value_property->NominalValue())
                        {
                           bHasIncomingStation = true;
                           incoming_station = *(single_value_property->NominalValue()->as<Ifc4x3_add2::IfcLengthMeasure>());
                        }
                     }
                  }

                  if (bHasStation && bHasIncomingStation)
                  {
                     vStationEquations.emplace_back(distance_along, incoming_station, station);
                  }
                  else if (bHasStation && !bHasIncomingStation)
                  {
                     vStations.emplace_back(distance_along, station);
                  }
                  else
                  {
                     ATLASSERT(false); // not expecting incoming station without station
                  }
               }
            }
         }
      }
   }
}

Float64 CIfcImporter::GetStartStation(Ifc4x3_add2::IfcAlignment* pAlignment)
{
   auto value = GetProperty<Ifc4x3_add2,Ifc4x3_add2::IfcReal>(pAlignment, "Pset_Stationing", "Station");
   if (value)
   {
      return (Float64)(*value);
   }
   else
   {
      return 0.0;
   }
}

Float64 CIfcImporter::LoadAlignment(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment)
{
    USES_CONVERSION;
    m_bAlignmentStarted = false; // the alignment data block has not yet been started

    m_AlignmentData.Name = A2T(pAlignment->Name() ? (*(pAlignment->Name())).c_str() : pAlignment->Description() ? (*(pAlignment->Description())).c_str() : "");

    // initialize the alignment data
    m_AlignmentData.Direction = 0.00;
    m_AlignmentData.xRefPoint = 0.00;
    m_AlignmentData.yRefPoint = 0.00;
    m_AlignmentData.RefStation = 0.00;
    m_AlignmentData.CompoundCurves.clear();

    std::vector<std::pair<Float64, Float64>> vStations; // distance along, station
    std::vector<std::tuple<Float64, Float64, Float64>> vStationEquations; // distance along, incoming station, station
    GetStations(pAlignment, vStations, vStationEquations);

    Float64 station_adjustment = 0;

    if (0 < vStations.size())
    {
       m_AlignmentData.RefStation = vStations.front().second; // .second is the station value
       station_adjustment = m_AlignmentData.RefStation - vStations.front().first; // .first is the distance along value
    }

    if (0 < vStationEquations.size())
    {
       // we don't handle equations yet, but the underlying COGO model does - need up update PGSuper to model equations
    }

    auto horizontal_alignment = GetAlignmentHorizontal(pAlignment);
    ATLASSERT(horizontal_alignment); // should have found one

    Float64 current_station = GetStartStation(pAlignment);
    current_station += station_adjustment;

    // alignment is made up of Line, Spiral, and/or Curve elements
    auto nested = horizontal_alignment->IsNestedBy();
    auto related_objects = (*nested->begin())->RelatedObjects();
    auto begin = related_objects->begin();
    auto iter = begin;
    auto end = related_objects->end();
    for (; iter != end; iter++)
    {
        auto alignment_segment = (*iter)->as<Ifc4x3_add2::IfcAlignmentSegment>();
        bool bIsThereANextSegment = ((iter + 1) != end);

        auto horizontal_alignment_segment = GetHorizontalAlignmentSegment(alignment_segment);

        auto predefined_type = horizontal_alignment_segment->PredefinedType();

        Ifc4x3_add2::IfcAlignmentHorizontalSegment* entrySpiral = nullptr;
        Ifc4x3_add2::IfcAlignmentHorizontalSegment* exitSpiral = nullptr;

        Float64 end_station = current_station;

        if (predefined_type == Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE)
        {
            end_station = OnLine(current_station, horizontal_alignment_segment);
        }
        else if (predefined_type == Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
        {
            // PGSuper can only handle
            // Spiral-Curve
            // Spiral-Curve-Spiral
            // Curve-Spiral
            //
            // Curve-Spiral-Curve, where Spiral is a transition spiral with the start and end radius equal
            // to the curve radii... PGSuper cannot do this case

            entrySpiral = horizontal_alignment_segment;
            if (bIsThereANextSegment)
            {
                // if there is a next segment, check to see if it is a curve
                iter++; // advance to next segment
                auto curve = GetHorizontalAlignmentSegment((*(iter))->as<Ifc4x3_add2::IfcAlignmentSegment>());
                if (curve->PredefinedType() != Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
                    curve = nullptr;

                if (curve)
                {
                    // it's a curve... is there an element that follows the curve?
                    bIsThereANextSegment = ((iter + 1) != end);
                    if (bIsThereANextSegment)
                    {
                        // if there is a next segment, see if it is a spiral
                        exitSpiral = GetHorizontalAlignmentSegment((*(iter + 1))->as<Ifc4x3_add2::IfcAlignmentSegment>());
                        if (exitSpiral->PredefinedType() != Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
                            exitSpiral = nullptr;

                        // if not a spiral, pExitSpiral will be nullptr
                        // this is OK, it just means we have a Spiral-Curve situation
                        if (exitSpiral)
                        {
                            // it is an exit spiral, so advance the iterator
                            iter++;
                            bIsThereANextSegment = ((iter + 1) != end);
                        }
                    }

                    end_station = OnCurve(current_station, entrySpiral, curve, exitSpiral);
                }
                else
                {
                    IFC_THROW(_T("A curve must follow a spiral")); // because PGSuper can't handle it otherwise
                    ATLASSERT(false); // a curve must follow a spiral
                }
            }
            else
            {
                m_Notes.push_back(std::_tstring(_T("Element ignored: The last element in the alignment cannot be a Spiral."))); // PGSuper can't model a lone spiral
            }
        }
        else if (predefined_type == Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
        {
            // looking for Curve-Spiral case
            if (bIsThereANextSegment)
            {
                // check to see if the next element is a spiral
                exitSpiral = GetHorizontalAlignmentSegment((*(iter + 1))->as<Ifc4x3_add2::IfcAlignmentSegment>());
                if (exitSpiral->PredefinedType() != Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
                    exitSpiral = nullptr;

                // if not a spiral, pExitSpiral will be nullptr
                // this is OK, it just means we have a Spiral-Curve situation
                if (exitSpiral)
                {
                    iter++;
                    bIsThereANextSegment = ((iter + 1) != end);
                }

                // Check if the next object is a Curve and if
                // the exit spiral and the curve are touching
                if (bIsThereANextSegment)
                {
                    auto next_curve = GetHorizontalAlignmentSegment((*(iter))->as<Ifc4x3_add2::IfcAlignmentSegment>());
                    if (next_curve->PredefinedType() != Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
                        next_curve = nullptr;

                    if (next_curve && exitSpiral)
                    {
                        CComPtr<IPoint2d> pntSpiralStart, pntSpiralPI, pntSpiralEnd;
                        GetSpiralPoints(exitSpiral, &pntSpiralStart, &pntSpiralPI, &pntSpiralEnd);

                        CComPtr<IPoint2d> pntNextCurveStart, pntNextCurvePI, pntNextCurveEnd, pntNextCurveCenter;
                        GetCurvePoints(next_curve, &pntNextCurveStart, &pntNextCurvePI, &pntNextCurveEnd, &pntNextCurveCenter);

                        if (SameLocation(pntSpiralEnd, pntNextCurveStart, m_Precision) == S_FALSE)
                        {
                            IFC_THROW(_T("PGSuper cannot model a transition spiral between two circular curves."));
                        }
                    }
                }
            }

            ATLASSERT(entrySpiral == nullptr); // this must be the case
            end_station = OnCurve(current_station, entrySpiral, horizontal_alignment_segment, exitSpiral);
        }
        else
        {
            ATLASSERT(false);
        }
        current_station = end_station;
    }

    return station_adjustment;
}

void CIfcImporter::LoadProfile(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment, Float64 stationAdjustment)
{
   m_ProfileState = PROFILE_NOT_STARTED;
   m_ProfileData.Station = 0;
   m_ProfileData.Elevation = 0;
   m_ProfileData.Grade = 0;
   m_ProfileData.VertCurves.clear();


   auto horizontal_alignment = GetAlignmentHorizontal(pAlignment);
   auto vertical_alignment = GetAlignmentVertical(pAlignment);

   Float64 start_station = GetStartStation(pAlignment);
   start_station += stationAdjustment;

   if(vertical_alignment)
   {
       auto nested = vertical_alignment->IsNestedBy();
       ATLASSERT((*nested).size() == 1);
       auto related_objects = (*nested->begin())->RelatedObjects();

      if (related_objects->size() == 0)
      {
         // the profile geometry list is empty so assume a flat grade
         m_Notes.push_back(std::_tstring(_T("A profile was not found or the profile does not contain segments. Assuming a default profile.")));
         m_ProfileData.Station = start_station;
         m_ProfileData.Elevation = 0.0;
         m_ProfileData.Grade = 0.0;
         return;
      }

      auto begin = related_objects->begin();
      auto iter = begin;
      auto end = related_objects->end();
      for (; iter != end; iter++)
      {
          auto alignment_segment = (*iter)->as<Ifc4x3_add2::IfcAlignmentSegment>();
          auto vertical_alignment_segment = GetVerticalAlignmentSegment(alignment_segment);
          auto predefined_type = vertical_alignment_segment->PredefinedType();

          if (predefined_type == Ifc4x3_add2::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT)
          {
              OnLinearSegment(start_station,vertical_alignment_segment);
          }
          else if (predefined_type == Ifc4x3_add2::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC)
          {
              OnParabolicSegment(start_station,vertical_alignment_segment);
          }
          else if (predefined_type == Ifc4x3_add2::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CIRCULARARC)
          {
#pragma Reminder("WORKING HERE - Need to deal with vertical circular arcs") // treat it as a parabola for now
             OnParabolicSegment(start_station, vertical_alignment_segment);
          }
          else if (predefined_type == Ifc4x3_add2::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CLOTHOID)
          {
#pragma Reminder("WORKING HERE - Need to deal with vertical clothoid arcs") // treat it as a parabola for now
             OnParabolicSegment(start_station, vertical_alignment_segment);
          }
          else
          {
             ATLASSERT(false); // is there a new type ???
                                // TODO: provide a better exception
             IFC_THROW(_T("An unknown profile element was encountered"));
          }
      }
   }

   if (m_ProfileState != PROFILE_ESTABLISHED)
   {
      // we are out of elements and the profile definition is not finished

      // THIS IS A GOOD PLACE TO GENERATE INFORMATION MESSAGES ABOUT ANY ASSUMPTIONS
      // THIS IMPORTER HAD TO MAKE
      if (m_ProfileData.VertCurves.size() == 0)
      {
         // A second element was not provided so the main grade
         // could not be established... use the default value of 0
         m_Notes.push_back(std::_tstring(_T("More elements are needed to determine the starting grade of the profile. Assuming a grade of 0.0%")));
      }
      else
      {
         // The exit profile of the last vertical curve could not be established...
         // We could either pop the last vertical curve out of the list or use the
         // default exit grade of 0.
         //
         // Use the default exit grade of 0.
         m_Notes.push_back(std::_tstring(_T("More elements are needed to determine the exit grade of the last vertical curve. Assuming a grade of 0.0%")));
      }
   }
}

Float64 CIfcImporter::OnLine(Float64 startStation, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pLine)
{
   Float64 sx, sy;
   GetPoint(pLine->StartPoint(), &sx, &sy);

   Float64 length = WBFL::Units::ConvertToSysUnits(pLine->SegmentLength(),*m_pLengthUnit);
   Float64 startDirection = WBFL::Units::ConvertToSysUnits(pLine->StartDirection(), *m_pAngleUnit);;
   return OnLine(sx, sy, startStation, startDirection, length);
}

Float64 CIfcImporter::OnLine(Float64 sx,Float64 sy,Float64 startStation,Float64 startDirection, Float64 length)
{
   Float64 end_station = startStation + length;

   if (!m_bAlignmentStarted)
   {
      // the bridge starts somewhere in this line segment
      m_AlignmentData.RefStation = startStation;
      m_AlignmentData.xRefPoint = sx;
      m_AlignmentData.yRefPoint = sy;
      m_AlignmentData.Direction = startDirection;

      if (IsZero(m_AlignmentData.Direction))
         m_AlignmentData.Direction = 0.0;
      else if (m_AlignmentData.Direction < 0)
         m_AlignmentData.Direction += TWO_PI;

      m_bAlignmentStarted = true;
   }
   else
   {
      if (m_LastAlignmentType != Curve)
      {
         // add an angle point
         CompoundCurveData hcData;
         hcData.PIStation = startStation;
         hcData.Radius = 0;
         hcData.EntrySpiral = 0;
         hcData.ExitSpiral = 0;
         hcData.bFwdTangent = true;
         hcData.FwdTangent = startDirection;

         m_AlignmentData.CompoundCurves.push_back(hcData);
      }
   }

   m_LastAlignmentType = Line;

   return end_station;
}

Float64 CIfcImporter::OnCurve(Float64 startStation, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pEntrySpiral, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pCurve, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pExitSpiral)
{
   ATLASSERT(pCurve != nullptr);

   Float64 radius = WBFL::Units::ConvertToSysUnits(pCurve->StartRadiusOfCurvature(), *m_pLengthUnit);

   // Get all the construction points
   CComPtr<IPoint2d> pntEntryStart, pntEntryPI, pntEntryEnd;
   CComPtr<IPoint2d> pntExitStart, pntExitPI, pntExitEnd;
   CComPtr<IPoint2d> pntCurveStart, pntCurvePI, pntCurveEnd, pntCurveCenter;

   Float64 entry_spiral_length = 0;
   Float64 exit_spiral_length = 0;
   if (pEntrySpiral)
   {
      GetSpiralPoints(pEntrySpiral, &pntEntryStart, &pntEntryPI, &pntEntryEnd);
      entry_spiral_length = WBFL::Units::ConvertToSysUnits(pEntrySpiral->SegmentLength(), *m_pLengthUnit);

      if (pntEntryStart == nullptr || pntEntryPI == nullptr || pntEntryEnd == nullptr)
      {
         m_Notes.push_back(std::_tstring(_T("Entry spiral ignored.")));
         pntEntryStart.Release();
         pntEntryPI.Release();
         pntEntryEnd.Release();
         pEntrySpiral = nullptr;
         entry_spiral_length = 0;
      }
      else
      {
         if (!IsZero(pEntrySpiral->StartRadiusOfCurvature()))
         {
            m_Notes.push_back(std::_tstring(_T("Start radius of entry spiral taken to be infinite")));
         }
         CheckSpiralType(pEntrySpiral);
      }
   }

   GetCurvePoints(pCurve, &pntCurveStart, &pntCurvePI, &pntCurveEnd, &pntCurveCenter);

   if (pntCurveStart == nullptr || pntCurvePI == nullptr || pntCurveEnd == nullptr || pntCurveCenter == nullptr)
   {
      m_Notes.push_back(std::_tstring(_T("Zero radius curve could not be constructed.")));
      return startStation;
   }

   if (pExitSpiral)
   {
      GetSpiralPoints(pExitSpiral, &pntExitStart, &pntExitPI, &pntExitEnd);
      exit_spiral_length = WBFL::Units::ConvertToSysUnits(pExitSpiral->SegmentLength(), *m_pLengthUnit);

      if (pntExitStart == nullptr || pntExitPI == nullptr || pntExitEnd == nullptr)
      {
         m_Notes.push_back(std::_tstring(_T("Exit spiral ignored.")));
         pntExitStart.Release();
         pntExitPI.Release();
         pntExitEnd.Release();
         pExitSpiral = nullptr;
         exit_spiral_length = 0;
      }
      else
      {
         if (!IsZero(pExitSpiral->EndRadiusOfCurvature()))
         {
            m_Notes.push_back(std::_tstring(_T("End radius of exit spiral taken to be infinite")));
         }
         CheckSpiralType(pExitSpiral);
      }
   }

   Float64 sx, sy;
   pntCurveStart->get_X(&sx);
   pntCurveStart->get_Y(&sy);

#if defined _DEBUG
   // check radius based on the circular curve parameters
   Float64 ex, ey;
   Float64 cx, cy;
   Float64 px, py;

   pntCurvePI->get_X(&px);
   pntCurvePI->get_Y(&py);

   pntCurveEnd->get_X(&ex);
   pntCurveEnd->get_Y(&ey);

   pntCurveCenter->get_X(&cx);
   pntCurveCenter->get_Y(&cy);

   Float64 dx = sx - cx;
   Float64 dy = sy - cy;
   ATLASSERT(IsEqual(fabs(radius), sqrt(dx*dx + dy*dy)));
#endif // _DEBUG

   // Determine the control points
   CComPtr<IPoint2d> pntStart, pntPI, pntEnd;
   if (pEntrySpiral && !pExitSpiral)
   {
      // Spiral-Curve
      pntStart = pntEntryStart;

      // PI is at the intersection of the forward and back tangents
      CComPtr<IIntersect2> intersect;
      m_CogoEngine->get_Intersect(&intersect);
      intersect->LinesByPoints(pntEntryStart, pntEntryPI, 0.0, pntCurvePI, pntCurveEnd, 0.0, &pntPI);

      pntEnd = pntCurveEnd;

      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pEntrySpiral->EndRadiusOfCurvature(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("End radius of the entry spiral does not match the radius of the circular curve. The entry spiral end radius will be ignored.")));
      }

      if (SameLocation(pntEntryEnd, pntCurveStart, m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The end of the entry spiral does not coincide with the start of the circular curve. The end of the entry spiral has been adjusted.")));
      }
   }
   else if (!pEntrySpiral && pExitSpiral)
   {
      // Curve-Spiral
      pntStart = pntCurveStart;

      // PI is at the intersection of the forward and back tangents
      CComPtr<IIntersect2> intersect;
      m_CogoEngine->get_Intersect(&intersect);
      intersect->LinesByPoints(pntStart, pntCurvePI, 0.0, pntExitPI, pntExitEnd, 0.0, &pntPI);

      pntEnd = pntExitEnd;

      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pExitSpiral->StartRadiusOfCurvature(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("Start radius of the exit spiral does not match the radius of the circular curve. The exit spiral start radius will be ignored.")));
      }

      if (SameLocation(pntCurveEnd, pntExitStart, m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The start of the exit spiral does not coincide with the end of the circular curve. The exit spiral has been adjusted.")));
      }
   }
   else if (pEntrySpiral && pExitSpiral)
   {
      // Spiral-Curve-Spiral
      pntStart = pntEntryStart;

      CComPtr<IIntersect2> intersect;
      m_CogoEngine->get_Intersect(&intersect);
      intersect->LinesByPoints(pntEntryStart, pntEntryPI, 0.0, pntExitPI, pntExitEnd, 0.0, &pntPI);

      ATLASSERT(pntPI->SameLocation(pntCurvePI));

      pntEnd = pntExitEnd;


      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pEntrySpiral->EndRadiusOfCurvature(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("End radius of the entry spiral does not match the radius of the circular curve. The entry spiral end radius will be ignored.")));
      }

      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pExitSpiral->StartRadiusOfCurvature(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("Start radius of the exit spiral does not match the radius of the circular curve. The exit spiral start radius will be ignored.")));
      }

      if (SameLocation(pntEntryEnd, pntCurveStart, m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The end of the entry spiral does not coincide with the start of the circular curve. The entry spiral has been adjusted.")));
      }

      if (SameLocation(pntCurveEnd, pntExitStart, m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The start of the exit spiral does not coincide with the end of the circular curve. The exit spiral has been adjusted.")));
      }
   }
   else
   {
      // Curve
      pntStart = pntCurveStart;
      pntPI = pntCurvePI;
      pntEnd = pntCurveEnd;
   }

   // create a horizontal curve object so that we can get some information from it
   CComPtr<ICompoundCurve> hc;
   hc.CoCreateInstance(CLSID_CompoundCurve);
   hc->put_PBT(pntStart);
   hc->put_PI(pntPI);
   hc->put_PFT(pntEnd);
   hc->put_Radius(fabs(radius));
   hc->put_SpiralLength(spEntry, entry_spiral_length);
   hc->put_SpiralLength(spExit, exit_spiral_length);

   Float64 length;
   hc->get_TotalLength(&length);

   CComPtr<IAngle> objAngle;
   hc->get_CircularCurveAngle(&objAngle);

   Float64 end_station = startStation + length;

   if (!m_bAlignmentStarted)
   {
      // this is the first element... start the alignment data
      m_AlignmentData.RefStation = startStation;
      m_AlignmentData.xRefPoint = sx;
      m_AlignmentData.yRefPoint = sy;

      CComPtr<IDirection> dirBkTangent;
      hc->get_BkTangentBrg(&dirBkTangent);
      dirBkTangent->get_Value(&m_AlignmentData.Direction);

      m_bAlignmentStarted = true;
   }

   CompoundCurveData hcData;
   hcData.EntrySpiral = entry_spiral_length;
   hcData.ExitSpiral = exit_spiral_length;
   hcData.Radius = fabs(radius);

   Float64 tangent;
   hc->get_BkTangentLength(&tangent);
   hcData.PIStation = startStation + tangent;

   CComPtr<IAngle> curve_angle;
   hc->get_CurveAngle(&curve_angle);

   Float64 delta;
   curve_angle->get_Value(&delta);

   CurveDirectionType dir;
   hc->get_Direction(&dir);

   if (dir == cdRight)
      delta *= -1;

   hcData.bFwdTangent = false;
   hcData.FwdTangent = delta;

   m_AlignmentData.CompoundCurves.push_back(hcData);

   m_LastAlignmentType = Curve;

   return end_station;
}

void CIfcImporter::GetCurvePoints(Ifc4x3_add2::IfcAlignmentHorizontalSegment* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter)
{
   auto pStart = pCurve->StartPoint();
   auto bkTangentBrg = WBFL::Units::ConvertToSysUnits(pCurve->StartDirection(), *m_pAngleUnit);
   auto L = WBFL::Units::ConvertToSysUnits(pCurve->SegmentLength(), *m_pLengthUnit);
   auto R = WBFL::Units::ConvertToSysUnits(pCurve->StartRadiusOfCurvature(), *m_pLengthUnit);
   bool bIsCCW = (R < 0 ? false : true);

   Float64 delta = fabs(L / R);
   Float64 T = R*tan(delta / 2);

   Float64 sx, sy;
   GetPoint(pStart, &sx, &sy);
   CComPtr<IPoint2d> pntStart;
   pntStart.CoCreateInstance(CLSID_Point2d);
   pntStart->Move(sx, sy);
   pntStart.CopyTo(ppStart);

   CComPtr<ILocate2> locate;
   m_CogoEngine->get_Locate(&locate);

   locate->ByDistDir(*ppStart, T, CComVariant(bkTangentBrg), 0.0, ppPI);

   Float64 fwdTangentBrg = bkTangentBrg + (bIsCCW ? 1 : -1)*delta;

   locate->ByDistDir(*ppPI, T, CComVariant(fwdTangentBrg), 0.0, ppEnd);

   locate->ByDistDir(*ppStart, R, CComVariant(bkTangentBrg + (bIsCCW ? 1 : -1)*PI_OVER_2), 0.0, ppCenter);
}

Float64 SpiralX(Float64 ls, Float64 angle)
{
   return ls*(1 - pow(angle, 2) / 10 + pow(angle, 4) / 216 - pow(angle, 6) / 9360);
}

Float64 SpiralY(Float64 ls, Float64 angle)
{
   return ls*(angle / 3 - pow(angle, 3) / 42 + pow(angle, 5) / 1320 - pow(angle, 7) / 75600);
}

void CIfcImporter::GetSpiralPoints(Ifc4x3_add2::IfcAlignmentHorizontalSegment* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd)
{
   auto pStart = pSpiral->StartPoint();
   auto bkTangentBrg = WBFL::Units::ConvertToSysUnits(pSpiral->StartDirection(), *m_pAngleUnit);
   auto L = WBFL::Units::ConvertToSysUnits(pSpiral->SegmentLength(), *m_pLengthUnit);
   auto Rstart = WBFL::Units::ConvertToSysUnits(pSpiral->StartRadiusOfCurvature(), *m_pLengthUnit);
   auto Rend = WBFL::Units::ConvertToSysUnits(pSpiral->EndRadiusOfCurvature(), *m_pLengthUnit);
   auto R = IsZero(Rstart) ? Rend : Rstart; // zero means infinite radius
   bool bIsCCW = (R < 0 ? true : false);

   Float64 sx, sy;
   GetPoint(pStart, &sx, &sy);
   CComPtr<IPoint2d> pntStart;
   pntStart.CoCreateInstance(CLSID_Point2d);
   pntStart->Move(sx, sy);
   pntStart.CopyTo(ppStart);

   Float64 DE = L / (2 * R); // deflection angle
   Float64 X = SpiralX(L, DE);
   Float64 Y = SpiralY(L, DE);
   Float64 v = Y / sin(DE); // short tangent
   Float64 u = X - Y / tan(DE); // long tangent

   CComPtr<ILocate2> locate;
   m_CogoEngine->get_Locate(&locate);

   if (IsZero(Rend)) // zero means infinite radius at end so working from start of spiral
   {
      locate->ByDistDir(*ppStart, v, CComVariant(bkTangentBrg), 0.0, ppPI);
      locate->ByDistDir(*ppPI, u, CComVariant(bkTangentBrg + (bIsCCW ? 1 : -1)*DE), 0.0, ppEnd);
   }
   else
   {
      locate->ByDistDir(*ppStart, u, CComVariant(bkTangentBrg), 0.0, ppPI);
      locate->ByDistDir(*ppPI, v, CComVariant(bkTangentBrg + (bIsCCW ? 1 : -1)*DE), 0.0, ppEnd);
   }
}

void CIfcImporter::OnLinearSegment(Float64 startStation,Ifc4x3_add2::IfcAlignmentVerticalSegment* pLinearSegment)
{
   Float64 length = WBFL::Units::ConvertToSysUnits(pLinearSegment->HorizontalLength(),*m_pLengthUnit);
   Float64 start_gradient = pLinearSegment->StartGradient();
   Float64 start_dist = WBFL::Units::ConvertToSysUnits(pLinearSegment->StartDistAlong(),*m_pLengthUnit);

   Float64 start_height = WBFL::Units::ConvertToSysUnits(pLinearSegment->StartHeight(),*m_pLengthUnit);

   if (m_ProfileState == PROFILE_NOT_STARTED)
   {
      m_ProfileData.Station = startStation + start_dist;
      m_ProfileData.Elevation = start_height;
      m_ProfileData.Grade = start_gradient;
      m_ProfileState = PROFILE_ESTABLISHED;
   }
   else
   {
      // PGSuper models linear segments as zero length vertical curves
      if (m_ProfileData.VertCurves.size() == 0 || !IsEqual(m_ProfileData.VertCurves.back().ExitGrade, start_gradient))
      {
         VertCurveData vcData;
         vcData.PVIStation = startStation + start_dist + length / 2;
         vcData.L1 = length;
         vcData.L2 = 0;
         vcData.ExitGrade = start_gradient;

         m_ProfileData.VertCurves.push_back(vcData);

         m_ProfileState = PROFILE_ESTABLISHED;
      }
   }
}

void CIfcImporter::OnParabolicSegment(Float64 startStation, Ifc4x3_add2::IfcAlignmentVerticalSegment* pParaCurve)
{
   // finish any open profile element
   Float64 start_gradient = pParaCurve->StartGradient();
   Float64 start_dist = WBFL::Units::ConvertToSysUnits(pParaCurve->StartDistAlong(), *m_pLengthUnit);
   Float64 start_height = WBFL::Units::ConvertToSysUnits(pParaCurve->StartHeight(), *m_pLengthUnit);
   Float64 length = WBFL::Units::ConvertToSysUnits(pParaCurve->HorizontalLength(), *m_pLengthUnit);
   Float64 exit_gradient = pParaCurve->EndGradient();

   if (m_ProfileState == PROFILE_NOT_STARTED)
   {
      m_ProfileData.Station = startStation + start_dist;
      m_ProfileData.Elevation = start_height;
      m_ProfileData.Grade = start_gradient;
      m_ProfileState = PROFILE_ESTABLISHED;
   }

   // add this vertical curve
   VertCurveData vcData;
   vcData.PVIStation = startStation + start_dist + length / 2;
   vcData.L1 = length;
   vcData.L2 = 0;
   vcData.ExitGrade = exit_gradient;

   m_ProfileData.VertCurves.push_back(vcData);
}

void CIfcImporter::CheckSpiralType(Ifc4x3_add2::IfcAlignmentHorizontalSegment* pSpiral)
{
   switch (pSpiral->PredefinedType())
   {
   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE:
   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE:
   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC:
   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE:
   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE:
   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND:
      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid.")));
      break;

   case Ifc4x3_add2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID:
      // this is ok... we were expecting clothoid
      break;

   default:
      ATLASSERT(false); // is there a new spiral type???
      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid.")));
      break;
   }
}