///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2023  Washington State Department of Transportation
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
#include "IfcAlignmentConverter.h"
#include "IfcAlignmentConverterException.h"

#include <MFCTools\Prompts.h>
#include <EAF\EAFAutoProgress.h>

#include <Units\Units.h>

// Constants for tracking the state of converting the profile data
#define PROFILE_NOT_STARTED -1
#define PROFILE_ESTABLISHED  1

// Use this throw macro when the data conversion cannot continue
// The catcher, or other, is responsible for deleting it
#define IFC_THROW(_s_) throw new CIfcAlignmentConverterException(_s_);

Float64 CIfcAlignmentConverter::m_Precision = 0.001;


HRESULT SameLocation(IPoint2d* pnt1, IPoint2d* pnt2,Float64 tolerance)
{
   Float64 x1, y1, x2, y2;
   pnt1->Location(&x1, &y1);
   pnt2->Location(&x2, &y2);

   return IsEqual(x1, x2, tolerance) && IsEqual(y1, y2, tolerance) ? S_OK : S_FALSE;
}

CIfcAlignmentConverter::CIfcAlignmentConverter(void)
{
   m_pLengthUnit = nullptr;
   m_pAngleUnit = nullptr;

   m_bAlignmentStarted = false;
   m_ProfileState = PROFILE_NOT_STARTED;

   m_CogoEngine.CoCreateInstance(CLSID_CogoEngine);
   m_GeomUtil.CoCreateInstance(CLSID_GeomUtil);

   m_LastAlignmentType = Unknown;
}

CIfcAlignmentConverter::~CIfcAlignmentConverter(void)
{
}

template <typename Schema>
void CIfcAlignmentConverter::InitUnits(IfcParse::IfcFile& file)
{
   auto geometric_representation_contexts = file.instances_by_type<Schema::IfcGeometricRepresentationContext>();
   auto geometric_representation_context = (0 < geometric_representation_contexts->size()) ? *(geometric_representation_contexts->begin()) : nullptr;
#pragma Reminder("WORKING HERE - There could be multiple geometric representation contexts, how do we know if we have the right one?")
   if (geometric_representation_context && geometric_representation_context->Precision() != boost::none)
   {
      m_Precision = *(geometric_representation_context->Precision());
   }

#pragma Reminder("WORKING HERE - UNITS - THERE ARE MANY CASES THIS DOESN'T DEAL WITH")
   auto unit_assignment_instances = file.instances_by_type<Schema::IfcUnitAssignment>();
   ATLASSERT(unit_assignment_instances->size() == 1);
   auto unit_assignment = *(unit_assignment_instances->begin());
   auto units = unit_assignment->Units();
   for (auto unit : *units)
   {
      auto derived_unit = unit->as<Schema::IfcDerivedUnit>();
      auto monitary_unit = unit->as<Schema::IfcMonetaryUnit>();
      auto si_unit = unit->as<Schema::IfcSIUnit>();
      auto conversion_based_unit = unit->as<Schema::IfcConversionBasedUnit>();
      auto conversion_based_unit_with_offset = unit->as<Schema::IfcConversionBasedUnitWithOffset>();

      if (si_unit)
      {
         if (si_unit->Name() == Schema::IfcSIUnitName::IfcSIUnitName_METRE)
         {
            if (si_unit->Prefix() != boost::none)
            {
               switch (*(si_unit->Prefix()))
               {
               case Schema::IfcSIPrefix::IfcSIPrefix_KILO:
                  m_pLengthUnit = &WBFL::Units::Measure::Kilometer;
                  break;

               case Schema::IfcSIPrefix::IfcSIPrefix_CENTI:
                  m_pLengthUnit = &WBFL::Units::Measure::Centimeter;
                  break;

               case Schema::IfcSIPrefix::IfcSIPrefix_MILLI:
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

         if (si_unit->Name() == Schema::IfcSIUnitName::IfcSIUnitName_RADIAN)
         {
            ATLASSERT(si_unit->Prefix() == boost::none); // not expecting anything like Kilo-radians
            m_pAngleUnit = &WBFL::Units::Measure::Radian;
            continue;
         }
      }

      if (conversion_based_unit)
      {
         if (conversion_based_unit->UnitType() == Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT)
         {
            if (conversion_based_unit->UnitType() == Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT)
            {
               auto measure_with_unit = conversion_based_unit->ConversionFactor();
               auto unit_component = measure_with_unit->UnitComponent()->as<Schema::IfcSIUnit>();
               ATLASSERT(unit_component->Name() == Schema::IfcSIUnitName::IfcSIUnitName_RADIAN);
               ATLASSERT(unit_component->Prefix() == boost::none); // not dealing with conversion factors to anything but meter
               Float64 conversion_factor;
               try
               {
                  auto value_component = measure_with_unit->ValueComponent();
                  ATLASSERT(value_component); // not dealing with anything but simple conversion factors
                   //auto value = *(value_component->as<Schema::IfcLengthMeasure>()); // as<> returns nullptr if the type is different so deferencing could lead to crash
                   //auto value = *(value_component->as<Schema::IfcRatio>()); // as<> returns nullptr if the type is different so deferencing could lead to crash
                  auto value = static_cast<Float64>(*value_component->data().getArgument(0));
                  conversion_factor = 1 / (value);
               }
               catch (IfcParse::IfcInvalidTokenException& e)
               {
                  // Was expecting something like 
                  // #15 = IFCMEASUREWITHUNIT(IFCLENGTHMEASURE(3.28083333333333), #16);
                  // where the expected token is IFCLENGTHMEASURE, but instead found something like
                  // #15=IFCMEASUREWITHUNIT(3.28083333333333,#16);
                  // we'll just get the value and keep going
                  TRACE(e.what());
                  Argument* pArgument = measure_with_unit->get("ValueComponent");
                  ATLASSERT(pArgument->type() == IfcUtil::Argument_DOUBLE);
                  double value = double(*pArgument);
                  conversion_factor = 1 / value;
               }

               if (IsEqual(conversion_factor, WBFL::Units::Measure::Degree.GetConvFactor()))
               {
                  m_pAngleUnit = &WBFL::Units::Measure::Degree;
               }
            }
            else if (conversion_based_unit->UnitType() == Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT)
            {
               auto measure_with_unit = conversion_based_unit->ConversionFactor();
               auto unit_component = measure_with_unit->UnitComponent()->as<Schema::IfcSIUnit>();
               ATLASSERT(unit_component->Name() == Schema::IfcSIUnitName::IfcSIUnitName_METRE);
               ATLASSERT(unit_component->Prefix() == boost::none); // not dealing with conversion factors to anything but meter

               Float64 conversion_factor;
               try
               {
                  auto value_component = measure_with_unit->ValueComponent();
                  ATLASSERT(value_component); // not dealing with anything but simple conversion factors
                  //auto value = *(value_component->as<Schema::IfcLengthMeasure>()); // as<> returns nullptr if the type is different so deferencing could lead to crash
                  //auto value = *(value_component->as<Schema::IfcRatio>()); // as<> returns nullptr if the type is different so deferencing could lead to crash
                  auto value = static_cast<Float64>(*value_component->data().getArgument(0));
                  conversion_factor = (value);
               }
               catch (IfcParse::IfcInvalidTokenException& e)
               {
                  // Was expecting something like 
                  // #15 = IFCMEASUREWITHUNIT(IFCLENGTHMEASURE(3.28083333333333), #16);
                  // where the expected token is IFCLENGTHMEASURE, but instead found something like
                  // #15=IFCMEASUREWITHUNIT(3.28083333333333,#16);
                  // we'll just get the value and keep going
                  TRACE(e.what());
                  Argument* pArgument = measure_with_unit->get("ValueComponent");
                  ATLASSERT(pArgument->type() == IfcUtil::Argument_DOUBLE);
                  double value = double(*pArgument);
                  conversion_factor = value;
               }

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
}

template <typename Schema>
bool  CIfcAlignmentConverter::IsValidAlignment(typename Schema::IfcAlignment* pAlignment)
{
   auto axis = pAlignment->Axis()->as<Schema::IfcAlignmentCurve>();
   auto horizontal = axis->Horizontal();
   auto segments = horizontal->Segments();
   auto nSegments = segments->size();
   if (nSegments == 0)
   {
      // alignment doesn't have any segments
      return false;
   }
   else if (nSegments == 1)
   {
      // our model doesn't support isolated transition segments
      // transition curves must be adjacent to circular curves
      auto segment = (*segments->begin());
      auto transition = segment->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
      return (transition == nullptr ? true : false);
   }
   else if (nSegments == 2)
   {
      // our model doesn't support isolated transition segments
      // transition curves must be adjacent to circular curves
      // can't have two transitions adjacent to each other either
      auto segment1 = (*segments->begin());
      auto transition1 = segment1->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();

      auto segment2 = (*segments->begin() + 1);
      auto transition2 = segment2->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
      return (transition1 == nullptr || transition2 == nullptr ? true : false);
   }
   else
   {
      // walk the segments - if we run into a transition curve we have to check the following
      // * transition is adjacent to a circular curve
      // * common radius with circular curve and transition curve are equal
      // * radius of transition curve away from circular curve is infinite
      auto begin = segments->begin();
      auto iter = begin;
      auto end = segments->end();
      for (; iter != end; iter++)
      {
         auto segment(*iter);
         auto transition = segment->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
         if (transition)
         {
            Float64 start_radius = (transition->hasStartRadius() ? transition->StartRadius() : 0.0);
            Float64 end_radius = (transition->hasEndRadius() ? transition->EndRadius() : 0.0);
            if (!IsZero(start_radius) && !IsZero(end_radius))
               return false; // one radius must be zero

            if ((iter == begin && !IsZero(start_radius)) || (iter == end - 1 && !IsZero(end_radius)))
               return false; // if starting with a transition curve, start radius must be zero or if ending with a transition curve, end radius must be zero

            if (iter != begin && !IsZero(start_radius))
            {
               // transition starts with a radius so a circular curve must preceed this transition curve
               auto arc = (*(iter - 1))->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
               if (!arc) return false; // previous is not an arc

               if (!IsEqual(arc->Radius(), start_radius)) return false; // common radii must be equal

               if (arc->IsCCW() != transition->IsStartRadiusCCW()) return false; // curves must be same direction
            }

            if (iter != end - 1 && !IsZero(end_radius))
            {
                // transition ends with a radius so a circular curve most come after this transition curve
                auto arc = (*(iter + 1))->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
                if (!arc) return false; // previous is not an arc

                if (!IsEqual(arc->Radius(), end_radius)) return false; // common radii must be equal

                if (arc->IsCCW() != transition->IsEndRadiusCCW()) return false; // curves must be same direction
            }
         }
      }
   }

   return true;
}

//bool IsTransitionCurve(Ifc4x3_rc2::IfcAlignmentHorizontalSegment* horizontal_segment)
//{
//    static std::set<Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::Value> transition_curve_types
//    {
//       Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE,
//       Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID,
//       Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE,
//       //Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC,
//       Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL,
//       //Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE,
//       Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE,
//       Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND
//    };
//
//    auto found = transition_curve_types.find(horizontal_segment->PredefinedType());
//    return (found == transition_curve_types.end() ? false : true);
//}
//
//bool IsTransitionCurve(Ifc4x3_rc3::IfcAlignmentHorizontalSegment* horizontal_segment)
//{
//    static std::set<Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::Value> transition_curve_types
//    {
//       Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE,
//       Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID,
//       Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE,
//       Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC,
//       //Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL,
//       //Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE,
//       Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE,
//       Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND
//    };
//
//    auto found = transition_curve_types.find(horizontal_segment->PredefinedType());
//    return (found == transition_curve_types.end() ? false : true);
//}
//
//bool IsTransitionCurve(Ifc4x3_rc4::IfcAlignmentHorizontalSegment* horizontal_segment)
//{
//    static std::set<Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::Value> transition_curve_types
//    {
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE,
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID,
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE,
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC,
//       //Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL,
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE,
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE,
//       Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND
//    };
//
//    auto found = transition_curve_types.find(horizontal_segment->PredefinedType());
//    return (found == transition_curve_types.end() ? false : true);
//}

bool IsTransitionCurve(Ifc4x3_tc1::IfcAlignmentHorizontalSegment* horizontal_segment)
{
   static std::set<Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::Value> transition_curve_types
   {
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE,
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID,
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE,
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC,
      //Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL,
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE,
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE,
      Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND
   };

   auto found = transition_curve_types.find(horizontal_segment->PredefinedType());
   return (found == transition_curve_types.end() ? false : true);
}

bool IsTransitionCurve(Ifc4x3_add1::IfcAlignmentHorizontalSegment* horizontal_segment)
{
   static std::set<Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::Value> transition_curve_types
   {
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE,
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID,
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE,
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC,
      //Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL,
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE,
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE,
      Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND
   };

   auto found = transition_curve_types.find(horizontal_segment->PredefinedType());
   return (found == transition_curve_types.end() ? false : true);
}



template <typename Schema>
bool IsCircularCurve(typename Schema::IfcAlignmentHorizontalSegment* horizontal_segment)
{
    return horizontal_segment->PredefinedType() == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC ? true : false;
}

//Ifc4x3_rc2::IfcAlignmentHorizontal* GetAlignmentHorizontal(Ifc4x3_rc2::IfcAlignment* pAlignment)
//{
//    Ifc4x3_rc2::IfcAlignmentHorizontal* horizontal_alignment = nullptr;
//    auto nested = pAlignment->IsDecomposedBy(); // these are the things that are nested by the alignment
//    for (auto rel_nests : *nested)
//    {
//        ATLASSERT(rel_nests->RelatingObject() == pAlignment);
//        auto related_objects = rel_nests->RelatedObjects();
//        for (auto related_object : *related_objects)
//        {
//            horizontal_alignment = related_object->as<Ifc4x3_rc2::IfcAlignmentHorizontal>();
//            if (horizontal_alignment) break;
//        }
//        if (horizontal_alignment) break;
//    }
//
//    return horizontal_alignment;
//}

template <typename Schema> 
typename Schema::IfcAlignmentHorizontal* GetAlignmentHorizontal(typename Schema::IfcAlignment* pAlignment)
{
    typename Schema::IfcAlignmentHorizontal* horizontal_alignment = nullptr;
    auto nested = pAlignment->IsNestedBy(); // these are the things that are nested by the alignment
    for (auto rel_nests : *nested)
    {
        ATLASSERT(rel_nests->RelatingObject() == pAlignment);
        auto related_objects = rel_nests->RelatedObjects();
        for (auto related_object : *related_objects)
        {
            horizontal_alignment = related_object->as<Schema::IfcAlignmentHorizontal>();
            if (horizontal_alignment) break;
        }
        if (horizontal_alignment) break;
    }

    return horizontal_alignment;
}

//Ifc4x3_rc3::IfcAlignmentHorizontal* GetAlignmentHorizontal(Ifc4x3_rc3::IfcAlignment* pAlignment)
//{
//    return GetAlignmentHorizontal<Ifc4x3_rc3>(pAlignment);
//}
//
//Ifc4x3_rc4::IfcAlignmentHorizontal* GetAlignmentHorizontal(Ifc4x3_rc4::IfcAlignment* pAlignment)
//{
//    return GetAlignmentHorizontal<Ifc4x3_rc4>(pAlignment);
//}

Ifc4x3_tc1::IfcAlignmentHorizontal* GetAlignmentHorizontal(Ifc4x3_tc1::IfcAlignment* pAlignment)
{
   return GetAlignmentHorizontal<Ifc4x3_tc1>(pAlignment);
}

Ifc4x3_add1::IfcAlignmentHorizontal* GetAlignmentHorizontal(Ifc4x3_add1::IfcAlignment* pAlignment)
{
   return GetAlignmentHorizontal<Ifc4x3_add1>(pAlignment);
}


//Ifc4x3_rc2::IfcAlignmentVertical* GetAlignmentVertical(Ifc4x3_rc2::IfcAlignment* pAlignment)
//{
//    Ifc4x3_rc2::IfcAlignmentVertical* vertical_alignment = nullptr;
//    auto nested = pAlignment->IsDecomposedBy(); // these are the things that are nested by the alignment
//    for (auto rel_nests : *nested)
//    {
//        ATLASSERT(rel_nests->RelatingObject() == pAlignment);
//        auto related_objects = rel_nests->RelatedObjects();
//        for (auto related_object : *related_objects)
//        {
//            vertical_alignment = related_object->as<Ifc4x3_rc2::IfcAlignmentVertical>();
//            if (vertical_alignment) break;
//        }
//        if (vertical_alignment) break;
//    }
//
//    return vertical_alignment;
//}

template <typename Schema>
typename Schema::IfcAlignmentVertical* GetAlignmentVertical(typename Schema::IfcAlignment* pAlignment)
{
    Schema::IfcAlignmentVertical* vertical_alignment = nullptr;
    auto nested = pAlignment->IsNestedBy(); // these are the things that are nested by the alignment
    for (auto rel_nests : *nested)
    {
        ATLASSERT(rel_nests->RelatingObject() == pAlignment);
        auto related_objects = rel_nests->RelatedObjects();
        for (auto related_object : *related_objects)
        {
            vertical_alignment = related_object->as<Schema::IfcAlignmentVertical>();
            if (vertical_alignment) break;
        }
        if (vertical_alignment) break;
    }

    return vertical_alignment;
}

//Ifc4x3_rc3::IfcAlignmentVertical* GetAlignmentVertical(Ifc4x3_rc3::IfcAlignment* pAlignment)
//{
//    return GetAlignmentVertical<Ifc4x3_rc3>(pAlignment);
//}
//
//Ifc4x3_rc4::IfcAlignmentVertical* GetAlignmentVertical(Ifc4x3_rc4::IfcAlignment* pAlignment)
//{
//    return GetAlignmentVertical<Ifc4x3_rc4>(pAlignment);
//}

Ifc4x3_tc1::IfcAlignmentVertical* GetAlignmentVertical(Ifc4x3_tc1::IfcAlignment* pAlignment)
{
   return GetAlignmentVertical<Ifc4x3_tc1>(pAlignment);
}

Ifc4x3_add1::IfcAlignmentVertical* GetAlignmentVertical(Ifc4x3_add1::IfcAlignment* pAlignment)
{
   return GetAlignmentVertical<Ifc4x3_add1>(pAlignment);
}

template <typename Schema>
bool CIfcAlignmentConverter::IsValidAlignment_4x3(IfcParse::IfcFile& file,typename Schema::IfcAlignment* pAlignment)
{
    Schema::IfcAlignmentHorizontal* horizontal_alignment = GetAlignmentHorizontal(pAlignment);
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
        auto alignment_segment = (*(related_objects->begin()))->as<Schema::IfcAlignmentSegment>();
        auto horizontal_segment = GetHorizontalAlignmentSegment(alignment_segment);
        return !IsTransitionCurve(horizontal_segment);
    }
    else if (nSegments == 2)
    {
        // our model doesn't support isolated transition segments
        // transition curves must be adjacent to circular curves
        // can't have two transitions adjacent to each other either
        auto alignment_segment1 = (*(related_objects->begin()))->as<Schema::IfcAlignmentSegment>();
        auto horizontal_segment1 = GetHorizontalAlignmentSegment(alignment_segment1);
        bool bIsTransitionCurve1 = IsTransitionCurve(horizontal_segment1);


        auto alignment_segment2 = (*(related_objects->begin() + 1))->as<Schema::IfcAlignmentSegment>();
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
            auto alignment_segment = (*iter)->as<Schema::IfcAlignmentSegment>();
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
                    // transition starts with a radius so a circular curve must preceed this transition curve
                    auto prev_alignment_segment = (*(iter - 1))->as<Schema::IfcAlignmentSegment>();
                    auto prev_horizontal_segment = GetHorizontalAlignmentSegment(prev_alignment_segment);
                    if (!IsCircularCurve<Schema>(prev_horizontal_segment)) return false; // previous is not a circular curve

                    Float64 circular_curve_radius = prev_horizontal_segment->StartRadiusOfCurvature();
                    ATLASSERT(IsEqual(circular_curve_radius, prev_horizontal_segment->EndRadiusOfCurvature()));
                    if (!IsEqual(circular_curve_radius, start_radius)) return false; // common radii must be equal

                    if (::BinarySign(start_radius) != ::BinarySign(circular_curve_radius)) return false; // curves must be same direction
                }

                if (iter != end - 1 && !IsZero(end_radius))
                {
                    // transition ends with a radius so a circular curve most come after this transition curve
                    auto next_alignment_segment = (*(iter + 1))->as<Schema::IfcAlignmentSegment>();
                    auto next_horizontal_segment = GetHorizontalAlignmentSegment(next_alignment_segment);
                    if (!IsCircularCurve<Schema>(next_horizontal_segment)) return false; // next is not a circular curve

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
    void SetProgress(IProgress* pProgress) { m_pProgress = pProgress; }
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
    IProgress* m_pProgress;
};

class ProgressStream : public std::ostream
{
public:
    ProgressStream() : std::ostream(&_progress) {};
    void SetProgress(IProgress* pProgress) { _progress.SetProgress(pProgress); }
private:
    ProgressStringBuf _progress;
};

HRESULT CIfcAlignmentConverter::ConvertToPGSuper(IBroker* pBroker, CString& strFilePath)
{
    USES_CONVERSION;

    std::unique_ptr<IfcParse::IfcFile> pFile = nullptr;

    { // scope the progress window so it closes automatically when we are done with it
        GET_IFACE2(pBroker, IProgress, pProgress);
        CEAFAutoProgress ap(pProgress);

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

   AlignmentData2 alignment_data;
   ProfileData2 profile_data;
   RoadwaySectionData section_data;

   m_Notes.clear();

   auto strSchemaName = pFile->schema()->name();
   bool bResult = false;
   /*if (strSchemaName == std::string("IFC4X1"))
   {
      bResult = ConvertToPGSuper<Ifc4x1>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X2"))
   {
      bResult = ConvertToPGSuper<Ifc4x2>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X3_RC1"))
   {
      bResult = ConvertToPGSuper<Ifc4x3_rc1>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X3_RC2"))
   {
      bResult = ConvertToPGSuper_4x3<Ifc4x3_rc2>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X3_RC3"))
   {
      bResult = ConvertToPGSuper_4x3<Ifc4x3_rc3>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X3_RC4"))
   {
       bResult = ConvertToPGSuper_4x3<Ifc4x3_rc4>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else*/ if (strSchemaName == std::string("IFC4X3_ADD1"))
   {
      bResult = ConvertToPGSuper_4x3<Ifc4x3_add1>(*pFile, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X3_TC1"))
   {
      bResult = ConvertToPGSuper_4x3<Ifc4x3_tc1>(*pFile, &alignment_data, &profile_data, &section_data);
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

   if (bResult)
   {
      GET_IFACE2(pBroker, IEvents, pEvents);
      pEvents->HoldEvents();

      GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
      pRoadwayData->SetAlignmentData2(alignment_data);
      pRoadwayData->SetProfileData2(profile_data);
      pRoadwayData->SetRoadwaySectionData(section_data);

      pEvents->FirePendingEvents();
   }

   return S_OK;
}

template <typename Schema>
bool CIfcAlignmentConverter::ConvertToPGSuper(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData)
{
   InitUnits<Schema>(file);

   auto alignment = GetAlignment<Schema>(file);
   if (alignment == nullptr)
      return false;

   LoadAlignment<Schema>(alignment);
   *pAlignmentData = m_AlignmentData;

   LoadProfile<Schema>(alignment);
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

template <typename Schema>
bool CIfcAlignmentConverter::ConvertToPGSuper_4x3(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData)
{
   InitUnits<Schema>(file);

   auto alignment = GetAlignment_4x3<Schema>(file);
   if (alignment == nullptr)
      return false;

   Float64 alignment_adjustment = LoadAlignment_4x3<Schema>(file,alignment);
   *pAlignmentData = m_AlignmentData;

   LoadProfile_4x3<Schema>(file,alignment, alignment_adjustment);
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

std::vector<std::_tstring> CIfcAlignmentConverter::GetNotes()
{
   return m_Notes;
}

template <typename Schema>
typename Schema::IfcAlignment* CIfcAlignmentConverter::GetAlignment(IfcParse::IfcFile& file)
{
   USES_CONVERSION;

   Schema::IfcAlignment::list::ptr alignments = file.instances_by_type<Schema::IfcAlignment>();
   std::vector<Schema::IfcAlignment*> valid_alignments;

   for (auto alignment : *alignments)
   {
      if (IsValidAlignment<Schema>(alignment))
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
         auto strLabel = (alignment->hasName() ? alignment->Name() : alignment->hasDescription() ? alignment->Description() : "Unnamed");
         os << strLabel << std::endl;
      }
      int result = AfxChoose(_T("Select Alignment"), _T("Select alignment to import"), A2T(os.str().c_str()), 0, TRUE);
      if (result < 0)
         return nullptr; // dialog was canceled
      else
         return valid_alignments[result];
   }


   return nullptr;
}

template <typename Schema>
typename Schema::IfcAlignment* CIfcAlignmentConverter::GetAlignment_4x3(IfcParse::IfcFile& file)
{
   USES_CONVERSION;

   Schema::IfcAlignment::list::ptr alignments = file.instances_by_type<Schema::IfcAlignment>();
   std::vector<Schema::IfcAlignment*> valid_alignments;

   for (auto alignment : *alignments)
   {
      if (IsValidAlignment_4x3<Schema>(file,alignment))
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
      int result = AfxChoose(_T("Select Alignment"), _T("Select alignment to import"), A2T(os.str().c_str()), 0, TRUE);
      if (result < 0)
         return nullptr; // dialog was canceled
      else
         return valid_alignments[result];
   }


   return nullptr;
}

template <typename Schema>
void CIfcAlignmentConverter::LoadAlignment(typename Schema::IfcAlignment* pAlignment)
{
   USES_CONVERSION;
   m_bAlignmentStarted = false; // the alignment data block has not yet been started
   
   m_AlignmentData.Name = A2T(pAlignment->hasName() ? pAlignment->Name().c_str() : pAlignment->hasDescription() ? pAlignment->Description().c_str() : "");

   // initialize the alignment data
   m_AlignmentData.Direction = 0.00;
   m_AlignmentData.xRefPoint = 0.00;
   m_AlignmentData.yRefPoint = 0.00;
   m_AlignmentData.CompoundCurves.clear();

   auto axis = pAlignment->Axis();
   auto curve = axis->as<Schema::IfcAlignmentCurve>();
   auto horizontal = curve->Horizontal();

   Float64 current_station; // station at the start of the current element
   if(horizontal->hasStartDistAlong())
   {
      // as I understand IFC 8.7.3.1, StartDistAlong is the value of the distance along at the start of the alignment... that seems like a starting station
      current_station = WBFL::Units::ConvertToSysUnits(horizontal->StartDistAlong(), *m_pLengthUnit);
   }
   else
   {
      current_station = 0.0;
   }

   // alignment is made up of Line, Spiral, and/or Curve elements
   auto segments = horizontal->Segments();
   auto begin = segments->begin();
   auto iter = begin;
   auto end = segments->end();
   for (; iter != end; iter++)
   {
      auto segment(*iter);
      bool bIsThereANextSegment = ((iter+1) != end);

      auto linear = segment->CurveGeometry()->as<Schema::IfcLineSegment2D>();
      auto transition = segment->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
      auto curve = segment->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();

      Schema::IfcTransitionCurveSegment2D* entrySpiral = nullptr;
      Schema::IfcTransitionCurveSegment2D* exitSpiral = nullptr;

      Float64 end_station = current_station;

      if (linear)
      {
         end_station = OnLine<Schema,Schema::IfcLineSegment2D>(current_station, linear);
      }
      else if (transition)
      {
         // PGSuper can only handle
         // Spiral-Curve
         // Spiral-Curve-Spiral
         // Curve-Spiral
         //
         // Curve-Spiral-Curve, where Spiral is a transition spiral with the start and end radius equal
         // to the curve radii... PGSuper cannot do this case


         entrySpiral = transition;
         if (bIsThereANextSegment)
         {
            // if there is a next segment, check to see if it is a curve
            iter++; // advance to next segment
            curve = (*iter)->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
            if (curve)
            {
               // it's a curve... is there an element that follows the curve?
               bIsThereANextSegment = ((iter+1) != end);
               if (bIsThereANextSegment)
               {
                  // if there is a next segment, see if it is a spiral
                  exitSpiral = (*(iter+1))->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();

                  // if not a spiral, pExitSpiral will be nullptr
                  // this is OK, it just means we have a Spiral-Curve situation
                  if (exitSpiral)
                  {
                     // it is an exit spiral, so advance the iterator
                     iter++;
                     bIsThereANextSegment = ((iter+1) != end);
                  }
               }

               end_station = OnCurve<Schema, Schema::IfcTransitionCurveSegment2D, Schema::IfcCircularArcSegment2D>(current_station, entrySpiral, curve, exitSpiral);
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
      else if (curve)
      {
         // looking for Curve-Spiral case
         if (bIsThereANextSegment)
         {
            // check to see if the next element is a spiral
            exitSpiral = (*(iter + 1))->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
            
            // if not a spiral, pExitSpiral will be nullptr
            // this is OK, it just means we have a Spiral-Curve situation
            if (exitSpiral)
            {
               iter++;
               bIsThereANextSegment = ((iter+1) != end);
            }

            // Check if the next object is a Curve and if
            // the exit spiral and the curve are touching
            if (bIsThereANextSegment)
            {
               auto next_curve = (*(iter))->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
               if (next_curve && exitSpiral)
               {
                  CComPtr<IPoint2d> pntSpiralStart, pntSpiralPI, pntSpiralEnd;
                  GetSpiralPoints<Schema>(exitSpiral, &pntSpiralStart, &pntSpiralPI, &pntSpiralEnd);

                  CComPtr<IPoint2d> pntNextCurveStart, pntNextCurvePI, pntNextCurveEnd, pntNextCurveCenter;
                  GetCurvePoints<Schema,Schema::IfcCircularArcSegment2D>(next_curve, &pntNextCurveStart, &pntNextCurvePI, &pntNextCurveEnd, &pntNextCurveCenter);

                  if (SameLocation(pntSpiralEnd,pntNextCurveStart,m_Precision) == S_FALSE)
                  {
                     IFC_THROW(_T("PGSuper cannot model a transition spiral between two circular curves."));
                  }
               }
            }
         }

         ATLASSERT(entrySpiral == nullptr); // this must be the case
         end_station = OnCurve<Schema, Schema::IfcTransitionCurveSegment2D, Schema::IfcCircularArcSegment2D>(current_station, entrySpiral, curve, exitSpiral);
      }
      else
      {
         ATLASSERT(false);
      }
      current_station = end_station;
   }
}

//Ifc4x3_rc2::IfcAlignmentHorizontalSegment* GetHorizontalAlignmentSegment(Ifc4x3_rc2::IfcAlignmentSegment* alignment_segment)
//{
//   return alignment_segment->GeometricParameters()->as<Ifc4x3_rc2::IfcAlignmentHorizontalSegment>();
//}
//
//Ifc4x3_rc3::IfcAlignmentHorizontalSegment* GetHorizontalAlignmentSegment(Ifc4x3_rc3::IfcAlignmentSegment* alignment_segment)
//{
//   return alignment_segment->DesignParameters()->as<Ifc4x3_rc3::IfcAlignmentHorizontalSegment>();
//}
//
//Ifc4x3_rc4::IfcAlignmentHorizontalSegment* GetHorizontalAlignmentSegment(Ifc4x3_rc4::IfcAlignmentSegment* alignment_segment)
//{
//    return alignment_segment->DesignParameters()->as<Ifc4x3_rc4::IfcAlignmentHorizontalSegment>();
//}

Ifc4x3_tc1::IfcAlignmentHorizontalSegment* GetHorizontalAlignmentSegment(Ifc4x3_tc1::IfcAlignmentSegment* alignment_segment)
{
   return alignment_segment->DesignParameters()->as<Ifc4x3_tc1::IfcAlignmentHorizontalSegment>();
}

Ifc4x3_add1::IfcAlignmentHorizontalSegment* GetHorizontalAlignmentSegment(Ifc4x3_add1::IfcAlignmentSegment* alignment_segment)
{
   return alignment_segment->DesignParameters()->as<Ifc4x3_add1::IfcAlignmentHorizontalSegment>();
}


//Ifc4x3_rc2::IfcAlignmentVerticalSegment* GetVerticalAlignmentSegment(Ifc4x3_rc2::IfcAlignmentSegment* alignment_segment)
//{
//   return alignment_segment->GeometricParameters()->as<Ifc4x3_rc2::IfcAlignmentVerticalSegment>();
//}
//
//Ifc4x3_rc3::IfcAlignmentVerticalSegment* GetVerticalAlignmentSegment(Ifc4x3_rc3::IfcAlignmentSegment* alignment_segment)
//{
//   return alignment_segment->DesignParameters()->as<Ifc4x3_rc3::IfcAlignmentVerticalSegment>();
//}
//
//Ifc4x3_rc4::IfcAlignmentVerticalSegment* GetVerticalAlignmentSegment(Ifc4x3_rc4::IfcAlignmentSegment* alignment_segment)
//{
//    return alignment_segment->DesignParameters()->as<Ifc4x3_rc4::IfcAlignmentVerticalSegment>();
//}

Ifc4x3_tc1::IfcAlignmentVerticalSegment* GetVerticalAlignmentSegment(Ifc4x3_tc1::IfcAlignmentSegment* alignment_segment)
{
   return alignment_segment->DesignParameters()->as<Ifc4x3_tc1::IfcAlignmentVerticalSegment>();
}

Ifc4x3_add1::IfcAlignmentVerticalSegment* GetVerticalAlignmentSegment(Ifc4x3_add1::IfcAlignmentSegment* alignment_segment)
{
   return alignment_segment->DesignParameters()->as<Ifc4x3_add1::IfcAlignmentVerticalSegment>();
}


//Float64 CIfcAlignmentConverter::GetStartDistAlong(Ifc4x3_rc2::IfcAlignmentHorizontal* pHorizontal)
//{
//    Float64 station = pHorizontal->hasStartDistAlong() ? pHorizontal->StartDistAlong() : 0.0;
//    return ::ConvertToSysUnits(station, *m_pLengthUnit);
//}
//
//Float64 CIfcAlignmentConverter::GetStartDistAlong(Ifc4x3_rc3::IfcAlignmentHorizontal* pHorizontal)
//{
//    //Float64 station = pHorizontal->hasStartDistAlong() ? pHorizontal->StartDistAlong() : 0.0;
//    //return ::ConvertToSysUnits(station, *m_pLengthUnit);
//    return 0.0; // not a property of IfcAlignmentHorizontal in rc3
//}
//
//Float64 CIfcAlignmentConverter::GetStartDistAlong(Ifc4x3_rc4::IfcAlignmentHorizontal* pHorizontal)
//{
//    return 0.0; // not a property of IfcAlignmentHorizontal in rc4
//}

Float64 CIfcAlignmentConverter::GetStartDistAlong(Ifc4x3_tc1::IfcAlignmentHorizontal* pHorizontal)
{
   return 0.0; // not a property of IfcAlignmentHorizontal in tc1
}

Float64 CIfcAlignmentConverter::GetStartDistAlong(Ifc4x3_add1::IfcAlignmentHorizontal* pHorizontal)
{
   return 0.0; // not a property of IfcAlignmentHorizontal in add1
}


template <typename Schema>
void CIfcAlignmentConverter::GetStations(typename Schema::IfcAlignment* pAlignment, std::vector<std::pair<Float64, Float64>>& vStations, std::vector<std::tuple<Float64, Float64, Float64>>& vStationEquations)
{
   auto nested = pAlignment->IsNestedBy();
   for (auto rel_nests : *nested)
   {
      ATLASSERT(rel_nests->RelatingObject() == pAlignment);
      auto related_objects = rel_nests->RelatedObjects();
      for (auto related_object : *related_objects)
      {
         auto referent = related_object->as<Schema::IfcReferent>();
         if (referent && referent->PredefinedType() && *(referent->PredefinedType()) == Schema::IfcReferentTypeEnum::IfcReferentType_STATION)
         {
            Float64 distance_along = 0;
            if (referent->ObjectPlacement())
            {
               auto object_placement = referent->ObjectPlacement();
               auto linear_placement = object_placement->as<Schema::IfcLinearPlacement>();
               if (linear_placement)
               {
                  // get the distance along the curve for the placement of the referent
                  auto axis2placementlinear = linear_placement->RelativePlacement();
                  auto location = axis2placementlinear->Location();
                  auto point_by_distance_expression = location->as<Schema::IfcPointByDistanceExpression>();
                  if (point_by_distance_expression)
                  {
                     distance_along = *(point_by_distance_expression->DistanceAlong()->as<Schema::IfcNonNegativeLengthMeasure>());
                  }
               }
            }

            auto rel_defines_by_properties = referent->IsDefinedBy();
            for (auto rel_defines_property : *rel_defines_by_properties)
            {
               auto property_set = rel_defines_property->RelatingPropertyDefinition()->as<Schema::IfcPropertySet>();
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
                        auto single_value_property = prop->as<Schema::IfcPropertySingleValue>();
                        if (single_value_property->NominalValue())
                        {
                           bHasStation = true;
                           station = *(single_value_property->NominalValue()->as<Schema::IfcLengthMeasure>());
                        }
                     }
                     else if (prop->Name() == "IncomingStation")
                     {
                        auto single_value_property = prop->as<Schema::IfcPropertySingleValue>();
                        if (single_value_property->NominalValue())
                        {
                           bHasIncomingStation = true;
                           incoming_station = *(single_value_property->NominalValue()->as<Schema::IfcLengthMeasure>());
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

template <typename Schema>
Float64 CIfcAlignmentConverter::LoadAlignment_4x3(IfcParse::IfcFile& file, typename Schema::IfcAlignment* pAlignment)
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
    GetStations<Schema>(pAlignment, vStations, vStationEquations);

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

    Schema::IfcAlignmentHorizontal* horizontal_alignment = GetAlignmentHorizontal(pAlignment);
    ATLASSERT(horizontal_alignment); // should have found one

    Float64 current_station = GetStartDistAlong(horizontal_alignment); // not part of rc4, so returns zero
    current_station += station_adjustment;

    // alignment is made up of Line, Spiral, and/or Curve elements
    auto nested = horizontal_alignment->IsNestedBy();
    auto related_objects = (*nested->begin())->RelatedObjects();
    auto begin = related_objects->begin();
    auto iter = begin;
    auto end = related_objects->end();
    for (; iter != end; iter++)
    {
        auto alignment_segment = (*iter)->as<Schema::IfcAlignmentSegment>();
        bool bIsThereANextSegment = ((iter + 1) != end);

        auto horizontal_alignment_segment = GetHorizontalAlignmentSegment(alignment_segment);

        auto predefined_type = horizontal_alignment_segment->PredefinedType();

        Schema::IfcAlignmentHorizontalSegment* entrySpiral = nullptr;
        Schema::IfcAlignmentHorizontalSegment* exitSpiral = nullptr;

        Float64 end_station = current_station;

        if (predefined_type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE)
        {
            end_station = OnLine<Schema, Schema::IfcAlignmentHorizontalSegment>(current_station, horizontal_alignment_segment);
        }
        else if (predefined_type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
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
                auto curve = GetHorizontalAlignmentSegment((*(iter))->as<Schema::IfcAlignmentSegment>());
                if (curve->PredefinedType() != Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
                    curve = nullptr;

                if (curve)
                {
                    // it's a curve... is there an element that follows the curve?
                    bIsThereANextSegment = ((iter + 1) != end);
                    if (bIsThereANextSegment)
                    {
                        // if there is a next segment, see if it is a spiral
                        exitSpiral = GetHorizontalAlignmentSegment((*(iter + 1))->as<Schema::IfcAlignmentSegment>());
                        if (exitSpiral->PredefinedType() != Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
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

                    end_station = OnCurve_4x3<Schema>(current_station, entrySpiral, curve, exitSpiral);
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
        else if (predefined_type == Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
        {
            // looking for Curve-Spiral case
            if (bIsThereANextSegment)
            {
                // check to see if the next element is a spiral
                exitSpiral = GetHorizontalAlignmentSegment((*(iter + 1))->as<Schema::IfcAlignmentSegment>());
                if (exitSpiral->PredefinedType() != Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID)
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
                    auto next_curve = GetHorizontalAlignmentSegment((*(iter))->as<Schema::IfcAlignmentSegment>());
                    if (next_curve->PredefinedType() != Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC)
                        next_curve = nullptr;

                    if (next_curve && exitSpiral)
                    {
                        CComPtr<IPoint2d> pntSpiralStart, pntSpiralPI, pntSpiralEnd;
                        GetSpiralPoints_4x3<Schema>(exitSpiral, &pntSpiralStart, &pntSpiralPI, &pntSpiralEnd);

                        CComPtr<IPoint2d> pntNextCurveStart, pntNextCurvePI, pntNextCurveEnd, pntNextCurveCenter;
                        GetCurvePoints_4x3<Schema>(next_curve, &pntNextCurveStart, &pntNextCurvePI, &pntNextCurveEnd, &pntNextCurveCenter);

                        if (SameLocation(pntSpiralEnd, pntNextCurveStart, m_Precision) == S_FALSE)
                        {
                            IFC_THROW(_T("PGSuper cannot model a transition spiral between two circular curves."));
                        }
                    }
                }
            }

            ATLASSERT(entrySpiral == nullptr); // this must be the case
            end_station = OnCurve_4x3<Schema>(current_station, entrySpiral, horizontal_alignment_segment, exitSpiral);
        }
        else
        {
            ATLASSERT(false);
        }
        current_station = end_station;
    }

    return station_adjustment;
}

template <typename Schema>
void CIfcAlignmentConverter::LoadProfile(typename Schema::IfcAlignment* pAlignment)
{
   m_ProfileState = PROFILE_NOT_STARTED;
   m_ProfileData.VertCurves.clear();


   auto axis = pAlignment->Axis();
   auto curve = axis->as<Schema::IfcAlignmentCurve>();
   auto vertical = curve->hasVertical() ? curve->Vertical() : nullptr;

   Float64 start_station = curve->Horizontal()->hasStartDistAlong() ? curve->Horizontal()->StartDistAlong() : 0.0;
   start_station = WBFL::Units::ConvertToSysUnits(start_station, *m_pLengthUnit);

   auto segments = vertical ? vertical->Segments() : nullptr;

  if (vertical == nullptr || segments->size() == 0)
   {
      // the profile geometry list is empty so assume a flat grade
      m_Notes.push_back(std::_tstring(_T("A profile was not found or the profile does not contain segments. Assuming a default profile.")));
      m_ProfileData.Station = start_station;
      m_ProfileData.Elevation = 0.0;
      m_ProfileData.Grade = 0.0;
      return;
   }

   for (auto segment : *segments)
   {
      auto linear_segment = segment->as<Schema::IfcAlignment2DVerSegLine>();
      auto parabolic_arc = segment->as<Schema::IfcAlignment2DVerSegParabolicArc>();
      auto circular_arc = segment->as<Schema::IfcAlignment2DVerSegCircularArc>();

      if (linear_segment)
      {
         LinearSegment<Schema, Schema::IfcAlignment2DVerSegLine>(start_station,linear_segment);
      }
      else if (parabolic_arc)
      {
         ParabolicSegment<Schema, Schema::IfcAlignment2DVerSegParabolicArc>(start_station, parabolic_arc);
      }
      else if (circular_arc)
      {
#pragma Reminder("WORKING HERE - Need to deal with circular arcs") // treat it as a parabola for now
         parabolic_arc = (Schema::IfcAlignment2DVerSegParabolicArc*)circular_arc;
         ParabolicSegment<Schema, Schema::IfcAlignment2DVerSegParabolicArc>(start_station,parabolic_arc);
         // TODO: provide a better exception
         //IFC_THROW(_T("Circular curve element was found in the profile definition. PGSuper does not support circular vertical curves."));
      }
      else
      {
         ATLASSERT(false); // is there a new type ???
                           // TODO: provide a better exception
         IFC_THROW(_T("An unknown profile element was encountered"));
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

template <typename Schema>
void CIfcAlignmentConverter::LoadProfile_4x3(IfcParse::IfcFile& file, typename Schema::IfcAlignment* pAlignment, Float64 stationAdjustment)
{
   m_ProfileState = PROFILE_NOT_STARTED;
   m_ProfileData.Station = 0;
   m_ProfileData.Elevation = 0;
   m_ProfileData.Grade = 0;
   m_ProfileData.VertCurves.clear();


   Schema::IfcAlignmentHorizontal* horizontal_alignment = GetAlignmentHorizontal(pAlignment);
   Schema::IfcAlignmentVertical* vertical_alignment = GetAlignmentVertical(pAlignment);

   Float64 start_station = GetStartDistAlong(horizontal_alignment);
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
          auto alignment_segment = (*iter)->as<Schema::IfcAlignmentSegment>();
          auto vertical_alignment_segment = GetVerticalAlignmentSegment(alignment_segment);
          auto predefined_type = vertical_alignment_segment->PredefinedType();

          if (predefined_type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT)
          {
              LinearSegment_4x3<Schema>(start_station,vertical_alignment_segment);
          }
          else if (predefined_type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC)
          {
              ParabolicSegment_4x3<Schema>(start_station,vertical_alignment_segment);
          }
          else if (predefined_type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CIRCULARARC)
          {
#pragma Reminder("WORKING HERE - Need to deal with vertical circular arcs") // treat it as a parabola for now
             ParabolicSegment_4x3<Schema>(start_station, vertical_alignment_segment);
          }
          else if (predefined_type == Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CLOTHOID)
          {
#pragma Reminder("WORKING HERE - Need to deal with vertical clothoid arcs") // treat it as a parabola for now
             ParabolicSegment_4x3<Schema>(start_station, vertical_alignment_segment);
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

template <typename Schema,typename Segment>
Float64 CIfcAlignmentConverter::OnLine(Float64 startStation, typename Segment* pLine)
{
   Float64 sx, sy;
   GetPoint<Schema>(pLine->StartPoint(), &sx, &sy);

   Float64 length = WBFL::Units::ConvertToSysUnits(pLine->SegmentLength(),*m_pLengthUnit);
   Float64 startDirection = WBFL::Units::ConvertToSysUnits(pLine->StartDirection(), *m_pAngleUnit);;
   return OnLine(sx, sy, startStation, startDirection, length);
}

Float64 CIfcAlignmentConverter::OnLine(Float64 sx,Float64 sy,Float64 startStation,Float64 startDirection, Float64 length)
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

template <typename Schema,typename SpiralType, typename CurveType>
Float64 CIfcAlignmentConverter::OnCurve(Float64 startStation, typename SpiralType* pEntrySpiral, typename CurveType* pCurve, typename SpiralType* pExitSpiral)
{
   ATLASSERT(pCurve != nullptr);

   Float64 radius = WBFL::Units::ConvertToSysUnits(pCurve->Radius(), *m_pLengthUnit);

   // Get all the construction points
   CComPtr<IPoint2d> pntEntryStart, pntEntryPI, pntEntryEnd;
   CComPtr<IPoint2d> pntExitStart, pntExitPI, pntExitEnd;
   CComPtr<IPoint2d> pntCurveStart, pntCurvePI, pntCurveEnd, pntCurveCenter;

   Float64 entry_spiral_length = 0;
   Float64 exit_spiral_length = 0;
   if (pEntrySpiral)
   {
      GetSpiralPoints<Schema>(pEntrySpiral, &pntEntryStart, &pntEntryPI, &pntEntryEnd);
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
         if(pEntrySpiral->hasStartRadius() && !IsZero(pEntrySpiral->StartRadius()))
         {
            m_Notes.push_back(std::_tstring(_T("Start radius of entry spiral taken to be infinite")));
         }
         CheckSpiralType<Schema,SpiralType>(pEntrySpiral);
      }
   }

   GetCurvePoints<Schema,CurveType>(pCurve, &pntCurveStart, &pntCurvePI, &pntCurveEnd, &pntCurveCenter);

   if (pntCurveStart == nullptr || pntCurvePI == nullptr || pntCurveEnd == nullptr || pntCurveCenter == nullptr)
   {
      m_Notes.push_back(std::_tstring(_T("Zero radius curve could not be constructed.")));
      return startStation;
   }

   if (pExitSpiral)
   {
      GetSpiralPoints<Schema>(pExitSpiral, &pntExitStart, &pntExitPI, &pntExitEnd);
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
         if (pExitSpiral->hasEndRadius() && !IsZero(pExitSpiral->EndRadius()))
         {
            m_Notes.push_back(std::_tstring(_T("End radius of exit spiral taken to be infinite")));
         }
         CheckSpiralType<Schema,SpiralType>(pExitSpiral);
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

      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pEntrySpiral->EndRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("End radius of the entry spiral does not match the radius of the circular curve. The entry spiral end radius will be ignored.")));
      }

      if (SameLocation(pntEntryEnd,pntCurveStart,m_Precision) == S_FALSE)
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

      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pExitSpiral->StartRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("Start radius of the exit spiral does not match the radius of the circular curve. The exit spiral start radius will be ignored.")));
      }

      if (SameLocation(pntCurveEnd,pntExitStart,m_Precision) == S_FALSE)
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


      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pEntrySpiral->EndRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("End radius of the entry spiral does not match the radius of the circular curve. The entry spiral end radius will be ignored.")));
      }

      if (!IsEqual(WBFL::Units::ConvertToSysUnits(pExitSpiral->StartRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("Start radius of the exit spiral does not match the radius of the circular curve. The exit spiral start radius will be ignored.")));
      }

      if (SameLocation(pntEntryEnd,pntCurveStart,m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The end of the entry spiral does not coincide with the start of the circular curve. The entry spiral has been adjusted.")));
      }

      if (SameLocation(pntCurveEnd,pntExitStart,m_Precision) == S_FALSE)
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
   hc->putref_PBT(pntStart);
   hc->putref_PI(pntPI);
   hc->putref_PFT(pntEnd);
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

template <typename Schema>
Float64 CIfcAlignmentConverter::OnCurve_4x3(Float64 startStation, typename Schema::IfcAlignmentHorizontalSegment* pEntrySpiral, typename Schema::IfcAlignmentHorizontalSegment* pCurve, typename Schema::IfcAlignmentHorizontalSegment* pExitSpiral)
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
      GetSpiralPoints_4x3<Schema>(pEntrySpiral, &pntEntryStart, &pntEntryPI, &pntEntryEnd);
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
         CheckSpiralType_4x3(pEntrySpiral);
      }
   }

   GetCurvePoints_4x3<Schema>(pCurve, &pntCurveStart, &pntCurvePI, &pntCurveEnd, &pntCurveCenter);

   if (pntCurveStart == nullptr || pntCurvePI == nullptr || pntCurveEnd == nullptr || pntCurveCenter == nullptr)
   {
      m_Notes.push_back(std::_tstring(_T("Zero radius curve could not be constructed.")));
      return startStation;
   }

   if (pExitSpiral)
   {
      GetSpiralPoints_4x3<Schema>(pExitSpiral, &pntExitStart, &pntExitPI, &pntExitEnd);
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
         CheckSpiralType_4x3(pExitSpiral);
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
   hc->putref_PBT(pntStart);
   hc->putref_PI(pntPI);
   hc->putref_PFT(pntEnd);
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

template <typename Schema, typename CurveType>
void CIfcAlignmentConverter::GetCurvePoints(typename CurveType* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter)
{
   auto pStart = pCurve->StartPoint();
   auto bkTangentBrg = WBFL::Units::ConvertToSysUnits(pCurve->StartDirection(), *m_pAngleUnit);
   auto L = WBFL::Units::ConvertToSysUnits(pCurve->SegmentLength(),*m_pLengthUnit);
   auto R = WBFL::Units::ConvertToSysUnits(pCurve->Radius(),*m_pLengthUnit);

   Float64 delta = L / R;
   Float64 T = R*tan(delta/2);

   Float64 sx, sy;
   GetPoint<Schema>(pStart, &sx, &sy);
   CComPtr<IPoint2d> pntStart;
   pntStart.CoCreateInstance(CLSID_Point2d);
   pntStart->Move(sx, sy);
   pntStart.CopyTo(ppStart);

   CComPtr<ILocate2> locate;
   m_CogoEngine->get_Locate(&locate);

   locate->ByDistDir(*ppStart, T, CComVariant(bkTangentBrg), 0.0, ppPI);

   Float64 fwdTangentBrg = bkTangentBrg + (pCurve->IsCCW() ? 1 : -1)*delta;

   locate->ByDistDir(*ppPI, T, CComVariant(fwdTangentBrg), 0.0, ppEnd);

   locate->ByDistDir(*ppStart, R, CComVariant(bkTangentBrg + (pCurve->IsCCW() ? 1 : -1)*PI_OVER_2), 0.0, ppCenter);
}

template <typename Schema>
void CIfcAlignmentConverter::GetCurvePoints_4x3(typename Schema::IfcAlignmentHorizontalSegment* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter)
{
   auto pStart = pCurve->StartPoint();
   auto bkTangentBrg = WBFL::Units::ConvertToSysUnits(pCurve->StartDirection(), *m_pAngleUnit);
   auto L = WBFL::Units::ConvertToSysUnits(pCurve->SegmentLength(), *m_pLengthUnit);
   auto R = WBFL::Units::ConvertToSysUnits(pCurve->StartRadiusOfCurvature(), *m_pLengthUnit);
   bool bIsCCW = (R < 0 ? true : false);

   Float64 delta = fabs(L / R);
   Float64 T = R*tan(delta / 2);

   Float64 sx, sy;
   GetPoint<Schema>(pStart, &sx, &sy);
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

template <typename Schema,typename SpiralType>
void CIfcAlignmentConverter::GetSpiralPoints(typename SpiralType* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd)
{
   auto pStart = pSpiral->StartPoint();
   auto bkTangentBrg = WBFL::Units::ConvertToSysUnits(pSpiral->StartDirection(), *m_pAngleUnit);
   auto L = WBFL::Units::ConvertToSysUnits(pSpiral->SegmentLength(),*m_pLengthUnit);
   auto R = WBFL::Units::ConvertToSysUnits((pSpiral->hasStartRadius() ? pSpiral->StartRadius() : pSpiral->EndRadius()),*m_pLengthUnit);
   bool bIsCCW = (pSpiral->hasStartRadius() ? pSpiral->IsStartRadiusCCW() : pSpiral->IsEndRadiusCCW());

   Float64 sx, sy;
   GetPoint<Schema>(pStart, &sx, &sy);
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

   if (pSpiral->hasStartRadius())
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

template <typename Schema>
void CIfcAlignmentConverter::GetSpiralPoints_4x3(typename Schema::IfcAlignmentHorizontalSegment* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd)
{
   auto pStart = pSpiral->StartPoint();
   auto bkTangentBrg = WBFL::Units::ConvertToSysUnits(pSpiral->StartDirection(), *m_pAngleUnit);
   auto L = WBFL::Units::ConvertToSysUnits(pSpiral->SegmentLength(), *m_pLengthUnit);
   auto Rstart = WBFL::Units::ConvertToSysUnits(pSpiral->StartRadiusOfCurvature(), *m_pLengthUnit);
   auto Rend = WBFL::Units::ConvertToSysUnits(pSpiral->EndRadiusOfCurvature(), *m_pLengthUnit);
   auto R = IsZero(Rstart) ? Rend : Rstart; // zero means infinite radius
   bool bIsCCW = (R < 0 ? true : false);

   Float64 sx, sy;
   GetPoint<Schema>(pStart, &sx, &sy);
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

template <typename Schema,typename LineSegmentType>
void CIfcAlignmentConverter::LinearSegment(Float64 startStation,typename LineSegmentType* pLinearSegment)
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

template <typename Schema>
void CIfcAlignmentConverter::LinearSegment_4x3(Float64 startStation, typename Schema::IfcAlignmentVerticalSegment* pLinearSegment)
{
   LinearSegment<Schema::IfcAlignmentVerticalSegment>(startStation, pLinearSegment);
}

template <typename Schema, typename ParabolicSegmentType>
void CIfcAlignmentConverter::ParabolicSegment(Float64 startStation, typename ParabolicSegmentType* pParaCurve)
{
   // finish any open profile element
   Float64 start_gradient = pParaCurve->StartGradient();
   Float64 start_dist = WBFL::Units::ConvertToSysUnits(pParaCurve->StartDistAlong(), *m_pLengthUnit);
   Float64 start_height = WBFL::Units::ConvertToSysUnits(pParaCurve->StartHeight(), *m_pLengthUnit);
   Float64 length = WBFL::Units::ConvertToSysUnits(pParaCurve->HorizontalLength(), *m_pLengthUnit);
   Float64 R = WBFL::Units::ConvertToSysUnits(pParaCurve->ParabolaConstant(), *m_pLengthUnit);
   if (pParaCurve->IsConvex())
      R *= -1;

   Float64 exit_gradient = length/R + start_gradient;

   if (m_ProfileState == PROFILE_NOT_STARTED)
   {
      m_ProfileData.Station = startStation + start_dist;
      m_ProfileData.Elevation = start_height;
      m_ProfileData.Grade = start_gradient;
      m_ProfileState = PROFILE_ESTABLISHED;
   }

   // add this vertical curve
   VertCurveData vcData;
   vcData.PVIStation = startStation + start_dist + length/2;
   vcData.L1 = length;
   vcData.L2 = 0;
   vcData.ExitGrade = exit_gradient;

   m_ProfileData.VertCurves.push_back(vcData);
}

template <typename Schema>
void CIfcAlignmentConverter::ParabolicSegment_4x3(Float64 startStation, typename Schema::IfcAlignmentVerticalSegment* pParaCurve)
{
   // finish any open profile element
   Float64 start_gradient = pParaCurve->StartGradient();
   Float64 start_dist = WBFL::Units::ConvertToSysUnits(pParaCurve->StartDistAlong(), *m_pLengthUnit);
   Float64 start_height = WBFL::Units::ConvertToSysUnits(pParaCurve->StartHeight(), *m_pLengthUnit);
   Float64 length = WBFL::Units::ConvertToSysUnits(pParaCurve->HorizontalLength(), *m_pLengthUnit);
   Float64 R = WBFL::Units::ConvertToSysUnits(pParaCurve->RadiusOfCurvature() ? *(pParaCurve->RadiusOfCurvature()) : 0.0, *m_pLengthUnit);

   Float64 exit_gradient = pParaCurve->EndGradient();

   // determine of the parabola is convex
   // if the second derivative is > 0 it is convex
   // y = (g2 - g1)/(2L)X^2 + g1*X + startElev
   // y' = (g2 - g1)/L * X + g1
   // y'' = (g2 - g1)/L
   if (0 < (exit_gradient - start_gradient)/length) // is convex
      R *= -1; // convex

   ATLASSERT(IsEqual(exit_gradient, start_gradient + length / R));

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

template <typename Schema, typename SpiralType>
void CIfcAlignmentConverter::CheckSpiralType(typename SpiralType* pSpiral)
{
   switch (pSpiral->TransitionCurveType())
   {
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_BIQUADRATICPARABOLA:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_BLOSSCURVE:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_COSINECURVE:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_CUBICPARABOLA:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_SINECURVE:
      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid")));
      break;

   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_CLOTHOIDCURVE:
      // this is ok... it is what we were expecting
      break;

   default:
      ATLASSERT(false); // is there a new spiral type???
      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid")));
      break;
   }
}

//void CIfcAlignmentConverter::CheckSpiralType_4x3(Ifc4x3_rc2::IfcAlignmentHorizontalSegment* pSpiral)
//{
//   switch (pSpiral->PredefinedType())
//   {
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE:
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE:
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL:
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BIQUADRATICPARABOLA:
//   //case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC:
//   //case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE:
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE:
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND:
//      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid.")));
//      break;
//
//   case Ifc4x3_rc2::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID:
//      // this is ok... we were expecting clothoid
//      break;
//
//   default:
//      ATLASSERT(false); // is there a new spiral type???
//      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid.")));
//      break;
//   }
//}
//
//void CIfcAlignmentConverter::CheckSpiralType_4x3(Ifc4x3_rc3::IfcAlignmentHorizontalSegment* pSpiral)
//{
//   switch (pSpiral->PredefinedType())
//   {
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE:
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE:
//   //case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL:
//   //case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BIQUADRATICPARABOLA:
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC:
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE:
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE:
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND:
//      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid.")));
//      break;
//
//   case Ifc4x3_rc3::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID:
//      // this is ok... we were expecting clothoid
//      break;
//
//   default:
//      ATLASSERT(false); // is there a new spiral type???
//      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid.")));
//      break;
//   }
//}
//
//void CIfcAlignmentConverter::CheckSpiralType_4x3(Ifc4x3_rc4::IfcAlignmentHorizontalSegment* pSpiral)
//{
//    switch (pSpiral->PredefinedType())
//    {
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE:
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE:
//    //case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL:
//    //case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BIQUADRATICPARABOLA:
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC:
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE:
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE:
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND:
//        m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid.")));
//        break;
//
//    case Ifc4x3_rc4::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID:
//        // this is ok... we were expecting clothoid
//        break;
//
//    default:
//        ATLASSERT(false); // is there a new spiral type???
//        m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid.")));
//        break;
//    }
//}

void CIfcAlignmentConverter::CheckSpiralType_4x3(Ifc4x3_tc1::IfcAlignmentHorizontalSegment* pSpiral)
{
   switch (pSpiral->PredefinedType())
   {
   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE:
   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE:
      //case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL:
      //case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BIQUADRATICPARABOLA:
   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC:
   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE:
   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE:
   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND:
      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid.")));
      break;

   case Ifc4x3_tc1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID:
      // this is ok... we were expecting clothoid
      break;

   default:
      ATLASSERT(false); // is there a new spiral type???
      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid.")));
      break;
   }
}

void CIfcAlignmentConverter::CheckSpiralType_4x3(Ifc4x3_add1::IfcAlignmentHorizontalSegment* pSpiral)
{
   switch (pSpiral->PredefinedType())
   {
   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BLOSSCURVE:
   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_COSINECURVE:
      //case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBICSPIRAL:
      //case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_BIQUADRATICPARABOLA:
   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CUBIC:
   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_HELMERTCURVE:
   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_SINECURVE:
   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_VIENNESEBEND:
      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid.")));
      break;

   case Ifc4x3_add1::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID:
      // this is ok... we were expecting clothoid
      break;

   default:
      ATLASSERT(false); // is there a new spiral type???
      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid.")));
      break;
   }
}

