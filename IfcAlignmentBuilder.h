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

#include "IfcModelBuilder.h"

#include <IFace/Tools.h>
#include <IFace\Project.h>
#include <IFace\Alignment.h>
#include <EAF\EAFDisplayUnits.h>


// creates geometry and business logic segments for horizontal alignment tangent runs
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_tangent(typename Schema::IfcCartesianPoint* p, double dir, double length, const CIfcModelBuilderOptions& options)
{
   // business logic
   auto design_parameters = new Schema::IfcAlignmentHorizontalSegment(
      boost::none, boost::none, p, dir, 0.0, 0.0, length, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);

   auto alignment_segment = new Schema::IfcAlignmentSegment(
      IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, design_parameters);

   // geometry
   typename Schema::IfcCurveSegment* curve_segment = nullptr;
   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      curve_segment = mapAlignmentSegment(alignment_segment).first;
   }

   return { curve_segment, alignment_segment };
}

// creates geometry and business logic segments for horizontal alignment horizonal curves
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_hcurve(typename Schema::IfcCartesianPoint* pc, double dir, double radius, double lc, const CIfcModelBuilderOptions& options)
{
   // business logic
   auto design_parameters = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, pc, dir, radius, radius, lc, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC);
   auto alignment_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, design_parameters);

   // geometry
   typename Schema::IfcCurveSegment* curve_segment = nullptr;
   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      curve_segment = mapAlignmentSegment(alignment_segment).first;
   }

   return { curve_segment, alignment_segment };
}

// creates geometry and business logic segments for horizontal alignment entry clothoid transition curve
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_entry_spiral(typename Schema::IfcCartesianPoint* pc, double dir, double radius, double ls, const CIfcModelBuilderOptions& options)
{
   // business logic
   auto design_parameters = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, pc, dir, 0.0, radius, ls, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
   auto alignment_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, design_parameters);

   // geometry
   typename Schema::IfcCurveSegment* curve_segment = nullptr;
   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      curve_segment = mapAlignmentSegment(alignment_segment).first;
   }

   return { curve_segment, alignment_segment };
}

// creates geometry and business logic segments for horizontal alignment exit clothoid transition curve
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_exit_spiral(typename Schema::IfcCartesianPoint* pc, double dir, double radius, double ls, const CIfcModelBuilderOptions& options)
{
   // business logic
   auto design_parameters = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, pc, dir, radius, 0.0, ls, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
   auto alignment_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, design_parameters);

   // geometry
   typename Schema::IfcCurveSegment* curve_segment = nullptr;
   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      curve_segment = mapAlignmentSegment(alignment_segment).first;
   }

   return { curve_segment, alignment_segment };
}

// creates geometry and business logic segments for vertical profile gradient runs
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_gradient(typename Schema::IfcCartesianPoint* p, double slope, double length, const CIfcModelBuilderOptions& options)
{
   CHECK(0 <= length);

   // business logic
   auto design_parameters = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, p->Coordinates()[0], length, p->Coordinates()[1], slope, slope, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
   auto alignment_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, design_parameters);

   // geometry
   typename Schema::IfcCurveSegment* curve_segment = nullptr;
   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      curve_segment = mapAlignmentSegment(alignment_segment).first;
   }

   return { curve_segment, alignment_segment };
}

// creates geometry and business logic segments for vertical profile parabolic vertical curves
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_vcurve(typename Schema::IfcCartesianPoint* p, double start_slope, double end_slope, double length, const CIfcModelBuilderOptions& options)
{
   CHECK(0 < length);

   if (IsEqual(start_slope, end_slope))
   {
      // this is actually a gradient line
      return create_gradient<Schema>(p, start_slope, length, options);
   }

   // business logic
   double R = length / (end_slope - start_slope);
   auto design_parameters = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, p->Coordinates()[0], length, p->Coordinates()[1], start_slope, end_slope, R, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
   auto alignment_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, design_parameters);

   // geometry
   typename Schema::IfcCurveSegment* curve_segment = nullptr;
   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      curve_segment = mapAlignmentSegment(alignment_segment).first;
   }

   return { curve_segment, alignment_segment };
}

CComPtr<IPoint2d> GetAlignmentStartPoint(std::shared_ptr<WBFL::EAF::Broker> pBroker,Float64* pStation,Float64* pElevation,Float64* pGrade)
{
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, pStation, pElevation,pGrade, &startPoint);

   return startPoint;
}

template <typename Schema>
void CreateHorizontalAlignment(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options, typename Schema::IfcAlignmentHorizontal** phorizontal_alignment, typename Schema::IfcRelNests** pnests_horizontal_segments, typename Schema::IfcCompositeCurve** phorizontal_geometry_base_curve)
{
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr alignment_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
   typename aggregate_of<typename Schema::IfcSegment>::ptr curve_segments(new aggregate_of<typename Schema::IfcSegment>());

   Float64 startStation, startElevation, startGrade;
   auto startPoint = GetAlignmentStartPoint(pBroker,&startStation,&startElevation,&startGrade);

   // create the start point
   auto ifc_start_point = ConvertPoint<Schema>(startPoint);

   // loop over all the horizontal curves
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> prevPoint = startPoint;
   auto ifc_prev_point = ifc_start_point;
   IndexType nHCurves = pAlignment->GetCurveCount();
   for (IndexType i = 0; i < nHCurves; i++)
   {
      CComPtr<ICompoundCurve> curve;
      pAlignment->GetCurve(i, pgsTypes::pcGlobal, &curve);

      CComPtr<IPoint2d> pntTS;
      curve->get_TS(&pntTS);

      // create a line segment from end of previous alignment segment to the start of this curve
      if (pntTS->SameLocation(prevPoint) == S_FALSE)
      {
         GET_IFACE2(pBroker, IGeometry, pGeometry);
         Float64 dist;
         CComPtr<IDirection> direction;
         pGeometry->Inverse(prevPoint, pntTS, &dist, &direction);
         Float64 angle;
         direction->get_Value(&angle);

         auto [geometry_segment, business_segment] = create_tangent<Schema>(ifc_prev_point, angle, dist, options);
         file.addEntity(business_segment);
         alignment_segments->push(business_segment);
         if (geometry_segment)
         {
            file.addEntity(geometry_segment);
            curve_segments->push(geometry_segment);
         }
      }

      // create this horizontal curve

      std::array<Float64, 2> Lspiral;
      curve->get_SpiralLength(spEntry, &Lspiral[spEntry]);
      curve->get_SpiralLength(spExit, &Lspiral[spExit]);

      Float64 Lc;
      curve->get_CurveLength(&Lc);

      Float64 R;
      curve->get_Radius(&R);

      ATLASSERT(0 < R); // need to deal with zero radius curves, which gives us an angle point in the alignment

      CurveDirectionType curve_direction;
      curve->get_Direction(&curve_direction);
      bool bIsCCW = (curve_direction == cdLeft) ? true : false;

      if (0.0 < Lspiral[spEntry])
      {
         // there is an entry spiral

         // spiral starts at the Tangent to Spiral point (TS)
         auto ifc_ts = ConvertPoint<Schema>(pntTS);

         // tangent to spiral is the back tangent of the full curve
         CComPtr<IDirection> bkTangentBrg;
         curve->get_BkTangentBrg(&bkTangentBrg);
         Float64 bk_tangent_direction;
         bkTangentBrg->get_Value(&bk_tangent_direction);

         auto [geometry_segment, business_segment] = create_entry_spiral<Schema>(ifc_ts, bk_tangent_direction, (bIsCCW ? 1.0 : -1.0) * R, Lspiral[spEntry], options);
         file.addEntity(business_segment);
         alignment_segments->push(business_segment);
         if (geometry_segment)
         {
            file.addEntity(geometry_segment);
            curve_segments->push(geometry_segment);
         }
      }

      //
      // build the horizontal curve
      //

      // curve starts at the Spiral-to-Curve point
      CComPtr<IPoint2d> sc;
      curve->get_SC(&sc);
      auto ifc_sc = ConvertPoint<Schema>(sc);

      // tanget at the start is for the circular curve, not the full curve
      CComPtr<IDirection> bkTangentBrgCurve;
      curve->get_CurveBkTangentBrg(&bkTangentBrgCurve);
      Float64 bk_tangent_direction_curve;
      bkTangentBrgCurve->get_Value(&bk_tangent_direction_curve);

      auto [geometry_segment, business_segment] = create_hcurve<Schema>(ifc_sc, bk_tangent_direction_curve, (bIsCCW ? 1.0 : -1.0) * R, Lc, options);
      file.addEntity(business_segment);
      alignment_segments->push(business_segment);
      if (geometry_segment)
      {
         file.addEntity(geometry_segment);
         curve_segments->push(geometry_segment);
      }

      if (0.0 < Lspiral[spExit])
      {
         // there is an exit spiral

         // spiral starts at the Curve to Spiral point (CS)
         CComPtr<IPoint2d> pntCS;
         curve->get_CS(&pntCS);
         auto ifc_cs = ConvertPoint<Schema>(pntCS);

         CComPtr<IDirection> fwdTangentBrgCurve;
         curve->get_CurveFwdTangentBrg(&fwdTangentBrgCurve); // forward tangent of curve is start tangent to exit spiral
         Float64 fwd_tangent_direction_curve;
         fwdTangentBrgCurve->get_Value(&fwd_tangent_direction_curve);

         auto [geometry_segment, business_segment] = create_exit_spiral<Schema>(ifc_cs, fwd_tangent_direction_curve, (bIsCCW ? 1.0 : -1.0) * R, Lspiral[spExit], options);
         file.addEntity(business_segment);
         alignment_segments->push(business_segment);
         if (geometry_segment)
         {
            file.addEntity(geometry_segment);
            curve_segments->push(geometry_segment);
         }
      }

      // end of this curve (Spiral to Tangent, ST) becomes previous point for next alignment segment
      prevPoint.Release();
      curve->get_ST(&prevPoint);
      ifc_prev_point = ConvertPoint<Schema>(prevPoint);
   }

   // build a linear segment from end of previous alignment segment to the end of the alignment
   Float64 endStation, endElevation, endGrade;
   CComPtr<IPoint2d> endPoint;
   pAlignment->GetEndPoint(2, &endStation, &endElevation, &endGrade, &endPoint);

   GET_IFACE2(pBroker, IGeometry, pGeometry);
   Float64 dist;
   CComPtr<IDirection> direction;
   pGeometry->Inverse(prevPoint, endPoint, &dist, &direction);
   auto ifc_end_point = ConvertPoint<Schema>(endPoint);
   Float64 angle;
   direction->get_Value(&angle);

   if (prevPoint->SameLocation(endPoint) == S_FALSE)
   {
      // end the alignment with a line segment
      auto [geometry_segment, business_segment] = create_tangent<Schema>(ifc_prev_point, angle, dist, options);
      file.addEntity(business_segment);
      alignment_segments->push(business_segment);
      if (geometry_segment)
      {
         file.addEntity(geometry_segment);
         curve_segments->push(geometry_segment);
      }
   }


   // Add terminator segment
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Product_Shape/Product_Geometric_Representation/Alignment_Geometry/Alignment_Geometry_-_Horizontal_and_Vertical/content.html
   // 4.1.7.1.1.2 Zero length segment shall be added at the end of the list of segments
   // 4.1.7.1.1.2 If the geometry definition is present, then a zero length curve segment must be provided as well
   auto [geometry_segment, business_segment] = create_tangent<Schema>(ifc_end_point, angle, 0.0, options);
   file.addEntity(business_segment);
   alignment_segments->push(business_segment);
   if (geometry_segment)
   {
      // the last segment must have a DISCONTINUOUS transition code... create_tangent assumes continuous segments
      geometry_segment->setTransition(Schema::IfcTransitionCode::IfcTransitionCode_DISCONTINUOUS);
      file.addEntity(geometry_segment);
      curve_segments->push(geometry_segment);
   }

   // create a horizontal alignment from all the alignment segments
   // position the alignment relative to the site localPlacement
   auto horizontal_alignment = new Schema::IfcAlignmentHorizontal(IfcParse::IfcGlobalId(), nullptr, std::string("Horizontal Alignment"), boost::none, boost::none, nullptr, nullptr/*representation*/);
   file.addEntity(horizontal_alignment);
   auto site = file.getSingle<typename Schema::IfcSite>();
   file.relatePlacements(site, horizontal_alignment);

   // name the segments
   IndexType idx = 1;
   for (auto& segment : *alignment_segments)
   {
      std::ostringstream os;
      os << "H" << idx++;
      segment->setName(os.str());
   }

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Nests horizontal alignment segments with horizontal alignment"), horizontal_alignment, alignment_segments);
   file.addEntity(nests);

   *phorizontal_alignment = horizontal_alignment;
   *pnests_horizontal_segments = nests;

   auto composite_curve = new Schema::IfcCompositeCurve(curve_segments, false/*not self-intersecting*/);
   file.addEntity(composite_curve);
   *phorizontal_geometry_base_curve = composite_curve;
}

template <typename Schema>
void CreateVerticalProfile(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcCompositeCurve* horizontal_geometry_base_curve, const CIfcModelBuilderOptions& options, typename Schema::IfcAlignmentVertical** pvertical_profile, typename Schema::IfcRelNests** pnests_vertical_segments, typename Schema::IfcGradientCurve** palignment_gradient_curve)
{
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr profile_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
   typename aggregate_of<typename Schema::IfcSegment>::ptr curve_segments(new aggregate_of<typename Schema::IfcSegment>());

   // Profile is defined by profile segments located at "distance from start" of the alignment and "length".
   // We can't use stations to define the profile.
   // Distance from start is taken to be Station - Start Station

   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   Float64 startStation, startElevation, startGrade;
   auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   Float64 prev_end_dist_along = 0; // startStation; // this is distance along alignment, not station
   Float64 prev_end_gradient = startGrade;
   Float64 prev_end_height = startElevation;

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   IndexType nVCurves = pAlignment->GetVertCurveCount();
   for (IndexType i = 0; i < nVCurves; i++)
   {
      CComPtr<IVerticalCurve> curve;
      pAlignment->GetVertCurve(i, &curve);

      CComPtr<IProfilePoint> startPoint;
      curve->get_BVC(&startPoint);
      CComPtr<IStation> station;
      startPoint->get_Station(&station);
      Float64 start_height;
      startPoint->get_Elevation(&start_height);
      ZoneIndexType zoneIdx;
      Float64 start_dist_along;
      station->GetStation(&zoneIdx, &start_dist_along);
      start_dist_along -= startStation;
#pragma Reminder("How to deal with station equations?") // see IfcReferent

      if (!IsEqual(prev_end_dist_along, start_dist_along))
      {
         // create a linear segment between the last profile element and this curve
         Float64 length = start_dist_along - prev_end_dist_along;
         auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{prev_end_dist_along, prev_end_height});
         auto [geometry_segment, business_segment] = create_gradient<Schema>(vertical_point, prev_end_gradient, length, options);
         file.addEntity(business_segment);
         profile_segments->push(business_segment);
         if (geometry_segment)
         {
            file.addEntity(geometry_segment);
            curve_segments->push(geometry_segment);
         }
      }

      Float64 l1, l2;
      curve->get_L1(&l1);
      curve->get_L2(&l2);
      if (!IsEqual(l1, l2) && !IsZero(l2))
      {
         // compound vertical curve
         CComPtr<IProfilePoint> pviPoint;
         curve->get_PVI(&pviPoint);
         CComPtr<IStation> pviStation;
         pviPoint->get_Station(&pviStation);
         Float64 pviElevation;
         curve->Elevation(CComVariant(pviStation), &pviElevation);
         Float64 pviGrade;
         curve->Grade(CComVariant(pviStation), &pviGrade);

         Float64 start_gradient, end_gradient;
         curve->get_EntryGrade(&start_gradient);
         curve->get_ExitGrade(&end_gradient);

         auto vertical_point1 = new Schema::IfcCartesianPoint(std::vector<double>{start_dist_along, start_height});
         auto [geometry_segment1, business_segment1] = create_vcurve<Schema>(vertical_point1, start_gradient, pviGrade, l1, options);
         file.addEntity(business_segment1);
         profile_segments->push(business_segment1);
         if (geometry_segment1)
         {
            file.addEntity(geometry_segment1);
            curve_segments->push(geometry_segment1);
         }

         auto vertical_point2 = new Schema::IfcCartesianPoint(std::vector<double>{start_dist_along + l1, pviElevation});
         auto [geometry_segment2, business_segment2] = create_vcurve<Schema>(vertical_point2, pviGrade, end_gradient, l2, options);
         file.addEntity(business_segment2);
         profile_segments->push(business_segment2);
         if (geometry_segment2)
         {
            file.addEntity(geometry_segment2);
            curve_segments->push(geometry_segment2);
         }
      }
      else
      {
         Float64 horizontal_length;
         CComQIPtr<IProfileElement> element(curve);
         element->GetLength(&horizontal_length);
         Float64 start_gradient, end_gradient;
         curve->get_EntryGrade(&start_gradient);
         curve->get_ExitGrade(&end_gradient);

         if (IsEqual(start_gradient, end_gradient))
         {
            // this is just a straight line
            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{prev_end_dist_along, prev_end_height});
            auto [geometry_segment, business_segment] = create_gradient<Schema>(vertical_point, prev_end_gradient, l1, options);
            file.addEntity(business_segment);
            profile_segments->push(business_segment);
            if (geometry_segment)
            {
               file.addEntity(geometry_segment);
               curve_segments->push(geometry_segment);
            }
         }
         else
         {
            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{start_dist_along, start_height});
            auto [geometry_segment, business_segment] = create_vcurve<Schema>(vertical_point, start_gradient, end_gradient, horizontal_length, options);
            file.addEntity(business_segment);
            profile_segments->push(business_segment);
            if (geometry_segment)
            {
               file.addEntity(geometry_segment);
               curve_segments->push(geometry_segment);
            }
         }
      }

      // setup parameters for next loop
      CComPtr<IProfilePoint> evc;
      curve->get_EVC(&evc);
      CComPtr<IStation> evcStation;
      evc->get_Station(&evcStation);
      evcStation->GetStation(&zoneIdx, &prev_end_dist_along);
      prev_end_dist_along -= startStation;
      evc->get_Elevation(&prev_end_height);
      curve->get_ExitGrade(&prev_end_gradient);
#pragma Reminder("How to deal with station equations?") // see IfcReferent
   }

   Float64 endStation, endElevation, endGrade;
   CComPtr<IPoint2d> endPoint;
   pAlignment->GetEndPoint(2, &endStation, &endElevation, &endGrade, &endPoint);
   if (!IsEqual(prev_end_dist_along, endStation))
   {
      // create a linear segment between the last profile element and the end of the alignment
      ATLASSERT(IsEqual(prev_end_gradient, endGrade));
      Float64 length = endStation - startStation - prev_end_dist_along;
      auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{prev_end_dist_along, prev_end_height});
      auto [geometry_segment, business_segment] = create_gradient<Schema>(vertical_point, prev_end_gradient, length, options);
      file.addEntity(business_segment);
      profile_segments->push(business_segment);
      if (geometry_segment)
      {
         file.addEntity(geometry_segment);
         curve_segments->push(geometry_segment);
      }

      // check elevation
      ATLASSERT(IsEqual(endElevation, prev_end_height + length * prev_end_gradient));

      // parameters at end of segment
      prev_end_dist_along += length;
      prev_end_height += length * prev_end_gradient;

   }

   // Add terminator segment
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Product_Shape/Product_Geometric_Representation/Alignment_Geometry/Alignment_Geometry_-_Horizontal_and_Vertical/content.html
   // 4.1.7.1.1.2 Zero length segment shall be added at the end of the list of segments
   // 4.1.7.1.1.2 If the geometry definition is present, then a zero length curve segment must be provided as well
   auto terminator_vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{prev_end_dist_along, prev_end_height});
   auto [geometry_segment, business_segment] = create_gradient<Schema>(terminator_vertical_point, prev_end_gradient, 0.0, options);
   file.addEntity(business_segment);
   profile_segments->push(business_segment);
   if (geometry_segment)
   {
      // the last segment must have a DISCONTINUOUS transition code... create_gradient assumes continuous segments
      geometry_segment->setTransition(Schema::IfcTransitionCode::IfcTransitionCode_DISCONTINUOUS);

      file.addEntity(geometry_segment);
      curve_segments->push(geometry_segment);
   }


   // name the segments
   IndexType idx = 1;
   for (auto& segment : *profile_segments)
   {
      std::ostringstream os;
      os << "V" << idx++;
      segment->setName(os.str());
   }

   auto vertical_profile = new Schema::IfcAlignmentVertical(IfcParse::IfcGlobalId(), nullptr, std::string("Vertical Alignment"), boost::none, boost::none, file.getSingle<typename Schema::IfcLocalPlacement>(), nullptr);
   file.addEntity(vertical_profile);

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Nests vertical alignment segments with vertical alignment"), vertical_profile, profile_segments);
   file.addEntity(nests);

   *pvertical_profile = vertical_profile;
   *pnests_vertical_segments = nests;

   // define the roadway surface geometric representation with IfcSectionedSurface
#pragma Reminder("This geometric construction does not take into account the different ways PGSuper defines slope of the roadway section")
   //// expect surfaces to come out wrong until this is made more robust. for now, stick with simple crowns
   //
   //typename aggregate_of<typename Schema::IfcProfileDef>::ptr cross_sections(new aggregate_of<typename Schema::IfcProfileDef>()); // container of roadway cross sections
   //typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions( new aggregate_of<typename Schema::IfcAxis2PlacementLinear>()); // container of positions where cross sections located

   //GET_IFACE2(pBroker, IRoadwayData, pRoadway);
   //const RoadwaySectionData& roadway_sections = pRoadway->GetRoadwaySectionData();
   //Float64 ref_station = pRoadway->GetAlignmentData2().RefStation;
   //for (const RoadwaySectionTemplate& section_template : roadway_sections.RoadwaySectionTemplates)
   //{
   //   // create a point on the alignment curve for this cross section template
   //   auto point_on_alignment = new Schema::IfcPointByDistanceExpression(
   //      new Schema::IfcLengthMeasure(section_template.Station - ref_station),  // distance from start of curve
   //      boost::none, // lateral offset
   //      boost::none, // vertical offset
   //      boost::none, // longitudinal offset
   //      horizontal_geometry_base_curve // the basis curve (eg, the alignment curve)
   //   );

   //   // create a linear placement object to position alignment point in space
   //   auto linear_placement_of_roadway_section = new Schema::IfcAxis2PlacementLinear(point_on_alignment, nullptr/*local z axis*/, nullptr/*ref direction to determine local x axis*/);

   //   // add the cross section placement into the container
   //   cross_section_positions->push(linear_placement_of_roadway_section);

   //   // start building the section profile by defining the widths and slopes
   //   std::vector<double> widths;
   //   std::vector<double> slopes;
   //   if (roadway_sections.NumberOfSegmentsPerSection == 2)
   //   {
   //      // this should be restricted to the bridge width instead of 100 m

   //      // slope measure type needs to be considered here
   //      widths.push_back(100);
   //      slopes.push_back(section_template.LeftSlope);
   //      widths.push_back(100);
   //      slopes.push_back(section_template.RightSlope);
   //   }
   //   else
   //   {
   //      for (const RoadwaySegmentData& segment_data : section_template.SegmentDataVec)
   //      {
   //         // slope measure type needs to be considered here
   //         widths.push_back(segment_data.Length);
   //         slopes.push_back(segment_data.Slope);
   //      }
   //   }

   //   // create a human readable label for this cross section profile
   //   std::ostringstream os;
   //   os << "Roadway Template at Station " << (LPCSTR)::FormatStation(pDisplayUnits->GetStationFormat(), section_template.Station).GetBuffer() << std::endl;

   //   // create the cross section profile
   //   auto cross_section = new Schema::IfcOpenCrossProfileDef(
   //      Schema::IfcProfileTypeEnum::IfcProfileType_CURVE, // Profile is treated as a curve that will be used in conjunction with a swept surface (otherwise, area makes a swept solid)
   //      os.str(), // optional profile name - this is just human readable information
   //      true, // widths are horizontal, not alone the slope
   //      widths, // left to right, widths of the profile line elements
   //      slopes, // left to right, slopes of the profile line elements
   //      boost::none, // optional list of tags. used to match points between sequential profile definitions when there are different number of points per profile
   //      nullptr); // point to designate as the start point of the profile. If nullptr, profile starts at the alignment
   //
   //   // add this cross section to the container of cross sections
   //   cross_sections->push(cross_section); 
   //}

   //// using the cross sections and their positions, create a sectioned surface attached to the alignment curve
   //auto sectioned_surface = new Schema::IfcSectionedSurface(horizontal_geometry_base_curve, cross_section_positions, cross_sections);
   //representation_items->push(sectioned_surface);

   auto gradient_curve = new Schema::IfcGradientCurve(curve_segments, false, horizontal_geometry_base_curve, nullptr);
   file.addEntity(gradient_curve);
   *palignment_gradient_curve = gradient_curve;
}

// creates representations for each IfcAlignmentSegment per CT 4.1.7.1.1.4
// https://standards.buildingsmart.org/IFC/RELEASE/IFC4_3/HTML/concepts/Product_Shape/Product_Geometric_Representation/Alignment_Geometry/Alignment_Geometry_-_Segments/content.html
template <typename Schema>
void CreateAlignmentSegmentRepresentations(IfcHierarchyHelper<typename Schema>& file, typename Schema::IfcLocalPlacement* global_placement, typename Schema::IfcGeometricRepresentationSubContext* segment_axis_subcontext, typename aggregate_of<typename Schema::IfcSegment>::ptr curve_segments, typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr segments)
{
   auto cs_iter = curve_segments->begin();
   auto s_iter = segments->begin();
   for (; cs_iter != curve_segments->end(); cs_iter++, s_iter++)
   {
      auto curve_segment = *cs_iter;
      auto alignment_segment = (*s_iter)->as<typename Schema::IfcAlignmentSegment>();

      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      representation_items->push(curve_segment);

      auto axis_representation = new Schema::IfcShapeRepresentation(segment_axis_subcontext, std::string("Axis"), std::string("Segment"), representation_items);
      file.addEntity(axis_representation);

      typename aggregate_of<typename Schema::IfcRepresentation>::ptr representations(new aggregate_of<typename Schema::IfcRepresentation>());
      representations->push(axis_representation);

      auto product = new Schema::IfcProductDefinitionShape(std::string("Product Definition of a Segment"), boost::none, representations);
      file.addEntity(product);

      alignment_segment->setObjectPlacement(global_placement);
      alignment_segment->setRepresentation(product);
   }
}

template <typename Schema>
void CreateAlignment(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   typename Schema::IfcProductDefinitionShape* alignment_representation = nullptr;

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);

   // Need Axis representation for Polyline, Gradient, and Segments
   auto axis_model_representation_subcontext = new Schema::IfcGeometricRepresentationSubContext(std::string("Axis"), std::string("Model"), geometric_representation_context, boost::none, Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_MODEL_VIEW, boost::none);
   file.addEntity(axis_model_representation_subcontext);

   typename Schema::IfcAlignmentHorizontal* horizontal_alignment_layout = nullptr;
   typename Schema::IfcAlignmentVertical* vertical_profile_layout = nullptr;
   typename Schema::IfcCompositeCurve* composite_curve = nullptr;
   typename Schema::IfcGradientCurve* gradient_curve = nullptr;
   typename Schema::IfcPolyline* polyline = nullptr;

   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      typename Schema::IfcRelNests* nests_horizontal_segments;
      CreateHorizontalAlignment<Schema>(file, pBroker, options, &horizontal_alignment_layout, &nests_horizontal_segments, &composite_curve);

      typename Schema::IfcRelNests* nests_vertical_segments;
      CreateVerticalProfile<Schema>(file, pBroker, composite_curve, options, &vertical_profile_layout, &nests_vertical_segments, &gradient_curve);

      // Need FootPrint representation for Horizontal+Vertical composite curve
      typename Schema::IfcGeometricRepresentationSubContext* footprint_model_representation_subcontext = nullptr;
      if (options.representations == CIfcModelBuilderOptions::Representations::Curve3dAndFootPrint)
      {
         footprint_model_representation_subcontext = new Schema::IfcGeometricRepresentationSubContext(std::string("FootPrint"), std::string("Model"), geometric_representation_context, boost::none, Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_MODEL_VIEW, boost::none);
         file.addEntity(footprint_model_representation_subcontext);
      }

      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr horizontal_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      horizontal_representation_items->push(composite_curve);

      typename Schema::IfcShapeRepresentation* footprint_curve2d_shape_representation = nullptr;
      if (options.representations == CIfcModelBuilderOptions::Representations::Curve3dAndFootPrint)
      {
         footprint_curve2d_shape_representation = new Schema::IfcShapeRepresentation(footprint_model_representation_subcontext, std::string("FootPrint"), std::string("Curve2D"), horizontal_representation_items);
         file.addEntity(footprint_curve2d_shape_representation);
      }

      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr vertical_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      vertical_representation_items->push(gradient_curve);

      auto curve3d_shape_representation = new Schema::IfcShapeRepresentation(axis_model_representation_subcontext, std::string("Axis"), std::string("Curve3D"), vertical_representation_items);
      file.addEntity(curve3d_shape_representation);

      typename aggregate_of<typename Schema::IfcRepresentation>::ptr representations(new aggregate_of<typename Schema::IfcRepresentation>());
      if (options.representations == CIfcModelBuilderOptions::Representations::Curve3dAndFootPrint)
      {
         representations->push(footprint_curve2d_shape_representation); // 2D alignment geometry (Horizontal + Vertical)
      }
      representations->push(curve3d_shape_representation); // 3D alignment geometry (Horizontal + Vertical)
      alignment_representation = new Schema::IfcProductDefinitionShape(std::string("Alignment Product Definition Shape"), boost::none, representations);
      // this alignment_representation will be assigned to the IfcAlignment when it is created a little further down.

      // loops over all the individual segments in the horizontal and vertical alignments setting up 'Axis' 'Segment' representations for each individual segment
      auto global_placement = file.addLocalPlacement();
      CreateAlignmentSegmentRepresentations(file, global_placement, axis_model_representation_subcontext, composite_curve->Segments(), nests_horizontal_segments->RelatedObjects());
      if (gradient_curve)
      {
         CreateAlignmentSegmentRepresentations(file, global_placement, axis_model_representation_subcontext, gradient_curve->Segments(), nests_vertical_segments->RelatedObjects());
      }
   }
   else
   {
      // Instead of IfcGradientCurve, we are using a generalized 3D polyline geometric representation of the alignment (a 3D wire)
      // This isn't as accurate, but some viewer may be able to deal with this better
      GET_IFACE2(pBroker, IRoadway, pAlignment);

      Float64 startStation, startElevation, startGrade;
      CComPtr<IPoint2d> startPoint;
      pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

      Float64 endStation, endElevation, endGrade;
      CComPtr<IPoint2d> endPoint;
      pAlignment->GetEndPoint(2, &endStation, &endElevation, &endGrade, &endPoint);

      IndexType nAlignmentPoints = 100;
      Float64 stationInc = (endStation - startStation) / (nAlignmentPoints + 1);
      typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());
      for (IndexType i = 0; i <= nAlignmentPoints; i++)
      {
         Float64 offset = 0.0;
         Float64 station = startStation + i * stationInc;
         CComPtr<IPoint2d> pnt;
         pAlignment->GetPoint(station, offset, nullptr /*normal offset*/, pgsTypes::pcGlobal, &pnt);
         Float64 x, y;
         pnt->Location(&x, &y);
         Float64 z = pAlignment->GetElevation(station, offset);

         points->push(new Schema::IfcCartesianPoint(std::vector<Float64>{x, y, z}));
      }
      polyline = new Schema::IfcPolyline(points);

      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr alignment_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      alignment_representation_items->push(polyline);

      auto curve3d_shape_representation = new Schema::IfcShapeRepresentation(axis_model_representation_subcontext, std::string("Axis"), std::string("Curve3D"), alignment_representation_items);
      file.addEntity(curve3d_shape_representation);

      typename aggregate_of<typename Schema::IfcRepresentation>::ptr representations(new aggregate_of<typename Schema::IfcRepresentation>());
      representations->push(curve3d_shape_representation); // 3D alignment geometry (Horizontal + Vertical)
      alignment_representation = new Schema::IfcProductDefinitionShape(std::string("Alignment Product Definition Shape"), boost::none, representations);
      // this alignment_representation will be assigned to the IfcAlignment when it is created a little further down.
   }

   // place the alignment relative to the site
   auto site = file.getSingle<typename Schema::IfcSite>();
   auto local_placement = site->ObjectPlacement();
   if (!local_placement)
   {
      local_placement = file.addLocalPlacement();
   }

   GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
   std::string strAlignmentName(T2A(pRoadwayData->GetAlignmentData2().Name.c_str()));
   if (strAlignmentName.empty()) strAlignmentName = "Unnamed alignment";
   auto alignment = new Schema::IfcAlignment(IfcParse::IfcGlobalId(), nullptr, strAlignmentName, boost::none, boost::none, local_placement, alignment_representation, boost::none);
   file.addEntity(alignment);

   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
   {
      // 4.1.4.4.1 Alignments nest horizontal and vertical layouts
      // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Object_Composition/Nesting/Alignment_Layouts/content.html
      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr alignment_layout_list(new aggregate_of<typename Schema::IfcObjectDefinition>());
      alignment_layout_list->push(horizontal_alignment_layout);
      alignment_layout_list->push(vertical_profile_layout);

      auto nests_alignment_layouts = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, std::string("Nest horizontal and vertical alignment layouts with the alignment"), boost::none, alignment, alignment_layout_list);
      file.addEntity(nests_alignment_layouts);
   }

   // IFC 4.1.4.1.1 "Every IfcAlignment must be related to IfcProject using the IfcRelAggregates relationship"
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Object_Composition/Aggregation/Alignment_Aggregation_To_Project/content.html
   // IfcProject <-> IfcRelAggregates <-> IfcAlignment
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_alignments_in_project(new aggregate_of<typename Schema::IfcObjectDefinition>());
   list_of_alignments_in_project->push(alignment);
   auto project = file.getSingle<typename Schema::IfcProject>();
   auto aggregate_alignments_with_project = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Alignments in project"), boost::none, project, list_of_alignments_in_project);
   file.addEntity(aggregate_alignments_with_project);

   // IFC 4.1.5.1 alignment is referenced in spatial structure of an IfcSpatialElement. In this case IfcSite is the highest level IfcSpatialElement
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Object_Connectivity/Alignment_Spatial_Reference/content.html
   // IfcSite <-> IfcRelReferencedInSpatialStructure <-> IfcAlignment
   // This means IfcAlignment is not part of the IfcSite (it is not an aggregate component) but instead IfcAlignment is used within
   // the IfcSite by reference. This implies an IfcAlignment can traverse many IfcSite instances within an IfcProject
   typename Schema::IfcSpatialReferenceSelect::list::ptr list_alignments_referenced_in_site(new Schema::IfcSpatialReferenceSelect::list);
   list_alignments_referenced_in_site->push(alignment);
   auto rel_referenced_in_spatial_structure = new Schema::IfcRelReferencedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, list_alignments_referenced_in_site, site);
   file.addEntity(rel_referenced_in_spatial_structure);
}


template <typename Schema>
typename Schema::IfcCurve* GetAlignmentDirectrix(IfcHierarchyHelper<Schema>& file, const CIfcModelBuilderOptions& options)
{
   // get the directrix line of the alignment
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();
   auto alignment_representation = alignment->Representation();
   auto alignment_representations = alignment_representation->Representations();
   for (auto& representation : *alignment_representations)
   {
      auto alignment_representation_items = representation->Items();
      for (auto& representation_item : *alignment_representation_items)
      {
         auto directrix = representation_item->as<typename Schema::IfcCurve>();
         if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
         {
            if (directrix->as<typename Schema::IfcGradientCurve>())
               return directrix;
         }
         else
         {
            // if we are using a polyline, then the directrix is the polyline curve (which is 3d, not 2d, in this case)
            if (directrix->as<typename Schema::IfcPolyline>())
               return directrix;
         }
      }
   }
   CHECK(false); // didn't find the alignment curve
   return nullptr;
}