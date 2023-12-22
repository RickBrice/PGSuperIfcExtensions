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
#include "IfcModelBuilder.h"

#include <IFace\Project.h>
#include <IFace\VersionInfo.h>
#include <IFace\Alignment.h>
#include <IFace\DocumentType.h>
#include <IFace\Bridge.h>
#include <IFace\Intervals.h>
#include <IFace\PrestressForce.h>
#include <EAF\EAFDisplayUnits.h>

#include <PgsExt\GirderLabel.h>
#include <PgsExt\PrecastSegmentData.h>

#include <WBFLCogo\CogoHelpers.h>

#define CLOCKWISE 0
#define COUNTERCLOCKWISE 1
int GetVertexOrdering(IShape* pShape)
{
    Float64 area;

    // Initialize and check for null polygon.
    area = 0;

    CComPtr<IPoint2dCollection> points;
    pShape->get_PolyPoints(&points);

    IndexType cPoints;
    points->get_Count(&cPoints);

    if (cPoints < 3)
    {
        return CLOCKWISE;
    }

    Float64 x0, y0;
    Float64 x1, y1;
    Float64 dy, dx;
    Float64 ar, at;

    // loop over all points - make sure of closure
    IndexType idx0, idx1;
    idx0 = 0;
    idx1 = 1;

    bool loop = true;
    bool last_round = false;

    while (loop)
    {
        CComPtr<IPoint2d> p0;
        CComPtr<IPoint2d> p1;

        points->get_Item(idx0, &p0);
        points->get_Item(idx1, &p1);

        p0->Location(&x0, &y0);
        p1->Location(&x1, &y1);

        dx = x1 - x0;
        dy = y1 - y0;

        ar = dx * y0;
        at = 0.5 * dy * dx;

        area += (ar + at);

        // loop termination test - need to go one more iteration if loop is not closed
        if (last_round)
        {
            // just finished closure loop. time to quit
            loop = false;
        }
        else
        {
            // increment for next go-around
            idx0++;
            idx1++;

            if (idx0 == cPoints - 1)
            {
                idx0 = cPoints - 1;
                idx1 = 0;

                // check if extra loop is required for closure
                CComPtr<IPoint2d> pStart;
                CComPtr<IPoint2d> pEnd;
                points->get_Item(idx1, &pStart);
                points->get_Item(idx0, &pEnd);

                if (pStart->SameLocation(pEnd) == S_FALSE)
                {
                    // one more loop to close poly
                    last_round = true;
                }
                else
                {
                    // loop is closed - just quit
                    loop = false;
                }
            }
        }
    }     // while

    if (area < 0)
    {
        return COUNTERCLOCKWISE;
    }
    else
    {
        return CLOCKWISE;
    }
}

template <typename Schema>
typename Schema::IfcCartesianPoint* ConvertPoint(IPoint2d* pPoint)
{
   Float64 x, y;
   pPoint->Location(&x, &y);
   return new Schema::IfcCartesianPoint(std::vector<double>{x, y});
}

template <typename Schema>
typename Schema::IfcPolyline* CreatePolyline(IShape* shape)
{
   CComPtr<IPoint2dCollection> polyPoints;
   shape->get_PolyPoints(&polyPoints);

   if (GetVertexOrdering(shape) == CLOCKWISE)
   {
      polyPoints->Reverse();
   }

   typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());
   IndexType nPoints;
   polyPoints->get_Count(&nPoints);

   // polygon must be closed and it must be closed by reference, not different points at same location
   CComPtr<IPoint2d> first, last;
   polyPoints->get_Item(0, &first);
   polyPoints->get_Item(nPoints - 1, &last);
   if (first->SameLocation(last) == S_OK)
   {
      // first and last point are at same location so polyPoints is a closed polygon
      // reduce the number of points traversed by 1 so points is open
      nPoints--;
   }

   for (IndexType idx = 0; idx < nPoints; idx++)
   {
      CComPtr<IPoint2d> point;
      polyPoints->get_Item(idx, &point);
      points->push(ConvertPoint<Schema>(point));
   }

   // we know that points is open (last point is not the same as the first)
   // the polygon must be closed by reference
   points->push(*(points->begin()));

   auto polyline = new Schema::IfcPolyline(points);

   return polyline;
}

template <typename Schema>
typename Schema::IfcCurve* GetAlignmentDirectrix(IfcHierarchyHelper<Schema>& file)
{
   // get the directrix line of the alignment
   auto alignment = file.getSingle<typename Schema::IfcAlignment>();
   auto alignment_representation = alignment->Representation();
   auto alignment_representations = alignment_representation->Representations();
   auto representation = *(alignment_representations->begin() + 1);
   auto alignment_representation_items = representation->Items();
   auto representation_item = *(alignment_representation_items->begin());
   auto directrix = representation_item->as<typename Schema::IfcCurve>();
   CHECK(directrix->as<typename Schema::IfcGradientCurve>());
   return directrix;
}


template <typename Schema>
void CreateHorizontalAlignment(IfcHierarchyHelper<Schema>& file,IBroker* pBroker, const CIfcModelBuilderOptions& options, typename Schema::IfcAlignmentHorizontal** phorizontal_alignment, typename Schema::IfcCompositeCurve** phorizontal_geometry_base_curve)
{
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr alignment_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
   typename aggregate_of<typename Schema::IfcSegment>::ptr curve_segments(new aggregate_of<typename Schema::IfcSegment>());

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   // create the start point
   auto ifc_start_point = ConvertPoint<Schema>(startPoint);

   // loop over all the horizontal curves
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
      if(pntTS->SameLocation(prevPoint) == S_FALSE)
      {
         GET_IFACE2(pBroker, IGeometry, pGeometry);
         Float64 dist;
         CComPtr<IDirection> direction;
         pGeometry->Inverse(prevPoint, pntTS, &dist, &direction);
         Float64 angle;
         direction->get_Value(&angle);

         auto ifc_line_segment = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_prev_point, angle, 0.0, 0.0, dist, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
         file.addEntity(ifc_line_segment);
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, ifc_line_segment);
         file.addEntity(segment);
         alignment_segments->push(segment);

         auto placement = new Schema::IfcAxis2Placement2D(ifc_prev_point, new Schema::IfcDirection(std::vector<double>{cos(angle), sin(angle)}));
         if (options.tangents == CIfcModelBuilderOptions::Tangents::Polyline)
         {
            typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());
            points->push(ifc_prev_point);
            points->push(ConvertPoint<Schema>(pntTS));
            auto parent_curve = new Schema::IfcPolyline(points);
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcLengthMeasure(0.0), new Schema::IfcLengthMeasure(dist), parent_curve);
            file.addEntity(curve_segment);
            curve_segments->push(curve_segment);
         }
         else
         {
            auto parent_curve = new Schema::IfcLine(
               new Schema::IfcCartesianPoint(std::vector<double>({ 0,0 })),
               new Schema::IfcVector(new Schema::IfcDirection(std::vector<double>{1.0, 0.0}), 1.0));

            auto curve_segment = new Schema::IfcCurveSegment(
               Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT,
               new Schema::IfcAxis2Placement2D(ifc_prev_point, new Schema::IfcDirection(std::vector<double>{cos(angle), sin(angle)})),
               new Schema::IfcLengthMeasure(0.0), // start
               new Schema::IfcLengthMeasure(dist),
               parent_curve);

            file.addEntity(curve_segment);
            curve_segments->push(curve_segment);
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

      // curve starts at the Tangent to Spiral point (TS)
      auto ifc_ts = ConvertPoint<Schema>(pntTS);

      CComPtr<IDirection> bkTangentBrg;
      curve->get_BkTangentBrg(&bkTangentBrg);
      Float64 bk_tangent_direction;
      bkTangentBrg->get_Value(&bk_tangent_direction);

      if (0.0 < Lspiral[spEntry])
      {
         // there is an entry spiral


         auto ifc_entry_spiral = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_ts, bk_tangent_direction, 0, (bIsCCW ? 1.0 : -1.0)*R, Lspiral[spEntry], boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, ifc_entry_spiral);
         file.addEntity(segment);
         alignment_segments->push(segment);

         Float64 A = (bIsCCW ? 1.0 : -1.0) * sqrt(Lspiral[spEntry] * R);
         
         auto parent_curve = new Schema::IfcClothoid(
            new Schema::IfcAxis2Placement2D(
               new Schema::IfcCartesianPoint(std::vector<double>({ 0,0 })), 
               new Schema::IfcDirection(std::vector<double>({ 1,0 }))),
            A);

         auto curve_segment = new Schema::IfcCurveSegment(
            Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, 
            new Schema::IfcAxis2Placement2D(ifc_ts, new Schema::IfcDirection(std::vector<double>({cos(bk_tangent_direction),sin(bk_tangent_direction)}))),
            new Schema::IfcLengthMeasure(0.0), 
            new Schema::IfcLengthMeasure(Lspiral[spEntry]), 
            parent_curve);

         file.addEntity(parent_curve);
         file.addEntity(curve_segment);
         curve_segments->push(curve_segment);
      }

      // build the horizontal curve
      // curve starts at the Spiral-to-Curve point
      CComPtr<IPoint2d> sc;
      curve->get_SC(&sc);
      auto ifc_sc = ConvertPoint<Schema>(sc);

      CComPtr<IDirection> bkTangentBrgCurve;
      curve->get_CurveBkTangentBrg(&bkTangentBrgCurve);
      Float64 bk_tangent_direction_curve;
      bkTangentBrgCurve->get_Value(&bk_tangent_direction_curve);
      auto dx = cos(bk_tangent_direction_curve);
      auto dy = sin(bk_tangent_direction_curve);

      auto ifc_hcurve = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_sc, bk_tangent_direction_curve, (bIsCCW ? 1.0 : -1.0)*R, (bIsCCW ? 1.0 : -1.0)*R, Lc, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC);
      auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, ifc_hcurve);
      file.addEntity(segment);
      alignment_segments->push(segment);


      auto parent_curve = new Schema::IfcCircle(
         new Schema::IfcAxis2Placement2D(
            new Schema::IfcCartesianPoint(std::vector<double>{0,0}),
            new Schema::IfcDirection(std::vector<double>({1,0}))),
         R);
      
      auto curve_segment = new Schema::IfcCurveSegment(
         Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, 
         new Schema::IfcAxis2Placement2D(ifc_sc, new Schema::IfcDirection(std::vector<double>({ dx, dy }))),
         new Schema::IfcLengthMeasure(0.0), 
         new Schema::IfcLengthMeasure((bIsCCW ? 1.0 : -1.0) * Lc), 
         parent_curve);

      file.addEntity(parent_curve);
      file.addEntity(curve_segment);
      curve_segments->push(curve_segment);

      if (0.0 < Lspiral[spExit])
      {
         // there is an exit spiral

         CComPtr<IDirection> fwdTangentBrgCurve;
         curve->get_CurveFwdTangentBrg(&fwdTangentBrgCurve); // forward tangent of curve is start tangent to exit spiral
         Float64 fwd_tangent_direction_curve;
         fwdTangentBrgCurve->get_Value(&fwd_tangent_direction_curve);

         CComPtr<IDirection> fwdTangentBrgSpiral;
         curve->get_FwdTangentBrg(&fwdTangentBrgSpiral); // forward tangent of curve at the end of the exit spiral
         Float64 fwd_tangent_direction_spiral; 
         fwdTangentBrgSpiral->get_Value(&fwd_tangent_direction_spiral);

         // spiral starts at the Curve to Spiral point (CS)
         CComPtr<IPoint2d> pntCS;
         curve->get_CS(&pntCS);
         auto ifc_cs = ConvertPoint<Schema>(pntCS);

         // The geometry of the spiral (point of infinite radius) begins at the Spiral to Tangent point (ST)
         CComPtr<IPoint2d> pntST;
         curve->get_ST(&pntST);
         auto ifc_st = ConvertPoint<Schema>(pntST);

         auto ifc_exit_spiral = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_cs, fwd_tangent_direction_curve, (bIsCCW ? 1.0 : -1.0)*R, 0, Lspiral[spExit], boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, ifc_exit_spiral);
         file.addEntity(segment);
         alignment_segments->push(segment);

         Float64 A = (bIsCCW ? -1.0 : 1.0) * sqrt(Lspiral[spExit] * R);

         auto parent_curve = new Schema::IfcClothoid(
            new Schema::IfcAxis2Placement2D(
               new Schema::IfcCartesianPoint(std::vector<double>({ 0,0 })), 
               new Schema::IfcDirection(std::vector<double>({ 1,0 }))),
            A);

         auto curve_segment = new Schema::IfcCurveSegment(
            Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT,
            new Schema::IfcAxis2Placement2D(ifc_cs, new Schema::IfcDirection(std::vector<double>({ cos(fwd_tangent_direction_curve),sin(fwd_tangent_direction_curve) }))),
            new Schema::IfcLengthMeasure(-1.0 * Lspiral[spExit]),
            new Schema::IfcLengthMeasure(1.0*Lspiral[spExit]),
            parent_curve);

         file.addEntity(parent_curve);
         file.addEntity(curve_segment);
         curve_segments->push(curve_segment);
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

   if(prevPoint->SameLocation(endPoint) == S_FALSE)
   {
      // end the alignment with a line segment
      auto ifc_line_segment = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_prev_point, angle, 0.0, 0.0, dist, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
      file.addEntity(ifc_line_segment);
      auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, ifc_line_segment);
      file.addEntity(segment);
      alignment_segments->push(segment);
   
      auto placement = new Schema::IfcAxis2Placement2D(ifc_prev_point, new Schema::IfcDirection(std::vector<double>{cos(angle), sin(angle)}));
      if (options.tangents == CIfcModelBuilderOptions::Tangents::Polyline)
      {
         typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());
         points->push(ifc_prev_point);
         points->push(ifc_end_point);
         auto parent_curve = new Schema::IfcPolyline(points);
         auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcLengthMeasure(0.0), new Schema::IfcLengthMeasure(dist), parent_curve);
         file.addEntity(curve_segment);
         curve_segments->push(curve_segment);
      }
      else
      {
         auto dir = new Schema::IfcDirection(std::vector<double>{1.0,0.0});
         double mag = 1.0;
         auto vector = new Schema::IfcVector(dir, mag);
         auto parent_curve = new Schema::IfcLine(new Schema::IfcCartesianPoint(std::vector<double>{0, 0}), vector);
         auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcLengthMeasure(0.0), new Schema::IfcLengthMeasure(dist), parent_curve);
         file.addEntity(curve_segment);
         curve_segments->push(curve_segment);
      }
   }


   // Add terminator segment
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Product_Shape/Product_Geometric_Representation/Alignment_Geometry/Alignment_Geometry_-_Horizontal_and_Vertical/content.html
   // 4.1.7.1.1.2 Zero length segment shall be added at the end of the list of segments
   auto terminator_line_segment = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_end_point, 0.0, 0.0, 0.0, 0.0, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
   file.addEntity(terminator_line_segment);
   auto terminator_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, terminator_line_segment);
   file.addEntity(terminator_segment);
   alignment_segments->push(terminator_segment);

   // 4.1.7.1.1.2 If the geometry definition is present, then a zero length curve segment must be provided as well
   auto segment_curve_placement = new Schema::IfcAxis2Placement2D(ifc_end_point, new Schema::IfcDirection(std::vector<double>{cos(angle), sin(angle)}));
   auto dir = new Schema::IfcDirection(std::vector<double>{1.0,0.0});
   double mag = 1.0;
   auto vector = new Schema::IfcVector(dir, mag);
   auto terminator_parent_curve = new Schema::IfcLine(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), vector);
   auto terminator_curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_DISCONTINUOUS, segment_curve_placement, new Schema::IfcLengthMeasure(0.0), new Schema::IfcLengthMeasure(0.0), terminator_parent_curve);
   file.addEntity(terminator_curve_segment);
   curve_segments->push(terminator_curve_segment);

   // create a horizontal alignment from all the alignment segments
   // position the alignment relative to the site localPlacement
   auto horizontal_alignment = new Schema::IfcAlignmentHorizontal(IfcParse::IfcGlobalId(), nullptr, std::string("Horizontal Alignment"), boost::none, boost::none, nullptr, nullptr/*representation*/);
   file.addEntity(horizontal_alignment);
   auto site = file.getSingle<typename Schema::IfcSite>();
   file.relatePlacements(site, horizontal_alignment);

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Nests horizontal alignment segments with horizontal alignment"), horizontal_alignment, alignment_segments);
   file.addEntity(nests);

   *phorizontal_alignment = horizontal_alignment;

   auto composite_curve = new Schema::IfcCompositeCurve(curve_segments, false/*not self-intersecting*/);
   file.addEntity(composite_curve);
   *phorizontal_geometry_base_curve = composite_curve;
}

template <typename Schema>
void CreateVerticalProfile(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcCompositeCurve* horizontal_geometry_base_curve, const CIfcModelBuilderOptions& options, typename Schema::IfcAlignmentVertical** pvertical_profile, typename Schema::IfcGradientCurve** palignment_gradient_curve)
{
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr profile_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
   typename aggregate_of<typename Schema::IfcSegment>::ptr curve_segments(new aggregate_of<typename Schema::IfcSegment>());

   // Profile is defined by profile segments located at "distance from start" of the alignment and "length".
   // We can't use stations to define the profile.
   // Distance from start is taken to be Station - Start Station

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   GET_IFACE2_NOCHECK(pBroker, IEAFDisplayUnits, pDisplayUnits);

   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   Float64 prev_end_dist_along = 0; // startStation; // this is distance along alignment, not station
   Float64 prev_end_gradient = startGrade;
   Float64 prev_end_height = startElevation;

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
         auto* linear_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, prev_end_dist_along, length, prev_end_height, prev_end_gradient, prev_end_gradient, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
         file.addEntity(linear_segment);

         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, linear_segment);
         file.addEntity(segment);
         profile_segments->push(segment);

         if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::GradientCurve)
         {
            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
            auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(prev_end_gradient), sin(prev_end_gradient)});
            auto parent_curve = new Schema::IfcLine(vertical_point, new Schema::IfcVector(gradient_direction, 1.0));
            auto segment_curve_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0,0.0}), new Schema::IfcDirection(std::vector<double>{1.0,0.0}));
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, segment_curve_placement, new Schema::IfcLengthMeasure(prev_end_dist_along), new Schema::IfcLengthMeasure(length), parent_curve);
            file.addEntity(curve_segment);
            curve_segments->push(curve_segment);
         }
      }

      Float64 l1, l2;
      curve->get_L1(&l1);
      curve->get_L2(&l2);
      if (!IsEqual(l1, l2) && !IsZero(l2))
      {
         // compound curve
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

         Float64 k1, k2;
         curve->get_K1(&k1);
         curve->get_K2(&k2);

         auto* parabolic_segment1 = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along, l1, start_height, start_gradient, pviGrade, IsZero(k1) ? Float64_Max : 1 / k1, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
         auto segment1 = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, parabolic_segment1);
         file.addEntity(segment1);
         profile_segments->push(segment1);

         Float64 A1 = (pviGrade - start_gradient) / (2 * l1);
         Float64 B1 = start_gradient;
         Float64 C1 = start_height;
         auto parent_curve1_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
         auto parent_curve1 = new Schema::IfcPolynomialCurve(parent_curve1_placement, std::vector<double>{0.0,1.0}, std::vector<double>{C1, B1, A1}, boost::none);
         auto segment_curve1_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
         auto curve_segment1 = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, segment_curve1_placement, new Schema::IfcLengthMeasure(start_dist_along), new Schema::IfcLengthMeasure(l1), parent_curve1);
         file.addEntity(curve_segment1);
         curve_segments->push(curve_segment1);

         auto* parabolic_segment2 = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along + l1, l2, pviElevation, pviGrade, end_gradient, IsZero(k2) ? Float64_Max : 1 / k2, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
         auto segment2 = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, parabolic_segment2);
         file.addEntity(segment2);
         profile_segments->push(segment2);

         Float64 A2 = (end_gradient - pviGrade) / (2 * l2);
         Float64 B2 = pviGrade;
         Float64 C2 = pviElevation;
         auto parent_curve2_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
         auto parent_curve2 = new Schema::IfcPolynomialCurve(parent_curve2_placement, boost::none, std::vector<double>{C2, B2, A2}, boost::none);
         auto segment_curve2_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
         auto curve_segment2 = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, segment_curve2_placement, new Schema::IfcLengthMeasure(start_dist_along + l1), new Schema::IfcLengthMeasure(l2), parent_curve2);
         file.addEntity(curve_segment2);
         curve_segments->push(curve_segment2);
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
            auto* linear_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, prev_end_dist_along, l1, prev_end_height, start_gradient, start_gradient, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
            auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, linear_segment);
            file.addEntity(segment);
            profile_segments->push(segment);

            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
            auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(start_gradient), sin(start_gradient)});
            auto parent_curve = new Schema::IfcLine(vertical_point, new Schema::IfcVector(gradient_direction, 1.0));
            auto segment_curve_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, segment_curve_placement, new Schema::IfcLengthMeasure(prev_end_dist_along), new Schema::IfcLengthMeasure(l1), parent_curve);
            file.addEntity(curve_segment);
            curve_segments->push(curve_segment);
         }
         else
         {
            Float64 k;
            curve->get_K1(&k);
            auto* parabolic_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along, horizontal_length, start_height, start_gradient, end_gradient, IsZero(k) ? Float64_Max : 1 / k, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
            auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, parabolic_segment);
            file.addEntity(segment);
            profile_segments->push(segment);

            // based on this post (https://forums.buildingsmart.org/t/ifcpolynomialcurve-clarification/4716), define the vertical curve as
            // CoefficientsX = {0.0, 1.0}
            // CoefficientsY = {C,B,A}
            Float64 A = (end_gradient - start_gradient) / (2 * horizontal_length);
            Float64 B = start_gradient;
            Float64 C = start_height;
            auto parent_curve_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
            auto parent_curve = new Schema::IfcPolynomialCurve(parent_curve_placement, std::vector<double>{0.0,1.0}, std::vector<double>{C, B, A}, boost::none);
            auto segment_curve_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, segment_curve_placement, 
               new Schema::IfcLengthMeasure(start_dist_along), 
               new Schema::IfcLengthMeasure(horizontal_length),
               parent_curve);
            file.addEntity(curve_segment);
            curve_segments->push(curve_segment);
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
      auto linear_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, prev_end_dist_along, length, prev_end_height, prev_end_gradient, prev_end_gradient, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
      auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, linear_segment);
      file.addEntity(segment);
      profile_segments->push(segment);

      auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
      auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(prev_end_gradient), sin(prev_end_gradient)});
      auto parent_curve = new Schema::IfcLine(vertical_point, new Schema::IfcVector(gradient_direction, 1.0));
      auto segment_curve_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
      auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, segment_curve_placement, new Schema::IfcLengthMeasure(prev_end_dist_along), new Schema::IfcLengthMeasure(length), parent_curve);
      file.addEntity(curve_segment);
      curve_segments->push(curve_segment);

      // check elevation
      ATLASSERT(IsEqual(endElevation, prev_end_height + length * prev_end_gradient));
   }

   // Add terminator segment
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Product_Shape/Product_Geometric_Representation/Alignment_Geometry/Alignment_Geometry_-_Horizontal_and_Vertical/content.html
   // 4.1.7.1.1.2 Zero length segment shall be added at the end of the list of segments
   auto terminator_line_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, 0.0, 0, 0, 0, 0, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
   auto terminator_segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, nullptr, nullptr, terminator_line_segment);
   file.addEntity(terminator_segment);
   profile_segments->push(terminator_segment);

   // 4.1.7.1.1.2 If the geometry definition is present, then a zero length curve segment must be provided as well
   auto terminator_vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
   auto terminator_gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(prev_end_gradient), sin(prev_end_gradient)});
   auto terminator_parent_curve = new Schema::IfcLine(terminator_vertical_point, new Schema::IfcVector(terminator_gradient_direction, 1.0));
   auto terminator_placement = new Schema::IfcAxis2Placement2D(new Schema::IfcCartesianPoint(std::vector<double>{0.0, 0.0}), new Schema::IfcDirection(std::vector<double>{1.0, 0.0}));
   auto terminator_curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_DISCONTINUOUS, terminator_placement, new Schema::IfcLengthMeasure(0.0), new Schema::IfcLengthMeasure(0.0), terminator_parent_curve);
   file.addEntity(terminator_curve_segment);
   curve_segments->push(terminator_curve_segment);

   auto vertical_profile = new Schema::IfcAlignmentVertical(IfcParse::IfcGlobalId(), nullptr, std::string("Vertical Alignment"), boost::none, boost::none, file.getSingle<typename Schema::IfcLocalPlacement>(), nullptr);
   file.addEntity(vertical_profile);

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Nests vertical alignment segments with vertical alignment"), vertical_profile, profile_segments);
   file.addEntity(nests);

   *pvertical_profile = vertical_profile;

   
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

   if (options.alignment_model == CIfcModelBuilderOptions::AlignmentModel::Polyline)
   {
      // Instead of IfcGradientCurve, we are using a generalized 3D polyline geometric representation of the alignment (a 3D wire)
      // This isn't as accurate, but some viewer may be able to deal with this better
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
      auto polyline = new Schema::IfcPolyline(points);
      *palignment_gradient_curve = nullptr;
   }
   else
   {
      auto gradient_curve = new Schema::IfcGradientCurve(curve_segments, false, horizontal_geometry_base_curve, nullptr);
      file.addEntity(gradient_curve);
      *palignment_gradient_curve = gradient_curve;
   }
}

template <typename Schema>
void CreateAlignment(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;
   typename Schema::IfcAlignmentHorizontal* horizontal_alignment_layout;
   typename Schema::IfcCompositeCurve* composite_curve;
   CreateHorizontalAlignment<Schema>(file, pBroker, options, &horizontal_alignment_layout, &composite_curve);

   typename Schema::IfcAlignmentVertical* vertical_profile_layout;
   typename Schema::IfcGradientCurve* gradient_curve;
   CreateVerticalProfile<Schema>(file, pBroker, composite_curve, options, &vertical_profile_layout, &gradient_curve);


   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);

   auto axis_model_representation_subcontext = new Schema::IfcGeometricRepresentationSubContext(std::string("Axis"), std::string("Model"), geometric_representation_context, boost::none, Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_MODEL_VIEW, boost::none);
   file.addEntity(axis_model_representation_subcontext);

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr horizontal_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   horizontal_representation_items->push(composite_curve);

   auto footprint_shape_representation = new Schema::IfcShapeRepresentation(axis_model_representation_subcontext, std::string("FootPrint"), std::string("Curve2D"), horizontal_representation_items);
   file.addEntity(footprint_shape_representation);

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr vertical_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   vertical_representation_items->push(gradient_curve);

   auto curve3d_shape_representation = new Schema::IfcShapeRepresentation(axis_model_representation_subcontext, std::string("Axis"), std::string("Curve3D"), vertical_representation_items);
   file.addEntity(curve3d_shape_representation);


   typename aggregate_of<typename Schema::IfcRepresentation>::ptr representations(new aggregate_of<typename Schema::IfcRepresentation>());
   representations->push(footprint_shape_representation); // 2D alignment geometry
   representations->push(curve3d_shape_representation); // 3D alignment geometry
   auto alignment_product = new Schema::IfcProductDefinitionShape(std::string("Alignment Product Definition Shape"), boost::none, representations);


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
   auto alignment = new Schema::IfcAlignment(IfcParse::IfcGlobalId(), nullptr, strAlignmentName, boost::none, boost::none, local_placement, alignment_product, boost::none);
   file.addEntity(alignment);

   // 4.1.4.4.1 Alignments nest horizontal and vertical layouts
   // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/concepts/Object_Composition/Nesting/Alignment_Layouts/content.html
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr alignment_layout_list(new aggregate_of<typename Schema::IfcObjectDefinition>());
   alignment_layout_list->push(horizontal_alignment_layout);
   alignment_layout_list->push(vertical_profile_layout);

   auto nests_alignment_layouts = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, std::string("Nest horizontal and vertical alignment layouts with the alignment"), boost::none, alignment, alignment_layout_list);
   file.addEntity(nests_alignment_layouts);

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


   // add stationing information
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   auto point_on_alignment = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(0.0), boost::none, boost::none, boost::none, gradient_curve ? gradient_curve : composite_curve);
   auto relative_placement = new Schema::IfcAxis2PlacementLinear(point_on_alignment, nullptr, nullptr);
   auto referent_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, nullptr);


   typename aggregate_of<typename Schema::IfcProperty>::ptr pset_station_properties(new aggregate_of<typename Schema::IfcProperty>());
   pset_station_properties->push(new Schema::IfcPropertySingleValue(std::string("Station"), boost::none, new Schema::IfcLengthMeasure(startStation), nullptr));

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_Stationing"), boost::none, pset_station_properties);
   file.addEntity(property_set);

   auto stationing_referent = new Schema::IfcReferent(IfcParse::IfcGlobalId(), nullptr, std::string("Start of alignment station"), boost::none, boost::none, referent_placement, nullptr, Schema::IfcReferentTypeEnum::IfcReferentType_STATION);
   file.addEntity(stationing_referent);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_stationing_objects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_stationing_objects->push(stationing_referent);

   auto nests_stationing = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), nullptr, std::string("Nests Referents with station information with alignment"), boost::none, alignment, related_stationing_objects);
   file.addEntity(nests_stationing);

   auto rel_defines_by_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Relates station properties to referent"), boost::none, related_stationing_objects, property_set);
   file.addEntity(rel_defines_by_properties);
}

template <typename Schema>
typename Schema::IfcProfileDef* CreateSectionProfile(IShapes* pShapes,const pgsPointOfInterest& poi,IntervalIndexType intervalIdx, const CIfcModelBuilderOptions& options)
{
   CComPtr<IShape> shape;
   IndexType gdrIdx, slabIdx;
   pShapes->GetSegmentShape(intervalIdx, poi, true/*orient*/, pgsTypes::scBridge, &shape, &gdrIdx, &slabIdx);

   CComQIPtr<ICompositeShape> composite(shape);
   CComPtr<ICompositeShapeItem> shapeItem;
   composite->get_Item(gdrIdx, &shapeItem);

   CComPtr<IShape> gdrShape;
   CComPtr<IShape> _gdrShape;
   shapeItem->get_Shape(&_gdrShape);
   CComQIPtr<ICompositeShape> compGdrShape(_gdrShape);
   if (compGdrShape)
   {
      CComPtr<ICompositeShapeItem> compItem;
      compGdrShape->get_Item(0, &compItem);
      compItem->get_Shape(&gdrShape);
   }
   else
   {
      gdrShape = _gdrShape;
   }

   auto polyline = CreatePolyline<Schema>(gdrShape);
   auto girder_section = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, std::string("CrossSectionProfile"), polyline);

   return girder_section;
}

template <typename Schema> 
typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr CreateStrands(IfcHierarchyHelper<Schema>& file, IBroker* pBroker,const pgsPointOfInterest& poiStart,const pgsPointOfInterest& poiEnd,typename Schema::IfcObjectPlacement* strand_placement)
{
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr strands(new aggregate_of<typename Schema::IfcObjectDefinition>());

   const CSegmentKey& segmentKey(poiStart.GetSegmentKey());

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   GET_IFACE2_NOCHECK(pBroker, IMaterials, pMaterials);

   PoiList vHP;
   pPoi->GetPointsOfInterest(segmentKey, POI_HARPINGPOINT, &vHP);
   std::array<std::string, 3> strStrandType{ "Straight","Harped","Temporary" };
   for (int i = 0; i < 3; i++)
   {
      pgsTypes::StrandType strandType = pgsTypes::StrandType(i);

      StrandIndexType nStrands = pStrandGeom->GetStrandCount(segmentKey, strandType);
      if (nStrands == 0) continue;

      const auto* pStrand = pMaterials->GetStrandMaterial(segmentKey, strandType);

      CComPtr<IPoint2dCollection> strand_points_start, strand_points_end;
      pStrandGeom->GetStrandPositions(poiStart, strandType, &strand_points_start);
      pStrandGeom->GetStrandPositions(poiEnd, strandType, &strand_points_end);

      std::vector<CComPtr<IPoint2dCollection>> strands_at_harp_points;
      if (strandType == pgsTypes::Harped)
      {
         for (const pgsPointOfInterest& poi : vHP)
         {
            CComPtr<IPoint2dCollection> points;
            pStrandGeom->GetStrandPositions(poi, strandType, &points);
            strands_at_harp_points.push_back(points);
         }
      }

      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr strand_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      for (StrandIndexType strandIdx = 0; strandIdx < nStrands; strandIdx++)
      {
         typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());

         CComPtr<IPoint2d> pntStart;
         strand_points_start->get_Item(strandIdx, &pntStart);

         Float64 X, Y;
         pntStart->Location(&X, &Y);
         auto start_point = new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, poiStart.GetDistFromStart()});

         points->push(start_point);

         auto begin = std::begin(strands_at_harp_points);
         auto end = std::end(strands_at_harp_points);
         for (auto iter = begin; iter != end; iter++)
         {
            auto harp_points(*iter);
            CComPtr<IPoint2d> point;
            harp_points->get_Item(strandIdx, &point);

            auto i = std::distance(begin, iter);
            const pgsPointOfInterest& poi = vHP[i];

            point->Location(&X, &Y);
            auto hp = new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, poi.GetDistFromStart()});
            points->push(hp);
         }

         CComPtr<IPoint2d> pntEnd;
         strand_points_end->get_Item(strandIdx, &pntEnd);

         pntEnd->Location(&X, &Y);
         auto end_point = new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, poiEnd.GetDistFromStart()});
         points->push(end_point);

         auto directrix = new Schema::IfcPolyline(points);
         file.addEntity(directrix);

         // NOTE: IfcSweptDiskSolid is not part of AbRV.
         auto swept_disk_solid = new Schema::IfcSweptDiskSolid(directrix, pStrand->GetNominalDiameter() / 2, boost::none, boost::none, boost::none);
         file.addEntity(swept_disk_solid);
         strand_representation_items->push(swept_disk_solid);
      }

      auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
      ATLASSERT(geometric_representation_context);
      auto strand_shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), strand_representation_items);
      typename aggregate_of<typename Schema::IfcRepresentation>::ptr strand_shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
      strand_shape_representation_list->push(strand_shape_representation);
      auto strand_product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, strand_shape_representation_list);

      auto strand = new Schema::IfcTendon(IfcParse::IfcGlobalId(), nullptr, strStrandType[strandType], boost::none, boost::none, strand_placement, strand_product_definition_shape, boost::none, boost::none,
         Schema::IfcTendonTypeEnum::IfcTendonType_STRAND,
         pStrand->GetNominalDiameter(),
         pStrand->GetNominalArea(),
         pStrandGeom->GetPjack(segmentKey, strandType),
         pStrandGeom->GetJackingStress(segmentKey, strandType),
         boost::none, boost::none, boost::none);
      file.addEntity(strand);

      strands->push(strand);
   }

   return strands;
}


CIfcModelBuilder::CIfcModelBuilder(void)
{
}

CIfcModelBuilder::~CIfcModelBuilder(void)
{
}

bool CIfcModelBuilder::BuildModel(IBroker* pBroker, const CString& strFilePath, const CIfcModelBuilderOptions& options)
{
   bool bResult = false;
   switch (options.schema)
   {
   //case Schema_4x3_rc3: bResult = BuildModel<Ifc4x3_rc3>(pBroker, strFilePath, bSimplifiedAlignment); break;
   //case Schema_4x3_rc4: bResult = BuildModel<Ifc4x3_rc4>(pBroker, strFilePath, bSimplifiedAlignment); break;
   //case CIfcModelBuilderOptions::Schema::Schema_4x3_tc1: bResult = BuildModel<Ifc4x3_tc1>(pBroker, strFilePath, options); break;
   //case CIfcModelBuilderOptions::Schema::Schema_4x3_add1: bResult = BuildModel<Ifc4x3_add1>(pBroker, strFilePath, options); break;
   case CIfcModelBuilderOptions::Schema::Schema_4x3_add2: bResult = BuildModel<Ifc4x3_add2>(pBroker, strFilePath, options); break;
   default:
      ATLASSERT(false); // is there a new schema type
   }

   return bResult;
}

template <typename Schema>
void InitializeFile(IfcHierarchyHelper<Schema>& file, IBroker* pBroker,const CString& strFilePath)
{
   USES_CONVERSION;

   int nPos = strFilePath.ReverseFind('\\');
   CString strFileName(strFilePath);
   if (nPos != -1)
   {
      strFileName = strFileName.Right(strFilePath.GetLength() - nPos - 1);
   }

   GET_IFACE2(pBroker, IProjectProperties, pProjectProperties);
   GET_IFACE2(pBroker, IVersionInfo, pVersionInfo);
   GET_IFACE2(pBroker, IDocumentType, pDocType);

   auto owner_history = file.addOwnerHistory();

   // See https://standards.buildingsmart.org/documents/Implementation/ImplementationGuide_IFCHeaderData_Version_1.0.2.pdf for details about required information
   file.header().file_name().name(T2A(strFileName)); // filename without path
   std::vector<std::string> authors;
   authors.push_back(T2A(pProjectProperties->GetEngineer()));
   file.header().file_name().author(authors);
   std::vector<std::string> organizations;
   organizations.push_back(T2A(pProjectProperties->GetCompany()));
   file.header().file_name().organization(organizations);
   //file.header().file_name().preprocessor_version(); // this is info about the toolkit we are using which is IfcOpenShell... this field is filled in by default

   std::vector<std::string> file_description;
   std::ostringstream os;
   os << "ViewDefinition[Alignment-basedReferenceView]" << std::ends;
   file_description.push_back(os.str().c_str());
   file.header().file_description().description(file_description);

   std::_tostringstream _os;
   _os << _T("BridgeLink:") << (pDocType->IsPGSuperDocument() ? _T("PGSuper") : _T("PGSplice")) << _T(" Version ") << pVersionInfo->GetVersion(true).GetBuffer() << std::ends;
   std::string strVersion(T2A(_os.str().c_str()));
   file.header().file_name().originating_system(strVersion);

   //auto project = file.addProject(); // Don't like the default units in IfcOpenShell so we have do build our own
   /////////////////////////// The following is copied from IfcHierarchyHelper<Schema>::addProject and tweaked
   typename Schema::IfcUnit::list::ptr units(new typename Schema::IfcUnit::list);

   auto* dimexp = new Schema::IfcDimensionalExponents(0, 0, 0, 0, 0, 0, 0);
   auto* unit1 = new Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
   auto* unit2 = new Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN);

   units->push(unit1);
   units->push(unit2);

   auto* unit_assignment = new Schema::IfcUnitAssignment(units);

   typename Schema::IfcRepresentationContext::list::ptr rep_contexts(new typename Schema::IfcRepresentationContext::list);
   auto* project = new Schema::IfcProject(IfcParse::IfcGlobalId(), owner_history, std::string("MyProject"), boost::none, boost::none, boost::none, boost::none, rep_contexts, unit_assignment);

   file.addEntity(dimexp);
   file.addEntity(unit1);
   file.addEntity(unit2);
   file.addEntity(unit_assignment);
   file.addEntity(project);
   ///////////////////////////////////////// end of copy from IfcHierarchyHelper<Schema>::addProject

   auto site = file.addSite(project);
   auto site_local_placement = file.getSingle<typename Schema::IfcLocalPlacement>(); // addSite creates a local placement so get it here
   
   std::string bridge_name(T2A(pProjectProperties->GetBridgeName()));
   if (bridge_name.empty()) bridge_name = "Unnamed Bridge";

   std::string project_name = bridge_name + std::string(" Project");
   project->setName(project_name);

   std::string site_name = std::string("Site of ") + bridge_name;
   site->setName(site_name);

   owner_history->OwningApplication()->setApplicationFullName(std::string(pDocType->IsPGSuperDocument() ? "BridgeLink:PGSuper" : "BridgeLink:PGSplice"));
   owner_history->OwningApplication()->setApplicationIdentifier(std::string(pDocType->IsPGSuperDocument() ? "PGSuper" : "PGSplice"));
   owner_history->OwningApplication()->setVersion(std::string(T2A(pVersionInfo->GetVersion(true))));
   owner_history->OwningApplication()->ApplicationDeveloper()->setIdentification(std::string("Washington State Department of Transportation, Bridge and Structures Office"));
   owner_history->OwningApplication()->ApplicationDeveloper()->setName(std::string("Richard Brice, PE"));
}


template <typename Schema>
void CreateGirderSegmentRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IIntervals, pIntervals);
   IntervalIndexType intervalIdx = pIntervals->GetErectSegmentInterval(segmentKey);

   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   const CPrecastSegmentData* pSegment = pIBridgeDesc->GetPrecastSegmentData(segmentKey);
   pgsTypes::SegmentVariationType variationType = pSegment->GetVariationType();

   GET_IFACE2(pBroker, IBridge, pBridge);

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   if (variationType == pgsTypes::svtParabolic)
   {
      // single parabola
      Float64 Ls = pBridge->GetSegmentLength(segmentKey);
      Float64 Lleft = pSegment->GetVariationLength(pgsTypes::sztLeftPrismatic);
      Float64 Lright = pSegment->GetVariationLength(pgsTypes::sztRightPrismatic);
      Float64 L = Ls - Lleft - Lright; // length of the non-prismatic portion of the segment
      IndexType nSections = 10; // break into nSections along the parabolic taper
      for (IndexType i = 0; i < nSections; i++)
      {
         Float64 X = Lleft + i * L / nSections;
         pgsPointOfInterest poi = pPoi->GetPointOfInterest(segmentKey, X);
         vPoi.push_back(poi);
      }
      pPoi->SortPoiList(&vPoi);
   }
   else if (variationType == pgsTypes::svtDoubleParabolic)
   {
      // double parabola
      IndexType nSections = 10; // break into nSections along the parabolic taper

      // left parabola
      Float64 Lleft = pSegment->GetVariationLength(pgsTypes::sztLeftPrismatic);
      Float64 Lt = pSegment->GetVariationLength(pgsTypes::sztLeftTapered);
      for (IndexType i = 0; i < nSections; i++)
      {
         Float64 X = Lleft + i * Lt / nSections;
         pgsPointOfInterest poi = pPoi->GetPointOfInterest(segmentKey, X);
         vPoi.push_back(poi);
      }

      // right parabola
      Float64 Ls = pBridge->GetSegmentLength(segmentKey);
      Float64 Lright = pSegment->GetVariationLength(pgsTypes::sztRightPrismatic);
      Float64 Lr = pSegment->GetVariationLength(pgsTypes::sztRightTapered);
      Lleft = Ls - Lright - Lr; // location of the left end of the right parabola
      for (IndexType i = 0; i < nSections; i++)
      {
         Float64 X = Lleft + i * Lr / nSections;
         pgsPointOfInterest poi = pPoi->GetPointOfInterest(segmentKey, X);
         vPoi.push_back(poi);
      }
      pPoi->SortPoiList(&vPoi);
   }


   const pgsPointOfInterest& poiStart(vPoi.front());
   const pgsPointOfInterest& poiEnd(vPoi.back());

   CComPtr<IDirection> segment_direction;
   pBridge->GetSegmentBearing(segmentKey, &segment_direction);
   Float64 direction;
   segment_direction->get_Value(&direction);

   Float64 segment_grade = pBridge->GetSegmentSlope(segmentKey);

   Float64 Ls = pBridge->GetSegmentLength(segmentKey);

   typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr girder_line_points(new aggregate_of<typename Schema::IfcCartesianPoint>());
   girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}));
   girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{Ls, 0, 0}));
   auto girder_line = new Schema::IfcPolyline(girder_line_points);
   file.addEntity(girder_line);


   typename aggregate_of<typename Schema::IfcProfileDef>::ptr cross_sections(new aggregate_of<typename Schema::IfcProfileDef>());

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   GET_IFACE2(pBroker, IShapes, pShapes);
   // extrusion (but this doesn't work for tapered sections like spliced girders with parabolic haunch)
   //auto iter = std::begin(vPoi);
   //pgsPointOfInterest prevPoi(*iter);
   //auto prev_girder_profile = CreateSectionProfile<Schema>(pShapes, prevPoi, intervalIdx);
   //iter++;
   //auto end = std::end(vPoi);
   //for(; iter != end; iter++)
   //{
   //   const pgsPointOfInterest& poi(*iter);
   //   auto girder_profile = CreateSectionProfile<Schema>(pShapes, poi, intervalIdx);

   //   // extrude the shape in the global X direction
   //   auto position = new Ifc4x3_rc3::IfcAxis2Placement3D(
   //      new Schema::IfcCartesianPoint(std::vector<Float64>{prevPoi.GetDistFromStart(), 0, 0}), // begin extrusion dist from start from (0,0,0)
   //      new Schema::IfcDirection(std::vector<Float64>{1, 0, 0}), // direction the Z-axis of the extrusion in the global X direction
   //      new Schema::IfcDirection(std::vector<Float64>{0, 1, 0}) // direction the X-axis of the cross section in the global Y direction
   //   );

   //   Float64 depth = poi.GetDistFromStart() - prevPoi.GetDistFromStart();
   //   auto extruded_direction = new Schema::IfcDirection(std::vector<Float64>{0, 0, 1}); // extrude in the Z-direction relative to the placement (the placement puts the Z-axis in the global X direction)
   //   //auto solid = new Schema::IfcExtrudedAreaSolid(prev_girder_profile, position, extruded_direction, depth);
   //   auto solid = new Schema::IfcExtrudedAreaSolidTapered(prev_girder_profile, position, extruded_direction, depth, girder_profile);
   //   representation_items->push(solid);

   //   prevPoi = poi;
   //   prev_girder_profile = girder_profile;
   //}

   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());

   for (const pgsPointOfInterest& poi : vPoi)
   {
      auto girder_perimeter = CreateSectionProfile<Schema>(pShapes, poi, intervalIdx, options);
      file.addEntity(girder_perimeter);
      cross_sections->push(girder_perimeter);

      auto pde = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(poi.GetDistFromStart()), boost::none, boost::none, boost::none, girder_line);
      file.addEntity(pde);

      auto lp = new Schema::IfcAxis2PlacementLinear(pde, nullptr, nullptr);
      file.addEntity(lp);

      cross_section_positions->push(lp);
   }

   auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions);
   file.addEntity(sectioned_solid);
   representation_items->push(sectioned_solid);

   Float64 startSegmentStation, startSegmentOffset;
   pBridge->GetStationAndOffset(poiStart, &startSegmentStation, &startSegmentOffset);

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file);

   // place the segment relative to the alignment
   auto distance_along = startSegmentStation - startStation;
   auto segment_origin_point = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(distance_along), 0.0, 0.0, 0.0, directrix);
   auto relative_placement = new Schema::IfcAxis2PlacementLinear(
      segment_origin_point,
      new Schema::IfcDirection(std::vector<double>{ -segment_grade * cos(direction), -segment_grade * sin(direction), 1.0}), // Exact Z direction
      new Schema::IfcDirection(std::vector<double>{cos(direction), sin(direction), 0}) // RefDirection
      );
   auto segment_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, nullptr);
   file.addEntity(segment_origin_point);
   file.addEntity(relative_placement);
   file.addEntity(segment_placement);

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);


   segment->setObjectPlacement(segment_placement);
   segment->setRepresentation(product_definition_shape);
}

template <typename Schema>
void CreateGirderSegmentMaterials(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment)
{
   USES_CONVERSION;
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   GET_IFACE2(pBroker, IMaterials, pMaterials);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   IntervalIndexType releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
   IntervalIndexType liftingIntervalIdx = pIntervals->GetLiftSegmentInterval(segmentKey);
   IntervalIndexType haulingIntervalIdx = pIntervals->GetHaulSegmentInterval(segmentKey);

   // create the material
   auto material = new Schema::IfcMaterial("Precast Segment Concrete", boost::none/*description*/, boost::none/*category*/);
   file.addEntity(material);

   // Pset_MaterialConcrete
   typename aggregate_of<typename Schema::IfcProperty>::ptr material_concrete_properties(new aggregate_of<typename Schema::IfcProperty>());
   material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("CompressiveStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc28(segmentKey)), nullptr));
   material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("MaxAggregateSize"), boost::none, new Schema::IfcPositiveLengthMeasure(pMaterials->GetSegmentMaxAggrSize(segmentKey)), nullptr));
   auto pset_material_concrete = new Schema::IfcMaterialProperties(std::string("Pset_MaterialConcrete"), boost::none/*description*/, material_concrete_properties, material);
   file.addEntity(pset_material_concrete);

   // Pset_PrecastConcreteElementGeneral
   typename aggregate_of<typename Schema::IfcProperty>::ptr precast_concrete_properties(new aggregate_of<typename Schema::IfcProperty>());
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("FormStrippingStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx)), nullptr));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("LiftingStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, liftingIntervalIdx)), nullptr));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("ReleaseStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx)), nullptr));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("TransportationStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, haulingIntervalIdx)), nullptr));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("InitialTension"), boost::none, new Schema::IfcPressureMeasure(pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent)), nullptr));
   auto pset_material_precast_concrete = new Schema::IfcMaterialProperties(std::string("Pset_PrecastConcreteElementGeneral"), boost::none, precast_concrete_properties, material);
   file.addEntity(pset_material_precast_concrete);

   // need a list of entities that are associated with this material
   // right now we are creating a unique material for each segment but we still need the list
   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr segments(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   segments->push(segment);

   // associate the material with the segment (ie segments collection)
   auto rel_associates_materials = new Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Concrete_To_Precast_Segment"), boost::none, segments, material);
   file.addEntity(rel_associates_materials);
}

//template <typename Schema>
//void CreateStrandRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, typename Schema::IfcTendon* tendon)
//{
//   // place strands relative to the segment
//   auto strand_placement = file.addLocalPlacement(segment->ObjectPlacement(),
//      0, 0, 0, // (0,0,0) of the strands is at (0,0,0) of the segment
//      1, 0, 0, // direction the Z-axis of the extrusion in the global X direction 
//      0, 1, 0 // direction the X-axis of the cross section in the global Y direction
//   );
//
//   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
//   PoiList vPoi;
//   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
//   ATLASSERT(2 <= vPoi.size());
//
//   const pgsPointOfInterest& poiStart(vPoi.front());
//   const pgsPointOfInterest& poiEnd(vPoi.back());
//
//   auto strands = CreateStrands<Schema>(file, pBroker, poiStart, poiEnd, strand_placement);
//
//   if (0 < strands->size())
//   {
//      auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Segment_Aggregates_Strands"), boost::none, segment, strands);
//      file.addEntity(rel_aggregates);
//   }
//}

template <typename Schema>
void CreateClosureJointRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CClosureKey& closureKey, typename Schema::IfcElementAssembly* closureJoint, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IIntervals, pIntervals);
   IntervalIndexType intervalIdx = pIntervals->GetCompositeClosureJointInterval(closureKey);

   CSegmentKey prevSegmentKey(closureKey);
   CSegmentKey nextSegmentKey(closureKey.groupIndex, closureKey.girderIndex, closureKey.segmentIndex + 1);

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPrevSegmentPoi;
   pPoi->GetPointsOfInterest(prevSegmentKey, POI_END_FACE, &vPrevSegmentPoi);
   ATLASSERT(vPrevSegmentPoi.size() == 1);

   PoiList vNextSegmentPoi;
   pPoi->GetPointsOfInterest(nextSegmentKey, POI_START_FACE, &vNextSegmentPoi);
   ATLASSERT(vNextSegmentPoi.size() == 1);

   const pgsPointOfInterest& poiStart(vPrevSegmentPoi.front());
   const pgsPointOfInterest& poiEnd(vNextSegmentPoi.front());

   GET_IFACE2(pBroker, IBridge, pBridge);
   CComPtr<IPoint2d> pntStart, pntEnd;
   pBridge->GetPoint(poiStart, pgsTypes::pcGlobal, &pntStart);
   pBridge->GetPoint(poiEnd, pgsTypes::pcGlobal, &pntEnd);

   GET_IFACE2(pBroker, IGirder, pGirder);
   Float64 Zs = pGirder->GetTopGirderChordElevation(poiStart);
   Float64 Ze = pGirder->GetTopGirderChordElevation(poiEnd);


   Float64 distance;
   CComPtr<IDirection> closure_direction;
   cogoUtil::Inverse(pntStart, pntEnd, &distance, &closure_direction);

   Float64 direction;
   closure_direction->get_Value(&direction);

   Float64 Lc = pBridge->GetClosureJointLength(closureKey);
   Float64 closure_grade = (Ze - Zs) / Lc;

   GET_IFACE2(pBroker, IShapes, pShapes);
   typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr girder_line_points(new aggregate_of<typename Schema::IfcCartesianPoint>());
   girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}));
   girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{Lc, 0, 0}));
   auto girder_line = new Schema::IfcPolyline(girder_line_points);
   file.addEntity(girder_line);

   typename aggregate_of<typename Schema::IfcProfileDef>::ptr cross_sections(new aggregate_of<typename Schema::IfcProfileDef>());
   cross_sections->push(CreateSectionProfile<Schema>(pShapes, poiStart, intervalIdx, options));
   cross_sections->push(CreateSectionProfile<Schema>(pShapes, poiEnd, intervalIdx, options));

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());

   std::string representation_type;
   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());
   auto pde_start = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(0.0), boost::none, boost::none, boost::none, girder_line);
   auto pde_end = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(Lc), boost::none, boost::none, boost::none, girder_line);
   auto start_section = new Schema::IfcAxis2PlacementLinear(pde_start, nullptr, nullptr);
   auto end_section = new Schema::IfcAxis2PlacementLinear(pde_end, nullptr, nullptr);
   cross_section_positions->push(start_section);
   cross_section_positions->push(end_section);

   representation_type = "AdvancedSweptSolid";
   auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions);

   representation_items->push(sectioned_solid);

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), representation_type, representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);


   Float64 startClosureStation, startClosureOffset;
   pBridge->GetStationAndOffset(poiStart, &startClosureStation, &startClosureOffset);

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startClosureElevation = pAlignment->GetElevation(startClosureStation, 0.0);

   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file);

   // place closure joint relative to alignment
   auto closure_origin_point = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(startClosureStation - startStation), -startClosureOffset, Zs - startClosureElevation, 0.0, directrix);
   auto relative_placement = new Schema::IfcAxis2PlacementLinear(closure_origin_point, new Schema::IfcDirection(std::vector<double>{ -closure_grade * cos(direction), -closure_grade * sin(direction), 1}), new Schema::IfcDirection(std::vector<double>{cos(direction), sin(direction), 0}));
   auto closure_placement = new Schema::IfcLinearPlacement(nullptr/*alignment->ObjectPlacement()*/, relative_placement, nullptr);

   closureJoint->setObjectPlacement(closure_placement);
   closureJoint->setRepresentation(product_definition_shape);
}

template <typename Schema>
void CreateDeckRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcBridgePart* deck, const CIfcModelBuilderOptions& options)
{
   GET_IFACE2(pBroker, IBridge, pBridge);

   if (pBridge->GetDeckType() == pgsTypes::sdtNone)
      return;

   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file);

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   IndexType nDeckSections = 10;

   GET_IFACE2(pBroker, IShapes, pShapes);
   typename aggregate_of<typename Schema::IfcProfileDef>::ptr cross_sections(new aggregate_of<typename Schema::IfcProfileDef>());
   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());
   for (IndexType i = 0; i <= nDeckSections; i++)
   {
      auto station = i * (endBrgStation - startBrgStation) / nDeckSections + startBrgStation;

      CComPtr<IShape> slab_shape;
      pShapes->GetSlabShape(station, nullptr, true/*include haunch*/, &slab_shape);

      auto polyline = CreatePolyline<Schema>(slab_shape);
      auto deck_perimeter = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, std::string("deck perimeter"), polyline);
      cross_sections->push(deck_perimeter);
      file.addEntity(deck_perimeter);

      auto pde = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(station - startStation), boost::none, boost::none, boost::none, directrix);
      auto deck_section_placement = new Schema::IfcAxis2PlacementLinear(pde, nullptr, nullptr);
      cross_section_positions->push(deck_section_placement);
      file.addEntity(pde);
      file.addEntity(deck_section_placement);
   }

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(directrix, cross_sections, cross_section_positions);
   representation_items->push(sectioned_solid);
   file.addEntity(sectioned_solid);

   auto site = file.getSingle<typename Schema::IfcSite>();
   auto deck_placement = site->ObjectPlacement();

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   file.addEntity(product_definition_shape);

   auto slab = new Schema::IfcSlab(IfcParse::IfcGlobalId(), nullptr, boost::none, std::string("Deck slab element"), boost::none, deck_placement, product_definition_shape, boost::none,
      Schema::IfcSlabTypeEnum::IfcSlabType_FLOOR); // see Ifc 4x3 6.1.2.19.2 (FLOOR represents a bridge deck)
   typename aggregate_of<typename Schema::IfcProduct>::ptr list_of_slabs(new aggregate_of<typename Schema::IfcProduct>());
   list_of_slabs->push(slab);

   auto rel_slab_contained_in_deck = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, std::string("Places slab into spatial structure of deck"), boost::none, list_of_slabs, deck);
   file.addEntity(rel_slab_contained_in_deck);
}


template <typename Schema>
void CreateRailingSystemRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, pgsTypes::TrafficBarrierOrientation tbOrientation, typename Schema::IfcRailing* railing, const CIfcModelBuilderOptions& options)
{
   GET_IFACE2(pBroker, IBarriers, pBarriers);

   bool bHasSidewalk = pBarriers->HasSidewalk(tbOrientation);
   bool bHasInteriorBarrier = pBarriers->HasInteriorBarrier(tbOrientation);

   IndexType nShapes = 1 + (bHasSidewalk ? 1 : 0) + (bHasInteriorBarrier ? 1 : 0);

   GET_IFACE2(pBroker, IBridge, pBridge);
   //Float64 startStation = pBridge->GetPierStation(0);
   //Float64 endStation = pBridge->GetPierStation(pBridge->GetPierCount() - 1);
   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);

   GET_IFACE2(pBroker, IShapes, pShapes);
   CComPtr<IShape> start_barrier_shape, end_barrier_shape;
   if (tbOrientation == pgsTypes::tboLeft)
   {
      pShapes->GetLeftTrafficBarrierShape(startBrgStation, nullptr, &start_barrier_shape);
      pShapes->GetLeftTrafficBarrierShape(endBrgStation, nullptr, &end_barrier_shape);
   }
   else
   {
      pShapes->GetRightTrafficBarrierShape(startBrgStation, nullptr, &start_barrier_shape);
      pShapes->GetRightTrafficBarrierShape(endBrgStation, nullptr, &end_barrier_shape);
   }

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file);

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());
   std::vector<typename aggregate_of<typename Schema::IfcProfileDef>::ptr> cross_sections;
   for (int i = 0; i < nShapes; i++)
   {
      typename aggregate_of<typename Schema::IfcProfileDef>::ptr ptr(new aggregate_of<typename Schema::IfcProfileDef>());
      cross_sections.push_back(ptr);
   }

   for (int i = 0; i < 2; i++)
   {
      CComPtr<IShape> barrier_shape = (i == 0) ? start_barrier_shape : end_barrier_shape;
      CComQIPtr<ICompositeShape> composite(barrier_shape);
      if (!composite)
      {
         CComPtr<ICompositeShape> compShape;
         compShape.CoCreateInstance(CLSID_CompositeShape);
         compShape->AddShape(barrier_shape, VARIANT_FALSE);
         composite = compShape;
      }

#if defined _DEBUG
      IndexType _nShapes;
      composite->get_Count(&_nShapes);
      ATLASSERT(nShapes == _nShapes); // if this fires the actual number of shapes is not the same as the expected number of shapes
#endif

      auto distance_along = (i == 0 ? startBrgStation - startStation : endBrgStation - startStation);
      auto pde = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(distance_along), boost::none, boost::none, boost::none, directrix);
      auto placement = new Schema::IfcAxis2PlacementLinear(pde, nullptr, nullptr);
      cross_section_positions->push(placement);
      file.addEntity(pde);
      file.addEntity(placement);

      for (IndexType shapeIdx = 0; shapeIdx < nShapes; shapeIdx++)
      {
         CComPtr<ICompositeShapeItem> shape_item;
         composite->get_Item(shapeIdx, &shape_item);

         CComPtr<IShape> shape;
         shape_item->get_Shape(&shape);

         auto polyline = CreatePolyline<Schema>(shape);
         auto shape_perimeter = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, boost::none, polyline);
         cross_sections[shapeIdx]->push(shape_perimeter);
      } // next shape
   } // next end (i)

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   for (IndexType shapeIdx = 0; shapeIdx < nShapes; shapeIdx++)
   {
      auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(directrix, cross_sections[shapeIdx], cross_section_positions);
      representation_items->push(sectioned_solid);
      file.addEntity(sectioned_solid);
   }

   auto site = file.getSingle<typename Schema::IfcSite>();
   auto railing_placement = site->ObjectPlacement();


   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   file.addEntity(product_definition_shape);

   railing->setObjectPlacement(railing_placement);
   railing->setRepresentation(product_definition_shape);
}

template <typename Schema>
typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr CreatePiers(IfcHierarchyHelper<Schema>& file, IBroker* pBroker)
{
   USES_CONVERSION;

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_piers(new aggregate_of<typename Schema::IfcObjectDefinition>());

   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nPiers = pBridge->GetPierCount();
   for (IndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      std::string pier_name(T2A(LABEL_PIER_EX(pBridge->IsAbutment(pierIdx), pierIdx)));
      auto pier = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, pier_name, boost::none, boost::none, nullptr, nullptr, boost::none,
         Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
         Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
         pBridge->IsAbutment(pierIdx) ? Schema::IfcBridgePartTypeEnum::IfcBridgePartType_ABUTMENT : Schema::IfcBridgePartTypeEnum::IfcBridgePartType_PIER);
      file.addEntity(pier);

      auto foundation = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none,nullptr,nullptr,boost::none,
         Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
         Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
         Schema::IfcBridgePartTypeEnum::IfcBridgePartType_FOUNDATION);
      file.addEntity(foundation);

      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_foundations(new aggregate_of<typename Schema::IfcObjectDefinition>());
      list_of_foundations->push(foundation);

      // This was the TPF concept, but it violates bSI rules
      //// IfcBridgePart::PIER <-> IfcRelContainedInSpatialStructure <-> IfcBridgePart::FOUNDATION
      //auto rel_contained_in_pier_spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, std::string("Foundations in pier spatial structure"), boost::none, list_of_foundations, pier);
      //file.addEntity(rel_contained_in_pier_spatial_structure);

      // IfcBridgePart::PIER <-> IfcRelAggregates <-> IfcBridgePart::FOUNDATION
      auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Foundation is an aggregate component of pier"), boost::none, pier, list_of_foundations);
      file.addEntity(rel_aggregates);


      list_of_piers->push(pier);
   }

   return list_of_piers;
}

template <typename Schema>
void CreateBridge(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   // Define spatial structure
   // From IfcSite https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/IfcSite.htm
   // IfcProject <-> IfcRelAggregates <-> IfcSite <-> IfcRelAggregates <-> IfcBridge
   // InitializeFile sets up IfcProject <-> IfcRelAggregates <-> IfcSite



   // Create IfcBridge, which is an IfcSpatialStructureElement
   GET_IFACE2(pBroker, IProjectProperties, pProjectProperties);
   std::string bridge_name(T2A(pProjectProperties->GetBridgeName()));
   if (bridge_name.empty()) bridge_name = "Unnamed Bridge";

   auto bridge = new Schema::IfcBridge(IfcParse::IfcGlobalId(), nullptr, bridge_name, boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcElementCompositionEnum::IfcElementComposition_COMPLEX, Schema::IfcBridgeTypeEnum::IfcBridgeType_GIRDER);
   file.addEntity(bridge);

   // create a list of bridges in the site
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_bridges_in_the_site(new aggregate_of<typename Schema::IfcObjectDefinition>());
   list_of_bridges_in_the_site->push(bridge); // add the bridge to the list

   // aggregate the bridges with the side
   // IfcSite <-> IfcRelAggregates <-> IfcBridge
   auto site = file.getSingle<typename Schema::IfcSite>();
   auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Bridges in the site"), boost::none, site, list_of_bridges_in_the_site);
   file.addEntity(rel_aggregates);



   // Create top level spatial structure of bridge
   auto superstructure = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Superstructure"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);
   file.addEntity(superstructure);

   auto substructure = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Substructure"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUBSTRUCTURE);
   file.addEntity(substructure);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_bridge_parts(new aggregate_of<typename Schema::IfcObjectDefinition>());
   list_of_bridge_parts->push(superstructure);
   list_of_bridge_parts->push(substructure);

   // IfcBridge <-> IfcRelAggregates <-> IfcBridgePart::SUPERSTRUCTURE, SUBSTRUCTURE
   auto bridge_spatial_elements = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in bridge spatial structure"), boost::none, bridge, list_of_bridge_parts);
   file.addEntity(bridge_spatial_elements);

   
   
   // Create spatial structure of substructure
   // IfcBridgePart::SUBSTRUCTURE <-> IfcRelAggregates <-> IfcBridgePart::ABUTMENT, PIER
   auto list_of_piers = CreatePiers(file, pBroker);
   auto substructure_spatial_elements = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in substructure spatial structure"), boost::none, substructure, list_of_piers);
   file.addEntity(substructure_spatial_elements);

   // Create spatial structure of superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelAggregates <-> IfcBridgePart::DECK, 
   // NOTE: Could also be DECK_SEGMENT if we looked at the spliced girder staged deck construction model
   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_superstructure_spatial_elements(new aggregate_of<typename Schema::IfcObjectDefinition>());

   auto deck = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Deck"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_DECK);
   CreateDeckRepresentation(file, pBroker, deck, options);
   file.addEntity(deck);
   list_of_superstructure_spatial_elements->push(deck);

   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelAggregates <-> IfcBridgePart::DECK
   auto rel_aggregates_elements_of_superstructure_spatial_structure = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in superstructure spatial structure"), boost::none, superstructure, list_of_superstructure_spatial_elements);
   file.addEntity(rel_aggregates_elements_of_superstructure_spatial_structure);



   // Add railings to the spatial structure of the superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcRailing
   typename aggregate_of<typename Schema::IfcProduct>::ptr list_of_superstructure_elements(new aggregate_of<typename Schema::IfcProduct>());
   auto left_railing = new Schema::IfcRailing(IfcParse::IfcGlobalId(), nullptr, std::string("Left Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_GUARDRAIL);
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboLeft, left_railing, options);
   file.addEntity(left_railing);
   list_of_superstructure_elements->push(left_railing);

   auto right_railing = new Schema::IfcRailing(IfcParse::IfcGlobalId(), nullptr, std::string("Right Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_GUARDRAIL);
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboRight, right_railing, options);
   file.addEntity(right_railing);
   list_of_superstructure_elements->push(right_railing);

   // Add girders to the spatial structure of the superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcElementAssembly::GIRDER
   GET_IFACE2(pBroker, IBridge, pBridge);
   GroupIndexType nGroups = pBridge->GetGirderGroupCount();
   for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
   {
      GirderIndexType nGirders = pBridge->GetGirderCount(grpIdx);
      for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
      {
         std::_tostringstream os;
         os << GIRDER_LABEL(CGirderKey(grpIdx, gdrIdx));
         std::string girder_name(T2A(os.str().c_str()));

         SegmentIndexType nSegments = pBridge->GetSegmentCount(grpIdx, gdrIdx);

         auto girder = new Schema::IfcElementAssembly(IfcParse::IfcGlobalId(), nullptr, girder_name, boost::none, boost::none, nullptr, nullptr, boost::none, boost::none, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
         file.addEntity(girder);
         list_of_superstructure_elements->push(girder);

         typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_girder_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
         for (SegmentIndexType segIdx = 0; segIdx < nSegments; segIdx++)
         {
            CSegmentKey segmentKey(grpIdx, gdrIdx, segIdx);
            std::ostringstream os_segment_name;
            os_segment_name << "Segment " << LABEL_SEGMENT(segIdx);
            auto segment_name = os_segment_name.str();
            auto segment = new Schema::IfcBeam(IfcParse::IfcGlobalId(), nullptr, segment_name, boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcBeamTypeEnum::IfcBeamType_GIRDER_SEGMENT);
            CreateGirderSegmentRepresentation<Schema>(file, pBroker, segmentKey, segment, options);
            CreateGirderSegmentMaterials<Schema>(file, pBroker, segmentKey, segment);
            file.addEntity(segment);
            list_of_girder_segments->push(segment);

            // build segment internals (strand, rebar, etc)... need to do this for each strand, bar, etc
            // Save this for later
            //typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_segment_parts(new aggregate_of<typename Schema::IfcObjectDefinition>());

            //auto strand = new Schema::IfcTendon(IfcParse::IfcGlobalId(), nullptr, std::string("strand"), boost::none, boost::none, nullptr, nullptr, boost::none, boost::none /*grade*/,
            //   Schema::IfcTendonTypeEnum::IfcTendonType_STRAND, boost::none /*nominal diameter*/, boost::none /*area*/, boost::none/*force*/, boost::none/*prestress*/,
            //   boost::none/*friction coefficient*/, boost::none /*anchorage slip*/, boost::none/*min curvature radius*/);
            //file.addEntity(strand);
            //list_of_segment_parts->push(strand);

            //auto longitudinal_rebar = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("longitudinal rebar"), boost::none,
            //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade (depreciated)*/, boost::none/*nominal diameter (depreciated)*/,
            //   boost::none /*nominal area*/, boost::none/*bar length (depreciated)*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_MAIN, boost::none/*bar surface (depreciated)*/);
            //file.addEntity(longitudinal_rebar);
            //list_of_segment_parts->push(longitudinal_rebar);

            //auto stirrups = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("stirrups"), boost::none,
            //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade (depreciated)*/, boost::none/*nominal diameter (depreciated)*/,
            //   boost::none /*nominal area*/, boost::none/*bar length (depreciated)*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_SHEAR, boost::none/*bar surface (depreciated)*/);
            //file.addEntity(stirrups);
            //list_of_segment_parts->push(stirrups);

            //auto lifting_loops = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("lifting loops"), boost::none,
            //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade (depreciated)*/, boost::none/*nominal diameter (depreciated)*/,
            //   boost::none /*nominal area*/, boost::none/*bar length (depreciated)*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_LIGATURE, boost::none/*bar surface (depreciated)*/);
            //file.addEntity(lifting_loops);
            //list_of_segment_parts->push(lifting_loops);

            //auto segment_parts_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("precast segment parts"), boost::none, segment, list_of_segment_parts);
            //file.addEntity(segment_parts_aggregates);


            if (segIdx < nSegments - 1)
            {
               std::ostringstream os_closure_name;
               os_closure_name << "Closure Joint " << LABEL_SEGMENT(segIdx);
               auto closure_joint_name = os_closure_name.str();
               auto closure_joint = new Schema::IfcElementAssembly(IfcParse::IfcGlobalId(), nullptr, closure_joint_name, boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_SITE, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_USERDEFINED);
               CreateClosureJointRepresentation<Schema>(file, pBroker, segmentKey, closure_joint, options);
               file.addEntity(closure_joint);
               list_of_girder_segments->push(closure_joint);

               // need to do this for each bar, stirrup, etc
               //typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_closure_joint_parts(new aggregate_of<typename Schema::IfcObjectDefinition>());

               //auto longitudinal_rebar = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("longitudinal rebar"), boost::none,
               //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade*/, boost::none/*nominal diameter*/,
               //   boost::none /*nominal area*/, boost::none/*bar length*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_MAIN, Schema::IfcReinforcingBarSurfaceEnum::IfcReinforcingBarSurface_TEXTURED);
               //file.addEntity(longitudinal_rebar);
               //list_of_closure_joint_parts->push(longitudinal_rebar);

               //auto stirrups = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("stirrups"), boost::none,
               //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade*/, boost::none/*nominal diameter*/,
               //   boost::none /*nominal area*/, boost::none/*bar length*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_SHEAR, Schema::IfcReinforcingBarSurfaceEnum::IfcReinforcingBarSurface_TEXTURED);
               //file.addEntity(stirrups);
               //list_of_closure_joint_parts->push(stirrups);

               //auto closure_joint_parts_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("cast in place closure joint parts"), boost::none, closure_joint, list_of_closure_joint_parts);
               //file.addEntity(closure_joint_parts_aggregates);
            }

            std::ostringstream os_relationship_name;
            os_relationship_name << "Elements of girder for Group " << LABEL_GROUP(grpIdx) << " Girder " << T2A(LABEL_GIRDER(gdrIdx));
            auto girder_aggregation_name = os_relationship_name.str();
            auto girder_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, girder_aggregation_name, boost::none, girder, list_of_girder_segments);
            file.addEntity(girder_aggregates);
         }
      }
   }

   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcRailing, IfcElementAssembly::GIRDER
   auto rel_contained_in_superstructure_spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in superstructure spatial structure"), boost::none, list_of_superstructure_elements, superstructure);
   file.addEntity(rel_contained_in_superstructure_spatial_structure);
}

template <typename Schema>
bool CIfcModelBuilder::BuildModel(IBroker* pBroker, const CString& strFilePath, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   IfcHierarchyHelper<Schema> file;
   InitializeFile<Schema>(file, pBroker, strFilePath); // creates project and site

   CreateAlignment<Schema>(file, pBroker, options); // creates alignment and aggregates with project, references into site spatial structure

   if (options.model_elements == CIfcModelBuilderOptions::ModelElements::AlignmentAndBridge)
   {
      CreateBridge<Schema>(file, pBroker, options); // creates bridge with site spatial structure
   }

   std::ofstream ofs(T2A(strFilePath));
   ofs << file;

   return true;
}
