#include "stdafx.h"
#include "IfcModelBuilder.h"

#include <IFace\Project.h>
#include <IFace\VersionInfo.h>
#include <IFace\Alignment.h>
#include <IFace\DocumentType.h>
#include <IFace\Bridge.h>
#include <IFace\Intervals.h>
#include <IFace\PrestressForce.h>

#include <PgsExt\GirderLabel.h>
#include <PgsExt\PrecastSegmentData.h>

#include <WBFLCogo\CogoHelpers.h>




#define CLOCKWISE 0
#define COUNTERCLOCKWISE 1
int GetVertexOrdering(IShape* pShape)
{
    Float64 area;

    // Intialize and check for null polygon.
    area = 0;

    CComPtr<IPoint2dCollection> points;
    pShape->get_PolyPoints(&points);

    CollectionIndexType cPoints;
    points->get_Count(&cPoints);

    if (cPoints == 0)
    {
        return CLOCKWISE;
    }

    Float64 x0, y0;
    Float64 x1, y1;
    Float64 dy, dx;
    Float64 ar, at;

    // loop over all points - make sure of closure
    CollectionIndexType idx0, idx1;
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
void CreateHorizontalAlignment_4x3(IfcHierarchyHelper<Schema>& file,IBroker* pBroker, typename Schema::IfcAlignmentHorizontal** phorizontal_alignment, typename Schema::IfcCompositeCurve** phorizontal_geometry_base_curve)
{
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> alignment_segments(new IfcTemplatedEntityList<Schema::IfcObjectDefinition>());
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcSegment>> curve_segments(new IfcTemplatedEntityList<Schema::IfcSegment>());

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
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, ifc_line_segment);
         file.addEntity(segment);
         alignment_segments->push(segment);

         boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
         points->push(ifc_prev_point);
         points->push(ConvertPoint<Schema>(pntTS));
         auto parent_curve = new Schema::IfcPolyline(points);
         auto placement = new Schema::IfcAxis2Placement2D(ifc_prev_point, new Schema::IfcDirection(std::vector<double>{cos(angle), sin(angle)}));
         auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(0.0), new Schema::IfcNonNegativeLengthMeasure(dist), parent_curve);
         curve_segments->push(curve_segment);
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
      auto ifc_start = ConvertPoint<Schema>(pntTS);

      if (0.0 < Lspiral[spEntry])
      {
         // there is an entry spiral
         CComPtr<IDirection> bkTangentBrg;
         curve->get_BkTangentBrg(&bkTangentBrg);
         Float64 bk_tangent_direction;
         bkTangentBrg->get_Value(&bk_tangent_direction);

         auto ifc_entry_spiral = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_start, bk_tangent_direction, 0, (bIsCCW ? 1.0 : -1.0)*R, Lspiral[spEntry], boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, ifc_entry_spiral);
         file.addEntity(segment);
         alignment_segments->push(segment);

         Float64 A = (bIsCCW ? 1.0 : -1.0) * sqrt(Lspiral[spEntry] / R);
         auto parent_curve = new Schema::IfcClothoid(ifc_start, A);
         auto placement = new Schema::IfcAxis2Placement2D(ifc_start, new Schema::IfcDirection(std::vector<double>{cos(bk_tangent_direction), sin(bk_tangent_direction)}));
         auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(0.0), new Schema::IfcNonNegativeLengthMeasure(Lspiral[spEntry]), parent_curve);
         curve_segments->push(curve_segment);

         // spiral ends at the Spiral to Curve point (CS)
         CComPtr<IPoint2d> pntSC;
         curve->get_SC(&pntSC);
         ifc_start = ConvertPoint<Schema>(pntSC);
      }

      // build the horizontal curve
      CComPtr<IDirection> bkTangentBrgCurve;
      curve->get_CurveBkTangentBrg(&bkTangentBrgCurve);
      Float64 bk_tangent_direction_curve;
      bkTangentBrgCurve->get_Value(&bk_tangent_direction_curve);
      auto ifc_hcurve = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_start, bk_tangent_direction_curve, (bIsCCW ? 1.0 : -1.0)*R, (bIsCCW ? 1.0 : -1.0)*R, Lc, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC);
      auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, ifc_hcurve);
      file.addEntity(segment);
      alignment_segments->push(segment);

      CComPtr<IPoint2d> sc, ccc, cs;
      curve->get_SC(&sc);
      curve->get_CCC(&ccc);
      curve->get_CS(&cs);
      auto ifc_ccc = ConvertPoint<Schema>(ccc);
      auto parent_curve = new Schema::IfcCircle(ifc_ccc, (bIsCCW ? 1.0 : -1.0) * R); // radius is always a positive value but CW curves don't render correctly in ACCA IFC BIM browser if the radius isn't negative
      auto placement = new Schema::IfcAxis2Placement2D(ifc_start, new Schema::IfcDirection(std::vector<double>{cos(bk_tangent_direction_curve), sin(bk_tangent_direction_curve)}));
      auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(0.0), new Schema::IfcNonNegativeLengthMeasure(Lc), parent_curve);
      curve_segments->push(curve_segment);

      if (0.0 < Lspiral[spExit])
      {
         // there is an exit spiral

         CComPtr<IDirection> fwdTangentBrgCurve;
         curve->get_CurveFwdTangentBrg(&fwdTangentBrgCurve); // forward tangent of curve is start tangent to exit spiral
         Float64 fwd_tangent_direction_curve;
         fwdTangentBrgCurve->get_Value(&fwd_tangent_direction_curve);

         // spiral starts at the Curve to Spiral point (CS)
         CComPtr<IPoint2d> pntCS;
         curve->get_CS(&pntCS);
         ifc_start = ConvertPoint<Schema>(pntCS);

         auto ifc_exit_spiral = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_start, fwd_tangent_direction_curve, (bIsCCW ? 1.0 : -1.0)*R, 0, Lspiral[spExit], boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, ifc_exit_spiral);
         file.addEntity(segment);
         alignment_segments->push(segment);

         Float64 A = (bIsCCW ? 1.0 : -1.0) * sqrt(Lspiral[spExit] / R);
         auto parent_curve = new Schema::IfcClothoid(ifc_start, A);
         auto placement = new Schema::IfcAxis2Placement2D(ifc_prev_point, new Schema::IfcDirection(std::vector<double>{cos(fwd_tangent_direction_curve), sin(fwd_tangent_direction_curve)}));
         auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(0.0), new Schema::IfcNonNegativeLengthMeasure(Lspiral[spExit]), parent_curve);
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

   if(prevPoint->SameLocation(endPoint) == S_FALSE)
   {
      // end the alignment with a line segment
      GET_IFACE2(pBroker, IGeometry, pGeometry);
      Float64 dist;
      CComPtr<IDirection> direction;
      pGeometry->Inverse(prevPoint, endPoint, &dist, &direction);
      Float64 angle;
      direction->get_Value(&angle);

      auto ifc_line_segment = new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_prev_point, angle, 0.0, 0.0, dist, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
      auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, ifc_line_segment);
      file.addEntity(segment);
      alignment_segments->push(segment);
   
      boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
      points->push(ifc_prev_point);
      points->push(ConvertPoint<Schema>(endPoint));
      auto parent_curve = new Schema::IfcPolyline(points);
      auto placement = new Schema::IfcAxis2Placement2D(ifc_prev_point, new Schema::IfcDirection(std::vector<double>{cos(angle), sin(angle)}));
      auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_DISCONTINUOUS, placement, new Schema::IfcNonNegativeLengthMeasure(0.0), new Schema::IfcNonNegativeLengthMeasure(dist), parent_curve);
      curve_segments->push(curve_segment);
   }

   // create a horizontal alignment from all the alignment segments
   auto horizontal_alignment = new Schema::IfcAlignmentHorizontal(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Horz Alignment"), boost::none, boost::none, file.getSingle<Schema::IfcLocalPlacement>(), nullptr/*representation*/);
   file.addEntity(horizontal_alignment);

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Horz Segments"), std::string("Nests horizontal alignment segments with horizontal alignment"), horizontal_alignment, alignment_segments);
   file.addEntity(nests);

   *phorizontal_alignment = horizontal_alignment;

   auto composite_curve = new Schema::IfcCompositeCurve(curve_segments, false);
   *phorizontal_geometry_base_curve = composite_curve;
}

template <typename Schema>
void CreateVerticalProfile_4x3(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcCompositeCurve* horizontal_geometry_base_curve, bool bSimplifiedAlignment, typename Schema::IfcAlignmentVertical** pvertical_profile, typename Schema::IfcProductDefinitionShape** palignment_product_definition_shape)
{
   boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcObjectDefinition>> profile_segments(new IfcTemplatedEntityList<typename Schema::IfcObjectDefinition>());
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcSegment>> curve_segments(new IfcTemplatedEntityList<Schema::IfcSegment>());

   // Profile is defined by profile segments located at "distance from start" of the alignment and "length".
   // We can't use stations to define the profile.
   // Distance from start is taken to be Station - Start Station

   GET_IFACE2(pBroker, IRoadway, pAlignment);

   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   Float64 prev_end_dist_along = 0; // startStation; // this is distance along alignment, not station
   Float64 prev_end_gradient = startGrade;
   Float64 prev_end_height = startElevation;

   IndexType nVCurves = pAlignment->GetVertCurveCount();
   for (IndexType i = 0; i < nVCurves; i++)
   {
      CComPtr<IVertCurve> curve;
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
#pragma Reminder("WORKING HERE - How to deal with station equations?") // see IfcReferent

      if (!IsEqual(prev_end_dist_along, start_dist_along))
      {
         // create a linear segment between the last profile element and this curve
         Float64 length = start_dist_along - prev_end_dist_along;
         Schema::IfcAlignmentVerticalSegment* linear_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, prev_end_dist_along, length, prev_end_height, prev_end_gradient, prev_end_gradient, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
         auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, linear_segment);
         file.addEntity(segment);
         profile_segments->push(segment);

         if (!bSimplifiedAlignment)
         {
            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
            auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(prev_end_gradient), sin(prev_end_gradient)});
            auto parent_curve = new Schema::IfcLine(vertical_point, new Schema::IfcVector(gradient_direction, 1.0));
            auto placement = new Schema::IfcAxis2Placement3D(vertical_point, new Schema::IfcDirection(std::vector<double>{0, 0, 1}), gradient_direction);
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(prev_end_dist_along), new Schema::IfcNonNegativeLengthMeasure(length), parent_curve);
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

         Schema::IfcAlignmentVerticalSegment* parabolic_segment1 = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along, l1, start_height, start_gradient, pviGrade, IsZero(k1) ? Float64_Max : 1 / k1, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
         auto segment1 = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, parabolic_segment1);
         file.addEntity(segment1);
         profile_segments->push(segment1);

         Float64 A1 = (pviGrade - start_gradient) / (2 * l1);
         Float64 B1 = start_gradient;
         Float64 C1 = start_height;
         auto vertical_point1 = new Schema::IfcCartesianPoint(std::vector<double>{0.0, start_height});
         auto gradient_direction1 = new Schema::IfcDirection(std::vector<double>{cos(start_gradient), sin(start_gradient)});
         auto placement1 = new Schema::IfcAxis2Placement3D(vertical_point1, new Schema::IfcDirection(std::vector<double>{0, 0, 1}), gradient_direction1);
         auto parent_curve1 = new Schema::IfcPolynomialCurve(placement1, boost::none, std::vector<double>{A1, B1, C1}, boost::none);
         auto curve_segment1 = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement1, new Schema::IfcNonNegativeLengthMeasure(start_dist_along), new Schema::IfcNonNegativeLengthMeasure(l1), parent_curve1);
         curve_segments->push(curve_segment1);

         Schema::IfcAlignmentVerticalSegment* parabolic_segment2 = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along + l1, l2, pviElevation, pviGrade, end_gradient, IsZero(k2) ? Float64_Max : 1 / k2, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
         auto segment2 = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, parabolic_segment2);
         file.addEntity(segment2);
         profile_segments->push(segment2);

         Float64 A2 = (end_gradient - pviGrade) / (2 * l2);
         Float64 B2 = pviGrade;
         Float64 C2 = pviElevation;
         auto vertical_point2 = new Schema::IfcCartesianPoint(std::vector<double>{0.0, pviElevation});
         auto gradient_direction2 = new Schema::IfcDirection(std::vector<double>{cos(pviGrade), sin(pviGrade)});
         auto placement2 = new Schema::IfcAxis2Placement3D(vertical_point2, new Schema::IfcDirection(std::vector<double>{0, 0, 1}), gradient_direction1);
         auto parent_curve2 = new Schema::IfcPolynomialCurve(placement2, boost::none, std::vector<double>{A2, B1, C2}, boost::none);
         auto curve_segment2 = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement2, new Schema::IfcNonNegativeLengthMeasure(start_dist_along + l1), new Schema::IfcNonNegativeLengthMeasure(l2), parent_curve2);
         curve_segments->push(curve_segment2);
      }
      else
      {
         Float64 horizontal_length;
         curve->get_Length(&horizontal_length);
         Float64 start_gradient, end_gradient;
         curve->get_EntryGrade(&start_gradient);
         curve->get_ExitGrade(&end_gradient);

         if (IsEqual(start_gradient, end_gradient))
         {
            // this is just a straight line
            Schema::IfcAlignmentVerticalSegment* linear_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, prev_end_dist_along, l1, prev_end_height, start_gradient, start_gradient, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
            auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, linear_segment);
            file.addEntity(segment);
            profile_segments->push(segment);

            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
            auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(start_gradient), sin(start_gradient)});
            auto parent_curve = new Schema::IfcLine(vertical_point, new Schema::IfcVector(gradient_direction, 1.0));
            auto placement = new Schema::IfcAxis2Placement3D(vertical_point, new Schema::IfcDirection(std::vector<double>{0, 0, 1}), gradient_direction);
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(prev_end_dist_along), new Schema::IfcNonNegativeLengthMeasure(l1), parent_curve);
            curve_segments->push(curve_segment);
         }
         else
         {
            Float64 k;
            curve->get_K1(&k);
            Schema::IfcAlignmentVerticalSegment* parabolic_segment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along, horizontal_length, start_height, start_gradient, end_gradient, IsZero(k) ? Float64_Max : 1 / k, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
            auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, parabolic_segment);
            file.addEntity(segment);
            profile_segments->push(segment);

            Float64 A = (end_gradient - start_gradient) / (2 * horizontal_length);
            Float64 B = start_gradient;
            Float64 C = start_height;
            auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, start_height});
            auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(start_gradient), sin(start_gradient)});
            auto placement = new Schema::IfcAxis2Placement3D(vertical_point, new Schema::IfcDirection(std::vector<double>{0, 0, 1}), gradient_direction);
            auto parent_curve = new Schema::IfcPolynomialCurve(placement, boost::none, std::vector<double>{A, B, C}, boost::none);
            auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_CONTSAMEGRADIENT, placement, new Schema::IfcNonNegativeLengthMeasure(start_dist_along), new Schema::IfcNonNegativeLengthMeasure(horizontal_length), parent_curve);
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
#pragma Reminder("WORKING HERE - How to deal with station equations?") // see IfcReferent
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
      auto segment = new Schema::IfcAlignmentSegment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, nullptr, nullptr, linear_segment);
      file.addEntity(segment);
      profile_segments->push(segment);

      auto vertical_point = new Schema::IfcCartesianPoint(std::vector<double>{0.0, prev_end_height});
      auto gradient_direction = new Schema::IfcDirection(std::vector<double>{cos(prev_end_gradient), sin(prev_end_gradient)});
      auto parent_curve = new Schema::IfcLine(vertical_point, new Schema::IfcVector(gradient_direction, 1.0));
      auto placement = new Schema::IfcAxis2Placement3D(vertical_point, new Schema::IfcDirection(std::vector<double>{0, 0, 1}), gradient_direction);
      auto curve_segment = new Schema::IfcCurveSegment(Schema::IfcTransitionCode::IfcTransitionCode_DISCONTINUOUS, placement, new Schema::IfcNonNegativeLengthMeasure(prev_end_dist_along), new Schema::IfcNonNegativeLengthMeasure(length), parent_curve);
      curve_segments->push(curve_segment);

      // check elevation
      ATLASSERT(IsEqual(endElevation, prev_end_height + length * prev_end_gradient));
   }

   auto vertical_profile = new Schema::IfcAlignmentVertical(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Vert Alignment"), boost::none, boost::none, file.getSingle<Schema::IfcLocalPlacement>(), nullptr);
   file.addEntity(vertical_profile);

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Vert Segments"), std::string("Nests vertical alignment segments with vertical alignment"), vertical_profile, profile_segments);
   file.addEntity(nests);

   *pvertical_profile = vertical_profile;


   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentationItem>> representation_items(new IfcTemplatedEntityList<Schema::IfcRepresentationItem>());
   
   // define the roadway surface geometric representation with IfcSectionedSurface
#pragma Reminder("WORKING HERE - This geometric construction does not take into account the different ways PGSuper defines slope of the roadway section")
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProfileDef>> cross_sections(new IfcTemplatedEntityList<Schema::IfcProfileDef>());
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcPointByDistanceExpression>> cross_section_positions(new IfcTemplatedEntityList<Schema::IfcPointByDistanceExpression>());

   GET_IFACE2(pBroker, IRoadwayData, pRoadway);
   const RoadwaySectionData& roadway_sections = pRoadway->GetRoadwaySectionData();
   Float64 ref_station = pRoadway->GetAlignmentData2().RefStation;
   for (const RoadwaySectionTemplate& section_template : roadway_sections.RoadwaySectionTemplates)
   {
      auto point_on_alignment = new Schema::IfcPointByDistanceExpression(new Schema::IfcNonNegativeLengthMeasure(section_template.Station - ref_station), 0.0, 0.0, 0.0, horizontal_geometry_base_curve);
      cross_section_positions->push(point_on_alignment);

      std::vector<double> widths;
      std::vector<double> slopes;
      if (roadway_sections.NumberOfSegmentsPerSection == 2)
      {
         widths.push_back(100);
         slopes.push_back(section_template.LeftSlope);
         widths.push_back(100);
         slopes.push_back(section_template.RightSlope);
      }
      else
      {
         for (const RoadwaySegmentData& segment_data : section_template.SegmentDataVec)
         {
            widths.push_back(segment_data.Length);
            slopes.push_back(segment_data.Slope);
         }
      }

      auto cross_section = new Schema::IfcOpenCrossProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_CURVE, boost::none, true, widths, slopes, boost::none);
      cross_sections->push(cross_section);
   }
   auto sectioned_surface = new Schema::IfcSectionedSurface(horizontal_geometry_base_curve, cross_section_positions, cross_sections, false);
   representation_items->push(sectioned_surface);

   if (bSimplifiedAlignment)
   {
      // Instead of IfcGradientCurve, we are using a generalized 3D polyline geometric representation of the alignment (a 3D wire)
      // This isn't as accurate, but some viewer may be able to deal with this better
      IndexType nAlignmentPoints = 100;
      Float64 stationInc = (endStation - startStation) / (nAlignmentPoints + 1);
      boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
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
      representation_items->push(polyline);
   }
   else
   {
      auto gradient_curve = new Schema::IfcGradientCurve(curve_segments, false, horizontal_geometry_base_curve, nullptr);
      representation_items->push(gradient_curve);
   }

   auto geometric_representation_context = file.getRepresentationContext(std::string("3D"));
   ATLASSERT(geometric_representation_context);
   auto representation_subcontext = new Schema::IfcGeometricRepresentationSubContext(std::string("Axis"), std::string("Model"), geometric_representation_context, boost::none, Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_GRAPH_VIEW, boost::none);
   auto shape_representation = new Schema::IfcShapeRepresentation(representation_subcontext, std::string("Axis"), std::string("Curve3D"), representation_items);
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentation>> representations(new IfcTemplatedEntityList<Schema::IfcRepresentation>());
   representations->push(shape_representation);

   auto alignment_product_definition_shape = new Schema::IfcProductDefinitionShape(std::string("Alignment_Profile_3D"), boost::none, representations);

   *palignment_product_definition_shape = alignment_product_definition_shape;
}

template <typename Schema>
typename Schema::IfcAlignment* CreateAlignment_4x3(IfcHierarchyHelper<Schema>& file, IBroker* pBroker,bool bSimplifiedAlignment)
{
   USES_CONVERSION;
   Schema::IfcAlignmentHorizontal* horizontal_alignment;
   Schema::IfcCompositeCurve* horizontal_geometry_base_curve;
   CreateHorizontalAlignment_4x3<Schema>(file, pBroker, &horizontal_alignment, &horizontal_geometry_base_curve);

   Schema::IfcAlignmentVertical* vertical_profile;
   Schema::IfcProductDefinitionShape* alignment_product_definition_shape;
   CreateVerticalProfile_4x3<Schema>(file, pBroker, horizontal_geometry_base_curve, bSimplifiedAlignment, &vertical_profile, &alignment_product_definition_shape);

   Schema::IfcObjectPlacement* local_placement;
   auto site = file.getSingle<Schema::IfcSite>();
   if (site->hasObjectPlacement())
   {
      local_placement = site->ObjectPlacement();
   }
   else
   {
      local_placement = file.addLocalPlacement();
   }

   GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
   std::string strAlignmentName(T2A(pRoadwayData->GetAlignmentData2().Name.c_str()));
   if (strAlignmentName.length() == 0) strAlignmentName = "Unnamed";
   auto alignment = new Schema::IfcAlignment(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), strAlignmentName, boost::none, boost::none, local_placement, alignment_product_definition_shape, boost::none);
   file.addEntity(alignment);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> related_objects(new IfcTemplatedEntityList<Schema::IfcObjectDefinition>());
   related_objects->push(horizontal_alignment);
   related_objects->push(vertical_profile);

   auto nests = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Nest Horz and Vert Alignment"), boost::none, alignment, related_objects);
   file.addEntity(nests);

   // add stationing information
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   auto point_on_alignment = new Schema::IfcPointByDistanceExpression(new Schema::IfcNonNegativeLengthMeasure(0.0), 0.0, 0.0, 0.0, horizontal_geometry_base_curve);
   auto relative_placement = new Schema::IfcAxis2PlacementLinear(point_on_alignment, new Schema::IfcDirection(std::vector<double>{ 0, 0, 1}), new Schema::IfcDirection(std::vector<double>{1,0,0}));
   auto cartesian_position = new Schema::IfcAxis2Placement3D(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}), new Schema::IfcDirection(std::vector<double>{ 0, 0, 1 }), new Schema::IfcDirection(std::vector<double>{1, 0, 0}));
   auto referent_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, cartesian_position);


   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProperty>> pset_station_properties(new IfcTemplatedEntityList<Schema::IfcProperty>());
   pset_station_properties->push(new Schema::IfcPropertySingleValue(std::string("Station"), boost::none, new Schema::IfcLengthMeasure(startStation), nullptr));
   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Pset_Stationing"), boost::none, pset_station_properties);
   file.addEntity(property_set);
   auto stationing_referent = new Schema::IfcReferent(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Start of alignment station"), boost::none, boost::none, referent_placement, nullptr, Schema::IfcReferentTypeEnum::IfcReferentType_STATION, boost::none);
   file.addEntity(stationing_referent);
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> related_stationing_objects(new IfcTemplatedEntityList<Schema::IfcObjectDefinition>());
   related_stationing_objects->push(stationing_referent);
   auto nests_stationing = new Schema::IfcRelNests(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Nests Referents with station information with alignment"), boost::none, alignment, related_stationing_objects);
   file.addEntity(nests_stationing);
   auto rel_defines_by_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("Relates station properties to referent"), boost::none, related_stationing_objects, property_set);
   file.addEntity(rel_defines_by_properties);

   return alignment;
}

template <typename Schema>
typename Schema::IfcProfileDef* CreateSectionProfile(IShapes* pShapes,const pgsPointOfInterest& poi,IntervalIndexType intervalIdx)
{
   CComPtr<IShape> shape;
   IndexType gdrIdx, slabIdx;
   pShapes->GetSegmentShape(intervalIdx, poi, true/*orient*/, pgsTypes::scGirder, &shape, &gdrIdx, &slabIdx);

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

   CComPtr<IPoint2dCollection> polyPoints;
   gdrShape->get_PolyPoints(&polyPoints);

   if (GetVertexOrdering(gdrShape) == CLOCKWISE)
   {
      polyPoints->Reverse();
   }

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
   IndexType nPoints;
   polyPoints->get_Count(&nPoints);
   for (IndexType idx = 0; idx < nPoints; idx++)
   {
      CComPtr<IPoint2d> point;
      polyPoints->get_Item(idx, &point);
      Float64 x, y;
      point->Location(&x, &y);

      // NOTE: Use -x because PGSuper has X > 0 to the right looking down station at the start of the girder and IFC has X < 0
      points->push(new Schema::IfcCartesianPoint(std::vector<double>{-x, y}));
   }

   // make sure the polygon is closed
   CComPtr<IPoint2d> first, last;
   polyPoints->get_Item(0, &first);
   polyPoints->get_Item(nPoints - 1, &last);
   if (first->SameLocation(last) == S_FALSE)
   {
      points->push(*(points->begin()));
   }

   auto polyline = new Schema::IfcPolyline(points);
   auto girder_section = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, std::string("CrossSectionProfile"), polyline);

   return girder_section;
}

template <typename Schema> 
boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcObjectDefinition>> CreateStrands(IfcHierarchyHelper<Schema>& file, IBroker* pBroker,const pgsPointOfInterest& poiStart,const pgsPointOfInterest& poiEnd,typename Schema::IfcObjectPlacement* strand_placement)
{
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> strands(new IfcTemplatedEntityList<Schema::IfcObjectDefinition>());

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

      boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentationItem>> strand_representation_items(new IfcTemplatedEntityList<Schema::IfcRepresentationItem>());
      for (StrandIndexType strandIdx = 0; strandIdx < nStrands; strandIdx++)
      {
         boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());

         CComPtr<IPoint2d> pntStart;
         strand_points_start->get_Item(strandIdx, &pntStart);

         Float64 X, Y;
         pntStart->Location(&X, &Y);
         auto start_point = new Schema::IfcCartesianPoint(std::vector<Float64>{-X, Y, poiStart.GetDistFromStart()});

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
            auto hp = new Schema::IfcCartesianPoint(std::vector<Float64>{-X, Y, poi.GetDistFromStart()});
            points->push(hp);
         }

         CComPtr<IPoint2d> pntEnd;
         strand_points_end->get_Item(strandIdx, &pntEnd);

         pntEnd->Location(&X, &Y);
         auto end_point = new Schema::IfcCartesianPoint(std::vector<Float64>{-X, Y, poiEnd.GetDistFromStart()});
         points->push(end_point);

         auto directrix = new Schema::IfcPolyline(points);
         file.addEntity(directrix);

         auto swept_disk_solid = new Schema::IfcSweptDiskSolid(directrix, pStrand->GetNominalDiameter() / 2, boost::none, boost::none, boost::none);
         file.addEntity(swept_disk_solid);
         strand_representation_items->push(swept_disk_solid);
      }

      auto geometric_representation_context = file.getRepresentationContext(std::string("3D"));
      ATLASSERT(geometric_representation_context);
      auto strand_shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), strand_representation_items);
      boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentation>> strand_shape_representation_list(new IfcTemplatedEntityList<Schema::IfcRepresentation>());
      strand_shape_representation_list->push(strand_shape_representation);
      auto strand_product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, strand_shape_representation_list);

      auto owner_history = file.getSingle<Schema::IfcOwnerHistory>();
      auto strand = new Schema::IfcTendon(IfcParse::IfcGlobalId(), owner_history, strStrandType[strandType], boost::none, boost::none, strand_placement, strand_product_definition_shape, boost::none, boost::none,
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

template <typename Schema>
typename Schema::IfcProduct* CreateGirderSegment_4x3(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcAlignment* alignment)
{
    USES_CONVERSION;

    auto owner_history = file.getSingle<Schema::IfcOwnerHistory>();
    ATLASSERT(owner_history);

    GET_IFACE2(pBroker, IIntervals, pIntervals);
    IntervalIndexType intervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);

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

    GET_IFACE2(pBroker, IGirder, pGirder);
    Float64 Zs = pGirder->GetTopGirderChordElevation(poiStart);

    CComPtr<IDirection> segment_direction;
    pBridge->GetSegmentBearing(segmentKey, &segment_direction);
    Float64 direction;
    segment_direction->get_Value(&direction);

    Float64 segment_grade = pBridge->GetSegmentSlope(segmentKey);

    Float64 Ls = pBridge->GetSegmentLength(segmentKey);

    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> girder_line_points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
    girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}));
    girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{Ls, 0, 0}));
    auto girder_line = new Schema::IfcPolyline(girder_line_points);
    file.addEntity(girder_line);


    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProfileDef>> cross_sections(new IfcTemplatedEntityList<Schema::IfcProfileDef>());
    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>> cross_section_positions(new IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>());

    //CComPtr<IAngle> objStartSkew, objEndSkew;
    //pBridge->GetSegmentSkewAngle(segmentKey, pgsTypes::metStart, &objStartSkew);
    //pBridge->GetSegmentSkewAngle(segmentKey, pgsTypes::metEnd, &objEndSkew);

    //Float64 start_skew, end_skew;
    //objStartSkew->get_Value(&start_skew);
    //objEndSkew->get_Value(&end_skew);

    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentationItem>> representation_items(new IfcTemplatedEntityList<Schema::IfcRepresentationItem>());
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
    //      new Schema::IfcDirection(std::vector<Float64>{1, 0, 0}), // direction the Z-axis of the extrsion in the global X direction
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

    for (const pgsPointOfInterest& poi : vPoi)
    {
       cross_section_positions->push(new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{poi.GetDistFromStart(), 0}), nullptr, nullptr));
       auto girder_section = CreateSectionProfile<Schema>(pShapes, poi, intervalIdx);
       cross_sections->push(girder_section);
    }

    auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions, true/*fixed vertical axis*/);
    representation_items->push(sectioned_solid);

    auto geometric_representation_context = file.getRepresentationContext(std::string("3D"));
    ATLASSERT(geometric_representation_context);
    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentation>> shape_representation_list(new IfcTemplatedEntityList<Schema::IfcRepresentation>());
    auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
    shape_representation_list->push(shape_representation);
    auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);


    std::_tostringstream os;
    os << SEGMENT_LABEL(segmentKey);
    std::string segment_name(T2A(os.str().c_str()));

    auto bridge = file.getSingle<Schema::IfcBridge>();
    ATLASSERT(bridge);
    auto bridge_placement = bridge->ObjectPlacement();
    ATLASSERT(bridge_placement);

    Float64 startSegmentStation, startSegmentOffset;
    pBridge->GetStationAndOffset(poiStart, &startSegmentStation, &startSegmentOffset);

    GET_IFACE2(pBroker, IRoadway, pAlignment);
    Float64 startSegmentElevation = pAlignment->GetElevation(startSegmentStation, 0.0);

    Float64 startStation, startElevation, startGrade;
    CComPtr<IPoint2d> startPoint;
    pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

    auto alignment_representation = alignment->Representation();
    auto representation = *(alignment_representation->Representations()->begin());
    auto representation_item = *(representation->Items()->begin());
    auto directrix = representation_item->as<Schema::IfcCurve>();


    auto point_on_alignment = new Schema::IfcPointByDistanceExpression(new Schema::IfcNonNegativeLengthMeasure(startSegmentStation - startStation), -startSegmentOffset, Zs - startSegmentElevation, 0.0, directrix);
    auto relative_placement = new Schema::IfcAxis2PlacementLinear(point_on_alignment, new Schema::IfcDirection(std::vector<double>{ -segment_grade*cos(direction), -segment_grade * sin(direction), 1}), new Schema::IfcDirection(std::vector<double>{cos(direction), sin(direction), 0}));
    auto cartesian_position = new Schema::IfcAxis2Placement3D(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}), new Schema::IfcDirection(std::vector<double>{ 0, 0, 1 }), new Schema::IfcDirection(std::vector<double>{1, 0, 0}));
    auto segment_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, cartesian_position);

    auto segment = new Schema::IfcBeam(IfcParse::IfcGlobalId(), owner_history, segment_name, boost::none, boost::none, segment_placement, product_definition_shape, boost::none, Schema::IfcBeamTypeEnum::IfcBeamType_GIRDER_SEGMENT);
    file.addEntity(segment);

    //
    // materials
    //

    GET_IFACE2(pBroker, IMaterials, pMaterials);
    GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
    IntervalIndexType releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
    IntervalIndexType liftingIntervalIdx = pIntervals->GetLiftSegmentInterval(segmentKey);
    IntervalIndexType haulingIntervalIdx = pIntervals->GetHaulSegmentInterval(segmentKey);

    // create the material
    auto material = new Schema::IfcMaterial("Precast Segment Concrete", boost::none/*description*/, boost::none/*category*/);
    file.addEntity(material);

    // Pset_MaterialConcrete
    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProperty>> material_concrete_properties(new IfcTemplatedEntityList<Schema::IfcProperty>());
    material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("CompressiveStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc28(segmentKey)), nullptr));
    material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("MaxAggregateSize"), boost::none, new Schema::IfcPositiveLengthMeasure(pMaterials->GetSegmentMaxAggrSize(segmentKey)), nullptr));
    auto pset_material_concrete = new Schema::IfcMaterialProperties(std::string("Pset_MaterialConcrete"), boost::none/*description*/, material_concrete_properties, material);
    file.addEntity(pset_material_concrete);

    // Pset_PrecastConcreteElementGeneral
    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProperty>> precast_concrete_properties(new IfcTemplatedEntityList<Schema::IfcProperty>());
    precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("FormStrippingStrength"),boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey,releaseIntervalIdx)),nullptr));
    precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("LiftingStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, liftingIntervalIdx)), nullptr));
    precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("ReleaseStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx)), nullptr));
    precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("TransportationStrength"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey, haulingIntervalIdx)), nullptr));
    precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("InitialTension"), boost::none, new Schema::IfcPressureMeasure(pStrandGeom->GetJackingStress(segmentKey,pgsTypes::Permanent)), nullptr));
    auto pset_material_precast_concrete = new Schema::IfcMaterialProperties(std::string("Pset_PrecastConcreteElementGeneral"), boost::none, precast_concrete_properties, material);
    file.addEntity(pset_material_precast_concrete);

    //// define the individual properties of the material
    //boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProperty>> individual_properties(new IfcTemplatedEntityList<Schema::IfcProperty>());
    //auto release_strength = new Schema::IfcPropertySingleValue(std::string("f'ci"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc(segmentKey,releaseIntervalIdx)), nullptr);
    //auto concrete_strength = new Schema::IfcPropertySingleValue(std::string("f'c"), boost::none, new Schema::IfcPressureMeasure(pMaterials->GetSegmentFc28(segmentKey)), nullptr);
    //individual_properties->push(release_strength);
    //individual_properties->push(concrete_strength);

    //// create the material properties for the material using the individual properties
    //auto material_properties = new Schema::IfcMaterialProperties(std::string("Material Properties for Precast Segment Concrete"), boost::none/*description*/, pset_material_concrete, material);
    //file.addEntity(material_properties);

    // need a list of entities that are associated with this material
    // right now we are creating a unique material for each segment but we still need the list
    IfcEntityList::ptr segments(new IfcEntityList);
    segments->push(segment);

    // associate the material with the segment (ie segments collection)
    auto rel_associates_materials = new Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), owner_history, std::string("Associates_Concrete_To_Precast_Segment"), boost::none, segments, material);
    file.addEntity(rel_associates_materials);

    //
    // begin modeling of prestressing strands
    //
    
    // place strands relative to the segment
    auto strand_placement = file.addLocalPlacement(segment_placement,
       0,0,0, // (0,0,0) of the strands is at (0,0,0) of the segment
       1,0,0, // direction the Z-axis of the extrsion in the global X direction 
       0,1,0 // direction the X-axis of the cross section in the global Y direction
    );

    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> strands = CreateStrands<Schema>(file, pBroker, poiStart, poiEnd, strand_placement);

    if (0 < strands->size())
    {
       auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), owner_history, std::string("Segment_Aggregates_Strands"), boost::none, segment, strands);
       file.addEntity(rel_aggregates);
    }


    return segment;
}


template <typename Schema>
typename Schema::IfcProduct* CreateClosureJoint_4x3(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CClosureKey& closureKey, typename Schema::IfcAlignment* alignment)
{
    USES_CONVERSION;

    auto owner_history = file.getSingle<Schema::IfcOwnerHistory>();
    ATLASSERT(owner_history);

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
    cogoUtil::Inverse(pntStart,pntEnd,&distance,&closure_direction);

    Float64 direction;
    closure_direction->get_Value(&direction);

    Float64 Lc = pBridge->GetClosureJointLength(closureKey);
    Float64 closure_grade = (Ze - Zs)/Lc;

    GET_IFACE2(pBroker, IShapes, pShapes);
    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> girder_line_points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
    girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}));
    girder_line_points->push(new Schema::IfcCartesianPoint(std::vector<double>{Lc, 0, 0}));
    auto girder_line = new Schema::IfcPolyline(girder_line_points);
    file.addEntity(girder_line);

    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProfileDef>> cross_sections(new IfcTemplatedEntityList<Schema::IfcProfileDef>());
    cross_sections->push(CreateSectionProfile<Schema>(pShapes, poiStart, intervalIdx));
    cross_sections->push(CreateSectionProfile<Schema>(pShapes, poiEnd, intervalIdx));

    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>> cross_section_positions(new IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>());
    auto start_section = new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{0, 0}), nullptr, nullptr);
    auto end_section = new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{Lc, 0}), nullptr, nullptr);
    cross_section_positions->push(start_section);
    cross_section_positions->push(end_section);

    auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions, true/*fixed vertical axis*/);

    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentationItem>> representation_items(new IfcTemplatedEntityList<Schema::IfcRepresentationItem>());
    representation_items->push(sectioned_solid);

    auto geometric_representation_context = file.getRepresentationContext(std::string("3D"));
    ATLASSERT(geometric_representation_context);
    boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentation>> shape_representation_list(new IfcTemplatedEntityList<Schema::IfcRepresentation>());
    auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
    shape_representation_list->push(shape_representation);
    auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);


    std::_tostringstream os;
    os << CLOSURE_LABEL(closureKey);
    std::string segment_name(T2A(os.str().c_str()));

    auto bridge = file.getSingle<Schema::IfcBridge>();
    ATLASSERT(bridge);
    auto bridge_placement = bridge->ObjectPlacement();
    ATLASSERT(bridge_placement);

    Float64 startClosureStation, startClosureOffset;
    pBridge->GetStationAndOffset(poiStart, &startClosureStation, &startClosureOffset);

    GET_IFACE2(pBroker, IRoadway, pAlignment);
    Float64 startClosureElevation = pAlignment->GetElevation(startClosureStation, 0.0);

    Float64 startStation, startElevation, startGrade;
    CComPtr<IPoint2d> startPoint;
    pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

    auto alignment_representation = alignment->Representation();
    auto representation = *(alignment_representation->Representations()->begin());
    auto representation_item = *(representation->Items()->begin());
    auto directrix = representation_item->as<Schema::IfcCurve>();


    auto point_on_alignment = new Schema::IfcPointByDistanceExpression(new Schema::IfcNonNegativeLengthMeasure(startClosureStation - startStation), -startClosureOffset, Zs - startClosureElevation, 0.0, directrix);
    auto relative_placement = new Schema::IfcAxis2PlacementLinear(point_on_alignment, new Schema::IfcDirection(std::vector<double>{ -closure_grade * cos(direction), -closure_grade * sin(direction), 1}), new Schema::IfcDirection(std::vector<double>{cos(direction), sin(direction), 0}));
    auto cartesian_position = new Schema::IfcAxis2Placement3D(new Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}), new Schema::IfcDirection(std::vector<double>{ 0, 0, 1 }), new Schema::IfcDirection(std::vector<double>{1, 0, 0}));
    auto closure_placement = new Schema::IfcLinearPlacement(nullptr, relative_placement, cartesian_position);

    auto segment = new Schema::IfcBeam(IfcParse::IfcGlobalId(), owner_history, segment_name, boost::none, std::string("Cast in Place Closure Joint"), closure_placement, product_definition_shape, boost::none, Schema::IfcBeamTypeEnum::IfcBeamType_USERDEFINED);
    file.addEntity(segment);
    return segment;
}

template <typename Schema>
typename Schema::IfcProduct* CreateDeck_4x3(IfcHierarchyHelper<Schema>& file, IBroker* pBroker,  typename Schema::IfcAlignment* alignment)
{
   GET_IFACE2(pBroker, IBridge, pBridge);

   //Float64 startStation = pBridge->GetPierStation(0);
   //Float64 endStation = pBridge->GetPierStation(pBridge->GetPierCount() - 1);
   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);

   GET_IFACE2(pBroker, IShapes, pShapes);
   CComPtr<IShape> slab_shape;
   pShapes->GetSlabShape(startBrgStation, nullptr, true, &slab_shape);

   // The slab shape is in Bridge Section Coordinates. For IFC, we need it relative to the alignment.
   // That is, the top of the slab that is concident with the alignment needs to be at elevation 0.0.
   // Get the alignment elevation at the slab cut station and lower the slab shape by this amount.
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 Z = pAlignment->GetElevation(startBrgStation, 0.0);

   CComPtr<IPoint2dCollection> polyPoints;
   slab_shape->get_PolyPoints(&polyPoints);

   if (GetVertexOrdering(slab_shape) == CLOCKWISE)
   {
      polyPoints->Reverse();
   }

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
   IndexType nPoints;
   polyPoints->get_Count(&nPoints);

   for (IndexType idx = 0; idx < nPoints; idx++)
   {
      CComPtr<IPoint2d> point;
      polyPoints->get_Item(idx, &point);
      Float64 x, y;
      point->Location(&x, &y);

      // NOTE: Use -x because PGSuper has X > 0 to the right looking down station at the start of the girder and IFC has X < 0
      points->push(new Schema::IfcCartesianPoint(std::vector<double>{-x, y - Z}));
   }

   // make sure the polygon is closed
   CComPtr<IPoint2d> first, last;
   polyPoints->get_Item(0, &first);
   polyPoints->get_Item(nPoints - 1, &last);
   if(first->SameLocation(last) == S_FALSE)
   {
      points->push(*(points->begin()));
   }

   auto polyline = new Schema::IfcPolyline(points);
   auto slab_section = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, boost::none, polyline);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProfileDef>> cross_sections(new IfcTemplatedEntityList<Schema::IfcProfileDef>());
   cross_sections->push(slab_section);
   cross_sections->push(slab_section);

   auto alignment_representation = alignment->Representation();
   auto representation = *(alignment_representation->Representations()->begin());
   auto representation_item = *(representation->Items()->begin());
   auto directrix = representation_item->as<Schema::IfcCurve>();

   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>> cross_section_positions(new IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>());
   auto start_section = new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{startBrgStation - startStation, 0}), nullptr, nullptr);
   auto end_section = new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{endBrgStation - startStation, 0}), nullptr, nullptr);
   cross_section_positions->push(start_section);
   cross_section_positions->push(end_section);

   auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(directrix, cross_sections, cross_section_positions, true/*Y is always up*/);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentationItem>> representation_items(new IfcTemplatedEntityList<Schema::IfcRepresentationItem>());
   representation_items->push(sectioned_solid);

   auto geometric_representation_context = file.getRepresentationContext(std::string("3D"));
   ATLASSERT(geometric_representation_context);
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentation>> shape_representation_list(new IfcTemplatedEntityList<Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

   auto deck = new Schema::IfcSlab(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), std::string("MyDeck"), boost::none, boost::none, nullptr/*placement*/, product_definition_shape, boost::none, Schema::IfcSlabTypeEnum::IfcSlabType_NOTDEFINED);
   file.addEntity(deck);
   return deck;
}

template <typename Schema>
typename Schema::IfcProduct* CreateRailingSystem_4x3(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, pgsTypes::TrafficBarrierOrientation tbOrientation, typename Schema::IfcAlignment* alignment)
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
   CComPtr<IShape> barrier_shape;
   if (tbOrientation == pgsTypes::tboLeft)
      pShapes->GetLeftTrafficBarrierShape(startBrgStation, nullptr, &barrier_shape);
   else
      pShapes->GetRightTrafficBarrierShape(startBrgStation, nullptr, &barrier_shape);

   // See note in CreateDeck<>
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 Z = pAlignment->GetElevation(startBrgStation, 0.0);

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

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentationItem>> representation_items(new IfcTemplatedEntityList<Schema::IfcRepresentationItem>());
   for (IndexType shapeIdx = 0; shapeIdx < nShapes; shapeIdx++)
   {
       CComPtr<ICompositeShapeItem> shape_item;
       composite->get_Item(shapeIdx, &shape_item);

       CComPtr<IShape> shape;
       shape_item->get_Shape(&shape);

       CComPtr<IPoint2dCollection> polyPoints;
       shape->get_PolyPoints(&polyPoints);

       if (GetVertexOrdering(shape) == CLOCKWISE)
       {
           polyPoints->Reverse();
       }

       boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcCartesianPoint>> points(new IfcTemplatedEntityList<Schema::IfcCartesianPoint>());
       IndexType nPoints;
       polyPoints->get_Count(&nPoints);

       // make sure the polygon is closed
       CComPtr<IPoint2d> first, last;
       polyPoints->get_Item(0, &first);
       polyPoints->get_Item(nPoints - 1, &last);
       if (first->SameLocation(last) == S_FALSE)
       {
           polyPoints->Add(first);
           nPoints++;
       }

       for (IndexType idx = 0; idx < nPoints; idx++)
       {
           CComPtr<IPoint2d> point;
           polyPoints->get_Item(idx, &point);
           Float64 x, y;
           point->Location(&x, &y);

           // NOTE: Use -x because PGSuper has X > 0 to the right looking down station at the start of the girder and IFC has X < 0
           points->push(new Schema::IfcCartesianPoint(std::vector<double>{-x, y - Z}));
       }

       auto polyline = new Schema::IfcPolyline(points);
       auto barrier_section = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, boost::none, polyline);

       boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProfileDef>> cross_sections(new IfcTemplatedEntityList<Schema::IfcProfileDef>());
       cross_sections->push(barrier_section);
       cross_sections->push(barrier_section);

       auto alignment_representation = alignment->Representation();
       auto representation = *(alignment_representation->Representations()->begin());
       auto representation_item = *(representation->Items()->begin());
       auto directrix = representation_item->as<Schema::IfcCurve>();

       Float64 startStation, startElevation, startGrade;
       CComPtr<IPoint2d> startPoint;
       pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

       boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>> cross_section_positions(new IfcTemplatedEntityList<Schema::IfcAxis2PlacementLinear>());
       auto start_section = new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{startBrgStation - startStation, 0}), nullptr, nullptr);
       auto end_section = new Schema::IfcAxis2PlacementLinear(new Schema::IfcCartesianPoint(std::vector<double>{endBrgStation - startStation, 0}), nullptr, nullptr);
       cross_section_positions->push(start_section);
       cross_section_positions->push(end_section);

       auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(directrix, cross_sections, cross_section_positions, true/*Y is always up*/);

       representation_items->push(sectioned_solid);
   }

   auto geometric_representation_context = file.getRepresentationContext(std::string("3D"));
   ATLASSERT(geometric_representation_context);
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcRepresentation>> shape_representation_list(new IfcTemplatedEntityList<Schema::IfcRepresentation>());
    auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
    shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

   std::string barrier_name(tbOrientation == pgsTypes::tboLeft ? "Left Barrier" : "Right Barrier");

   auto railing = new Schema::IfcRailing(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), barrier_name, boost::none, boost::none, nullptr/*placement*/, product_definition_shape, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_GUARDRAIL);
   file.addEntity(railing);
   return railing;
}

template <typename Schema>
typename Schema::IfcBridge* CreateBridge_4x3(IfcHierarchyHelper<Schema>& file,IBroker* pBroker,  typename Schema::IfcAlignment* alignment)
{
   USES_CONVERSION;

   auto site = file.getSingle<Schema::IfcSite>();
   ATLASSERT(site);

   auto owner_history = file.getSingle<Schema::IfcOwnerHistory>();
   ATLASSERT(owner_history);

   // place the bridge relative to the alignment
   ATLASSERT(alignment->hasRepresentation());
   auto representation = alignment->Representation();
   auto representations = representation->Representations();
   auto items = (*(representations->begin()))->Items();
   auto representation_item = *(items->begin());
   auto directrix = representation_item->as<Schema::IfcCurve>();

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 startStation = pBridge->GetPierStation(0);
   auto placement_point = new Schema::IfcPointByDistanceExpression(new Schema::IfcNonNegativeLengthMeasure(startStation), boost::none, boost::none, boost::none, directrix);
   auto relative_placement = new Schema::IfcAxis2PlacementLinear(placement_point, nullptr/*Z axis is up*/, nullptr /*RefDirection - makes X be normal to curve*/);
   auto bridge_placement = new Schema::IfcLinearPlacement(nullptr/*placement is relative to the origin of the alignment directrix*/, relative_placement, nullptr);

   auto bridge = new Schema::IfcBridge(IfcParse::IfcGlobalId(), owner_history, std::string("My Bridge"), std::string("Sample Bridge"), boost::none, bridge_placement, nullptr, boost::none, Schema::IfcElementCompositionEnum::IfcElementComposition_COMPLEX, Schema::IfcBridgeTypeEnum::IfcBridgeType_GIRDER);
   file.addEntity(bridge);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProduct>> related_elements(new IfcTemplatedEntityList<Schema::IfcProduct>());

   // Bug in IfcOpenShell is preventing the use of IfcFacilityPart to create a sub-spatial structure for the superstructure
   auto bridge_part_type = new Schema::IfcBridgePartTypeEnum(Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);
   //IfcWrite::IfcWriteArgument* attr = new IfcWrite::IfcWriteArgument();
   //attr->set(IfcWrite::IfcWriteArgument::EnumerationReference(Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE, Schema::IfcBridgePartTypeEnum::ToString(Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE)));
   //auto data_ = new IfcEntityInstanceData(&Schema::IfcBridgePartTypeEnum::Class());
   //data_->setArgument(0, attr);
   //auto bridge_part_type = new Schema::IfcBridgePartTypeEnum(data_);

   auto superstructure = new Schema::IfcFacilityPart(IfcParse::IfcGlobalId(), owner_history, std::string("Superstructure"), boost::none, boost::none, bridge_placement, nullptr, boost::none, 
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL, 
      bridge_part_type, 
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL);
   file.addEntity(superstructure);
   related_elements->push(superstructure);

   GroupIndexType nGroups = pBridge->GetGirderGroupCount();
   for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
   {
      GirderIndexType nGirders = pBridge->GetGirderCount(grpIdx);
      for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
      {
         std::_tostringstream os;
         os << GIRDER_LABEL(CGirderKey(grpIdx, gdrIdx));
         std::string girder_name(T2A(os.str().c_str()));

         auto girder = new Schema::IfcElementAssembly(IfcParse::IfcGlobalId(), owner_history, girder_name, boost::none, boost::none, bridge_placement, nullptr, boost::none, Schema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_FACTORY, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
         file.addEntity(girder);

         related_elements->push(girder);

         boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> segments(new IfcTemplatedEntityList<Schema::IfcObjectDefinition>());

         SegmentIndexType nSegments = pBridge->GetSegmentCount(grpIdx, gdrIdx);
         for (SegmentIndexType segIdx = 0; segIdx < nSegments; segIdx++)
         {
            auto segment = CreateGirderSegment_4x3<Schema>(file, pBroker, CSegmentKey(grpIdx, gdrIdx, segIdx), alignment);
            segments->push(segment);

            if (segIdx < nSegments - 1)
            {
                auto closure_joint = CreateClosureJoint_4x3<Schema>(file, pBroker, CClosureKey(grpIdx, gdrIdx, segIdx), alignment);
                segments->push(closure_joint);
            }
         }

         auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), owner_history, std::string("Girder_Aggregates_Segments"), boost::none, girder, segments);
         file.addEntity(rel_aggregates);
      }
   }

   auto deck = CreateDeck_4x3<Schema>(file, pBroker, alignment);
   related_elements->push(deck);

   auto left_barrier = CreateRailingSystem_4x3<Schema>(file, pBroker, pgsTypes::tboLeft, alignment);
   related_elements->push(left_barrier);

   auto right_barrier = CreateRailingSystem_4x3<Schema>(file, pBroker, pgsTypes::tboRight, alignment);
   related_elements->push(right_barrier);

   auto spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), owner_history, boost::none, boost::none, related_elements, bridge);
   file.addEntity(spatial_structure);
   bridge->ContainsElements()->push(spatial_structure);

   return bridge;
}






CIfcModelBuilder::CIfcModelBuilder(void)
{
}

CIfcModelBuilder::~CIfcModelBuilder(void)
{
}

void CIfcModelBuilder::BuildModel(IBroker* pBroker, const CString& strFilePath, CIfcModelBuilder::SchemaType schemaType, bool bSimplifiedAlignment)
{
   switch (schemaType)
   {
   case Schema_4x3_rc3: BuildModel<Ifc4x3_rc3>(pBroker, strFilePath, bSimplifiedAlignment); break;
   case Schema_4x3_rc4: BuildModel<Ifc4x3_rc4>(pBroker, strFilePath, bSimplifiedAlignment); break;
   default:
      ATLASSERT(false); // is there a new schema type
   }
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
   os << "ViewDefinition[MyTestViews]" << std::ends;// << file.schema()->name() << " Issue date: 17 February 2021]" << std::ends;
   file_description.push_back(os.str().c_str());
   file.header().file_description().description(file_description);

   std::_tostringstream _os;
   _os << _T("BridgeLink:") << (pDocType->IsPGSuperDocument() ? _T("PGSuper") : _T("PGSplice")) << _T(" Version ") << pVersionInfo->GetVersion(true).GetBuffer() << std::ends;
   std::string strVersion(T2A(_os.str().c_str()));
   file.header().file_name().originating_system(strVersion);

   //auto project = file.addProject(); // Don't like the default units in IfcOpenShell so we have do build our own
   /////////////////////////// The following is copied from addProject and tweeked
   IfcEntityList::ptr units(new IfcEntityList);
   Schema::IfcDimensionalExponents* dimexp = new Schema::IfcDimensionalExponents(0, 0, 0, 0, 0, 0, 0);
   Schema::IfcSIUnit* unit1 = new Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
   Schema::IfcSIUnit* unit2 = new Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN);

   units->push(unit1);
   units->push(unit2);

   Schema::IfcUnitAssignment* unit_assignment = new Schema::IfcUnitAssignment(units);

   Schema::IfcRepresentationContext::list::ptr rep_contexts(new Schema::IfcRepresentationContext::list);
   Schema::IfcProject* project = new Schema::IfcProject(IfcParse::IfcGlobalId(), owner_history, std::string("MyProject"), boost::none, boost::none, boost::none, boost::none, rep_contexts, unit_assignment);

   file.addEntity(dimexp);
   file.addEntity(unit1);
   file.addEntity(unit2);
   file.addEntity(unit_assignment);
   file.addEntity(project);
   ///////////////////////////////////////// end of copy from addProject

   auto site = file.addSite(project);
   auto local_placement = file.getSingle<Schema::IfcLocalPlacement>(); // addSite creates a local placement so get it here

   project->setName(T2A(pProjectProperties->GetBridgeName()));
   site->setName(T2A(pProjectProperties->GetBridgeName()));

#pragma Reminder("WORKING HERE - Need to set up all owner history information")
   owner_history->OwningApplication()->setApplicationFullName(pDocType->IsPGSuperDocument() ? "BridgeLink:PGSuper" : "BridgeLink:PGSplice");
   owner_history->OwningApplication()->setApplicationIdentifier(pDocType->IsPGSuperDocument() ? "PGSuper" : "PGSplice");
   owner_history->OwningApplication()->setVersion(T2A(pVersionInfo->GetVersion(true)));
   owner_history->OwningApplication()->ApplicationDeveloper()->setIdentification("Washington State Department of Transportation, Bridge and Structures Office");
   owner_history->OwningApplication()->ApplicationDeveloper()->setName("Richard Brice, PE");

   auto world_coordinate_system = file.addPlacement3d();
   auto geometric_representation_context = new Schema::IfcGeometricRepresentationContext(std::string("3D"), std::string("Model"), 3, boost::none, world_coordinate_system, nullptr/*true north*/);
   file.addEntity(geometric_representation_context); // ADD THE CONTEXT TO THE FILE!!!!
   auto contexts = project->RepresentationContexts(); // get the old context
   contexts->push(geometric_representation_context); // add the new context
   project->setRepresentationContexts(contexts); // set the context back into the project
}

template <typename Schema>
void CIfcModelBuilder::BuildModel(IBroker* pBroker, const CString& strFilePath, bool bSimplifiedAlignment)
{
   USES_CONVERSION;

   IfcHierarchyHelper<Schema> file;
   InitializeFile<Schema>(file, pBroker, strFilePath);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProduct>> related_elements(new IfcTemplatedEntityList<Schema::IfcProduct>());

   auto alignment = CreateAlignment_4x3<Schema>(file, pBroker, bSimplifiedAlignment);
   related_elements->push(alignment);

   auto bridge = CreateBridge_4x3<Schema>(file, pBroker, alignment);


   auto owner_history = file.getSingle<Schema::IfcOwnerHistory>();
   auto spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), owner_history, boost::none, boost::none, related_elements, bridge);
   file.addEntity(spatial_structure);

   auto site = file.getSingle<Schema::IfcSite>();

   site->ContainsElements()->push(spatial_structure);

   auto road = new Schema::IfcRoad(IfcParse::IfcGlobalId(), owner_history, std::string("Some road"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcElementCompositionEnum::IfcElementComposition_COMPLEX, Schema::IfcRoadTypeEnum::IfcRoadType_NOTDEFINED);
   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcObjectDefinition>> related_objects(new IfcTemplatedEntityList<Schema::IfcObjectDefinition>());
   related_objects->push(road);
   auto agg = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), owner_history, boost::none, boost::none, site, related_objects);
   file.addEntity(road);
   file.addEntity(agg);

   std::ofstream ofs(T2A(strFilePath));
   ofs << file;
}
