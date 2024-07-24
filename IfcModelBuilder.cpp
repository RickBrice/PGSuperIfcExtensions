///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2024  Washington State Department of Transportation
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
#include <IFace\AnalysisResults.h>
#include <EAF\EAFDisplayUnits.h>
#include <EAF\EAFAutoProgress.h>
#include <PgsExt\GirderLabel.h>
#include <PgsExt\PrecastSegmentData.h>

#include <WBFLCogo\CogoHelpers.h>


constexpr IndexType NUM_DECK_SECTIONS = 10;


#pragma Reminder("TODO - generalize the property enum methods and move to IfcHierarchyHelper")
// Need to cache the IfcPropertyEnumeration for lookup - it can be used multiple times by reference
// Need to have a getPropertyEnumeration method
// Need to generalize the enumValues from strings to IfcValue
// createPropertyEnumeratedValue needs two forms, a single value and a vector of values
template <typename Schema>
typename Schema::IfcPropertyEnumeration* createPropertyEnumeration(const std::string& name, std::vector<std::string>& enumValues, typename Schema::IfcUnit* unit = nullptr)
{
   typename aggregate_of<typename Schema::IfcValue>::ptr enum_values(new aggregate_of<typename Schema::IfcValue>());
   for (const auto& value : enumValues)
   {
      enum_values->push(new Schema::IfcLabel(value));
   }

   auto property_enum = new Schema::IfcPropertyEnumeration(name, enum_values, unit);
   return property_enum;
}

template <typename Schema>
typename Schema::IfcPropertyEnumeratedValue* createPropertyEnumeratedValue(const std::string& property_name,typename Schema::IfcPropertyEnumeration* enumeration,const std::string& value)
{
   typename aggregate_of<typename Schema::IfcValue>::ptr list_of_selected_enum_values(new aggregate_of<typename Schema::IfcValue>());
   list_of_selected_enum_values->push(new Schema::IfcLabel(value));
   auto property_enum_value = new Schema::IfcPropertyEnumeratedValue(property_name, boost::none, list_of_selected_enum_values, enumeration);
   return property_enum_value;
}

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
typename Schema::IfcCartesianPoint* ConvertPoint(IPoint2d* pPoint,bool bMirror = false)
{
   Float64 x, y;
   pPoint->Location(&x, &y);
   x = IsZero(x) ? 0.0 : x;
   y = IsZero(y) ? 0.0 : y;
   return new Schema::IfcCartesianPoint(std::vector<double>{bMirror ? -x : x, y});
}

template <typename Schema>
typename Schema::IfcCurve* CreatePolyline(IShape* shape, const CIfcModelBuilderOptions& options)
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
      points->push(ConvertPoint<Schema>(point,true/*mirror about Y axis*/));
   }

   // we know that points is open (last point is not the same as the first)
   // the polygon must be closed by reference
   points->push(*(points->begin()));

   typename Schema::IfcCurve* curve = nullptr;
   if (options.sweep_profile == CIfcModelBuilderOptions::SweepProfile::IndexPolyCurve)
   {
      std::vector<std::vector<double>> vpoints;
      for (auto point : *points)
      {
         auto coord = point->Coordinates();
         vpoints.emplace_back(coord);
      }
      auto point_list = new Schema::IfcCartesianPointList2D(vpoints,boost::none);
      curve = new Schema::IfcIndexedPolyCurve(point_list,boost::none,false);
   }
   else
   {
      curve = new Schema::IfcPolyline(points);
   }

   return curve;
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


// creates geometry and business logic segments for horizontal alignment tangent runs
template <typename Schema>
std::pair<typename Schema::IfcCurveSegment*, typename Schema::IfcAlignmentSegment*> create_tangent(typename Schema::IfcCartesianPoint* p, double dir, double length,const CIfcModelBuilderOptions& options)
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


template <typename Schema>
void CreateHorizontalAlignment(IfcHierarchyHelper<Schema>& file,IBroker* pBroker, const CIfcModelBuilderOptions& options, typename Schema::IfcAlignmentHorizontal** phorizontal_alignment, typename Schema::IfcRelNests** pnests_horizontal_segments,typename Schema::IfcCompositeCurve** phorizontal_geometry_base_curve)
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

   if(prevPoint->SameLocation(endPoint) == S_FALSE)
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
void CreateVerticalProfile(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcCompositeCurve* horizontal_geometry_base_curve, const CIfcModelBuilderOptions& options, typename Schema::IfcAlignmentVertical** pvertical_profile, typename Schema::IfcRelNests** pnests_vertical_segments, typename Schema::IfcGradientCurve** palignment_gradient_curve)
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
            auto [geometry_segment, business_segment] = create_vcurve<Schema>(vertical_point, start_gradient, end_gradient, horizontal_length,options);
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
void CreateAlignment(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
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


   // add stationing information
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   typename Schema::IfcCurve* curve = nullptr;
   if (gradient_curve) curve = gradient_curve;
   else if (composite_curve) curve = composite_curve;
   else curve = polyline;
   auto point_on_alignment = new Schema::IfcPointByDistanceExpression(
      new Schema::IfcLengthMeasure(0.0), 
      boost::none, boost::none, boost::none, 
      curve);
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

   auto polyline = CreatePolyline<Schema>(gdrShape, options);
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

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 slope = pBridge->GetSegmentSlope(segmentKey);

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

         Float64 X, Y, Z; // X = distance along beam, Z = vertical distance in beam section, Y = horizontal distance in beam section = Z.cross(X)
         pntStart->Location(&Y, &Z);
         X = poiStart.GetDistFromStart() * sqrt(1 + slope * slope); // adjust distance along plan length to distance along girder

         auto start_point = new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z});

         points->push(start_point);

         auto begin = std::begin(strands_at_harp_points);
         auto end = std::end(strands_at_harp_points);
         for (auto iter = begin; iter != end; iter++)
         {
            auto harp_points(*iter);
            CComPtr<IPoint2d> point;
            harp_points->get_Item(strandIdx, &point);

            auto i = std::distance(begin, iter);
            const pgsPointOfInterest& poiHP = vHP[i];

            point->Location(&Y, &Z);
            X = poiHP.GetDistFromStart() * sqrt(1 + slope * slope);
            auto hp = new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z});
            points->push(hp);
         }

         CComPtr<IPoint2d> pntEnd;
         strand_points_end->get_Item(strandIdx, &pntEnd);

         pntEnd->Location(&Y, &Z);
         X = poiEnd.GetDistFromStart() * sqrt(1 + slope * slope);
         auto end_point = new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z});
         points->push(end_point);

         auto directrix = new Schema::IfcPolyline(points);
         file.addEntity(directrix);

         // NOTE: IfcSweptDiskSolid is not part of AbV.
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

      // 6.3.4.9 Pset_ElementComponentCommon
      typename aggregate_of<typename Schema::IfcProperty>::ptr element_component_common_properties(new aggregate_of<typename Schema::IfcProperty>());

      if (strandType == pgsTypes::Temporary)
      {
         // 6.1.8.8 PEnum_ElementStatus
         // This is the only PSet I could find with TEMPORARY so use it for temporary strands
         std::vector<std::string> enum_values{ "DEMOLISH","EXISTING","NEW","TEMPORARY","OTHER","NOTKNOWN","UNSET" };
         auto element_status_enum = createPropertyEnumeration<Schema>("PEnum_ElementStatus", enum_values);
         auto enum_value = createPropertyEnumeratedValue<Schema>("Status", element_status_enum, "TEMPORARY");
         element_component_common_properties->push(enum_value);
      }

      // 6.3.8.1 PEnum_ElementComponentCorrosionTreatment
      std::vector<std::string> enum_values{ "EPOXYCOATED","GALVANISED","NONE","PAINTED","STAINLESS","NOTDEFINED" };
      auto corrosion_treatment_enum_values = createPropertyEnumeration<Schema>("PEnum_ElementComponentCorrosionTreatment", enum_values);
      auto corrosion_treatment_type = createPropertyEnumeratedValue<Schema>("CorrosionTreatment", corrosion_treatment_enum_values, pStrand->GetCoating() == WBFL::Materials::PsStrand::Coating::None ? "NONE" : "EPOXYCOATED");
      element_component_common_properties->push(corrosion_treatment_type);

      auto pset_element_component_common = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ElementComponentCommon"), boost::none, element_component_common_properties);
      file.addEntity(pset_element_component_common);


      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_strands(new aggregate_of<typename Schema::IfcObjectDefinition>());
      related_strands->push(strand);

      auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_strands, pset_element_component_common);
      file.addEntity(related_properties);
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
   os << "ViewDefinition[Alignment-basedView]" << std::ends;
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
   // owner_history->OwningApplication()->ApplicationDeveloper() is IfcOrganization
   auto organization = owner_history->OwningApplication()->ApplicationDeveloper();
   organization->setIdentification(std::string("Washington State Department of Transportation, Bridge and Structures Office"));
   organization->setName(std::string("Richard Brice, PE"));

   // this is an optional parameter, but the AASHTO IDS requires it
   typename aggregate_of<typename Schema::IfcActorRole>::ptr roles(new aggregate_of<typename Schema::IfcActorRole>());
   auto role = new Schema::IfcActorRole(Schema::IfcRoleEnum::IfcRole_CIVILENGINEER, boost::none, boost::none);
   file.addEntity(role);
   roles->push(role);
   organization->setRoles(roles);

   owner_history->OwningApplication()->setApplicationDeveloper(organization);
}


template <typename Schema>
void CreateGirderSegmentRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, const CIfcModelBuilderOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
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

   Float64 Ls = pBridge->GetSegmentLength(segmentKey);
   Float64 Lg = pBridge->GetSegmentPlanLength(segmentKey);
   Float64 slope = pBridge->GetSegmentSlope(segmentKey);

   if (variationType == pgsTypes::svtParabolic)
   {
      // single parabola
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

   CComPtr<IPoint2d> pntStart, pntEnd;
   pBridge->GetPoint(poiStart, pgsTypes::pcGlobal, &pntStart);
   pBridge->GetPoint(poiEnd, pgsTypes::pcGlobal, &pntEnd);

   GET_IFACE2(pBroker, IGirder, pGirder);

   Float64 sx, sy;
   pntStart->Location(&sx, &sy);
   Float64 sz = pGirder->GetTopGirderChordElevation(poiStart);


   Float64 ex, ey;
   pntEnd->Location(&ex, &ey);
   Float64 ez = pGirder->GetTopGirderChordElevation(poiEnd);

   typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr girder_line_points(new aggregate_of<typename Schema::IfcCartesianPoint>());
   // build the girder model in a simple coordinate system, then use ObjectPlacement to local in space
   girder_line_points->push(new Schema::IfcCartesianPoint({ 0,0,0 }));
   girder_line_points->push(new Schema::IfcCartesianPoint({ Lg,0,0 })); // due East from origin, length is plan length, not basic segment length
   auto girder_line = new Schema::IfcPolyline(girder_line_points);
   file.addEntity(girder_line);


   typename aggregate_of<typename Schema::IfcProfileDef>::ptr cross_sections(new aggregate_of<typename Schema::IfcProfileDef>());

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   GET_IFACE2(pBroker, IShapes, pShapes);

   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());

   for (const pgsPointOfInterest& poi : vPoi)
   {
      auto girder_perimeter = CreateSectionProfile<Schema>(pShapes, poi, intervalIdx, options);
      file.addEntity(girder_perimeter);
      cross_sections->push(girder_perimeter);

      Float64 x = poi.GetDistFromStart() * sqrt(1 + slope * slope); // adjust distance along plan length to distance along girder
      auto pde = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(x), boost::none, boost::none, boost::none, girder_line);
      file.addEntity(pde);

      auto lp = new Schema::IfcAxis2PlacementLinear(pde, nullptr, nullptr);
      file.addEntity(lp);

      cross_section_positions->push(lp);
   }

   auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions);
   file.addEntity(sectioned_solid);
   representation_items->push(sectioned_solid);

   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

   // Place the segment in 3D space
   WBFL::Geometry::Vector3d ref_direction(ex - sx, ey - sy, ez - sz); // along the length of the girder
   ref_direction.Normalize();
   WBFL::Geometry::Vector3d z(0, 0, 1); // true up direction
   WBFL::Geometry::Vector3d y = z.Cross(ref_direction); // cross product gives Y axis perpendicular to ref_direction and up
   WBFL::Geometry::Vector3d axis = ref_direction.Cross(y); // cross product gives Z axis of the girder
   auto segment_placement = file.addLocalPlacement(nullptr,
      sx, sy, sz,
      axis.X(), axis.Y(), axis.Z(),
      ref_direction.X(), ref_direction.Y(), ref_direction.Z());
   segment->setObjectPlacement(segment_placement);
   segment->setRepresentation(product_definition_shape);
}

template <typename Schema>
void CreateGirderSegmentMaterials(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, const CIfcModelBuilderOptions& options, typename Schema::IfcStyledRepresentation* styled_representation)
{
   USES_CONVERSION;
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   GET_IFACE2(pBroker, IMaterials, pMaterials);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   IntervalIndexType releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
   IntervalIndexType liftingIntervalIdx = pIntervals->GetLiftSegmentInterval(segmentKey);
   IntervalIndexType haulingIntervalIdx = pIntervals->GetHaulSegmentInterval(segmentKey);

   Float64 camber_ratio = 0.0;
   if (options.include_camber)
   {
      // Compute the camber ratio
      // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_PrecastConcreteElementGeneral.htm
      // The camber deflection, measured from the midpoint of a cambered face of a piece to the midpoint of the chord joining the ends of the same face, 
      // as shown in the figure below (figure not provided), divided by the original (nominal) straight length of the face of the piece.
      GET_IFACE2(pBroker, IPointOfInterest, pPoi);
      PoiList vPoi;
      pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
      CHECK(vPoi.size() == 1);
      const pgsPointOfInterest& poiMS = vPoi.front();
      GET_IFACE2(pBroker, ICamber, pCamber);
      Float64 D = pCamber->GetDCamberForGirderSchedule(poiMS, pgsTypes::CreepTime::Max);
      GET_IFACE2(pBroker, IBridge, pBridge);
      Float64 Ls = pBridge->GetSegmentPlanLength(segmentKey);
      camber_ratio = D / Ls;
   }

   // create the material
   auto material = new Schema::IfcMaterial("Precast Segment Concrete", boost::none/*description*/, boost::none/*category*/);
   file.addEntity(material);

   // assigns the presentation styles to the material
   typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
   list_of_representations->push(styled_representation);
   auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, material);
   file.addEntity(material_defintion_representation);

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
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("BatterAtStart"), boost::none, new Schema::IfcPlaneAngleMeasure(0.0), nullptr));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("BatterAtEnd"), boost::none, new Schema::IfcPlaneAngleMeasure(0.0), nullptr));
   if (options.include_camber) {
      precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("CamberAtMidspan"), boost::none, new Schema::IfcRatioMeasure(camber_ratio), nullptr));
   }
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("DesignLocationNumber"), boost::none, new Schema::IfcLabel(T2A(SEGMENT_LABEL(segmentKey))), nullptr));
   auto pset_material_precast_concrete = new Schema::IfcMaterialProperties(std::string("Pset_PrecastConcreteElementGeneral"), boost::none, precast_concrete_properties, material);
   file.addEntity(pset_material_precast_concrete);

   // need a list of entities that are associated with this material
   // right now we are creating a unique material for each segment but we still need the list
   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr segments(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   segments->push(segment);

   // associate the material with the segment (ie segments collection)
   auto rel_associates_materials = new Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Concrete_to_Precast_Segment"), boost::none, segments, material);
   file.addEntity(rel_associates_materials);
}

template <typename Schema>
void CreateStrandRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, typename Schema::IfcStyledRepresentation* styled_representation)
{
   USES_CONVERSION;

   // place strands relative to the segment origin
   auto strand_placement = file.addLocalPlacement(segment->ObjectPlacement());

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   const pgsPointOfInterest& poiStart(vPoi.front());
   const pgsPointOfInterest& poiEnd(vPoi.back());

   auto strands = CreateStrands<Schema>(file, pBroker, poiStart, poiEnd, strand_placement);

   if (0 < strands->size())
   {
      // create the material
      auto strand_material = new Schema::IfcMaterial("Prestressing Strand", boost::none/*description*/, boost::none/*category*/);
      file.addEntity(strand_material);

      // assigns the presentation styles to the material
      typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
      list_of_representations->push(styled_representation);
      auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, strand_material);
      file.addEntity(material_defintion_representation);


      pgsTypes::StrandType strandType = pgsTypes::Straight;
      GET_IFACE2(pBroker, IMaterials, pMaterials);
      const auto* pStrand = pMaterials->GetStrandMaterial(segmentKey, strandType);
      auto fy = pStrand->GetYieldStrength();
      auto fpu = pStrand->GetUltimateStrength();
      auto eu = 0.035; // from ASTM A416 spec

#pragma Reminder("WORKING HERE - Define strand material")
      // Need to clean this up
      // ASTM A416 is for low relaxation strand... PGSuper does low relaxation and stress relieved
      // ASTM A416 is for Grade 250 and Grade 270... PGSuper does grade 300 as well, but there doesn't seem to be an ASTM
      // We are assuming same material for all strands, but that is not the case in the PGSuper data model
      // straight, harped, and temporary can be different - Grade 250, Grade 270, Grade 300
      // Strand size/diameter is a property on IfcTendon
      std::ostringstream os;
      os << "ASTM A416 Grade " << T2A(WBFL::Materials::PsStrand::GetGrade(pStrand->GetGrade(), true/*US units*/).c_str());
      auto grade = os.str();

      // Pset_MaterialSteel
      typename aggregate_of<typename Schema::IfcProperty>::ptr material_steel_properties(new aggregate_of<typename Schema::IfcProperty>());
      //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("YieldStress"), boost::none, new Schema::IfcPressureMeasure(fy), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStress"), boost::none, new Schema::IfcPressureMeasure(fpu), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStrain"), boost::none, new Schema::IfcPositiveRatioMeasure(eu), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("StructuralGrade"), boost::none, new Schema::IfcLabel(grade.c_str()), nullptr));
      auto pset_material_steel = new Schema::IfcMaterialProperties(std::string("Pset_MaterialSteel"), boost::none/*description*/, material_steel_properties, strand_material);
      file.addEntity(pset_material_steel);

      // need a list of entities that are associated with this material
      // right now we are creating a unique material for each strand but we still need the list
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr strands_for_material(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      for (auto& strand : *strands)
      {
         strands_for_material->push(strand);
      }

      // associate the material with the segment (ie segments collection)
      auto rel_associates_materials = new Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Steel_to_Strand"), boost::none, strands_for_material, strand_material);
      file.addEntity(rel_associates_materials);

      auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Segment_Aggregates_Strands"), boost::none, segment, strands);
      file.addEntity(rel_aggregates);
   }
}

template <typename Schema>
void CreateClosureJointRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CClosureKey& closureKey, typename Schema::IfcElementAssembly* closureJoint, const CIfcModelBuilderOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
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
   Float64 sx, sy;
   pntStart->Location(&sx, &sy);
   Float64 sz = pGirder->GetTopGirderChordElevation(poiStart);


   Float64 ex, ey;
   pntEnd->Location(&ex, &ey);
   Float64 ez = pGirder->GetTopGirderChordElevation(poiEnd);

   Float64 Lc = pBridge->GetClosureJointLength(closureKey); // this is a plan length distance, we need the length along the grade
   Float64 slope = (ez - sz) / Lc;

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

   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), representation_type, representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

   // Place the segment in 3D space
   WBFL::Geometry::Vector3d ref_direction(ex - sx, ey - sy, ez - sz); // along the length of the girder
   ref_direction.Normalize();
   WBFL::Geometry::Vector3d z(0, 0, 1); // true up direction
   WBFL::Geometry::Vector3d y = z.Cross(ref_direction); // cross product gives Y axis perpendicular to ref_direction and up
   WBFL::Geometry::Vector3d axis = ref_direction.Cross(y); // cross product gives Z axis of the girder
   auto closure_placement = file.addLocalPlacement(nullptr,
      sx, sy, sz,
      axis.X(), axis.Y(), axis.Z(),
      ref_direction.X(), ref_direction.Y(), ref_direction.Z());
   closureJoint->setObjectPlacement(closure_placement);
   closureJoint->setRepresentation(product_definition_shape);
}

template <typename Schema>
void CreateDeckRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcBridgePart* deck, const CIfcModelBuilderOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
#pragma Reminder("WORKING HERE - Deck Model - need to re-think this approach")
   // Consider modeling the slab separately from the haunch. The basic slab is the same everywhere.
   // Each girder has it's own haunch.
   // This should eliminate the problem with different number of points in the cross section profile.
   // The deck representation would be a composite of the main slab and each haunch

   // This is not a good model of the deck. This model just creates NUM_DECK_SECTIONS cross sections and extrudes between them.
   GET_IFACE2(pBroker, IBridge, pBridge);
   USES_CONVERSION;

   if (pBridge->GetDeckType() == pgsTypes::sdtNone)
      return;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);

   CComPtr<IDirection> objDir;
   pBridge->GetPierDirection(0, &objDir);
   Float64 startDir;
   objDir->get_Value(&startDir);
   objDir.Release();
   pBridge->GetPierDirection(pBridge->GetPierCount() - 1, &objDir);
   Float64 endDir;
   objDir->get_Value(&endDir);

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file,options);

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   IndexType nDeckSections = NUM_DECK_SECTIONS;

   objDir.Release();
   objDir.CoCreateInstance(CLSID_Direction);

   IndexType point_count = 0;

   GET_IFACE2(pBroker, IShapes, pShapes);
   typename aggregate_of<typename Schema::IfcProfileDef>::ptr cross_sections(new aggregate_of<typename Schema::IfcProfileDef>());
   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());
   for (IndexType i = 0; i <= nDeckSections; i++)
   {
      auto station = i * (endBrgStation - startBrgStation) / nDeckSections + startBrgStation;

      // This code, and the objDir in GetSlabShape, are trying to account for skew by sweeping the cut line angle
      // between the start and end of the bridge. However, there appears to be an issue in WBFL::CoordinateGeometry 
      // that causes the top of deck section from the roadway to have duplicate points. Duplicate points are
      // not valid for the IFC section shape so we will just use a normal section cut for now
      //auto dir = i * (endDir - startDir) / nDeckSections + startDir;
      //objDir->put_Value(dir);

      CComPtr<IShape> slab_shape;
      pShapes->GetSlabShape(station, nullptr/*objDir*/, true/*include haunch*/, &slab_shape);

      // All of the deck cross sections must have exactly the same number of points or it is an invalid IFC representation
      // Capture the number of points for the first deck section, then compare all other deck sections
      // If the point count is different, just skip it.
      // It is typically different between spans at continuous piers... this is a hack, see note above about
      // modeling skews
      CComPtr<IPoint2dCollection> polyPoints;
      slab_shape->get_PolyPoints(&polyPoints);
      IndexType nPoints;
      polyPoints->get_Count(&nPoints);
      if (i == 0)
      {
         point_count = nPoints;
      }
      else
      {
         if (point_count != nPoints)
            continue;
      }

      double elev = pAlignment->GetElevation(station, 0.0);
      CComQIPtr<IXYPosition> pos(slab_shape);
      pos->Offset(0.0, -elev);

      auto polyline = CreatePolyline<Schema>(slab_shape, options);
      std::ostringstream os;
      os << "Deck Section at Station " << T2A(WBFL::COGO::Station(station).AsString(station_format).c_str());
      auto deck_perimeter = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, os.str(), polyline);
      cross_sections->push(deck_perimeter);
      file.addEntity(deck_perimeter);

      auto pde = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(station - startStation), boost::none, boost::none, boost::none, directrix);
      auto deck_section_placement = new Schema::IfcAxis2PlacementLinear(pde, nullptr/*axis*/, nullptr/*ref_direction*/);
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

   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   file.addEntity(product_definition_shape);

   auto slab = new Schema::IfcSlab(IfcParse::IfcGlobalId(), nullptr, std::string("Deck Slab"), boost::none, boost::none, deck_placement, product_definition_shape, boost::none,
      Schema::IfcSlabTypeEnum::IfcSlabType_FLOOR); // see Ifc 4x3 6.1.2.19.2 (FLOOR represents a bridge deck), name is option but AASHTO IDS requires it
   typename aggregate_of<typename Schema::IfcProduct>::ptr list_of_slabs(new aggregate_of<typename Schema::IfcProduct>());
   list_of_slabs->push(slab);

   auto rel_slab_contained_in_deck = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, std::string("Places slab into spatial structure of deck"), boost::none, list_of_slabs, deck);
   file.addEntity(rel_slab_contained_in_deck);
}


template <typename Schema>
void CreateRailingSystemRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, pgsTypes::TrafficBarrierOrientation tbOrientation, typename Schema::IfcProduct* railing, const CIfcModelBuilderOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
   GET_IFACE2(pBroker, IBarriers, pBarriers);

   bool bHasSidewalk = pBarriers->HasSidewalk(tbOrientation);
   bool bHasInteriorBarrier = pBarriers->HasInteriorBarrier(tbOrientation);

   IndexType nShapesPerBarrier = 1 + (bHasSidewalk ? 1 : 0) + (bHasInteriorBarrier ? 1 : 0);

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);

   IndexType nSections = NUM_DECK_SECTIONS;
   std::vector<std::pair<Float64,CComPtr<IShape>>> barrier_shapes;

   GET_IFACE2(pBroker, IShapes, pShapes);
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   if (tbOrientation == pgsTypes::tboLeft)
   {
      for (IndexType i = 0; i <= nSections; i++)
      {
         auto station = i * (endBrgStation - startBrgStation) / nSections + startBrgStation;
         CComPtr<IShape> shape;
         pShapes->GetLeftTrafficBarrierShape(station, nullptr, &shape);

         double elev = pAlignment->GetElevation(station, 0.0);
         CComQIPtr<IXYPosition> pos(shape);
         pos->Offset(0.0, -elev);

         barrier_shapes.emplace_back(station,shape);
      }
   }
   else
   {
      for (IndexType i = 0; i <= nSections; i++)
      {
         auto station = i * (endBrgStation - startBrgStation) / nSections + startBrgStation;
         CComPtr<IShape> shape;
         pShapes->GetRightTrafficBarrierShape(station, nullptr, &shape);

         double elev = pAlignment->GetElevation(station, 0.0);
         CComQIPtr<IXYPosition> pos(shape);
         pos->Offset(0.0, -elev);

         barrier_shapes.emplace_back(station, shape);
      }
   }

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file, options);

   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   typename aggregate_of<typename Schema::IfcAxis2PlacementLinear>::ptr cross_section_positions(new aggregate_of<typename Schema::IfcAxis2PlacementLinear>());
   std::vector<typename aggregate_of<typename Schema::IfcProfileDef>::ptr> cross_sections;
   for (int i = 0; i < nShapesPerBarrier; i++)
   {
      typename aggregate_of<typename Schema::IfcProfileDef>::ptr ptr(new aggregate_of<typename Schema::IfcProfileDef>());
      cross_sections.push_back(ptr);
   }

   for(auto [station,barrier_shape] : barrier_shapes)
   {
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
      ATLASSERT(nShapesPerBarrier == _nShapes); // if this fires the actual number of shapes is not the same as the expected number of shapes
#endif

      auto distance_along = station - startStation;
      auto pde = new Schema::IfcPointByDistanceExpression(new Schema::IfcLengthMeasure(distance_along), boost::none, boost::none, boost::none, directrix);
      auto placement = new Schema::IfcAxis2PlacementLinear(pde, nullptr, nullptr);
      cross_section_positions->push(placement);
      file.addEntity(pde);
      file.addEntity(placement);

      for (IndexType shapeIdx = 0; shapeIdx < nShapesPerBarrier; shapeIdx++)
      {
         CComPtr<ICompositeShapeItem> shape_item;
         composite->get_Item(shapeIdx, &shape_item);

         CComPtr<IShape> shape;
         shape_item->get_Shape(&shape);

         auto polyline = CreatePolyline<Schema>(shape, options);
         std::ostringstream os;
         os << (tbOrientation == pgsTypes::tboLeft ? "Left" : "Right") << " Barrier";
         if (1 < nShapesPerBarrier)
            os << " Shape " << shapeIdx;

         auto shape_perimeter = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, os.str(), polyline);
         cross_sections[shapeIdx]->push(shape_perimeter);
      } // next shape
   } // next section

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   for (IndexType shapeIdx = 0; shapeIdx < nShapesPerBarrier; shapeIdx++)
   {
      auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(directrix, cross_sections[shapeIdx], cross_section_positions);
      representation_items->push(sectioned_solid);
      file.addEntity(sectioned_solid);
   }

   auto site = file.getSingle<typename Schema::IfcSite>();
   auto railing_placement = site->ObjectPlacement();


   typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   auto shape_representation = new Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
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

      if (pBridge->IsAbutment(pierIdx))
         Classify_TPFAbutment<Schema>(file,pier);
      else
         Classify_TPFPier<Schema>(file,pier);

      std::ostringstream os;
      os << "Foundation at " << pier_name; // name is not required, but is specified in AASHTO IDS
      auto foundation = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, nullptr, nullptr, boost::none,
         Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
         Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
         Schema::IfcBridgePartTypeEnum::IfcBridgePartType_FOUNDATION);
      file.addEntity(foundation);
      Classify_TPFFoundation<Schema>(file, foundation);

      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_foundations(new aggregate_of<typename Schema::IfcObjectDefinition>());
      list_of_foundations->push(foundation);

      // IfcBridgePart::PIER <-> IfcRelAggregates <-> IfcBridgePart::FOUNDATION
      auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Foundation is an aggregate component of pier"), boost::none, pier, list_of_foundations);
      file.addEntity(rel_aggregates);


      list_of_piers->push(pier);
   }

   return list_of_piers;
}

template <typename Schema>
typename Schema::IfcStyledRepresentation* CreateMaterialRepresentation(std::string name,double r,double g,double b,typename Schema::IfcGeometricRepresentationContext* geometric_representation_context)
{
   auto color = new Schema::IfcColourRgb(name,r,g,b);
   auto ssr = new Schema::IfcSurfaceStyleRendering(color, boost::none, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, Schema::IfcReflectanceMethodEnum::IfcReflectanceMethod_NOTDEFINED);
   typename aggregate_of<typename Schema::IfcSurfaceStyleElementSelect>::ptr list_of_surface_styles(new aggregate_of<typename Schema::IfcSurfaceStyleElementSelect>());
   list_of_surface_styles->push(ssr);
   auto ss = new Schema::IfcSurfaceStyle(name, Schema::IfcSurfaceSide::IfcSurfaceSide_BOTH, list_of_surface_styles);
   typename aggregate_of<typename Schema::IfcPresentationStyle>::ptr list_of_presentation_styles(new aggregate_of<typename Schema::IfcPresentationStyle>());
   list_of_presentation_styles->push(ss);
   auto styled_item = new Schema::IfcStyledItem(nullptr, list_of_presentation_styles, boost::none);
   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr styled_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   styled_items->push(styled_item);
   auto styled_representation = new Schema::IfcStyledRepresentation(geometric_representation_context, boost::none, boost::none, styled_items);
   return styled_representation;
}

template <typename Schema>
void CreateBridge(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   // we are using the TPFBridge bSDD for classifications
   auto classification = new Schema::IfcClassification(std::string("TPF Bridge")/*Source*/, std::string("1") /*Edition*/, std::string("2024-06-08") /*EditionDate*/, std::string("TPFBridge (USA)"), boost::none /*Description*/, std::string("https://search.bsdd.buildingsmart.org/uri/aashto") /*Specification*/, boost::none /*ReferenceTokens*/);
   file.addEntity(classification);


   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);

   auto body_model_representation_subcontext = new Schema::IfcGeometricRepresentationSubContext(std::string("Body"), std::string("Model"), geometric_representation_context, boost::none, Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_MODEL_VIEW, boost::none);
   file.addEntity(body_model_representation_subcontext);


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
   Classify_TPFSuperstructure<Schema>(file, superstructure);

   auto substructure = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Substructure"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUBSTRUCTURE);
   file.addEntity(substructure);
   Classify_TPFSubstructure<Schema>(file, substructure);

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
   CreateDeckRepresentation(file, pBroker, deck, options, body_model_representation_subcontext);
   file.addEntity(deck);
   list_of_superstructure_spatial_elements->push(deck);

   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelAggregates <-> IfcBridgePart::DECK
   auto rel_aggregates_elements_of_superstructure_spatial_structure = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in superstructure spatial structure"), boost::none, superstructure, list_of_superstructure_spatial_elements);
   file.addEntity(rel_aggregates_elements_of_superstructure_spatial_structure);



   // Add railings to the spatial structure of the superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcRailing
   typename aggregate_of<typename Schema::IfcProduct>::ptr list_of_superstructure_elements(new aggregate_of<typename Schema::IfcProduct>());
   auto left_railing = new Schema::IfcWall(IfcParse::IfcGlobalId(), nullptr, std::string("Left Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcWallTypeEnum::IfcWallType_PARAPET);
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboLeft, left_railing, options, body_model_representation_subcontext);
   file.addEntity(left_railing);
   list_of_superstructure_elements->push(left_railing);

   auto right_railing = new Schema::IfcWall(IfcParse::IfcGlobalId(), nullptr, std::string("Right Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcWallTypeEnum::IfcWallType_PARAPET);
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboRight, right_railing, options, body_model_representation_subcontext);
   file.addEntity(right_railing);
   list_of_superstructure_elements->push(right_railing);

   // Add girders to the spatial structure of the superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcElementAssembly::GIRDER

   // first create a representation object for the girder material. this is one way to add presentation information such as color.
   // this helper function just sets color and hard codes all the other parameters - this can be expanded in the future
   auto girder_material_representation = CreateMaterialRepresentation<Schema>("Girder", 7.6078431372549E-1, 7.72549019607843E-1, 8.E-1, geometric_representation_context);
   auto strand_material_representation = CreateMaterialRepresentation<Schema>("Strand", 1, 0, 0, geometric_representation_context);

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

         auto girder = new Schema::IfcElementAssembly(IfcParse::IfcGlobalId(), nullptr, girder_name, boost::none, boost::none, 
            nullptr/*ObjectPlacement - to be set in CreateGirderSegmentRepresentation*/, 
            nullptr/*Representation - to be set in CreateGirderSegmentRepresentation*/,
            boost::none, boost::none, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
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
            CreateGirderSegmentRepresentation<Schema>(file, pBroker, segmentKey, segment, options, body_model_representation_subcontext);
            CreateGirderSegmentMaterials<Schema>(file, pBroker, segmentKey, segment, options, girder_material_representation);
            
            CreateStrandRepresentation<Schema>(file, pBroker, segmentKey, segment, strand_material_representation);

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
               auto closure_joint = new Schema::IfcElementAssembly(IfcParse::IfcGlobalId(), nullptr, closure_joint_name, std::string("Cast in place concrete closure joint"), std::string("CLOSUREJOINT"), nullptr, nullptr, boost::none, Schema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_SITE, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_USERDEFINED);
               CreateClosureJointRepresentation<Schema>(file, pBroker, segmentKey, closure_joint, options, body_model_representation_subcontext);
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
         } // next segment

         std::ostringstream os_relationship_name;
         os_relationship_name << "Elements of girder for Group " << LABEL_GROUP(grpIdx) << " Girder " << T2A(LABEL_GIRDER(gdrIdx));
         auto girder_aggregation_name = os_relationship_name.str();
         auto girder_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, girder_aggregation_name, boost::none, girder, list_of_girder_segments);
         file.addEntity(girder_aggregates);
      } // next girder
   } // next group

   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcRailing, IfcElementAssembly::GIRDER
   auto rel_contained_in_superstructure_spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in superstructure spatial structure"), boost::none, list_of_superstructure_elements, superstructure);
   file.addEntity(rel_contained_in_superstructure_spatial_structure);

   Create_Pset_BridgeCommon<Schema>(file,bridge);
   Create_TPFBridge_BridgeCommon<Schema>(file, pBroker, bridge);
   Classify_TPFBridge<Schema>(file, bridge);
 }

template <typename Schema>
bool CIfcModelBuilder::BuildModel(IBroker* pBroker, const CString& strFilePath, const CIfcModelBuilderOptions& options)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IProgress, pProgress);
   CEAFAutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Exporting IFC model"));

   IfcHierarchyHelper<Schema> file;
   InitializeFile<Schema>(file, pBroker, strFilePath); // creates project and site
   Create_Pset_ProjectCommon<Schema>(file);
   Create_TPFBridge_ProjectCommon<Schema>(file);

   CreateAlignment<Schema>(file, pBroker, options); // creates alignment and aggregates with project, references into site spatial structure

   if (options.model_elements == CIfcModelBuilderOptions::ModelElements::AlignmentAndBridge)
   {
      CreateBridge<Schema>(file, pBroker, options); // creates bridge with site spatial structure
   }


   std::ofstream ofs(T2A(strFilePath));
   ofs << file;

   return true;
}

template <typename Schema>
void Create_Pset_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   auto project = file.getSingle<typename Schema::IfcProject>();

   // 5.1.8.1 PEnum_ProjectType
   std::vector<std::string> enum_values{ "MODIFICAITON","NEWBUILD","OPERATIONMAINTENANCE","RENOVATION","REPAIR" };
   auto project_type_enum = createPropertyEnumeration<Schema>("PEnum_ProjectType", enum_values);
   auto project_type_property = createPropertyEnumeratedValue<Schema>("ProjectType", project_type_enum, "NEWBUILD");

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(project_type_property);

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ProjectCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_projects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_projects->push(project);

   auto project_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_projects, property_set);
   file.addEntity(project_properties);
}

template <typename Schema>
void Create_TPFBridge_ProjectCommon(IfcHierarchyHelper<Schema>& file)
{
   auto project = file.getSingle<typename Schema::IfcProject>();

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("ContractNumber"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("DesignNumber"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("ProjectNumber"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("ProjectWebsite"), boost::none, new Schema::IfcLabel(std::string("Unknown")), nullptr));

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_ProjectCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_projects(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_projects->push(project);

   auto project_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_projects, property_set);
   file.addEntity(project_properties);
}

template <typename Schema>
void Create_Pset_BridgeCommon(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
{
   // 5.4.8.4 PEnum_StructureIndicator
   std::vector<std::string> enum_values{ "COATED","COMPOSITE","HOMOGENEOUS" };
   auto penum = createPropertyEnumeration<Schema>("PEnum_StructureIndicator", enum_values);
   auto property = createPropertyEnumeratedValue<Schema>("StructureIndicator", penum, "COMPOSITE");

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(property);

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_BridgeCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_bridges(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_bridges->push(bridge);

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Create_TPFBridge_BridgeCommon(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, typename Schema::IfcBridge* bridge)
{
   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nSpans = pBridge->GetSpanCount();
   auto nPiers = pBridge->GetPierCount();

   typename aggregate_of<typename Schema::IfcProperty>::ptr list_of_properties(new aggregate_of<typename Schema::IfcProperty>());
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("tpfBridge_NumberOfSpans"), boost::none, new Schema::IfcInteger((int)nSpans), nullptr));
   list_of_properties->push(new Schema::IfcPropertySingleValue(std::string("tpfBridge_NumberOfSupports"), boost::none, new Schema::IfcInteger((int)nPiers), nullptr));

   auto property_set = new Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("TPFBridge_BridgeCommon"), boost::none, list_of_properties);

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_bridges(new aggregate_of<typename Schema::IfcObjectDefinition>());
   related_bridges->push(bridge);

   auto related_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, property_set);
   file.addEntity(related_properties);
}

template <typename Schema>
void Classify_TPFBridge(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridge* bridge)
{
   auto classifications = file.instances_by_type<typename Schema::IfcClassification>();
   auto classification = (*classifications->begin())->as<typename Schema::IfcClassification>();

   auto classification_reference = new Schema::IfcClassificationReference(std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/1/class/tpfBridge_Bridge"),
      boost::none /*Identification*/, std::string("Bridge") /*Name*/, 
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_bridges(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   related_bridges->push(bridge);

   auto related_classes = new Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_bridges, classification_reference);

   file.addEntity(related_classes);
}

template <typename Schema>
void Classify_TPFBridgePart(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* part,const std::string& uri,const std::string& name, const std::string& part_type)
{
   auto classifications = file.instances_by_type<typename Schema::IfcClassification>();
   auto classification = (*classifications->begin())->as<typename Schema::IfcClassification>();

   auto classification_reference = new Schema::IfcClassificationReference(uri,
      boost::none /*Identification*/, name,
      classification,
      boost::none /*Description*/, boost::none /*Sort*/);

   typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_parts(new aggregate_of<typename Schema::IfcDefinitionSelect>());
   related_parts->push(part);

   auto related_classes = new Schema::IfcRelAssociatesClassification(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_parts, classification_reference);

   file.addEntity(related_classes);
}

template <typename Schema>
void Classify_TPFSuperstructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* superstructure)
{
   Classify_TPFBridgePart(file, superstructure, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/1/class/tpfBridge_BridgeSuperstructure"), std::string("Bridge Superstructure"), std::string("IfcBridgePart.SUPERSTRUCTURE"));
}

template <typename Schema>
void Classify_TPFSubstructure(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* substructure)
{
   Classify_TPFBridgePart(file, substructure, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/1/class/tpfBridge_BridgeSubstructure"), std::string("Bridge Substructure"), std::string("IfcBridgePart.SUBSTRUCTURE"));
}

template <typename Schema>
void Classify_TPFAbutment(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* abutment)
{
   Classify_TPFBridgePart(file, abutment, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/1/class/tpfBridge_AbutmentSpatial"), std::string("Abutment (Spatial)"), std::string("IfcBridgePart.ABUTMENT"));
}

template <typename Schema>
void Classify_TPFPier(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* pier)
{
   Classify_TPFBridgePart(file, pier, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/1/class/tpfBridge_PierSpatial"), std::string("Pier (Spatial)"), std::string("IfcBridgePart.PIER"));
}

template <typename Schema>
void Classify_TPFFoundation(IfcHierarchyHelper<Schema>& file, typename Schema::IfcBridgePart* foundation)
{
   Classify_TPFBridgePart(file, foundation, std::string("https://identifier.buildingsmart.org/uri/aashto/tpfBridge/1/class/tpfBridge_Foundation"), std::string("Foundation"), std::string("IfcBridgePart.FOUNDATION"));
}