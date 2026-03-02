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
#include "IfcExporter.h"
#include "IfcAlignmentBuilder.h"
#include "Referents.h"
#include "USBridge_Classifications.h"
#include "Units.h"
#include "Rebar.h"
#include "GirderSheets.h"
#include "Materials.h"
#include "ConstructionSequence.h"

#include <IFace/Tools.h>
#include <IFace\VersionInfo.h>
#include <IFace\DocumentType.h>
#include <IFace\PrestressForce.h>


#include <EAF/AutoProgress.h>
#include <PsgLib\PrecastSegmentData.h>
#include <psgLib/GirderLibraryEntry.h>
#include <WBFLGenericBridgeTools.h>
#include <PsgLib\BridgeDescription2.h>
#include <Plugins\BeamFamilyCLSID.h>
#include <GeomModel/GeomModel.h>

//namespace std
//{
//   double lerp(double a, double b, double t)
//   {
//      return a + t * (b - a);
//   }
//}

constexpr IndexType NUM_DECK_SECTIONS = 10;

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
typename Schema::IfcCurve* CreatePolyline(IPoint2dCollection* polyPoints)
{
   typename Schema::IfcCartesianPoint::list::ptr points(new typename Schema::IfcCartesianPoint::list);
   IndexType nPoints;
   polyPoints->get_Count(&nPoints);

   if (nPoints < 3)
      return nullptr; // there must be at least 3 points in the cross section or this isn't a polygon cross section

   for (IndexType idx = 0; idx < nPoints; idx++)
   {
      CComPtr<IPoint2d> point;
      polyPoints->get_Item(idx, &point);
      points->push(ConvertPoint<Schema>(point, true/*mirror about Y axis*/));
   }

   return new typename Schema::IfcPolyline(points);
}

template <typename Schema>
typename Schema::IfcCurve* CreatePolyline(IShape* shape, const CIfcExportOptions& options,double cut_angle = PI_OVER_2)
{
   CComPtr<IPoint2dCollection> polyPoints;
   shape->get_PolyPoints(&polyPoints);

   if (GetVertexOrdering(shape) == CLOCKWISE)
   {
      polyPoints->Reverse();
   }

   typename Schema::IfcCartesianPoint::list::ptr points(new typename Schema::IfcCartesianPoint::list);
   IndexType nPoints;
   polyPoints->get_Count(&nPoints);

   if (nPoints < 3)
      return nullptr; // there must be at least 3 points in the cross section or this isn't a polygon cross section

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

      double x, y;
      point->Location(&x, &y);
      x /= sin(cut_angle);
      CComPtr<IPoint2d> new_point;
      new_point.CoCreateInstance(CLSID_Point2d);
      new_point->Move(x, y);

      points->push(ConvertPoint<Schema>(new_point,true/*mirror about Y axis*/));
   }

   // we know that points is open (last point is not the same as the first)
   // the polygon must be closed by reference
   points->push(*(points->begin()));

   typename Schema::IfcCurve* curve = nullptr;
   if (options.sweep_profile == CIfcExportOptions::SweepProfile::IndexedPolyCurve)
   {
      std::vector<std::vector<double>> vpoints;
      for (auto point : *points)
      {
         auto coord = point->Coordinates();
         vpoints.emplace_back(coord);
      }
      auto point_list = new typename Schema::IfcCartesianPointList2D(vpoints,boost::none);
      curve = new typename Schema::IfcIndexedPolyCurve(point_list,boost::none,false);
   }
   else
   {
      curve = new typename Schema::IfcPolyline(points);
   }

   return curve;
}


template <typename Schema>
typename Schema::IfcProfileDef* CreateSectionProfile(std::shared_ptr<IShapes> pShapes,const pgsPointOfInterest& poi,IntervalIndexType intervalIdx, const CIfcExportOptions& options,double cut_angle=PI_OVER_2)
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

   auto polyline = CreatePolyline<Schema>(gdrShape, options, cut_angle);
   auto girder_section = new typename Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, std::string("CrossSectionProfile"), polyline);

   return girder_section;
}

template <typename Schema>
typename Schema::IfcTendonType* GetTendonType(IfcHierarchyHelper<Schema>& file, const WBFL::Materials::PsStrand* pStrand)
{
   USES_CONVERSION;
   std::string name(T2A(pStrand->GetName().c_str()));

   // search to see if an IfcTendonType has already been created
   auto project = file.getSingle<typename Schema::IfcProject>();
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   for (auto& rel_declares : *rel_declares_instances)
   {
      if (rel_declares->RelatingContext()->as<typename Schema::IfcProject>())
      {
         auto related_definitions = rel_declares->RelatedDefinitions();
         for (auto& reldef : *related_definitions)
         {
            auto tendon_type = reldef->as<typename Schema::IfcTendonType>();
            if (tendon_type && tendon_type->Name() == name)
            {
               return tendon_type;
            }
         }
      }
   }

   // if we get this far, we need a new IfcTendonType
   auto tendon_type = new typename Schema::IfcTendonType(
      IfcParse::IfcGlobalId(),
      nullptr,
      name, /*Name*/
      boost::none, /*Description*/
      boost::none, /*ApplicableOccurrence*/
      boost::none, /*HasPropertySets*/
      boost::none, /*RepresentationMaps*/
      boost::none, /*Tag*/
      boost::none, /*ElementType*/
      Schema::IfcTendonTypeEnum::IfcTendonType_STRAND, /*PredefinedType*/
      pStrand->GetNominalDiameter(), /*NominalDiameter*/
      pStrand->GetNominalArea(), /*CrossSectionArea*/
      boost::none /*SheathDiameter*/
   );

   Classify_Prestressing<Schema>(file, tendon_type);

   file.addEntity(tendon_type);

   // add the new definition to the project
   if (rel_declares_instances->size() == 0)
   {
      typename Schema::IfcDefinitionSelect::list::ptr related_definitions(new typename Schema::IfcDefinitionSelect::list);
      related_definitions->push(tendon_type);

      auto rel_declares = new typename Schema::IfcRelDeclares(
         IfcParse::IfcGlobalId(),
         nullptr,
         boost::none,
         boost::none,
         project,
         related_definitions);

      file.addEntity(rel_declares);
   }
   else
   {
      for (auto& rel_declares : *rel_declares_instances)
      {
         if (rel_declares->RelatingContext()->as<typename Schema::IfcProject>())
         {
            auto related_definitions = rel_declares->RelatedDefinitions();
            related_definitions->push(tendon_type);
            rel_declares->setRelatedDefinitions(related_definitions);
            break;
         }
      }
   }

   return tendon_type;
}

template <typename Schema> 
typename Schema::IfcObjectDefinition::list::ptr CreateStrands(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker,const pgsPointOfInterest& poiStart,const pgsPointOfInterest& poiEnd,typename Schema::IfcBeam* beam)
{
   typename Schema::IfcObjectDefinition::list::ptr strands(new typename Schema::IfcObjectDefinition::list);

   const CSegmentKey& segmentKey(poiStart.GetSegmentKey());

   GET_IFACE2(pBroker, IBridge, pBridge);
   CComPtr<IAngle> angle_start_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metStart, &angle_start_face);
   Float64 start_face_angle;
   angle_start_face->get_Value(&start_face_angle);


   CComPtr<IAngle> angle_end_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metEnd, &angle_end_face);
   Float64 end_face_angle;
   angle_end_face->get_Value(&end_face_angle);

   Float64 segment_length = pBridge->GetSegmentPlanLength(segmentKey);
   Float64 slope = pBridge->GetSegmentSlope(segmentKey); // need slope to adjust poi (which is a plan view measure) to an along the girder distance

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   GET_IFACE2_NOCHECK(pBroker, IMaterials, pMaterials);

   // place strands relative to the segment origin
   typename Schema::IfcLocalPlacement* strand_placement = nullptr;


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

      strand_placement = (strand_placement == nullptr ? file.addLocalPlacement(beam->ObjectPlacement()) : strand_placement);
      typename Schema::IfcRepresentationItem::list::ptr strand_representation_items(new typename Schema::IfcRepresentationItem::list);
      for (StrandIndexType strandIdx = 0; strandIdx < nStrands; strandIdx++)
      {
         typename Schema::IfcCartesianPoint::list::ptr points(new typename Schema::IfcCartesianPoint::list);

         CComPtr<IPoint2d> pntStart;
         strand_points_start->get_Item(strandIdx, &pntStart);

         Float64 X, Y, Z; // X = distance along beam, Z = vertical distance in beam section, Y = horizontal distance in beam section = Z.cross(X)
         pntStart->Location(&Y, &Z);
         Y *= -1.0;

         X = poiStart.GetDistFromStart() * sqrt(1 + slope * slope); // adjust distance along plan length to distance along girder

         auto start_offset = Y / tan(start_face_angle);
         
         X += start_offset;

         auto start_point = new typename Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z});

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
            Y *= -1.0;

            X = poiHP.GetDistFromStart() * sqrt(1 + slope * slope);
            auto hp = new typename Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z});
            points->push(hp);

            if (iter == begin)
            {
               // if this is the first harp point, the strand elevation at the start face of the beam
               // needs to be adjusted for the start_offset distance.
               auto dx = X - (start_point->Coordinates()[0] - start_offset);
               auto dz = Z - start_point->Coordinates()[2];
               auto harped_strand_slope = dz / dx;
               start_point->Coordinates()[2] += start_offset * harped_strand_slope;
            }
         }

         CComPtr<IPoint2d> pntEnd;
         strand_points_end->get_Item(strandIdx, &pntEnd);

         pntEnd->Location(&Y, &Z);
         Y *= -1.0;

         X = poiEnd.GetDistFromStart() * sqrt(1 + slope * slope);
         auto end_offset = Y / tan(end_face_angle);
         X += end_offset;

         if (strandType == pgsTypes::Harped)
         { 
            auto last_point = *(points->end()-1);
            auto dx = (X - end_offset) - last_point->Coordinates()[0];
            auto dz = Z - last_point->Coordinates()[2];
            auto harped_strand_slope = dz / dx;
            Z += end_offset * harped_strand_slope;
         }

         auto end_point = new typename Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z});
         points->push(end_point);

         auto directrix = new typename Schema::IfcPolyline(points);
         file.addEntity(directrix);

         // NOTE: IfcSweptDiskSolid is not part of AbV.
         auto swept_disk_solid = new typename Schema::IfcSweptDiskSolid(directrix, pStrand->GetNominalDiameter() / 2, boost::none, boost::none, boost::none);
         file.addEntity(swept_disk_solid);
         strand_representation_items->push(swept_disk_solid);

         auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
         ATLASSERT(geometric_representation_context);
         auto strand_shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), strand_representation_items);
         typename Schema::IfcRepresentation::list::ptr strand_shape_representation_list(new typename Schema::IfcRepresentation::list);
         strand_shape_representation_list->push(strand_shape_representation);
         auto strand_product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, strand_shape_representation_list);

         std::ostringstream os;
         os << strStrandType[strandType] << ":" << strandIdx + 1;
         auto strand = new typename Schema::IfcTendon(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, strand_placement, strand_product_definition_shape, boost::none, boost::none,
            boost::none, /*Schema::IfcTendonTypeEnum::IfcTendonType_STRAND,*/ // per 4.1.3.2, this must not be used unless PredefinedType at the ObjectType level is set to NOTDEFINED
            boost::none, /*pStrand->GetNominalDiameter() depreciated*/
            boost::none, /*pStrand->GetNominalArea() depreciated*/
            pStrandGeom->GetPjack(segmentKey, strandType),
            pStrandGeom->GetJackingStress(segmentKey, strandType),
            boost::none, boost::none, boost::none);

         file.addEntity(strand);

         auto* tendon_type = GetTendonType<Schema>(file,pStrand);
         file.addRelatedObject<typename Schema::IfcRelDefinesByType>(tendon_type, strand);


         Float64 Pjack = pStrandGeom->GetPjack(segmentKey, strandType);
         Float64 db_start, db_end;
         bool bDebonded = pStrandGeom->IsStrandDebonded(segmentKey, strandIdx, strandType, nullptr, &db_start, &db_end);
         Create_Pset_usBridge_ReinforcementCommon(file, strand, Pjack, bDebonded, db_start); // assumes symmetric debonding since classification can't handle unsymmetric

         // 6.3.4.9 Pset_ElementComponentCommon
         typename Schema::IfcProperty::list::ptr element_component_common_properties(new typename Schema::IfcProperty::list);

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

         auto pset_element_component_common = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ElementComponentCommon"), boost::none, element_component_common_properties);
         file.addEntity(pset_element_component_common);


         typename Schema::IfcObjectDefinition::list::ptr related_strands(new typename Schema::IfcObjectDefinition::list);
         related_strands->push(strand);

         auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_strands, pset_element_component_common);
         file.addEntity(related_properties);

         strands->push(strand);
      }
   }

   return strands;
}

template <typename Schema>
typename Schema::IfcObjectDefinition::list::ptr CreateRebars(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const pgsPointOfInterest& poiStart, const pgsPointOfInterest& poiEnd, typename Schema::IfcBeam* beam)
{
   USES_CONVERSION;

   const CSegmentKey& segmentKey(poiStart.GetSegmentKey());

   GET_IFACE2(pBroker, IBridge, pBridge);
   CComPtr<IAngle> angle_start_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metStart, &angle_start_face);
   Float64 start_face_angle;
   angle_start_face->get_Value(&start_face_angle);


   CComPtr<IAngle> angle_end_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metEnd, &angle_end_face);
   Float64 end_face_angle;
   angle_end_face->get_Value(&end_face_angle);

   Float64 segment_length = pBridge->GetSegmentPlanLength(segmentKey); // length along grade
   Float64 slope = pBridge->GetSegmentSlope(segmentKey); // need slope to adjust bar start/end distance (which is a plan view measure) to an along the girder distance

   GET_IFACE2(pBroker, ILongitudinalRebar, pLongRebar);
   const CLongitudinalRebarData* pLRD = pLongRebar->GetSegmentLongitudinalRebarData(segmentKey);

   GET_IFACE2(pBroker, ILongRebarGeometry, pLongRebarGeom);
   CComPtr<IRebarLayout> rebar_layout;
   pLongRebarGeom->GetRebarLayout(segmentKey, &rebar_layout);

   IndexType nRebars;
   rebar_layout->get_Count(&nRebars);

   typename Schema::IfcObjectDefinition::list::ptr rebars(new typename Schema::IfcObjectDefinition::list);
   if (nRebars == 0)
      return rebars;


   // place rebar relative to the segment origin
   auto segment_origin = file.addLocalPlacement(beam->ObjectPlacement());

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   CComPtr<IEnumRebarLayoutItems> enum_items;
   rebar_layout->get__EnumRebarLayoutItems(&enum_items);

   IndexType layout_item_idx = 0;
   CComPtr<IRebarLayoutItem> rebar_layout_item;
   while (enum_items->Next(1, &rebar_layout_item, nullptr) != S_FALSE)
   {
      Float64 start, centerline_bar_length;
      rebar_layout_item->get_Start(&start);
      rebar_layout_item->get_Length(&centerline_bar_length); // length of the bar measured along CL Girder

      start *= sqrt(1 + slope * slope);
      centerline_bar_length *= sqrt(1 + slope * slope);

      CComPtr<IEnumRebarPatterns> enum_patterns;
      rebar_layout_item->get__EnumRebarPatterns(&enum_patterns);
      CComPtr<IRebarPattern> rebar_pattern;
      while (enum_patterns->Next(1, &rebar_pattern, nullptr) != S_FALSE)
      {
         typename Schema::IfcReinforcingBarType* rebar_type = nullptr;
         CComPtr<IRebar> rb;
         rebar_pattern->get_Rebar(&rb);

         CComBSTR bar_name;
         rb->get_Name(&bar_name);
         WBFL::Materials::Rebar::Size bar_size = WBFL::LRFD::RebarPool::GetBarSize(OLE2CT(bar_name));
         const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(pLRD->BarType, pLRD->BarGrade, bar_size);

         if (rebar_type == nullptr)
         {
            Float64 db;
            rb->get_NominalDiameter(&db);
            
            // create a basic representation of the bar based on the bar's length at the centerline of the girder
            // This will be used in mapped representations and the bar length will be scaled to the actual bar length
            // accounting for the actual bar's offset from centerline of beam as well as the effect of girder end face skew
            typename Schema::IfcCartesianPoint::list::ptr points(new typename Schema::IfcCartesianPoint::list);
            points->push(new typename Schema::IfcCartesianPoint(std::vector<double>{0., 0., 0.}));
            points->push(new typename Schema::IfcCartesianPoint(std::vector<double>{centerline_bar_length, 0., 0.}));
            auto directrix = new typename Schema::IfcPolyline(points);
            auto swept_disk_solid = new typename Schema::IfcSweptDiskSolid(directrix, db / 2, boost::none, boost::none, boost::none);
            typename Schema::IfcRepresentationItem::list::ptr representation_items(new Schema::IfcRepresentationItem::list);
            representation_items->push(swept_disk_solid);
            typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
            auto shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
            std::ostringstream os;
            os << "Girder_Longitudinal_Bar_" << OLE2A(bar_name);
            rebar_type = CreateReinforcingBarType<Schema>(file, os.str(), false, pRebar, shape_representation);
         }


         typename Schema::IfcRepresentationItem::list::ptr mapped_representation_items(new Schema::IfcRepresentationItem::list);
         IndexType nBars;
         rebar_pattern->get_Count(&nBars);
         for (IndexType barIdx = 0; barIdx < nBars; barIdx++)
         {
            CComPtr<IPoint2d> p1, p2;
            rebar_pattern->get_Location(0.0, barIdx, &p1);

            Float64 X, Y, Z;
            p1->Location(&Y, &Z);
            // p1.X = horizontal from CL beam, p1.Y is distance from top of beam
            // p1.X = -Y in IFC beam coordinates
            // p1.Y = Z in IFC beam coordinates
            Y *= -1.0;

            // X is distance along CL beam in IFC beam coordinates
            X = start;

            Float64 start_offset = 0.0;
            Float64 end_offset = 0.0;

            if (IsEqual(segment_length, centerline_bar_length))
            {
               // Bar runs full length of segment so adjust it's length based on rebar offset from CL of beam and skews
               // No adjustments are made for partial length bars
               // TODO: will need to update this and adjust for partial length bars that are tied to the end faces of the beam

               // adjust start position based on girder start face skew
               start_offset = Y / tan(start_face_angle);
               X += start_offset;

               // length adjustment based on girder end face skew
               end_offset = Y / tan(end_face_angle);
            }

            Float64 actual_bar_length = -start_offset + centerline_bar_length + end_offset;
            Float64 scaleX = actual_bar_length / centerline_bar_length;

            auto rebar_type_representation_maps = rebar_type->RepresentationMaps();
            auto mapping_source = *((*rebar_type_representation_maps)->begin());

            // Use a nonUniform transformation so we can scale only the length of the bar (Uniform transformation scales in all directions which would increase the diameter of the bar - we don't want that)
            auto mapping_target = new typename Schema::IfcCartesianTransformationOperator3DnonUniform(new typename Schema::IfcDirection({ 1.0,0.0,0.0 }), nullptr, new typename Schema::IfcCartesianPoint({ X,Y,Z }), scaleX, nullptr, boost::none, boost::none);
            auto mapped_item = new typename Schema::IfcMappedItem(mapping_source, mapping_target);
            mapped_representation_items->push(mapped_item);
         }

         typename Schema::IfcRepresentation::list::ptr shape_representation_list(new Schema::IfcRepresentation::list);
         auto shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), mapped_representation_items);
         shape_representation_list->push(shape_representation);
         auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

         std::ostringstream os;
         os << "Rebar Row " << (layout_item_idx+1) << " " << T2A(WBFL::LRFD::RebarPool::GetBarSize(bar_size).c_str());

         auto rebar = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, segment_origin, product_definition_shape, boost::none,
            boost::none, // steel grade: depreciated
            boost::none, // nominal diameter: depreciated
            boost::none, // cross section area: depreciated
            boost::none, // bar length: depreciated
            boost::none, // predefined type: depreciated
            boost::none  // predefined type: depreciated
         );
         file.addEntity(rebar);

         file.addRelatedObject<typename Schema::IfcRelDefinesByType>(rebar_type, rebar);

         rebars->push(rebar);
         rebar_pattern.Release();
      }
      rebar_layout_item.Release();
      layout_item_idx++;
   }
   return rebars;
}

template <typename Schema>
typename Schema::IfcObjectDefinition::list::ptr CreateStirrups(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* beam)
{
   // WORKING HERE - The idea is to check to see if the beam is of the IBeam family, otherwise, don't model stirrups (Already doing this step in the calling function)
   // For I-beams, start with WSDOT G2 bars, then change to G1 bars (but there are 2 bars, not 1)... then add the G3 bar in the top flange
   // This is just an experiment for how to model stirrups and a rebar cage.
   // When this is re-built as an extension agent, bar shape will be an input as part of the girder definition

   USES_CONVERSION;

   typename Schema::IfcObjectDefinition::list::ptr rebars(new typename Schema::IfcObjectDefinition::list);
   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   typename Schema::IfcLocalPlacement* segment_origin = nullptr; // only create if needed

   // Assume the beam is constant depth
   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   auto poiStart = pPoi->GetPointOfInterest(segmentKey, 0.0);

   Float64 cover = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);

   GET_IFACE2(pBroker, IMaterials, pMaterials);
   WBFL::Materials::Rebar::Type bar_type;
   WBFL::Materials::Rebar::Grade bar_grade;
   pMaterials->GetSegmentTransverseRebarMaterial(segmentKey, &bar_type, &bar_grade);

   GET_IFACE2(pBroker, IGirder, pGirder);
   Float64 wbf = pGirder->GetBottomFlangeWidth(poiStart, 0);
   Float64 hbf = pGirder->GetBottomFlangeThickness(poiStart, 0) - WBFL::Units::ConvertToSysUnits(4.5,WBFL::Units::Measure::Inch);

   // Create G3 #5 bar type and representation (this is a dummy bar only for WSDOT girders)
   const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(bar_type, bar_grade, WBFL::Materials::Rebar::Size::bs5);
   std::ostringstream os;
   os << "G3 Top Bars";
   auto g3_rebar_type = GetReinforcingBarType<Schema>(file, os.str(), false, pRebar);
   Float64 wtf = pGirder->GetTopFlangeWidth(poiStart);
   Float64 g3_bar_length = wtf - 2*cover;
   if (g3_rebar_type == nullptr)
   {
      Float64 db = pRebar->GetNominalDimension();
      typename Schema::IfcCartesianPoint::list::ptr points(new typename Schema::IfcCartesianPoint::list);
      points->push(new typename Schema::IfcCartesianPoint(std::vector<double>{0., -g3_bar_length / 2, 0.}));
      points->push(new typename Schema::IfcCartesianPoint(std::vector<double>{0., g3_bar_length / 2, 0.}));
      auto directrix = new typename Schema::IfcPolyline(points);
      auto swept_disk_solid = new typename Schema::IfcSweptDiskSolid(directrix, db / 2, boost::none, boost::none, boost::none);
      typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);
      representation_items->push(swept_disk_solid);
      typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
      auto shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);

      g3_rebar_type = CreateReinforcingBarType<Schema>(file, os.str(), false, pRebar, shape_representation);
   }

   // G9 bars
   pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(bar_type, bar_grade, WBFL::Materials::Rebar::Size::bs3);
   os.str("");
   os.clear();
   os << "G9 Bottom Confinement Bars";
   auto g9_rebar_type = GetReinforcingBarType<Schema>(file, os.str(), false, pRebar);

   if (g9_rebar_type == nullptr)
   {
      auto db = pRebar->GetNominalDimension();

      // This is a totally hard coded G9 bar without curves - need to update this later
      std::vector<std::vector<double>> g9_point_list;
      g9_point_list.push_back({ 0.0, (wbf - 2 * cover) / 2, 0.0 });
      g9_point_list.push_back({ 0.0, (wbf - 2 * cover) / 2, hbf - 2 * cover });
      g9_point_list.push_back({ 0.0, 0.0, WBFL::Units::ConvertToSysUnits(9.125, WBFL::Units::Measure::Inch) }); // no way to get height of bottom bulb, this is WSDOT's G9 bar dimension
      g9_point_list.push_back({ 0.0, -(wbf - 2 * cover) / 2, hbf - 2 * cover });
      g9_point_list.push_back({ 0.0, -(wbf - 2 * cover) / 2, 0.0 });

      typename Schema::IfcSegmentIndexSelect::list::ptr g9_segments(new typename Schema::IfcSegmentIndexSelect::list);
      g9_segments->push(new typename Schema::IfcLineIndex({ 1,2 }));
      g9_segments->push(new typename Schema::IfcLineIndex({ 2,3 }));
      g9_segments->push(new typename Schema::IfcLineIndex({ 3,4 }));
      g9_segments->push(new typename Schema::IfcLineIndex({ 4,5 }));

      auto g9_directrix = new typename Schema::IfcIndexedPolyCurve(new typename Schema::IfcCartesianPointList3D(g9_point_list, boost::none), g9_segments, boost::none);

      auto swept_disk_solid = new typename Schema::IfcSweptDiskSolid(g9_directrix, db / 2, boost::none, boost::none, boost::none);
      auto representation_items = Schema::IfcRepresentationItem::list::ptr(new typename Schema::IfcRepresentationItem::list);
      representation_items->push(swept_disk_solid);
      auto shape_representation_list = typename Schema::IfcRepresentation::list::ptr(new typename Schema::IfcRepresentation::list);
      auto shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);

      g9_rebar_type = CreateReinforcingBarType<Schema>(file, os.str(), false, pRebar, shape_representation);
   }

   // G10 bars
   os.str("");
   os.clear();
   os << "G10 Bottom Confinement Bars";
   auto g10_rebar_type = GetReinforcingBarType<Schema>(file, os.str(), false, pRebar);

   if (g10_rebar_type == nullptr)
   {
      auto db = pRebar->GetNominalDimension();

      auto three_inch = WBFL::Units::ConvertToSysUnits(3.0, WBFL::Units::Measure::Inch);
      Float64 r = 4.5 * db;
      Float64 h = three_inch - 5. * db;
      Float64 d = 0.5 * (wbf - 2 * cover - db - 2 * r);
      Float64 delta = PI_OVER_2;

      std::vector<std::vector<double>> g10_point_list;
      g10_point_list.push_back({ 0., (d + r), h + r });
      g10_point_list.push_back({ 0., (d + r), r });
      g10_point_list.push_back({ 0., (d + r * sin(delta / 2)), r * cos(delta / 2) });
      g10_point_list.push_back({ 0., d, 0.0 });
      g10_point_list.push_back({ 0., -d, 0.0 });
      g10_point_list.push_back({ 0., -(d + r * sin(delta / 2)), r * cos(delta / 2) });
      g10_point_list.push_back({ 0., -(d + r), r });
      g10_point_list.push_back({ 0., -(d + r), h + r });

      typename Schema::IfcSegmentIndexSelect::list::ptr g10_segments(new typename Schema::IfcSegmentIndexSelect::list);
      g10_segments->push(new typename Schema::IfcLineIndex({ 1,2 }));
      g10_segments->push(new typename Schema::IfcArcIndex({ 2,3,4 }));
      g10_segments->push(new typename Schema::IfcLineIndex({ 4,5 }));
      g10_segments->push(new typename Schema::IfcArcIndex({ 5,6,7 }));
      g10_segments->push(new typename Schema::IfcLineIndex({ 7,8 }));

      auto g10_directrix = new typename Schema::IfcIndexedPolyCurve(new typename Schema::IfcCartesianPointList3D(g10_point_list, boost::none), g10_segments, boost::none);

      auto swept_disk_solid = new typename Schema::IfcSweptDiskSolid(g10_directrix, db / 2, boost::none, boost::none, boost::none);
      auto representation_items = typename Schema::IfcRepresentationItem::list::ptr(new typename Schema::IfcRepresentationItem::list);
      representation_items->push(swept_disk_solid);
      auto shape_representation_list = typename Schema::IfcRepresentation::list::ptr(new typename Schema::IfcRepresentation::list);
      auto shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);

      g10_rebar_type = CreateReinforcingBarType<Schema>(file, os.str(), false, pRebar, shape_representation);
   }


   // Get some basic geometry for G2 stirrups
   Float64 Hg = pGirder->GetHeight(poiStart);
   Float64 t = pGirder->GetWebThickness(poiStart, 0);

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 Lg = pBridge->GetSegmentPlanLength(segmentKey);
   Float64 A = pBridge->GetSlabOffset(segmentKey, pgsTypes::metStart);
   Float64 H1 = Hg + A + WBFL::Units::ConvertToSysUnits(3.0, WBFL::Units::Measure::Inch); // h1 = Hg + "A" + 3"


   CComPtr<IAngle> angle_start_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metStart, &angle_start_face);
   Float64 start_face_angle;
   angle_start_face->get_Value(&start_face_angle);


   CComPtr<IAngle> angle_end_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metEnd, &angle_end_face);
   Float64 end_face_angle;
   angle_end_face->get_Value(&end_face_angle);

   // Interpolation function of bar angle relative to CL beam.
   // This can be any function, but for now, we hard code it to look like
   // the splay layout on the WSDOT WF Girder 4 of 5 Details sheet
   // except that the length of Zone 1 is taken to be wtf/tan(angle)/2
   // and the length of Zone 2 is taken to be 10 ft. These parameters
   // can be updated in the future based on the stirrup zone layouts
   auto fn_bar_angle = [start_face_angle, 
                        end_face_angle, 
                        Lstart = std::max(wtf,wbf) / fabs(tan(start_face_angle)) / 2,
                        Lsplay = WBFL::Units::ConvertToSysUnits(10.0,WBFL::Units::Measure::Feet), 
                        Lend   = std::max(wtf, wbf) / fabs(tan(end_face_angle)) / 2,
                        Lg](Float64 x)->Float64 {
      if (x < Lstart)
         return start_face_angle;
      else if (x < Lstart + Lsplay)
         return std::lerp(start_face_angle, PI_OVER_2, (x - Lstart) / Lsplay);
      else if (Lg - Lend - Lsplay < x && x < Lg - Lend)
         return std::lerp(PI_OVER_2, end_face_angle, (x - (Lg - Lsplay - Lend)) / Lsplay);
      else if (Lg - Lend < x)
         return end_face_angle;
      else
         return PI_OVER_2;
      };

   GET_IFACE2(pBroker, IStirrupGeometry, pStirrupGeometry);
   ZoneIndexType nZones = pStirrupGeometry->GetPrimaryZoneCount(segmentKey);

   for (ZoneIndexType zoneIdx = 0; zoneIdx < nZones; zoneIdx++)
   {
      Float64 start, end;
      pStirrupGeometry->GetPrimaryZoneBounds(segmentKey, zoneIdx, &start, &end);

      WBFL::Materials::Rebar::Size bar_size;
      Float64 nLegs, spacing;
      pStirrupGeometry->GetPrimaryVertStirrupBarInfo(segmentKey, zoneIdx, &bar_size, &nLegs, &spacing);

      const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(bar_type, bar_grade, bar_size);
      Float64 db = pRebar->GetNominalDimension();

     
      std::ostringstream os;
      os << "G2 bars";
      auto* g2_rebar_type = GetReinforcingBarType<Schema>(file, os.str(), true, pRebar);
      if (g2_rebar_type == nullptr)
      {
         // Bar type doesn't exist, create it
         // Create geometry of a "G2" bar
         Float64 dl = Hg - cover - db / 2; // distance from top of beam to center of hair-pin bend
         Float64 du = H1 - dl; // distance from top of beam upwards to the end of the bar
         Float64 dx = t / 2 - cover - db / 2; // horizontal distance from CL Beam to CL bar (this is basically the bend radius)

         std::vector<std::vector<double>> point_list;
         // X = longitudinal axis of beam (use 0.0 for start face of beam)
         // Y = horizontal distance relative to start face of beam, positive values to the left
         // Z = vertical elevation. From PGSuper, elevation is 0.0 at top of beam
         point_list.push_back({ 0.0,dx,du }); // top left of bar
         point_list.push_back({ 0.0,dx,-(dl - dx) }); // left side of bar at start of bend
         point_list.push_back({ 0.0,0.0,-dl }); // low point at center of hair-pin bend
         point_list.push_back({ 0.0,-dx,-(dl - dx) }); // right side of bar at end of bend
         point_list.push_back({ 0.0,-dx,du }); // top right of bar

         typename Schema::IfcSegmentIndexSelect::list::ptr segments(new typename Schema::IfcSegmentIndexSelect::list);
         segments->push(new typename Schema::IfcLineIndex({ 1,2 }));
         segments->push(new typename Schema::IfcArcIndex({ 2,3,4 }));
         segments->push(new typename Schema::IfcLineIndex({ 4,5 }));

         auto directrix = new typename Schema::IfcIndexedPolyCurve(new typename Schema::IfcCartesianPointList3D(point_list, boost::none), segments, boost::none);
         file.addEntity(directrix);

         auto swept_disk_solid = new typename Schema::IfcSweptDiskSolid(directrix, db / 2, boost::none, boost::none, boost::none);
         file.addEntity(swept_disk_solid);

         typename Schema::IfcRepresentationItem::list::ptr rebar_representation_items(new typename Schema::IfcRepresentationItem::list);
         rebar_representation_items->push(swept_disk_solid);

         auto rebar_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), rebar_representation_items);

         g2_rebar_type = CreateReinforcingBarType<Schema>(file, os.str(), true, pRebar, rebar_representation);
      }

      typename Schema::IfcRepresentationItem::list::ptr g2_mapped_representation_items(new typename Schema::IfcRepresentationItem::list);
      typename Schema::IfcRepresentationItem::list::ptr g3_mapped_representation_items(new typename Schema::IfcRepresentationItem::list);
      typename Schema::IfcRepresentationItem::list::ptr g9_mapped_representation_items(new typename Schema::IfcRepresentationItem::list);
      typename Schema::IfcRepresentationItem::list::ptr g10_mapped_representation_items(new typename Schema::IfcRepresentationItem::list);

      Float64 offset = start + (start < Lg/2.0 ? 1.0 : -1.0)*spacing;
      Float64 sign = (offset < Lg / 2.0 ? 1.0 : -1.0);
      IndexType nBars = (IndexType)((end - start) / spacing);
      for (IndexType barIdx = 0; barIdx < nBars; barIdx++, offset += spacing)
      {
         Float64 bar_angle = fn_bar_angle(offset); // angle of bar in plan view, measured from horizontal
         auto X_direction = new typename Schema::IfcDirection({ sin(bar_angle),-cos(bar_angle),0.0 });
         auto Y_direction = new typename Schema::IfcDirection({ cos(bar_angle),sin(bar_angle),0.0 });
         Float64 scaleY = fabs(1 / sin(bar_angle));

         auto g2_mapping_target = new typename Schema::IfcCartesianTransformationOperator3DnonUniform(X_direction, Y_direction, new typename Schema::IfcCartesianPoint({ offset, 0., -cover }), 1.0, nullptr, scaleY, boost::none);
         auto g2_rebar_type_representation_maps = g2_rebar_type->RepresentationMaps();
         auto g2_mapping_source = *((*g2_rebar_type_representation_maps)->begin());
         auto g2_mapped_item = new typename Schema::IfcMappedItem(g2_mapping_source, g2_mapping_target);
         g2_mapped_representation_items->push(g2_mapped_item);

         // G3 and G2 bars can't have the same offset, otherwise they will conflict with each other
         // Offset the G3 bars 1-db towards the center of the beam (+1 db in left half, and -1 db in right half)
         auto g3_mapping_target = new typename Schema::IfcCartesianTransformationOperator3DnonUniform(X_direction, Y_direction, new typename Schema::IfcCartesianPoint({ offset + sign * db, 0., -cover }), 1.0, nullptr, scaleY, boost::none);
         auto g3_rebar_type_representation_maps = g3_rebar_type->RepresentationMaps();
         auto g3_mapping_source = *((*g3_rebar_type_representation_maps)->begin());
         auto g3_mapped_item = new typename Schema::IfcMappedItem(g3_mapping_source, g3_mapping_target);
         g3_mapped_representation_items->push(g3_mapped_item);

         auto g9_mapping_target = new typename Schema::IfcCartesianTransformationOperator3DnonUniform(X_direction, Y_direction, new typename Schema::IfcCartesianPoint({ offset + sign * db, 0., -(Hg - cover/* - db#3*/) }), 1.0, nullptr, scaleY, boost::none);
         auto g9_rebar_type_representation_maps = g9_rebar_type->RepresentationMaps();
         auto g9_mapping_source = *((*g9_rebar_type_representation_maps)->begin());
         auto g9_mapped_item = new typename Schema::IfcMappedItem(g9_mapping_source, g9_mapping_target);
         g9_mapped_representation_items->push(g9_mapped_item);

         auto g10_mapping_target = new typename Schema::IfcCartesianTransformationOperator3DnonUniform(X_direction, Y_direction, new typename Schema::IfcCartesianPoint({ offset + sign * db, 0., -(Hg - cover/* - db#3*/) }), 1.0, nullptr, scaleY, boost::none);
         auto g10_rebar_type_representation_maps = g10_rebar_type->RepresentationMaps();
         auto g10_mapping_source = *((*g10_rebar_type_representation_maps)->begin());
         auto g10_mapped_item = new typename Schema::IfcMappedItem(g10_mapping_source, g10_mapping_target);
         g10_mapped_representation_items->push(g10_mapped_item);
      }

      typename Schema::IfcRepresentation::list::ptr g2_shape_representation_list(new typename Schema::IfcRepresentation::list);
      auto g2_shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), g2_mapped_representation_items);
      g2_shape_representation_list->push(g2_shape_representation);

      auto g2_product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, g2_shape_representation_list);

      segment_origin = (segment_origin == nullptr ? file.addLocalPlacement(beam->ObjectPlacement()) : segment_origin);

      os.str("");
      os.clear();
      os << "Zone " << LABEL_STIRRUP_ZONE(zoneIdx) << " Stirrups";
      auto g2_rebar = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, segment_origin, g2_product_definition_shape, boost::none,
         boost::none, // steel grade: depreciated
         boost::none, // nominal diameter: depreciated
         boost::none, // cross section area: depreciated
         boost::none, // bar length: depreciated
         boost::none, // predefined type: depreciated
         boost::none  // predefined type: depreciated
      );
      file.addEntity(g2_rebar);

      file.addRelatedObject<typename Schema::IfcRelDefinesByType>(g2_rebar_type, g2_rebar);

      rebars->push(g2_rebar);


      typename Schema::IfcRepresentation::list::ptr g3_shape_representation_list(new Schema::IfcRepresentation::list);
      auto g3_shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), g3_mapped_representation_items);
      g3_shape_representation_list->push(g3_shape_representation);

      auto g3_product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, g3_shape_representation_list);

      segment_origin = (segment_origin == nullptr ? file.addLocalPlacement(beam->ObjectPlacement()) : segment_origin);
      os.str("");
      os.clear();
      os << "Zone " << LABEL_STIRRUP_ZONE(zoneIdx) << " Top Bars";
      auto g3_rebar = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, segment_origin, g3_product_definition_shape, boost::none,
         boost::none, // steel grade: depreciated
         boost::none, // nominal diameter: depreciated
         boost::none, // cross section area: depreciated
         boost::none, // bar length: depreciated
         boost::none, // predefined type: depreciated
         boost::none  // predefined type: depreciated
      );
      file.addEntity(g3_rebar);

      file.addRelatedObject<typename Schema::IfcRelDefinesByType>(g3_rebar_type, g3_rebar);

      rebars->push(g3_rebar);

      typename Schema::IfcRepresentation::list::ptr g9_shape_representation_list(new Schema::IfcRepresentation::list);
      auto g9_shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), g9_mapped_representation_items);
      g9_shape_representation_list->push(g9_shape_representation);

      auto g9_product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, g9_shape_representation_list);

      segment_origin = (segment_origin == nullptr ? file.addLocalPlacement(beam->ObjectPlacement()) : segment_origin);
      os.str("");
      os.clear();
      os << "Zone " << LABEL_STIRRUP_ZONE(zoneIdx) << " G9 Bottom Confinement Bars";
      auto g9_rebar = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, segment_origin, g9_product_definition_shape, boost::none,
         boost::none, // steel grade: depreciated
         boost::none, // nominal diameter: depreciated
         boost::none, // cross section area: depreciated
         boost::none, // bar length: depreciated
         boost::none, // predefined type: depreciated
         boost::none  // predefined type: depreciated
      );
      file.addEntity(g9_rebar);

      file.addRelatedObject<typename Schema::IfcRelDefinesByType>(g9_rebar_type, g9_rebar);

      rebars->push(g9_rebar);


      typename Schema::IfcRepresentation::list::ptr g10_shape_representation_list(new Schema::IfcRepresentation::list);
      auto g10_shape_representation = new typename Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), g10_mapped_representation_items);
      g10_shape_representation_list->push(g10_shape_representation);

      auto g10_product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, g10_shape_representation_list);

      segment_origin = (segment_origin == nullptr ? file.addLocalPlacement(beam->ObjectPlacement()) : segment_origin);
      os.str("");
      os.clear();
      os << "Zone " << LABEL_STIRRUP_ZONE(zoneIdx) << " G10 Bottom Confinement Bars";
      auto g10_rebar = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, segment_origin, g10_product_definition_shape, boost::none,
         boost::none, // steel grade: depreciated
         boost::none, // nominal diameter: depreciated
         boost::none, // cross section area: depreciated
         boost::none, // bar length: depreciated
         boost::none, // predefined type: depreciated
         boost::none  // predefined type: depreciated
      );
      file.addEntity(g10_rebar);

      file.addRelatedObject<typename Schema::IfcRelDefinesByType>(g10_rebar_type, g10_rebar);

      rebars->push(g10_rebar);
   }
   return rebars;
}

template <typename Schema>
void GirderSegment_SectionedSolidHorizontal(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker,const CSegmentKey& segmentKey,const PoiList& vPoi,std::function<Float64(Float64)>fn_cut_angle, typename Schema::IfcBeam* segment, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 Lg = pBridge->GetSegmentPlanLength(segmentKey);
   typename Schema::IfcCartesianPoint::list::ptr girder_line_points(new Schema::IfcCartesianPoint::list);
   girder_line_points->push(new typename Schema::IfcCartesianPoint({ 0,0,0 }));
   girder_line_points->push(new typename Schema::IfcCartesianPoint({ Lg,0,0 })); // due East from origin, length is plan length, not basic segment length
   auto girder_line = new typename Schema::IfcPolyline(girder_line_points);
   file.addEntity(girder_line);


   typename Schema::IfcProfileDef::list::ptr cross_sections(new typename Schema::IfcProfileDef::list);

   typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);
   GET_IFACE2(pBroker, IShapes, pShapes);

   typename Schema::IfcAxis2PlacementLinear::list::ptr cross_section_positions(new typename Schema::IfcAxis2PlacementLinear::list);

   GET_IFACE2(pBroker, IIntervals, pIntervals);
   IntervalIndexType intervalIdx = pIntervals->GetErectSegmentInterval(segmentKey);

   Float64 slope = pBridge->GetSegmentSlope(segmentKey);

   for (const pgsPointOfInterest& poi : vPoi)
   {
      auto cut_angle = fn_cut_angle(poi.GetDistFromStart());
      auto girder_profile = CreateSectionProfile<Schema>(pShapes, poi, intervalIdx, options, cut_angle);
      file.addEntity(girder_profile);
      cross_sections->push(girder_profile);

      Float64 x = poi.GetDistFromStart() * sqrt(1 + slope * slope); // adjust distance along plan length to distance along girder
      auto pde = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(x), boost::none, boost::none, boost::none, girder_line);
      file.addEntity(pde);

      // contrary to the IFC documentation, the RefDirection is normal to the plane of the cross section (at least that is how many of implemented it)

      WBFL::Geometry::Vector3d up(options.batter_ends ? slope : 0, 0, 1); // along the length of the girder
      up.Normalize();

      auto rd = new typename Schema::IfcDirection({ sin(cut_angle),-cos(cut_angle),0.0 }); // normal to the plane of the cross section
      auto axis = new typename Schema::IfcDirection({ up.X(),up.Y(),up.Z() }); // up
      auto lp = new typename Schema::IfcAxis2PlacementLinear(pde, axis, rd);
      file.addEntity(lp);

      cross_section_positions->push(lp);
   }

   auto sectioned_solid = new typename Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions);
   file.addEntity(sectioned_solid);
   representation_items->push(sectioned_solid);

   typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
   auto shape_representation = new typename Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   segment->setRepresentation(product_definition_shape);
}

std::pair<IndexType,std::vector<std::vector<double>>> generate_point_list(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, const PoiList& vPoi, std::function<Float64(Float64)>fn_cut_angle)
{
   std::vector<std::vector<double>> point_list;
   IndexType nPointsPerProfile = 0;

   GET_IFACE2(pBroker, IIntervals, pIntervals);
   auto intervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 slope = pBridge->GetSegmentSlope(segmentKey);

   GET_IFACE2(pBroker, IShapes, pShapes);
   
   for (const pgsPointOfInterest& poi : vPoi)
   {
      auto z = poi.GetDistFromStart() * sqrt(1 + slope * slope);  // adjust distance along plan length to distance along girder

      CComPtr<IShape> shape;
      pShapes->GetSegmentShape(intervalIdx, poi, false, pgsTypes::scGirder, &shape);

      CComPtr<IPoint2dCollection> shape_points;
      shape->get_PolyPoints(&shape_points);

      IndexType nPoints;
      shape_points->get_Count(&nPoints);

      CComPtr<IPoint2d> p0, pn;
      shape_points->get_Item(0, &p0);
      shape_points->get_Item(nPoints - 1, &pn);
      Float64 dist;
      p0->DistanceEx(pn, &dist);
      if (IsZero(dist))
         nPoints--; // profile is closed so skip the last point

      nPointsPerProfile = nPoints;

      for (auto i = 0; i < nPoints; i++)
      {
         CComPtr<IPoint2d> pnt;
         shape_points->get_Item(i, &pnt);

         Float64 x, y;
         pnt->Location(&x, &y);

         // the polygonal face set points are in global X,Y,Z
         // Distance along beam (z) = Global X
         // Vertical Beam Dimension (y) = Global Z
         // Horizontal Beam Dimension (x) = Global Y (negative x because of difference in PGSuper coordinates and IFC global coordinates)
#pragma Reminder("WORKING HERE - polygons - need to skew adjust the local x coordinate")
         point_list.push_back({z, -x, y});
      }
   }

   return { nPointsPerProfile,point_list };
}

std::vector<std::vector<int>> build_faces(IndexType nPointsPerProfile, const std::vector<std::vector<double>>& vPoints)
{
   // outer face vertices must connect counter-clockwise when seen from the outside
   auto nSideFaces = nPointsPerProfile - 1;

   std::vector<std::vector<int>> face_indices;
   auto get_point_index = [nPointsPerProfile](IndexType profileIdx, IndexType pntIdx)->int {return int(profileIdx * nPointsPerProfile + pntIdx) + 1; };
   auto nProfiles = vPoints.size() / nPointsPerProfile;
   for (auto profileIdx = 0; profileIdx < nProfiles - 1; profileIdx++)
   {
      for (auto profilePointIdx = 0; profilePointIdx < nSideFaces; profilePointIdx++)
      {
         auto idx0 = get_point_index(profileIdx, profilePointIdx + 1);
         auto idx1 = get_point_index(profileIdx, profilePointIdx);
         auto idx2 = get_point_index(profileIdx + 1, profilePointIdx);
         auto idx3 = get_point_index(profileIdx + 1, profilePointIdx + 1);
         face_indices.push_back({ idx0,idx1,idx2,idx3 });
      }

      // last face on perimeter of girder
      auto idx0 = get_point_index(profileIdx, 0);
      auto idx1 = get_point_index(profileIdx, nSideFaces);
      auto idx2 = get_point_index(profileIdx + 1, nSideFaces);
      auto idx3 = get_point_index(profileIdx + 1, 0);
      face_indices.push_back({ idx0,idx1,idx2,idx3 });
   }

   // start face
   std::vector<int> start_point_indices;
   for (auto profilePointIdx = 0; profilePointIdx < nPointsPerProfile; profilePointIdx++)
   {
      start_point_indices.push_back(get_point_index(0, profilePointIdx));
   }
   face_indices.push_back(start_point_indices);


   // end face
   std::vector<int> end_point_indices;
   for (auto profilePointIdx = nSideFaces; profilePointIdx != INVALID_INDEX; profilePointIdx--)
   {
      end_point_indices.push_back(get_point_index(nProfiles-1, profilePointIdx));
   }
   face_indices.push_back(end_point_indices);

   return face_indices;
}

template <typename Schema>
void GirderSegment_PolygonalFaceSet(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, const PoiList& vPoi, std::function<Float64(Float64)>fn_cut_angle, typename Schema::IfcBeam* segment, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
   //std::vector<std::vector<double>> point_list;
   //point_list.push_back({ 0.0, 0.0, 0.0 });
   //point_list.push_back({ 1.0, 0.0, 0.0 });
   //point_list.push_back({ 1.0, 1.0, 0.0 });
   //point_list.push_back({ 0.0, 1.0, 0.0 });
   //point_list.push_back({ 0.0, 0.0, 2.0 });
   //point_list.push_back({ 1.0, 0.0, 2.0 });
   //point_list.push_back({ 1.0, 1.0, 2.0 });
   //point_list.push_back({ 0.0, 1.0, 2.0 });
   //auto coordinates = new typename Schema::IfcCartesianPointList3D(point_list, boost::none);
   //auto nPointsPerProfile = 4;

   //std::vector<std::vector<int>> face_indices_list;
   //face_indices_list.push_back({ 1,2,6,5 });
   //face_indices_list.push_back({ 6,2,3,7 });
   //face_indices_list.push_back({ 7,3,4,8 });
   //face_indices_list.push_back({ 8,4,1,5 });
   //face_indices_list.push_back({ 1,4,3,2 });
   //face_indices_list.push_back({ 6,7,8,5 });

   auto [nPointsPerProfile, point_list] = generate_point_list(pBroker,segmentKey,vPoi,fn_cut_angle);
   auto coordinates = new typename Schema::IfcCartesianPointList3D(point_list, boost::none);
   auto face_indices_list = build_faces(nPointsPerProfile, point_list);

   typename Schema::IfcIndexedPolygonalFace::list::ptr faces(new typename Schema::IfcIndexedPolygonalFace::list());
   for (auto face_indices : face_indices_list)
   {
      auto face = new typename Schema::IfcIndexedPolygonalFace(face_indices);
      faces->push(face);
   }

   auto faceset = new typename Schema::IfcPolygonalFaceSet(coordinates,boost::none,faces,boost::none);

   typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);
   representation_items->push(faceset);

   typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
   auto shape_representation = new typename Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("Tessellation"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   segment->setRepresentation(product_definition_shape);
}

template <typename Schema>
void GirderSegment_FacetedBrep(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, const PoiList& vPoi, std::function<Float64(Float64)>fn_cut_angle, typename Schema::IfcBeam* segment, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
   auto [nPointsPerProfile, point_list] = generate_point_list(pBroker, segmentKey, vPoi, fn_cut_angle);
   auto face_indices_list = build_faces(nPointsPerProfile, point_list);

   typename Schema::IfcFace::list::ptr faces(new typename Schema::IfcFace::list);
   for (auto& face_indices : face_indices_list)
   {
      // create IfcPolyLoop with each point being point_list[face_indices[i]-1];
      typename Schema::IfcCartesianPoint::list::ptr polygon(new typename Schema::IfcCartesianPoint::list);
      for (auto idx : face_indices)
      {
         auto& p = point_list[idx - 1];
         polygon->push(new typename Schema::IfcCartesianPoint({ p[0],p[1],p[2] }));
      }
      auto polyloop = new typename Schema::IfcPolyLoop(polygon);
      auto face_bound = new typename Schema::IfcFaceOuterBound(polyloop, new typename Schema::IfcBoolean(true));
      typename Schema::IfcFaceBound::list::ptr face_bounds(new typename Schema::IfcFaceBound::list);
      face_bounds->push(face_bound);
      auto face = new typename Schema::IfcFace(face_bounds);
      faces->push(face);
   }
   auto shell = new typename Schema::IfcClosedShell(faces);
   auto faceted_brep = new typename Schema::IfcFacetedBrep(shell);

   typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);
   representation_items->push(faceted_brep);

   typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
   auto shape_representation = new typename Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("Brep"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   segment->setRepresentation(product_definition_shape);
}

template <typename Schema>
void CreateGirderSegmentRepresentation(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   const CPrecastSegmentData* pSegment = pIBridgeDesc->GetPrecastSegmentData(segmentKey);
   pgsTypes::SegmentVariationType variationType = pSegment->GetVariationType();

   GET_IFACE2(pBroker, IBridge, pBridge);

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   Float64 Ls = pBridge->GetSegmentLength(segmentKey);

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



   CComPtr<IAngle> angle_start_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metStart, &angle_start_face);
   CComPtr<IAngle> angle_end_face;
   pBridge->GetSegmentAngle(segmentKey, pgsTypes::metEnd, &angle_end_face);

   Float64 start_face_angle, end_face_angle;
   angle_start_face->get_Value(&start_face_angle);
   angle_end_face->get_Value(&end_face_angle);

   // linear interpolation of skew angle along length of segment
   auto fn_cut_angle = [start_face_angle, end_face_angle, Ls](Float64 x)->Float64 {
      return std::lerp(start_face_angle, end_face_angle, x / Ls);
      };


   // create the solid model
   // build the girder model in a simple coordinate system, then use segment->setObjectPlacement(segment_placement) to position beam in space
   switch (options.beam_model)
   {
   case CIfcExportOptions::BeamModel::SectionedSolidHorizontal:
      GirderSegment_SectionedSolidHorizontal(file, pBroker, segmentKey, vPoi, fn_cut_angle, segment, options, pGeometricRepresentationSubContext);
      break;
   case CIfcExportOptions::BeamModel::PolygonalFaceSet:
      GirderSegment_PolygonalFaceSet(file, pBroker, segmentKey, vPoi, fn_cut_angle, segment, options, pGeometricRepresentationSubContext);
      break;
   case CIfcExportOptions::BeamModel::FacetedBrep:
      GirderSegment_FacetedBrep(file, pBroker, segmentKey, vPoi, fn_cut_angle, segment, options, pGeometricRepresentationSubContext);
      break;

   default:
      ASSERT(false);
   }



   // Place the segment in 3D space

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

   WBFL::Geometry::Vector3d ref_direction(ex - sx, ey - sy, ez - sz); // along the length of the girder
   ref_direction.Normalize();
   WBFL::Geometry::Vector3d z(0, 0, 1); // true up direction
   WBFL::Geometry::Vector3d y = z.Cross(ref_direction); // cross product gives Y axis perpendicular to ref_direction and up
   WBFL::Geometry::Vector3d axis = ref_direction.Cross(y); // cross product gives Z axis of the girder
   
   
   typename Schema::IfcObjectPlacement* segment_placement = nullptr;
   if (options.beam_placement == CIfcExportOptions::BeamPlacement::Local)
   {
      segment_placement = file.addLocalPlacement(nullptr,
         sx, sy, sz,
         axis.X(), axis.Y(), axis.Z(),
         ref_direction.X(), ref_direction.Y(), ref_direction.Z());
   }
   else
   {
      auto directrix = GetAlignmentDirectrix(file, options);
      typename Schema::IfcCurve* basis_curve = nullptr;
      if (auto gc = directrix->as<typename Schema::IfcGradientCurve>())
      {
         basis_curve = gc->BaseCurve();
      }
      else
      {
         basis_curve = directrix;
      }

      Float64 startStation, startElevation, startGrade;
      auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

      Float64 station, offset;
      pBridge->GetStationAndOffset(poiStart, &station, &offset);
      // per PGSuper, positive offset is to the right, per IFC, positive value is to the left.... use -offset
      auto pde = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(station - startStation), -offset, sz, boost::none, basis_curve);
      auto a2pl = new typename Schema::IfcAxis2PlacementLinear(pde,
         new typename Schema::IfcDirection({axis.X(), axis.Y(), axis.Z()}),
         new typename Schema::IfcDirection({ ref_direction.X(),ref_direction.Y(),ref_direction.Z() })
         );

      auto fallback_placement = file.addPlacement3d(sx, sy, sz,
                                                    axis.X(), axis.Y(), axis.Z(),
                                                    ref_direction.X(), ref_direction.Y(), ref_direction.Z());

      segment_placement = new typename Schema::IfcLinearPlacement(nullptr, a2pl, fallback_placement);
      file.addEntity(segment_placement);
   }

   segment->setObjectPlacement(segment_placement);
}

template <typename Schema>
void CreateGirderSegmentMaterials(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, const CIfcExportOptions& options)
{
   USES_CONVERSION;
   GET_IFACE2(pBroker, IIntervals, pIntervals);
   GET_IFACE2(pBroker, IMaterials, pMaterials);
   GET_IFACE2(pBroker, IStrandGeometry, pStrandGeom);
   IntervalIndexType releaseIntervalIdx = pIntervals->GetPrestressReleaseInterval(segmentKey);
   IntervalIndexType liftingIntervalIdx = pIntervals->GetLiftSegmentInterval(segmentKey);
   IntervalIndexType haulingIntervalIdx = pIntervals->GetHaulSegmentInterval(segmentKey);

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
   CHECK(vPoi.size() == 1);
   const pgsPointOfInterest& poiMS = vPoi.front();

   Float64 fc = pMaterials->GetSegmentFc28(segmentKey);
   Float64 max_agg_size = pMaterials->GetSegmentMaxAggrSize(segmentKey);

   auto material = GetConcreteMaterial<Schema>(file, pBroker, fc, max_agg_size, "Precast Concrete", SEGMENT_BORDER_COLOR);

   // need a list of entities that are associated with this material
   // right now we are creating a unique material for each segment but we still need the list
   typename Schema::IfcDefinitionSelect::list::ptr segments(new typename Schema::IfcDefinitionSelect::list);
   segments->push(segment);

   // associate the material with the segment (ie segments collection)
   auto rel_associates_materials = new typename Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Concrete_to_Precast_Segment"), boost::none, segments, material);
   file.addEntity(rel_associates_materials);

   // Gather data for Pset_PrecastConcreteElementGeneral
   //GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* displacement_unit = nullptr;

   Float64 fci = pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx);
   Float64 fcl = pMaterials->GetSegmentFc(segmentKey, liftingIntervalIdx);
   Float64 fch = pMaterials->GetSegmentFc(segmentKey, haulingIntervalIdx);
   Float64 fpj = pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent);
   //if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
   //{
   //   stress_unit = GetStressUnit<Schema>(file,pBroker);
   //   displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);

   //   fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //   fci = WBFL::Units::ConvertFromSysUnits(fci, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //   fcl = WBFL::Units::ConvertFromSysUnits(fcl, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //   fch = WBFL::Units::ConvertFromSysUnits(fch, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //   fpj = WBFL::Units::ConvertFromSysUnits(fpj, pDisplayUnits->GetStressUnit().UnitOfMeasure);
   //}

   Float64 batter = 0.0;
   if (options.batter_ends)
   {
      GET_IFACE2(pBroker, IBridge, pBridge);
      Float64 slope = pBridge->GetSegmentSlope(segmentKey);
      batter = atan(slope);
   }

   // Pset_PrecastConcreteElementGeneral
   typename Schema::IfcProperty::list::ptr precast_concrete_properties(new typename Schema::IfcProperty::list);
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("FormStrippingStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LiftingStrength"), boost::none, new typename Schema::IfcPressureMeasure(fcl), stress_unit));
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("ReleaseStrength"), boost::none, new typename Schema::IfcPressureMeasure(fci), stress_unit));
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("TransportationStrength"), boost::none, new typename Schema::IfcPressureMeasure(fch), stress_unit));
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("InitialTension"), boost::none, new typename Schema::IfcPressureMeasure(fpj), stress_unit));
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BatterAtStart"), boost::none, new typename Schema::IfcPlaneAngleMeasure(batter), nullptr));
   precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("BatterAtEnd"), boost::none, new typename Schema::IfcPlaneAngleMeasure(batter), nullptr));
   if (options.include_camber) {
      // Compute the camber ratio
      // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_PrecastConcreteElementGeneral.htm
      // The camber deflection, measured from the midpoint of a cambered face of a piece to the midpoint of the chord joining the ends of the same face, 
      // as shown in the figure below (figure not provided), divided by the original (nominal) straight length of the face of the piece.
      GET_IFACE2(pBroker, ICamber, pCamber);
      Float64 D = pCamber->GetDCamberForGirderSchedule(poiMS, pgsTypes::CreepTime::Max);
      GET_IFACE2(pBroker, IBridge, pBridge);
      Float64 Ls = pBridge->GetSegmentPlanLength(segmentKey);
      Float64 camber_ratio = D / Ls;

      precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("CamberAtMidspan"), boost::none, new typename Schema::IfcRatioMeasure(camber_ratio), nullptr));
   }

   {
      pgsAutoGirderLabel autoLabel;
      pgsGirderLabel::UseAlphaLabel(false);
      precast_concrete_properties->push(new typename Schema::IfcPropertySingleValue(std::string("DesignLocationNumber"), boost::none, new typename Schema::IfcLabel(T2A(SEGMENT_LABEL(segmentKey))), nullptr));
   }
   auto pset_precast_concrete_element_general = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_PrecastConcreteElementGeneral"), boost::none, precast_concrete_properties);
   file.addEntity(pset_precast_concrete_element_general);

   typename Schema::IfcObjectDefinition::list::ptr related_segments(new typename Schema::IfcObjectDefinition::list);
   related_segments->push(segment);

   auto related_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, pset_precast_concrete_element_general);
   file.addEntity(related_properties);

   if (options.include_quantities)
   {
      // Qto_BeamBaseQuantities
#pragma Reminder("NOTE: These are a little bit dummy quantities - updated in the future")
   // assuming simple sections (no change in cross section or depth like end blocks are variable depth hammerhead segments)
   // need to update PGSuper so we can get the different surface areas directly instead of having to compute them here
      GET_IFACE2(pBroker, IBridge, pBridge);
      GET_IFACE2(pBroker, ISectionProperties, pSectProps);
      auto L = pBridge->GetSegmentPlanLength(segmentKey);
      auto A = pSectProps->GetAg(releaseIntervalIdx, poiMS);
      auto P = pSectProps->GetPerimeter(poiMS);
      auto OSA = L * P;
      auto GSA = OSA + 2 * A;
      auto GV = L * A;
      auto W = pSectProps->GetSegmentWeight(segmentKey);
      auto g = WBFL::Units::System::GetGravitationalAcceleration();
      auto Mass = W / g; // this is a unit of mass

      typename Schema::IfcConversionBasedUnit* big_area_unit = nullptr;
      typename Schema::IfcConversionBasedUnit* small_area_unit = nullptr;
      typename Schema::IfcConversionBasedUnit* volume_unit = nullptr;
      typename Schema::IfcConversionBasedUnit* mass_unit = nullptr;
      typename Schema::IfcConversionBasedUnit* length_unit = nullptr;

      //if (pDisplayUnits->GetUnitMode() == WBFL::EAF::UnitMode::US)
      //{
      //   big_area_unit = GetBigAreaUnit<Schema>(file, pBroker);
      //   small_area_unit = GetSmallAreaUnit<Schema>(file, pBroker);
      //   volume_unit = GetVolumeUnit<Schema>(file, pBroker);
      //   mass_unit = GetMassUnit<Schema>(file, pBroker);
      //   length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

      //   GSA = WBFL::Units::ConvertFromSysUnits(GSA, WBFL::Units::Measure::Feet2);
      //   GV = WBFL::Units::ConvertFromSysUnits(GV, WBFL::Units::Measure::Feet3);
      //   L = WBFL::Units::ConvertFromSysUnits(L, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
      //   A = WBFL::Units::ConvertFromSysUnits(A, pDisplayUnits->GetAreaUnit().UnitOfMeasure);
      //   Mass = WBFL::Units::ConvertFromSysUnits(Mass, WBFL::Units::Measure::PoundMass);
      //}


      typename Schema::IfcPhysicalQuantity::list::ptr beam_quantities(new typename Schema::IfcPhysicalQuantity::list);
      beam_quantities->push(new typename Schema::IfcQuantityArea(std::string("GrossSurfaceArea"), boost::none, big_area_unit, GSA, boost::none));
      beam_quantities->push(new typename Schema::IfcQuantityVolume(std::string("GrossVolume"), boost::none, volume_unit, GV, boost::none));

      auto qto_bodygeometryvalidation = new typename Schema::IfcElementQuantity(IfcParse::IfcGlobalId(), nullptr, std::string("Qto_BodyGeometryValidation"), boost::none, std::string("BaseQuantities"), beam_quantities);
      file.addEntity(qto_bodygeometryvalidation);

      beam_quantities->push(new typename Schema::IfcQuantityLength(std::string("Length"), boost::none, length_unit, L, boost::none));
      beam_quantities->push(new typename Schema::IfcQuantityArea(std::string("CrossSectionArea"), boost::none, small_area_unit, A, boost::none));
      beam_quantities->push(new typename Schema::IfcQuantityArea(std::string("OuterSurfaceArea"), boost::none, big_area_unit, OSA, boost::none));
      beam_quantities->push(new typename Schema::IfcQuantityWeight(std::string("GrossWeight"), boost::none, mass_unit, Mass, boost::none));

      auto qto_beambasequantities = new typename Schema::IfcElementQuantity(IfcParse::IfcGlobalId(), nullptr, std::string("Qto_BeamBaseQuantities"), boost::none, std::string("BaseQuantities"), beam_quantities);
      file.addEntity(qto_beambasequantities);

      auto rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, qto_bodygeometryvalidation);
      file.addEntity(rel_defines_by_properties);

      rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, qto_beambasequantities);
      file.addEntity(rel_defines_by_properties);
   }
}

template <typename Schema>
void CreateStrandRepresentation(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* beam, const CIfcExportOptions& options)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   const pgsPointOfInterest& poiStart(vPoi.front());
   const pgsPointOfInterest& poiEnd(vPoi.back());

   auto strands = CreateStrands<Schema>(file, pBroker, poiStart, poiEnd, beam);

   if (0 < strands->size())
   {
      pgsTypes::StrandType strandType = pgsTypes::Straight;
      GET_IFACE2(pBroker, IMaterials, pMaterials);
      const auto* pStrand = pMaterials->GetStrandMaterial(segmentKey, strandType);
      auto strand_material = GetStrandMaterial(file, pStrand);

      typename Schema::IfcDefinitionSelect::list::ptr strands_for_material(new typename Schema::IfcDefinitionSelect::list);
      for (auto& strand : *strands)
      {
         strands_for_material->push(strand);
      }

      // associate the material with the strand
      auto rel_associates_materials = new typename Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Steel_to_Strand"), boost::none, strands_for_material, strand_material);
      file.addEntity(rel_associates_materials);

      // strands are a aggregate part of a beam
      auto rel_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Segment_Aggregates_Strands"), boost::none, beam, strands);
      file.addEntity(rel_aggregates);
   }
}

template <typename Schema>
void CreateLongitudinalRebarRepresentation(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* beam,typename Schema::IfcElementAssembly* rebar_assembly, const CIfcExportOptions& options)
{
   if (!options.include_rebar)
      return;

   USES_CONVERSION;

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   const pgsPointOfInterest& poiStart(vPoi.front());
   const pgsPointOfInterest& poiEnd(vPoi.back());


   auto rebars = CreateRebars<Schema>(file, pBroker, poiStart, poiEnd, beam);

   if (0 < rebars->size())
   {
      GET_IFACE2(pBroker, IMaterials, pMaterials);
      WBFL::Materials::Rebar::Type rebar_type;
      WBFL::Materials::Rebar::Grade rebar_grade;
      pMaterials->GetSegmentLongitudinalRebarMaterial(segmentKey, &rebar_type, &rebar_grade);
      const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(rebar_type, rebar_grade, WBFL::Materials::Rebar::Size::bs3);
      auto rebar_material = GetRebarMaterial(file,pRebar, "Rebar", REBAR_COLOR);

      // need a list of entities that are associated with this material
      // right now we are creating a unique material for each strand but we still need the list
      typename Schema::IfcDefinitionSelect::list::ptr rebars_for_material(new typename Schema::IfcDefinitionSelect::list);
      for (auto& rebar : *rebars)
      {
         rebars_for_material->push(rebar);
      }

      // associate the material with the rebar
      auto rel_associates_materials = new typename Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Steel_to_Rebar"), boost::none, rebars_for_material, rebar_material);
      file.addEntity(rel_associates_materials);

      // aggregate the rebar with its assembly
      auto rel_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Rebar_Assembly_Aggregates_Rebars"), boost::none, rebar_assembly, rebars);
      file.addEntity(rel_aggregates);
   }
}


template <typename Schema>
void CreateStirrupRepresentation(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* beam, typename Schema::IfcElementAssembly* rebar_assembly, const CIfcExportOptions& options)
{
   if (!options.include_rebar)
      return;

   USES_CONVERSION;

   // For now, we only do stirrups for WF-Beams (that's because stirrups are dummy rebars)
   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   const CBridgeDescription2* pBridgeDesc = pIBridgeDesc->GetBridgeDescription();
   const CGirderGroupData* pGroup = pBridgeDesc->GetGirderGroup(segmentKey.groupIndex);
   const GirderLibraryEntry* pGdrEntry = pGroup->GetGirderLibraryEntry(segmentKey.girderIndex);
   auto beam_factory = pGdrEntry->GetBeamFactory();
   if (!::IsEqualGUID(beam_factory->GetFamilyCLSID(), CLSID_WFBeamFamily))
      return;


   auto rebars = CreateStirrups<Schema>(file, pBroker, segmentKey, beam);

   if (0 < rebars->size())
   {
      GET_IFACE2(pBroker, IMaterials, pMaterials);
      WBFL::Materials::Rebar::Type rebar_type;
      WBFL::Materials::Rebar::Grade rebar_grade;
      pMaterials->GetSegmentLongitudinalRebarMaterial(segmentKey, &rebar_type, &rebar_grade);
      const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(rebar_type, rebar_grade, WBFL::Materials::Rebar::Size::bs3);

      auto rebar_material = GetRebarMaterial(file, pRebar, "Stirrup", STIRRUP_COLOR);

      // need a list of entities that are associated with this material
      // right now we are creating a unique material for each strand but we still need the list
      typename Schema::IfcDefinitionSelect::list::ptr rebars_for_material(new typename Schema::IfcDefinitionSelect::list);
      for (auto& rebar : *rebars)
      {
         rebars_for_material->push(rebar);
      }

      // associate the material with the segment (ie segments collection)
      auto rel_associates_materials = new typename Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Steel_to_Rebar"), boost::none, rebars_for_material, rebar_material);
      file.addEntity(rel_associates_materials);

      // aggregate the rebar with its assembly
      auto rel_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Rebar_Assembly_Aggregates_Rebars"), boost::none, rebar_assembly, rebars);
      file.addEntity(rel_aggregates);
   }
}


template <typename Schema>
void CreateClosureJointRepresentation(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CClosureKey& closureKey, typename Schema::IfcBeam* closureJoint, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
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
   typename Schema::IfcCartesianPoint::list::ptr girder_line_points(new typename Schema::IfcCartesianPoint::list);
   girder_line_points->push(new typename Schema::IfcCartesianPoint(std::vector<double>{0, 0, 0}));
   girder_line_points->push(new typename Schema::IfcCartesianPoint(std::vector<double>{Lc, 0, 0}));
   auto girder_line = new typename Schema::IfcPolyline(girder_line_points);
   file.addEntity(girder_line);

   typename Schema::IfcProfileDef::list::ptr cross_sections(new typename Schema::IfcProfileDef::list);
   cross_sections->push(CreateSectionProfile<Schema>(pShapes, poiStart, intervalIdx, options));
   cross_sections->push(CreateSectionProfile<Schema>(pShapes, poiEnd, intervalIdx, options));

   typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);

   std::string representation_type;
   typename Schema::IfcAxis2PlacementLinear::list::ptr cross_section_positions(new typename Schema::IfcAxis2PlacementLinear::list);
   auto pde_start = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(0.0), boost::none, boost::none, boost::none, girder_line);
   auto pde_end = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(Lc), boost::none, boost::none, boost::none, girder_line);
   auto start_section = new typename Schema::IfcAxis2PlacementLinear(pde_start, nullptr, nullptr);
   auto end_section = new typename Schema::IfcAxis2PlacementLinear(pde_end, nullptr, nullptr);
   cross_section_positions->push(start_section);
   cross_section_positions->push(end_section);

   auto sectioned_solid = new typename Schema::IfcSectionedSolidHorizontal(girder_line, cross_sections, cross_section_positions);
   file.addEntity(sectioned_solid);
   representation_items->push(sectioned_solid);

   typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
   auto shape_representation = new typename Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

   // Place the segment in 3D space
   WBFL::Geometry::Vector3d ref_direction(ex - sx, ey - sy, ez - sz); // along the length of the girder
   ref_direction.Normalize();
   WBFL::Geometry::Vector3d z(0, 0, 1); // true up direction
   WBFL::Geometry::Vector3d y = z.Cross(ref_direction); // cross product gives Y axis perpendicular to ref_direction and up
   WBFL::Geometry::Vector3d axis = ref_direction.Cross(y); // cross product gives Z axis of the girder

   typename Schema::IfcObjectPlacement* closure_placement = nullptr;
   if (options.beam_placement == CIfcExportOptions::BeamPlacement::Local)
   {
      closure_placement = file.addLocalPlacement(nullptr,
         sx, sy, sz,
         axis.X(), axis.Y(), axis.Z(),
         ref_direction.X(), ref_direction.Y(), ref_direction.Z());
   }
   else
   {
      auto directrix = GetAlignmentDirectrix(file, options);
      typename Schema::IfcCurve* basis_curve = nullptr;
      if (auto gc = directrix->as<typename Schema::IfcGradientCurve>())
      {
         basis_curve = gc->BaseCurve();
      }
      else
      {
         basis_curve = directrix;
      }

      Float64 startStation, startElevation, startGrade;
      auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

      Float64 station, offset;
      pBridge->GetStationAndOffset(poiStart, &station, &offset);
      auto pde = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(station - startStation), offset, sz, boost::none, basis_curve);
      auto a2pl = new typename Schema::IfcAxis2PlacementLinear(pde,
         new typename Schema::IfcDirection({ axis.X(), axis.Y(), axis.Z() }),
         new typename Schema::IfcDirection({ ref_direction.X(),ref_direction.Y(),ref_direction.Z() })
      );

      auto fallback_placement = file.addPlacement3d(sx, sy, sz,
         axis.X(), axis.Y(), axis.Z(),
         ref_direction.X(), ref_direction.Y(), ref_direction.Z());

      closure_placement = new typename Schema::IfcLinearPlacement(nullptr, a2pl, fallback_placement);
      file.addEntity(closure_placement);
   }
   closureJoint->setObjectPlacement(closure_placement);
   closureJoint->setRepresentation(product_definition_shape);
}

template <typename Schema>
void CreateSlab(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcBridgePart* deck, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
#pragma Reminder("WORKING HERE - Deck Model - need to re-think this approach")
   // Consider modeling the slab separately from the haunch. The basic slab is the same everywhere.
   // Each girder has it's own haunch.
   // This should eliminate the problem with different number of points in the cross section profile.
   // The deck representation would be a composite of the main slab and each haunch (This is what Representation.Items is for. Deck + each haunch are items)

   // This is not a good model of the deck. This model just creates NUM_DECK_SECTIONS cross sections and extrudes between them.
   GET_IFACE2(pBroker, IBridge, pBridge);
   USES_CONVERSION;

   if (pBridge->GetDeckType() == pgsTypes::sdtNone)
      return;

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);

   CComPtr<IAngle> objAngle;
   pBridge->GetPierSkew(0, &objAngle);
   Float64 start_skew;
   objAngle->get_Value(&start_skew);
   objAngle.Release();
   pBridge->GetPierSkew(pBridge->GetPierCount() - 1, &objAngle);
   Float64 end_skew;
   objAngle->get_Value(&end_skew);

   auto fn_cut_angle = [start_skew, end_skew, L = endBrgStation - startBrgStation](Float64 x)->Float64 { return PI_OVER_2 + std::lerp(start_skew, end_skew, x / L); };

   // get the directrix line of the alignment
   auto directrix = GetAlignmentDirectrix(file,options);

   Float64 startStation, startElevation, startGrade;
   auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   IndexType nDeckSections = NUM_DECK_SECTIONS;

   IndexType point_count = 0;

   bool bIncludeHaunch = IsZero(start_skew) && IsZero(end_skew) ? true : false;

   GET_IFACE2(pBroker, IShapes, pShapes);
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   typename Schema::IfcProfileDef::list::ptr cross_sections(new typename Schema::IfcProfileDef::list);
   typename Schema::IfcAxis2PlacementLinear::list::ptr cross_section_positions(new typename Schema::IfcAxis2PlacementLinear::list);
   for (IndexType i = 0; i <= nDeckSections; i++)
   {
      auto station = i * (endBrgStation - startBrgStation) / nDeckSections + startBrgStation;

      // This code, and the objDir in GetSlabShape, are trying to account for skew by sweeping the cut line angle
      // between the start and end of the bridge. However, there appears to be an issue in WBFL::CoordinateGeometry 
      // that causes the top of deck section from the roadway to have duplicate points. Duplicate points are
      // not valid for the IFC section shape so we will just use a normal section cut for now.
      // NOTE: Below, we skew the normal cut section polygon similar to what we do for girders.
      //auto dir = i * (endDir - startDir) / nDeckSections + startDir;
      //objDir->put_Value(dir);



      CComPtr<IShape> slab_shape;
      pShapes->GetSlabShape(station, nullptr/*objDir*/, bIncludeHaunch, &slab_shape);

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

      auto cut_angle = fn_cut_angle(station);
      auto polyline = CreatePolyline<Schema>(slab_shape, options, cut_angle);

      CComPtr<IDirection> objDir;
      pAlignment->GetBearingNormal(station, &objDir);
      Float64 normal_angle;
      objDir->get_Value(&normal_angle);

      std::ostringstream os;
      os << "Deck Section at Station " << T2A(WBFL::COGO::Station(station).AsString(station_format).c_str());
      auto deck_perimeter = new typename Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, os.str(), polyline);
      cross_sections->push(deck_perimeter);
      file.addEntity(deck_perimeter);

      auto distance_along = station - startStation;
      auto pde = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(distance_along), boost::none, boost::none, boost::none, directrix);
      auto rd = new typename Schema::IfcDirection({ cos(normal_angle + cut_angle),sin(normal_angle + cut_angle),0.0 }); // normal to the plane of the cross section
      auto axis = new typename Schema::IfcDirection({ 0,0,1 }); // up
      auto deck_section_placement = new typename Schema::IfcAxis2PlacementLinear(pde, axis, rd);
      cross_section_positions->push(deck_section_placement);
      file.addEntity(pde);
      file.addEntity(deck_section_placement);
   }

   typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);
   auto sectioned_solid = new typename Schema::IfcSectionedSolidHorizontal(directrix, cross_sections, cross_section_positions);
   representation_items->push(sectioned_solid);
   file.addEntity(sectioned_solid);

   auto site = file.getSingle<typename Schema::IfcSite>();
   auto deck_placement = site->ObjectPlacement();

   typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
   auto shape_representation = new typename Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   shape_representation_list->push(shape_representation);
   auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
   file.addEntity(product_definition_shape);

   auto slab = new typename Schema::IfcSlab(IfcParse::IfcGlobalId(), nullptr, std::string("Deck Slab"), boost::none, boost::none, deck_placement, product_definition_shape, boost::none,
      Schema::IfcSlabTypeEnum::IfcSlabType_FLOOR); // see Ifc 4x3 6.1.2.19.2 (FLOOR represents a bridge deck), name is option but AASHTO IDS requires it
   file.addEntity(slab);

   if (options.classify)
   {
      Classify_Slab(file, slab);
   }

   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   const auto* deck_desc = pIBridgeDesc->GetDeckDescription();

   typename Schema::IfcProperty::list::ptr deck_properties(new typename Schema::IfcProperty::list);
   deck_properties->push(new typename Schema::IfcPropertySingleValue(std::string("GrossDepth"), boost::none, new typename Schema::IfcLengthMeasure(deck_desc->GrossDepth), nullptr));
   deck_properties->push(new typename Schema::IfcPropertySingleValue(std::string("LeftEdgeDepth"), boost::none, new typename Schema::IfcLengthMeasure(deck_desc->OverhangEdgeDepth[pgsTypes::stLeft]), nullptr));
   deck_properties->push(new typename Schema::IfcPropertySingleValue(std::string("RightEdgeDepth"), boost::none, new typename Schema::IfcLengthMeasure(deck_desc->OverhangEdgeDepth[pgsTypes::stRight]), nullptr));

   typename Schema::IfcValue::list::ptr station_list(new typename Schema::IfcValue::list);
   typename Schema::IfcValue::list::ptr left_list(new typename Schema::IfcValue::list);
   typename Schema::IfcValue::list::ptr right_list(new typename Schema::IfcValue::list);
   for (const auto& deck_point : deck_desc->DeckEdgePoints)
   {
      station_list->push(new typename Schema::IfcLengthMeasure(deck_point.Station));
      left_list->push(new typename Schema::IfcLengthMeasure(deck_point.LeftEdge));   
      right_list->push(new typename Schema::IfcLengthMeasure(deck_point.RightEdge));
   }
   deck_properties->push(new typename Schema::IfcPropertyListValue(std::string("Stations"), boost::none, station_list, nullptr));
   deck_properties->push(new typename Schema::IfcPropertyListValue(std::string("LeftEdges"), boost::none, left_list, nullptr));
   deck_properties->push(new typename Schema::IfcPropertyListValue(std::string("RightEdges"), boost::none, right_list, nullptr));

   auto pset_deck = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("pgsDeck"), boost::none, deck_properties);
   file.addEntity(pset_deck);
   typename Schema::IfcObjectDefinition::list::ptr decks(new typename Schema::IfcObjectDefinition::list);
   decks->push(slab);
   auto rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Defines slab edge offset geometry"), boost::none, decks, pset_deck);
   file.addEntity(rel_defines_by_properties);

   file.addRelatedObject<typename Schema::IfcRelContainedInSpatialStructure>(deck, slab);
}


template <typename Schema>
void CreateRailingSystemRepresentation(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, pgsTypes::TrafficBarrierOrientation tbOrientation, typename Schema::IfcProduct* railing, const CIfcExportOptions& options, typename Schema::IfcGeometricRepresentationSubContext* pGeometricRepresentationSubContext)
{
   GET_IFACE2(pBroker, IBarriers, pBarriers);

   bool bHasSidewalk = pBarriers->HasSidewalk(tbOrientation);
   bool bHasInteriorBarrier = pBarriers->HasInteriorBarrier(tbOrientation);

   IndexType nShapesPerBarrier = 1 + (bHasSidewalk ? 1 : 0) + (bHasInteriorBarrier ? 1 : 0);

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 startBrgStation = pBridge->GetBearingStation(0, pgsTypes::Ahead);
   Float64 endBrgStation = pBridge->GetBearingStation(pBridge->GetPierCount() - 1, pgsTypes::Back);


   CComPtr<IAngle> objAngle;
   pBridge->GetPierSkew(0, &objAngle);
   Float64 start_skew;
   objAngle->get_Value(&start_skew);
   objAngle.Release();
   pBridge->GetPierSkew(pBridge->GetPierCount() - 1, &objAngle);
   Float64 end_skew;
   objAngle->get_Value(&end_skew);

   auto fn_cut_angle = [start_skew, end_skew, L = endBrgStation - startBrgStation](Float64 x)->Float64 { return PI_OVER_2 + std::lerp(start_skew, end_skew, x / L); };

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
   auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   typename Schema::IfcAxis2PlacementLinear::list::ptr cross_section_positions(new typename Schema::IfcAxis2PlacementLinear::list);
   std::vector<typename Schema::IfcProfileDef::list::ptr> cross_sections;
   for (int i = 0; i < nShapesPerBarrier; i++)
   {
      typename Schema::IfcProfileDef::list::ptr ptr(new typename Schema::IfcProfileDef::list);
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

      auto cut_angle = fn_cut_angle(station);


      CComPtr<IDirection> objDir;
      pAlignment->GetBearingNormal(station, &objDir);
      Float64 normal_angle;
      objDir->get_Value(&normal_angle);

      auto distance_along = station - startStation;
      auto pde = new typename Schema::IfcPointByDistanceExpression(new typename Schema::IfcLengthMeasure(distance_along), boost::none, boost::none, boost::none, directrix);
      auto rd = new typename Schema::IfcDirection({ cos(normal_angle + cut_angle),sin(normal_angle + cut_angle),0.0 }); // normal to the plane of the cross section
      auto axis = new typename Schema::IfcDirection({ 0,0,1 }); // up
      auto placement = new typename Schema::IfcAxis2PlacementLinear(pde, axis, rd);
      cross_section_positions->push(placement);
      file.addEntity(pde);
      file.addEntity(placement);

      for (IndexType shapeIdx = 0; shapeIdx < nShapesPerBarrier; shapeIdx++)
      {
         CComPtr<ICompositeShapeItem> shape_item;
         composite->get_Item(shapeIdx, &shape_item);

         CComPtr<IShape> shape;
         shape_item->get_Shape(&shape);

         auto polyline = CreatePolyline<Schema>(shape, options, cut_angle);
         if (polyline)
         {
            std::ostringstream os;
            os << (tbOrientation == pgsTypes::tboLeft ? "Left" : "Right") << " Barrier";
            if (1 < nShapesPerBarrier)
               os << " Shape " << shapeIdx;

            auto shape_perimeter = new typename Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, os.str(), polyline);
            cross_sections[shapeIdx]->push(shape_perimeter);
         }
      } // next shape
   } // next section

   typename Schema::IfcRepresentationItem::list::ptr representation_items(new typename Schema::IfcRepresentationItem::list);
   for (IndexType shapeIdx = 0; shapeIdx < nShapesPerBarrier; shapeIdx++)
   {
      if (0 < cross_sections[shapeIdx]->size())
      {
         auto sectioned_solid = new typename Schema::IfcSectionedSolidHorizontal(directrix, cross_sections[shapeIdx], cross_section_positions);
         representation_items->push(sectioned_solid);
         file.addEntity(sectioned_solid);
      }
   }

   if (0 < representation_items->size())
   {
      auto site = file.getSingle<typename Schema::IfcSite>();
      auto railing_placement = site->ObjectPlacement();


      typename Schema::IfcRepresentation::list::ptr shape_representation_list(new typename Schema::IfcRepresentation::list);
      auto shape_representation = new typename Schema::IfcShapeRepresentation(pGeometricRepresentationSubContext, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
      shape_representation_list->push(shape_representation);
      auto product_definition_shape = new typename Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);
      file.addEntity(product_definition_shape);

      railing->setObjectPlacement(railing_placement);
      railing->setRepresentation(product_definition_shape);
   }
}

template <typename Schema>
typename Schema::IfcObjectDefinition::list::ptr CreatePiers(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   USES_CONVERSION;

   typename Schema::IfcObjectDefinition::list::ptr list_of_piers(new typename Schema::IfcObjectDefinition::list);

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   auto station_format = pDisplayUnits->GetStationFormat();

   // get stationing information
   Float64 startStation, startElevation, startGrade;
   auto startPoint = GetAlignmentStartPoint(pBroker, &startStation, &startElevation, &startGrade);

   auto directrix = GetAlignmentDirectrix(file, options);

   GET_IFACE2(pBroker, IBridge, pBridge);
   auto nPiers = pBridge->GetPierCount();
   for (IndexType pierIdx = 0; pierIdx < nPiers; pierIdx++)
   {
      std::string pier_name(T2A(LABEL_PIER_EX(pBridge->IsAbutment(pierIdx), pierIdx)));
      auto pier = new typename Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, pier_name, boost::none, boost::none, nullptr, nullptr, boost::none,
         Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
         Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
         pBridge->IsAbutment(pierIdx) ? Schema::IfcBridgePartTypeEnum::IfcBridgePartType_ABUTMENT : Schema::IfcBridgePartTypeEnum::IfcBridgePartType_PIER);
      file.addEntity(pier);

      if (options.classify)
      {
         if (pBridge->IsAbutment(pierIdx))
            Classify_Abutment<Schema>(file, pier);
         else
            Classify_Pier<Schema>(file, pier);
      }

      std::ostringstream os;
      os << "Foundation at " << pier_name; // name is not required, but is specified in AASHTO IDS
      auto foundation = new typename Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, nullptr, nullptr, boost::none,
         Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
         Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
         Schema::IfcBridgePartTypeEnum::IfcBridgePartType_FOUNDATION);
      file.addEntity(foundation);
      if (options.classify)
      {
         Classify_Foundation<Schema>(file, foundation);
      }

      typename Schema::IfcObjectDefinition::list::ptr list_of_foundations(new typename Schema::IfcObjectDefinition::list);
      list_of_foundations->push(foundation);

      // IfcBridgePart::PIER <-> IfcRelAggregates <-> IfcBridgePart::FOUNDATION
      auto rel_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Foundation is an aggregate component of pier"), boost::none, pier, list_of_foundations);
      file.addEntity(rel_aggregates);

      if (options.include_work_plan)
      {
         // related the Stage 1 construction task to the pier product
         auto task = GetStage1Task(file, pBroker, options);
         AssignTaskToProduct(task, pier, file, pBroker, options);
      }

      // referent position
      auto pierStation = pBridge->GetPierStation(pierIdx);

      auto point_on_alignment = new typename Schema::IfcPointByDistanceExpression(
         new typename Schema::IfcLengthMeasure(pierStation - startStation),
         boost::none, boost::none, boost::none,
         directrix);
      auto relative_placement = new typename Schema::IfcAxis2PlacementLinear(point_on_alignment, nullptr, nullptr);
      auto referent_placement = new typename Schema::IfcLinearPlacement(nullptr, relative_placement, nullptr);

      // create referent
      std::ostringstream os2;
      os2 << "Station " << T2A(WBFL::COGO::Station(pierStation).AsString(station_format).c_str()) << " " << T2A(LABEL_PIER_EX(pBridge->IsAbutment(pierIdx), pierIdx));
      auto referent = new typename Schema::IfcReferent(IfcParse::IfcGlobalId(), nullptr, os2.str(), boost::none, boost::none, referent_placement, nullptr, Schema::IfcReferentTypeEnum::IfcReferentType_POSITION);
      file.addEntity(referent);

      // create and assign Pset_Stationing
      typename Schema::IfcProperty::list::ptr pset_station_properties(new typename Schema::IfcProperty::list);
      pset_station_properties->push(new typename Schema::IfcPropertySingleValue(std::string("Station"), boost::none, new typename Schema::IfcLengthMeasure(pierStation), nullptr));

      auto property_set = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_Stationing"), boost::none, pset_station_properties);
      file.addEntity(property_set);

      typename Schema::IfcObjectDefinition::list::ptr referents(new typename Schema::IfcObjectDefinition::list);
      referents->push(referent);

      auto rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Relates pier station properties to referent"), boost::none, referents, property_set);
      file.addEntity(rel_defines_by_properties);

      // IfcReferent <-> IfcRelPositions <-> IfcBridgePart::PIER,FOUNDATION
      // Without providing geometry, this is how the pier and foundation are positions - referent informs on the position of the products it positions
      std::string strPositions("Positions pier and foundation");
      typename Schema::IfcProduct::list::ptr related_products(new typename Schema::IfcProduct::list);
      related_products->push(pier);
      related_products->push(foundation);
      if (pierIdx == 0 || pierIdx == nPiers - 1)
      {
         auto bridge = file.getSingle<typename Schema::IfcBridge>();
         related_products->push(bridge);
         strPositions = "Positions ends of bridge";
      }

      auto rel_positions = new typename Schema::IfcRelPositions(IfcParse::IfcGlobalId(), nullptr, strPositions, boost::none, referent, related_products);
      file.addEntity(rel_positions);

      auto alignment = file.getSingle<typename Schema::IfcAlignment>();
      AddReferent(file, alignment, referent);

      // WORKING HERE - experimenting with custom property sets for girder spacing
      typename Schema::IfcProperty::list::ptr spacing_property(new typename Schema::IfcProperty::list);
      if ( 0 < pierIdx )
      {
         typename Schema::IfcValue::list::ptr list(new typename Schema::IfcValue::list);

         auto spacing = pBridge->GetGirderSpacing(pierIdx, pgsTypes::PierFaceType::Back, pgsTypes::MeasurementLocation::AtCenterlineBearing, pgsTypes::MeasurementType::NormalToItem);
         for (auto s : spacing)
         {
            list->push(new typename Schema::IfcLengthMeasure(s));
         }
         spacing_property->push(new typename Schema::IfcPropertyListValue(std::string("Back_Spacing"), boost::none, list, nullptr));
      }

      if (pierIdx < nPiers - 1)
      {
         typename Schema::IfcValue::list::ptr list(new typename Schema::IfcValue::list);

         auto spacing = pBridge->GetGirderSpacing(pierIdx, pgsTypes::PierFaceType::Ahead, pgsTypes::MeasurementLocation::AtCenterlineBearing, pgsTypes::MeasurementType::NormalToItem);
         for (auto s : spacing)
         {
            list->push(new typename Schema::IfcLengthMeasure(s));
         }
         spacing_property->push(new typename Schema::IfcPropertyListValue(std::string("Ahead_Spacing"), boost::none, list, nullptr));
      }
      auto pset_spacing = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("pgsSpacing"), boost::none, spacing_property);
      file.addEntity(pset_spacing);
      typename Schema::IfcObjectDefinition::list::ptr piers(new typename Schema::IfcObjectDefinition::list);
      piers->push(pier);
      rel_defines_by_properties = new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, std::string("Relates girder spacing property set to pier"), boost::none, piers, pset_spacing);
      file.addEntity(rel_defines_by_properties);


      list_of_piers->push(pier);
   }

   return list_of_piers;
}

template <typename Schema>
typename Schema::IfcBeam* CreatePrecastSegment(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const std::string& name,const CSegmentKey& segmentKey,typename Schema::IfcBeamType* beam_type,const CIfcExportOptions& options)
{
   auto beam = new IfcSchema::IfcBeam(
      IfcParse::IfcGlobalId(),
      nullptr, // OwnerHistory
      name,  // Name
      boost::none, // Description
      boost::none, // ObjectType
      nullptr,  // ObjectPlacement
      nullptr,  // Representation
      boost::none, // Tag
      boost::none // PredefinedType (must not be used if defined in IfcBeamType)
   );

   auto body_model_representation_subcontext = file.getRepresentationSubContext(std::string("Body"), std::string("Model"));

   CreateGirderSegmentRepresentation<Schema>(file, pBroker, segmentKey, beam, options, body_model_representation_subcontext);
   CreateGirderSegmentMaterials<Schema>(file, pBroker, segmentKey, beam, options);
   CreateStrandRepresentation<Schema>(file, pBroker, segmentKey, beam, options);

   auto rebar_assembly = new typename Schema::IfcElementAssembly(
      IfcParse::IfcGlobalId(),
      nullptr, // OwnerHistory
      std::string("Girder Rebar"), // Name
      boost::none, // Description
      boost::none, // ObjectType
      nullptr, // ObjectPlacement
      nullptr, // Representation
      boost::none, // Tag
      Schema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_FACTORY, // AssemblyPlace
      Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_REINFORCEMENT_UNIT // PredefinedType
   );

   // aggregate the rebar assembly with its girder segment
   typename Schema::IfcObjectDefinition::list::ptr beam_aggregate_elements(new typename Schema::IfcObjectDefinition::list);
   beam_aggregate_elements->push(rebar_assembly);
   auto rel_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Girder_Segment_Aggregates_Rebar_Assembly"), boost::none, beam, beam_aggregate_elements);
   file.addEntity(rel_aggregates);

   CreateLongitudinalRebarRepresentation<Schema>(file, pBroker, segmentKey, beam, rebar_assembly, options);
   CreateStirrupRepresentation<Schema>(file, pBroker, segmentKey, beam, rebar_assembly, options);

   // do this last so the beam is complete defined before it is put into the model
   file.addEntity(beam);
   file.addRelatedObject<typename Schema::IfcRelDefinesByType>(beam_type, beam);
   return beam;
}

template <typename Schema>
void CreateBridge(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CIfcExportOptions& options)
{
   USES_CONVERSION;

   if (options.classify)
   {
      Add_USBridge_Classification<Schema>(file);
   }


   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   ATLASSERT(geometric_representation_context);

   auto body_model_representation_subcontext = new typename Schema::IfcGeometricRepresentationSubContext(std::string("Body"), std::string("Model"), geometric_representation_context, boost::none, Schema::IfcGeometricProjectionEnum::IfcGeometricProjection_MODEL_VIEW, boost::none);
   file.addEntity(body_model_representation_subcontext);


   // Define spatial structure
   // From IfcSite https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/IfcSite.htm
   // IfcProject <-> IfcRelAggregates <-> IfcSite <-> IfcRelAggregates <-> IfcBridge
   // InitializeFile sets up IfcProject <-> IfcRelAggregates <-> IfcSite



   // Create IfcBridge, which is an IfcSpatialStructureElement
   GET_IFACE2(pBroker, IProjectProperties, pProjectProperties);
   std::string bridge_name(T2A(pProjectProperties->GetBridgeName()));
   if (bridge_name.empty()) bridge_name = "Unnamed Bridge";

   auto bridge = new typename Schema::IfcBridge(IfcParse::IfcGlobalId(), nullptr, bridge_name, boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcElementCompositionEnum::IfcElementComposition_COMPLEX, Schema::IfcBridgeTypeEnum::IfcBridgeType_GIRDER);
   file.addEntity(bridge);

   // create the assumed construction sequence and related it to the bridge
   CreateAssumedConstructionSequence(file, pBroker, options); // must come after IfcBridge is created because this function looks up the bridge

   // create a list of bridges in the site
   typename Schema::IfcObjectDefinition::list::ptr list_of_bridges_in_the_site(new typename Schema::IfcObjectDefinition::list);
   list_of_bridges_in_the_site->push(bridge); // add the bridge to the list

   // aggregate the bridges with the side
   // IfcSite <-> IfcRelAggregates <-> IfcBridge
   auto site = file.getSingle<typename Schema::IfcSite>();
   auto rel_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Bridges in the site"), boost::none, site, list_of_bridges_in_the_site);
   file.addEntity(rel_aggregates);



   // Create top level spatial structure of bridge
   auto superstructure = new typename Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Superstructure"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);
   file.addEntity(superstructure);
   if (options.classify)
   {
      Classify_Superstructure<Schema>(file, superstructure);
   }

   auto substructure = new typename Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Substructure"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUBSTRUCTURE);
   file.addEntity(substructure);
   if (options.classify)
   {
      Classify_Substructure<Schema>(file, substructure);
   }

   auto deck = new typename Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Deck"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_DECK);
   CreateSlab(file, pBroker, deck, options, body_model_representation_subcontext);
   file.addEntity(deck);
   if (options.classify)
   {
      Classify_Deck<Schema>(file, deck);
   }

   typename Schema::IfcObjectDefinition::list::ptr list_of_bridge_parts(new typename Schema::IfcObjectDefinition::list);
   list_of_bridge_parts->push(deck);
   list_of_bridge_parts->push(superstructure);
   list_of_bridge_parts->push(substructure);

   // IfcBridge <-> IfcRelAggregates <-> IfcBridgePart::DECK, SUPERSTRUCTURE, SUBSTRUCTURE
   auto bridge_spatial_elements = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in bridge spatial structure"), boost::none, bridge, list_of_bridge_parts);
   file.addEntity(bridge_spatial_elements);


   // Create spatial structure of substructure
   // IfcBridgePart::SUBSTRUCTURE <-> IfcRelAggregates <-> IfcBridgePart::ABUTMENT, PIER
   auto list_of_piers = CreatePiers(file, pBroker, options);
   auto substructure_spatial_elements = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in substructure spatial structure"), boost::none, substructure, list_of_piers);
   file.addEntity(substructure_spatial_elements);

   if (options.include_work_plan)
   {
      // related the Stage 1 construction task to the pier product
      auto task = GetStage3Task(file, pBroker, options);
      auto rel_contained_in_spatial_structures = deck->ContainsElements();
      for (auto& rel_contained_in_spatial_structure : *rel_contained_in_spatial_structures)
      {
         auto related_elements = rel_contained_in_spatial_structure->RelatedElements();
         for (auto& related_element : *related_elements)
         {
            AssignTaskToProduct(task, related_element, file, pBroker, options);
         }
      }
   }


   // Add railings to the spatial structure of the superstructure
   // IfcBridgePart::DECK <-> IfcRelContainedInSpatialStructure <-> IfcRailing

#pragma Reminder("PGSUPER needs the concept of No Railing - here we define a railing product with no representation which says there is a railing")
   typename Schema::IfcProduct* left_railing;
   if (options.railings == CIfcExportOptions::Railings::Parapet)
   {
      left_railing = new typename Schema::IfcWall(IfcParse::IfcGlobalId(), nullptr, std::string("Left Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcWallTypeEnum::IfcWallType_PARAPET);
   }
   else
   {
      left_railing = new typename Schema::IfcRailing(IfcParse::IfcGlobalId(), nullptr, std::string("Left Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_BALUSTRADE);
   }
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboLeft, left_railing, options, body_model_representation_subcontext);
   file.addEntity(left_railing);
   file.addRelatedObject<typename Schema::IfcRelContainedInSpatialStructure>(deck, left_railing);

   typename Schema::IfcProduct* right_railing;
   if (options.railings == CIfcExportOptions::Railings::Parapet)
   {
      right_railing = new typename Schema::IfcWall(IfcParse::IfcGlobalId(), nullptr, std::string("Right Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcWallTypeEnum::IfcWallType_PARAPET);
   }
   else
   {
      right_railing = new typename Schema::IfcRailing(IfcParse::IfcGlobalId(), nullptr, std::string("Right Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_BALUSTRADE);
   }
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboRight, right_railing, options, body_model_representation_subcontext);
   file.addEntity(right_railing);
   file.addRelatedObject<typename Schema::IfcRelContainedInSpatialStructure>(deck, right_railing);

   std::vector<typename Schema::IfcProduct*> railings{ left_railing,right_railing };
   if (options.classify)
   {
      Classify_Railings(file, railings);
      Create_Pset_usBridge_RailingCommon(file, railings);
   }

   if (options.include_work_plan)
   {
      // related the Stage 1 construction task to the pier product
      auto task = GetStage4Task(file, pBroker, options);
      AssignTaskToProduct(task, left_railing, file, pBroker, options);
      AssignTaskToProduct(task, right_railing, file, pBroker, options);
   }

   // Add girders to the spatial structure of the superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcElementAssembly::GIRDER

   // define the common properties of the precast girder type
   // specifically, precast concrete is defined by Pset_ConcreteElementGeneral with AssemblyPlace=FACTORY, CastingMethod=PRECAST
   typename Schema::IfcPropertySetDefinition::list::ptr property_sets(new typename Schema::IfcPropertySetDefinition::list);

   // Pset_ConcreteElementGeneral
   typename Schema::IfcProperty::list::ptr concrete_element_general_properties(new typename Schema::IfcProperty::list);
   // PEnum_AssemblyPlace
   std::vector<std::string> assembly_place_enum_values{ "FACTORY","OFFSITE","SITE","OTHER","UNKNOWN","UNSET" };
   auto assembly_place_property_enum_values = createPropertyEnumeration<Schema>("PEnum_AssemblyPlace", assembly_place_enum_values);
   auto assembly_place = createPropertyEnumeratedValue<Schema>("AssemblyPlace", assembly_place_property_enum_values, "FACTORY");
   concrete_element_general_properties->push(assembly_place);
   // PEnum_ConcreteCastingMethod
   std::vector<std::string> casting_method_enum_values{ "INSITU","MIXED","PRECAST","PRINTED","OTHER","UNKNOWN","UNSET" };
   auto casting_method_property_enum_values = createPropertyEnumeration<Schema>("PEnum_ConcreteCastingMethod", casting_method_enum_values);
   auto casting_method = createPropertyEnumeratedValue<Schema>("CastingMethod", casting_method_property_enum_values, "PRECAST");
   concrete_element_general_properties->push(casting_method);
   // create Pset_ConcreteElementGeneral
   auto pset_concrete_element_general = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ConcreteElementGeneral"), boost::none, concrete_element_general_properties);
   property_sets->push(pset_concrete_element_general);

   GET_IFACE2(pBroker, IDocumentType, pDocType);
   bool bIsPGSplice = pDocType->IsPGSpliceDocument();

   std::set<std::_tstring> beam_names;
   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   GroupIndexType nGroups = pIBridgeDesc->GetGirderGroupCount();
   for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
   {
      GirderIndexType nGirders = pIBridgeDesc->GetGirderGroup(grpIdx)->GetGirderCount();
      for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
      {
         auto pGirder = pIBridgeDesc->GetGirder(CGirderKey(grpIdx, gdrIdx));
         beam_names.insert(pGirder->GetGirderName());
      }
   }

   std::map<std::_tstring, typename Schema::IfcBeamType*> beam_types;
   if (bIsPGSplice)
   {
      auto beam_type = new typename Schema::IfcBeamType(
         IfcParse::IfcGlobalId(),
         nullptr, // OwnerHistory
         std::string("Precast Girder Type"), // Name
         boost::none, // Description
         std::string("IfcBeam/BEAM"), // ApplicableOccurrence 
         property_sets, // HasPropertySets (properties common to all beams of this type)
         boost::none, // RepresentationMaps (representations common to all beams of this type)
         boost::none, // Tag
         boost::none, // ElementType (type name if PredefinedType is USERDEFINED)
         Schema::IfcBeamTypeEnum::IfcBeamType_GIRDER_SEGMENT
      );
      file.addEntity(beam_type);
      beam_types.insert(std::make_pair(_T("Spliced_Girder_Type"), beam_type));
   }
   else
   {
      typename Schema::IfcObjectDefinition::list::ptr beam_object_definitions(new typename Schema::IfcObjectDefinition::list);
      for (auto beam_name : beam_names)
      {
         auto beam_type = new typename Schema::IfcBeamType(
            IfcParse::IfcGlobalId(),
            nullptr, // OwnerHistory
            std::string(T2A(beam_name.c_str())), // Name
            boost::none, // Description
            std::string("IfcBeam/BEAM"), // ApplicableOccurrence 
            property_sets, // HasPropertySets (properties common to all beams of this type)
            boost::none, // RepresentationMaps (representations common to all beams of this type)
            boost::none, // Tag
            boost::none, // ElementType (type name if PredefinedType is USERDEFINED)
            Schema::IfcBeamTypeEnum::IfcBeamType_BEAM
         );
         file.addEntity(beam_type);
         beam_object_definitions->push(beam_type);
         beam_types.insert(std::make_pair(beam_name, beam_type));
      }
      AssociateDocuments<Schema>(file, beam_object_definitions);
   }


   // build the beams
   typename Schema::IfcObject::list::ptr beam_objects(new typename Schema::IfcObject::list);


   std::vector<typename Schema::IfcProduct*> girders;
   for (GroupIndexType grpIdx = 0; grpIdx < nGroups; grpIdx++)
   {
      GirderIndexType nGirders = pIBridgeDesc->GetGirderGroup(grpIdx)->GetGirderCount();
      for (GirderIndexType gdrIdx = 0; gdrIdx < nGirders; gdrIdx++)
      {
         std::_tostringstream os;
         os << GIRDER_LABEL(CGirderKey(grpIdx, gdrIdx));
         std::string girder_name(T2A(os.str().c_str()));

         SegmentIndexType nSegments = pIBridgeDesc->GetGirderGroup(grpIdx)->GetGirder(gdrIdx)->GetSegmentCount();

         typename Schema::IfcElementAssembly* girder = nullptr;
         if (1 < nSegments)
         {
            girder = new typename Schema::IfcElementAssembly(IfcParse::IfcGlobalId(), nullptr, girder_name, boost::none, boost::none,
               nullptr/*ObjectPlacement - to be set in CreateGirderSegmentRepresentation*/,
               nullptr/*Representation - to be set in CreateGirderSegmentRepresentation*/,
               boost::none, Schema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_SITE, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
            file.addEntity(girder);
            file.addRelatedObject<typename Schema::IfcRelContainedInSpatialStructure>(superstructure, girder);
            girders.push_back(girder);
         }


         if (1 < nSegments && options.include_work_plan)
         {
            // related the Stage 1 construction task to the girder product
            auto task = GetStage1Task(file, pBroker, options);
            AssignTaskToProduct(task, girder, file, pBroker, options);
         }

         if (1 < nSegments)
         {
            typename Schema::IfcObjectDefinition::list::ptr list_of_girder_segments(new typename Schema::IfcObjectDefinition::list);
            for (SegmentIndexType segIdx = 0; segIdx < nSegments; segIdx++)
            {
               CSegmentKey segmentKey(grpIdx, gdrIdx, segIdx);
               std::ostringstream os_segment_name;
               os_segment_name << "Segment " << LABEL_SEGMENT(segIdx);
               auto segment_name = os_segment_name.str();

               auto beam_type = beam_types[_T("Spliced_Girder_Type")];
               auto beam = CreatePrecastSegment(file, pBroker, segment_name, segmentKey, beam_type, options);

               list_of_girder_segments->push(beam); // beams in this girder
               beam_objects->push(beam); // all beams

               if (segIdx < nSegments - 1)
               {
                  std::ostringstream os_closure_name;
                  os_closure_name << "Closure Joint " << LABEL_SEGMENT(segIdx);
                  auto closure_joint_name = os_closure_name.str();
                  auto closure_joint = new IfcSchema::IfcBeam(
                     IfcParse::IfcGlobalId(),
                     nullptr, // OwnerHistory
                     closure_joint_name,
                     std::string("Cast in place concrete closure joint"),
                     std::string("CLOSUREJOINT"),
                     nullptr,  // ObjectPlacement
                     nullptr,  // Representation
                     boost::none, // Tag
                     boost::none // PredefinedType (must not be used if defined in IfcBeamType)
                  );

                  CreateClosureJointRepresentation<Schema>(file, pBroker, segmentKey, closure_joint, options, body_model_representation_subcontext);

                  // need to do this for each bar, stirrup, etc
                  //typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_closure_joint_parts(new aggregate_of<typename Schema::IfcObjectDefinition>());

                  //auto longitudinal_rebar = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("longitudinal rebar"), boost::none,
                  //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade*/, boost::none/*nominal diameter*/,
                  //   boost::none /*nominal area*/, boost::none/*bar length*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_MAIN, Schema::IfcReinforcingBarSurfaceEnum::IfcReinforcingBarSurface_TEXTURED);
                  //file.addEntity(longitudinal_rebar);
                  //list_of_closure_joint_parts->push(longitudinal_rebar);

                  //auto stirrups = new typename Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, std::string("stirrups"), boost::none,
                  //   boost::none, nullptr, nullptr, boost::none, boost::none/*steel grade*/, boost::none/*nominal diameter*/,
                  //   boost::none /*nominal area*/, boost::none/*bar length*/, Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_SHEAR, Schema::IfcReinforcingBarSurfaceEnum::IfcReinforcingBarSurface_TEXTURED);
                  //file.addEntity(stirrups);
                  //list_of_closure_joint_parts->push(stirrups);

                  //auto closure_joint_parts_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("cast in place closure joint parts"), boost::none, closure_joint, list_of_closure_joint_parts);
                  //file.addEntity(closure_joint_parts_aggregates);


                  // Should do this in a IfcBeamType 
                  typename Schema::IfcObjectDefinition::list::ptr related_segments(new typename Schema::IfcObjectDefinition::list);
                  related_segments->push(closure_joint);

                  // Pset_ConcreteElementGeneral
                  typename Schema::IfcProperty::list::ptr concrete_element_general_properties(new typename Schema::IfcProperty::list);
                  // PEnum_AssemblyPlace
                  std::vector<std::string> assembly_place_enum_values{ "FACTORY","OFFSITE","SITE","OTHER","UNKNOWN","UNSET" };
                  auto assembly_place_property_enum_values = createPropertyEnumeration<Schema>("PEnum_AssemblyPlace", assembly_place_enum_values);
                  auto assembly_place = createPropertyEnumeratedValue<Schema>("AssemblyPlace", assembly_place_property_enum_values, "SITE");
                  concrete_element_general_properties->push(assembly_place);

                  std::vector<std::string> casting_method_enum_values{ "INSITU","MIXED","PRECAST","PRINTED","OTHER","UNKNOWN","UNSET" };
                  auto casting_method_property_enum_values = createPropertyEnumeration<Schema>("PEnum_ConcreteCastingMethod", casting_method_enum_values);
                  auto casting_method = createPropertyEnumeratedValue<Schema>("CastingMethod", casting_method_property_enum_values, "INSITU");
                  concrete_element_general_properties->push(casting_method);
                  auto pset_concrete_element_general = new typename Schema::IfcPropertySet(IfcParse::IfcGlobalId(), nullptr, std::string("Pset_ConcreteElementGeneral"), boost::none, concrete_element_general_properties);
                  file.addEntity(new typename Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, pset_concrete_element_general));

                  file.addEntity(closure_joint);
                  list_of_girder_segments->push(closure_joint);
               }
            } // next segment

            std::ostringstream os_relationship_name;
            os_relationship_name << "Elements of girder for Group " << LABEL_GROUP(grpIdx) << " Girder " << T2A(LABEL_GIRDER(gdrIdx));
            auto girder_aggregation_name = os_relationship_name.str();
            auto girder_aggregates = new typename Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, girder_aggregation_name, boost::none, girder, list_of_girder_segments);
            file.addEntity(girder_aggregates);
         }
         else
         {
            CSegmentKey segmentKey(grpIdx, gdrIdx, 0);
            std::_tostringstream os;
            os << GIRDER_LABEL(CGirderKey(grpIdx, gdrIdx));
            std::string girder_name(T2A(os.str().c_str()));

            auto beam_type = beam_types[pIBridgeDesc->GetGirder(segmentKey)->GetGirderName()];
            auto beam = CreatePrecastSegment(file, pBroker, girder_name, segmentKey, beam_type, options);

            beam_objects->push(beam); // all beams
            file.addRelatedObject<typename Schema::IfcRelContainedInSpatialStructure>(superstructure, beam);
            if (options.classify)
            {
               Create_Pset_usBridge_GirderCommon(file, pBroker, options, segmentKey, beam);
            }
         }
      } // next girder
   } // next group

   if (options.classify)
   {
      std::vector<typename Schema::IfcProduct*> beams;
      for (auto beam_object : *beam_objects)
      {
         beams.push_back(beam_object->as<typename Schema::IfcProduct>(beam_object));
      }
      Classify_PrecastGirderElements(file, beams);

      Classify_Girders(file, girders);
   }

   Create_Pset_BridgeCommon<Schema>(file,bridge);
   if (options.classify)
   {
      Create_Pset_usBridge_BridgeCommon<Schema>(file, pBroker, bridge);
      Classify_Bridge<Schema>(file, bridge);
   }
 }

  template <typename Schema>
 void InitializeFile(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, const CString& strFilePath)
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
    file.header().file_name()->setname(T2A(strFileName)); // filename without path
    //file.header().file_name().name(T2A(strFileName)); // filename without path
    std::vector<std::string> authors;
    authors.push_back(T2A(pProjectProperties->GetEngineer()));
    file.header().file_name()->setauthor(authors);
    //file.header().file_name().author(authors);
    std::vector<std::string> organizations;
    organizations.push_back(T2A(pProjectProperties->GetCompany()));
    file.header().file_name()->setorganization(organizations);
    //file.header().file_name().organization(organizations);
    //file.header().file_name().preprocessor_version(); // this is info about the toolkit we are using which is IfcOpenShell... this field is filled in by default

    std::vector<std::string> file_description;
    std::ostringstream os;
    os << "ViewDefinition [Alignment-basedView]" << std::ends;
    file_description.push_back(os.str().c_str());
    file.header().file_description()->setdescription(file_description);
    //file.header().file_description().description(file_description);

    std::_tostringstream _os;
    // https://github.com/buildingSMART/IFC4.x-IF/tree/header-policy/docs/IFC-file-header#originating_system
    _os << _T("Washington State Department of Transportation") << _T(" - ") << _T("BridgeLink:") << (pDocType->IsPGSuperDocument() ? _T("PGSuper") : _T("PGSplice")) << _T(" - ") << pVersionInfo->GetVersion(true).GetBuffer() << std::ends;
    std::string strVersion(T2A(_os.str().c_str()));
    file.header().file_name()->setoriginating_system(strVersion);
    //file.header().file_name().originating_system(strVersion);

    //auto project = file.addProject(); // Don't like the default units in IfcOpenShell so we have do build our own
    /////////////////////////// The following is copied from IfcHierarchyHelper<Schema>::addProject and tweaked
    typename Schema::IfcUnit::list::ptr units(new typename Schema::IfcUnit::list);

    auto* unit1 = new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
    auto* unit2 = new typename Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN);

    units->push(unit1);
    units->push(unit2);

    auto* unit_assignment = new typename Schema::IfcUnitAssignment(units);

    typename Schema::IfcRepresentationContext::list::ptr rep_contexts(new typename Schema::IfcRepresentationContext::list);
    auto* project = new typename Schema::IfcProject(IfcParse::IfcGlobalId(), owner_history, std::string("MyProject"), boost::none, boost::none, boost::none, boost::none, rep_contexts, unit_assignment);

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
    typename Schema::IfcActorRole::list::ptr roles(new typename Schema::IfcActorRole::list);
    auto role = new typename Schema::IfcActorRole(Schema::IfcRoleEnum::IfcRole_CIVILENGINEER, boost::none, boost::none);
    file.addEntity(role);
    roles->push(role);
    organization->setRoles(roles);

    owner_history->OwningApplication()->setApplicationDeveloper(organization);
 }

template <typename Schema>
bool CIfcExporter::BuildModel(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CString& strFilePath, const CIfcExportOptions& options)
{
   USES_CONVERSION;

   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Exporting IFC model"));

   IfcHierarchyHelper<Schema> file;
   InitializeFile<Schema>(file, pBroker, strFilePath); // creates project and site
   Create_Pset_ProjectCommon<Schema>(file);
   if (options.classify)
   {
      Create_Pset_usBridge_ProjectCommon<Schema>(file);
   }

   CreateAlignment<Schema>(file, pBroker, options); // creates alignment and aggregates with project, references into site spatial structure

   if (options.model_elements == CIfcExportOptions::ModelElements::AlignmentAndBridge)
   {
      CreateBridge<Schema>(file, pBroker, options); // creates bridge with site spatial structure
   }


   std::ofstream ofs(T2A(strFilePath));
   ofs << file;

   return true;
}


CIfcExporter::CIfcExporter(void)
{
}

CIfcExporter::~CIfcExporter(void)
{
}

bool CIfcExporter::BuildModel(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CString& strFilePath, const CIfcExportOptions& options)
{
   bool bResult = false;
   switch (options.schema)
   {
   case CIfcExportOptions::Schema::Schema_4x3_add2: bResult = BuildModel<IfcSchema>(pBroker, strFilePath, options); break;
   default:
      ATLASSERT(false); // is there a new typename Schema type
   }

   return bResult;
}
