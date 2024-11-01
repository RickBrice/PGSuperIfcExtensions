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
#include "IfcAlignmentBuilder.h"
#include "Referents.h"
#include "USBridge_Classifications.h"
#include "Units.h"

#include <IFace\VersionInfo.h>
#include <IFace\DocumentType.h>
#include <IFace\PrestressForce.h>

#include <EAF\EAFAutoProgress.h>
#include <PgsExt\PrecastSegmentData.h>
#include <WBFLGenericBridgeTools.h>
#include <PgsExt\BridgeDescription2.h>
#include <Plugins\BeamFamilyCLSID.h>


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
typename Schema::IfcCartesianPoint* ConvertPoint(IPoint2d* pPoint,bool bMirror = false)
{
   Float64 x, y;
   pPoint->Location(&x, &y);
   x = IsZero(x) ? 0.0 : x;
   y = IsZero(y) ? 0.0 : y;
   return new Schema::IfcCartesianPoint(std::vector<double>{bMirror ? -x : x, y});
}

template <typename Schema>
typename Schema::IfcCurve* CreatePolyline(IPoint2dCollection* polyPoints)
{
   typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());
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

   return new Schema::IfcPolyline(points);
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
      points->push(ConvertPoint<Schema>(point,true/*mirror about Y axis*/));
   }

   // we know that points is open (last point is not the same as the first)
   // the polygon must be closed by reference
   points->push(*(points->begin()));

   typename Schema::IfcCurve* curve = nullptr;
   if (options.sweep_profile == CIfcModelBuilderOptions::SweepProfile::IndexedPolyCurve)
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
   auto tendon_type = new Schema::IfcTendonType(
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

   file.addEntity(tendon_type);

   // add the new definition to the project
   if (rel_declares_instances->size() == 0)
   {
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_definitions(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      related_definitions->push(tendon_type);

      auto rel_declares = new Schema::IfcRelDeclares(
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
typename Schema::IfcReinforcingBarType* GetReinforcingBarType(IfcHierarchyHelper<Schema>& file, const WBFL::Materials::Rebar* pRebar)
{
   USES_CONVERSION;
   std::string name(T2A(pRebar->GetName().c_str()));

   // search to see if an IfcReinforcingBarType has already been created
   auto project = file.getSingle<typename Schema::IfcProject>();
   auto rel_declares_instances = file.instances_by_type<typename Schema::IfcRelDeclares>();
   for (auto& rel_declares : *rel_declares_instances)
   {
      if (rel_declares->RelatingContext()->as<typename Schema::IfcProject>())
      {
         auto related_definitions = rel_declares->RelatedDefinitions();
         for (auto& reldef : *related_definitions)
         {
            auto rebar_type = reldef->as<typename Schema::IfcReinforcingBarType>();
            if (rebar_type && rebar_type->Name() == name)
            {
               return rebar_type;
            }
         }
      }
   }

   // if we get this far, we need a new IfcReinforcingBarType
   auto rebar_type = new Schema::IfcReinforcingBarType(
      IfcParse::IfcGlobalId(),
      nullptr,
      name, /*Name*/
      boost::none, /*Description*/
      boost::none, /*ApplicableOccurrence*/
      boost::none, /*HasPropertySets*/
      boost::none, /*RepresentationMaps*/
      boost::none, /*Tag*/
      boost::none, /*ElementType*/
      Schema::IfcReinforcingBarTypeEnum::IfcReinforcingBarType_MAIN, /*PredefinedType*/
      pRebar->GetNominalDimension(), /*NominalDiameter*/
      pRebar->GetNominalArea(), /*CrossSectionArea*/
      boost::none, /*BarLength*/
      boost::none, /*BarSurface*/
      boost::none, /*BendingShapeCode*/
      boost::none /*BendingParameters*/
   );

   file.addEntity(rebar_type);

   // add the new definition to the project
   if (rel_declares_instances->size() == 0)
   {
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr related_definitions(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      related_definitions->push(rebar_type);

      auto rel_declares = new Schema::IfcRelDeclares(
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
            related_definitions->push(rebar_type);
            rel_declares->setRelatedDefinitions(related_definitions);
            break;
         }
      }
   }

   return rebar_type;
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

         auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
         ATLASSERT(geometric_representation_context);
         auto strand_shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), strand_representation_items);
         typename aggregate_of<typename Schema::IfcRepresentation>::ptr strand_shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
         strand_shape_representation_list->push(strand_shape_representation);
         auto strand_product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, strand_shape_representation_list);

         std::ostringstream os;
         os << strStrandType[strandType] << ":" << strandIdx + 1;
         auto strand = new Schema::IfcTendon(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, strand_placement, strand_product_definition_shape, boost::none, boost::none,
            boost::none, /*Schema::IfcTendonTypeEnum::IfcTendonType_STRAND,*/ // per 4.1.3.2, this must not be used unless PredefinedType at the ObjecType level is set to NOTDEFINED
            boost::none, /*pStrand->GetNominalDiameter() depreciated*/
            boost::none, /*pStrand->GetNominalArea() depreciated*/
            pStrandGeom->GetPjack(segmentKey, strandType),
            pStrandGeom->GetJackingStress(segmentKey, strandType),
            boost::none, boost::none, boost::none);
         file.addEntity(strand);

         auto* tendon_type = GetTendonType<Schema>(file,pStrand);

         if (tendon_type->Types()->size() == 0)
         {
            typename aggregate_of<typename Schema::IfcObject>::ptr related_objects(new aggregate_of<typename Schema::IfcObject>());
            related_objects->push(strand);

            auto rel_defines_by_type = new Schema::IfcRelDefinesByType(
               IfcParse::IfcGlobalId(),
               nullptr,
               std::string("strand defined by IfcTendonType"),
               boost::none,
               related_objects,
               tendon_type);

            file.addEntity(rel_defines_by_type);
         }
         else
         {
            auto rel_defines_set = tendon_type->Types();
            auto rel_defines = *(rel_defines_set->begin());
            auto rel_objects = rel_defines->RelatedObjects();
            rel_objects->push(strand);
            rel_defines->setRelatedObjects(rel_objects);
         }


         Classify_TPFPrestressing<Schema>(file, strand);

         Float64 db_start, db_end;
         bool bDebonded = pStrandGeom->IsStrandDebonded(segmentKey, strandIdx, strandType, nullptr, &db_start, &db_end);
         Create_Pset_TPFBridge_ReinforcementCommon(file, strand, bDebonded, db_start); // assumes symmetric debonding since classification can't handle unsymmetric

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

         strands->push(strand);
      }
   }

   return strands;
}

template <typename Schema>
typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr CreateRebars(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const pgsPointOfInterest& poiStart, const pgsPointOfInterest& poiEnd, typename Schema::IfcObjectPlacement* rebar_placement)
{
   USES_CONVERSION;

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr rebars(new aggregate_of<typename Schema::IfcObjectDefinition>());

   const CSegmentKey& segmentKey(poiStart.GetSegmentKey());

   GET_IFACE2(pBroker, IBridge, pBridge);
   Float64 Ls = pBridge->GetSegmentLength(segmentKey);
   Float64 Lg = pBridge->GetSegmentPlanLength(segmentKey);
   Float64 slope = pBridge->GetSegmentSlope(segmentKey);

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   GET_IFACE2(pBroker, ILongitudinalRebar, pLongRebar);
   const CLongitudinalRebarData* pLRD = pLongRebar->GetSegmentLongitudinalRebarData(segmentKey);

   GET_IFACE2(pBroker, ILongRebarGeometry, pLongRebarGeom);
   CComPtr<IRebarLayout> rebar_layout;
   pLongRebarGeom->GetRebarLayout(segmentKey, &rebar_layout);

   CComPtr<IEnumRebarLayoutItems> enum_items;
   rebar_layout->get__EnumRebarLayoutItems(&enum_items);

   IndexType layout_item_idx = 0;
   CComPtr<IRebarLayoutItem> rebar_layout_item;
   while (enum_items->Next(1, &rebar_layout_item, nullptr) != S_FALSE)
   {
      Float64 start, length;
      rebar_layout_item->get_Start(&start);
      rebar_layout_item->get_Length(&length);

      Float64 end = start + length;

      CComPtr<IEnumRebarPatterns> enum_patterns;
      rebar_layout_item->get__EnumRebarPatterns(&enum_patterns);
      CComPtr<IRebarPattern> rebar_pattern;
      while (enum_patterns->Next(1, &rebar_pattern, nullptr) != S_FALSE)
      {
         IndexType nBars;
         rebar_pattern->get_Count(&nBars);
         for (IndexType barIdx = 0; barIdx < nBars; barIdx++)
         {
            CComPtr<IPoint2d> p1, p2;
            rebar_pattern->get_Location(0.0, barIdx, &p1);
            rebar_pattern->get_Location(length, barIdx, &p2);

            typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());

            Float64 X, Y, Z;
            p1->Location(&Y, &Z);
            X = start * sqrt(1 + slope * slope);
            points->push(new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z}));

            p2->Location(&Y, &Z);
            X = end * sqrt(1 + slope * slope);
            points->push(new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z}));

            auto directrix = new Schema::IfcPolyline(points);
            file.addEntity(directrix);

            typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());

            CComPtr<IRebar> rb;
            rebar_pattern->get_Rebar(&rb);
            Float64 db;
            rb->get_NominalDiameter(&db);

            CComBSTR bar_name;
            rb->get_Name(&bar_name);
            WBFL::Materials::Rebar::Size bar_size = WBFL::LRFD::RebarPool::GetBarSize(OLE2CT(bar_name));

            auto swept_disk_solid = new Schema::IfcSweptDiskSolid(directrix, db / 2, boost::none, boost::none, boost::none);
            file.addEntity(swept_disk_solid);
            representation_items->push(swept_disk_solid);

            typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
            ATLASSERT(geometric_representation_context);
            auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
            shape_representation_list->push(shape_representation);
            auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

            const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(pLRD->BarType, pLRD->BarGrade, bar_size);
            std::ostringstream os;
            os << "Rebar Row " << (layout_item_idx+1) << " " << T2A(WBFL::LRFD::RebarPool::GetBarSize(bar_size).c_str());

            auto rebar = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, rebar_placement, product_definition_shape, boost::none,
               boost::none, // steel grade: depreciated
               boost::none, // nominal diameter: depreciated
               boost::none, // cross section area: depreciated
               boost::none, // bar length: depreciated
               boost::none, // predefined type: depreciated
               boost::none  // predefined type: depreciated
            );
            file.addEntity(rebar);

            auto* rebar_type = GetReinforcingBarType<Schema>(file, pRebar);

            if (rebar_type->Types()->size() == 0)
            {
               typename aggregate_of<typename Schema::IfcObject>::ptr related_objects(new aggregate_of<typename Schema::IfcObject>());
               related_objects->push(rebar);

               auto rel_defines_by_type = new Schema::IfcRelDefinesByType(
                  IfcParse::IfcGlobalId(),
                  nullptr,
                  std::string("rebar defined by IfcReinforcingBarType"),
                  boost::none,
                  related_objects,
                  rebar_type);

               file.addEntity(rel_defines_by_type);
            }
            else
            {
               auto rel_defines_set = rebar_type->Types();
               auto rel_defines = *(rel_defines_set->begin());
               auto rel_objects = rel_defines->RelatedObjects();
               rel_objects->push(rebar);
               rel_defines->setRelatedObjects(rel_objects);
            }
            rebars->push(rebar);
         }
         rebar_pattern.Release();
      }
      rebar_layout_item.Release();
      layout_item_idx++;
   }
   return rebars;
}

template <typename Schema>
typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr CreateStirrups(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const pgsPointOfInterest& poiStart, const pgsPointOfInterest& poiEnd, typename Schema::IfcObjectPlacement* segment_origin)
{
   USES_CONVERSION;

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr rebars(new aggregate_of<typename Schema::IfcObjectDefinition>());
   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   const auto& segmentKey = poiStart.GetSegmentKey();

   GET_IFACE2(pBroker, IMaterials, pMaterials);
   WBFL::Materials::Rebar::Type bar_type;
   WBFL::Materials::Rebar::Grade bar_grade;
   pMaterials->GetSegmentTransverseRebarMaterial(segmentKey, &bar_type, &bar_grade);

   WBFL::Materials::Rebar::Size bar_size = WBFL::Materials::Rebar::Size::bs4;

   const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(bar_type, bar_grade, bar_size);

   auto* rebar_type = GetReinforcingBarType<Schema>(file, pRebar);
   auto representation_maps = rebar_type->RepresentationMaps();
   if (!representation_maps || (*representation_maps)->size() == 0)
   {
      // WORKING HERE - The idea is to check to see if the beam is of the IBeam family, otherwise, don't model stirrups (Already doing this step in the calling function)
      // For Ibeams, start with WSDOT G2 bars, then change to G1 bars (but there are 2 bars, not 1)... then add the G3 bar in the top flange
      // This is just an experiment for how to model stirrups and a rebar cage.
      // When this is re-built as an extension agent, bar shape will be an input

      GET_IFACE2(pBroker, IGirder, pGirder);
      Float64 Hg = pGirder->GetHeight(poiStart);
      Float64 t = pGirder->GetWebThickness(poiStart, 0);

      GET_IFACE2(pBroker, IBridge, pBridge);
      Float64 A = pBridge->GetSlabOffset(segmentKey, pgsTypes::metStart);

      Float64 H1 = Hg + A + WBFL::Units::ConvertToSysUnits(3.0, WBFL::Units::Measure::Inch); // h1 = Hg + "A" + 3"

      Float64 cover = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);

      Float64 db = WBFL::Units::ConvertToSysUnits(0.5, WBFL::Units::Measure::Inch); // Assuming #4 bar, which is a dummy value

      Float64 dl = Hg - cover - db / 2;
      Float64 du = H1 - dl;
      Float64 dx = t / 2 - cover - db / 2;

      std::vector<std::vector<double>> point_list;
      point_list.push_back({ 0.0,-dx,du });
      point_list.push_back({ 0.0,-dx,-(dl - dx) });
      point_list.push_back({ 0.0,0.0,-dl });
      point_list.push_back({ 0.0,dx,-(dl - dx) });
      point_list.push_back({ 0.0,dx,du });

      typename aggregate_of<typename Schema::IfcSegmentIndexSelect>::ptr segments(new aggregate_of<typename Schema::IfcSegmentIndexSelect>());
      segments->push(new Schema::IfcLineIndex({ 1,2 }));
      segments->push(new Schema::IfcArcIndex({ 2,3,4 }));
      segments->push(new Schema::IfcLineIndex({ 4,5 }));

      auto directrix = new Schema::IfcIndexedPolyCurve(new Schema::IfcCartesianPointList3D(point_list, boost::none), segments, boost::none);
      file.addEntity(directrix);

      auto swept_disk_solid = new Schema::IfcSweptDiskSolid(directrix, db / 2, boost::none, boost::none, boost::none);
      file.addEntity(swept_disk_solid);

      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr rebar_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      rebar_representation_items->push(swept_disk_solid);

      auto rebar_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), rebar_representation_items);

      auto representation_map = new Schema::IfcRepresentationMap(file.addPlacement3d(), rebar_representation);

      if (!representation_maps) {
         typename aggregate_of<typename Schema::IfcRepresentationMap>::ptr rm(new aggregate_of<typename Schema::IfcRepresentationMap>());
         representation_maps = rm;
      }
      (*representation_maps)->push(representation_map);
      rebar_type->setRepresentationMaps(representation_maps);
   }

   representation_maps = rebar_type->RepresentationMaps();
   auto rebar_representation = (*((*representation_maps)->begin()))->MappedRepresentation();

   auto start = WBFL::Units::ConvertToSysUnits(1.5, WBFL::Units::Measure::Inch); // stirrups start 1.5" from face of beam
   Float64 offset = start;
   for (IndexType barIdx = 0; barIdx < 2; barIdx++, offset += 1.0)
   {
      typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());

      auto placement = segment_origin->as<typename Schema::IfcLocalPlacement>()->RelativePlacement()->as<typename Schema::IfcAxis2Placement3D>();
      auto mapping_source = new Schema::IfcRepresentationMap(placement, rebar_representation);

      auto mapping_target = new Schema::IfcCartesianTransformationOperator3D(nullptr, nullptr, new Schema::IfcCartesianPoint({ offset,0.,0. }), 1.0, nullptr);
      auto mapped_item = new Schema::IfcMappedItem(mapping_source, mapping_target);
      typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr mapped_representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
      mapped_representation_items->push(mapped_item);

      auto mapped_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("MappedRepresentation"), mapped_representation_items);

      shape_representation_list->push(mapped_representation);

      auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

      std::ostringstream os;
      os << "Dummy Stirrup " << barIdx;
      auto rebar = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, segment_origin, product_definition_shape, boost::none,
         boost::none, // steel grade: depreciated
         boost::none, // nominal diameter: depreciated
         boost::none, // cross section area: depreciated
         boost::none, // bar length: depreciated
         boost::none, // predefined type: depreciated
         boost::none  // predefined type: depreciated
      );
      file.addEntity(rebar);

      if (rebar_type->Types()->size() == 0)
      {
         typename aggregate_of<typename Schema::IfcObject>::ptr related_objects(new aggregate_of<typename Schema::IfcObject>());
         related_objects->push(rebar);

         auto rel_defines_by_type = new Schema::IfcRelDefinesByType(
            IfcParse::IfcGlobalId(),
            nullptr,
            std::string("rebar defined by IfcReinforcingBarType"),
            boost::none,
            related_objects,
            rebar_type);

         file.addEntity(rel_defines_by_type);
      }
      else
      {
         auto rel_defines_set = rebar_type->Types();
         auto rel_defines = *(rel_defines_set->begin());
         auto rel_objects = rel_defines->RelatedObjects();
         rel_objects->push(rebar);
         rel_defines->setRelatedObjects(rel_objects);
      }

      rebars->push(rebar);
   }
   return rebars;

   ///////////////////

   //const CSegmentKey& segmentKey(poiStart.GetSegmentKey());

   //GET_IFACE2(pBroker, IBridge, pBridge);
   //Float64 Ls = pBridge->GetSegmentLength(segmentKey);
   //Float64 Lg = pBridge->GetSegmentPlanLength(segmentKey);
   //Float64 slope = pBridge->GetSegmentSlope(segmentKey);

   //auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist

   //GET_IFACE2(pBroker, ILongitudinalRebar, pLongRebar);
   //const CLongitudinalRebarData* pLRD = pLongRebar->GetSegmentLongitudinalRebarData(segmentKey);

   //GET_IFACE2(pBroker, ILongRebarGeometry, pLongRebarGeom);
   //CComPtr<IRebarLayout> rebar_layout;
   //pLongRebarGeom->GetRebarLayout(segmentKey, &rebar_layout);

   //CComPtr<IEnumRebarLayoutItems> enum_items;
   //rebar_layout->get__EnumRebarLayoutItems(&enum_items);

   //IndexType layout_item_idx = 0;
   //CComPtr<IRebarLayoutItem> rebar_layout_item;
   //while (enum_items->Next(1, &rebar_layout_item, nullptr) != S_FALSE)
   //{
   //   Float64 start, length;
   //   rebar_layout_item->get_Start(&start);
   //   rebar_layout_item->get_Length(&length);

   //   Float64 end = start + length;

   //   CComPtr<IEnumRebarPatterns> enum_patterns;
   //   rebar_layout_item->get__EnumRebarPatterns(&enum_patterns);
   //   CComPtr<IRebarPattern> rebar_pattern;
   //   while (enum_patterns->Next(1, &rebar_pattern, nullptr) != S_FALSE)
   //   {
   //      IndexType nBars;
   //      rebar_pattern->get_Count(&nBars);
   //      for (IndexType barIdx = 0; barIdx < nBars; barIdx++)
   //      {
   //         CComPtr<IPoint2d> p1, p2;
   //         rebar_pattern->get_Location(0.0, barIdx, &p1);
   //         rebar_pattern->get_Location(length, barIdx, &p2);

   //         typename aggregate_of<typename Schema::IfcCartesianPoint>::ptr points(new aggregate_of<typename Schema::IfcCartesianPoint>());

   //         Float64 X, Y, Z;
   //         p1->Location(&Y, &Z);
   //         X = start * sqrt(1 + slope * slope);
   //         points->push(new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z}));

   //         p2->Location(&Y, &Z);
   //         X = end * sqrt(1 + slope * slope);
   //         points->push(new Schema::IfcCartesianPoint(std::vector<Float64>{X, Y, Z}));

   //         auto directrix = new Schema::IfcPolyline(points);
   //         file.addEntity(directrix);

   //         typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());

   //         CComPtr<IRebar> rb;
   //         rebar_pattern->get_Rebar(&rb);
   //         Float64 db;
   //         rb->get_NominalDiameter(&db);

   //         CComBSTR bar_name;
   //         rb->get_Name(&bar_name);
   //         WBFL::Materials::Rebar::Size bar_size = WBFL::LRFD::RebarPool::GetBarSize(OLE2CT(bar_name));

   //         auto swept_disk_solid = new Schema::IfcSweptDiskSolid(directrix, db / 2, boost::none, boost::none, boost::none);
   //         file.addEntity(swept_disk_solid);
   //         representation_items->push(swept_disk_solid);

   //         typename aggregate_of<typename Schema::IfcRepresentation>::ptr shape_representation_list(new aggregate_of<typename Schema::IfcRepresentation>());
   //         ATLASSERT(geometric_representation_context);
   //         auto shape_representation = new Schema::IfcShapeRepresentation(geometric_representation_context, std::string("Body"), std::string("AdvancedSweptSolid"), representation_items);
   //         shape_representation_list->push(shape_representation);
   //         auto product_definition_shape = new Schema::IfcProductDefinitionShape(boost::none, boost::none, shape_representation_list);

   //         const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(pLRD->BarType, pLRD->BarGrade, bar_size);
   //         std::ostringstream os;
   //         os << "Rebar Row " << (layout_item_idx + 1) << " " << T2A(WBFL::LRFD::RebarPool::GetBarSize(bar_size).c_str());

   //         auto rebar = new Schema::IfcReinforcingBar(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, rebar_placement, product_definition_shape, boost::none,
   //            boost::none, // steel grade: depreciated
   //            boost::none, // nominal diameter: depreciated
   //            boost::none, // cross section area: depreciated
   //            boost::none, // bar length: depreciated
   //            boost::none, // predefined type: depreciated
   //            boost::none  // predefined type: depreciated
   //         );
   //         file.addEntity(rebar);

   //         auto* rebar_type = GetReinforcingBarType<Schema>(file, pRebar);

   //         if (rebar_type->Types()->size() == 0)
   //         {
   //            typename aggregate_of<typename Schema::IfcObject>::ptr related_objects(new aggregate_of<typename Schema::IfcObject>());
   //            related_objects->push(rebar);

   //            auto rel_defines_by_type = new Schema::IfcRelDefinesByType(
   //               IfcParse::IfcGlobalId(),
   //               nullptr,
   //               std::string("rebar defined by IfcReinforcingBarType"),
   //               boost::none,
   //               related_objects,
   //               rebar_type);

   //            file.addEntity(rel_defines_by_type);
   //         }
   //         else
   //         {
   //            auto rel_defines_set = rebar_type->Types();
   //            auto rel_defines = *(rel_defines_set->begin());
   //            auto rel_objects = rel_defines->RelatedObjects();
   //            rel_objects->push(rebar);
   //            rel_defines->setRelatedObjects(rel_objects);
   //         }
   //         rebars->push(rebar);
   //      }
   //      rebar_pattern.Release();
   //   }
   //   rebar_layout_item.Release();
   //   layout_item_idx++;
   //}
   //return rebars;
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

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_RELEASED_SEGMENT | POI_5L, &vPoi);
   CHECK(vPoi.size() == 1);
   const pgsPointOfInterest& poiMS = vPoi.front();

   Float64 camber_ratio = 0.0;
   if (options.include_camber)
   {
      // Compute the camber ratio
      // https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_PrecastConcreteElementGeneral.htm
      // The camber deflection, measured from the midpoint of a cambered face of a piece to the midpoint of the chord joining the ends of the same face, 
      // as shown in the figure below (figure not provided), divided by the original (nominal) straight length of the face of the piece.
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

   GET_IFACE2(pBroker, IEAFDisplayUnits, pDisplayUnits);
   typename Schema::IfcConversionBasedUnit* stress_unit = nullptr;
   typename Schema::IfcConversionBasedUnit* displacement_unit = nullptr;

   Float64 fc = pMaterials->GetSegmentFc28(segmentKey);
   Float64 fci = pMaterials->GetSegmentFc(segmentKey, releaseIntervalIdx);
   Float64 fcl = pMaterials->GetSegmentFc(segmentKey, liftingIntervalIdx);
   Float64 fch = pMaterials->GetSegmentFc(segmentKey, haulingIntervalIdx);
   Float64 fpj = pStrandGeom->GetJackingStress(segmentKey, pgsTypes::Permanent);
   Float64 max_agg_size = pMaterials->GetSegmentMaxAggrSize(segmentKey);
   if (pDisplayUnits->GetUnitMode() == eafTypes::umUS)
   {
      stress_unit = GetStressUnit<Schema>(file,pBroker);
      displacement_unit = GetDisplacementUnit<Schema>(file, pBroker);

      fc = WBFL::Units::ConvertFromSysUnits(fc, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fci = WBFL::Units::ConvertFromSysUnits(fci, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fcl = WBFL::Units::ConvertFromSysUnits(fcl, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fch = WBFL::Units::ConvertFromSysUnits(fch, pDisplayUnits->GetStressUnit().UnitOfMeasure);
      fpj = WBFL::Units::ConvertFromSysUnits(fpj, pDisplayUnits->GetStressUnit().UnitOfMeasure);

      max_agg_size = WBFL::Units::ConvertFromSysUnits(max_agg_size, pDisplayUnits->GetDeflectionUnit().UnitOfMeasure);
   }

   // Pset_MaterialConcrete
   typename aggregate_of<typename Schema::IfcProperty>::ptr material_concrete_properties(new aggregate_of<typename Schema::IfcProperty>());
   material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("CompressiveStrength"), boost::none, new Schema::IfcPressureMeasure(fc), stress_unit));
   material_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("MaxAggregateSize"), boost::none, new Schema::IfcPositiveLengthMeasure(max_agg_size), displacement_unit));
   auto pset_material_concrete = new Schema::IfcMaterialProperties(std::string("Pset_MaterialConcrete"), boost::none/*description*/, material_concrete_properties, material);
   file.addEntity(pset_material_concrete);

   // Pset_PrecastConcreteElementGeneral
   typename aggregate_of<typename Schema::IfcProperty>::ptr precast_concrete_properties(new aggregate_of<typename Schema::IfcProperty>());
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("FormStrippingStrength"), boost::none, new Schema::IfcPressureMeasure(fci), stress_unit));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("LiftingStrength"), boost::none, new Schema::IfcPressureMeasure(fcl), stress_unit));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("ReleaseStrength"), boost::none, new Schema::IfcPressureMeasure(fci), stress_unit));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("TransportationStrength"), boost::none, new Schema::IfcPressureMeasure(fch), stress_unit));
   precast_concrete_properties->push(new Schema::IfcPropertySingleValue(std::string("InitialTension"), boost::none, new Schema::IfcPressureMeasure(fpj), stress_unit));
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

   if (options.include_quantities)
   {
      // Qto_BeamBaseQuantities
#pragma Reminder("NOTE: These are a little bit dummy quantities - updated in the future")
   // assuming simple sections (no change in cross section or depth like end blocks are variable depth hammerhead segments)
   // need to update pgsuper so we can get the different surface areas directly instead of having to compute them here
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

      if (pDisplayUnits->GetUnitMode() == eafTypes::umUS)
      {
         big_area_unit = GetBigAreaUnit<Schema>(file, pBroker);
         small_area_unit = GetSmallAreaUnit<Schema>(file, pBroker);
         volume_unit = GetVolumeUnit<Schema>(file, pBroker);
         mass_unit = GetMassUnit<Schema>(file, pBroker);
         length_unit = GetSpanLengthUnit<Schema>(file, pBroker);

         GSA = WBFL::Units::ConvertFromSysUnits(GSA, WBFL::Units::Measure::Feet2);
         GV = WBFL::Units::ConvertFromSysUnits(GV, WBFL::Units::Measure::Feet3);
         L = WBFL::Units::ConvertFromSysUnits(L, pDisplayUnits->GetSpanLengthUnit().UnitOfMeasure);
         A = WBFL::Units::ConvertFromSysUnits(A, pDisplayUnits->GetAreaUnit().UnitOfMeasure);
         Mass = WBFL::Units::ConvertFromSysUnits(Mass, WBFL::Units::Measure::PoundMass);
      }


      typename aggregate_of<typename Schema::IfcPhysicalQuantity>::ptr beam_quantities(new aggregate_of<typename Schema::IfcPhysicalQuantity>());
      beam_quantities->push(new Schema::IfcQuantityArea(std::string("GrossSurfaceArea"), boost::none, big_area_unit, GSA, boost::none));
      beam_quantities->push(new Schema::IfcQuantityVolume(std::string("GrossVolume"), boost::none, volume_unit, GV, boost::none));

      auto qto_bodygeometryvalidation = new Schema::IfcElementQuantity(IfcParse::IfcGlobalId(), nullptr, std::string("Qto_BodyGeometryValidation"), boost::none, boost::none, beam_quantities);
      file.addEntity(qto_bodygeometryvalidation);

      beam_quantities->push(new Schema::IfcQuantityLength(std::string("Length"), boost::none, length_unit, L, boost::none));
      beam_quantities->push(new Schema::IfcQuantityArea(std::string("CrossSectionArea"), boost::none, small_area_unit, A, boost::none));
      beam_quantities->push(new Schema::IfcQuantityArea(std::string("OuterSurfaceArea"), boost::none, big_area_unit, OSA, boost::none));
      beam_quantities->push(new Schema::IfcQuantityWeight(std::string("GrossWeight"), boost::none, mass_unit, Mass, boost::none));

      auto qto_beambasequantities = new Schema::IfcElementQuantity(IfcParse::IfcGlobalId(), nullptr, std::string("Qto_BeamBaseQuantities"), boost::none, boost::none, beam_quantities);
      file.addEntity(qto_beambasequantities);

      typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr related_segments(new aggregate_of<typename Schema::IfcObjectDefinition>());
      related_segments->push(segment);

      auto rel_defines_by_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, qto_bodygeometryvalidation);
      file.addEntity(rel_defines_by_properties);

      rel_defines_by_properties = new Schema::IfcRelDefinesByProperties(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, related_segments, qto_beambasequantities);
      file.addEntity(rel_defines_by_properties);
   }
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
void CreateLongitudinalRebarRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, typename Schema::IfcStyledRepresentation* styled_representation)
{
   USES_CONVERSION;

   // place rebar relative to the segment origin
   auto rebar_placement = file.addLocalPlacement(segment->ObjectPlacement());

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   const pgsPointOfInterest& poiStart(vPoi.front());
   const pgsPointOfInterest& poiEnd(vPoi.back());

   auto rebars = CreateRebars<Schema>(file, pBroker, poiStart, poiEnd, rebar_placement);

   if (0 < rebars->size())
   {
      // create the material
      auto rebar_material = new Schema::IfcMaterial("Reinforcement", boost::none/*description*/, boost::none/*category*/);
      file.addEntity(rebar_material);

      // assigns the presentation styles to the material
      typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
      list_of_representations->push(styled_representation);
      auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, rebar_material);
      file.addEntity(material_defintion_representation);

      GET_IFACE2(pBroker, IMaterials, pMaterials);
      WBFL::Materials::Rebar::Type rebar_type;
      WBFL::Materials::Rebar::Grade rebar_grade;
      pMaterials->GetSegmentLongitudinalRebarMaterial(segmentKey, &rebar_type, &rebar_grade);
      const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(rebar_type, rebar_grade, WBFL::Materials::Rebar::Size::bs3);
      auto fy = pRebar->GetYieldStrength();
      auto fpu = pRebar->GetUltimateStrength();
      auto eu = pRebar->GetElongation(); // depends on bar size and we are using a dummy #3 bar

      std::ostringstream os;
      os << T2A(pRebar->GetName().c_str());
      auto grade = os.str();

      // Pset_MaterialSteel
      typename aggregate_of<typename Schema::IfcProperty>::ptr material_steel_properties(new aggregate_of<typename Schema::IfcProperty>());
      //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("YieldStress"), boost::none, new Schema::IfcPressureMeasure(fy), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStress"), boost::none, new Schema::IfcPressureMeasure(fpu), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStrain"), boost::none, new Schema::IfcPositiveRatioMeasure(eu), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("StructuralGrade"), boost::none, new Schema::IfcLabel(grade.c_str()), nullptr));
      auto pset_material_steel = new Schema::IfcMaterialProperties(std::string("Pset_MaterialSteel"), boost::none/*description*/, material_steel_properties, rebar_material);
      file.addEntity(pset_material_steel);

      // need a list of entities that are associated with this material
      // right now we are creating a unique material for each strand but we still need the list
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr rebars_for_material(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      for (auto& rebar : *rebars)
      {
         rebars_for_material->push(rebar);
      }

      // associate the material with the segment (ie segments collection)
      auto rel_associates_materials = new Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Steel_to_Rebar"), boost::none, rebars_for_material, rebar_material);
      file.addEntity(rel_associates_materials);

      // aggregate the rebar with the segment
      auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Segment_Aggregates_Rebars"), boost::none, segment, rebars);
      file.addEntity(rel_aggregates);
   }
}


template <typename Schema>
void CreateStirrupRepresentation(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CSegmentKey& segmentKey, typename Schema::IfcBeam* segment, typename Schema::IfcStyledRepresentation* styled_representation)
{
   USES_CONVERSION;

   // For now, we only do stirrups for WF-Beams (that's because stirrups are dummy rebars)
   GET_IFACE2(pBroker, IBridgeDescription, pIBridgeDesc);
   const CBridgeDescription2* pBridgeDesc = pIBridgeDesc->GetBridgeDescription();
   const CGirderGroupData* pGroup = pBridgeDesc->GetGirderGroup(segmentKey.groupIndex);
   const GirderLibraryEntry* pGdrEntry = pGroup->GetGirderLibraryEntry(segmentKey.girderIndex);
   CComPtr<IBeamFactory> beam_factory;
   pGdrEntry->GetBeamFactory(&beam_factory);
   if (!::IsEqualGUID(beam_factory->GetFamilyCLSID(), CLSID_WFBeamFamily))
      return;


   // place rebar relative to the segment origin
   auto rebar_placement = file.addLocalPlacement(segment->ObjectPlacement());

   GET_IFACE2(pBroker, IPointOfInterest, pPoi);
   PoiList vPoi;
   pPoi->GetPointsOfInterest(segmentKey, POI_START_FACE | POI_END_FACE | POI_SECTCHANGE, &vPoi, POIFIND_OR);
   ATLASSERT(2 <= vPoi.size());

   const pgsPointOfInterest& poiStart(vPoi.front());
   const pgsPointOfInterest& poiEnd(vPoi.back());

   auto rebars = CreateStirrups<Schema>(file, pBroker, poiStart, poiEnd, rebar_placement);

   if (0 < rebars->size())
   {
      // create the material
      auto rebar_material = new Schema::IfcMaterial("Reinforcement", boost::none/*description*/, boost::none/*category*/);
      file.addEntity(rebar_material);

      // assigns the presentation styles to the material
      typename aggregate_of<typename Schema::IfcRepresentation>::ptr list_of_representations(new aggregate_of<typename Schema::IfcRepresentation>());
      list_of_representations->push(styled_representation);
      auto material_defintion_representation = new Schema::IfcMaterialDefinitionRepresentation(boost::none, boost::none, list_of_representations, rebar_material);
      file.addEntity(material_defintion_representation);

      GET_IFACE2(pBroker, IMaterials, pMaterials);
      WBFL::Materials::Rebar::Type rebar_type;
      WBFL::Materials::Rebar::Grade rebar_grade;
      pMaterials->GetSegmentLongitudinalRebarMaterial(segmentKey, &rebar_type, &rebar_grade);
      const auto* pRebar = WBFL::LRFD::RebarPool::GetInstance()->GetRebar(rebar_type, rebar_grade, WBFL::Materials::Rebar::Size::bs3);
      auto fy = pRebar->GetYieldStrength();
      auto fpu = pRebar->GetUltimateStrength();
      auto eu = pRebar->GetElongation(); // depends on bar size and we are using a dummy #3 bar

      std::ostringstream os;
      os << T2A(pRebar->GetName().c_str());
      auto grade = os.str();

      // Pset_MaterialSteel
      typename aggregate_of<typename Schema::IfcProperty>::ptr material_steel_properties(new aggregate_of<typename Schema::IfcProperty>());
      //https://ifc43-docs.standards.buildingsmart.org/IFC/RELEASE/IFC4x3/HTML/lexical/Pset_MaterialSteel.htm
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("YieldStress"), boost::none, new Schema::IfcPressureMeasure(fy), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStress"), boost::none, new Schema::IfcPressureMeasure(fpu), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("UltimateStrain"), boost::none, new Schema::IfcPositiveRatioMeasure(eu), nullptr));
      material_steel_properties->push(new Schema::IfcPropertySingleValue(std::string("StructuralGrade"), boost::none, new Schema::IfcLabel(grade.c_str()), nullptr));
      auto pset_material_steel = new Schema::IfcMaterialProperties(std::string("Pset_MaterialSteel"), boost::none/*description*/, material_steel_properties, rebar_material);
      file.addEntity(pset_material_steel);

      // need a list of entities that are associated with this material
      // right now we are creating a unique material for each strand but we still need the list
      typename aggregate_of<typename Schema::IfcDefinitionSelect>::ptr rebars_for_material(new aggregate_of<typename Schema::IfcDefinitionSelect>());
      for (auto& rebar : *rebars)
      {
         rebars_for_material->push(rebar);
      }

      // associate the material with the segment (ie segments collection)
      auto rel_associates_materials = new Schema::IfcRelAssociatesMaterial(IfcParse::IfcGlobalId(), nullptr, std::string("Associates_Steel_to_Rebar"), boost::none, rebars_for_material, rebar_material);
      file.addEntity(rel_associates_materials);

      // aggregate the rebar with the segment
      auto rel_aggregates = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Segment_Aggregates_Rebars"), boost::none, segment, rebars);
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
         if (polyline)
         {
            std::ostringstream os;
            os << (tbOrientation == pgsTypes::tboLeft ? "Left" : "Right") << " Barrier";
            if (1 < nShapesPerBarrier)
               os << " Shape " << shapeIdx;

            auto shape_perimeter = new Schema::IfcArbitraryClosedProfileDef(Schema::IfcProfileTypeEnum::IfcProfileType_AREA, os.str(), polyline);
            cross_sections[shapeIdx]->push(shape_perimeter);
         }
      } // next shape
   } // next section

   typename aggregate_of<typename Schema::IfcRepresentationItem>::ptr representation_items(new aggregate_of<typename Schema::IfcRepresentationItem>());
   for (IndexType shapeIdx = 0; shapeIdx < nShapesPerBarrier; shapeIdx++)
   {
      if (0 < cross_sections[shapeIdx]->size())
      {
         auto sectioned_solid = new Schema::IfcSectionedSolidHorizontal(directrix, cross_sections[shapeIdx], cross_section_positions);
         representation_items->push(sectioned_solid);
         file.addEntity(sectioned_solid);
      }
   }

   if (0 < representation_items->size())
   {
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
}

template <typename Schema>
typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr CreatePiers(IfcHierarchyHelper<Schema>& file, IBroker* pBroker, const CIfcModelBuilderOptions& options)
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

      if (options.classify)
      {
         if (pBridge->IsAbutment(pierIdx))
            Classify_TPFAbutment<Schema>(file, pier);
         else
            Classify_TPFPier<Schema>(file, pier);
      }

      std::ostringstream os;
      os << "Foundation at " << pier_name; // name is not required, but is specified in AASHTO IDS
      auto foundation = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, os.str(), boost::none, boost::none, nullptr, nullptr, boost::none,
         Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
         Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
         Schema::IfcBridgePartTypeEnum::IfcBridgePartType_FOUNDATION);
      file.addEntity(foundation);
      if (options.classify)
      {
         Classify_TPFFoundation<Schema>(file, foundation);
      }

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

   if (options.classify)
   {
      Add_TPF_Classification<Schema>(file);
   }


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
   if (options.classify)
   {
      Classify_TPFSuperstructure<Schema>(file, superstructure);
   }

   auto substructure = new Schema::IfcBridgePart(IfcParse::IfcGlobalId(), nullptr, std::string("Substructure"), boost::none, boost::none, nullptr, nullptr, boost::none,
      Schema::IfcElementCompositionEnum::IfcElementComposition_PARTIAL,
      Schema::IfcFacilityUsageEnum::IfcFacilityUsage_LONGITUDINAL,
      Schema::IfcBridgePartTypeEnum::IfcBridgePartType_SUBSTRUCTURE);
   file.addEntity(substructure);
   if (options.classify)
   {
      Classify_TPFSubstructure<Schema>(file, substructure);
   }

   typename aggregate_of<typename Schema::IfcObjectDefinition>::ptr list_of_bridge_parts(new aggregate_of<typename Schema::IfcObjectDefinition>());
   list_of_bridge_parts->push(superstructure);
   list_of_bridge_parts->push(substructure);

   // IfcBridge <-> IfcRelAggregates <-> IfcBridgePart::SUPERSTRUCTURE, SUBSTRUCTURE
   auto bridge_spatial_elements = new Schema::IfcRelAggregates(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in bridge spatial structure"), boost::none, bridge, list_of_bridge_parts);
   file.addEntity(bridge_spatial_elements);



   // Create spatial structure of substructure
   // IfcBridgePart::SUBSTRUCTURE <-> IfcRelAggregates <-> IfcBridgePart::ABUTMENT, PIER
   auto list_of_piers = CreatePiers(file, pBroker, options);
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
#pragma Reminder("PGSUPER needs the concept of No Railing - here we define a railing product with no representation which says there is a railing")
   typename Schema::IfcProduct* left_railing;
   if (options.railings == CIfcModelBuilderOptions::Railings::Parapet)
   {
      left_railing = new Schema::IfcWall(IfcParse::IfcGlobalId(), nullptr, std::string("Left Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcWallTypeEnum::IfcWallType_PARAPET);
   }
   else
   {
      left_railing = new Schema::IfcRailing(IfcParse::IfcGlobalId(), nullptr, std::string("Left Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_BALUSTRADE);
   }
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboLeft, left_railing, options, body_model_representation_subcontext);
   file.addEntity(left_railing);
   list_of_superstructure_elements->push(left_railing);

   typename Schema::IfcProduct* right_railing;
   if (options.railings == CIfcModelBuilderOptions::Railings::Parapet)
   {
      right_railing = new Schema::IfcWall(IfcParse::IfcGlobalId(), nullptr, std::string("Right Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcWallTypeEnum::IfcWallType_PARAPET);
   }
   else
   {
      right_railing = new Schema::IfcRailing(IfcParse::IfcGlobalId(), nullptr, std::string("Right Railing"), boost::none, boost::none, nullptr, nullptr, boost::none, Schema::IfcRailingTypeEnum::IfcRailingType_BALUSTRADE);
   }
   CreateRailingSystemRepresentation(file, pBroker, pgsTypes::tboRight, right_railing, options, body_model_representation_subcontext);
   file.addEntity(right_railing);
   list_of_superstructure_elements->push(right_railing);

   std::vector<typename Schema::IfcProduct*> railings{ left_railing,right_railing };
   if (options.classify)
   {
      Classify_TPFRailings(file, railings);
      Create_Pset_TPFBridge_RailingCommon(file, railings);
   }

   // Add girders to the spatial structure of the superstructure
   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcElementAssembly::GIRDER

   // first create a representation object for the girder material. this is one way to add presentation information such as color.
   // this helper function just sets color and hard codes all the other parameters - this can be expanded in the future
   auto girder_material_representation = CreateMaterialRepresentation<Schema>("Girder", 7.6078431372549E-1, 7.72549019607843E-1, 8.E-1, geometric_representation_context);
   auto strand_material_representation = CreateMaterialRepresentation<Schema>("Strand", 1, 0, 0, geometric_representation_context);
   auto rebar_material_representation = CreateMaterialRepresentation<Schema>("Rebar", 0, 1, 0, geometric_representation_context);

   std::vector<typename Schema::IfcProduct*> girders;
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
            boost::none, Schema::IfcAssemblyPlaceEnum::IfcAssemblyPlace_FACTORY, Schema::IfcElementAssemblyTypeEnum::IfcElementAssemblyType_GIRDER);
         file.addEntity(girder);
         list_of_superstructure_elements->push(girder);
         girders.push_back(girder);

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

            CreateLongitudinalRebarRepresentation<Schema>(file, pBroker, segmentKey, segment, rebar_material_representation);

            CreateStirrupRepresentation<Schema>(file, pBroker, segmentKey, segment, rebar_material_representation);

            file.addEntity(segment);
            list_of_girder_segments->push(segment);

            if (nSegments == 1 && options.classify)
            {
               // TPFBridge_GirderCommon doesn't work for spliced girders so only do this for precast
               Create_Pset_TPFBridge_GirderCommon(file, pBroker, options, segmentKey, girder);
            }

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

   if (options.classify)
   {
      Classify_TPFGirders(file, girders);
   }

   // IfcBridgePart::SUPERSTRUCTURE <-> IfcRelContainedInSpatialStructure <-> IfcRailing, IfcElementAssembly::GIRDER
   auto rel_contained_in_superstructure_spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), nullptr, std::string("Elements in superstructure spatial structure"), boost::none, list_of_superstructure_elements, superstructure);
   file.addEntity(rel_contained_in_superstructure_spatial_structure);

   Create_Pset_BridgeCommon<Schema>(file,bridge);
   if (options.classify)
   {
      Create_Pset_TPFBridge_BridgeCommon<Schema>(file, pBroker, bridge);
      Classify_TPFBridge<Schema>(file, bridge);
   }
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
   if (options.classify)
   {
      Create_Pset_TPFBridge_ProjectCommon<Schema>(file);
   }

   CreateAlignment<Schema>(file, pBroker, options); // creates alignment and aggregates with project, references into site spatial structure

   if (options.model_elements == CIfcModelBuilderOptions::ModelElements::AlignmentAndBridge)
   {
      CreateBridge<Schema>(file, pBroker, options); // creates bridge with site spatial structure
   }

   CreateReferents(file, pBroker, options);


   std::ofstream ofs(T2A(strFilePath));
   ofs << file;

   return true;
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
