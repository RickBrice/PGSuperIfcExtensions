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
#include "DeckSlab.h"
#include "Geometry.h"

#include <map>

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/ifcgeomelement.h>
#include <ifcgeom/kernels/opencascade/OpenCascadeKernel.h>

#include <IFace\Project.h>
#include <IFace/Alignment.h>
#include <EAF/EAFProgress.h>
#include <EAF/AutoProgress.h>


#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file


// returns the id of the IfcSlab... assumes there is only one
int get_slab_id(IfcParse::IfcFile& file)
{
   auto slabs = file.instances_by_type<IfcSchema::IfcSlab>();
   auto it = std::find_if(slabs->begin(), slabs->end(), [](const auto& slab) {return slab->PredefinedType() == IfcSchema::IfcSlabTypeEnum::IfcSlabType_FLOOR; });
   if (it == slabs->end())
      return 0; // no slabs

   auto slab = *it;
   return slab->id();
}

template <typename T,typename U>
std::pair<typename std::multimap<T,U>::const_iterator,
          typename std::multimap<T,U>::const_iterator>
   equal_range_tolerance(const std::multimap<T, U>& mm, T key, T tol)
{
   auto lower = mm.lower_bound(key - tol);
   auto upper = mm.upper_bound(key + tol);

   return { lower,upper };
}

std::map<double, std::pair<double, double>> condense(const std::multimap<double, double>& mm)
{
   std::map<double, std::pair<double, double>> result;
   for (auto it = mm.begin(); it != mm.end();)
   {
      double key = it->first;

      // range of equal keys
      //auto range = mm.equal_range(key);
      auto range = equal_range_tolerance(mm, key, 1e-6);

      // initialize min/max
      double minv = range.first->second;
      double maxv = range.first->second;

      int count = 0;
      for (auto jt = range.first; jt != range.second; ++jt) {
         minv = std::min(minv, jt->second);
         maxv = std::max(maxv, jt->second);
         ++count;
      }

      if (count == 1) {
         result[key] = { minv,std::numeric_limits<double>::infinity() };
      }
      else {
         result[key] = { minv,maxv };
      }

      it = range.second; // move to the next key group
   }

   // remove all entries where second == std::numeric_limits<double>::infinity()
   for (auto it = result.begin(); it != result.end(); ) {
      if (it->second.second == std::numeric_limits<double>::infinity())
         it = result.erase(it);
      else
         ++it;
   }

   return result;
}

std::map<double, std::pair<double, double>> get_deck_slab(std::shared_ptr<WBFL::EAF::Broker> pBroker, IfcParse::IfcFile& file)
{
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Progressing deck geometry"));

   auto slab_id = get_slab_id(file);
   ASSERT(slab_id != 0); // should not be calling into this function if there isn't a deck slab

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   std::multimap<double, double> deck_points;
   CComPtr<IPoint2d> point;
   point.CoCreateInstance(CLSID_Point2d);

   ifcopenshell::geometry::Settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   // set up filter for the geometry iterator
   IfcGeom::instance_id_filter filter(true, false, { slab_id });
   std::vector<IfcGeom::filter_t> filters({ filter });

   std::unique_ptr<IfcGeom::OpenCascadeKernel> kernel(std::make_unique<IfcGeom::OpenCascadeKernel>(settings));
   IfcGeom::Iterator iterator(std::move(kernel), settings, &file, filters,1);
   bool bResult = iterator.initialize();
   do
   {
      auto element = iterator.get();

      auto triangulation = dynamic_cast<IfcGeom::TriangulationElement*>(element);
      auto geometry = triangulation->geometry_pointer();
      const auto& verts = geometry->verts();
      const auto& faces = geometry->faces();

      std::vector<Eigen::Vector3d> v;
      v.reserve(verts.size() / 3);
      for (auto i = 0; i < verts.size(); i += 3)
      {
         v.emplace_back(verts[i], verts[i + 1], verts[i + 2]);
      }

      std::vector<std::array<int, 3>> f;
      f.reserve(faces.size() / 3);
      for (auto i = 0; i < faces.size(); i += 3)
      {
         // This line of code should be correct, however....
         f.emplace_back(std::array<int, 3>({ faces[i],faces[i + 1],faces[i + 2] }));

         // There appears to be a bug in the IfcOpenShell processing of IfcSectionedSolidHorizontal
         // The vertex indices of the triangular faces are in the wrong direction resulting in an
         // inward surface normal. Working around this problem by manually reversing the order
         //f.emplace_back(std::array<int, 3>({ faces[i+2],faces[i + 1],faces[i] }));
      }

      // 1. Build mesh
      Mesh m(v, f);

      // 2. Decompose mesh into smooth components
      std::vector<Mesh> decomposed = m.decompose();

      // 3. Select top component (max Z of vertex + avg normal)
      Mesh top_component = get_top_mesh(decomposed);

      // 4. Extract boundary wire
      Wire boundary = top_component.boundary();

      // 5. Decompose boundary into smooth segments
      std::vector<Wire> segments = boundary.decompose();

      // 6. Sort segments by length
      std::sort(segments.begin(), segments.end(),
         [](const Wire& a, const Wire& b)
         {
            return a.plan_length() < b.plan_length();
         });

      // Assume two longest segments are the slab edges
      for (auto iter = segments.end() - 2; iter != segments.end(); iter++)
      {
         Wire side = *iter;

         // first point of first edge
         auto v1 = side.verts[side.edges[0][0]];
         point->Move(v1.x(), v1.y());
         Float64 station, offset;
         pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &station, &offset);
         deck_points.insert(std::make_pair(station, offset));

         // second point of all edges (assuming end of prev = start of next)
         for (auto& edge : side.edges)
         {
            auto v2 = side.verts[edge[1]];
            point->Move(v2.x(), v2.y());
            Float64 station, offset;
            pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &station, &offset);
            deck_points.insert(std::make_pair(station, offset));
         }
      }
   } while (iterator.next());

   auto results = condense(deck_points);

   return results;
}

bool create_alignment_from_deck(std::shared_ptr<WBFL::EAF::Broker> pBroker, IfcParse::IfcFile& file)
{
   // The bridge has a deck, but the model does not have an alignment.
   // We can't compute station and offset of the deck edge points without an alignment.
   // Assume the alignment is straight. Create the alignment using the center point of the deck ends to establish the start point and direction.
   // This is a best effort approach to get an alignment that is usable for stationing the deck geometry. 
   // It will not be correct, but it will at least allow the user to see the deck geometry in the correct location and station the deck geometry.

   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Progressing deck geometry to establish an alignment"));

   auto slab_id = get_slab_id(file);
   ASSERT(slab_id != 0); // should not be calling into this function if there isn't a deck slab

   ifcopenshell::geometry::Settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   // set up filter for the geometry iterator
   IfcGeom::instance_id_filter filter(true, false, { slab_id });
   std::vector<IfcGeom::filter_t> filters({ filter });

   std::unique_ptr<IfcGeom::OpenCascadeKernel> kernel(std::make_unique<IfcGeom::OpenCascadeKernel>(settings));
   IfcGeom::Iterator iterator(std::move(kernel), settings, &file, filters, 1);
   bool bResult = iterator.initialize();
   do
   {
      auto element = iterator.get();

      auto triangulation = dynamic_cast<IfcGeom::TriangulationElement*>(element);
      auto geometry = triangulation->geometry_pointer();
      const auto& verts = geometry->verts();
      const auto& faces = geometry->faces();

      std::vector<Eigen::Vector3d> v;
      v.reserve(verts.size() / 3);
      for (auto i = 0; i < verts.size(); i += 3)
      {
         v.emplace_back(verts[i], verts[i + 1], verts[i + 2]);
      }

      std::vector<std::array<int, 3>> f;
      f.reserve(faces.size() / 3);
      for (auto i = 0; i < faces.size(); i += 3)
      {
         // This line of code should be correct, however....
         f.emplace_back(std::array<int, 3>({ faces[i],faces[i + 1],faces[i + 2] }));

         // There appears to be a bug in the IfcOpenShell processing of IfcSectionedSolidHorizontal
         // The vertex indices of the triangular faces are in the wrong direction resulting in an
         // inward surface normal. Working around this problem by manually reversing the order
         //f.emplace_back(std::array<int, 3>({ faces[i+2],faces[i + 1],faces[i] }));
      }

      // 1. Build mesh
      Mesh m(v, f);

      // 2. Decompose mesh into smooth components
      std::vector<Mesh> decomposed = m.decompose();

      // 3. Select top component (max Z of vertex + avg normal)
      Mesh top_component = get_top_mesh(decomposed);

      // 4. Extract boundary wire
      Wire boundary = top_component.boundary();

      // 5. Decompose boundary into smooth segments
      std::vector<Wire> segments = boundary.decompose();

      // 6. Sort segments by length
      std::sort(segments.begin(), segments.end(),
         [](const Wire& a, const Wire& b)
         {
            return a.plan_length() < b.plan_length();
         });

      // 7. Assume two shortest segments are the ends of the slab
      Wire start_seg = segments[0];
      Wire end_seg = segments[1];

      // 8. Get center point of edge wire bounding-box
      Eigen::Vector3d c1 = start_seg.bbox_center();
      Eigen::Vector3d c2 = end_seg.bbox_center();

      // 9. Create alignment between those two points
      GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
      AlignmentData2 alignmentData = pRoadwayData->GetAlignmentData2();
      alignmentData.xRefPoint = c1.x();
      alignmentData.yRefPoint = c1.y();
      alignmentData.Direction = atan2(c2.y() - c1.y(), c2.x() - c1.x());
      pRoadwayData->SetAlignmentData2(alignmentData);

      return true;

   } while (iterator.next());

   return false;
}