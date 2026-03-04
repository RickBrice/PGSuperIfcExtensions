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
#include "BeamSpacing.h"
#include "Geometry.h"
#include "Utilities.h"
#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/ifcgeomelement.h>
#include <ifcgeom/kernels/opencascade/OpenCascadeKernel.h>

// gets all the IfcBeams
// This needs to be updated so we get only the superstructure beams
std::set<int> get_beam_ids(IfcParse::IfcFile& file)
{
   auto superstructure = GetBridgePart(file, IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   // lambda function to filter beams that are contained in the superstructure
   auto filter = [&superstructure](auto beam) {
         auto related_elements = beam->ContainedInStructure();
         if (related_elements)
         {
            for (auto related_element : *related_elements)
            {
               if (related_element->RelatingStructure() == superstructure)
                  return true;
            }
            return false;
         }
         else
         {
            return false; // not in a spatial structure
         }
      };


   auto beams = file.instances_by_type<IfcSchema::IfcBeam>();

   std::set<int> beam_ids;
   for (int id : *beams
      //| std::views::filter([&os](auto beam) {return beam->Name() && beam->Name()->starts_with(os.str()); }) // filter all beams that start with "Span n"
      | std::views::filter(filter)
      | std::views::transform([](auto beam) {return beam->id(); })) // transform the beam to its id
   {
      beam_ids.insert(id);
   }

   return beam_ids;
}

// this function attempts to get the beam spacing based on the geometry (not using a property from a pset)
std::pair<Spacing, Spacing> get_beam_spacing(IfcParse::IfcFile& file)
{
   WBFL::System::Logger::Debug(_T("Beam Spacing"));

   auto beam_ids = get_beam_ids(file);

   ifcopenshell::geometry::Settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   // set up filter for the geometry iterator
   IfcGeom::instance_id_filter filter(true, false, beam_ids);
   std::vector<IfcGeom::filter_t> filters({ filter });

   std::unique_ptr<IfcGeom::OpenCascadeKernel> kernel(std::make_unique<IfcGeom::OpenCascadeKernel>(settings));
   int num_threads = std::thread::hardware_concurrency();
   IfcGeom::Iterator iterator(std::move(kernel), settings, &file, filters, num_threads);

   std::map<GroupIndexType, std::vector<Eigen::Vector3d>> start_points;
   std::map<GroupIndexType, std::vector<Eigen::Vector3d>> end_points;
   bool bResult = iterator.initialize();
   // NOTE: need to deal with (bResult == false)

   // This do loop can be multi-threaded
   do
   {
      auto element = iterator.get();
      //std::cout << element->name() << std::endl;;
      auto beam = element->product()->as<IfcSchema::IfcBeam>();
      auto girder_key = get_girder_key(beam);
      if (girder_key == CGirderKey())
      {
         continue;
      }

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
         f.emplace_back(std::array<int, 3>({ faces[i],faces[i + 1],faces[i + 2] }));
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
      std::vector<Wire> segments = boundary.decompose(1.0/*one degree difference in direction*/);

      // 6. Sort segments by length
      std::sort(segments.begin(), segments.end(),
         [](const Wire& a, const Wire& b)
         {
            return a.length() < b.length();
         });

      // 7. Take the two shortest segments
      Wire start_seg = segments[0];
      Wire end_seg = segments[1];

      // 8. Get center point of edge wire bounding-box
      Eigen::Vector3d c1 = start_seg.bbox_center();
      Eigen::Vector3d c2 = end_seg.bbox_center();

      // since the short length of each end is essentially the same,
      // the start and end can be out of order. this is a bad way to test.
      // closest to origin is assumed to be at the start end
      // for beams going NE to SW from the NE quadrant, this would not be true
      if (c2.norm() < c1.norm()) std::swap(c1, c2);

      if (start_points[girder_key.groupIndex].size() <= girder_key.girderIndex)
         start_points[girder_key.groupIndex].resize(girder_key.girderIndex+1);

      if (end_points[girder_key.groupIndex].size() <= girder_key.girderIndex)
         end_points[girder_key.groupIndex].resize(girder_key.girderIndex+1);

      start_points[girder_key.groupIndex][girder_key.girderIndex] = c1;
      end_points[girder_key.groupIndex][girder_key.girderIndex] = c2;

      std::ostringstream os;
      os << "Group " << girder_key.groupIndex << ", " << "Girder " << girder_key.girderIndex << " " << c1.transpose() << " -> " << c2.transpose();
      WBFL::System::Logger::Debug(os.str().c_str());
   } while (iterator.next());

   Spacing ss;
   for (auto [grpIdx, points] : start_points)
   {
      ss[grpIdx].reserve(points.size() - 1);
      std::ranges::transform(std::views::iota(size_t{ 1 }, points.size()),
         std::back_inserter(ss[grpIdx]),
         [&](size_t i) {return (points[i] - points[i - 1]).norm(); });
   }

   Spacing es;
   for (auto [grpIdx, points] : end_points)
   {
      es[grpIdx].reserve(points.size() - 1);
      std::ranges::transform(std::views::iota(size_t{ 1 }, points.size()),
         std::back_inserter(es[grpIdx]),
         [&](size_t i) {return (points[i] - points[i - 1]).norm(); });
   }

   //std::cout << "Start Spacing, ";
   //for (auto s : ss) std::cout << s / 0.3048 << ", ";
   //std::cout << std::endl;

   return { ss, es };
}
