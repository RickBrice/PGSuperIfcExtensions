///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2026  Washington State Department of Transportation
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
#include "Piers.h"

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/element.h>
#include <ifcgeom/kernels/opencascade/opencascade_kernel.h>

#include <EAF/EAFProgress.h>
#include <EAF/AutoProgress.h>
#include <psgLib/GirderLabel.h>
#include <IFace/Alignment.h>

// gets all the IfcBeams
// This needs to be updated so we get only the superstructure beams
std::set<int> get_beam_ids(ifcopenshell::file& file)
{
   auto superstructure = GetBridgePart(file, IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_SUPERSTRUCTURE);

   // lambda function to filter beams that are contained in the superstructure
   auto filter = [&superstructure](auto beam) {
         auto related_elements = beam.ContainedInStructure();
         for (auto& related_element : related_elements)
         {
            if (related_element.RelatingStructure() == superstructure && GetPredefinedType<IfcSchema::IfcBeam, IfcSchema::IfcBeamType, IfcSchema::IfcBeamTypeEnum::Value>(beam) == IfcSchema::IfcBeamTypeEnum::IfcBeamType_BEAM)
               return true;
         }
         return false; // not in the superstructure's spatial structure
      };


   auto beams = file.instances_by_type<IfcSchema::IfcBeam>();

   std::set<int> beam_ids;
#if _HAS_CXX20
   for (int id : beams
      //| std::views::filter([&os](auto beam) {return beam->Name() && beam->Name()->starts_with(os.str()); }) // filter all beams that start with "Span n"
      | std::views::filter(filter)
      | std::views::transform([](auto beam) {return beam.id(); })) // transform the beam to its id
   {
      beam_ids.insert(id);
   }
#else
   for (auto beam : beams)
   {
      if (filter(beam))
      {
         beam_ids.insert(beam.id());
      }
   }
#endif

   return beam_ids;
}

// Locates the girders from the IfcBeam geometry
GirderLayout get_girder_layout(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const std::vector<Float64>& pier_stations)
{
   USES_CONVERSION;
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Locating girders"));

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> point;
   point.CoCreateInstance(CLSID_Point2d);
   auto station_and_offset = [&](const Eigen::Vector2d& p)
      {
         point->Move(p.x(), p.y());
         Float64 station, offset;
         pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &station, &offset);
         return std::make_pair(station, offset);
      };

   SpanIndexType nSpans = pier_stations.size() - 1;
   GirderLayout layout(nSpans);

   auto beam_ids = get_beam_ids(file);

   ifcopenshell::geom::settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);

   // set up filter for the geometry iterator
   ifcopenshell::geom::instance_id_filter filter(true, false, beam_ids);
   std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

   std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
   int num_threads = 1;// std::thread::hardware_concurrency();
   ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, num_threads);
   if (!iterator.initialize())
   {
      WBFL::System::Logger::Info(_T("Unable to process beam geometry. Girders will not be located."));
      return layout;
   }

   // This do loop can be multi-threaded
   do
   {
      auto element = iterator.get();
      auto beam = element->product().as<IfcSchema::IfcBeam>();

      std::_tostringstream os;
      os << _T("Processing geometry for ") << A2T(beam.Name().value_or("unnamed beam").c_str()) << std::endl;
      pProgress->UpdateMessage(os.str().c_str());

      auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
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
      if (segments.size() < 2)
      {
         std::ostringstream msg;
         msg << "Unable to identify the ends of " << beam.Name().value_or("unnamed beam") << " from the beam geometry.";
         WBFL::System::Logger::Info(msg.str().c_str());
         continue;
      }

      // 6. Sort segments by length
      std::sort(segments.begin(), segments.end(),
         [](const Wire& a, const Wire& b)
         {
            return a.plan_length() < b.plan_length();
         });

      // 7. The two shortest segments are the ends of the top of the girder. Their centers are on the CL girder
      Eigen::Vector3d c1 = segments[0].bbox_center();
      Eigen::Vector3d c2 = segments[1].bbox_center();

      GirderCenterline girder;
      girder.beam_id = beam.id();
      girder.mesh = std::make_shared<Mesh>(std::move(m));
      girder.start = c1.head<2>();
      girder.end = c2.head<2>();
      std::tie(girder.start_station, girder.start_offset) = station_and_offset(girder.start);
      std::tie(girder.end_station, girder.end_offset) = station_and_offset(girder.end);
      if (girder.end_station < girder.start_station)
      {
         std::swap(girder.start, girder.end);
         std::swap(girder.start_station, girder.end_station);
         std::swap(girder.start_offset, girder.end_offset);
      }

      // 8. The girder is in the span that contains its mid-point
      Float64 mid_station = (girder.start_station + girder.end_station) / 2;
      auto next_pier = std::upper_bound(pier_stations.begin(), pier_stations.end(), mid_station);
      SpanIndexType spanIdx = (SpanIndexType)std::clamp<std::ptrdiff_t>(std::distance(pier_stations.begin(), next_pier) - 1, 0, (std::ptrdiff_t)nSpans - 1);
      layout[spanIdx].push_back(girder);
   } while (iterator.next());

   // 9. Order girders left to right (offsets are positive to the right) and assign girder keys
   for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
   {
      auto& girders = layout[spanIdx];
      std::sort(girders.begin(), girders.end(), [](const auto& a, const auto& b) {return a.start_offset + a.end_offset < b.start_offset + b.end_offset; });

      for (GirderIndexType gdrIdx = 0; gdrIdx < girders.size(); gdrIdx++)
      {
         auto& girder = girders[gdrIdx];
         girder.girder_key = CGirderKey(spanIdx, gdrIdx);

         // Cross check with the girder designation in the model
         auto beam = file.instance_by_id(girder.beam_id).as<IfcSchema::IfcBeam>();
         auto named_key = get_girder_key(beam);
         if (named_key != CGirderKey() && named_key != girder.girder_key)
         {
            std::_tostringstream msg;
            msg << A2T(beam.Name().value_or("unnamed beam").c_str()) << _T(" is designated ") << LABEL_GIRDER(named_key)
                << _T(" in the model, but its location is ") << LABEL_GIRDER(girder.girder_key) << _T(". The location is used.");
            WBFL::System::Logger::Info(msg.str().c_str());
         }
      }
   }

   return layout;
}

// Girder spacing measured along a line through the alignment at the pier stations
std::pair<Spacing, Spacing> get_beam_spacing(std::shared_ptr<WBFL::EAF::Broker> pBroker, const GirderLayout& layout, const std::vector<Float64>& pier_stations, const std::vector<Eigen::Vector2d>* pLineDirections)
{
   GET_IFACE2(pBroker, IRoadway, pAlignment);

   // distance along the line through the alignment at station, positive to the right, to where each girder CL crosses it.
   // The line is normal to the alignment unless a direction is given.
   auto offsets_at = [&](const std::vector<GirderCenterline>& girders, PierIndexType pierIdx)
      {
         Float64 station = pier_stations[pierIdx];
         CComPtr<IPoint2d> alignment_point;
         pAlignment->GetPoint(station, 0.0, nullptr, pgsTypes::PlanCoordinateType::pcGlobal, &alignment_point);
         Float64 x, y;
         alignment_point->Location(&x, &y);
         Eigen::Vector2d a(x, y);

         CComPtr<IDirection> normal;
         pAlignment->GetBearingNormal(station, &normal); // normal to the right
         Float64 angle;
         normal->get_Value(&angle);
         Eigen::Vector2d right(cos(angle), sin(angle));

         Eigen::Vector2d n = (pLineDirections ? (*pLineDirections)[pierIdx].normalized() : right);
         if (n.dot(right) < 0)
            n = -n;

         std::vector<Float64> offsets;
         for (const auto& girder : girders)
         {
            // intersect start + t*(end - start) with a + u*n
            Eigen::Vector2d d = girder.end - girder.start;
            Eigen::Matrix2d A;
            A << d.x(), -n.x(),
                 d.y(), -n.y();
            Eigen::Vector2d tu = A.colPivHouseholderQr().solve(a - girder.start);
            offsets.push_back(tu(1));
         }
         return offsets;
      };

   auto spacing = [](const std::vector<Float64>& offsets)
      {
         std::vector<Float64> s;
         for (size_t i = 1; i < offsets.size(); i++)
            s.push_back(offsets[i] - offsets[i - 1]);
         return s;
      };

   Spacing start_spacing, end_spacing;
   for (SpanIndexType spanIdx = 0; spanIdx < layout.size(); spanIdx++)
   {
      start_spacing[spanIdx] = spacing(offsets_at(layout[spanIdx], spanIdx));
      end_spacing[spanIdx] = spacing(offsets_at(layout[spanIdx], spanIdx + 1));
   }

   return { start_spacing, end_spacing };
}
