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
#include "TrafficBarriers.h"
#include "GirderSection.h"
#include "IfcImporter.h"
#include "IfcTargetReader.h"

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/element.h>
#include <ifcgeom/kernels/opencascade/opencascade_kernel.h>

#include <IFace\Project.h>
#include <IFace/Alignment.h>
#include <psgLib/TrafficBarrierEntry.h>
#include <psgLib/LibraryManager.h>
#include <EAF/EAFProgress.h>
#include <EAF/AutoProgress.h>

#include <iomanip>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   using Loop = SectionOutline::Loop;

   double inches(double length) { return WBFL::Units::ConvertFromSysUnits(length, WBFL::Units::Measure::Inch); }

   double loop_area(const Loop& loop)
   {
      double a = 0;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& p = loop[i];
         const auto& q = loop[(i + 1) % loop.size()];
         a += p.x() * q.y() - q.x() * p.y();
      }
      return a / 2;
   }

   Eigen::Vector2d loop_center(const Loop& loop)
   {
      Eigen::Vector2d c(0, 0);
      for (const auto& p : loop)
         c += p;
      return loop.empty() ? c : Eigen::Vector2d(c / (double)loop.size());
   }

   double sign_of(pgsTypes::TrafficBarrierOrientation side) { return side == pgsTypes::tboLeft ? 1.0 : -1.0; }

   // The edge of the deck on one side, in the cut's coordinates (distance to the left of the alignment, elevation): the
   // outermost point of the deck section, at the top of the deck there
   Eigen::Vector2d deck_edge(const Loop& deck, pgsTypes::TrafficBarrierOrientation side)
   {
      double sign = sign_of(side);
      double edge = -DBL_MAX;
      for (const auto& p : deck)
         edge = std::max(edge, sign * p.x());
      const double tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 16.0, WBFL::Units::Measure::Inch);
      double top = -DBL_MAX;
      for (const auto& p : deck)
      {
         if (edge - tolerance <= sign * p.x())
            top = std::max(top, p.y());
      }
      return Eigen::Vector2d(sign * edge, top);
   }

   // The points of a loop moved to the edge of the deck at the top of the deck, positive toward the roadway
   Loop from_edge(const Loop& loop, const Eigen::Vector2d& edge, pgsTypes::TrafficBarrierOrientation side)
   {
      double sign = sign_of(side);
      Loop points;
      for (const auto& p : loop)
         points.emplace_back(sign * (edge.x() - p.x()), p.y() - edge.y());
      return points;
   }

   // The highest point of a loop at x, or nullopt if the loop doesn't reach x
   std::optional<double> top_at(const Loop& loop, double x)
   {
      std::optional<double> top;
      for (size_t i = 0; i < loop.size(); i++)
      {
         const auto& a = loop[i];
         const auto& b = loop[(i + 1) % loop.size()];
         if ((a.x() <= x) != (b.x() <= x))
            top = std::max(top.value_or(-DBL_MAX), a.y() + (x - a.x()) * (b.y() - a.y()) / (b.x() - a.x()));
      }
      return top;
   }

   // A barrier's section in traffic barrier coordinates: from the edge of the deck at the top of the deck, positive
   // toward the roadway, counterclockwise. PGSuper turns a barrier to the cross slope of the roadway at the edge of the
   // deck (WBFL CSectionCutTool::CreateBarrierShape), so the section is turned back by that slope (rise toward the roadway)
   Loop barrier_coordinates(const Loop& loop, const Eigen::Vector2d& edge, pgsTypes::TrafficBarrierOrientation side, double slope)
   {
      auto points = from_edge(loop, edge, side);
      Eigen::Rotation2Dd rotation(-atan(slope));
      for (auto& p : points)
         p = rotation * p;
      if (loop_area(points) < 0)
         std::reverse(points.begin(), points.end());
      return points;
   }

   Loop entry_points(const TrafficBarrierEntry* pEntry)
   {
      Loop loop;
      CComPtr<IPoint2dCollection> points;
      pEntry->GetBarrierPoints(&points);
      IndexType n = 0;
      if (points)
         points->get_Count(&n);
      for (IndexType i = 0; i < n; i++)
      {
         CComPtr<IPoint2d> point;
         points->get_Item(i, &point);
         Float64 x, y;
         point->Location(&x, &y);
         if (loop.empty() || 1e-9 < (loop.back() - Eigen::Vector2d(x, y)).norm())
            loop.emplace_back(x, y);
      }
      if (1 < loop.size() && (loop.back() - loop.front()).norm() < 1e-9)
         loop.pop_back(); // closed polygons may repeat the first point
      return loop;
   }

   SectionOutline outline_of(const Loop& loop)
   {
      SectionOutline outline;
      outline.solids.push_back(loop);
      return outline;
   }

   struct BarrierCut
   {
      IfcSchema::IfcObject element;
      std::string name;
      Loop loop; // the largest loop of the cut, in the cut's coordinates
   };
}

std::vector<ElementMesh> get_barrier_meshes(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file)
{
   const auto& reader = CIfcImporter::GetTargetReader();
   if (!reader.HasSelector(ElementKind::Barrier))
      return {};

   auto barriers = reader.Select(ElementKind::Barrier, file);
   if (barriers.empty())
   {
      WBFL::System::Logger::Info("The model has no traffic barriers (" + reader.GetTable().GetSelector(ElementKind::Barrier)->Describe() + ").");
      return {};
   }

   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Processing traffic barrier geometry"));

   std::set<int> ids;
   for (const auto& barrier : barriers)
      ids.insert(barrier.id());

   ifcopenshell::geom::settings settings;
   settings.set("use-world-coords", true);
   settings.set("weld-vertices", true);
   settings.set("disable-opening-subtractions", true);
   ifcopenshell::geom::instance_id_filter filter(true, false, ids);
   std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });
   std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
   ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, 1);
   if (!iterator.initialize())
   {
      WBFL::System::Logger::Info("Unable to process the traffic barrier geometry.");
      return {};
   }

   std::vector<ElementMesh> meshes;
   do
   {
      auto element = iterator.get();
      auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
      if (!triangulation)
         continue;
      auto geometry = triangulation->geometry_pointer();
      const auto& verts = geometry->verts();
      const auto& faces = geometry->faces();
      std::vector<Eigen::Vector3d> v;
      v.reserve(verts.size() / 3);
      for (size_t i = 0; i + 2 < verts.size(); i += 3)
         v.emplace_back(verts[i], verts[i + 1], verts[i + 2]);
      std::vector<std::array<int, 3>> f;
      f.reserve(faces.size() / 3);
      for (size_t i = 0; i + 2 < faces.size(); i += 3)
         f.emplace_back(std::array<int, 3>({ faces[i], faces[i + 1], faces[i + 2] }));
      meshes.push_back({ element->product().as<IfcSchema::IfcObject>(), Mesh(v, f) });
   } while (iterator.next());
   return meshes;
}

AlignmentCut get_alignment_cut(std::shared_ptr<WBFL::EAF::Broker> pBroker, Float64 station)
{
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> alignment_point;
   pAlignment->GetPoint(station, 0.0, nullptr, pgsTypes::PlanCoordinateType::pcGlobal, &alignment_point);
   Float64 x, y;
   alignment_point->Location(&x, &y);
   CComPtr<IDirection> bearing;
   pAlignment->GetBearing(station, &bearing);
   Float64 angle;
   bearing->get_Value(&angle);
   return { station, Eigen::Vector2d(x, y), Eigen::Vector2d(cos(angle), sin(angle)) };
}

void import_traffic_barriers(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, const std::vector<Float64>& pier_stations, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;
   if (barriers.empty() || pier_stations.size() < 2)
      return;

   const auto& reader = CIfcImporter::GetTargetReader();
   GET_IFACE2(pBroker, IRoadway, pAlignment);

   // the cut: normal to the alignment at the middle of the bridge (barrier end sections and transitions are at the ends)
   auto cut_plane = get_alignment_cut(pBroker, (pier_stations.front() + pier_stations.back()) / 2);
   Float64 station = cut_plane.station;
   auto deck_loops = CutMesh(deck, cut_plane.point, cut_plane.axis);
   if (deck_loops.empty())
   {
      WBFL::System::Logger::Info("The deck couldn't be cut at the middle of the bridge, so the traffic barriers can't be located. The railing systems are not changed.");
      return;
   }

   // the barrier sections at the cut, by side
   std::array<std::vector<BarrierCut>, 2> sides; // pgsTypes::tboLeft, tboRight
   for (const auto& barrier : barriers)
   {
      auto loops = CutMesh(barrier.mesh, cut_plane.point, cut_plane.axis);
      if (loops.empty())
         continue; // e.g. an end section, which doesn't reach the middle of the bridge

      auto side = (0 < loop_center(loops.front()).x() ? pgsTypes::tboLeft : pgsTypes::tboRight);
      sides[side].push_back({ barrier.element, barrier.element.Name().value_or("unnamed barrier"), loops.front() });
   }

   GET_IFACE2(pBroker, ILibrary, pLibrary);
   auto& library = pLibrary->GetTrafficBarrierLibrary();
   GET_IFACE2(pBroker, ILibraryNames, pLibraryNames);
   std::vector<std::_tstring> entry_names;
   pLibraryNames->EnumTrafficBarrierNames(&entry_names);

   const double tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 8.0, WBFL::Units::Measure::Inch);
   for (auto side : { pgsTypes::tboLeft, pgsTypes::tboRight })
   {
      std::string side_name(side == pgsTypes::tboLeft ? "left" : "right");
      auto& cuts = sides[side];
      if (cuts.empty())
      {
         WBFL::System::Logger::Info("No traffic barrier is on the " + side_name + " side at the middle of the bridge. The " + side_name + " railing system is not changed.");
         continue;
      }

      // the exterior barrier is the outermost
      double sign = sign_of(side);
      std::sort(cuts.begin(), cuts.end(), [sign](const BarrierCut& a, const BarrierCut& b) {return sign * loop_center(a.loop).x() > sign * loop_center(b.loop).x(); });
      const auto& cut = cuts.front();

      auto edge = deck_edge(deck_loops.front(), side);
      // the roadway's cross slope at the edge of the deck (PGSuper's, so it's undone exactly), as a rise toward the
      // roadway. Offsets are positive to the right
      double slope = sign * pAlignment->GetSlope(station, -edge.x());

      // The slope of the top of the deck just inside the barrier (6 to 18 in from it) is a check: if it isn't the
      // roadway's, PGSuper's roadway section doesn't match the model, and the barrier is turned by the wrong slope
      std::optional<double> deck_slope;
      {
         auto barrier_points = from_edge(cut.loop, edge, side);
         double inside = -DBL_MAX;
         for (const auto& p : barrier_points)
            inside = std::max(inside, p.x());
         const double inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
         auto deck_points = from_edge(deck_loops.front(), edge, side);
         auto y1 = top_at(deck_points, inside + 6 * inch);
         auto y2 = top_at(deck_points, inside + 18 * inch);
         if (y1 && y2)
            deck_slope = (*y2 - *y1) / (12 * inch);
      }
      auto model = SimplifyLoop(barrier_coordinates(cut.loop, edge, side, slope), WBFL::Units::ConvertToSysUnits(1.0 / 64.0, WBFL::Units::Measure::Inch));
      auto model_outline = outline_of(model);

      std::ostringstream os;
      os << "The " << side_name << " exterior traffic barrier is " << cut.name << " (#" << cut.element.id() << "): area "
         << std::fixed << std::setprecision(1) << WBFL::Units::ConvertFromSysUnits(fabs(loop_area(model)), WBFL::Units::Measure::Inch2) << " in^2, "
         << model.size() << " corners, turned back " << std::setprecision(2) << 100 * slope << "% (the roadway cross slope at the edge of the deck, rising toward the roadway).";
      if (deck_slope && 0.0025 < fabs(*deck_slope - slope))
         os << " The top of the deck inside the barrier slopes " << std::setprecision(2) << 100 * *deck_slope
            << "%, so the roadway section doesn't match the model (it isn't imported yet), and the barrier is turned by the wrong slope.";
      for (size_t i = 1; i < cuts.size(); i++)
         os << " " << cuts[i].name << " (#" << cuts[i].element.id() << ") is also on the " << side_name << " side (sidewalks and interior barriers aren't imported yet).";

      // The library entry with the same section (outlines within the tolerance). The barrier type named in the model is
      // a clue: that entry is used if its section is the same, otherwise the closest
      std::optional<std::string> type_name;
      if (auto reading = reader.Read("barrier.type", cut.element))
         type_name = std::get<std::string>(reading->value);

      const TrafficBarrierEntry* pEntry = nullptr;
      double entry_distance = DBL_MAX;
      std::string closest;
      double closest_distance = DBL_MAX;
      for (const auto& name : entry_names)
      {
         const auto* pCandidate = pLibrary->GetTrafficBarrierEntry(name.c_str());
         auto points = entry_points(pCandidate);
         if (points.size() < 3)
            continue;
         double distance = OutlineDistance(model_outline, outline_of(points));
         if (distance < closest_distance)
         {
            closest_distance = distance;
            closest = T2A(name.c_str());
         }
      }
      if (closest_distance <= tolerance)
      {
         pEntry = pLibrary->GetTrafficBarrierEntry(A2T(closest.c_str()));
         entry_distance = closest_distance;
         if (type_name && *type_name != closest)
         {
            if (const auto* pNamed = pLibrary->GetTrafficBarrierEntry(A2T(type_name->c_str())))
            {
               auto points = entry_points(pNamed);
               double distance = (points.size() < 3 ? DBL_MAX : OutlineDistance(model_outline, outline_of(points)));
               if (distance <= tolerance)
               {
                  pEntry = pNamed;
                  entry_distance = distance;
               }
            }
         }
      }

      if (pEntry)
      {
         os << " Its section is the library's " << T2A(pEntry->GetName().c_str()) << " (outlines within " << std::setprecision(3) << inches(entry_distance) << " in).";
      }
      else
      {
         if (!closest.empty())
            os << " No library traffic barrier has its section (the closest is " << closest << ", outlines " << std::setprecision(3) << inches(closest_distance) << " in apart).";

         // a new entry, named for the barrier type in the model, else the element's type, else the element
         std::string base = type_name.value_or("");
         if (base.empty())
         {
            for (const auto& rel : cut.element.IsTypedBy())
            {
               auto type_object = rel.RelatingType();
               if (type_object && type_object.Name())
                  base = *type_object.Name();
            }
         }
         if (base.empty())
            base = cut.name;
         std::_tstring name(A2T(base.c_str()));
         for (int n = 2; pLibrary->GetTrafficBarrierEntry(name.c_str()); n++)
         {
            std::_tostringstream os_name;
            os_name << A2T(base.c_str()) << _T(" (") << n << _T(")");
            name = os_name.str();
         }

         CComPtr<IPoint2dCollection> points;
         points.CoCreateInstance(CLSID_Point2dCollection);
         for (const auto& p : model)
         {
            CComPtr<IPoint2d> pt;
            pt.CoCreateInstance(CLSID_Point2d);
            pt->Move(p.x(), p.y());
            points->Add(pt);
         }
         TrafficBarrierEntry entry;
         entry.SetBarrierPoints(points);
         entry.SetWeightMethod(TrafficBarrierEntry::Compute);
         entry.SetCurbOffset(0);
         entry.IsBarrierStructurallyContinuous(false);
         library.AddEntry(entry, name.c_str());
         library.EnableEditing(name.c_str(), true);
         pEntry = pLibrary->GetTrafficBarrierEntry(name.c_str());
         entry_names.push_back(name);
         os << " A new library entry, " << T2A(name.c_str()) << ", is made from it (its weight is computed from the section). Its points (x, y in inches):";
         for (const auto& p : model)
            os << " (" << std::setprecision(3) << inches(p.x()) << ", " << inches(p.y()) << ")";
         os << ".";
      }
      WBFL::System::Logger::Info(os.str());

      auto* pRailing = (side == pgsTypes::tboLeft ? bridge_desc.GetLeftRailingSystem() : bridge_desc.GetRightRailingSystem());
      pRailing->strExteriorRailing = pEntry->GetName();
      pRailing->SetExteriorRailing(pEntry);
      pRailing->bUseSidewalk = false;
      pRailing->bUseInteriorRailing = false;
   }
}
