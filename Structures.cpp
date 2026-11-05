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
#include "Structures.h"
#include "DeckSlab.h"
#include "IfcImporter.h"
#include "IfcTargetReader.h"
#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/element.h>
#include <ifcgeom/kernels/opencascade/opencascade_kernel.h>
#include <IFace/Alignment.h>
#include <EAF/EAFProgress.h>
#include <EAF/AutoProgress.h>
#include <functional>
#include <iomanip>
#include <numeric>

#undef min
#undef max

namespace
{
   // Calls f(id, vertices) for each element whose geometry can be processed. vertices are x, y, z, x, y, z, ... in world coordinates
   void for_each_geometry(ifcopenshell::file& file, const std::set<int>& ids, const std::function<void(int, const std::vector<double>&)>& f)
   {
      if (ids.empty())
         return;

      ifcopenshell::geom::settings settings;
      settings.set("use-world-coords", true);
      settings.set("weld-vertices", true);
      settings.set("disable-opening-subtractions", true);

      ifcopenshell::geom::instance_id_filter filter(true, false, ids);
      std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

      std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
      ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, 1);
      if (!iterator.initialize())
         return;

      do
      {
         auto element = iterator.get();
         auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
         if (!triangulation)
            continue;
         const auto& verts = triangulation->geometry_pointer()->verts();
         if (!verts.empty())
            f(element->product().id(), verts);
      } while (iterator.next());
   }
}

std::map<int, Eigen::Vector2d> get_plan_centers(ifcopenshell::file& file, const std::set<int>& ids)
{
   std::map<int, Eigen::Vector2d> centers;
   for_each_geometry(file, ids, [&centers](int id, const std::vector<double>& verts)
      {
         Eigen::Vector2d vmin(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
         Eigen::Vector2d vmax(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
         for (size_t i = 0; i < verts.size(); i += 3)
         {
            Eigen::Vector2d v(verts[i], verts[i + 1]);
            vmin = vmin.cwiseMin(v);
            vmax = vmax.cwiseMax(v);
         }
         centers[id] = 0.5 * (vmin + vmax);
      });
   return centers;
}

namespace
{
   // A deck slab's station and offset ranges (offsets are positive to the right of the alignment)
   struct SlabExtent
   {
      int id;
      Float64 start = DBL_MAX, end = -DBL_MAX;
      Float64 left = DBL_MAX, right = -DBL_MAX;
      size_t structure = 0;
   };

   std::string feet(Float64 value)
   {
      std::ostringstream os;
      os << std::fixed << std::setprecision(2) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Feet) << " ft";
      return os.str();
   }

   // the element roles whose elements are in one structure or another. The deck slabs are assigned by themselves
   const std::array<ElementKind, 11> structure_roles{ ElementKind::Girder, ElementKind::GirderAssembly, ElementKind::ClosureJoint, ElementKind::Tendon,
      ElementKind::TendonBundle, ElementKind::Rebar, ElementKind::ReinforcementCage, ElementKind::Haunch, ElementKind::Bearing, ElementKind::Diaphragm, ElementKind::Barrier };
}

void select_structure(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file)
{
   auto slab_ids = get_slab_ids(file);
   if (slab_ids.size() < 2)
      return; // one deck slab, one structure

   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Looking for structures side by side"));

   GET_IFACE2(pBroker, IRoadway, pAlignment);
   CComPtr<IPoint2d> point;
   point.CoCreateInstance(CLSID_Point2d);
   auto station_and_offset = [&](Float64 x, Float64 y)
      {
         point->Move(x, y);
         Float64 station, offset;
         pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, point, &station, &offset);
         return std::make_pair(station, offset);
      };

   // 1. The station and offset range of each slab
   std::vector<SlabExtent> slabs;
   for_each_geometry(file, std::set<int>(slab_ids.begin(), slab_ids.end()), [&](int id, const std::vector<double>& verts)
      {
         SlabExtent slab{ id };
         for (size_t i = 0; i < verts.size(); i += 3)
         {
            auto [station, offset] = station_and_offset(verts[i], verts[i + 1]);
            slab.start = std::min(slab.start, station);
            slab.end = std::max(slab.end, station);
            slab.left = std::min(slab.left, offset);
            slab.right = std::max(slab.right, offset);
         }
         slabs.push_back(slab);
      });
   if (slabs.size() < 2)
      return;

   // 2. Structures: slabs whose offset ranges overlap are one deck (e.g. the segments between deck joints). Slabs side by
   // side, with a joint along the bridge and no slab across it, are different structures
   const Float64 inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
   std::vector<size_t> parent(slabs.size());
   std::iota(parent.begin(), parent.end(), 0);
   std::function<size_t(size_t)> find = [&](size_t i) { return parent[i] == i ? i : parent[i] = find(parent[i]); };
   for (size_t i = 0; i < slabs.size(); i++)
   {
      for (size_t j = i + 1; j < slabs.size(); j++)
      {
         Float64 overlap = std::min(slabs[i].right, slabs[j].right) - std::max(slabs[i].left, slabs[j].left);
         if (IsGT(0.0, overlap, inch / 16))
            parent[find(i)] = find(j);
      }
   }

   std::map<size_t, size_t> structure_index; // root -> structure number, in order of offset
   {
      std::map<size_t, Float64> left_of; // root -> smallest offset
      for (size_t i = 0; i < slabs.size(); i++)
      {
         auto root = find(i);
         auto [it, bNew] = left_of.emplace(root, slabs[i].left);
         if (!bNew)
            it->second = std::min(it->second, slabs[i].left);
      }
      if (left_of.size() < 2)
         return; // one structure

      std::vector<std::pair<Float64, size_t>> order;
      for (const auto& [root, left] : left_of)
         order.emplace_back(left, root);
      std::sort(order.begin(), order.end());
      for (size_t k = 0; k < order.size(); k++)
         structure_index[order[k].second] = k;
   }
   for (size_t i = 0; i < slabs.size(); i++)
      slabs[i].structure = structure_index[find(i)];
   size_t nStructures = structure_index.size();

   struct Structure
   {
      Float64 left = DBL_MAX, right = -DBL_MAX;
      std::vector<int> slab_ids;
   };
   std::vector<Structure> structures(nStructures);
   for (const auto& slab : slabs)
   {
      auto& s = structures[slab.structure];
      s.left = std::min(s.left, slab.left);
      s.right = std::max(s.right, slab.right);
      s.slab_ids.push_back(slab.id);
   }

   auto describe_structures = [&]()
      {
         std::ostringstream os;
         for (size_t k = 0; k < nStructures; k++)
         {
            os << (k == 0 ? "" : (k + 1 == nStructures ? ", and " : ", ")) << feet(structures[k].left) << " to " << feet(structures[k].right)
               << " (" << structures[k].slab_ids.size() << (structures[k].slab_ids.size() == 1 ? " slab)" : " slabs)");
         }
         return os.str();
      };

   // 3. The structure on the alignment
   std::optional<size_t> imported;
   for (size_t k = 0; k < nStructures; k++)
   {
      if (InRange(structures[k].left, 0.0, structures[k].right))
         imported = k;
   }
   if (!imported)
   {
      std::ostringstream os;
      os << "The deck is " << nStructures << " structures side by side, with joints along the bridge (deck offsets " << describe_structures()
         << "). None of them is on the alignment, so they're imported as one bridge. PGSuper models one structure: the roadway section and barriers may not be right.";
      WBFL::System::Logger::Warning(os.str());
      return;
   }

   // 4. The elements of the other structures: each element goes with the structure whose deck is over its plan center
   // (a slab over its station and offset), else with the structure nearest in offset
   const auto& reader = CIfcImporter::GetTargetReader();
   std::map<int, ElementKind> role_of;
   for (auto role : structure_roles)
   {
      for (const auto& element : reader.Select(role, file))
         role_of.emplace(element.id(), role);
   }
   std::set<int> element_ids;
   for (const auto& [id, role] : role_of)
      element_ids.insert(id);
   auto centers = get_plan_centers(file, element_ids);

   std::set<int> excluded;
   IndexType nExcludedSlabs = 0;
   for (const auto& slab : slabs)
   {
      if (slab.structure != *imported)
      {
         excluded.insert(slab.id);
         nExcludedSlabs++;
      }
   }
   std::map<ElementKind, IndexType> excluded_by_role;
   for (const auto& [id, center] : centers)
   {
      auto [station, offset] = station_and_offset(center.x(), center.y());
      std::optional<size_t> structure;
      for (const auto& slab : slabs)
      {
         if (InRange(slab.start, station, slab.end) && InRange(slab.left, offset, slab.right))
         {
            structure = slab.structure;
            break;
         }
      }
      if (!structure)
      {
         Float64 nearest = DBL_MAX;
         for (size_t k = 0; k < nStructures; k++)
         {
            Float64 distance = offset < structures[k].left ? structures[k].left - offset : (structures[k].right < offset ? offset - structures[k].right : 0.0);
            if (distance < nearest)
            {
               nearest = distance;
               structure = k;
            }
         }
      }
      if (*structure != *imported)
      {
         excluded.insert(id);
         excluded_by_role[role_of[id]]++;
      }
   }

   CIfcImporter::ExcludeElements(excluded);

   std::ostringstream os;
   os << "The deck is " << nStructures << " structures side by side, with joints along the bridge and no deck slab across them (deck offsets "
      << describe_structures() << "). PGSuper models one structure, so the one on the alignment (" << feet(structures[*imported].left) << " to "
      << feet(structures[*imported].right) << ") is imported. The others' elements aren't: " << nExcludedSlabs << " deck";
   for (const auto& [role, count] : excluded_by_role)
      os << ", " << count << " " << GetElementRoleName(role);
   os << " elements.";
   WBFL::System::Logger::Warning(os.str());
}
