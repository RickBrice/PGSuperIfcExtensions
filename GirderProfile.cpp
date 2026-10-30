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
#include "GirderProfile.h"

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/conversion_settings.h>
#include <ifcgeom/taxonomy.h>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   namespace taxonomy = ifcopenshell::geom::taxonomy;

   // The 4x4 transformation of a placement (object placement, axis placement, or cartesian transformation operator), in
   // metres (as in GirderReinforcement.cpp)
   template <typename Placement>
   Eigen::Matrix4d placement_matrix(ifcopenshell::geom::abstract_mapping& mapping, Placement placement)
   {
      if (!placement)
         return Eigen::Matrix4d::Identity();
      auto matrix = taxonomy::dcast<taxonomy::matrix4>(mapping.map(placement));
      return matrix ? matrix->ccomponents() : Eigen::Matrix4d::Identity();
   }

   Eigen::Vector3d transform(const Eigen::Matrix4d& m, const Eigen::Vector3d& p)
   {
      return (m * Eigen::Vector4d(p.x(), p.y(), p.z(), 1.0)).head<3>();
   }

   Eigen::Vector3d direction(IfcSchema::IfcDirection direction)
   {
      auto ratios = direction.DirectionRatios();
      return Eigen::Vector3d(ratios[0], 1 < ratios.size() ? ratios[1] : 0.0, 2 < ratios.size() ? ratios[2] : 0.0).normalized();
   }

   // The loops of a profile placed in 3D (global coordinates, system units)
   struct ProfileLoops
   {
      std::vector<std::vector<Eigen::Vector3d>> solids, voids;
   };

   Eigen::Matrix4d item_matrix(const taxonomy::geom_item& item)
   {
      return item.matrix ? item.matrix->ccomponents() : Eigen::Matrix4d::Identity();
   }

   void add_loop(const taxonomy::loop& loop, const Eigen::Matrix4d& m, bool external, const ifcopenshell::geom::settings& settings, ProfileLoops& loops)
   {
      std::vector<Eigen::Vector3d> points;
      for (const auto& point : loop.to_polygon(settings)) // arcs as short straight pieces, the loop's own placement applied
      {
         Eigen::Vector3d p = transform(m, point.ccomponents());
         if (points.empty() || 1e-9 < (points.back() - p).norm())
            points.push_back(p);
      }
      if (1 < points.size() && IsZero((points.back() - points.front()).norm(), 1e-9))
         points.pop_back(); // closed loops repeat the first point
      if (3 <= points.size())
         (external ? loops.solids : loops.voids).push_back(std::move(points));
   }

   // The loops of a mapped profile: a loop (parameterized profiles), a face (arbitrary profiles, the first loop the
   // outside unless it says otherwise), or a collection (composite profiles)
   void collect_loops(const taxonomy::ptr& item, const Eigen::Matrix4d& m, const ifcopenshell::geom::settings& settings, ProfileLoops& loops)
   {
      if (auto loop = taxonomy::dcast<taxonomy::loop>(item))
      {
         add_loop(*loop, m, loop->external.value_or(true), settings, loops);
      }
      else if (auto face = taxonomy::dcast<taxonomy::face>(item))
      {
         Eigen::Matrix4d face_m = m * item_matrix(*face);
         for (size_t i = 0; i < face->children.size(); i++)
            add_loop(*face->children[i], face_m, face->children[i]->external.value_or(i == 0), settings, loops);
      }
      else if (auto collection = taxonomy::dcast<taxonomy::collection>(item))
      {
         Eigen::Matrix4d collection_m = m * item_matrix(*collection);
         for (const auto& child : collection->children)
            collect_loops(child, collection_m, settings, loops);
      }
   }

   // The section between two sections with the same loops (point by point, as a loft between them), at a fraction t
   // from a to b. nullopt if their loops don't correspond
   std::optional<ProfileLoops> interpolate(const ProfileLoops& a, const ProfileLoops& b, double t)
   {
      auto blend = [t](const std::vector<std::vector<Eigen::Vector3d>>& la, const std::vector<std::vector<Eigen::Vector3d>>& lb, std::vector<std::vector<Eigen::Vector3d>>& result)
         {
            if (la.size() != lb.size())
               return false;
            for (size_t i = 0; i < la.size(); i++)
            {
               if (la[i].size() != lb[i].size())
                  return false;
               std::vector<Eigen::Vector3d> loop;
               for (size_t j = 0; j < la[i].size(); j++)
                  loop.push_back((1 - t) * la[i][j] + t * lb[i][j]);
               result.push_back(std::move(loop));
            }
            return true;
         };
      ProfileLoops loops;
      if (!blend(a.solids, b.solids, loops.solids) || !blend(a.voids, b.voids, loops.voids))
         return std::nullopt;
      return loops;
   }

   // A swept solid of the body and its section at mid-length (global coordinates)
   struct Sweep
   {
      std::string type;                  // its IFC class
      std::optional<ProfileLoops> loops; // the section at mid-length, or nullopt (why_not says why)
      std::string why_not;
      Eigen::Vector3d direction;         // of the sweep
      Eigen::Vector3d normal;            // of the profile's plane
   };

   class SweepCollector
   {
   public:
      SweepCollector(ifcopenshell::geom::abstract_mapping& mapping, const ifcopenshell::geom::settings& settings, const Eigen::Vector3d& girder_axis) :
         m_Mapping(mapping), m_Settings(settings), m_GirderAxis(girder_axis)
      {
      }

      std::vector<Sweep> sweeps;
      std::set<std::string> others; // the IFC classes of the other items of the body

      // the swept solids of a representation item (m: from the item's coordinates to global). Mapped items and the
      // solids of clipping results (the girder ends cut by half spaces) are followed
      void Collect(IfcSchema::IfcRepresentationItem item, const Eigen::Matrix4d& m)
      {
         if (auto mapped = item.as<IfcSchema::IfcMappedItem>())
         {
            auto source = mapped.MappingSource();
            Eigen::Matrix4d to_global = m * placement_matrix(m_Mapping, mapped.MappingTarget()) * placement_matrix(m_Mapping, source.MappingOrigin());
            for (const auto& mapped_item : source.MappedRepresentation().Items())
               Collect(mapped_item, to_global);
         }
         else if (auto clipping = item.as<IfcSchema::IfcBooleanClippingResult>())
         {
            if (auto first = clipping.FirstOperand().as<IfcSchema::IfcRepresentationItem>())
               Collect(first, m);
            else
               others.insert(clipping.declaration().name());
         }
         else if (auto extrusion = item.as<IfcSchema::IfcExtrudedAreaSolid>())
         {
            CollectExtrusion(extrusion, m);
         }
         else if (auto sectioned = item.as<IfcSchema::IfcSectionedSolidHorizontal>())
         {
            CollectSectionedSolid(sectioned, m);
         }
         else
         {
            others.insert(item.declaration().name());
         }
      }

   private:
      ifcopenshell::geom::abstract_mapping& m_Mapping;
      const ifcopenshell::geom::settings& m_Settings;
      Eigen::Vector3d m_GirderAxis;

      std::optional<ProfileLoops> GetLoops(IfcSchema::IfcProfileDef profile, const Eigen::Matrix4d& m)
      {
         try
         {
            auto item = m_Mapping.map(profile);
            if (!item)
               return std::nullopt;
            ProfileLoops loops;
            collect_loops(item, m, m_Settings, loops);
            if (loops.solids.empty())
               return std::nullopt;
            return loops;
         }
         catch (...)
         {
            return std::nullopt;
         }
      }

      // An extrusion, or a tapered extrusion with its section halfway between its start and end profiles
      void CollectExtrusion(IfcSchema::IfcExtrudedAreaSolid extrusion, const Eigen::Matrix4d& m)
      {
         Sweep sweep;
         sweep.type = extrusion.declaration().name();
         Eigen::Matrix4d to_global = m * placement_matrix(m_Mapping, extrusion.Position());
         Eigen::Matrix3d rotation = to_global.block<3, 3>(0, 0);
         sweep.direction = (rotation * direction(extrusion.ExtrudedDirection())).normalized();
         sweep.normal = (rotation * Eigen::Vector3d::UnitZ()).normalized();
         sweep.loops = GetLoops(extrusion.SweptArea(), to_global);
         if (!sweep.loops)
         {
            sweep.why_not = "its profile can't be read";
         }
         else if (auto tapered = extrusion.as<IfcSchema::IfcExtrudedAreaSolidTapered>())
         {
            auto end = GetLoops(tapered.EndSweptArea(), to_global);
            sweep.loops = end ? interpolate(*sweep.loops, *end, 0.5) : std::nullopt;
            if (!sweep.loops)
               sweep.why_not = "its end profile can't be read or doesn't correspond to its start profile";
         }
         sweeps.push_back(std::move(sweep));
      }

      // A sectioned solid: its cross sections before and after mid-length, interpolated. Each cross section is placed as
      // IfcOpenShell places it (make_loft in infra_sweep_helper.cpp): the profile's y axis is the placement's Axis (else
      // up), its plane is normal to the RefDirection (else to the directrix), and its x axis is y cross the normal (to the
      // left of the directrix). The directrix is taken to run along the girder
      void CollectSectionedSolid(IfcSchema::IfcSectionedSolidHorizontal solid, const Eigen::Matrix4d& m)
      {
         Sweep sweep;
         sweep.type = solid.declaration().name();
         sweep.direction = m_GirderAxis;
         std::vector<IfcSchema::IfcProfileDef> profiles = solid.CrossSections();
         std::vector<IfcSchema::IfcAxis2PlacementLinear> positions = solid.CrossSectionPositions();
         if (profiles.size() != positions.size() || profiles.size() < 2)
         {
            sweep.why_not = "its cross sections and their positions don't correspond";
            sweeps.push_back(std::move(sweep));
            return;
         }

         std::vector<std::pair<double, size_t>> stations; // distance along the directrix, cross section
         for (size_t i = 0; i < positions.size(); i++)
         {
            auto location = positions[i].Location().as<IfcSchema::IfcPointByDistanceExpression>();
            auto distance = location ? location.DistanceAlong().as<IfcSchema::IfcLengthMeasure>() : IfcSchema::IfcLengthMeasure();
            if (!distance)
            {
               sweep.why_not = "a cross section isn't positioned by a length along the directrix";
               sweeps.push_back(std::move(sweep));
               return;
            }
            stations.emplace_back((double)distance, i);
         }
         std::sort(stations.begin(), stations.end());
         double middle = (stations.front().first + stations.back().first) / 2;
         size_t k = 1;
         while (k + 1 < stations.size() && stations[k].first < middle)
            k++;
         const auto& [before, i] = stations[k - 1];
         const auto& [after, j] = stations[k];
         double t = IsEqual(before, after) ? 0.0 : (middle - before) / (after - before);

         Eigen::Matrix3d rotation = m.block<3, 3>(0, 0);
         Eigen::Vector3d lateral = Eigen::Vector3d::UnitZ().cross(m_GirderAxis).normalized();
         auto place = [&](size_t index, Eigen::Vector3d& normal)
            {
               const auto& position = positions[index];
               auto axis = position.Axis();
               auto ref_direction = position.RefDirection();
               Eigen::Vector3d y = rotation * (axis ? direction(axis) : Eigen::Vector3d(Eigen::Vector3d::UnitZ()));
               y.normalize();
               Eigen::Vector3d n = ref_direction ? (rotation * direction(ref_direction)).normalized() : m_GirderAxis;
               Eigen::Vector3d x = (axis || ref_direction) ? y.cross(n) : lateral;
               if (x.norm() < 1e-9)
                  x = lateral - lateral.dot(y) * y; // the axis is normal to the profile's plane
               x.normalize();
               normal = x.cross(y);
               Eigen::Matrix4d to_global = Eigen::Matrix4d::Identity();
               to_global.block<3, 1>(0, 0) = x;
               to_global.block<3, 1>(0, 1) = y;
               to_global.block<3, 1>(0, 2) = normal;
               return GetLoops(profiles[index], to_global);
            };
         Eigen::Vector3d normal_a, normal_b;
         auto a = place(i, normal_a);
         auto b = place(j, normal_b);
         sweep.normal = (1 - t) * normal_a + t * normal_b;
         if (!a || !b)
            sweep.why_not = "a cross section's profile can't be read";
         else if (!(sweep.loops = interpolate(*a, *b, t)))
            sweep.why_not = "the cross sections around mid-length don't correspond";
         sweeps.push_back(std::move(sweep));
      }
   };
}

GirderProfileReader::GirderProfileReader(ifcopenshell::file& file)
{
   m_Settings = std::make_unique<ifcopenshell::geom::settings>();
   m_Settings->set("mesher-linear-deflection", 0.0001); // arcs within 0.1 mm (a curve's deviation from its straight pieces)
   m_Mapping.reset(ifcopenshell::geom::impl::mapping_implementations().construct(&file, *m_Settings));
}

GirderProfileReader::~GirderProfileReader() = default;

std::optional<SectionOutline> GirderProfileReader::Read(IfcSchema::IfcBeam beam, const GirderCenterline& girder, std::string& source, std::string& why_not) const
{
   if (!beam.Representation())
   {
      why_not = "it has no representation";
      return std::nullopt;
   }

   // the girder axis, and the plane CutGirderSection cuts: vertical, normal to the axis in plan
   Eigen::Vector3d plan_axis(girder.end.x() - girder.start.x(), girder.end.y() - girder.start.y(), 0.0);
   if (IsZero(plan_axis.norm()))
   {
      why_not = "the girder has no length";
      return std::nullopt;
   }
   plan_axis.normalize();
   Eigen::Vector3d girder_axis = (Eigen::Vector3d(girder.end.x(), girder.end.y(), girder.end_top) - Eigen::Vector3d(girder.start.x(), girder.start.y(), girder.start_top)).normalized();
   Eigen::Vector3d left = Eigen::Vector3d::UnitZ().cross(plan_axis);

   SweepCollector collector(*m_Mapping, *m_Settings, girder_axis);
   Eigen::Matrix4d to_global = placement_matrix(*m_Mapping, beam.ObjectPlacement());
   for (const auto& representation : beam.Representation().Representations())
   {
      auto identifier = representation.RepresentationIdentifier();
      if (identifier && *identifier != "Body")
         continue;
      for (const auto& item : representation.Items())
         collector.Collect(item, to_global);
   }

   if (!collector.others.empty())
   {
      std::string names;
      for (const auto& name : collector.others)
         names += (names.empty() ? "" : ", ") + name;
      why_not = "its body isn't swept solids only (" + names + ")";
      return std::nullopt;
   }
   if (collector.sweeps.size() != 1)
   {
      why_not = collector.sweeps.empty() ? "its body has no swept solid" : "its body has " + std::to_string(collector.sweeps.size()) + " swept solids";
      return std::nullopt;
   }

   const auto& sweep = collector.sweeps.front();
   source = sweep.type;
   if (!sweep.loops)
   {
      why_not = "the " + sweep.type + ": " + sweep.why_not;
      return std::nullopt;
   }
   if (fabs(sweep.direction.dot(girder_axis)) < cos(M_PI / 12))
   {
      why_not = "the " + sweep.type + " isn't swept along the girder";
      return std::nullopt;
   }
   if (fabs(sweep.normal.normalized().dot(sweep.direction)) < 0.25)
   {
      why_not = "the " + sweep.type + "'s profile is nearly parallel to the girder";
      return std::nullopt;
   }

   // Each point is moved along the sweep onto the cut plane: the section the plane cuts from the solid, so a profile in a
   // skewed plane gives the normal section
   auto project = [&](const std::vector<Eigen::Vector3d>& loop)
      {
         SectionOutline::Loop section;
         for (const auto& p : loop)
         {
            Eigen::Vector3d q = p - sweep.direction * (p.dot(plan_axis) / sweep.direction.dot(plan_axis));
            section.emplace_back(q.dot(left), q.z());
         }
         return section;
      };
   SectionOutline outline;
   for (const auto& loop : sweep.loops->solids)
      outline.solids.push_back(project(loop));
   for (const auto& loop : sweep.loops->voids)
      outline.voids.push_back(project(loop));
   return UprightSection(outline);
}
