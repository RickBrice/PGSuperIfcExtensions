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
#include "Piers.h"
#include "Properties.h"
#include "Geometry.h"
#include "IfcImporter.h"

#include <psgLib\GirderLabel.h>
#include <IFace/Alignment.h>

#include <psgLib/BridgeDescription2.h>
#include <EAF/EAFProgress.h>
#include <EAF/AutoProgress.h>

#include <limits>
#include <numeric>

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/element.h>
#include <ifcgeom/kernels/opencascade/opencascade_kernel.h>


PierIndexType get_pier_count(ifcopenshell::file& file)
{
   auto parts = file.instances_by_type<IfcSchema::IfcBridgePart>();
   PierIndexType nPiers = 0;
   for (auto& part : parts)
   {
      if (part.PredefinedType().has_value() && (part.PredefinedType().value() == IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_ABUTMENT || part.PredefinedType().value() == IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_PIER))
         nPiers++;
   }
   return nPiers;
}

std::set<int> get_bearing_ids(IfcSchema::IfcBridgePart pier)
{
   std::set<int> bearing_ids;
   auto rels = pier.ContainsElements();
   for (auto& rel : rels)
   {
      auto related_elements = rel.RelatedElements();
      for (auto& related_element : related_elements)
      {
         auto bearing = related_element.as<IfcSchema::IfcBearing>();
         if (bearing)
         {
            bearing_ids.insert(bearing.id());
         }
      }
   }

   return bearing_ids;
}

double get_pier_station(std::shared_ptr<WBFL::EAF::Broker> pBroker,ifcopenshell::file& file, PierIndexType pierIdx, IfcSchema::IfcBridgePart pier)
{
   auto rel_positions = pier.PositionedRelativeTo();
   IfcSchema::IfcPositioningElement positioning_element;
   if (!rel_positions.empty())
      positioning_element = rel_positions.front().RelatingPositioningElement();
   auto referent = positioning_element ? positioning_element.as<IfcSchema::IfcReferent>() : IfcSchema::IfcReferent{};
   auto station = referent ? GetMeasureProperty<IfcSchema::IfcLengthMeasure>(CIfcImporter::GetUnits(), referent, "Pset_Stationing", "Station") : std::nullopt;
   if (referent && station)
   {
      return *station;
   }
   else
   {
      USES_CONVERSION;
      std::ostringstream os;
      os << "Expected Pier " << T2A(LABEL_PIER(pierIdx)) << " ";
      pier.to_string(os);
      os << " to be positioned with an IfcReferent and have stationing defined with Pset_Stationing. Attempting to estimate station from all IfcBearing in the spatial structure of the IfcBridgePart.PIER";
      WBFL::System::Logger::Info(os.str().c_str());

      ifcopenshell::geom::settings settings;
      settings.set("use-world-coords", true);
      settings.set("weld-vertices", true);
      settings.set("disable-opening-subtractions", true);

      auto bearing_ids = get_bearing_ids(pier);

      if (bearing_ids.empty())
      {
         WBFL::System::Logger::Info("IfcBearing not found in the spatial structure of this pier. Assuming piers to be at stations on 100ft increment.");
         double station = WBFL::Units::ConvertToSysUnits(pierIdx * 100.0, WBFL::Units::Measure::Feet);
         return station;
      }
      else
      {

         // set up filter for the geometry iterator
         ifcopenshell::geom::instance_id_filter filter(true, false, bearing_ids);
         std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

         std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
         int num_threads = 1;// std::thread::hardware_concurrency();
         ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, num_threads);


         Eigen::Vector3d vmax(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
         Eigen::Vector3d vmin(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
         if (!iterator.initialize())
         {
            WBFL::System::Logger::Info("Unable to process bearing geometry. Assuming piers to be at stations on 100ft increment.");
            return WBFL::Units::ConvertToSysUnits(pierIdx * 100.0, WBFL::Units::Measure::Feet);
         }

         do
         {
            auto element = iterator.get();

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

            for (auto& v_ : v)
            {
               vmin = vmin.cwiseMin(v_);
               vmax = vmax.cwiseMax(v_);
            }
         } while (iterator.next());

         Eigen::Vector3d center_point = 0.5 * (vmin + vmax);

         CComPtr<IPoint2d> pnt;
         pnt.CoCreateInstance(CLSID_Point2d);
         pnt->Move(center_point.x(), center_point.y());

         GET_IFACE2(pBroker, IRoadway, pAlignment);
         Float64 station, offset;
         pAlignment->GetStationAndOffset(pgsTypes::PlanCoordinateType::pcGlobal, pnt, &station, &offset);
         return station;
      }
   }
}


namespace
{
   Float64 median(std::vector<Float64> values)
   {
      std::sort(values.begin(), values.end());
      size_t n = values.size();
      return n % 2 ? values[n / 2] : (values[n / 2 - 1] + values[n / 2]) / 2;
   }

   // plan outline, center (of the bounding box), and elevation range of each IfcBearing
   std::vector<BearingGeometry> get_bearing_geometry(ifcopenshell::file& file)
   {
      std::set<int> bearing_ids;
      for (auto& bearing : file.instances_by_type<IfcSchema::IfcBearing>())
         bearing_ids.insert(bearing.id());

      std::vector<BearingGeometry> bearings;
      if (bearing_ids.empty())
         return bearings;

      ifcopenshell::geom::settings settings;
      settings.set("use-world-coords", true);
      settings.set("weld-vertices", true);
      settings.set("disable-opening-subtractions", true);

      ifcopenshell::geom::instance_id_filter filter(true, false, bearing_ids);
      std::vector<ifcopenshell::geom::filter_function> filters({ std::ref(filter) });

      std::unique_ptr<ifcopenshell::geom::kernels::abstract_kernel> kernel(std::make_unique<ifcopenshell::geom::open_cascade_kernel>(settings));
      ifcopenshell::geom::iterator iterator(std::move(kernel), settings, &file, filters, 1);
      if (!iterator.initialize())
         return bearings;

      do
      {
         auto element = iterator.get();
         auto triangulation = dynamic_cast<ifcopenshell::geom::triangulation_element*>(element.get());
         const auto& verts = triangulation->geometry_pointer()->verts();
         if (verts.empty())
            continue;

         BearingGeometry bearing;
         bearing.id = element->product().id();
         bearing.zmin = std::numeric_limits<double>::infinity();
         bearing.zmax = -std::numeric_limits<double>::infinity();
         Eigen::Vector2d vmin(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
         Eigen::Vector2d vmax(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
         for (size_t i = 0; i < verts.size(); i += 3)
         {
            Eigen::Vector2d v(verts[i], verts[i + 1]);
            bearing.plan.push_back(v);
            vmin = vmin.cwiseMin(v);
            vmax = vmax.cwiseMax(v);
            bearing.zmin = std::min(bearing.zmin, (double)verts[i + 2]);
            bearing.zmax = std::max(bearing.zmax, (double)verts[i + 2]);
         }
         bearing.center = 0.5 * (vmin + vmax);
         bearings.push_back(std::move(bearing));
      } while (iterator.next());

      return bearings;
   }
}

void locate_bearings(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, GirderLayout& layout)
{
   GET_IFACE2(pBroker, IEAFProgress, pProgress);
   WBFL::EAF::AutoProgress ap(pProgress);
   pProgress->UpdateMessage(_T("Locating bearings"));

   auto bearings = get_bearing_geometry(file);
   if (bearings.empty())
   {
      WBFL::System::Logger::Info(_T("No IfcBearing geometry found. Girder end distances and bearing offsets will not be imported."));
      return;
   }

   const Float64 feet = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Feet);
   const Float64 inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);

   // Each bearing supports the girder end whose CL is nearest to it, if the bearing is under the end of that girder:
   // from the end of the girder to a few feet in from the end. Bearings beyond the end of a girder support the girder
   // on the other side of the pier.
   struct GirderEnd { GirderCenterline* girder; pgsTypes::MemberEndType end; Eigen::Vector2d point; Eigen::Vector2d inward; };
   std::vector<GirderEnd> girder_ends;
   for (auto& girders : layout)
   {
      for (auto& girder : girders)
      {
         Eigen::Vector2d direction = (girder.end - girder.start).normalized();
         girder_ends.push_back({ &girder, pgsTypes::metStart, girder.start, direction });
         girder_ends.push_back({ &girder, pgsTypes::metEnd, girder.end, -direction });
      }
   }

   for (const auto& bearing : bearings)
   {
      GirderEnd* nearest = nullptr;
      Float64 nearest_lateral = 3 * feet;
      for (auto& girder_end : girder_ends)
      {
         Eigen::Vector2d v = bearing.center - girder_end.point;
         Float64 along = v.dot(girder_end.inward);
         Float64 lateral = fabs(v.x() * girder_end.inward.y() - v.y() * girder_end.inward.x());
         if (-inch / 2 < along && along < 6 * feet && lateral < nearest_lateral)
         {
            nearest = &girder_end;
            nearest_lateral = lateral;
         }
      }
      if (nearest)
         nearest->girder->bearings[nearest->end].push_back(bearing);
   }

   // CL bearing: the average distance of the bearings from the end of the girder, on the CL girder
   for (auto& girder_end : girder_ends)
   {
      const auto& end_bearings = girder_end.girder->bearings[girder_end.end];
      auto& cl_bearing = (girder_end.end == pgsTypes::metStart ? girder_end.girder->start_bearing : girder_end.girder->end_bearing);
      if (end_bearings.empty())
      {
         std::_tostringstream os;
         os << _T("No IfcBearing found at the ") << (girder_end.end == pgsTypes::metStart ? _T("start") : _T("end")) << _T(" of ") << LABEL_GIRDER(girder_end.girder->girder_key) << _T(".");
         WBFL::System::Logger::Info(os.str().c_str());
         continue;
      }

      Float64 distance = 0;
      for (const auto& bearing : end_bearings)
         distance += (bearing.center - girder_end.point).dot(girder_end.inward);
      distance /= end_bearings.size();
      cl_bearing = girder_end.point + distance * girder_end.inward;
   }
}

namespace
{
   // fixity of a bearing (fixed longitudinally, fixed transversely) from Pset_BearingCommon.DisplacementAccommodated,
   // or from a text property named Fixity (e.g. "Fixed" or "Expansion") that some authoring tools use
   std::optional<std::pair<bool, bool>> get_bearing_fixity(ifcopenshell::file& file, int bearing_id)
   {
      auto bearing = file.instance_by_id(bearing_id).as<IfcSchema::IfcBearing>();
      if (!bearing)
         return std::nullopt;

      auto lower = [](std::string s) {std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {return (char)std::tolower(c); }); return s; };

      std::optional<std::pair<bool, bool>> fixity;
      for (auto& rel : bearing.IsDefinedBy())
      {
         auto pset = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySet>();
         if (!pset)
            continue;

         for (auto& property : pset.HasProperties())
         {
            std::string name = lower(property.Name());
            if (pset.Name() == std::string("Pset_BearingCommon") && name == "displacementaccommodated")
            {
               if (auto list = property.as<IfcSchema::IfcPropertyListValue>(); list && list.ListValues())
               {
                  std::vector<bool> accommodated;
                  for (auto& value : *list.ListValues())
                  {
                     if (auto b = value.as<IfcSchema::IfcBoolean>())
                        accommodated.push_back((bool)b);
                  }
                  if (2 <= accommodated.size())
                     return std::make_pair(!accommodated[0], !accommodated[1]); // the standard property set wins
               }
            }
            else if (!fixity && name.ends_with("fixity"))
            {
               auto value = property.as<IfcSchema::IfcPropertySingleValue>();
               if (!value || !value.NominalValue())
                  continue;
               std::string text;
               if (auto label = value.NominalValue().as<IfcSchema::IfcLabel>())
                  text = lower(label);
               else if (auto t = value.NominalValue().as<IfcSchema::IfcText>())
                  text = lower(t);

               if (text.find("fixed") != std::string::npos)
                  fixity = std::make_pair(true, true);
               else if (text.find("expansion") != std::string::npos || text.find("free") != std::string::npos)
                  fixity = std::make_pair(false, false);
            }
         }
      }
      return fixity;
   }

   bool same_bearing(const CBearingData2& a, const CBearingData2& b)
   {
      const Float64 tolerance = WBFL::Units::ConvertToSysUnits(1.0 / 16.0, WBFL::Units::Measure::Inch);
      return a.Shape == b.Shape && a.BearingCount == b.BearingCount && a.FixedX == b.FixedX && a.FixedY == b.FixedY &&
         fabs(a.Length - b.Length) <= tolerance && fabs(a.Width - b.Width) <= tolerance && fabs(a.Height - b.Height) <= tolerance && fabs(a.Spacing - b.Spacing) <= tolerance;
   }
}

void set_bearing_data(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const GirderLayout& layout, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   // bearing data for each girder at each pier face, from the bearing geometry
   std::map<std::tuple<PierIndexType, pgsTypes::PierFaceType, GirderIndexType>, CBearingData2> bearing_data;
   for (SpanIndexType spanIdx = 0; spanIdx < layout.size(); spanIdx++)
   {
      for (GirderIndexType gdrIdx = 0; gdrIdx < layout[spanIdx].size(); gdrIdx++)
      {
         const auto& girder = layout[spanIdx][gdrIdx];
         Eigen::Vector2d direction = (girder.end - girder.start).normalized();
         Eigen::Vector2d normal(-direction.y(), direction.x());

         for (auto end : { pgsTypes::metStart, pgsTypes::metEnd })
         {
            const auto& bearings = girder.bearings[end];
            if (bearings.empty())
               continue;

            // start from the bearing in the template so values that can't be measured keep their defaults
            CBearingData2 bd(*bridge_desc.GetBearingData());
            bd.DefinitionType = btBasic;

            // size of the first bearing, along and normal to the girder
            const auto& bearing = bearings.front();
            Float64 min_along = Float64_Max, max_along = -Float64_Max, min_across = Float64_Max, max_across = -Float64_Max;
            std::vector<Eigen::Vector2d> plan_points; // distinct points in plan (the top and bottom of the bearing coincide in plan)
            for (const auto& p : bearing.plan)
            {
               Eigen::Vector2d v = p - bearing.center;
               min_along = std::min(min_along, v.dot(direction)); max_along = std::max(max_along, v.dot(direction));
               min_across = std::min(min_across, v.dot(normal)); max_across = std::max(max_across, v.dot(normal));
               if (std::none_of(plan_points.begin(), plan_points.end(), [&v](const auto& q) {return (q - v).norm() < 1.0e-4; }))
                  plan_points.push_back(v);
            }

            // round if many points in plan, all at the same distance from the center (the corners of a rectangle are also at the same distance)
            std::vector<Float64> radii;
            for (const auto& v : plan_points)
               radii.push_back(v.norm());
            auto [min_radius, max_radius] = std::minmax_element(radii.begin(), radii.end());
            Float64 mean_radius = std::accumulate(radii.begin(), radii.end(), 0.0) / radii.size();
            if (8 <= plan_points.size() && *max_radius - *min_radius < 0.02 * mean_radius)
            {
               bd.Shape = bsRound;
               bd.Length = bd.Width = 2 * mean_radius;
            }
            else
            {
               bd.Shape = bsRectangular;
               bd.Length = max_along - min_along;
               bd.Width = max_across - min_across;
            }
            bd.Height = bearing.zmax - bearing.zmin;

            // number of bearings and their spacing across the girder
            bd.BearingCount = bearings.size();
            bd.Spacing = 0;
            if (1 < bearings.size())
            {
               std::vector<Float64> offsets;
               for (const auto& b : bearings)
                  offsets.push_back(b.center.dot(normal));
               std::sort(offsets.begin(), offsets.end());
               bd.Spacing = (offsets.back() - offsets.front()) / (offsets.size() - 1);
            }

            if (auto fixity = get_bearing_fixity(file, bearing.id))
            {
               bd.FixedX = fixity->first;
               bd.FixedY = fixity->second;
            }

            PierIndexType pierIdx = (end == pgsTypes::metStart ? spanIdx : spanIdx + 1);
            pgsTypes::PierFaceType face = (end == pgsTypes::metStart ? pgsTypes::Ahead : pgsTypes::Back);
            bearing_data[{pierIdx, face, gdrIdx}] = bd;
         }
      }
   }

   if (bearing_data.empty())
      return;

   // the same bearing everywhere, the same bearing at each pier face, or bearings that vary by girder
   const auto& first = bearing_data.begin()->second;
   bool bSameEverywhere = std::all_of(bearing_data.begin(), bearing_data.end(), [&first](const auto& item) {return same_bearing(item.second, first); });

   std::map<std::pair<PierIndexType, pgsTypes::PierFaceType>, std::vector<const CBearingData2*>> faces;
   for (const auto& [key, bd] : bearing_data)
      faces[{std::get<0>(key), std::get<1>(key)}].push_back(&bd);
   bool bSameAtEachFace = std::all_of(faces.begin(), faces.end(), [](const auto& face) {return std::all_of(face.second.begin(), face.second.end(), [&face](const auto* bd) {return same_bearing(*bd, *face.second.front()); }); });

   auto describe = [](const CBearingData2& bd)
      {
         std::_tostringstream os;
         os << std::fixed << std::setprecision(3) << (bd.Shape == bsRound ? _T("round, diameter ") : _T("rectangular, "))
            << WBFL::Units::ConvertFromSysUnits(bd.Length, WBFL::Units::Measure::Inch) << _T(" in");
         if (bd.Shape != bsRound)
            os << _T(" x ") << WBFL::Units::ConvertFromSysUnits(bd.Width, WBFL::Units::Measure::Inch) << _T(" in");
         os << _T(", height ") << WBFL::Units::ConvertFromSysUnits(bd.Height, WBFL::Units::Measure::Inch) << _T(" in");
         if (1 < bd.BearingCount)
            os << _T(", ") << bd.BearingCount << _T(" bearings at ") << WBFL::Units::ConvertFromSysUnits(bd.Spacing, WBFL::Units::Measure::Inch) << _T(" in");
         os << (bd.FixedX || bd.FixedY ? _T(", fixed") : _T(""));
         return os.str();
      };

   if (bSameEverywhere)
   {
      bridge_desc.SetBearingType(pgsTypes::brtBridge);
      bridge_desc.SetBearingData(first);
      std::_tostringstream os;
      os << _T("Bearings: the same bearing everywhere, ") << describe(first);
      WBFL::System::Logger::Info(os.str().c_str());
   }
   else if (bSameAtEachFace)
   {
      bridge_desc.SetBearingType(pgsTypes::brtPier);
      for (const auto& [key, bds] : faces)
      {
         auto [pierIdx, face] = key;
         auto* pPier = bridge_desc.GetPier(pierIdx);
         // an abutment has girders on one face only. the other face gets the same bearing.
         for (auto f : { pgsTypes::Back, pgsTypes::Ahead })
         {
            if (f == face || pPier->IsAbutment())
               pPier->SetBearingData(f, *bds.front());
         }
         std::_tostringstream os;
         os << _T("Bearings at ") << LABEL_PIER_EX(pPier->IsAbutment(), pierIdx) << (face == pgsTypes::Ahead ? _T(" ahead") : _T(" back")) << _T(" face: ") << describe(*bds.front());
         WBFL::System::Logger::Info(os.str().c_str());
      }
   }
   else
   {
      bridge_desc.SetBearingType(pgsTypes::brtGirder);
      for (const auto& [key, bd] : bearing_data)
      {
         auto [pierIdx, face, gdrIdx] = key;
         auto* pPier = bridge_desc.GetPier(pierIdx);
         for (auto f : { pgsTypes::Back, pgsTypes::Ahead })
         {
            if (f == face || pPier->IsAbutment())
               pPier->SetBearingData(gdrIdx, f, bd);
         }
      }
      WBFL::System::Logger::Info(_T("Bearings: the bearings vary by girder. Bearing data is set for each girder."));
   }
}

std::optional<Eigen::Vector2d> get_referent_direction(IfcSchema::IfcBridgePart pier)
{
   // an explicit RefDirection of the linear placement of the referent that positions the pier is the direction of the CL pier.
   // Without one, the placement follows the alignment and says nothing about the pier direction.
   auto rel_positions = pier.PositionedRelativeTo();
   if (rel_positions.empty())
      return std::nullopt;

   auto referent = rel_positions.front().RelatingPositioningElement().as<IfcSchema::IfcReferent>();
   if (!referent)
      return std::nullopt;

   auto linear_placement = referent.ObjectPlacement().as<IfcSchema::IfcLinearPlacement>();
   if (!linear_placement)
      return std::nullopt;

   auto ref_direction = linear_placement.RelativePlacement().RefDirection();
   if (!ref_direction)
      return std::nullopt;

   auto ratios = ref_direction.DirectionRatios();
   Eigen::Vector2d direction(ratios[0], ratios[1]);
   if (direction.norm() < 1.0e-9)
      return std::nullopt; // vertical
   return direction.normalized();
}

namespace
{
   // direction of the best fit line through the points (principal direction), if there are at least two distinct points
   std::optional<Eigen::Vector2d> fit_line_direction(const std::vector<Eigen::Vector2d>& points)
   {
      if (points.size() < 2)
         return std::nullopt;

      Eigen::Vector2d mean = Eigen::Vector2d::Zero();
      for (const auto& p : points) mean += p;
      mean /= (double)points.size();
      Eigen::Matrix2d covariance = Eigen::Matrix2d::Zero();
      for (const auto& p : points) covariance += (p - mean) * (p - mean).transpose();
      if (covariance.trace() < 1.0e-12)
         return std::nullopt;

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> solver(covariance);
      return Eigen::Vector2d(solver.eigenvectors().col(1)); // eigenvector of the largest eigenvalue
   }

   // angle between two lines, 0 to 90 degrees
   Float64 angle_between_lines(const Eigen::Vector2d& a, const Eigen::Vector2d& b)
   {
      return acos(std::clamp(fabs(a.dot(b)), 0.0, 1.0));
   }

   // the ends of the girders at each pier face
   struct GirderEnd
   {
      Eigen::Vector2d end;                     // end of the girder
      std::optional<Eigen::Vector2d> bearing;  // CL bearing
      Eigen::Vector2d inward;                  // direction from the end of the girder into the span
   };
   using PierFaces = std::map<std::pair<PierIndexType, pgsTypes::PierFaceType>, std::vector<GirderEnd>>;

   PierFaces get_pier_faces(const GirderLayout& layout)
   {
      PierFaces pier_faces;
      for (SpanIndexType spanIdx = 0; spanIdx < layout.size(); spanIdx++)
      {
         for (const auto& girder : layout[spanIdx])
         {
            Eigen::Vector2d direction = (girder.end - girder.start).normalized();
            pier_faces[{spanIdx, pgsTypes::Ahead}].push_back({ girder.start, girder.start_bearing, direction });
            pier_faces[{spanIdx + 1, pgsTypes::Back}].push_back({ girder.end, girder.end_bearing, -direction });
         }
      }
      return pier_faces;
   }

   // PGSuper orientation string for a CL pier direction. The skew is measured from the left normal to the
   // alignment (the direction of a NORMAL pier line), positive counterclockwise (left)
   std::_tstring orientation_string(std::shared_ptr<WBFL::EAF::Broker> pBroker, Float64 station, const Eigen::Vector2d& direction)
   {
      GET_IFACE2(pBroker, IRoadway, pAlignment);
      CComPtr<IDirection> normal;
      pAlignment->GetBearingNormal(station, &normal); // normal to the right
      Float64 angle;
      normal->get_Value(&angle);
      Float64 left_normal = angle + M_PI;

      Float64 skew = atan2(direction.y(), direction.x()) - left_normal;
      skew = atan2(sin(skew), cos(skew));
      // the pier line has no direction, so the skew is between -90 and 90 degrees
      if (skew <= -PI_OVER_2)
         skew += M_PI;
      else if (PI_OVER_2 < skew)
         skew -= M_PI;

      Float64 skew_deg = RoundOff(WBFL::Units::ConvertFromSysUnits(skew, WBFL::Units::Measure::Degree), 1.0e-6);
      if (IsZero(skew_deg))
         return _T("NORMAL");

      std::_tostringstream os;
      os << std::setprecision(10) << fabs(skew_deg) << (0 < skew_deg ? _T(" L") : _T(" R"));
      return os.str();
   }
}

std::vector<Eigen::Vector2d> set_pier_orientation(std::shared_ptr<WBFL::EAF::Broker> pBroker, const GirderLayout& layout, const std::vector<Float64>& pier_stations, const std::vector<std::optional<Eigen::Vector2d>>& referent_directions, CBridgeDescription2& bridge_desc)
{
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   auto pier_faces = get_pier_faces(layout);

   const Float64 tolerance = WBFL::Units::ConvertToSysUnits(0.5, WBFL::Units::Measure::Degree);

   std::vector<Eigen::Vector2d> pier_directions;
   for (PierIndexType pierIdx = 0; pierIdx < pier_stations.size(); pierIdx++)
   {
      // directions of the CL bearing lines and the girder end lines on the faces of the pier. Both parallel the CL pier
      std::vector<Eigen::Vector2d> bearing_lines, girder_end_lines;
      for (auto face : { pgsTypes::Back, pgsTypes::Ahead })
      {
         auto found = pier_faces.find({ pierIdx, face });
         if (found == pier_faces.end())
            continue;

         std::vector<Eigen::Vector2d> bearings, ends;
         for (const auto& e : found->second)
         {
            ends.push_back(e.end);
            if (e.bearing)
               bearings.push_back(*e.bearing);
         }
         if (bearings.size() == found->second.size())
         {
            if (auto d = fit_line_direction(bearings)) bearing_lines.push_back(*d);
         }
         if (auto d = fit_line_direction(ends)) girder_end_lines.push_back(*d);
      }

      // the CL pier direction from the first available source
      std::vector<std::pair<LPCTSTR, Eigen::Vector2d>> sources;
      if (referent_directions[pierIdx]) sources.emplace_back(_T("positioning referent"), *referent_directions[pierIdx]);
      if (!bearing_lines.empty()) sources.emplace_back(_T("CL bearing line"), bearing_lines.front());
      if (!girder_end_lines.empty()) sources.emplace_back(_T("girder ends"), girder_end_lines.front());

      Eigen::Vector2d direction;
      std::_tostringstream os;
      os << LABEL_PIER_EX(bridge_desc.GetPier(pierIdx)->IsAbutment(), pierIdx) << _T(" orientation ");
      if (sources.empty())
      {
         CComPtr<IDirection> normal;
         pAlignment->GetBearingNormal(pier_stations[pierIdx], &normal);
         Float64 angle;
         normal->get_Value(&angle);
         direction = Eigen::Vector2d(cos(angle), sin(angle));
         bridge_desc.GetPier(pierIdx)->SetOrientation(_T("NORMAL"));
         os << _T("NORMAL (assumed, no referent direction, bearings, or girder ends to determine it)");
      }
      else
      {
         direction = sources.front().second;
         auto orientation = orientation_string(pBroker, pier_stations[pierIdx], direction);
         bridge_desc.GetPier(pierIdx)->SetOrientation(orientation.c_str());
         os << orientation << _T(" (from ") << sources.front().first << _T(")");

         // cross check with the other sources
         for (size_t i = 1; i < sources.size(); i++)
         {
            Float64 difference = angle_between_lines(direction, sources[i].second);
            if (tolerance < difference)
            {
               os << _T(". The ") << sources[i].first << _T(" differs by ") << WBFL::Units::ConvertFromSysUnits(difference, WBFL::Units::Measure::Degree) << _T(" deg");
            }
         }
      }
      WBFL::System::Logger::Info(os.str().c_str());
      pier_directions.push_back(direction);
   }

   return pier_directions;
}

void set_end_distance_and_bearing_offset(std::shared_ptr<WBFL::EAF::Broker> pBroker, GirderLayout& layout, const std::vector<Float64>& pier_stations, const std::vector<Eigen::Vector2d>& pier_directions, CBridgeDescription2& bridge_desc)
{
   GET_IFACE2(pBroker, IRoadway, pAlignment);

   auto pier_faces = get_pier_faces(layout);

   const Float64 inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);

   for (const auto& [key, ends] : pier_faces)
   {
      auto [pierIdx, face] = key;
      Float64 station = pier_stations[pierIdx];
      auto* pPier = bridge_desc.GetPier(pierIdx);

      // CL pier: through the alignment at the pier station
      CComPtr<IPoint2d> point;
      pAlignment->GetPoint(station, 0.0, nullptr, pgsTypes::PlanCoordinateType::pcGlobal, &point);
      Float64 x, y;
      point->Location(&x, &y);
      Eigen::Vector2d pier_point(x, y);

      // normal to the CL pier, pointing into the span
      const auto& pier_line = pier_directions[pierIdx];
      Eigen::Vector2d normal(-pier_line.y(), pier_line.x());
      if (normal.dot(ends.front().inward) < 0)
         normal = -normal;

      bool bBearings = std::all_of(ends.begin(), ends.end(), [](const auto& e) {return e.bearing.has_value(); });

      std::vector<Float64> end_distances, bearing_offsets;
      if (bBearings)
      {
         for (const auto& e : ends)
         {
            end_distances.push_back((*e.bearing - e.end).dot(normal));
            bearing_offsets.push_back((*e.bearing - pier_point).dot(normal));
         }
      }
      else
      {
         // Without bearings, only the distance from the CL pier to the end of the girder is known. Keep the ends of the
         // girders where they are in the model with the end distance from the template: the bearing offset is the rest.
         auto [template_end_distance, measure] = pPier->GetGirderEndDistance(face, true);
         if (measure != ConnectionLibraryEntry::EndDistanceMeasurementType::FromBearingNormalToPier)
            template_end_distance = pPier->GetGirderEndDistance(face == pgsTypes::Back ? pgsTypes::Ahead : pgsTypes::Back, true).first; // best we can do
         for (const auto& e : ends)
         {
            end_distances.push_back(template_end_distance);
            bearing_offsets.push_back((e.end - pier_point).dot(normal) + template_end_distance);
         }
      }

      Float64 end_distance = median(end_distances);
      Float64 bearing_offset = median(bearing_offsets);

      if (!bBearings)
      {
         // CL bearing on the CL girder, end distance from the end of the girder normal to the CL pier
         SpanIndexType spanIdx = (face == pgsTypes::Ahead ? pierIdx : pierIdx - 1);
         for (auto& girder : layout[spanIdx])
         {
            Eigen::Vector2d direction = (girder.end - girder.start).normalized();
            if (face == pgsTypes::Ahead)
               girder.start_bearing = girder.start + direction * end_distance / direction.dot(normal);
            else
               girder.end_bearing = girder.end - direction * end_distance / (-direction).dot(normal);
         }
      }

      // an abutment has girders on one face only. the other face gets the same values.
      for (auto f : { pgsTypes::Back, pgsTypes::Ahead })
      {
         if (f == face || pPier->IsAbutment())
         {
            pPier->SetGirderEndDistance(f, end_distance, ConnectionLibraryEntry::EndDistanceMeasurementType::FromBearingNormalToPier);
            pPier->SetBearingOffset(f, bearing_offset, ConnectionLibraryEntry::BearingOffsetMeasurementType::NormalToPier);
         }
      }

      std::_tostringstream os;
      os << LABEL_PIER_EX(pPier->IsAbutment(), pierIdx) << (face == pgsTypes::Ahead ? _T(" ahead") : _T(" back")) << _T(" face: end distance ")
         << WBFL::Units::ConvertFromSysUnits(end_distance, WBFL::Units::Measure::Inch) << _T(" in, bearing offset ")
         << WBFL::Units::ConvertFromSysUnits(bearing_offset, WBFL::Units::Measure::Inch) << _T(" in")
         << (bBearings ? _T(" (from the CL bearings)") : _T(" (no bearings: the end distance is from the template and the bearing offset locates the ends of the girders)"));

      auto [min_ed, max_ed] = std::minmax_element(end_distances.begin(), end_distances.end());
      auto [min_bo, max_bo] = std::minmax_element(bearing_offsets.begin(), bearing_offsets.end());
      if (inch / 4 < *max_ed - *min_ed || inch / 4 < *max_bo - *min_bo)
      {
         os << _T(". The girders have different end distances (")
            << WBFL::Units::ConvertFromSysUnits(*min_ed, WBFL::Units::Measure::Inch) << _T(" to ") << WBFL::Units::ConvertFromSysUnits(*max_ed, WBFL::Units::Measure::Inch)
            << _T(" in) or bearing offsets (") << WBFL::Units::ConvertFromSysUnits(*min_bo, WBFL::Units::Measure::Inch) << _T(" to ") << WBFL::Units::ConvertFromSysUnits(*max_bo, WBFL::Units::Measure::Inch)
            << _T(" in). The median is used");
      }
      WBFL::System::Logger::Info(os.str().c_str());
   }
}
