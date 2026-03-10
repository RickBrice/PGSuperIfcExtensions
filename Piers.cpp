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
#include "Piers.h"
#include "Properties.h"
#include "Geometry.h"

#include <psgLib\GirderLabel.h>
#include <IFace/Alignment.h>

#include <limits>

#include <ifcgeom/abstract_mapping.h>
#include <ifcgeom/iterator.h>
#include <ifcgeom/ifcgeomelement.h>
#include <ifcgeom/kernels/opencascade/OpenCascadeKernel.h>


PierIndexType get_pier_count(IfcParse::IfcFile& file)
{
   auto parts = file.instances_by_type<IfcSchema::IfcBridgePart>();
   PierIndexType nPiers = 0;
   for (auto& part : *parts)
   {
      if (part->PredefinedType().has_value() && (part->PredefinedType().get() == IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_ABUTMENT || part->PredefinedType().get() == IfcSchema::IfcBridgePartTypeEnum::IfcBridgePartType_PIER))
         nPiers++;
   }
   return nPiers;
}

std::set<int> get_bearing_ids(IfcSchema::IfcBridgePart* pier)
{
   std::set<int> bearing_ids;
   auto rels = pier->ContainsElements();
   for (auto& rel : *rels)
   {
      auto related_elements = rel->RelatedElements();
      for (auto& related_element : *related_elements)
      {
         auto bearing = related_element->as<IfcSchema::IfcBearing>();
         if (bearing)
         {
            bearing_ids.insert(bearing->id());
         }
      }
   }

   return bearing_ids;
}

double get_pier_station(std::shared_ptr<WBFL::EAF::Broker> pBroker,IfcParse::IfcFile& file, PierIndexType pierIdx, IfcSchema::IfcBridgePart* pier)
{
   auto rel_positions = pier->PositionedRelativeTo();
   auto positioning_element = (rel_positions && 0 < rel_positions->size() ? (*rel_positions->begin())->RelatingPositioningElement() : nullptr);
   auto referent = (positioning_element ? positioning_element->as<IfcSchema::IfcReferent>() : nullptr);
   auto station = (referent ? GetProperty<IfcSchema, IfcSchema::IfcLengthMeasure>(referent, "Pset_Stationing", "Station") : nullptr);
   if (referent && station)
   {
      return *station;
   }
   else
   {
      USES_CONVERSION;
      std::ostringstream os;
      os << "Expected Pier " << T2A(LABEL_PIER(pierIdx)) << " ";
      pier->toString(os);
      os << " to be positioned with an IfcReferent and have stationing defined with Pset_Stationing. Attempting to estimate station from all IfcBearing in the spatial structure of the IfcBridgePart.PIER";
      WBFL::System::Logger::Info(os.str().c_str());

      ifcopenshell::geometry::Settings settings;
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
         IfcGeom::instance_id_filter filter(true, false, bearing_ids);
         std::vector<IfcGeom::filter_t> filters({ filter });

         std::unique_ptr<IfcGeom::OpenCascadeKernel> kernel(std::make_unique<IfcGeom::OpenCascadeKernel>(settings));
         int num_threads = 1;// std::thread::hardware_concurrency();
         IfcGeom::Iterator iterator(std::move(kernel), settings, &file, filters, num_threads);


         Eigen::Vector3d vmax(-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity());
         Eigen::Vector3d vmin(std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
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
