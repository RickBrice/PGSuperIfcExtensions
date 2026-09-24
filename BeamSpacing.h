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

#pragma once

#include <Eigen/Dense>
#include "Geometry.h"
#include <psgLib/Keys.h>

namespace WBFL { namespace EAF { class Broker; }; };

using Spacing = std::map<GroupIndexType, std::vector<double>>;
std::set<int> get_beam_ids(ifcopenshell::file& file);
// An IfcBearing: plan outline, center, and elevation range (global coordinates)
struct BearingGeometry
{
   int id;                           // IfcBearing
   Eigen::Vector2d center;           // center of the plan bounding box
   std::vector<Eigen::Vector2d> plan; // vertices in plan
   Float64 zmin, zmax;
};

// The centerline of a girder in plan, from the girder geometry
struct GirderCenterline
{
   int beam_id;                     // IfcBeam
   CGirderKey girder_key;           // group = span, girder index from left to right looking ahead on station
   Eigen::Vector2d start, end;      // CL girder at the ends of the girder (global coordinates), start is at the lesser station
   Float64 start_station, start_offset;
   Float64 end_station, end_offset;
   std::shared_ptr<Mesh> mesh;      // girder geometry (global coordinates)
   std::optional<Eigen::Vector2d> start_bearing, end_bearing; // CL bearing on the CL girder at the start and end, from the IfcBearing geometry
   std::array<std::vector<BearingGeometry>, 2> bearings;      // bearings at the start and end (pgsTypes::metStart, metEnd)
};

// [span][girder], girders ordered left to right
using GirderLayout = std::vector<std::vector<GirderCenterline>>;

// Locates the girders from the IfcBeam geometry. Each girder is assigned to the span that contains its
// mid-point and girders are ordered left to right within the span by offset from the alignment.
GirderLayout get_girder_layout(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const std::vector<Float64>& pier_stations);

// Girder spacing at the start (ahead side of the previous pier) and end (back side of the next pier) of each span,
// measured along a line through the alignment at each pier station: normal to the alignment, or in the direction
// given for each pier (e.g. along the CL pier)
std::pair<Spacing,Spacing> get_beam_spacing(std::shared_ptr<WBFL::EAF::Broker> pBroker, const GirderLayout& layout, const std::vector<Float64>& pier_stations, const std::vector<Eigen::Vector2d>* pLineDirections = nullptr);
