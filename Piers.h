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

#include "BeamSpacing.h"

class CBridgeDescription2;

PierIndexType get_pier_count(ifcopenshell::file& file);

double get_pier_station(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, PierIndexType pierIdx, IfcSchema::IfcBridgePart pier);

// Locates the CL bearing at each end of each girder from the IfcBearing geometry. Each bearing belongs to the girder end
// whose CL is nearest to it. The CL bearing is the center of the bearings under the end of the girder, projected onto the CL girder.
void locate_bearings(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, GirderLayout& layout);

// Sets the bearing data (shape, size, number and spacing of bearings, fixity) from the IfcBearing geometry and properties.
// The bearing data is for the whole bridge, each pier face, or each girder, depending on how much the bearings vary.
void set_bearing_data(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const GirderLayout& layout, CBridgeDescription2& bridge_desc);

// Direction of the CL pier from an explicit RefDirection of the linear placement of the IfcReferent that positions the pier
std::optional<Eigen::Vector2d> get_referent_direction(IfcSchema::IfcBridgePart pier);

// Sets the orientation of each pier and returns the direction of each CL pier. The direction is from the first of
// the positioning referent, the CL bearing lines, and the line through the girder ends. The other sources are cross checked.
// Without any of them, the pier is normal to the alignment.
std::vector<Eigen::Vector2d> set_pier_orientation(std::shared_ptr<WBFL::EAF::Broker> pBroker, const GirderLayout& layout, const std::vector<Float64>& pier_stations, const std::vector<std::optional<Eigen::Vector2d>>& referent_directions, CBridgeDescription2& bridge_desc);

// Sets the girder end distance and bearing offset for each pier face, measured normal to the CL pier, from the CL bearings.
// Without bearings, the end distance is kept and the bearing offset locates the ends of the girders. The CL bearing
// of each girder is then located from the end distance.
void set_end_distance_and_bearing_offset(std::shared_ptr<WBFL::EAF::Broker> pBroker, GirderLayout& layout, const std::vector<Float64>& pier_stations, const std::vector<Eigen::Vector2d>& pier_directions, CBridgeDescription2& bridge_desc);
