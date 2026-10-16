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

namespace WBFL { namespace EAF { class Broker; }; };

#include "Geometry.h"
#include <psgLib/BridgeDescription2.h>

// An element's geometry
struct ElementMesh
{
   IfcSchema::IfcObject element;
   Mesh mesh;
};

// The geometry of the elements of the barrier role of the mapping table (empty if there are none)
std::vector<ElementMesh> get_barrier_meshes(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file);

// The vertical plane normal to the alignment at a station, for CutMesh: loop points are (distance to the left of the
// alignment, elevation)
struct AlignmentCut
{
   Float64 station;
   Eigen::Vector2d point; // on the alignment
   Eigen::Vector2d axis;  // the alignment's direction
};
AlignmentCut get_alignment_cut(std::shared_ptr<WBFL::EAF::Broker> pBroker, Float64 station);

// Imports the exterior traffic barriers (devdocs/IfcImportPlan.md, Phase 7). The barriers are cut at the middle of the
// bridge, normal to the alignment, with the deck. A barrier's side is the side of the alignment it's on, and its section
// is put in traffic barrier coordinates: from the edge of the deck at the top of the deck, positive toward the roadway
// (the right barrier is a mirror of its library entry), turned back by the roadway cross slope there (PGSuper turns
// barriers to it, so the roadway section must be imported first). The section selects the traffic barrier library entry
// with the same section, or a new entry is created from it. A side with no barrier keeps its railing system
void import_traffic_barriers(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const std::vector<ElementMesh>& barriers, const std::vector<Float64>& pier_stations, CBridgeDescription2& bridge_desc);
