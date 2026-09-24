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

#include <psgLib/DeckPoint.h>

#include "BeamSpacing.h"
#include <psgLib/BridgeDescription2.h>

// Geometry of the deck slab (IfcSlab.FLOOR), or nullopt if it can't be processed
std::optional<Mesh> get_deck_mesh(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file);

// Geometry of the deck slab and the haunches that are separate elements, in one mesh: parts aggregated under the deck,
// and slabs and building element parts identified as haunches by ObjectType or Name. Returns nullopt if there are no
// separate haunches; the deck slab geometry (get_deck_mesh) then has all of the deck concrete, including haunches that
// are part of the deck solid or other items of its representation.
std::optional<Mesh> get_deck_concrete_mesh(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file);

// Deck edge points, relative to the alignment, that reproduce the edges of the deck slab geometry
std::vector<CDeckPoint> get_deck_edge_points(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck);

// Cast-in-place deck dimensions measured from the deck and girder geometry.
// Values that couldn't be measured are nullopt
struct DeckSectionData
{
   std::optional<Float64> gross_depth;
   std::array<std::optional<Float64>, 2> edge_depth; // pgsTypes::stLeft, stRight
   std::array<std::optional<pgsTypes::DeckOverhangTaper>, 2> overhang_taper;
   std::optional<pgsTypes::HaunchShapeType> haunch_shape; // only if the haunches are part of the deck solid
   std::optional<Float64> fillet;
   std::vector<std::vector<std::array<std::optional<Float64>, 2>>> slab_offset; // [span][girder][pgsTypes::metStart/metEnd], at CL bearing
   bool slab_offset_extrapolated = false; // true if the deck doesn't extend over some bearings and the slab offset was extrapolated
};

// Measures the deck across lines normal to the alignment at the quarter points of each span and
// the slab offset ("A" dimension) at the CL bearings. The CL bearings are located with the girder
// end distances in bridge_desc.
DeckSectionData analyze_deck_section(std::shared_ptr<WBFL::EAF::Broker> pBroker, const Mesh& deck, const GirderLayout& layout, const std::vector<Float64>& pier_stations, const CBridgeDescription2& bridge_desc);

bool create_alignment_from_deck(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file);