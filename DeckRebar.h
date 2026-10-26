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

#include "Geometry.h"

namespace WBFL { namespace EAF { class Broker; }; };
class CBridgeDescription2;

// Imports the deck reinforcement (CDeckRebarData) as lump sums, the area of bars per unit width of deck. PGSuper's deck
// reinforcement is simplified: each mat has primary bars that run the length of the bridge and supplemental bars at the
// piers, and only their area per width, their depth, and the length of the supplemental bars matter. Lump sums are what
// a model's bars give exactly, whatever their sizes and spacings; bar size and spacing input waits for a more detailed
// deck reinforcement model in PGSuper.
//
// The deck bars are the elements of the mapping table's rebar role whose centerlines are inside the deck slab (90% of
// their length). Longitudinal bars (more along the alignment than across it) are used; transverse bars aren't PGSuper data
// and are counted in the log. Bars at the same offset and depth (within 2 in and 1 in) are a line, e.g. bars lapped end to
// end. Lines that cover at least 75% of the bridge length are the primary reinforcement; the others are supplemental.
//
// Mats: the depths of the primary lines below the top of the deck, split at their largest gap (more than 1 in) into the top
// and bottom mats; with one group, it's the top mat if it's in the upper half of the deck. Every other bar is in the mat
// whose primary bars are nearer in depth.
//
// For each mat, the lump sum is the area of one bar of each primary line (a lap splice doesn't add) over the width the
// lines cover (from the first line to the last, plus one spacing). PGSuper puts a lump sum at its mat's cover, so the
// cover is the depth of the primary bars' centroid below the top of the deck (top mat) or its height above the bottom
// (bottom mat). The material is the one most of the mat's bars have, else the existing one.
//
// Supplemental bars belong to the pier they cross nearest to their middle (bars that cross no pier are counted and
// skipped). For each pier and mat, the lump sum is the area of one bar of each of their lines over the width they cover,
// and the cutoffs are the median distances from the CL pier to the bars' ends.
//
// A model without longitudinal deck bars keeps the existing deck reinforcement (the template's), and the log says so.
void import_deck_rebar(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const Mesh& deck, const std::vector<Float64>& pier_stations, CBridgeDescription2& bridge_desc);
