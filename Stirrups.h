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

#include "GirderBars.h"

// Imports the stirrups of the girders (G3 of devdocs/GirderImportPlan.md): the primary shear zones of CShearData2.
//
// A stirrup is a bar of a girder (GirderBars.h) with a vertical leg: a straight piece within 15 degrees of vertical at
// least half the girder's depth long. Bars within 1 in of each other along the girder are at one station (e.g. a pair of
// single-leg bars), whose number of legs is the sum of theirs; a leg whose top is more than 1 in above the girder's top
// engages the deck (a horizontal interface shear bar).
//
// The zones follow from the stations: the space before each station is a spacing, starting from the girder's start, and
// consecutive spaces that are equal (within 1/4 in) with the same bar size, number of legs, and legs engaging the deck
// make a zone. The zones are symmetric when every station has one at the same distance from the other end (within 1 in)
// with the same bars: then they're made from the stations in the first half, and the last zone goes to midspan, with the
// odd space at the center left out. Zones that are physically the same (e.g. two zones of one spacing) come out as one.
// The spacings are rounded to 1/8 in (1 mm with SI display units), as a designer gives them, and a zone is its number of
// spaces at that spacing long (zone lengths rounded on their own would move the zones after them, e.g. PGSuper's 1.5 in
// first zone).
//
// Confinement bars give the zones' confinement bar size (PGSuper's confinement bars are at the stirrups, at the zone's
// spacing). A confinement bar is a bar of the girder that isn't a stirrup, lies below 40% of the girder's depth, and whose
// path in the girder's section, its ends joined, encloses the bottom row of strands where it is (strands within 1 in of
// the lowest; without strands in the model, the CL girder 1 in above the bar's lowest point). Nothing in it is
// agency-specific: bars around the bottom flange (e.g. Iowa's IA56 bars, the exporter's G9 bars) and closed bottom flange
// ties (Iowa's S10 bars, the exporter's G10 bars) both enclose the strands. A confinement bar is at the station nearest to
// it within 3 in and half the spacing there; the station's confinement bar size is the bar that encloses the most (the
// others are reported), and bars that aren't at a station are reported and skipped. The confinement bar size is part of a
// station's bars, so a zone ends where the confinement starts or stops.
//
// The bar material (one per girder) is the most common of the stirrups', else the template's. Each girder gets exactly
// the model's stirrups, none if it has none: the girder library entry's transverse reinforcement (and the template's,
// which has it) is a seed value the user copies into a girder. Splitting reinforcement isn't imported yet; the girders
// have none.
void import_stirrups(std::shared_ptr<WBFL::EAF::Broker> pBroker, const GirderBars& girder_bars, CBridgeDescription2& bridge_desc);
