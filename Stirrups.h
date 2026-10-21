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
//
// The bar material (one per girder) is the most common of the stirrups', else the template's. Each girder gets exactly
// the model's stirrups, none if it has none: the girder library entry's transverse reinforcement (and the template's,
// which has it) is a seed value the user copies into a girder. Confinement and splitting reinforcement aren't imported
// yet; the girders have none.
void import_stirrups(const GirderBars& girder_bars, CBridgeDescription2& bridge_desc);
