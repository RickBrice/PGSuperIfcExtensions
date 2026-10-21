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
#include <psgLib/BridgeDescription2.h>

// Imports the longitudinal reinforcement of the girders (G3 of devdocs/GirderImportPlan.md, CLongitudinalRebarData).
//
// The bars are the elements of the rebar role of the mapping table, assigned to girders and read as strands are (see
// GirderReinforcement.h). A longitudinal bar is a bar whose longest straight piece in the girder runs along it (within 2
// degrees) and is at least 2 ft long. Bars of the same size that lap (collinear within 1.5 bar diameters, overlapping
// along the girder) are one bar, since PGSuper doesn't model laps.
//
// Each bar is referenced to the nearer face of the girder (top or bottom, from the girder's depth in the model), with its
// clear cover. Its layout: full length when both ends are within 3 in of the girder's end faces (measured along the bar,
// so skewed ends are allowed for), from the left or right end when one is, centered on the girder when it's centered, else
// from the nearer end. A bar that projects beyond an end is extended there (PGSuper keeps that a bar is extended, not the
// shape of the extension), when it's full length or starts at that end. Bars with the same face, cover, size, and layout
// make rows: PGSuper centers the bars of a row on the CL girder at a common spacing, so bars that aren't equally spaced
// are split into pairs at +/-z (and a bar at the center); a bar without a mirror bar can't be input and is skipped with a
// warning.
//
// The bar size is from rebar.bar_size (e.g. "#5", "5"), else the bar type's nominal diameter, else the geometry. The bar
// material (one per girder) is from rebar.grade (ASTM A615, A706, or A1035 in the text) and rebar.fy, else the template's.
// Each girder gets exactly the model's longitudinal bars, none if it has none: the bars of the girder library entry (and
// of the template, which has them) are seed values the user copies into a girder, not the girder's reinforcement.
void import_longitudinal_rebar(ifcopenshell::file& file, const GirderLayout& layout, CBridgeDescription2& bridge_desc);
