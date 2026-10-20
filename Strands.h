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

#include "BeamSpacing.h"
#include <psgLib/BridgeDescription2.h>

// Imports the prestressing strands of the girders (G2 of devdocs/GirderImportPlan.md).
//
// The strands are the elements of the tendon role of the mapping table. A strand belongs to the girder it's a part of
// (through its assemblies), else to the girder whose geometry contains its mid-point. Its centerline is read from its
// body geometry, following mapped items: the straight pieces (circular profiles extruded along a line, swept disks along
// polylines) chained end to end, with revolved solids as the bends between them. The centerline is clipped to the girder's
// plan footprint, and tendons that don't run the length of the girder (e.g. lifting loops) are skipped. A straight strand
// that projects beyond an end of the girder (e.g. bent up into the diaphragm) is an extended strand at that end, unless it's
// debonded there; PGSuper keeps only that it's extended, not the shape of the extension, and extends only straight strands. The centerline is put in
// girder coordinates: the distance from the start of the
// girder in plan, the offset from the CL girder (positive to the right looking ahead on station), and the depth below the
// top of the girder, measured normal to the girder's top line from the CL girder at the ends.
//
// A strand is temporary if its tendon.status is TEMPORARY, harped if its centerline bends, and straight otherwise. When
// every strand is at a position of the girder library entry's strand grid (harped strands with a common vertical offset at
// each end and harp point), the grid is filled (pgsTypes::sdtDirectSelection) with the harped strand offsets, also when the
// import created the entry with a generated grid. Otherwise the strands are input as individual strands
// (pgsTypes::sdtDirectStrandInput) with the harped strands' common harp points.
// Both have the debonding of each strand (tendon.debond_start and tendon.debond_end), the strand material of each strand type (from
// tendon.diameter, else the type's nominal diameter, else the geometry, and tendon.fpu, tendon.grade, and
// tendon.corrosion_treatment), and the jacking force (IfcTendon.PreStress, else IfcTendon.TensionForce, else
// tendon.tension_force). Girders without strands keep the template's strands.
void import_strands(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file, const GirderLayout& layout, CBridgeDescription2& bridge_desc);
