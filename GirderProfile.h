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

// Girder sections read from the swept profiles of a girder's body representation (devdocs/GirderImportPlan.md, G1
// step 1), exact where a section cut from the triangulated geometry (CutGirderSection) is faceted

#include "BeamSpacing.h"
#include "GirderSection.h"

namespace ifcopenshell { namespace geom { class abstract_mapping; class settings; } }

class GirderProfileReader
{
public:
   GirderProfileReader(ifcopenshell::file& file);
   ~GirderProfileReader();

   // The section of a girder at mid-length from the profile of its body's swept solid: an extrusion (tapered ones
   // interpolated), or a sectioned solid (the cross sections before and after mid-length interpolated). Mapped items
   // and clipping results are followed. The profile is projected along the sweep onto the vertical plane normal to the
   // girder axis in plan, the plane CutGirderSection cuts, so a profile in a skewed plane (e.g. at a skewed end) gives
   // the girder's normal section. Returns nullopt, with the reason in why_not, if the body isn't one swept solid along
   // the girder or its profile can't be read
   std::optional<SectionOutline> Read(IfcSchema::IfcBeam beam, const GirderCenterline& girder, std::string& source, std::string& why_not) const;

private:
   std::unique_ptr<ifcopenshell::geom::settings> m_Settings;
   std::unique_ptr<ifcopenshell::geom::abstract_mapping> m_Mapping;
};
