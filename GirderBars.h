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

#include "GirderReinforcement.h"
#include <psgLib/BridgeDescription2.h>
#include <Materials/Rebar.h>

// The reinforcing bars of the girders, read once for the bar imports (LongitudinalRebar.h, Stirrups.h).

using BarMaterial = std::pair<WBFL::Materials::Rebar::Type, WBFL::Materials::Rebar::Grade>;

// A bar in its girder (one bar of a group of bars written as one element)
struct GirderBar
{
   IfcSchema::IfcReinforcingBar element;
   GirderReinforcement centerline;
   WBFL::Materials::Rebar::Size size = WBFL::Materials::Rebar::Size::bsNone; // bsNone if it isn't known
   std::optional<BarMaterial> material; // nullopt if the model doesn't say
};

class GirderBars
{
public:
   // Reads the elements of the mapping table's rebar role that are in imported girders (see GirderReinforcementReader::Read).
   // The bar size is from rebar.bar_size (e.g. "#5", "5"), else the bar type's nominal diameter, else the geometry. The
   // material is from rebar.grade (ASTM A615, A706, or A1035 in the text) and rebar.fy
   GirderBars(ifcopenshell::file& file, const GirderLayout& layout);

   // the mapping table has a rebar role and the model has elements in it
   bool HasBars() const { return m_bHasBars; }

   // why there are no bars (when HasBars is false)
   const std::string& GetReason() const { return m_Reason; }

   const GirderReinforcementReader& GetGeometry() const { return m_Geometry; }

   // the bars of a girder (IfcBeam id), empty if it has none
   const std::vector<GirderBar>& GetBars(int beam_id) const;

   // The centerlines of a girder's strands (IfcBeam id): the elements of the tendon role that run its length (as in
   // Strands.h), for the bars that enclose them (confinement). Empty if it has none
   const std::vector<std::vector<GirderPoint>>& GetStrands(int beam_id) const;

   // the name of a girder (IfcBeam)
   std::string GetGirderName(int beam_id) const;

   // the segment of an imported girder in the bridge description
   static CPrecastSegmentData* GetSegment(CBridgeDescription2& bridge_desc, const GirderCenterline& girder);

private:
   ifcopenshell::file& m_File;
   GirderReinforcementReader m_Geometry;
   bool m_bHasBars = false;
   std::string m_Reason;
   std::map<int, std::vector<GirderBar>> m_Bars;
   std::map<int, std::vector<std::vector<GirderPoint>>> m_Strands;
};

// the nominal diameter of a bar size (sizes are the same for every type and grade), 0 if the pool doesn't have it
Float64 nominal_diameter(WBFL::Materials::Rebar::Size size);

// the material most of the bars have, or nullopt if no bar has one; bMixed is set if they have different materials
std::optional<BarMaterial> most_common_material(const std::map<BarMaterial, IndexType>& materials, bool& bMixed);

// calls fn for every segment of the bridge
void for_each_segment(CBridgeDescription2& bridge_desc, const std::function<void(CPrecastSegmentData*)>& fn);
