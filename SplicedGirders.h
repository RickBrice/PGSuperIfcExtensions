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

// Spliced girders (PGSplice): girder lines of precast segments and closure joints, the girder groups, segment
// connections at the piers, and temporary supports they make, and the layouts the rest of the import works with
// (devdocs/PGSplicePlan.md, S3)

#pragma once

#include "BeamSpacing.h"
#include <PsgLib/BridgeDescription2.h>
#include <map>
#include <optional>
#include <vector>

class GirderProfileReader;

// A girder line of the model: an IfcElementAssembly (the mapping table's girder_assembly role) of precast segments (its
// beams that aren't closure joints) and the closure joints (closure_joint) between them, or a lone girder (one segment)
struct SplicedGirder
{
   int id = 0;                                   // the IfcElementAssembly, or the IfcBeam of a lone girder
   std::string name;
   std::vector<GirderCenterline> segments;       // in order by station
   std::vector<std::optional<int>> closure_joints; // between segments i and i+1: the closure joint IfcBeam, or nullopt
};

// The number of spliced girders: girder assemblies that aggregate more than one beam, or that have a closure joint, and
// closure joints outside an assembly
IndexType count_spliced_girders(ifcopenshell::file& file);

// The girder lines: the girder assemblies, with their segments located from the geometry (get_beam_centerlines) and their
// closure joints put in the gaps between segments by where they are, and the lone girders (get_beam_ids)
std::vector<SplicedGirder> get_spliced_girders(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file);

// How the girders are framed: the girder groups, the segment connections at the piers inside the groups, and the
// temporary supports
struct SplicedFraming
{
   std::vector<bool> group_boundary; // by pier: the girders end there (abutments always)
   std::map<PierIndexType, pgsTypes::PierSegmentConnectionType> connections; // the piers inside a group

   struct TemporarySupport
   {
      Float64 station;
      pgsTypes::TemporarySupportType type;
      pgsTypes::TempSupportSegmentConnectionType connection;
      Float64 closure_length = 0; // along the girder, for a closure joint
   };
   std::vector<TemporarySupport> temporary_supports; // in order by station

   GroupIndexType GetGroup(Float64 station, const std::vector<Float64>& pier_stations) const; // the group a station is in
   GroupIndexType GetGroupCount() const;
};

// The framing from the girder lines and the erection towers (the table's erection_tower role). At each pier inside the
// bridge: a segment over it is a continuous segment, a closure joint over it is a closure joint connection, and girders
// ending at it make it a group boundary. A closure joint between piers is at a temporary support: an erection tower if
// one is there, else a strongback (the user, 2026-11-10). An erection tower under a segment is a temporary support with a
// continuous segment. The leftmost girder decides; girders that disagree are logged
SplicedFraming get_spliced_framing(ifcopenshell::file& file, const std::vector<SplicedGirder>& girders, const std::vector<Float64>& pier_stations);

// The timeline events and segment end geometry for new temporary supports and closure joints, from the template's temporary
// support (before ReduceToSingleSegmentGirders removes it), else from its timeline: erect with the piers, cast the closure
// joints and remove the temporary supports with the deck
struct TemporarySupportDefaults
{
   EventIndexType erection = INVALID_INDEX;
   EventIndexType removal = INVALID_INDEX;
   EventIndexType cast_closure = INVALID_INDEX;
   Float64 end_distance = 0; // from the CL bearing to the end of the segment, along the girder
};
TemporarySupportDefaults get_temporary_support_defaults(const CBridgeDescription2& bridge_desc);

// Sets the segment connections at the piers inside the groups and adds the temporary supports (the bridge description
// has its spans, groups, pier stations, and girder counts)
void apply_spliced_framing(const SplicedFraming& framing, const TemporarySupportDefaults& defaults, CBridgeDescription2& bridge_desc);

// The layouts the rest of the import works with. Each sets the girder keys and segment indices of the girders' segments
// (girder index from left to right in a group).
// The span layout [span][girder]: the piece of each girder line in each span, for the pier, spacing, bearing, and deck
// functions. A piece ends where a segment ends at the pier (a group boundary, or a closure joint over the pier), or at the
// pier on a segment that continues over it (continuous[end])
GirderLayout get_spliced_span_layout(std::vector<SplicedGirder>& girders, const SplicedFraming& framing, const std::vector<Float64>& pier_stations);

// The segment layout [group][segments of the group's girders], for the strands, bars, materials, and depth variation
GirderLayout get_spliced_segment_layout(const std::vector<SplicedGirder>& girders, const SplicedFraming& framing);

// One segment of each girder [group][girder] for the girder's library entry: the segment closest to constant depth
GirderLayout get_spliced_girder_layout(ifcopenshell::file& file, const GirderProfileReader& reader, const std::vector<SplicedGirder>& girders, const SplicedFraming& framing);

// Sets each segment's depth variation (none, linear, parabolic, double linear, double parabolic) from the heights and
// bottom flange depths of its cross sections (GirderProfileReader::ReadCrossSections). PGSplice measures the zones over
// the segment's framing length, from the CL of the closure joint (or the pier at the end of a group), so the prismatic
// zones at the ends include the distance from there to the end of the segment
void set_segment_variations(ifcopenshell::file& file, const GirderProfileReader& reader, const GirderLayout& segment_layout, const std::vector<Float64>& pier_stations, CBridgeDescription2& bridge_desc);
