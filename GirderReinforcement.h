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

namespace ifcopenshell { namespace geom { class abstract_mapping; } }

// Reading the centerlines of a girder's reinforcement (strands and bars) from the model, in girder coordinates. Shared by
// the strand import (Strands.h) and the bar import (LongitudinalRebar.h).

// positions closer than this are the same (1/8 in)
extern const Float64 reinforcement_position_tolerance;

// end distances farther than this aren't measured (2 ft)
extern const Float64 end_distance_limit;

// A point of a centerline in girder coordinates: x is the distance from the start of the girder in plan, z the offset
// from the CL girder (positive to the right looking ahead), and y the depth below the top of the girder
struct GirderPoint
{
   Float64 x, z, y;
};

// the path at x, by linear interpolation along its path and extrapolation of its first and last segments
GirderPoint point_at(const std::vector<GirderPoint>& path, Float64 x);
inline Float64 depth_at(const std::vector<GirderPoint>& path, Float64 x) { return point_at(path, x).y; }

// The centerline of a piece of reinforcement in its girder
struct GirderReinforcement
{
   int beam_id;                             // the girder's IfcBeam
   std::vector<GirderPoint> path;           // inside the girder, from the end at the lesser x, with its ends and bends only
   std::array<bool, 2> projects{ false, false }; // it projects beyond the girder at the start and end (pgsTypes::MemberEndType)
   std::optional<Float64> diameter;         // as drawn (swept disk or circular profile)
   std::array<std::optional<Float64>, 2> end_distance; // from each end to the girder's end face along the centerline (0 if it
                                            // projects), or nullopt if it's farther than end_distance_limit
};

class GirderReinforcementReader
{
public:
   GirderReinforcementReader(ifcopenshell::file& file, const GirderLayout& layout);
   ~GirderReinforcementReader();

   enum class Status
   {
      Read,       // the reinforcement is in a girder
      NoGeometry, // its body geometry has no straight pieces that make a centerline
      NoGirder,   // it isn't part of an imported girder and isn't inside one
      Outside     // its centerline isn't inside its girder's plan footprint
   };

   // The centerline of a strand or bar: the straight pieces of its body geometry (circular profiles extruded along a line,
   // swept disks along polylines and indexed poly curves; mapped items are followed) chained end to end through the bends
   // between them, clipped to the plan footprint of its girder, and put in girder coordinates. Its girder is the one it's
   // part of (through its assemblies), else the one whose geometry contains the point halfway along it, else the one that
   // contains the most of it. An element with more than one representation item whose pieces make a line is a group (e.g.
   // bars of one type written as one element with a mapped item for each bar): each of those items is a strand or bar.
   // Read if at least one is in a girder
   Status Read(IfcSchema::IfcElement element, std::vector<GirderReinforcement>& reinforcement);

   // the imported girders, by IfcBeam id
   const std::map<int, const GirderCenterline*>& GetGirders() const { return m_Girders; }

   // the length of a girder in plan
   Float64 GetLength(int beam_id) const;

   // the depth of a girder at x (top to bottom), from its geometry, or nullopt if it isn't known there
   std::optional<Float64> GetDepth(int beam_id, Float64 x) const;

private:
   Status ReadPiece(std::vector<Eigen::Vector3d> centerline, std::optional<int> beam_id, GirderReinforcement& reinforcement);

   class Profile;
   const Profile& GetProfile(int beam_id) const;
   const Eigen::AlignedBox3d& GetBounds(int beam_id) const; // of a girder's geometry (it must have a mesh)

   std::unique_ptr<ifcopenshell::geom::abstract_mapping> m_Mapping;
   std::map<int, const GirderCenterline*> m_Girders;
   std::set<int> m_BeamIds;
   mutable std::map<int, std::unique_ptr<Profile>> m_Profiles; // made when first needed
   mutable std::map<int, Eigen::AlignedBox3d> m_Bounds;          // made when first needed
};
