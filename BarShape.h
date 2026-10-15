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
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
// the Alternate Route Open Source License for more details.
//
// You should have received a copy of the Alternate Route Open Source
// License along with this program; if not, write to the Washington
// State Department of Transportation, Bridge and Structures Office,
// P.O. Box  47340, Olympia, WA 98503, USA or e-mail
// Bridge_Support@wsdot.wa.gov
///////////////////////////////////////////////////////////////////////

#pragma once

// Reinforcing bar shapes, lengths, and masses.
//
// A bar is described by the centerline path of its shape, in the coordinates of its bar type. The bar type's
// geometry, BarLength, and BarMass all come from the path, so a new or more detailed bar shape (e.g. when
// PGSuper's rebar model describes bar bends) only needs its path: the IfcReinforcingBarType and its properties
// follow without other changes (see CreateBarType in Rebar.h).

#include <array>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace WBFL { namespace Materials { class Rebar; } }

// The centerline of a bar: points joined by straight lines and circular arcs, as an IfcIndexedPolyCurve
class BarPath
{
public:
   using Point = std::array<Float64, 3>;

   // a straight line from points[a] to points[b], or a circular arc from points[a] through points[b] to points[c]
   struct Segment
   {
      IndexType a, b, c = INVALID_INDEX;
      bool IsArc() const { return c != INVALID_INDEX; }
   };

   // a straight bar
   static BarPath Straight(const Point& start, const Point& end);

   // adds a point and returns its index
   IndexType AddPoint(const Point& p);
   void AddLine(IndexType a, IndexType b);
   void AddArc(IndexType a, IndexType b, IndexType c);

   const std::vector<Point>& GetPoints() const { return m_Points; }
   const std::vector<Segment>& GetSegments() const { return m_Segments; }

   // the centerline length: the straight lines plus the arcs
   Float64 GetLength() const;

   // the path with every point scaled (e.g. across the girder for a bar that isn't normal to the girder axis).
   // An arc stays the circular arc through its scaled points
   BarPath Scaled(Float64 sx, Float64 sy, Float64 sz) const;

   // text that's the same for two paths exactly when they have the same geometry (to 0.01 mm)
   std::string GetSignature() const;

private:
   std::vector<Point> m_Points;
   std::vector<Segment> m_Segments;
};

// The mass of a bar of the given length. The ASTM A615/A706 (AASHTO M31) nominal masses are based on the nominal
// diameter and 490 lb/ft^3, so the mass is computed the same way and matches the tables (e.g. #4 is 0.668 lb/ft).
// The nominal areas in the tables are rounded, so they aren't used
Float64 GetBarMass(const WBFL::Materials::Rebar* pRebar, Float64 length);

// The detailing of a bar: everything its bar type is made from. Bars with the same definition share a type
struct BarDefinition
{
   std::string name;       // name of the bar type, e.g. "G2 bars - WF100G"
   std::string mark;       // bar mark, e.g. "G2"
   const WBFL::Materials::Rebar* rebar = nullptr;
   BarPath path;           // centerline, in the bar type's coordinates

   // usBrPset_ACI_ReinforcingBar
   std::string element;    // e.g. "BEAM"
   std::string use;        // e.g. "STIRRUP"
   std::string position;   // e.g. "CENTER"

   // usBrPset_ReinforcingCover
   std::optional<Float64> cover_top, cover_side, cover_bottom, cover_end;

   // usBrPset_ACI_BarShape: the ACI bend shape and its dimensions
   std::string shape_code; // e.g. "STRAIGHT", "S11"
   Float64 bend_radius = 0.0;
   std::unordered_map<std::string, Float64> shape_dimensions;

   Float64 GetLength() const { return path.GetLength(); }
   Float64 GetMass() const { return GetBarMass(rebar, GetLength()); }
};
