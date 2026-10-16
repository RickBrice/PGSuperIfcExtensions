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

#include "stdafx.h"
#include "BarShape.h"

#include <Materials/Rebar.h>

#include <cmath>
#include <iomanip>
#include <sstream>

namespace
{
   using Point = BarPath::Point;

   Point minus(const Point& a, const Point& b) { return { a[0] - b[0], a[1] - b[1], a[2] - b[2] }; }
   Float64 dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
   Float64 norm(const Point& a) { return sqrt(dot(a, a)); }

   // length of the circular arc from a through b to c
   Float64 arc_length(const Point& a, const Point& b, const Point& c)
   {
      Float64 ab = norm(minus(b, a)), bc = norm(minus(c, b)), ac = norm(minus(c, a));
      if (ab == 0 || bc == 0)
         return ac;

      // the inscribed angle at b subtends the arc from a to c that doesn't contain b, so the arc through b
      // has the central angle 2*pi - 2*angle_b
      Float64 cos_b = std::clamp(dot(minus(a, b), minus(c, b)) / (ab * bc), -1.0, 1.0);
      Float64 angle_b = acos(cos_b);
      Float64 sin_b = sin(angle_b);
      if (IsZero(sin_b, 1e-12))
         return ab + bc; // collinear points: a straight line
      Float64 radius = ac / (2 * sin_b);
      return radius * (TWO_PI - 2 * angle_b);
   }
}

BarPath BarPath::Straight(const Point& start, const Point& end)
{
   BarPath path;
   path.AddLine(path.AddPoint(start), path.AddPoint(end));
   return path;
}

IndexType BarPath::AddPoint(const Point& p)
{
   m_Points.push_back(p);
   return m_Points.size() - 1;
}

void BarPath::AddLine(IndexType a, IndexType b)
{
   m_Segments.push_back({ a, b });
}

void BarPath::AddArc(IndexType a, IndexType b, IndexType c)
{
   m_Segments.push_back({ a, b, c });
}

Float64 BarPath::GetLength() const
{
   Float64 length = 0;
   for (const auto& s : m_Segments)
   {
      if (s.IsArc())
         length += arc_length(m_Points[s.a], m_Points[s.b], m_Points[s.c]);
      else
         length += norm(minus(m_Points[s.b], m_Points[s.a]));
   }
   return length;
}

BarPath BarPath::Scaled(Float64 sx, Float64 sy, Float64 sz) const
{
   BarPath path(*this);
   for (auto& p : path.m_Points)
      p = { p[0] * sx, p[1] * sy, p[2] * sz };
   return path;
}

std::string BarPath::GetSignature() const
{
   std::ostringstream os;
   os << std::fixed << std::setprecision(5);
   for (const auto& p : m_Points)
      os << "(" << p[0] << "," << p[1] << "," << p[2] << ")";
   for (const auto& s : m_Segments)
      os << (s.IsArc() ? "a" : "l") << s.a << "," << s.b << (s.IsArc() ? "," + std::to_string(s.c) : std::string());
   return os.str();
}

Float64 GetBarMass(const WBFL::Materials::Rebar* pRebar, Float64 length)
{
   static const Float64 density = WBFL::Units::ConvertToSysUnits(490.0, WBFL::Units::Measure::LbmPerFeet3);
   Float64 db = pRebar->GetNominalDimension();
   return length * M_PI * db * db / 4 * density;
}
