///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2025  Washington State Department of Transportation
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

#include "IfcExporter.h"
#include <IFace\Project.h>

class CIfcImporter;

class CIfcAlignmentImporter
{
public:
   CIfcAlignmentImporter(CIfcImporter& importer);
   bool Import(IfcParse::IfcFile& file);

private:
   CIfcImporter& m_Importer;
   AlignmentData2 m_AlignmentData;
   ProfileData2   m_ProfileData;
   RoadwaySectionData m_RoadwaySectionData;
   bool m_bAlignmentStarted;
   int m_ProfileState; // -1 = not yet started, 0 = started, but grade not determined, 1 = first point established
   CComPtr<ICogoEngine> m_CogoEngine;
   CComPtr<IGeomUtil2d> m_GeomUtil;

   bool InitAlignmentParameters(IfcParse::IfcFile& file);
   Ifc4x3_add2::IfcAlignment* GetAlignment(IfcParse::IfcFile& file);
   Float64 LoadAlignment(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment);
   void LoadProfile(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment, Float64 stationAdjustment);
   bool IsValidAlignment(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment);
   void GetStations(Ifc4x3_add2::IfcAlignment* pAlignment, std::vector<std::pair<Float64, Float64>>& vStations, std::vector<std::tuple<Float64, Float64, Float64>>& vStationEquations);
   Float64 GetStartStation(Ifc4x3_add2::IfcAlignment* pAlignment);

   enum LastAlignmentType { Unknown, Line, Curve } m_LastAlignmentType;

   // adds a line to the alignment. returns the station at the end of the line
   Float64 OnLine(Float64 startStation, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pLine);
   Float64 OnLine(Float64 sx, Float64 sy, Float64 startStation, Float64 startDirection, Float64 length);

   // adds a curve to the alignment. returns the station at the end of the curve
   Float64 OnCurve(Float64 startStation, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pEntrySpiral, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pCurve, Ifc4x3_add2::IfcAlignmentHorizontalSegment* pExitSpiral);

   // adds linear segment to the profile
   void OnLinearSegment(Float64 startStation, Ifc4x3_add2::IfcAlignmentVerticalSegment* pLinearSegment);

   // adds a parabolic curve to the profile
   void OnParabolicSegment(Float64 startStation, Ifc4x3_add2::IfcAlignmentVerticalSegment* pParaCurve);

   void GetCurvePoints(Ifc4x3_add2::IfcAlignmentHorizontalSegment* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter);

   void GetSpiralPoints(Ifc4x3_add2::IfcAlignmentHorizontalSegment* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd);

   void CheckSpiralType(Ifc4x3_add2::IfcAlignmentHorizontalSegment* pSpiral);

   void GetPoint(Ifc4x3_add2::IfcCartesianPoint* pPoint, Float64* pX, Float64* pY);
};