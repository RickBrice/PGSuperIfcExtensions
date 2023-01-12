///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright © 1999-2023  Washington State Department of Transportation
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

#include <IFace\Project.h>

///////////////////////////////////////////////////////////////////////////
// CIfcAlignmentConverter
//
// Converts data between IFC and PGSuper data structures
class CIfcAlignmentConverter
{
public:
   CIfcAlignmentConverter(void);
   ~CIfcAlignmentConverter(void);

   // Converts Ifc data to PGSuper data
   HRESULT ConvertToPGSuper(IBroker* pBroker, CString& strFilePath);

   // Returns a list of notes that were generated during the IFC to PGSuper conversion process
   std::vector<std::_tstring> GetNotes();

   static Float64 GetPrecision() { return m_Precision; }

private:
   static Float64 m_Precision;
   const WBFL::Units::Length* m_pLengthUnit;
   const WBFL::Units::Angle* m_pAngleUnit;
   std::vector<std::_tstring> m_Notes;
   CComPtr<ICogoEngine> m_CogoEngine;
   CComPtr<IGeomUtil2d> m_GeomUtil;

   AlignmentData2 m_AlignmentData;
   ProfileData2   m_ProfileData;
   RoadwaySectionData m_RoadwaySectionData;

   bool m_bAlignmentStarted;
   int m_ProfileState; // -1 = not yet started, 0 = started, but grade not determined, 1 = first point established

   //LX::CrossSects* CreateCrossSections(IBroker* pBroker, LX::IFactory* pFactory);
   //LX::Roadway*    CreateRoadway(IBroker* pBroker, LX::IFactory* pFactory);


   template <typename Schema>
   bool ConvertToPGSuper(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData);

   template <typename Schema>
   bool ConvertToPGSuper_4x3(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData);

   template <typename Schema>
   typename Schema::IfcAlignment* GetAlignment(IfcParse::IfcFile& file);

   template <typename Schema>
   typename Schema::IfcAlignment* GetAlignment_4x3(IfcParse::IfcFile& file);

   template <typename Schema>
   void InitUnits(IfcParse::IfcFile& file);

   Float64 GetStartDistAlong(Ifc4x3_rc2::IfcAlignmentHorizontal* pHorizontal);
   Float64 GetStartDistAlong(Ifc4x3_rc3::IfcAlignmentHorizontal* pHorizontal);
   Float64 GetStartDistAlong(Ifc4x3_rc4::IfcAlignmentHorizontal* pHorizontal);

   template <typename Schema>
   void GetStations(typename Schema::IfcAlignment* pAlignment, std::vector<std::pair<Float64, Float64>>& vStations, std::vector<std::tuple<Float64, Float64, Float64>>& vStationEquations);

   template <typename Schema>
   void LoadAlignment(typename Schema::IfcAlignment* pAlignment);

   template <typename Schema>
   void LoadProfile(typename Schema::IfcAlignment* pAlignment);

   template <typename Schema>
   Float64 LoadAlignment_4x3(IfcParse::IfcFile& file,typename Schema::IfcAlignment* pAlignment);

   template <typename Schema>
   void LoadProfile_4x3(IfcParse::IfcFile& file, typename Schema::IfcAlignment* pAlignment,Float64 stationAdjustment);

   //void LoadCrossSections(LX::CrossSects* pCrossSects, LX::String& strSurfaceName);

   template <typename Schema>
   void GetPoint(typename Schema::IfcCartesianPoint* pPoint, Float64* pX, Float64* pY)
   {
      auto coordinates = pPoint->Coordinates();
      ATLASSERT(2 <= coordinates.size());
      *pX = WBFL::Units::ConvertToSysUnits(coordinates[0],*m_pLengthUnit);
      *pY = WBFL::Units::ConvertToSysUnits(coordinates[1],*m_pLengthUnit);
   }

   enum LastAlignmentType { Unknown, Line, Curve } m_LastAlignmentType;

   // adds a line to the alignment. returns the station at the end of the line
   template <typename Schema, typename Segment> Float64 OnLine(Float64 startStation, typename Segment* pLine);
   Float64 OnLine(Float64 sx, Float64 sy, Float64 startStation, Float64 startDirection, Float64 length);

   // adds a curve to the alignment. returns the station at the end of the curve
   template <typename Schema,typename SpiralType,typename CurveType>
   Float64 OnCurve(Float64 startStation, typename SpiralType* pEntrySpiral, typename CurveType* pCurve, typename SpiralType* pExitSpiral);

   template <typename Schema>
   Float64 OnCurve_4x3(Float64 startStation, typename Schema::IfcAlignmentHorizontalSegment* pEntrySpiral, typename Schema::IfcAlignmentHorizontalSegment* pCurve, typename Schema::IfcAlignmentHorizontalSegment* pExitSpiral);

   // adds linear segment to the profile
   template <typename Schema, typename LineSegmentType>
   void LinearSegment(Float64 startStation, typename LineSegmentType* pLinearSegment);

   template <typename Schema>
   void LinearSegment_4x3(Float64 startStation, typename Schema::IfcAlignmentVerticalSegment* pLinearSegment);

   // adds a parabolic curve to the proflie
   template <typename Schema, typename ParabolicSegmentType>
   void ParabolicSegment(Float64 startStation, typename ParabolicSegmentType* pParaCurve);

   template <typename Schema>
   void ParabolicSegment_4x3(Float64 startStation,typename Schema::IfcAlignmentVerticalSegment* pParaCurve);

   template <typename Schema,typename CurveType>
   void GetCurvePoints(typename CurveType* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter);

   template <typename Schema>
   void GetCurvePoints_4x3(typename Schema::IfcAlignmentHorizontalSegment* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter);

   template <typename Schema,typename SpiralType>
   void GetSpiralPoints(typename SpiralType* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd);

   template <typename Schema>
   void GetSpiralPoints_4x3(typename Schema::IfcAlignmentHorizontalSegment* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd);

   template <typename Schema, typename SpiralType>
   void CheckSpiralType(typename SpiralType* pSpiral);

   void CheckSpiralType_4x3(Ifc4x3_rc2::IfcAlignmentHorizontalSegment* pSpiral);
   void CheckSpiralType_4x3(Ifc4x3_rc3::IfcAlignmentHorizontalSegment* pSpiral);
   void CheckSpiralType_4x3(Ifc4x3_rc4::IfcAlignmentHorizontalSegment* pSpiral);

   //// NOTE: Not converting LandXML to cross slope data. This struct and the functions that follow are
   //// trial implementations, but they don't work.
   //struct DesignCrossSectData
   //{
   //   Float64 Station;
   //   LX::EnumSideofRoadType::Values Side;

   //   Float64 CrownPointOffset;
   //   Float64 LeftSlope;
   //   Float64 RightSlope;
   //};
   //void GetSlopes(LX::DesignCrossSectSurf* pDesignSurf, DesignCrossSectData* pSectionData);
   //static bool Compare(const DesignCrossSectData& a, const DesignCrossSectData& b);
   //void MergeSections(std::vector<DesignCrossSectData>& vSectionData);


   template <typename Schema>
   bool IsValidAlignment(typename Schema::IfcAlignment* pAlignment);

   template <typename Schema>
   bool IsValidAlignment_4x3(IfcParse::IfcFile& file, typename Schema::IfcAlignment* pAlignment);
};

