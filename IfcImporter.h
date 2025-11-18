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

#include <IFace/Tools.h>
#include <IFace\Project.h>

class CIfcImportOptions
{
public:
   enum class ModelElements
   {
      AlignmentOnly,
      AlignmentAndBridge
   };

   ModelElements model_elements = ModelElements::AlignmentAndBridge;
};


///////////////////////////////////////////////////////////////////////////
// CIfcImporter
//
// Converts data between IFC and PGSuper data structures
class CIfcImporter
{
public:
   CIfcImporter(void);
   ~CIfcImporter(void);

   // Converts Ifc data to PGSuper data
   HRESULT ImportFromIFC(std::shared_ptr<WBFL::EAF::Broker> pBroker, CString& strFilePath, CIfcImportOptions options);

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

   //LX::CrossSects* CreateCrossSections(std::shared_ptr<WBFL::EAF::Broker> pBroker, LX::IFactory* pFactory);
   //LX::Roadway*    CreateRoadway(std::shared_ptr<WBFL::EAF::Broker> pBroker, LX::IFactory* pFactory);

   bool ImportAlignment(std::shared_ptr<WBFL::EAF::Broker> pBroker, IfcParse::IfcFile& file);
   bool ImportBridge(std::shared_ptr<WBFL::EAF::Broker> pBroker, IfcParse::IfcFile& file);

   bool GetAlignmentParameters(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData);

   Ifc4x3_add2::IfcAlignment* GetAlignment(IfcParse::IfcFile& file);

   void InitUnits(IfcParse::IfcFile& file);

   void GetStations(Ifc4x3_add2::IfcAlignment* pAlignment, std::vector<std::pair<Float64, Float64>>& vStations, std::vector<std::tuple<Float64, Float64, Float64>>& vStationEquations);

   Float64 LoadAlignment(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment);

   void LoadProfile(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment,Float64 stationAdjustment);

   Float64 GetStartStation(Ifc4x3_add2::IfcAlignment* pAlignment);

   //void LoadCrossSections(LX::CrossSects* pCrossSects, LX::String& strSurfaceName);

   void GetPoint(Ifc4x3_add2::IfcCartesianPoint* pPoint, Float64* pX, Float64* pY)
   {
      auto coordinates = pPoint->Coordinates();
      ATLASSERT(2 <= coordinates.size());
      *pX = WBFL::Units::ConvertToSysUnits(coordinates[0],*m_pLengthUnit);
      *pY = WBFL::Units::ConvertToSysUnits(coordinates[1],*m_pLengthUnit);
   }

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

   // returns true if the alignment is a valid PGSuper alignment
   bool IsValidAlignment(IfcParse::IfcFile& file, Ifc4x3_add2::IfcAlignment* pAlignment);

   void SetGirderProperties(IfcParse::IfcFile& file, CBridgeDescription2& bridge_desc);
};

