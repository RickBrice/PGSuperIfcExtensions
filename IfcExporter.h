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

class CIfcExportOptions
{
public:
   enum class Schema
   {
      Schema_4x3_add2
   };

   enum class ModelElements
   {
      AlignmentOnly,
      AlignmentAndBridge
   };

   enum class AlignmentModel
   {
      Polyline, // IfcPolyline (3D wire)
      GradientCurve // IfcGradientCurve
   };

   enum class Tangents
   {
      Polyline, // use IfcPolyline to model tangents
      Line // use IfcLine to model tangents
   };

   enum class Representations
   {
      Curve3dOnly,
      Curve3dAndFootPrint
   };

   enum class SweepProfile
   {
      Polyline, // use IfcPolyline for sweep profile
      IndexedPolyCurve // use IfcIndexedPolyCurve for sweep profile
   };

   enum class Railings
   {
      Parapet, // use IfcWall.PARAPET
      Balustrade // Use IfcRailing.BALUSTRADE
   };

   enum class BeamPlacement
   {
      Linear,
      Local
   };

   enum class BeamModel
   {
      SectionedSolidHorizontal,
      PolygonalFaceSet,
      FacetedBrep
   };

   Schema schema = Schema::Schema_4x3_add2;
   ModelElements model_elements = ModelElements::AlignmentAndBridge;
   bool classify = true;
   AlignmentModel alignment_model = AlignmentModel::GradientCurve;
   Tangents tangents = Tangents::Line;
   Representations representations = Representations::Curve3dOnly;
   SweepProfile sweep_profile = SweepProfile::Polyline;
   Railings railings = Railings::Balustrade;
   bool include_work_plan = true;
   bool include_rebar = true;
   bool rebar_per_aci131 = true;
   bool include_camber = true;
   bool include_quantities = true;
   BeamPlacement beam_placement = BeamPlacement::Linear;
   BeamModel beam_model = BeamModel::SectionedSolidHorizontal;

   // Sometimes when girders are installed at a very steep angle, the ends of the girders
   // are battered so that the end faces are vertical when the beam is erected.
   // This parameter controls batter. It does two things:
   // 1) It batters the end faces of the girder, BUT DOES NOT ADJUST THE REINFORCEMENT (@todo - update reinforcement for batter)
   // 2) Puts the batter angle in the Pset_PrecastConcreteElementGeneral
   bool batter_ends = false; // set to false because PGSuper doesn't support batter
};

///////////////////////////////////////////////////////////////////////////
// CIfcExporter
class CIfcExporter
{
public:
    CIfcExporter(void);
    ~CIfcExporter(void);

    bool BuildModel(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CString& strFilePath, const CIfcExportOptions& options);

private:
    template <typename Schema>
    bool BuildModel(std::shared_ptr<WBFL::EAF::Broker> pBroker, const CString& strFilePath, const CIfcExportOptions& options);
};


template <typename Schema>
typename Schema::IfcCartesianPoint* ConvertPoint(IPoint2d* pPoint, bool bMirror = false)
{
   Float64 x, y;
   pPoint->Location(&x, &y);
   x = IsZero(x) ? 0.0 : x;
   y = IsZero(y) ? 0.0 : y;
   return new typename Schema::IfcCartesianPoint(std::vector<double>{bMirror ? -x : x, y});
}
