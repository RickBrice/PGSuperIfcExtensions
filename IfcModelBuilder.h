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

class CIfcModelBuilderOptions
{
public:
   enum class Schema
   {
      //Schema_4x3_rc3,
      //Schema_4x3_rc4
      //Schema_4x3_tc1,
      //Schema_4x3_add1,
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

   enum class ObjectPlacement
   {
      Global, // IfcAxis2Placement3D
      Linear // IfcAxis2PlacementLinear
   };

   Schema schema = Schema::Schema_4x3_add2;
   ModelElements model_elements = ModelElements::AlignmentAndBridge;
   AlignmentModel alignment_model = AlignmentModel::GradientCurve;
   ObjectPlacement object_placement = ObjectPlacement::Linear;
   Tangents tangents = Tangents::Line;
};

///////////////////////////////////////////////////////////////////////////
// CIfcModelBuilder
class CIfcModelBuilder
{
public:
    CIfcModelBuilder(void);
    ~CIfcModelBuilder(void);

    bool BuildModel(IBroker* pBroker, const CString& strFilePath, const CIfcModelBuilderOptions& options);

private:
    template <typename Schema>
    bool BuildModel(IBroker* pBroker, const CString& strFilePath, const CIfcModelBuilderOptions& options);
};

