///////////////////////////////////////////////////////////////////////
// IFC Extension for PGSuper
// Copyright � 1999-2026  Washington State Department of Transportation
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

namespace WBFL
{
   namespace System
   {
      class IStructuredSave;
      class IStructuredLoad;
   };
};

struct GeoreferencingData
{
   // For IfcProjectedCRS
   CString Name = _T("2927"); // horizontal datum EPSG code, e.g. "EPSG:4326"
   CString Description= _T("Washington South (ftUS)");
   CString GeodeticDatum = _T("NAD83(HARN)");
   CString VerticalDatum = _T("5703");
   CString MapProjection = _T("Lambert Conformal Conic 2SP");

   // For IfcMapConversion
   Float64 Eastings = 1041929.;
   Float64 Northings = 630714.;
   Float64 OrthogonalHeight = 0.;
   Float64 XAxisAbscissa = 1.;
   Float64 XAxisOrdinate = 0.;
   Float64 Scale = 3937. / 1200.; // us survey foot to meter conversion factor

   void Save(WBFL::System::IStructuredSave* pSave) const;
   void Load(WBFL::System::IStructuredLoad* pLoad);
};

// {03917488-9929-41F1-AB85-8FED05E73009}
DEFINE_GUID(IID_IGeoreferencing,
   0x3917488, 0x9929, 0x41f1, 0xab, 0x85, 0x8f, 0xed, 0x5, 0xe7, 0x30, 0x9);
class IGeoreferencing
{
public:
   virtual void SetGeoreferencingData(const GeoreferencingData& data) = 0;
   virtual const GeoreferencingData& GetGeoreferencingData() const = 0;
};
