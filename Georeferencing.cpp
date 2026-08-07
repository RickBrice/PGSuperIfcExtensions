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

#include "stdafx.h"
#include "Georeferencing.h"

static LPCTSTR UnitName = _T("Georeferencing");
static LPCTSTR ProjectedCRSUnitName = _T("ProjectedCRS");
static LPCTSTR MapConversionUnitName = _T("MapConversion");

void GeoreferencingData::Save(WBFL::System::IStructuredSave* pSave) const
{
   pSave->BeginUnit(UnitName, 2.0);

   pSave->BeginUnit(ProjectedCRSUnitName, 2.0);
   pSave->Property(_T("Name"), Name);
   pSave->Property(_T("Description"), Description);
   pSave->Property(_T("GeodeticDatum"), GeodeticDatum);
   pSave->Property(_T("VerticalDatum"), VerticalDatum);
   pSave->Property(_T("MapProjection"), MapProjection);
   pSave->Property(_T("IsCRSValid"), IsCRSValid);
   pSave->Property(_T("MapZone"), MapZone);
   pSave->Property(_T("IsMapUnitSI"), IsMapUnitSI);
   pSave->Property(_T("MapUnitName"), MapUnitName);
   pSave->Property(_T("MapUnitToMeters"), MapUnitToMeters);
   pSave->EndUnit(); // ProjectedCRS

   pSave->BeginUnit(MapConversionUnitName, 1.0);
   pSave->Property(_T("IsMapConversionValid"), IsMapConversionValid);
   pSave->Property(_T("Eastings"), Eastings);
   pSave->Property(_T("Northings"), Northings);
   pSave->Property(_T("OrthogonalHeight"), OrthogonalHeight);
   pSave->Property(_T("XAxisAbscissa"), XAxisAbscissa);
   pSave->Property(_T("XAxisOrdinate"), XAxisOrdinate);
   pSave->Property(_T("Scale"), Scale);
   pSave->EndUnit(); // MapConversion

   pSave->EndUnit(); // Georeferencing
}

void GeoreferencingData::Load(WBFL::System::IStructuredLoad* pLoad)
{
   if (!pLoad->BeginUnit(UnitName)) THROW_LOAD(InvalidFileFormat, pLoad);
   Float64 georef_version = pLoad->GetVersion();

   if (!pLoad->BeginUnit(ProjectedCRSUnitName)) THROW_LOAD(InvalidFileFormat, pLoad);
   Float64 crs_version = pLoad->GetVersion();
   std::_tstring value;
   if (!pLoad->Property(_T("Name"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
   Name = value.c_str();
   if (!pLoad->Property(_T("Description"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
   Description = value.c_str();
   if (!pLoad->Property(_T("GeodeticDatum"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
   GeodeticDatum = value.c_str();
   if (!pLoad->Property(_T("VerticalDatum"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
   VerticalDatum = value.c_str();
   if (!pLoad->Property(_T("MapProjection"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
   MapProjection = value.c_str();

   if (2.0 <= crs_version)
   {
      bool bValue;
      if (!pLoad->Property(_T("IsCRSValid"), &bValue)) THROW_LOAD(InvalidFileFormat, pLoad);
      IsCRSValid = bValue;
      if (!pLoad->Property(_T("MapZone"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
      MapZone = value.c_str();
      if (!pLoad->Property(_T("IsMapUnitSI"), &bValue)) THROW_LOAD(InvalidFileFormat, pLoad);
      IsMapUnitSI = bValue;
      if (!pLoad->Property(_T("MapUnitName"), &value)) THROW_LOAD(InvalidFileFormat, pLoad);
      MapUnitName = value.c_str();
      if (!pLoad->Property(_T("MapUnitToMeters"), &MapUnitToMeters)) THROW_LOAD(InvalidFileFormat, pLoad);
   }
   // else: older file format predates these properties - leave the in-class default values in place
   if (!pLoad->EndUnit()) THROW_LOAD(InvalidFileFormat, pLoad); // ProjectedCRS

   if (2.0 <= georef_version)
   {
      if (!pLoad->BeginUnit(MapConversionUnitName)) THROW_LOAD(InvalidFileFormat, pLoad);
      bool bValue;
      if (!pLoad->Property(_T("IsMapConversionValid"), &bValue)) THROW_LOAD(InvalidFileFormat, pLoad);
      IsMapConversionValid = bValue;
      if (!pLoad->Property(_T("Eastings"), &Eastings)) THROW_LOAD(InvalidFileFormat, pLoad);
      if (!pLoad->Property(_T("Northings"), &Northings)) THROW_LOAD(InvalidFileFormat, pLoad);
      if (!pLoad->Property(_T("OrthogonalHeight"), &OrthogonalHeight)) THROW_LOAD(InvalidFileFormat, pLoad);
      if (!pLoad->Property(_T("XAxisAbscissa"), &XAxisAbscissa)) THROW_LOAD(InvalidFileFormat, pLoad);
      if (!pLoad->Property(_T("XAxisOrdinate"), &XAxisOrdinate)) THROW_LOAD(InvalidFileFormat, pLoad);
      if (!pLoad->Property(_T("Scale"), &Scale)) THROW_LOAD(InvalidFileFormat, pLoad);
      if (!pLoad->EndUnit()) THROW_LOAD(InvalidFileFormat, pLoad); // MapConversion
   }
   // else: older file format never wrote the MapConversion unit - leave the in-class default values in place

   if (!pLoad->EndUnit()) THROW_LOAD(InvalidFileFormat, pLoad); // Georeferencing
}
