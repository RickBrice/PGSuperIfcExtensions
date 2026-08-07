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
#include "IfcGeoreferencingImporter.h"
#include "IfcImporter.h"
#include <IFace/Tools.h>
#include "Georeferencing.h"

namespace
{
   // Strips a leading "EPSG:" (case-insensitive), if present, so GeoreferencingData::Name and
   // ::VerticalDatum consistently never carry the prefix internally - the exporter always adds
   // it back unconditionally, so leaving it in here would double it up on a re-export.
   CString StripEpsgPrefix(const CString& value)
   {
      CString v(value);
      v.TrimLeft();
      if (v.Left(5).CompareNoCase(_T("EPSG:")) == 0)
         v = v.Mid(5);
      return v;
   }
}

CIfcGeoreferencingImporter::CIfcGeoreferencingImporter(CIfcImporter& importer) :
   m_Importer(importer)
{
}

IfcSchema::IfcProjectedCRS* CIfcGeoreferencingImporter::GetProjectedCRS(IfcParse::IfcFile& file)
{
   auto projected_crs_list = file.instances_by_type<IfcSchema::IfcProjectedCRS>();
   if (projected_crs_list->size() < 1)
   {
      WBFL::System::Logger::Info(_T("IFC model does not contain an IfcProjectedCRS. Georeferencing will not be imported."));
      return nullptr;
   }

   return *(projected_crs_list->begin());
}

IfcSchema::IfcMapConversion* CIfcGeoreferencingImporter::GetMapConversion(IfcParse::IfcFile& file)
{
   auto map_conversion_list = file.instances_by_type<IfcSchema::IfcMapConversion>();
   if (map_conversion_list->size() < 1)
   {
      WBFL::System::Logger::Info(_T("IFC model does not contain an IfcMapConversion. Map conversion will be computed for this model."));
      return nullptr;
   }

   return *(map_conversion_list->begin());
}

CIfcImporter::ImportResult CIfcGeoreferencingImporter::Import(IfcParse::IfcFile& file)
{
   GET_IFACE2(m_Importer.GetBroker(), IGeoreferencing, pGeoRef);
   auto georefdata = pGeoRef->GetGeoreferencingData();

   auto projected_crs = GetProjectedCRS(file);
   if (!projected_crs)
   {
      return CIfcImporter::ImportResult::NotFound;
   }

   georefdata.IsCRSValid = true;
   georefdata.Name = StripEpsgPrefix(CString(projected_crs->Name().value_or("").c_str()));
   georefdata.Description = CString(projected_crs->Description().value_or("").c_str());
   georefdata.GeodeticDatum = CString(projected_crs->GeodeticDatum().value_or("").c_str());
   georefdata.VerticalDatum = StripEpsgPrefix(CString(projected_crs->VerticalDatum().value_or("").c_str()));
   georefdata.MapProjection = CString(projected_crs->MapProjection().value_or("").c_str());
   georefdata.MapZone = CString(projected_crs->MapZone().value_or("").c_str());
   ImportMapUnit<IfcSchema>(projected_crs->MapUnit(), georefdata);

   auto map_conversion = GetMapConversion(file);
   if (map_conversion)
   {
      georefdata.IsMapConversionValid = true;
      georefdata.Eastings = map_conversion->Eastings();
      georefdata.Northings = map_conversion->Northings();
      georefdata.OrthogonalHeight = map_conversion->OrthogonalHeight();
      georefdata.XAxisAbscissa = map_conversion->XAxisAbscissa().value_or(1.0); // absent means "no rotation"
      georefdata.XAxisOrdinate = map_conversion->XAxisOrdinate().value_or(0.0);
      georefdata.Scale = map_conversion->Scale().value_or(1.0);
   }
   else
   {
      georefdata.IsMapConversionValid = false;
   }

   pGeoRef->SetGeoreferencingData(georefdata);

   return CIfcImporter::ImportResult::Success;
}
