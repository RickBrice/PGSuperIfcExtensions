#include "stdafx.h"
#include "EpsgCatalog.h"
#include "MapConversion.h"

#include <algorithm>

namespace
{
   // Splits a west/east longitude range into one or two [west,east] pieces with west <= east,
   // to handle ranges that cross the antimeridian (where PROJ reports west > east).
   std::vector<std::pair<double, double>> SplitLonRange(double west, double east)
   {
      if (west <= east)
         return { { west, east } };
      return { { west, 180.0 }, { -180.0, east } };
   }

   bool LonRangesOverlap(double aWest, double aEast, double bWest, double bEast)
   {
      for (const auto& a : SplitLonRange(aWest, aEast))
      {
         for (const auto& b : SplitLonRange(bWest, bEast))
         {
            if (a.first <= b.second && b.first <= a.second)
               return true;
         }
      }
      return false;
   }

   std::vector<EpsgCrsInfo> FetchCrsList(PJ_TYPE type)
   {
      std::vector<EpsgCrsInfo> result;

      PJ_CONTEXT* C = CreatePjContext();

      PROJ_CRS_LIST_PARAMETERS* params = proj_get_crs_list_parameters_create();
      params->types = &type;
      params->typesCount = 1;
      params->allow_deprecated = FALSE;

      int count = 0;
      PROJ_CRS_INFO** list = proj_get_crs_info_list_from_database(C, "EPSG", params, &count);
      if (list)
      {
         result.reserve(count);
         for (int i = 0; i < count; i++)
         {
            EpsgCrsInfo info;
            info.Code = CString(list[i]->code);
            info.Name = CString(list[i]->name);
            info.AreaName = list[i]->area_name ? CString(list[i]->area_name) : CString();
            info.ProjectionMethod = list[i]->projection_method_name ? CString(list[i]->projection_method_name) : CString();
            info.Deprecated = (list[i]->deprecated != 0);
            info.BboxValid = (list[i]->bbox_valid != 0);
            info.WestLon = list[i]->west_lon_degree;
            info.SouthLat = list[i]->south_lat_degree;
            info.EastLon = list[i]->east_lon_degree;
            info.NorthLat = list[i]->north_lat_degree;
            result.push_back(info);
         }
         proj_crs_info_list_destroy(list);
      }
      else
      {
         WBFL::System::Logger::Info(_T("Failed to retrieve EPSG CRS list from PROJ database."));
      }

      proj_get_crs_list_parameters_destroy(params);
      proj_context_destroy(C);

      std::sort(result.begin(), result.end(), [](const EpsgCrsInfo& a, const EpsgCrsInfo& b)
         {
            return a.Code.CompareNoCase(b.Code) < 0;
         });

      return result;
   }
}

const std::vector<EpsgCrsInfo>& GetHorizontalCrsList()
{
   static std::vector<EpsgCrsInfo> list = FetchCrsList(PJ_TYPE_PROJECTED_CRS);
   return list;
}

const std::vector<EpsgCrsInfo>& GetVerticalCrsList()
{
   static std::vector<EpsgCrsInfo> list = FetchCrsList(PJ_TYPE_VERTICAL_CRS);
   return list;
}

CString GetGeodeticDatumName(const CString& horizontalEpsgCode)
{
   USES_CONVERSION;
   CString result;

   PJ_CONTEXT* C = CreatePjContext();

   std::string code = T2A(horizontalEpsgCode);
   PJ* crs = proj_create_from_database(C, "EPSG", code.c_str(), PJ_CATEGORY_CRS, false, nullptr);
   if (crs)
   {
      PJ* datum = proj_crs_get_datum(C, crs);
      if (!datum)
         datum = proj_crs_get_datum_ensemble(C, crs);

      if (datum)
      {
         const char* name = proj_get_name(datum);
         if (name)
            result = CString(name);
         proj_destroy(datum);
      }
      else
      {
         WBFL::System::Logger::Info(_T("EPSG code has no associated geodetic datum or datum ensemble."));
      }

      proj_destroy(crs);
   }
   else
   {
      WBFL::System::Logger::Info(_T("Failed to look up EPSG code in PROJ database."));
   }

   proj_context_destroy(C);

   return result;
}

bool EpsgAreaOfUseOverlaps(const EpsgCrsInfo& a, const EpsgCrsInfo& b)
{
   if (!a.BboxValid || !b.BboxValid)
      return false;

   if (a.NorthLat < b.SouthLat || a.SouthLat > b.NorthLat)
      return false;

   return LonRangesOverlap(a.WestLon, a.EastLon, b.WestLon, b.EastLon);
}
