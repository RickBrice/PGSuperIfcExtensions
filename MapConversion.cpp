#include "stdafx.h"
#include "MapConversion.h"
#include <EAF\EAFApp.h>

double get_map_unit_to_meters(PJ_CONTEXT* C, const char* projected_epsg)
{
   USES_CONVERSION;
   CString strApp = EAFGetApp()->GetAppLocation();
   std::string db_path = T2A(strApp + _T("proj\\proj.db"));
   proj_context_set_database_path(C, db_path.c_str(), nullptr, nullptr);

   PJ* crs = proj_create(C, projected_epsg);
   if (crs == nullptr) {
      WATCH(_T("Failed to create CRS: ") << proj_context_errno_string(C,proj_context_errno(C)));
      return 1.0;
   }
   PJ* cs = proj_crs_get_coordinate_system(C, crs);
   if (cs == nullptr) {
      WATCH(_T("Failed to get coordinate system"));
      proj_destroy(crs);
      return 1.0;
   }
   double unit_conv_factor = 1.0;
   const char* unit_name = nullptr;
   proj_cs_get_axis_info(C, cs, 0, nullptr, nullptr, nullptr,
      &unit_conv_factor, &unit_name, nullptr, nullptr);
   WATCH(_T("Target CRS native unit: ") << (unit_name ? unit_name : "?") << _T(" (") << unit_conv_factor << _T(" m per unit)"));

   proj_destroy(cs);
   proj_destroy(crs);
   return unit_conv_factor;
}

std::string get_projected_epsg(std::shared_ptr<WBFL::EAF::Broker> pBroker)
{
   USES_CONVERSION;
   GET_IFACE2(pBroker, IGeoreferencing, pGeoRef);
   const auto& georef = pGeoRef->GetGeoreferencingData();

   CString projectedEpsg;
   projectedEpsg.Format(_T("EPSG:%s\n"), georef.Name);

   return T2A(projectedEpsg.LockBuffer());
}

std::pair<double, double> map_to_lonlat(std::shared_ptr<WBFL::EAF::Broker> pBroker,double easting, double norting)
{
   PJ_CONTEXT* C = proj_context_create();
   auto projected_epsg = get_projected_epsg(pBroker);
   auto lonlat = map_to_lonlat(C, projected_epsg.c_str(), easting, norting);
   proj_context_destroy(C);
   return lonlat;
}

std::pair<double, double> map_to_lonlat(PJ_CONTEXT* C, const char* projected_crs, double easting, double northing)
{
   PJ* P = proj_create_crs_to_crs(C, projected_crs, "EPSG:4326", nullptr);
   if (P == nullptr) {
      WATCH(_T("Failed to create transform to EPSG:4326"));
      return { 0, 0 };
   }
   PJ* norm = proj_normalize_for_visualization(C, P);
   proj_destroy(P);
   if (norm == nullptr) {
      WATCH(_T("Failed to normalize transform"));
      return { 0, 0 };
   }
   PJ_COORD in = proj_coord(easting, northing, 0, 0);
   PJ_COORD out = proj_trans(norm, PJ_FWD, in);
   proj_destroy(norm);
   return { out.xy.x, out.xy.y };
}

std::pair<double, double> lonlat_to_map(std::shared_ptr<WBFL::EAF::Broker> pBroker, double lon, double lat)
{
   PJ_CONTEXT* C = proj_context_create();
   auto projected_epsg = get_projected_epsg(pBroker);
   auto lonlat = lonlat_to_map(C, projected_epsg.c_str(), lon, lat);
   proj_context_destroy(C);
   return lonlat;
}

std::pair<double, double> lonlat_to_map(PJ_CONTEXT* C, const char* projected_crs, double lon, double lat)
{
   PJ* P = proj_create_crs_to_crs(C, "EPSG:4326", projected_crs, nullptr);
   if (P == nullptr) {
      WATCH(_T("Failed to create transform: ") << proj_context_errno_string(C, proj_context_errno(C)));
      return { 0, 0 };
   }

   // Forces conventional (lon, lat) input / (easting, northing) output
   // order, rather than each CRS's authority-defined order - same
   // reasoning as in the earlier EPSG:2927 -> WGS84 example.
   PJ* norm = proj_normalize_for_visualization(C, P);
   proj_destroy(P);
   if (norm == nullptr) {
      WATCH(_T("Failed to normalize transform"));
      return { 0, 0 };
   }

   PJ_COORD in = proj_coord(lon, lat, 0, 0);
   PJ_COORD out = proj_trans(norm, PJ_FWD, in);
   proj_destroy(norm);

   return { out.xy.x, out.xy.y };
}

double get_grid_scale_factor(PJ_CONTEXT* C, const char* projected_epsg,  double lon_deg, double lat_deg) {
   PJ* P = proj_create(C, projected_epsg);
   if (P == nullptr) {
      WATCH(_T("Failed to create CRS object for factors"));
      return 1.0;
   }
   PJ_COORD lp = proj_coord(proj_torad(lon_deg), proj_torad(lat_deg), 0, 0);
   PJ_FACTORS factors = proj_factors(P, lp);
   proj_destroy(P);

   WATCH(_T("Grid scale factor: parallel=") << factors.parallel_scale << _T(" meridional=") << factors.meridional_scale 
      << _T(" (should closely agree - conformal projection)"));
   return factors.parallel_scale;
}