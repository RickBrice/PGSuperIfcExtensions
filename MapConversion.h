#pragma once

#include <proj.h>
#include <utility>

#include "Georeferencing.h"
#include <IFace/Project.h>

CString GetProjDbPath(); // path to proj.db, deployed alongside this plugin's exe
PJ_CONTEXT* CreatePjContext(); // proj_context_create() with the proj.db path already configured

double get_map_unit_to_meters(PJ_CONTEXT* C, const char* projected_epsg);
std::string get_projected_epsg(std::shared_ptr<WBFL::EAF::Broker> pBroker);
std::pair<double, double> map_to_lonlat(std::shared_ptr<WBFL::EAF::Broker> pBroker,double easting, double norting);
std::pair<double, double> map_to_lonlat(PJ_CONTEXT* C, const char* projected_crs, double easting, double norting);

std::pair<double, double> lonlat_to_map(std::shared_ptr<WBFL::EAF::Broker> pBroker, double lon, double lat);
std::pair<double, double> lonlat_to_map(PJ_CONTEXT* C, const char* projected_crs, double lon, double lat);

double get_grid_scale_factor(PJ_CONTEXT* C, const char* projected_epsg, double lon_deg, double lat_deg);

template<typename Schema>
typename Schema::IfcMapConversion* create_map_conversion(IfcHierarchyHelper<Schema>& file, std::shared_ptr<WBFL::EAF::Broker> pBroker, typename Schema::IfcProjectedCRS* projected_crs)
{
   auto projected_epsg = get_projected_epsg(pBroker);

   GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
   auto alignment_data = pRoadwayData->GetAlignmentData2();

   PJ_CONTEXT* C = proj_context_create();

   auto map_unit_to_meters = get_map_unit_to_meters(C, projected_epsg.c_str());

   auto anchor = map_to_lonlat(C, projected_epsg.c_str(), alignment_data.xRefPoint, alignment_data.yRefPoint);

   auto grid_scale_factor = get_grid_scale_factor(C, projected_epsg.c_str(), anchor.first, anchor.second);

   double project_unit_to_meters = 1.0; // this project's unit is meter
   auto scale = (project_unit_to_meters / map_unit_to_meters) * grid_scale_factor;

   proj_context_destroy(C);

   auto geometric_representation_context = file.getRepresentationContext(std::string("Model")); // creates the representation context if it doesn't already exist
   auto map_conversion = new typename Schema::IfcMapConversion(
      geometric_representation_context,
      projected_crs,
      0.0, 0.0, 0.0, // Eastings, Northings, OrthogonalHeight,
      1.0, // XAxisAbscissa,
      0.0, // XAxisOrdinate,
      scale);

   return map_conversion;
};