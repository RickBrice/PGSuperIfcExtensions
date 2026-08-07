#pragma once

#include <proj.h>
#include <vector>

// A trimmed-down, CString-based mirror of PROJ_CRS_INFO for use by the EPSG picker dialog.
struct EpsgCrsInfo
{
   CString Code;              // e.g. "2927"
   CString Name;              // e.g. "NAD83(HARN) / Washington South (ftUS)"
   CString AreaName;          // e.g. "United States (USA) - Washington - counties of ..."
   CString ProjectionMethod;  // e.g. "Lambert Conic Conformal (2SP)" - empty for vertical CRS
   bool Deprecated = false;

   // Area of use, in degrees. BboxValid is false if PROJ didn't have an extent on file for this CRS.
   bool BboxValid = false;
   double WestLon = 0.0;
   double SouthLat = 0.0;
   double EastLon = 0.0;
   double NorthLat = 0.0;
};

// Both lists are lazily populated on first call and cached for the lifetime of the process.
const std::vector<EpsgCrsInfo>& GetHorizontalCrsList(); // EPSG projected coordinate reference systems
const std::vector<EpsgCrsInfo>& GetVerticalCrsList();   // EPSG vertical coordinate reference systems

// Looks up the geodetic datum name for a horizontal (projected) EPSG code, e.g. "4326" -> "WGS 84".
// Returns an empty string if the code isn't found or has no associated datum/datum ensemble.
CString GetGeodeticDatumName(const CString& horizontalEpsgCode);

// True if the two CRS's areas of use overlap (accounting for antimeridian-crossing bboxes).
// Always false if either bbox is unknown - callers should treat "unknown" as "can't confirm
// relevance" rather than a positive or negative match.
bool EpsgAreaOfUseOverlaps(const EpsgCrsInfo& a, const EpsgCrsInfo& b);
