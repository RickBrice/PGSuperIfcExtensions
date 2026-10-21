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

#include "stdafx.h"
#include "Stirrups.h"

#include <psgLib/ShearData.h>
#include <LRFD/RebarPool.h>

#include <iomanip>

#undef min // undef our version of min in MathEx.h so std::min is used correctly in this file
#undef max

namespace
{
   using WBFL::Materials::Rebar;

   const Float64 max_leg_slope = tan(WBFL::Units::ConvertToSysUnits(15.0, WBFL::Units::Measure::Degree)); // horizontal per vertical
   const Float64 station_tolerance = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch); // bars this close along the girder are at one station
   const Float64 spacing_tolerance = WBFL::Units::ConvertToSysUnits(0.25, WBFL::Units::Measure::Inch);
   const Float64 symmetry_tolerance = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
   const Float64 projection_tolerance = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch); // a leg's top this far above the girder engages the deck
   const Float64 default_min_leg = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Feet); // when the girder's depth isn't known

   // A stirrup: its position along the girder (the middle of its legs) and its legs
   struct Stirrup
   {
      Float64 x;
      Rebar::Size size;
      IndexType nLegs = 0, nInterfaceLegs = 0;
   };

   // The bar as a stirrup, or nullopt if it has no vertical legs
   std::optional<Stirrup> stirrup_of(const GirderBar& bar, const GirderReinforcementReader& geometry)
   {
      const auto& path = bar.centerline.path;
      Stirrup stirrup{ 0.0, bar.size };
      Float64 sum_x = 0;
      for (size_t i = 1; i < path.size(); i++)
      {
         const auto& a = path[i - 1];
         const auto& b = path[i];
         Float64 dy = fabs(b.y - a.y);
         if (max_leg_slope * dy < sqrt(pow(b.x - a.x, 2) + pow(b.z - a.z, 2)))
            continue;
         Float64 x = (a.x + b.x) / 2;
         auto depth = geometry.GetDepth(bar.centerline.beam_id, x);
         if (dy < (depth ? *depth / 2 : default_min_leg))
            continue;
         stirrup.nLegs++;
         sum_x += x;
         if (std::min(a.y, b.y) < -projection_tolerance) // y is the depth below the girder's top
            stirrup.nInterfaceLegs++;
      }
      if (stirrup.nLegs == 0)
         return std::nullopt;
      stirrup.x = sum_x / stirrup.nLegs;
      return stirrup;
   }

   // The stirrups at one place along the girder
   struct Station
   {
      Float64 x;
      Rebar::Size size;
      IndexType nLegs = 0, nInterfaceLegs = 0;
   };

   bool same_bars(const Station& a, const Station& b)
   {
      return a.size == b.size && a.nLegs == b.nLegs && a.nInterfaceLegs == b.nInterfaceLegs;
   }

   // The stations of the stirrups (sorted by x): stirrups within station_tolerance of the first of a station are at it. Its
   // size is the one with the most legs; nMixed counts the stations with more than one size
   std::vector<Station> get_stations(std::vector<Stirrup> stirrups, IndexType& nMixed)
   {
      std::sort(stirrups.begin(), stirrups.end(), [](const auto& a, const auto& b) {return a.x < b.x; });
      std::vector<Station> stations;
      for (size_t i = 0; i < stirrups.size();)
      {
         size_t j = i;
         Station station{ 0.0, Rebar::Size::bsNone };
         std::map<Rebar::Size, IndexType> legs;
         while (j < stirrups.size() && stirrups[j].x - stirrups[i].x <= station_tolerance)
         {
            station.x += stirrups[j].x;
            station.nLegs += stirrups[j].nLegs;
            station.nInterfaceLegs += stirrups[j].nInterfaceLegs;
            legs[stirrups[j].size] += stirrups[j].nLegs;
            j++;
         }
         station.x /= (j - i);
         station.size = std::max_element(legs.begin(), legs.end(), [](const auto& a, const auto& b) {return a.second < b.second; })->first;
         if (1 < legs.size())
            nMixed++;
         stations.push_back(station);
         i = j;
      }
      return stations;
   }

   // the stations have mirror stations, at the same distance from the other end with the same bars
   bool is_symmetric(const std::vector<Station>& stations, Float64 Lg)
   {
      return std::all_of(stations.begin(), stations.end(), [&](const auto& s)
         {
            return std::any_of(stations.begin(), stations.end(), [&](const auto& m) {return same_bars(s, m) && IsEqual(m.x, Lg - s.x, symmetry_tolerance); });
         });
   }

   // The zones of the stations: the space before each station (from the girder's start) is a spacing, and consecutive equal
   // spaces with the same bars make a zone. The last zone is infinitely long, as PGSuper's last zone is
   std::vector<CShearZoneData2> make_zones(const std::vector<Station>& stations)
   {
      std::vector<CShearZoneData2> zones;
      std::vector<IndexType> nSpaces;
      Float64 previous = 0;
      const Station* pLast = nullptr;
      for (const auto& station : stations)
      {
         Float64 space = station.x - previous;
         previous = station.x;
         if (pLast && same_bars(*pLast, station) && IsEqual(space, zones.back().BarSpacing, spacing_tolerance))
         {
            auto& zone = zones.back();
            zone.ZoneLength += space;
            nSpaces.back()++;
            zone.BarSpacing = zone.ZoneLength / nSpaces.back();
         }
         else
         {
            CShearZoneData2 zone;
            zone.ZoneNum = zones.size();
            zone.ZoneLength = space;
            zone.BarSpacing = space;
            zone.VertBarSize = station.size;
            zone.nVertBars = (Float64)station.nLegs;
            zone.nHorzInterfaceBars = (Float64)station.nInterfaceLegs;
            zone.ConfinementBarSize = Rebar::Size::bsNone;
            zones.push_back(zone);
            nSpaces.push_back(1);
         }
         pLast = &station;
      }
      if (!zones.empty())
         zones.back().ZoneLength = Float64_Max;
      return zones;
   }
}

void import_stirrups(const GirderBars& girder_bars, CBridgeDescription2& bridge_desc)
{
   USES_CONVERSION;

   // The girders' stirrups come from the template, which has the library entries' transverse reinforcement. That's a seed
   // value, not the girder's reinforcement, so a girder gets only the stirrups in the model.
   for_each_segment(bridge_desc, [](CPrecastSegmentData* pSegment)
      {
         CShearData2 shear_data;
         shear_data.ShearBarType = pSegment->ShearData.ShearBarType; // the material if the model doesn't say
         shear_data.ShearBarGrade = pSegment->ShearData.ShearBarGrade;
         pSegment->ShearData = shear_data;
      });
   if (!girder_bars.HasBars())
   {
      WBFL::System::Logger::Info("Stirrups aren't imported: " + girder_bars.GetReason() + ". The girders have no transverse reinforcement.");
      return;
   }

   const auto& geometry = girder_bars.GetGeometry();
   IndexType nWithoutStirrups = 0, nWithoutSize = 0;
   for (const auto& [beam_id, girder] : geometry.GetGirders())
   {
      std::vector<Stirrup> stirrups;
      std::map<BarMaterial, IndexType> materials;
      for (const auto& bar : girder_bars.GetBars(beam_id))
      {
         auto stirrup = stirrup_of(bar, geometry);
         if (!stirrup)
            continue;
         if (bar.size == Rebar::Size::bsNone)
         {
            nWithoutSize++;
            continue;
         }
         if (bar.material)
            materials[*bar.material]++;
         stirrups.push_back(*stirrup);
      }
      if (stirrups.empty())
      {
         nWithoutStirrups++;
         continue;
      }

      std::string girder_name = girder_bars.GetGirderName(beam_id);
      auto* pSegment = GirderBars::GetSegment(bridge_desc, *girder);
      auto& shear_data = pSegment->ShearData;

      bool bMixed;
      if (auto material = most_common_material(materials, bMixed))
      {
         shear_data.ShearBarType = material->first;
         shear_data.ShearBarGrade = material->second;
      }
      if (bMixed)
         WBFL::System::Logger::Warning("Stirrups of " + girder_name + ": the bars have different materials. PGSuper has one for the girder's transverse reinforcement, so the most common is used.");

      IndexType nMixedSizes = 0;
      auto stations = get_stations(stirrups, nMixedSizes);
      if (0 < nMixedSizes)
         WBFL::System::Logger::Warning("Stirrups of " + girder_name + ": " + std::to_string(nMixedSizes) + " places along the girder have stirrups of different sizes. PGSuper has one size for a zone, so the size with the most legs is used.");
      IndexType nStations = stations.size();

      Float64 Lg = geometry.GetLength(beam_id);
      shear_data.bAreZonesSymmetrical = is_symmetric(stations, Lg);
      if (shear_data.bAreZonesSymmetrical)
      {
         // the first half (the stirrups at the center are in the last zone, which goes to midspan)
         stations.erase(std::remove_if(stations.begin(), stations.end(), [Lg](const auto& s) {return Lg / 2 - symmetry_tolerance < s.x; }), stations.end());
      }
      shear_data.ShearZones = make_zones(stations);
      if (shear_data.ShearZones.empty())
         shear_data.ShearZones = CShearData2().ShearZones; // only a stirrup at the center: PGSuper needs a zone

      std::ostringstream os;
      os << "Stirrups of " << girder_name << ": " << stirrups.size() << " bars at " << nStations << " places in " << shear_data.ShearZones.size()
         << (shear_data.bAreZonesSymmetrical ? " symmetric" : "") << " zones (";
      for (const auto& zone : shear_data.ShearZones)
      {
         if (zone.ZoneNum != 0)
            os << "; ";
         os << T2A(WBFL::LRFD::RebarPool::GetBarSize(zone.VertBarSize).c_str()) << " " << zone.nVertBars << " legs";
         if (0 < zone.nHorzInterfaceBars)
            os << " (" << zone.nHorzInterfaceBars << " into the deck)";
         os << " @ " << std::fixed << std::setprecision(2) << WBFL::Units::ConvertFromSysUnits(zone.BarSpacing, WBFL::Units::Measure::Inch) << " in";
         if (zone.ZoneLength < Float64_Max)
            os << " for " << WBFL::Units::ConvertFromSysUnits(zone.ZoneLength, WBFL::Units::Measure::Inch) << " in";
         os << std::defaultfloat;
      }
      os << "), " << T2A(WBFL::LRFD::RebarPool::GetMaterialName(shear_data.ShearBarType, shear_data.ShearBarGrade).c_str()) << ".";
      WBFL::System::Logger::Info(os.str());
   }

   if (0 < nWithoutSize)
      WBFL::System::Logger::Warning(std::to_string(nWithoutSize) + " stirrups were skipped: their size isn't known (rebar.bar_size, the nominal diameter, or a drawn diameter of a bar size).");
   if (0 < nWithoutStirrups)
      WBFL::System::Logger::Info(std::to_string(nWithoutStirrups) + " girders have no stirrups in the model, so they have no transverse reinforcement.");
   WBFL::System::Logger::Info("Confinement and splitting reinforcement aren't imported yet: the girders have none.");
}
