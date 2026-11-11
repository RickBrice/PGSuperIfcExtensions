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
#include "SplicedGirders.h"
#include "GirderProfile.h"
#include "Piers.h"
#include "IfcImporter.h"
#include "IfcTargetReader.h"

#include <PsgLib/GirderLabel.h>
#include <PsgLib/TemporarySupportData.h>
#include <PsgLib/ClosureJointData.h>
#include <PsgLib/TimelineManager.h>

#include <algorithm>
#include <numeric>

namespace
{
   const Float64 inch = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Inch);
   const Float64 feet = WBFL::Units::ConvertToSysUnits(1.0, WBFL::Units::Measure::Feet);

   // stations within this of a pier are at the pier (a segment end, the edge of a gap)
   const Float64 station_tolerance = inch;

   std::_tstring ft(Float64 value)
   {
      std::_tostringstream os;
      os << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Feet) << _T(" ft");
      return os.str();
   }

   std::_tstring in(Float64 value)
   {
      std::_tostringstream os;
      os << std::fixed << std::setprecision(3) << WBFL::Units::ConvertFromSysUnits(value, WBFL::Units::Measure::Inch) << _T(" in");
      return os.str();
   }

   Float64 mid_station(const GirderCenterline& g) { return (g.start_station + g.end_station) / 2; }

   // the point, top elevation, and offset on a segment's centerline at a station (linear between its ends)
   void at_station(const GirderCenterline& segment, Float64 station, Eigen::Vector2d& point, Float64& top, Float64& offset)
   {
      Float64 length = segment.end_station - segment.start_station;
      Float64 t = IsZero(length) ? 0.0 : (station - segment.start_station) / length;
      point = segment.start + t * (segment.end - segment.start);
      top = std::lerp(segment.start_top, segment.end_top, t);
      offset = std::lerp(segment.start_offset, segment.end_offset, t);
   }

   // What a girder line does at a pier
   enum class AtPier { NotThere, Ends, Continuous, ClosureJoint };

   AtPier at_pier(const SplicedGirder& girder, Float64 station)
   {
      const auto& segments = girder.segments;
      if (station < segments.front().start_station - station_tolerance || segments.back().end_station + station_tolerance < station)
         return AtPier::NotThere;

      for (const auto& segment : segments)
      {
         if (segment.start_station + station_tolerance < station && station < segment.end_station - station_tolerance)
            return AtPier::Continuous;
      }

      for (size_t i = 0; i + 1 < segments.size(); i++)
      {
         if (segments[i].end_station - station_tolerance <= station && station <= segments[i + 1].start_station + station_tolerance)
            return girder.closure_joints[i] ? AtPier::ClosureJoint : AtPier::Ends;
      }

      return AtPier::Ends; // at the start or end of the girder line
   }

   LPCTSTR describe(AtPier at)
   {
      switch (at)
      {
      case AtPier::Ends: return _T("ends");
      case AtPier::Continuous: return _T("is continuous");
      case AtPier::ClosureJoint: return _T("has a closure joint");
      default: return _T("isn't there");
      }
   }

   // the stations of the erection towers (the table's erection_tower role), from the referents that position them
   std::vector<Float64> get_erection_tower_stations(ifcopenshell::file& file)
   {
      std::vector<Float64> stations;
      for (auto object : CIfcImporter::GetTargetReader().Select(ElementKind::ErectionTower, file))
      {
         auto part = object.as<IfcSchema::IfcBridgePart>();
         if (!part)
            continue;
         if (auto station = get_referent_station(part))
         {
            stations.push_back(*station);
         }
         else
         {
            std::ostringstream os;
            os << "The erection tower " << part.Name().value_or("unnamed") << " (#" << part.id() << ") isn't positioned by a referent with a station (Pset_Stationing), so it isn't imported.";
            WBFL::System::Logger::Warning(os.str());
         }
      }
      std::sort(stations.begin(), stations.end());
      return stations;
   }
}

IndexType count_spliced_girders(ifcopenshell::file& file)
{
   const auto& reader = CIfcImporter::GetTargetReader();
   IndexType nSpliced = 0;
   std::set<int> assembled;
   for (auto assembly : reader.Select(ElementKind::GirderAssembly, file))
   {
      IndexType nBeams = 0;
      bool bClosureJoint = false;
      for (auto rel : assembly.IsDecomposedBy())
      {
         for (auto part : rel.RelatedObjects())
         {
            if (auto beam = part.as<IfcSchema::IfcBeam>())
            {
               nBeams++;
               assembled.insert(beam.id());
               if (reader.Matches(ElementKind::ClosureJoint, beam))
                  bClosureJoint = true;
            }
         }
      }
      if (1 < nBeams || bClosureJoint)
         nSpliced++;
   }

   for (auto closure_joint : reader.Select(ElementKind::ClosureJoint, file))
   {
      if (!assembled.contains(closure_joint.id()))
         nSpliced++;
   }
   return nSpliced;
}

std::vector<SplicedGirder> get_spliced_girders(std::shared_ptr<WBFL::EAF::Broker> pBroker, ifcopenshell::file& file)
{
   const auto& reader = CIfcImporter::GetTargetReader();

   // the parts of each assembly
   struct Parts
   {
      SplicedGirder girder;
      std::vector<int> segment_ids, closure_ids;
   };
   std::vector<Parts> assemblies;
   std::set<int> segment_ids, closure_ids;
   for (auto assembly : reader.Select(ElementKind::GirderAssembly, file))
   {
      Parts parts;
      parts.girder.id = assembly.id();
      parts.girder.name = assembly.Name().value_or("unnamed girder");
      for (auto rel : assembly.IsDecomposedBy())
      {
         for (auto part : rel.RelatedObjects())
         {
            auto beam = part.as<IfcSchema::IfcBeam>();
            if (!beam)
               continue;
            if (reader.Matches(ElementKind::ClosureJoint, beam))
            {
               parts.closure_ids.push_back(beam.id());
               closure_ids.insert(beam.id());
            }
            else
            {
               parts.segment_ids.push_back(beam.id());
               segment_ids.insert(beam.id());
            }
         }
      }
      if (!parts.segment_ids.empty())
         assemblies.push_back(std::move(parts));
   }

   // the lone girders: girders in the superstructure that aren't in an assembly
   std::set<int> lone_ids;
   for (auto id : get_beam_ids(file))
   {
      if (!segment_ids.contains(id))
         lone_ids.insert(id);
   }

   // the geometry of every segment and closure joint at once
   std::set<int> ids(segment_ids);
   ids.insert(closure_ids.begin(), closure_ids.end());
   ids.insert(lone_ids.begin(), lone_ids.end());
   std::map<int, GirderCenterline> located;
   for (auto& centerline : get_beam_centerlines(pBroker, file, ids))
      located.emplace(centerline.beam_id, std::move(centerline));

   std::vector<SplicedGirder> girders;
   for (auto& parts : assemblies)
   {
      auto& girder = parts.girder;
      for (auto id : parts.segment_ids)
      {
         if (auto found = located.find(id); found != located.end())
            girder.segments.push_back(found->second);
         else
            WBFL::System::Logger::Warning("A segment of " + girder.name + " (#" + std::to_string(id) + ") can't be located from its geometry. It isn't imported.");
      }
      if (girder.segments.empty())
         continue;
      std::sort(girder.segments.begin(), girder.segments.end(), [](const auto& a, const auto& b) {return a.start_station < b.start_station; });

      // a closure joint is in the gap its middle is in. The middle of its centerline is the middle of the joint, whichever
      // of its edges get_beam_centerlines took for its ends
      girder.closure_joints.resize(girder.segments.size() - 1);
      for (auto id : parts.closure_ids)
      {
         auto found = located.find(id);
         if (found == located.end())
         {
            WBFL::System::Logger::Warning("A closure joint of " + girder.name + " (#" + std::to_string(id) + ") can't be located from its geometry. The gap it fills is taken as the end of a girder group.");
            continue;
         }
         Float64 station = mid_station(found->second);
         bool bPlaced = false;
         for (size_t i = 0; i + 1 < girder.segments.size(); i++)
         {
            if (girder.segments[i].end_station - station_tolerance <= station && station <= girder.segments[i + 1].start_station + station_tolerance)
            {
               girder.closure_joints[i] = id;
               bPlaced = true;
               break;
            }
         }
         if (!bPlaced)
            WBFL::System::Logger::Warning("A closure joint of " + girder.name + " (#" + std::to_string(id) + ") isn't between two of its segments. It isn't imported.");
      }
      girders.push_back(std::move(girder));
   }

   for (auto id : lone_ids)
   {
      auto found = located.find(id);
      if (found == located.end())
         continue;
      SplicedGirder girder;
      girder.id = id;
      girder.name = file.instance_by_id(id).as<IfcSchema::IfcBeam>().Name().value_or("unnamed beam");
      girder.segments.push_back(found->second);
      girders.push_back(std::move(girder));
   }

   // left to right, then by station
   std::sort(girders.begin(), girders.end(), [](const auto& a, const auto& b)
      {
         Float64 sa = a.segments.front().start_station, sb = b.segments.front().start_station;
         if (station_tolerance < fabs(sa - sb))
            return sa < sb;
         return a.segments.front().start_offset < b.segments.front().start_offset;
      });

   for (const auto& girder : girders)
   {
      IndexType nClosures = std::count_if(girder.closure_joints.begin(), girder.closure_joints.end(), [](const auto& c) {return c.has_value(); });
      std::_tostringstream os;
      os << girder.name.c_str() << _T(": ") << girder.segments.size() << _T(" segments and ") << nClosures << _T(" closure joints, station ")
         << ft(girder.segments.front().start_station) << _T(" to ") << ft(girder.segments.back().end_station);
      WBFL::System::Logger::Info(os.str().c_str());
   }

   return girders;
}

GroupIndexType SplicedFraming::GetGroup(Float64 station, const std::vector<Float64>& pier_stations) const
{
   GroupIndexType grpIdx = 0;
   for (PierIndexType pierIdx = 1; pierIdx + 1 < pier_stations.size(); pierIdx++)
   {
      if (group_boundary[pierIdx] && pier_stations[pierIdx] < station)
         grpIdx++;
   }
   return grpIdx;
}

GroupIndexType SplicedFraming::GetGroupCount() const
{
   return (GroupIndexType)std::count(group_boundary.begin() + 1, group_boundary.end(), true);
}

SplicedFraming get_spliced_framing(ifcopenshell::file& file, const std::vector<SplicedGirder>& girders, const std::vector<Float64>& pier_stations)
{
   SplicedFraming framing;
   PierIndexType nPiers = pier_stations.size();
   framing.group_boundary.resize(nPiers, true);

   // the piers inside the bridge
   for (PierIndexType pierIdx = 1; pierIdx + 1 < nPiers; pierIdx++)
   {
      Float64 station = pier_stations[pierIdx];
      std::optional<AtPier> decided;
      std::_tstring decided_by;
      IndexType nDisagree = 0;
      for (const auto& girder : girders)
      {
         auto at = at_pier(girder, station);
         if (at == AtPier::NotThere)
            continue;
         if (!decided)
         {
            decided = at;
            decided_by = CString(girder.name.c_str()).GetString();
         }
         else if (at != *decided)
         {
            nDisagree++;
         }
      }

      std::_tostringstream os;
      os << _T("Pier ") << LABEL_PIER(pierIdx) << _T(": ");
      if (!decided || *decided == AtPier::Ends)
      {
         os << _T("the girders end here, so it's the end of a girder group");
      }
      else
      {
         framing.group_boundary[pierIdx] = false;
         auto connection = (*decided == AtPier::Continuous ? pgsTypes::psctContinuousSegment : pgsTypes::psctContinousClosureJoint);
         framing.connections[pierIdx] = connection;
         os << _T("the girders ") << describe(*decided) << _T(" over it, so it's inside a girder group with ") << CPierData2::AsString(connection);
      }
      if (decided)
         os << _T(" (") << decided_by << _T(")");
      if (0 < nDisagree)
         os << _T(". ") << nDisagree << _T(" girders don't agree; the leftmost girder is used");
      WBFL::System::Logger::Info(os.str().c_str());
   }

   // temporary supports: the closure joints between the piers (from the leftmost girder of each group), and the erection
   // towers
   auto towers = get_erection_tower_stations(file);
   std::vector<bool> tower_used(towers.size(), false);
   auto find_tower = [&](Float64 from, Float64 to) -> std::optional<size_t>
      {
         for (size_t i = 0; i < towers.size(); i++)
         {
            if (!tower_used[i] && from - station_tolerance <= towers[i] && towers[i] <= to + station_tolerance)
               return i;
         }
         return std::nullopt;
      };

   std::set<GroupIndexType> groups_done;
   for (const auto& girder : girders)
   {
      GroupIndexType grpIdx = framing.GetGroup(mid_station(girder.segments.front()), pier_stations);
      if (!groups_done.insert(grpIdx).second)
         continue; // the leftmost girder of the group decided

      for (size_t i = 0; i + 1 < girder.segments.size(); i++)
      {
         const auto& left = girder.segments[i];
         const auto& right = girder.segments[i + 1];
         Float64 from = left.end_station, to = right.start_station;
         bool bAtPier = std::any_of(pier_stations.begin(), pier_stations.end(), [&](Float64 p) {return from - station_tolerance <= p && p <= to + station_tolerance; });
         if (bAtPier)
            continue; // a pier connection

         if (!girder.closure_joints[i])
         {
            std::_tostringstream os;
            os << CString(girder.name.c_str()).GetString() << _T(": segments ") << (i + 1) << _T(" and ") << (i + 2) << _T(" end between the piers (station ") << ft(from) << _T(" to ") << ft(to)
               << _T(") with no closure joint between them. PGSplice needs a closure joint there; a strongback with a closure joint is used.");
            WBFL::System::Logger::Warning(os.str().c_str());
         }

         SplicedFraming::TemporarySupport ts;
         ts.station = (from + to) / 2;
         ts.connection = pgsTypes::tsctClosureJoint;
         ts.closure_length = (right.start - left.end).norm();
         if (auto tower = find_tower(from, to))
         {
            tower_used[*tower] = true;
            ts.type = pgsTypes::ErectionTower;
         }
         else
         {
            ts.type = pgsTypes::StrongBack;
         }
         framing.temporary_supports.push_back(ts);
      }
   }

   for (size_t i = 0; i < towers.size(); i++)
   {
      if (tower_used[i])
         continue;
      bool bAtPier = std::any_of(pier_stations.begin(), pier_stations.end(), [&](Float64 p) {return fabs(p - towers[i]) <= station_tolerance; });
      if (bAtPier || towers[i] <= pier_stations.front() || pier_stations.back() <= towers[i])
      {
         std::_tostringstream os;
         os << _T("The erection tower at station ") << ft(towers[i]) << _T(" is at a pier or outside the bridge. It isn't imported.");
         WBFL::System::Logger::Warning(os.str().c_str());
         continue;
      }
      framing.temporary_supports.push_back({ towers[i], pgsTypes::ErectionTower, pgsTypes::tsctContinuousSegment, 0.0 });
   }

   std::sort(framing.temporary_supports.begin(), framing.temporary_supports.end(), [](const auto& a, const auto& b) {return a.station < b.station; });
   for (const auto& ts : framing.temporary_supports)
   {
      std::_tostringstream os;
      os << _T("Temporary support at station ") << ft(ts.station) << _T(": ") << (ts.type == pgsTypes::ErectionTower ? _T("erection tower") : _T("strongback"));
      if (ts.connection == pgsTypes::tsctClosureJoint)
         os << _T(" with a closure joint ") << in(ts.closure_length) << _T(" long") << (ts.type == pgsTypes::StrongBack ? _T(" (no erection tower is there, so a strongback is assumed)") : _T(""));
      else
         os << _T(" under a continuous segment");
      WBFL::System::Logger::Info(os.str().c_str());
   }

   return framing;
}

TemporarySupportDefaults get_temporary_support_defaults(const CBridgeDescription2& bridge_desc)
{
   TemporarySupportDefaults defaults;
   const auto* pTimeline = bridge_desc.GetTimelineManager();
   if (0 < bridge_desc.GetTemporarySupportCount())
   {
      // the first temporary support with a closure joint, else the first
      const CTemporarySupportData* pTS = bridge_desc.GetTemporarySupport(0);
      for (SupportIndexType tsIdx = 0; tsIdx < bridge_desc.GetTemporarySupportCount(); tsIdx++)
      {
         if (bridge_desc.GetTemporarySupport(tsIdx)->GetConnectionType() == pgsTypes::tsctClosureJoint)
         {
            pTS = bridge_desc.GetTemporarySupport(tsIdx);
            break;
         }
      }
      pTimeline->GetTempSupportEvents(pTS->GetID(), &defaults.erection, &defaults.removal);
      if (const auto* pClosure = pTS->GetClosureJoint(0))
         defaults.cast_closure = pTimeline->GetCastClosureJointEventIndex(pClosure);
      auto [end_distance, measure] = pTS->GetGirderEndDistance();
      defaults.end_distance = end_distance;
   }

   if (defaults.erection == INVALID_INDEX)
      defaults.erection = pTimeline->GetPierErectionEventIndex(bridge_desc.GetPier(0)->GetID());
   if (defaults.cast_closure == INVALID_INDEX)
      defaults.cast_closure = pTimeline->GetCastDeckEventIndex();
   if (defaults.removal == INVALID_INDEX)
      defaults.removal = pTimeline->GetCastDeckEventIndex();
   return defaults;
}

void apply_spliced_framing(const SplicedFraming& framing, const TemporarySupportDefaults& defaults, CBridgeDescription2& bridge_desc)
{
   for (const auto& [pierIdx, connection] : framing.connections)
   {
      auto* pPier = bridge_desc.GetPier(pierIdx);
      ASSERT(pPier->IsInteriorPier());
      pPier->SetSegmentConnectionType(connection, defaults.cast_closure);
   }

   for (const auto& ts : framing.temporary_supports)
   {
      auto* pNew = new CTemporarySupportData;
      pNew->SetStation(ts.station);
      pNew->SetSupportType(ts.type);
      SupportIndexType tsIdx = bridge_desc.AddTemporarySupport(pNew, defaults.erection, defaults.removal, defaults.cast_closure);
      auto* pTS = bridge_desc.GetTemporarySupport(tsIdx);

      // the girder spacing at the temporary support, as at the start of its group (as PGSplice's temporary support dialog does)
      const auto* pGroup = bridge_desc.GetGirderGroup(pTS->GetSpan());
      const auto* pSpacing = pGroup->GetPier(pgsTypes::metStart)->GetGirderSpacing(pgsTypes::Ahead);
      pTS->GetSegmentSpacing()->SetGirderCount(pGroup->GetGirderCount());
      for (GroupIndexType spacingGrpIdx = 0; spacingGrpIdx < pSpacing->GetSpacingGroupCount(); spacingGrpIdx++)
         pTS->GetSegmentSpacing()->SetGirderSpacing(spacingGrpIdx, pSpacing->GetGirderSpacing(spacingGrpIdx));

      if (ts.connection == pgsTypes::tsctClosureJoint)
      {
         pTS->SetConnectionType(pgsTypes::tsctClosureJoint, defaults.cast_closure);

         // the segment ends are half the closure joint from the CL temporary support, along the girder. The end distance
         // is the template's (where the segment rests on the support), the bearing offset is the rest
         Float64 end_distance = defaults.end_distance;
         pTS->SetGirderEndDistance(end_distance, ConnectionLibraryEntry::EndDistanceMeasurementType::FromBearingAlongGirder);
         pTS->SetBearingOffset(ts.closure_length / 2 + end_distance, ConnectionLibraryEntry::BearingOffsetMeasurementType::AlongGirder);
      }
   }
}

namespace
{
   // girder keys of the girder lines: the group, and the girder index from left to right in the group
   std::vector<CGirderKey> get_girder_keys(const std::vector<SplicedGirder>& girders, const SplicedFraming& framing, const std::vector<Float64>* pPierStations)
   {
      std::vector<CGirderKey> keys(girders.size());
      std::map<GroupIndexType, std::vector<size_t>> by_group;
      for (size_t i = 0; i < girders.size(); i++)
      {
         Float64 station = mid_station(girders[i].segments.front());
         GroupIndexType grpIdx = pPierStations ? framing.GetGroup(station, *pPierStations) : girders[i].segments.front().girder_key.groupIndex;
         by_group[grpIdx].push_back(i);
      }
      for (auto& [grpIdx, indices] : by_group)
      {
         // left to right by the offset of the middle of the girder line
         auto offset = [&](size_t i)
            {
               const auto& s = girders[i].segments;
               return (s.front().start_offset + s.back().end_offset) / 2;
            };
         std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {return offset(a) < offset(b); });
         for (GirderIndexType gdrIdx = 0; gdrIdx < indices.size(); gdrIdx++)
            keys[indices[gdrIdx]] = CGirderKey(grpIdx, gdrIdx);
      }
      return keys;
   }
}

GirderLayout get_spliced_span_layout(std::vector<SplicedGirder>& girders, const SplicedFraming& framing, const std::vector<Float64>& pier_stations)
{
   // the girder keys and segment indices of the segments
   auto keys = get_girder_keys(girders, framing, &pier_stations);
   for (size_t i = 0; i < girders.size(); i++)
   {
      for (SegmentIndexType segIdx = 0; segIdx < girders[i].segments.size(); segIdx++)
      {
         girders[i].segments[segIdx].girder_key = keys[i];
         girders[i].segments[segIdx].segment_index = segIdx;
      }
   }

   SpanIndexType nSpans = pier_stations.size() - 1;
   GirderLayout layout(nSpans);
   for (SpanIndexType spanIdx = 0; spanIdx < nSpans; spanIdx++)
   {
      Float64 start = pier_stations[spanIdx], end = pier_stations[spanIdx + 1];
      std::vector<std::pair<CGirderKey, GirderCenterline>> pieces;
      for (size_t i = 0; i < girders.size(); i++)
      {
         const auto& girder = girders[i];
         const auto& segments = girder.segments;
         if (segments.back().end_station < start + station_tolerance || end - station_tolerance < segments.front().start_station)
            continue; // not in this span

         // the segments in the span
         std::vector<const GirderCenterline*> in_span;
         for (const auto& segment : segments)
         {
            if (start - station_tolerance < segment.end_station && segment.start_station < end + station_tolerance &&
                station_tolerance < std::min(segment.end_station, end) - std::max(segment.start_station, start))
               in_span.push_back(&segment);
         }
         if (in_span.empty())
            continue;

         // the piece of the span takes the segment with the most of it (its beam, mesh, and segment index)
         auto overlap = [&](const GirderCenterline* s) {return std::min(s->end_station, end) - std::max(s->start_station, start); };
         const auto* main = *std::max_element(in_span.begin(), in_span.end(), [&](auto a, auto b) {return overlap(a) < overlap(b); });
         GirderCenterline piece = *main;
         piece.bearings = {};
         piece.start_bearing.reset();
         piece.end_bearing.reset();

         // the start: on the first segment in the span, at the pier if it continues over it
         const auto* first = in_span.front();
         if (first->start_station < start - station_tolerance)
         {
            at_station(*first, start, piece.start, piece.start_top, piece.start_offset);
            piece.start_station = start;
            piece.continuous[pgsTypes::metStart] = true;
         }
         else
         {
            piece.start = first->start;
            piece.start_top = first->start_top;
            piece.start_station = first->start_station;
            piece.start_offset = first->start_offset;
            piece.continuous[pgsTypes::metStart] = false;
         }

         const auto* last = in_span.back();
         if (end + station_tolerance < last->end_station)
         {
            at_station(*last, end, piece.end, piece.end_top, piece.end_offset);
            piece.end_station = end;
            piece.continuous[pgsTypes::metEnd] = true;
         }
         else
         {
            piece.end = last->end;
            piece.end_top = last->end_top;
            piece.end_station = last->end_station;
            piece.end_offset = last->end_offset;
            piece.continuous[pgsTypes::metEnd] = false;
         }

         pieces.emplace_back(keys[i], std::move(piece));
      }

      std::sort(pieces.begin(), pieces.end(), [](const auto& a, const auto& b) {return a.first.girderIndex < b.first.girderIndex; });
      for (auto& [key, piece] : pieces)
         layout[spanIdx].push_back(std::move(piece));
   }

   return layout;
}

GirderLayout get_spliced_segment_layout(const std::vector<SplicedGirder>& girders, const SplicedFraming& framing)
{
   GirderLayout layout(framing.GetGroupCount());
   for (const auto& girder : girders)
   {
      for (const auto& segment : girder.segments)
      {
         if (segment.girder_key.groupIndex < layout.size())
            layout[segment.girder_key.groupIndex].push_back(segment);
      }
   }
   return layout;
}

GirderLayout get_spliced_girder_layout(ifcopenshell::file& file, const GirderProfileReader& reader, const std::vector<SplicedGirder>& girders, const SplicedFraming& framing)
{
   GirderLayout layout(framing.GetGroupCount());
   for (const auto& girder : girders)
   {
      // the segment whose height changes the least
      const GirderCenterline* best = &girder.segments.front();
      Float64 best_range = DBL_MAX;
      for (const auto& segment : girder.segments)
      {
         auto sections = reader.ReadCrossSections(file.instance_by_id(segment.beam_id).as<IfcSchema::IfcBeam>());
         if (sections.empty())
            continue;
         auto [lo, hi] = std::minmax_element(sections.begin(), sections.end(), [](const auto& a, const auto& b) {return a.height < b.height; });
         Float64 range = hi->height - lo->height;
         if (range < best_range - WIDTH_PROFILE_STEP)
         {
            best_range = range;
            best = &segment;
         }
      }
      auto key = best->girder_key;
      if (key.groupIndex < layout.size())
      {
         auto& girders_in_group = layout[key.groupIndex];
         if (girders_in_group.size() <= key.girderIndex)
            girders_in_group.resize(key.girderIndex + 1);
         girders_in_group[key.girderIndex] = *best;
      }
   }
   return layout;
}

void set_segment_variations(ifcopenshell::file& file, const GirderProfileReader& reader, const GirderLayout& segment_layout, const std::vector<Float64>& pier_stations, CBridgeDescription2& bridge_desc)
{
   // the segments of each girder in order, for the distance from each segment end to its framing line: half the gap to the
   // next segment (the CL of the closure joint), or the distance to the nearest pier at the end of the girder
   std::map<CGirderKey, std::map<SegmentIndexType, const GirderCenterline*>> by_girder;
   for (const auto& segments : segment_layout)
      for (const auto& segment : segments)
         by_girder[segment.girder_key][segment.segment_index] = &segment;
   auto framing_offsets = [&](const GirderCenterline& segment)
      {
         const auto& girder = by_girder[segment.girder_key];
         auto nearest_pier = [&](Float64 station)
            {
               Float64 d = DBL_MAX;
               for (auto p : pier_stations)
                  d = std::min(d, fabs(p - station));
               return d;
            };
         auto prev = girder.find(segment.segment_index - 1);
         auto next = girder.find(segment.segment_index + 1);
         Float64 start = (0 < segment.segment_index && prev != girder.end()) ? (segment.start - prev->second->end).norm() / 2 : nearest_pier(segment.start_station);
         Float64 end = (next != girder.end()) ? (next->second->start - segment.end).norm() / 2 : nearest_pier(segment.end_station);
         return std::make_pair(start, end);
      };

   const Float64 tolerance = WIDTH_PROFILE_STEP; // 1/8 in
   std::map<pgsTypes::SegmentVariationType, IndexType> counts;
   IndexType nNotRead = 0;
   for (const auto& segments : segment_layout)
   {
      for (const auto& segment : segments)
      {
         auto beam = file.instance_by_id(segment.beam_id).as<IfcSchema::IfcBeam>();
         auto sections = reader.ReadCrossSections(beam);
         auto* pSegment = bridge_desc.GetGirderGroup(segment.girder_key.groupIndex)->GetGirder(segment.girder_key.girderIndex)->GetSegment(segment.segment_index);
         if (sections.size() < 2)
         {
            nNotRead++;
            continue;
         }

         // distances along the segment: the cross sections are along the directrix (the alignment), so they're scaled to
         // the segment's length
         Float64 L = (segment.end - segment.start).norm();
         Float64 scale = IsZero(sections.back().distance) ? 1.0 : L / sections.back().distance;
         auto [start_offset, end_offset] = framing_offsets(segment);
         size_t n = sections.size();
         auto x = [&](size_t i) {return sections[i].distance * scale; };
         auto h = [&](size_t i) {return sections[i].height; };
         auto bf = [&](size_t i) {return sections[i].bottom_flange_depth; };

         bool bConstant = std::all_of(sections.begin(), sections.end(), [&](const auto& s) {return fabs(s.height - h(0)) <= tolerance && fabs(s.bottom_flange_depth - bf(0)) <= tolerance; });
         if (bConstant)
         {
            pSegment->SetVariationType(pgsTypes::svtNone);
            counts[pgsTypes::svtNone]++;
            continue;
         }

         // the prismatic ends
         size_t i0 = 0;
         while (i0 + 1 < n && fabs(h(i0 + 1) - h(0)) <= tolerance)
            i0++;
         size_t i1 = n - 1;
         while (0 < i1 && fabs(h(i1 - 1) - h(n - 1)) <= tolerance)
            i1--;

         // a deeper (or shallower) part inside the segment makes it a double variation (e.g. a pier segment)
         size_t k = i0;
         for (size_t i = i0; i <= i1; i++)
         {
            if (h(k) < h(i))
               k = i;
         }
         Float64 hmax = h(k);
         bool bDouble = (std::max(h(0), h(n - 1)) + tolerance < hmax) && i0 < k && k < i1;

         bool bVariableBottomFlange = std::any_of(sections.begin(), sections.end(), [&](const auto& s) {return tolerance < fabs(s.bottom_flange_depth - bf(0)); });
         pgsTypes::SegmentVariationType type;
         if (bDouble)
         {
            // the plateau at the deepest part
            size_t a = k, b = k;
            while (i0 < a - 1 && fabs(h(a - 1) - hmax) <= tolerance)
               a--;
            while (b + 1 < i1 && fabs(h(b + 1) - hmax) <= tolerance)
               b++;
            bool bParabolic = 1 < (a - i0) || 1 < (i1 - b); // sections inside a taper: a curve
            type = bParabolic ? pgsTypes::svtDoubleParabolic : pgsTypes::svtDoubleLinear;
            pSegment->SetVariationType(type);
            pSegment->SetVariationParameters(pgsTypes::sztLeftPrismatic, start_offset + x(i0), h(0), bf(0));
            pSegment->SetVariationParameters(pgsTypes::sztLeftTapered, x(a) - x(i0), hmax, bf(a));
            pSegment->SetVariationParameters(pgsTypes::sztRightTapered, x(i1) - x(b), hmax, bf(b));
            pSegment->SetVariationParameters(pgsTypes::sztRightPrismatic, L - x(i1) + end_offset, h(n - 1), bf(n - 1));
         }
         else
         {
            bool bParabolic = 1 < (i1 - i0);
            type = bParabolic ? pgsTypes::svtParabolic : pgsTypes::svtLinear;
            pSegment->SetVariationType(type);
            pSegment->SetVariationParameters(pgsTypes::sztLeftPrismatic, start_offset + x(i0), h(0), bf(0));
            pSegment->SetVariationParameters(pgsTypes::sztRightPrismatic, L - x(i1) + end_offset, h(n - 1), bf(n - 1));
         }
         pSegment->EnableVariableBottomFlangeDepth(bVariableBottomFlange);
         counts[type]++;

         std::_tostringstream os;
         os << CString(beam.Name().value_or("unnamed segment").c_str()).GetString() << _T(" of ") << LABEL_GIRDER(segment.girder_key) << _T(", group ") << LABEL_GROUP(segment.girder_key.groupIndex)
            << _T(": ") << CPrecastSegmentData::GetSegmentVariation(type) << _T(" depth, ") << in(h(0)) << _T(" to ") << in(hmax) << (bDouble ? _T(" to ") + in(h(n - 1)) : _T(""))
            << (bVariableBottomFlange ? _T(", with a variable bottom flange") : _T("")) << _T(" (") << n << _T(" cross sections)");
         WBFL::System::Logger::Info(os.str().c_str());
      }
   }

   std::_tostringstream os;
   os << _T("Segment depth variation:");
   for (const auto& [type, count] : counts)
      os << _T(" ") << count << _T(" ") << CPrecastSegmentData::GetSegmentVariation(type);
   if (0 < nNotRead)
      os << _T(", ") << nNotRead << _T(" segments without cross sections to read (the template's variation is kept)");
   WBFL::System::Logger::Info(os.str().c_str());
}
