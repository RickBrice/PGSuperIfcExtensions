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

#pragma once

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace WBFL { namespace EAF { class Broker; }; };
namespace WBFL { namespace System { class IStructuredSave; class IStructuredLoad; }; };
class CBridgeDescription2;

// Round trips (devdocs/RoundTripDesign.md): an IFC model imported into PGSuper is exported again as a copy of the model
// (the source model) with PGSuper's design in it. The source model isn't changed.
//
// The element registry ties the two together. The import records which elements of the source model became which
// PGSuper items (by PGSuper's stable IDs, so adding or removing a span or girder doesn't shift the others), and the
// project saves it. The export names the same items (IfcExportIdentities), and the merge (IfcRoundTripMerge.h) updates
// the source's elements in place and adds what PGSuper designs that the source doesn't have.

// The keys of the registry: an element role and PGSuper's ID of the item
namespace roundtrip_key
{
   std::string Project();
   std::string Site();
   std::string Bridge();
   std::string Superstructure();
   std::string Substructure();
   std::string DeckPart();
   std::string Deck();
   std::string Barrier(pgsTypes::TrafficBarrierOrientation side);
   std::string Pier(PierIDType pierID);                                     // IfcBridgePart of an abutment or pier
   std::string Girder(SegmentIDType segmentID);                             // IfcBeam of a segment
   std::string Bearings(SegmentIDType segmentID, pgsTypes::MemberEndType end); // IfcBearing at a segment end, in order
   std::string Diaphragms(PierIDType pierID);                               // pier diaphragms
   std::string Tendons(SegmentIDType segmentID);                            // the source has the segment's strands
   std::string Bars(SegmentIDType segmentID);                               // the source has the segment's reinforcing bars
};

// The model an imported project came from
struct SourceModel
{
   std::wstring path;
   Uint64 size = 0;
   std::wstring sha256;

   bool IsEmpty() const { return path.empty(); }
};

// The SHA-256 of a file, in lower case hex. Empty if the file can't be read
std::wstring ComputeSha256(const std::filesystem::path& path);

// The source model information of a file (path, size, SHA-256)
SourceModel DescribeSourceModel(const std::filesystem::path& path);

// Whether a file is the source model: the same size and SHA-256
enum class SourceModelState { Missing, Changed, Unchanged };
SourceModelState CheckSourceModel(const SourceModel& source, const std::filesystem::path& path);

class CRoundTripRegistry
{
public:
   SourceModel source;
   std::map<std::string, std::vector<std::string>> elements; // registry key -> GlobalIds of the source model's elements

   bool IsEmpty() const { return source.IsEmpty(); }

   // the GlobalIds of a key (empty if the source has none)
   const std::vector<std::string>& Find(const std::string& key) const;

   void Save(WBFL::System::IStructuredSave* pStrSave) const;
   void Load(WBFL::System::IStructuredLoad* pStrLoad); // throws if the data can't be read
};

// The registry of a project, saved with it (CIfcExtensionAgent)
// {9B1E5C47-3D7A-4F2B-8E61-2C4A7D90B3F5}
DEFINE_GUID(IID_IIfcRoundTrip,
   0x9b1e5c47, 0x3d7a, 0x4f2b, 0x8e, 0x61, 0x2c, 0x4a, 0x7d, 0x90, 0xb3, 0xf5);
class IIfcRoundTrip
{
public:
   virtual const CRoundTripRegistry& GetRoundTripRegistry() const = 0;
   virtual void SetRoundTripRegistry(const CRoundTripRegistry& registry) = 0;
};

// Collects the elements of the source model the import uses, by PGSuper index (IDs are known once the bridge description is
// set). The import adds to it as it reads the model (CIfcImporter::GetRoundTripRecorder) and makes the registry at the end
class CRoundTripRecorder
{
public:
   void Clear() { m_Entries.clear(); }

   void Record(const std::string& key, int entity_id); // an element with no PGSuper index (project, bridge, deck, ...)
   void RecordPier(PierIndexType pierIdx, int entity_id);
   void RecordGirder(const CSegmentKey& segmentKey, int entity_id);
   void RecordBearing(const CSegmentKey& segmentKey, pgsTypes::MemberEndType end, int entity_id);
   void RecordDiaphragm(PierIndexType pierIdx, int entity_id);
   void RecordTendons(const CSegmentKey& segmentKey, int entity_id); // entity_id: the girder
   void RecordBars(const CSegmentKey& segmentKey, int entity_id);

   // The registry: the keys with PGSuper's IDs, and the GlobalIds of the elements
   CRoundTripRegistry MakeRegistry(ifcopenshell::file& file, const CBridgeDescription2& bridge_desc, const SourceModel& source) const;

private:
   enum class Kind { Key, Pier, Girder, Bearing, Diaphragm, Tendons, Bars };
   struct Entry
   {
      Kind kind;
      std::string key;
      CSegmentKey segment;
      PierIndexType pier = INVALID_INDEX;
      pgsTypes::MemberEndType end = pgsTypes::metStart;
      int entity_id;
   };
   std::vector<Entry> m_Entries;
};

// The elements of a PGSuper export, by registry key (entity ids in the exported model). The exporter fills it when
// CIfcExportOptions::export_identities is set
using IfcExportIdentities = std::map<std::string, std::vector<int>>;
