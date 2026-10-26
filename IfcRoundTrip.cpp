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
#include "IfcRoundTrip.h"

#include <PsgLib/BridgeDescription2.h>
#include <System/IStructuredSave.h>
#include <System/IStructuredLoad.h>

#include <array>
#include <bcrypt.h>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace roundtrip_key
{
   std::string Project() { return "project"; }
   std::string Site() { return "site"; }
   std::string Bridge() { return "bridge"; }
   std::string Superstructure() { return "superstructure"; }
   std::string Substructure() { return "substructure"; }
   std::string DeckPart() { return "deck_part"; }
   std::string Deck() { return "deck"; }
   std::string Barrier(pgsTypes::TrafficBarrierOrientation side) { return side == pgsTypes::tboLeft ? "barrier:left" : "barrier:right"; }
   std::string Pier(PierIDType pierID) { return "pier:" + std::to_string(pierID); }
   std::string Girder(SegmentIDType segmentID) { return "girder:" + std::to_string(segmentID); }
   std::string Bearings(SegmentIDType segmentID, pgsTypes::MemberEndType end) { return "bearings:" + std::to_string(segmentID) + (end == pgsTypes::metStart ? ":start" : ":end"); }
   std::string Diaphragms(PierIDType pierID) { return "diaphragms:" + std::to_string(pierID); }
   std::string Tendons(SegmentIDType segmentID) { return "tendons:" + std::to_string(segmentID); }
   std::string Bars(SegmentIDType segmentID) { return "bars:" + std::to_string(segmentID); }
};

std::wstring ComputeSha256(const std::filesystem::path& path)
{
   std::ifstream file(path, std::ios::binary);
   if (!file)
      return {};

   BCRYPT_ALG_HANDLE hAlg = nullptr;
   BCRYPT_HASH_HANDLE hHash = nullptr;
   if (!BCRYPT_SUCCESS(BCryptOpenAlgorithmProvider(&hAlg, BCRYPT_SHA256_ALGORITHM, nullptr, 0)))
      return {};

   std::wstring result;
   if (BCRYPT_SUCCESS(BCryptCreateHash(hAlg, &hHash, nullptr, 0, nullptr, 0, 0)))
   {
      std::vector<char> buffer(1 << 20);
      bool bOK = true;
      while (bOK && file)
      {
         file.read(buffer.data(), buffer.size());
         auto n = file.gcount();
         if (0 < n)
            bOK = BCRYPT_SUCCESS(BCryptHashData(hHash, (PUCHAR)buffer.data(), (ULONG)n, 0));
      }

      std::array<UCHAR, 32> digest;
      if (bOK && BCRYPT_SUCCESS(BCryptFinishHash(hHash, digest.data(), (ULONG)digest.size(), 0)))
      {
         std::wostringstream os;
         for (auto byte : digest)
            os << std::hex << std::setw(2) << std::setfill(L'0') << (int)byte;
         result = os.str();
      }
      BCryptDestroyHash(hHash);
   }
   BCryptCloseAlgorithmProvider(hAlg, 0);
   return result;
}

SourceModel DescribeSourceModel(const std::filesystem::path& path)
{
   SourceModel source;
   std::error_code ec;
   auto absolute = std::filesystem::absolute(path, ec);
   source.path = (ec ? path : absolute).wstring();
   source.size = std::filesystem::file_size(path, ec);
   if (ec)
      source.size = 0;
   source.sha256 = ComputeSha256(path);
   return source;
}

SourceModelState CheckSourceModel(const SourceModel& source, const std::filesystem::path& path)
{
   std::error_code ec;
   if (!std::filesystem::exists(path, ec))
      return SourceModelState::Missing;

   auto size = std::filesystem::file_size(path, ec);
   if (ec || size != source.size)
      return SourceModelState::Changed;

   return ComputeSha256(path) == source.sha256 ? SourceModelState::Unchanged : SourceModelState::Changed;
}

const std::vector<std::string>& CRoundTripRegistry::Find(const std::string& key) const
{
   static const std::vector<std::string> none;
   auto found = elements.find(key);
   return found == elements.end() ? none : found->second;
}

void CRoundTripRegistry::Save(WBFL::System::IStructuredSave* pStrSave) const
{
   USES_CONVERSION;
   pStrSave->BeginUnit(_T("RoundTrip"), 1.0);
   pStrSave->Property(_T("SourceModel"), source.path.c_str());
   pStrSave->Property(_T("SourceModelSize"), source.size);
   pStrSave->Property(_T("SourceModelSha256"), source.sha256.c_str());
   pStrSave->Property(_T("ElementCount"), (Uint64)elements.size());
   for (const auto& [key, global_ids] : elements)
   {
      pStrSave->BeginUnit(_T("Element"), 1.0);
      pStrSave->Property(_T("Key"), A2T(key.c_str()));
      pStrSave->Property(_T("Count"), (Uint64)global_ids.size());
      for (const auto& global_id : global_ids)
         pStrSave->Property(_T("GlobalId"), A2T(global_id.c_str()));
      pStrSave->EndUnit();
   }
   pStrSave->EndUnit();
}

void CRoundTripRegistry::Load(WBFL::System::IStructuredLoad* pStrLoad)
{
   USES_CONVERSION;
   auto fail = []() { throw std::runtime_error("The IFC round trip data of the project can't be read."); };

   elements.clear();
   if (!pStrLoad->BeginUnit(_T("RoundTrip"))) fail();
   std::_tstring text;
   if (!pStrLoad->Property(_T("SourceModel"), &text)) fail();
   source.path = text;
   if (!pStrLoad->Property(_T("SourceModelSize"), &source.size)) fail();
   if (!pStrLoad->Property(_T("SourceModelSha256"), &text)) fail();
   source.sha256 = text;
   Uint64 nElements;
   if (!pStrLoad->Property(_T("ElementCount"), &nElements)) fail();
   for (Uint64 i = 0; i < nElements; i++)
   {
      if (!pStrLoad->BeginUnit(_T("Element"))) fail();
      if (!pStrLoad->Property(_T("Key"), &text)) fail();
      auto& global_ids = elements[std::string(T2A(text.c_str()))];
      Uint64 nGlobalIds;
      if (!pStrLoad->Property(_T("Count"), &nGlobalIds)) fail();
      for (Uint64 j = 0; j < nGlobalIds; j++)
      {
         if (!pStrLoad->Property(_T("GlobalId"), &text)) fail();
         global_ids.push_back(std::string(T2A(text.c_str())));
      }
      if (!pStrLoad->EndUnit()) fail();
   }
   if (!pStrLoad->EndUnit()) fail();
}

void CRoundTripRecorder::Record(const std::string& key, int entity_id)
{
   m_Entries.push_back({ Kind::Key, key, CSegmentKey(), INVALID_INDEX, pgsTypes::metStart, entity_id });
}

void CRoundTripRecorder::RecordPier(PierIndexType pierIdx, int entity_id)
{
   m_Entries.push_back({ Kind::Pier, {}, CSegmentKey(), pierIdx, pgsTypes::metStart, entity_id });
}

void CRoundTripRecorder::RecordGirder(const CSegmentKey& segmentKey, int entity_id)
{
   m_Entries.push_back({ Kind::Girder, {}, segmentKey, INVALID_INDEX, pgsTypes::metStart, entity_id });
}

void CRoundTripRecorder::RecordBearing(const CSegmentKey& segmentKey, pgsTypes::MemberEndType end, int entity_id)
{
   m_Entries.push_back({ Kind::Bearing, {}, segmentKey, INVALID_INDEX, end, entity_id });
}

void CRoundTripRecorder::RecordDiaphragm(PierIndexType pierIdx, int entity_id)
{
   m_Entries.push_back({ Kind::Diaphragm, {}, CSegmentKey(), pierIdx, pgsTypes::metStart, entity_id });
}

void CRoundTripRecorder::RecordTendons(const CSegmentKey& segmentKey, int entity_id)
{
   m_Entries.push_back({ Kind::Tendons, {}, segmentKey, INVALID_INDEX, pgsTypes::metStart, entity_id });
}

void CRoundTripRecorder::RecordBars(const CSegmentKey& segmentKey, int entity_id)
{
   m_Entries.push_back({ Kind::Bars, {}, segmentKey, INVALID_INDEX, pgsTypes::metStart, entity_id });
}

CRoundTripRegistry CRoundTripRecorder::MakeRegistry(ifcopenshell::file& file, const CBridgeDescription2& bridge_desc, const SourceModel& source) const
{
   CRoundTripRegistry registry;
   registry.source = source;

   auto segment_id = [&bridge_desc](const CSegmentKey& segmentKey) -> std::optional<SegmentIDType>
      {
         if (bridge_desc.GetGirderGroupCount() <= segmentKey.groupIndex)
            return std::nullopt;
         const auto* pGroup = bridge_desc.GetGirderGroup(segmentKey.groupIndex);
         if (pGroup->GetGirderCount() <= segmentKey.girderIndex)
            return std::nullopt;
         const auto* pGirder = pGroup->GetGirder(segmentKey.girderIndex);
         if (pGirder->GetSegmentCount() <= segmentKey.segmentIndex)
            return std::nullopt;
         return pGirder->GetSegment(segmentKey.segmentIndex)->GetID();
      };
   auto pier_id = [&bridge_desc](PierIndexType pierIdx) -> std::optional<PierIDType>
      {
         if (bridge_desc.GetPierCount() <= pierIdx)
            return std::nullopt;
         return bridge_desc.GetPier(pierIdx)->GetID();
      };

   for (const auto& entry : m_Entries)
   {
      std::string key;
      switch (entry.kind)
      {
      case Kind::Key: key = entry.key; break;
      case Kind::Pier: if (auto id = pier_id(entry.pier)) key = roundtrip_key::Pier(*id); break;
      case Kind::Diaphragm: if (auto id = pier_id(entry.pier)) key = roundtrip_key::Diaphragms(*id); break;
      case Kind::Girder: if (auto id = segment_id(entry.segment)) key = roundtrip_key::Girder(*id); break;
      case Kind::Bearing: if (auto id = segment_id(entry.segment)) key = roundtrip_key::Bearings(*id, entry.end); break;
      case Kind::Tendons: if (auto id = segment_id(entry.segment)) key = roundtrip_key::Tendons(*id); break;
      case Kind::Bars: if (auto id = segment_id(entry.segment)) key = roundtrip_key::Bars(*id); break;
      }
      if (key.empty())
         continue; // the item isn't in the bridge (e.g. a span with too few girders kept the default)

      auto instance = file.instance_by_id(entry.entity_id);
      auto root = instance.as<IfcSchema::IfcRoot>();
      if (!root)
         continue;

      auto& global_ids = registry.elements[key];
      if (std::find(global_ids.begin(), global_ids.end(), root.GlobalId()) == global_ids.end())
         global_ids.push_back(root.GlobalId());
   }
   return registry;
}
