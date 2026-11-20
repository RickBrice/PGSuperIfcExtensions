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
#include "TestExportExtension.h"

#include <sstream>

namespace
{
   constexpr const char* MarkerRole = "ifctest.marker";

   // Records every hook call, in order. The marker element's property set has them (ifctest.marker.hooks)
   class CTestHooks : public IfcExport::IHooks
   {
   public:
      explicit CTestHooks(CTestExportExtension::Mode mode) : m_Mode(mode) {}

      void OnExportBegin(IfcExport::HookContext& context) override
      {
         Record("OnExportBegin");
         Require(context.elements.project && context.elements.site, "OnExportBegin: the project and the site exist");
      }

      void OnDeck(IfcExport::HookContext& context, IfcExport::Schema::IfcSlab slab) override
      {
         Record("OnDeck");
         Require((bool)slab, "OnDeck: the slab exists");
      }

      void OnBridgeParts(IfcExport::HookContext& context) override
      {
         Record("OnBridgeParts");
         Require(context.elements.bridge && context.elements.superstructure && context.elements.substructure, "OnBridgeParts: the bridge and its parts exist");
      }

      void OnPier(IfcExport::HookContext& context, PierIndexType pierIdx, IfcExport::Schema::IfcBridgePart pier) override
      {
         if (m_Mode == CTestExportExtension::Mode::Throw)
            throw std::runtime_error("the test extension fails here on purpose");

         Record("OnPier");
         Require((bool)pier, "OnPier: the pier exists");
      }

      void OnBearing(IfcExport::HookContext& context, PierIndexType pierIdx, pgsTypes::PierFaceType face, const CGirderKey& girderKey, IfcExport::Schema::IfcBearing bearing) override
      {
         Record("OnBearing");
      }

      void OnBarrier(IfcExport::HookContext& context, pgsTypes::TrafficBarrierOrientation side, IfcExport::Schema::IfcProduct barrier) override
      {
         Record("OnBarrier");
      }

      void OnSegment(IfcExport::HookContext& context, const CSegmentKey& segmentKey, IfcExport::Schema::IfcBeam beam) override
      {
         Record("OnSegment");
         Require(context.elements.segments.find(segmentKey) != context.elements.segments.end(), "OnSegment: the segment is in the exported elements");
      }

      void OnSuperstructureComplete(IfcExport::HookContext& context) override
      {
         Record("OnSuperstructureComplete");
         std::ostringstream os;
         os << "Test extension: " << context.elements.segments.size() << " segments, " << context.elements.piers.size() << " piers and abutments, "
            << context.elements.barriers.size() << " barriers in the exported elements";
         context.services.Log(os.str());
      }

      void OnBridgeComplete(IfcExport::HookContext& context) override
      {
         Record("OnBridgeComplete");
      }

      void OnModelComplete(IfcExport::HookContext& context) override
      {
         Record("OnModelComplete");

         // the marker element, with the calls in its property set
         using Schema = IfcExport::Schema;
         auto role = context.services.Role(MarkerRole);
         auto marker = context.file.create<Schema::IfcBuildingElementProxy>().initialize(
            ifcopenshell::global_id(),
            {}, // OwnerHistory
            std::string("Test extension marker"), // Name
            std::string("Written by the IFC extension's test export extension (/IfcTestExtension)"), // Description
            std::string("TEST_MARKER"), // ObjectType
            {}, // ObjectPlacement
            {}, // Representation
            std::nullopt, // Tag
            Schema::IfcBuildingElementProxyTypeEnum::IfcBuildingElementProxyType_USERDEFINED);
         context.file.addRelatedObject<Schema::IfcRelContainedInSpatialStructure>(context.elements.site, marker);

         std::string calls;
         for (const auto& call : m_Calls)
            calls += (calls.empty() ? "" : ";") + call;

         IfcExport::ElementKeys keys;
         keys.index = m_Calls.size();
         keys.extension = calls;
         context.services.Classify(role, marker);
         context.services.WritePropertySets(role, { marker }, keys);

         context.services.Log("Test extension: " + std::to_string(m_Calls.size()) + " hook calls before OnExportEnd");
         for (const auto& problem : m_Problems)
            context.services.Log("Test extension problem: " + problem);
      }

   private:
      CTestExportExtension::Mode m_Mode;
      std::vector<std::string> m_Calls;
      std::vector<std::string> m_Problems;

      void Record(const char* call) { m_Calls.push_back(call); }

      void Require(bool bCondition, const char* what)
      {
         if (!bCondition)
            m_Problems.push_back(std::string("expected ") + what);
      }
   };
}

bool CTestExportExtension::GetMode(const CString& name, Mode& mode)
{
   if (name.IsEmpty() || name.CompareNoCase(_T("hooks")) == 0)
      mode = Mode::Hooks;
   else if (name.CompareNoCase(_T("throw")) == 0)
      mode = Mode::Throw;
   else if (name.CompareNoCase(_T("api")) == 0)
      mode = Mode::ApiVersion;
   else if (name.CompareNoCase(_T("fingerprint")) == 0)
      mode = Mode::Fingerprint;
   else
      return false;
   return true;
}

std::string CTestExportExtension::GetId() const
{
   return "PGSuperIfc.TestExtension";
}

std::string CTestExportExtension::GetVersion() const
{
   return "1";
}

std::string CTestExportExtension::GetPrefix() const
{
   return "ifctest";
}

int CTestExportExtension::GetApiVersion() const
{
   return m_Mode == Mode::ApiVersion ? IFC_EXPORT_EXTENSION_API_VERSION + 1 : IIfcExportExtension::GetApiVersion();
}

std::string CTestExportExtension::GetBuildFingerprint() const
{
   return m_Mode == Mode::Fingerprint ? IIfcExportExtension::GetBuildFingerprint() + ";test" : IIfcExportExtension::GetBuildFingerprint();
}

void CTestExportExtension::RegisterContent(IfcExport::IRegistry& registry)
{
   registry.RegisterRole({ MarkerRole, "the test extension's marker element" });

   auto add = [&registry](const char* name, const char* role, ValueKind kind, const char* description, IfcExport::Getter get)
      {
         IfcExport::TargetSpec target;
         target.name = name;
         target.role = role;
         target.kind = kind;
         target.description = description;
         target.get = std::move(get);
         registry.RegisterTarget(target);
      };

   // the extension's own element: the keys its hook gave
   add("ifctest.marker.hooks", MarkerRole, ValueKind::Text, "the hook calls, in order",
      [](const IfcExport::Element& e) { return ExportValue(std::any_cast<std::string>(e.Extension())); });
   add("ifctest.marker.hook_count", MarkerRole, ValueKind::Count, "the number of hook calls",
      [](const IfcExport::Element& e) { return ExportValue((Int64)e.Index()); });

   // the IFC extension's elements: the keys the IFC extension gives
   add("ifctest.girder.segment", "girder", ValueKind::Text, "the girder's segment key",
      [](const IfcExport::Element& e)
      {
         const auto& key = e.Segment();
         return ExportValue(std::to_string(key.groupIndex) + "," + std::to_string(key.girderIndex) + "," + std::to_string(key.segmentIndex));
      });
   add("ifctest.pier.index", "pier", ValueKind::Count, "the pier's index",
      [](const IfcExport::Element& e) { return ExportValue((Int64)e.Pier()); });
}

std::unique_ptr<IfcExport::IHooks> CTestExportExtension::CreateHooks()
{
   return std::make_unique<CTestHooks>(m_Mode);
}
