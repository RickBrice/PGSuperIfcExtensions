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

/*****************************************************************************
   The IFC extension's model builder, for other plug-ins (devdocs/ExtensionsDesign.md)

   Another plug-in (e.g. an exporter for a particular receiving application)
   builds the IFC model the way the IFC export does, and adds to it, without the
   IFC extension knowing about the additions:

   - its own mapping table, which extends the standard table, declares the
     property sets it adds
   - its own targets (IIfcExtensionRegistry::RegisterTarget) supply their values.
     They can belong to the IFC extension's element roles (e.g. "girder") or to
     element roles of its own (IIfcExtensionRegistry::RegisterRole)
   - its hooks (IIfcExportHooks) add elements at fixed points of the build, and
     write their property sets through IIfcExportServices

   The IFC extension builds the model in a fixed order and calls the hooks at
   each step (a template method). The standard IFC export never loads an
   extension: its targets, roles, and tables are only known during
   IIfcModelBuilder::BuildModel calls that name it.

   The extension shares C++ objects (e.g. the IfcOpenShell file) with the IFC
   extension, so both must be built with the same compiler, run time library,
   and IfcOpenShell build. Import IfcExtensionSdk.props in both projects. The
   IFC extension checks the API version and build fingerprint before it uses an
   extension.
*****************************************************************************/

#include "IfcExtensionTypes.h"
#include "IfcExportOptions.h"

#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include <ifcparse/hierarchy_helper.h>
#include <ifcparse/schemas/Ifc4x3_add2.h>

class IIfcExportExtension;

namespace IfcExport
{
   // The schema of the models extensions add to (the IFC extension's only export schema)
   using Schema = ::Ifc4x3_add2;
   using File = hierarchy_helper<Schema>;

   // Identifies the build of a module that uses the SDK. Computed where it's compiled, so the IFC extension and an export
   // extension each give their own
   inline std::string BuildFingerprint()
   {
      std::ostringstream os;
      os << "api=" << IFC_EXPORT_EXTENSION_API_VERSION
         << ";msc=" << _MSC_VER
#if defined _ITERATOR_DEBUG_LEVEL
         << ";idl=" << _ITERATOR_DEBUG_LEVEL
#endif
#if defined _DEBUG
         << ";debug"
#endif
         << ";file=" << sizeof(File)
         << ";context=" << sizeof(ExportContext)
         << ";target=" << sizeof(TargetDef)
         << ";options=" << sizeof(CIfcExportOptions);
      return os.str();
   }

   // An element role of an export extension. Mapping tables name it in "applies_to" and "elements"
   struct RoleSpec
   {
      std::string name;        // starts with the extension's prefix (e.g. "brdr.diaphragm")
      std::string description; // for messages
   };

   // A target of an export extension. Export only: tables can't name it in "targets" (import locations)
   struct TargetSpec
   {
      std::string name;        // starts with the extension's prefix (e.g. "brdr.girder.k1")
      std::string role;        // the element role whose elements it belongs to: an IFC extension role (e.g. "girder") or one the extension registered
      ValueKind kind = ValueKind::Real;
      std::string description; // for messages
      ExportUnit display_unit = ExportUnit::None; // when properties are exported in display units
      Float64 display_round = 0;
      ExportGetter get;
   };

   // Receives an export extension's element roles and targets. Implemented by the IFC extension
   class IRegistry
   {
   public:
      // Registers an element role and returns its ElementKind (FirstExtension or more)
      virtual ElementKind RegisterRole(const RoleSpec& role) = 0;
      virtual void RegisterTarget(const TargetSpec& target) = 0;

      virtual ~IRegistry() = default;
   };

   // The elements the IFC extension has created so far
   struct ExportedElements
   {
      Schema::IfcProject project;
      Schema::IfcSite site;
      Schema::IfcBridge bridge;
      Schema::IfcBridgePart superstructure;
      Schema::IfcBridgePart substructure;
      Schema::IfcBridgePart deck_part;
      Schema::IfcSlab slab;
      std::map<PierIndexType, Schema::IfcBridgePart> piers; // and abutments
      std::map<pgsTypes::TrafficBarrierOrientation, Schema::IfcProduct> barriers;
      std::map<CSegmentKey, Schema::IfcBeam> segments;
   };

   // What the IFC extension does for export extensions: writes the property sets their mapping table declares for an
   // element role, and gives them what they need to create elements. Implemented by the IFC extension
   class IServices
   {
   public:
      // The ElementKind of an element role, by name (e.g. "girder", "brdr.diaphragm"). Throws if there isn't one
      virtual ElementKind Role(const std::string& name) const = 0;

      // True if the mapping table declares property sets (or classifications) for the element role
      virtual bool HasPropertySets(ElementKind role) const = 0;

      // Writes the property sets and quantity sets the mapping table declares for the element role, on the objects
      virtual void WritePropertySets(ElementKind role, const std::vector<Schema::IfcObjectDefinition>& objects, const ExportContext& context) = 0;

      // The property sets the mapping table declares for the element role that attach to type objects
      virtual std::vector<Schema::IfcPropertySetDefinition> CreateTypePropertySets(ElementKind role, const ExportContext& context) = 0;

      // Classifies an object with the classification references the mapping table declares for the element role
      virtual void Classify(ElementKind role, Schema::IfcObjectDefinition object) = 0;

      // The representation subcontext for body geometry ("Body", "Model")
      virtual Schema::IfcGeometricRepresentationSubContext BodySubContext() = 0;

      // Adds a line to the export log the caller of BuildModel receives
      virtual void Log(const std::string& line) = 0;

      virtual ~IServices() = default;
   };

   // What hooks receive
   struct HookContext
   {
      File& file;
      std::shared_ptr<WBFL::EAF::Broker> broker;
      const CIfcExportOptions& options;
      IServices& services;
      const ExportedElements& elements;
   };

   // The steps of the model build where an export extension adds to the model. Called in this order:
   //    OnExportBegin
   //    OnDeck, OnBridgeParts, OnPier (each), OnBearing (each), OnBarrier (each), OnSegment (each),
   //    OnSuperstructureComplete, OnBridgeComplete
   //    OnModelComplete
   //    OnExportEnd
   // The bridge steps are skipped when the options don't export the bridge (e.g. alignment only). A girder only export
   // calls OnSegment for its girder.
   // Property sets are only written when the options classify the model (CIfcExportOptions::classify), as in the
   // standard export, so the hooks' elements are usually created then too.
   // An exception thrown by a hook stops the export; no file is written
   class IHooks
   {
   public:
      virtual ~IHooks() = default;

      // after the project, site, and georeferencing
      virtual void OnExportBegin(HookContext& context) {}

      // after the deck slab and its property sets
      virtual void OnDeck(HookContext& context, Schema::IfcSlab slab) {}

      // after the bridge and its superstructure, substructure, and deck parts
      virtual void OnBridgeParts(HookContext& context) {}

      // after a pier or abutment, its foundation, and their property sets
      virtual void OnPier(HookContext& context, PierIndexType pierIdx, Schema::IfcBridgePart pier) {}

      // after a bearing
      virtual void OnBearing(HookContext& context, PierIndexType pierIdx, pgsTypes::PierFaceType face, const CGirderKey& girderKey, Schema::IfcBearing bearing) {}

      // after a barrier and its property sets
      virtual void OnBarrier(HookContext& context, pgsTypes::TrafficBarrierOrientation side, Schema::IfcProduct barrier) {}

      // after a girder segment, its reinforcement, and its property sets
      virtual void OnSegment(HookContext& context, const CSegmentKey& segmentKey, Schema::IfcBeam beam) {}

      // after all the girders (e.g. to add the diaphragms between them)
      virtual void OnSuperstructureComplete(HookContext& context) {}

      // after the bridge's own property sets
      virtual void OnBridgeComplete(HookContext& context) {}

      // after everything else, before the project units are completed and the model is written
      virtual void OnModelComplete(HookContext& context) {}

      // at the end of the export, whether it succeeded or not
      virtual void OnExportEnd(bool bSucceeded) noexcept {}
   };
}

// An export extension. Implemented by the plug-in that extends the export
class IIfcExportExtension
{
public:
   // e.g. "BridgeSight.BrDrExchange". Mapping tables that need the extension name it in "requires"
   virtual std::string GetId() const = 0;

   // For logs
   virtual std::string GetVersion() const = 0;

   // The prefix of the names of the extension's targets and element roles (e.g. "brdr")
   virtual std::string GetPrefix() const = 0;

   // The SDK the extension was built with. Don't override: they're computed where the extension is compiled
   virtual int GetApiVersion() const { return IFC_EXPORT_EXTENSION_API_VERSION; }
   virtual std::string GetBuildFingerprint() const { return IfcExport::BuildFingerprint(); }

   // Registers the extension's element roles and targets. Called before each export (the ElementKind of a role can
   // differ from one export to the next)
   virtual void RegisterContent(IfcExport::IRegistry& registry) = 0;

   // The hooks for one export, or nullptr
   virtual std::unique_ptr<IfcExport::IHooks> CreateHooks() = 0;

protected:
   virtual ~IIfcExportExtension() = default;
};

// What BuildModel did
struct IfcBuildResult
{
   bool bSucceeded = false;
   std::string error;                    // why it failed, for users
   std::vector<std::string> table_files; // the mapping table files, the selected table first
   std::vector<std::string> log;         // lines the extension's hooks logged
};

// {5E0A1C3B-7B6F-4C61-9D3A-2A8E4F0B7C11}
inline constexpr GUID IID_IIfcModelBuilder = { 0x5e0a1c3b, 0x7b6f, 0x4c61, { 0x9d, 0x3a, 0x2a, 0x8e, 0x4f, 0x0b, 0x7c, 0x11 } };

// The IFC extension's model builder. A broker interface of the IFC extension's agent, available in PGSuper documents
class IIfcModelBuilder
{
public:
   // The SDK the IFC extension was built with
   virtual int GetApiVersion() const = 0;
   virtual std::string GetBuildFingerprint() const = 0;

   // Builds an IFC model with the mapping table and the export extension (which may be nullptr), and writes it to the file.
   // The mapping table is options.mapping_file; when it's empty the mapping table setting or the installed standard table is used.
   // Doesn't throw: problems are reported in the result
   virtual IfcBuildResult BuildModel(const CIfcExportOptions& options, const CString& strFilePath, IIfcExportExtension* extension) = 0;

   // Writes an IDS with a specification for each property set the mapping table exports (as /IfcTableToIds), with the
   // extension's targets and element roles available to the table. Doesn't throw: problems are reported in the result
   virtual IfcBuildResult WriteTableIds(const CString& strMappingTable, const CString& strIdsPath, IIfcExportExtension* extension) = 0;

protected:
   virtual ~IIfcModelBuilder() = default;
};
