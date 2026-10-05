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
#include "IfcExtensionRegistry.h"
#include "Include\IfcExport\IfcModelBuilder.h"

#include <map>

namespace
{
   // the content of every export of an extension in this process. Never destroyed: the getters are code in the extension's
   // DLL, which may be unloaded before this one is
   std::deque<std::unique_ptr<CIfcExtensionContent>>& registered()
   {
      static auto* pContent = new std::deque<std::unique_ptr<CIfcExtensionContent>>;
      return *pContent;
   }

   // the next ElementKind for an extension's element role, unique in the process
   int g_NextRoleKind = (int)ElementKind::FirstExtension;

   const CIfcExtensionContent* g_pActive = nullptr;

   bool starts_with_prefix(const std::string& name, const std::string& prefix)
   {
      return prefix.size() + 1 < name.size() && name.compare(0, prefix.size(), prefix) == 0 && name[prefix.size()] == '.';
   }
}

// Fills in an extension's content as it registers it
class CRegistryImpl : public IfcExport::IRegistry
{
public:
   CRegistryImpl(CIfcExtensionContent& content) : m_Content(content) {}

   ElementKind RegisterRole(const IfcExport::RoleSpec& role) override
   {
      if (!starts_with_prefix(role.name, m_Content.prefix))
         Error("the element role \"" + role.name + "\" doesn't start with \"" + m_Content.prefix + ".\"");

      ElementKind existing;
      if (GetElementKind(role.name, existing) || m_Content.FindRole(role.name))
         Error("the element role \"" + role.name + "\" is registered twice");

      auto kind = (ElementKind)g_NextRoleKind++;
      m_Content.m_Roles.push_back({ role.name, role.description, kind });
      return kind;
   }

   void RegisterTarget(const IfcExport::TargetSpec& target) override
   {
      if (!starts_with_prefix(target.name, m_Content.prefix))
         Error("the target \"" + target.name + "\" doesn't start with \"" + m_Content.prefix + ".\"");

      if (FindTargetDef(target.name) || m_Content.FindTarget(target.name))
         Error("the target \"" + target.name + "\" is registered twice");

      if (!target.get)
         Error("the target \"" + target.name + "\" has no value (export extension targets are export only)");

      ElementKind element = ElementKind::Bridge;
      if (const auto* role = m_Content.FindRole(target.role))
         element = role->kind;
      else if (GetElementKind(target.role, element))
         element = GetTargetElement(element); // e.g. "abutment" targets are pier targets
      else
         Error("the target \"" + target.name + "\" belongs to the unknown element role \"" + target.role + "\"");

      auto& stored = m_Content.m_Targets.emplace_back();
      stored.name = target.name;
      stored.description = target.description;
      stored.def.name = stored.name;
      stored.def.element = element;
      stored.def.kind = target.kind;
      stored.def.description = stored.description;
      stored.def.display_unit = target.display_unit;
      stored.def.display_round = target.display_round;
      stored.def.get = target.get;
      stored.def.provider = m_Content.id;
   }

private:
   CIfcExtensionContent& m_Content;

   [[noreturn]] void Error(const std::string& problem)
   {
      throw CIfcExtensionException("The export extension \"" + m_Content.id + "\" can't be used: " + problem + ".");
   }
};

const TargetDef* CIfcExtensionContent::FindTarget(std::string_view name) const
{
   for (const auto& target : m_Targets)
   {
      if (target.def.name == name)
         return &target.def;
   }
   return nullptr;
}

const CIfcExtensionContent::Role* CIfcExtensionContent::FindRole(std::string_view name) const
{
   for (const auto& role : m_Roles)
   {
      if (role.name == name)
         return &role;
   }
   return nullptr;
}

const CIfcExtensionContent::Role* CIfcExtensionContent::FindRole(ElementKind kind) const
{
   for (const auto& role : m_Roles)
   {
      if (role.kind == kind)
         return &role;
   }
   return nullptr;
}

const CIfcExtensionContent& GetExtensionContent(IIfcExportExtension& extension)
{
   auto id = extension.GetId();

   if (extension.GetApiVersion() != IFC_EXPORT_EXTENSION_API_VERSION)
   {
      throw CIfcExtensionException("The export extension \"" + id + "\" was built for version " + std::to_string(extension.GetApiVersion()) +
         " of the IFC extension's programming interface, but the installed IFC extension has version " + std::to_string(IFC_EXPORT_EXTENSION_API_VERSION) +
         ". Install versions of the two that were released together (reinstalling BridgeLink fixes this).");
   }

   auto theirs = extension.GetBuildFingerprint();
   auto ours = IfcExport::BuildFingerprint();
   if (theirs != ours)
   {
      throw CIfcExtensionException("The export extension \"" + id + "\" wasn't built with the same tools as the installed IFC extension (" + theirs + " vs " + ours +
         "). Install versions of the two that were released together (reinstalling BridgeLink fixes this).");
   }

   auto pContent = std::make_unique<CIfcExtensionContent>();
   pContent->id = id;
   pContent->version = extension.GetVersion();
   pContent->prefix = extension.GetPrefix();
   if (pContent->prefix.empty() || pContent->prefix.find('.') != std::string::npos)
      throw CIfcExtensionException("The export extension \"" + id + "\" can't be used: its prefix must be a name without a period.");

   // the extension's targets refer to its own roles, so they're found while it registers them
   CRegistryImpl registry(*pContent);
   CIfcExtensionScope scope(pContent.get());
   extension.RegisterContent(registry);

   const auto& result = *pContent;
   registered().push_back(std::move(pContent));
   return result;
}

const CIfcExtensionContent* GetActiveExtensionContent()
{
   return g_pActive;
}

CIfcExtensionScope::CIfcExtensionScope(const CIfcExtensionContent* pContent) :
   m_pPrevious(g_pActive)
{
   g_pActive = pContent;
}

CIfcExtensionScope::~CIfcExtensionScope()
{
   g_pActive = m_pPrevious;
}
