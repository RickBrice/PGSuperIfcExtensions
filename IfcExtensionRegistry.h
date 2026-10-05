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
   The element roles and targets of export extensions (devdocs/ExtensionsDesign.md)

   An extension's content is registered for each export that uses it, and kept
   for the life of the process (mapping tables hold pointers to its targets, and
   the extension's DLL may be unloaded and loaded again between exports).
   It's only visible while a CIfcExtensionScope for the extension exists (during
   IIfcModelBuilder calls that name the extension): FindTargetDef,
   GetElementKind, and GetElementRoleName find it then, and never during the
   standard export, import, or in the mapping table editor.
*****************************************************************************/

#include "IfcTargets.h"

#include <deque>
#include <stdexcept>

class IIfcExportExtension;

// An export extension can't be used, or one of its hooks failed. The message is for users
class CIfcExtensionException : public std::runtime_error
{
public:
   using std::runtime_error::runtime_error;
};

// The element roles and targets an export extension registered
class CIfcExtensionContent
{
public:
   struct Role
   {
      std::string name;
      std::string description;
      ElementKind kind;
   };

   std::string id;
   std::string version;
   std::string prefix;

   const TargetDef* FindTarget(std::string_view name) const;
   const Role* FindRole(std::string_view name) const;
   const Role* FindRole(ElementKind kind) const;

   const std::deque<Role>& GetRoles() const { return m_Roles; }

private:
   friend class CRegistryImpl;

   struct Target
   {
      std::string name;
      std::string description;
      TargetDef def; // its name and description refer to the strings above
   };

   // deques, so the elements don't move as more are added
   std::deque<Role> m_Roles;
   std::deque<Target> m_Targets;
};

// The content of an export extension, registered for an export. Checks that the extension was built with this
// version of the SDK. Throws CIfcExtensionException if the extension can't be used
const CIfcExtensionContent& GetExtensionContent(IIfcExportExtension& extension);

// The content of the extension whose export is in progress, or nullptr
const CIfcExtensionContent* GetActiveExtensionContent();

// Makes an extension's content visible while it exists. With nullptr, nothing is made visible
class CIfcExtensionScope
{
public:
   explicit CIfcExtensionScope(const CIfcExtensionContent* pContent);
   ~CIfcExtensionScope();

   CIfcExtensionScope(const CIfcExtensionScope&) = delete;
   CIfcExtensionScope& operator=(const CIfcExtensionScope&) = delete;

private:
   const CIfcExtensionContent* m_pPrevious;
};
