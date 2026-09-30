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
#include "IfcTargetHints.h"
#include "Utilities.h"

namespace
{
   std::string lower(std::string s)
   {
      std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {return (char)std::tolower(c); });
      return s;
   }

   // normalized name: upper case letters and digits only (e.g. "2_Type" -> "2TYPE")
   std::string normalize(const std::string& name)
   {
      std::string n;
      for (unsigned char c : name)
      {
         if (std::isalnum(c))
            n += (char)std::toupper(c);
      }
      return n;
   }

   // true for values that say nothing (empty, "<null>", "BEAM", ...), which aren't worth a hint
   bool is_generic(const std::string& value)
   {
      static const std::set<std::string> generic{ "", "NULL", "NONE", "BEAM", "GIRDER", "NOTDEFINED", "USERDEFINED" };
      return generic.find(normalize(value)) != generic.end();
   }

   std::optional<std::string> text_value(IfcSchema::IfcPropertySingleValue value)
   {
      if (!value || !value.NominalValue())
         return std::nullopt;

      if (auto label = value.NominalValue().as<IfcSchema::IfcLabel>())
         return std::string(label);
      if (auto text = value.NominalValue().as<IfcSchema::IfcText>())
         return std::string(text);
      return std::nullopt;
   }

   // "pset.property = "value"" for each text property of the element's property sets whose name satisfies the predicate
   template <typename Predicate>
   void find_text_properties(IfcSchema::IfcObject object, Predicate is_candidate, std::vector<std::string>& hints)
   {
      for (auto& rel : object.IsDefinedBy())
      {
         auto pset = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySet>();
         if (!pset)
            continue;

         for (auto& property : pset.HasProperties())
         {
            std::string name = property.Name();
            if (!is_candidate(name))
               continue;

            if (auto text = text_value(property.as<IfcSchema::IfcPropertySingleValue>()); text && !is_generic(*text))
               hints.push_back(pset.Name().value_or("") + "." + name + " = \"" + *text + "\"");
         }
      }
   }

   bool contains_haunch(const std::optional<std::string>& text)
   {
      return text && lower(*text).find("haunch") != std::string::npos;
   }
}

std::vector<std::string> CIfcTargetHints::Find(const TargetDef& target, IfcSchema::IfcObject object)
{
   std::vector<std::string> hints;
   if (!object)
      return hints;

   if (target.name == "girder.type_names")
   {
      // the rules the importer used to guess the girder type, before mapping tables
      if (object.ObjectType() && !is_generic(*object.ObjectType()))
         hints.push_back("attribute ObjectType = \"" + *object.ObjectType() + "\"");

      find_text_properties(object, [](const std::string& name)
         {
            auto n = normalize(name);
            return n.ends_with("TYPE") || n.find("SHAPE") != std::string::npos;
         }, hints);

      for (auto& rel : object.HasAssociations())
      {
         auto rel_classification = rel.as<IfcSchema::IfcRelAssociatesClassification>();
         auto reference = rel_classification ? rel_classification.RelatingClassification().as<IfcSchema::IfcClassificationReference>() : IfcSchema::IfcClassificationReference{};
         // usBridge classifications identify the kind of element, not the girder type
         if (reference && reference.Name() && !is_generic(*reference.Name()) && !reference.Identification().value_or("").starts_with("usBridge_"))
            hints.push_back("classification " + reference.Identification().value_or("") + " Name = \"" + *reference.Name() + "\"");
      }
   }
   else if (target.name == "bearing.fixed_x" || target.name == "bearing.fixed_y")
   {
      find_text_properties(object, [](const std::string& name) {return normalize(name).ends_with("FIXITY"); }, hints);
   }

   return hints;
}

void CIfcTargetHints::LogHaunchHints(ifcopenshell::file& file, const std::set<int>& excluded_ids)
{
   std::vector<std::string> candidates;
   auto check = [&](auto element)
      {
         if (excluded_ids.find(element.id()) != excluded_ids.end())
            return;
         if (contains_haunch(element.ObjectType()) || contains_haunch(element.Name()))
         {
            std::ostringstream os;
            os << element.declaration().name() << " #" << element.id() << " (Name \"" << element.Name().value_or("") << "\", ObjectType \"" << element.ObjectType().value_or("") << "\")";
            candidates.push_back(os.str());
         }
      };

   for (auto& slab : file.instances_by_type<IfcSchema::IfcSlab>())
      check(slab);
   for (auto& part : file.instances_by_type<IfcSchema::IfcBuildingElementPart>())
      check(part);

   if (candidates.empty())
      return;

   std::ostringstream os;
   os << "Hint: " << candidates.size() << " elements may be haunches, e.g. " << candidates.front()
      << ". They aren't used. To measure them with the deck, add a \"haunch\" element role to the mapping table.";
   WBFL::System::Logger::Info(os.str());
}
