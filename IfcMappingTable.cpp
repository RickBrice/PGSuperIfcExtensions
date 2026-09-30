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
#include "IfcMappingTable.h"
#include "bSDD.h"

#include <nlohmann/json.hpp>

#include <cstdio>
#include <cstring>

EXTERN_C IMAGE_DOS_HEADER __ImageBase; // this DLL, for finding the installed standard table

using json = nlohmann::json;

namespace
{
   constexpr const char* table_format = "PGSuperIfcMapping";
   constexpr int table_version = 1;

   std::string upper(std::string s)
   {
      std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {return (char)std::toupper(c); });
      return s;
   }

   std::string lower(std::string s)
   {
      std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {return (char)std::tolower(c); });
      return s;
   }

   std::string join(const std::vector<std::string>& items, const char* separator = ", ")
   {
      std::string result;
      for (const auto& item : items)
      {
         if (!result.empty())
            result += separator;
         result += item;
      }
      return result;
   }

   std::filesystem::path path_from_utf8(const std::string& text)
   {
      return std::filesystem::path(std::u8string(text.begin(), text.end()));
   }

   // The IFC schema declaration with the given name, or nullptr
   const ifcopenshell::declaration* find_declaration(const std::string& name)
   {
      try
      {
         return IfcSchema::get_schema().declaration_by_name(name);
      }
      catch (...)
      {
         return nullptr;
      }
   }

   std::string source_description(MappingTableSource source)
   {
      switch (source)
      {
      case MappingTableSource::CommandLine: return "the command line (/IfcMapping)";
      case MappingTableSource::ConfigurationSetting: return "the BridgeLink configuration setting";
      case MappingTableSource::InstalledStandard: return "the standard table installed with the IFC extension";
      case MappingTableSource::Extends: return "\"extends\" in another table";
      }
      ASSERT(false);
      return "";
   }

   // Where the table files came from, for error messages. chain is the tables that extend the failing one, the selected table first
   struct LoadContext
   {
      MappingTableSource root_source = MappingTableSource::InstalledStandard;
      std::vector<std::filesystem::path> chain;
   };

   [[noreturn]] void fail(const std::filesystem::path& path, MappingTableSource source, const LoadContext& context, const std::string& problem, bool bInstallationProblem)
   {
      std::ostringstream os;
      os << "IFC mapping table error." << std::endl;
      os << "File: " << PathToString(path) << std::endl;
      os << "The file was chosen by " << source_description(source);
      for (auto it = context.chain.rbegin(); it != context.chain.rend(); it++)
         os << ", extended by " << PathToString(*it);
      os << "." << std::endl;
      os << problem << std::endl;

      os << "What to do: ";
      if (bInstallationProblem && source == MappingTableSource::InstalledStandard)
      {
         os << "The IFC extension installation is incomplete or damaged. Reinstall the IFC extension.";
      }
      else if (source == MappingTableSource::Extends)
      {
         os << "Correct the table, or the \"extends\" of " << PathToString(context.chain.back()) << ".";
      }
      else
      {
         switch (context.root_source)
         {
         case MappingTableSource::CommandLine:
            os << "Correct the table, or give a different table with /IfcMapping=<file>.";
            break;
         case MappingTableSource::ConfigurationSetting:
            os << "Correct the table, or choose a different table in the BridgeLink configuration.";
            break;
         default:
            os << "Correct the table. It may have been edited after the IFC extension was installed; reinstalling the extension restores it.";
            break;
         }
      }

      throw CIfcMappingTableException(os.str());
   }

   // Reads the JSON of one table file, validates it, and fills in a MappingTableFile
   class CTableParser
   {
   public:
      CTableParser(MappingTableFile& file) : m_File(file) {}

      void Parse(const json& j)
      {
         if (!j.is_object())
         {
            Error("", "a mapping table must be a JSON object");
            return;
         }

         CheckKeys(j, "", { "format", "version", "name", "extends", "elements", "targets", "property_sets", "quantity_sets", "classification_systems", "classifications" });

         auto format = GetString(j, "format", "", true);
         if (format && *format != table_format)
            Error("format", "must be \"" + std::string(table_format) + "\". This file isn't an IFC mapping table");

         if (!j.contains("version") || !j["version"].is_number_integer())
         {
            Error("version", "missing, or not an integer");
         }
         else
         {
            m_File.version = j["version"].get<int>();
            if (m_File.version != table_version)
               Error("version", "the table has version " + std::to_string(m_File.version) + ". This version of the IFC extension reads version " + std::to_string(table_version));
         }

         if (!m_Errors.empty())
            return; // not a table we can read, don't report the rest

         m_File.name = GetString(j, "name", "", true).value_or("");
         m_File.extends = GetString(j, "extends", "", false).value_or("");

         if (j.contains("elements"))
            ParseElements(j["elements"]);

         if (j.contains("property_sets"))
            ParsePropertySets(j["property_sets"], "property_sets", false);

         if (j.contains("quantity_sets"))
            ParsePropertySets(j["quantity_sets"], "quantity_sets", true);

         if (j.contains("classification_systems"))
            ParseClassificationSystems(j["classification_systems"]);

         if (j.contains("classifications"))
            ParseClassifications(j["classifications"]);

         if (j.contains("targets"))
            ParseTargets(j["targets"]);
      }

      const std::vector<std::string>& GetErrors() const { return m_Errors; }

   private:
      MappingTableFile& m_File;
      std::vector<std::string> m_Errors;

      void Error(const std::string& path, const std::string& message)
      {
         m_Errors.push_back((path.empty() ? std::string("(top level)") : path) + ": " + message);
      }

      static std::string Path(const std::string& parent, const std::string& key)
      {
         return parent.empty() ? key : parent + "." + key;
      }

      static std::string Path(const std::string& parent, size_t index)
      {
         return parent + "[" + std::to_string(index) + "]";
      }

      // Every object may have a "comment". Any other key must be in allowed
      void CheckKeys(const json& j, const std::string& path, std::initializer_list<std::string_view> allowed)
      {
         for (const auto& [key, value] : j.items())
         {
            if (key == "comment")
               continue;

            if (std::find(allowed.begin(), allowed.end(), key) == allowed.end())
            {
               std::vector<std::string> names(allowed.begin(), allowed.end());
               Error(Path(path, key), "unknown key. Expected one of: " + join(names) + ", comment");
            }
         }
      }

      std::optional<std::string> GetString(const json& j, const char* key, const std::string& path, bool bRequired)
      {
         if (!j.contains(key))
         {
            if (bRequired)
               Error(Path(path, key), "missing");
            return std::nullopt;
         }

         const auto& value = j[key];
         if (!value.is_string() || value.get<std::string>().empty())
         {
            Error(Path(path, key), "must be a text value that isn't empty");
            return std::nullopt;
         }

         return value.get<std::string>();
      }

      std::optional<bool> GetBool(const json& j, const char* key, const std::string& path)
      {
         if (!j.contains(key))
            return std::nullopt;

         if (!j[key].is_boolean())
         {
            Error(Path(path, key), "must be true or false");
            return std::nullopt;
         }

         return j[key].get<bool>();
      }

      // an element role, or a list of element roles
      bool GetRoles(const json& j, const char* key, const std::string& path, std::vector<ElementKind>& roles)
      {
         if (!j.contains(key))
         {
            Error(Path(path, key), "missing");
            return false;
         }

         const auto& value = j[key];
         std::vector<std::string> names;
         if (value.is_string())
            names.push_back(value.get<std::string>());
         else if (value.is_array() && !value.empty() && std::all_of(value.begin(), value.end(), [](const auto& v) {return v.is_string(); }))
            names = value.get<std::vector<std::string>>();
         else
         {
            Error(Path(path, key), "must be an element role or a list of element roles");
            return false;
         }

         bool bOk = true;
         for (const auto& name : names)
         {
            ElementKind kind;
            if (GetElementKind(name, kind))
            {
               roles.push_back(kind);
            }
            else
            {
               Error(Path(path, key), "unknown element role \"" + name + "\"");
               bOk = false;
            }
         }
         return bOk;
      }

      PropertySetDeclaration::Condition GetCondition(const json& j, const std::string& path)
      {
         auto condition = GetString(j, "condition", path, false).value_or("classify");
         if (condition == "quantities")
            return PropertySetDeclaration::Condition::Quantities;
         if (condition == "always")
            return PropertySetDeclaration::Condition::Always;
         if (condition != "classify")
            Error(Path(path, "condition"), "must be \"classify\", \"quantities\", or \"always\"");
         return PropertySetDeclaration::Condition::Classify;
      }

      bool GetOwner(const json& j, const char* key, const std::string& path, PropertyOwner& owner)
      {
         auto text = GetString(j, key, path, false);
         if (!text)
            return true;

         if (*text == "occurrence")
            owner = PropertyOwner::Occurrence;
         else if (*text == "type")
            owner = PropertyOwner::Type;
         else if (*text == "material")
            owner = PropertyOwner::Material;
         else
         {
            Error(Path(path, key), "must be \"occurrence\", \"type\", or \"material\"");
            return false;
         }
         return true;
      }

      void ParseElements(const json& j)
      {
         if (!j.is_object())
         {
            Error("elements", "must be an object of element roles and selectors");
            return;
         }

         for (const auto& [role, value] : j.items())
         {
            auto path = Path("elements", role);
            if (role == "comment")
               continue;

            ElementKind kind;
            if (!GetElementKind(role, kind))
            {
               Error(path, "unknown element role");
               continue;
            }

            auto selector = ParseSelector(value, path);
            selector.table = m_File.name;
            m_File.elements[kind] = selector;
         }
      }

      ElementSelector ParseSelector(const json& j, const std::string& path)
      {
         ElementSelector selector;
         if (!j.is_object())
         {
            Error(path, "must be an object");
            return selector;
         }

         if (j.contains("any_of"))
         {
            CheckKeys(j, path, { "any_of" });
            const auto& alternatives = j["any_of"];
            if (!alternatives.is_array() || alternatives.empty())
            {
               Error(Path(path, "any_of"), "must be a list of selectors");
               return selector;
            }

            for (size_t i = 0; i < alternatives.size(); i++)
               selector.any_of.push_back(ParseSelector(alternatives[i], Path(Path(path, "any_of"), i)));

            return selector;
         }

         CheckKeys(j, path, { "entity", "predefined_type", "attributes", "classification" });

         const ifcopenshell::entity* entity = nullptr;
         if (auto name = GetString(j, "entity", path, true))
         {
            auto declaration = find_declaration(*name);
            entity = declaration ? declaration->as_entity() : nullptr;
            if (entity)
               selector.entity = *name;
            else
               Error(Path(path, "entity"), "\"" + *name + "\" isn't an IFC entity");
         }

         if (auto predefined_type = GetString(j, "predefined_type", path, false))
            selector.predefined_type = upper(*predefined_type);

         if (j.contains("attributes"))
         {
            const auto& attributes = j["attributes"];
            if (!attributes.is_object())
            {
               Error(Path(path, "attributes"), "must be an object of attribute names and values");
            }
            else
            {
               for (const auto& [name, value] : attributes.items())
               {
                  auto attribute_path = Path(Path(path, "attributes"), name);
                  if (!value.is_string())
                  {
                     Error(attribute_path, "must be a text value");
                     continue;
                  }

                  if (entity)
                  {
                     const auto& all_attributes = entity->all_attributes();
                     bool bFound = std::any_of(all_attributes.begin(), all_attributes.end(), [&name](const auto* attribute) {return attribute->name() == name; });
                     if (!bFound)
                     {
                        Error(attribute_path, "\"" + name + "\" isn't an attribute of " + selector.entity);
                        continue;
                     }
                  }

                  selector.attributes.emplace_back(name, value.get<std::string>());
               }
            }
         }

         if (auto classification = GetString(j, "classification", path, false))
            selector.classification = *classification;

         return selector;
      }

      // a bSDD URI ("bsdd" for the declaration's own name, "bsdd:Name" for another name), an explicit URI, or empty
      static std::string ResolveUri(const std::string& uri, const std::string& name, const char* kind)
      {
         if (uri.empty())
            return "";
         if (uri == "bsdd")
            return BSDD_URI + kind + "/" + name;
         if (uri.starts_with("bsdd:"))
            return BSDD_URI + kind + "/" + uri.substr(5);
         return uri;
      }

      // a constant property value. The JSON type must fit the IFC value type
      std::optional<TargetValue> GetConstant(const json& j, const std::string& path, const std::string& type)
      {
         bool bText = (type == "IFCLABEL" || type == "IFCTEXT" || type == "IFCIDENTIFIER");
         bool bBoolean = (type == "IFCBOOLEAN");
         bool bInteger = (type == "IFCINTEGER" || type == "IFCCOUNTMEASURE");
         if (bText && j.is_string())
            return j.get<std::string>();
         if (bBoolean && j.is_boolean())
            return j.get<bool>();
         if (bInteger && j.is_number_integer())
            return (Int64)j.get<int64_t>();
         if (!bText && !bBoolean && !bInteger && j.is_number())
            return (Float64)j.get<double>();

         Error(path, "doesn't fit the property's type " + type);
         return std::nullopt;
      }

      // property sets, or quantity sets (bQuantities). A declaration for a list of element roles becomes one declaration per role
      void ParsePropertySets(const json& j, const char* section, bool bQuantities)
      {
         if (!j.is_array())
         {
            Error(section, bQuantities ? "must be a list of quantity sets" : "must be a list of property sets");
            return;
         }

         for (size_t i = 0; i < j.size(); i++)
         {
            auto path = Path(section, i);
            const auto& jpset = j[i];
            if (!jpset.is_object())
            {
               Error(path, "must be an object");
               continue;
            }

            if (bQuantities)
               CheckKeys(jpset, path, { "name", "applies_to", "condition", "method", "uri", "remove", "quantities" });
            else
               CheckKeys(jpset, path, { "name", "applies_to", "attach", "condition", "uri", "shared", "remove", "properties" });

            PropertySetDeclaration pset;
            pset.name = GetString(jpset, "name", path, true).value_or("");
            std::vector<ElementKind> roles;
            bool bRoles = GetRoles(jpset, "applies_to", path, roles);
            if (!bQuantities)
               GetOwner(jpset, "attach", path, pset.attach);
            pset.condition = GetCondition(jpset, path);
            pset.uri = ResolveUri(GetString(jpset, "uri", path, false).value_or(""), pset.name, "class");
            pset.quantities = bQuantities;
            pset.method = GetString(jpset, "method", path, false).value_or("");
            pset.shared = GetBool(jpset, "shared", path).value_or(false);
            pset.remove = GetBool(jpset, "remove", path).value_or(false);

            const char* items_key = bQuantities ? "quantities" : "properties";
            if (!pset.remove)
            {
               if (!jpset.contains(items_key) || !jpset[items_key].is_array() || jpset[items_key].empty())
               {
                  Error(Path(path, items_key), bQuantities ? "missing, or not a list of quantities" : "missing, or not a list of properties");
                  continue;
               }

               const auto& jproperties = jpset[items_key];
               for (size_t k = 0; k < jproperties.size(); k++)
               {
                  auto property_path = Path(Path(path, items_key), k);
                  if (auto property = ParseProperty(jproperties[k], property_path, pset, bRoles ? roles : std::vector<ElementKind>{}))
                     pset.properties.push_back(*property);
               }
            }

            for (auto role : roles)
            {
               pset.applies_to = role;
               m_File.property_sets.push_back(pset);
            }
         }
      }

      std::optional<PropertyDeclaration> ParseProperty(const json& jproperty, const std::string& property_path, const PropertySetDeclaration& pset, const std::vector<ElementKind>& roles)
      {
         if (!jproperty.is_object())
         {
            Error(property_path, "must be an object");
            return std::nullopt;
         }

         if (pset.quantities)
            CheckKeys(jproperty, property_path, { "name", "type", "target", "value", "import" });
         else
            CheckKeys(jproperty, property_path, { "name", "type", "target", "value", "import", "uri", "enumeration" });

         PropertyDeclaration property;
         property.name = GetString(jproperty, "name", property_path, true).value_or("");

         if (auto type = GetString(jproperty, "type", property_path, true))
         {
            if (pset.quantities)
            {
               if (!IsExportQuantityType(upper(*type)))
                  Error(Path(property_path, "type"), "\"" + *type + "\" isn't a quantity type the exporter can write (IfcQuantityLength, IfcQuantityArea, IfcQuantityVolume, IfcQuantityWeight, IfcQuantityCount)");
               else
                  property.type = upper(*type);
            }
            else
            {
               auto declaration = find_declaration(*type);
               if (!declaration || declaration->as_entity())
                  Error(Path(property_path, "type"), "\"" + *type + "\" isn't an IFC value type");
               else
                  property.type = upper(*type);
            }
         }

         if (auto target_name = GetString(jproperty, "target", property_path, false))
         {
            property.target = FindTargetDef(*target_name);
            if (!property.target)
            {
               Error(Path(property_path, "target"), UnknownTarget(*target_name));
            }
            else if (pset.shared)
            {
               Error(Path(property_path, "target"), "a shared property set can't have properties with targets, because target values belong to one element");
            }
            else
            {
               for (auto role : roles)
               {
                  if (property.target->element != GetTargetElement(role))
                  {
                     Error(Path(property_path, "target"), "\"" + *target_name + "\" belongs to the " + std::string(GetElementRoleName(property.target->element)) +
                        ", but the property set applies to the " + std::string(GetElementRoleName(role)));
                     break;
                  }
               }
            }
         }

         if (jproperty.contains("value"))
         {
            if (property.target)
               Error(Path(property_path, "value"), "a property can have a target or a value, not both");
            else if (!property.type.empty())
               property.value = GetConstant(jproperty["value"], Path(property_path, "value"), pset.quantities ? std::string("IFCREAL") : property.type);
         }

         // a property without a value only declares the property, so any value type will do
         if (!pset.quantities && (property.target || property.value) && !property.type.empty() && !IsExportValueType(property.type))
            Error(Path(property_path, "type"), "\"" + property.type + "\" isn't a value type the exporter can write");

         property.import = GetBool(jproperty, "import", property_path).value_or(!pset.quantities);
         property.uri = ResolveUri(GetString(jproperty, "uri", property_path, false).value_or(""), property.name, "prop");

         if (jproperty.contains("enumeration"))
         {
            const auto& jenum = jproperty["enumeration"];
            auto enum_path = Path(property_path, "enumeration");
            if (!jenum.is_object())
            {
               Error(enum_path, "must be an object with a name and values");
            }
            else
            {
               CheckKeys(jenum, enum_path, { "name", "values" });
               property.enumeration_name = GetString(jenum, "name", enum_path, true).value_or("");
               if (!jenum.contains("values") || !jenum["values"].is_array() || jenum["values"].empty() || !std::all_of(jenum["values"].begin(), jenum["values"].end(), [](const auto& v) {return v.is_string(); }))
                  Error(Path(enum_path, "values"), "missing, or not a list of text values");
               else
                  property.enumeration_values = jenum["values"].get<std::vector<std::string>>();

               if (property.target && property.target->kind != ValueKind::Text)
                  Error(Path(property_path, "target"), "an enumerated property needs a text target");

               if (property.value && std::holds_alternative<std::string>(*property.value) &&
                  std::find(property.enumeration_values.begin(), property.enumeration_values.end(), std::get<std::string>(*property.value)) == property.enumeration_values.end())
                  Error(Path(property_path, "value"), "isn't one of the enumeration values");
            }
         }

         return property;
      }

      void ParseClassificationSystems(const json& j)
      {
         if (!j.is_array())
         {
            Error("classification_systems", "must be a list of classification systems");
            return;
         }

         for (size_t i = 0; i < j.size(); i++)
         {
            auto path = Path("classification_systems", i);
            const auto& jsystem = j[i];
            if (!jsystem.is_object())
            {
               Error(path, "must be an object");
               continue;
            }

            CheckKeys(jsystem, path, { "name", "source", "edition", "edition_date", "specification" });

            ClassificationSystemDeclaration system;
            system.name = GetString(jsystem, "name", path, true).value_or("");
            system.source = GetString(jsystem, "source", path, false).value_or("");
            system.edition = GetString(jsystem, "edition", path, false).value_or("");
            system.edition_date = GetString(jsystem, "edition_date", path, false).value_or("");
            auto specification = GetString(jsystem, "specification", path, false).value_or("");
            system.specification = (specification == "bsdd" ? BSDD_URI : specification);
            m_File.classification_systems.push_back(system);
         }
      }

      void ParseClassifications(const json& j)
      {
         if (!j.is_array())
         {
            Error("classifications", "must be a list of classifications");
            return;
         }

         for (size_t i = 0; i < j.size(); i++)
         {
            auto path = Path("classifications", i);
            const auto& jclassification = j[i];
            if (!jclassification.is_object())
            {
               Error(path, "must be an object");
               continue;
            }

            CheckKeys(jclassification, path, { "applies_to", "condition", "system", "identification", "name", "location", "remove" });

            ClassificationDeclaration classification;
            std::vector<ElementKind> roles;
            GetRoles(jclassification, "applies_to", path, roles);
            classification.condition = GetCondition(jclassification, path);
            classification.identification = GetString(jclassification, "identification", path, true).value_or("");
            classification.remove = GetBool(jclassification, "remove", path).value_or(false);
            if (!classification.remove)
            {
               classification.system = GetString(jclassification, "system", path, true).value_or("");
               classification.name = GetString(jclassification, "name", path, false).value_or("");
               auto location = GetString(jclassification, "location", path, false).value_or("");
               classification.location = (location == "bsdd" ? BSDD_URI + "class/" + classification.identification : location);
            }

            for (auto role : roles)
            {
               classification.applies_to = role;
               m_File.classifications.push_back(classification);
            }
         }
      }

      static std::string UnknownTarget(const std::string& name)
      {
         std::vector<std::string> names;
         for (const auto& target : GetTargetDefs())
            names.emplace_back(target.name);
         return "unknown target \"" + name + "\". Targets: " + join(names);
      }

      void ParseTargets(const json& j)
      {
         if (!j.is_object())
         {
            Error("targets", "must be an object of targets and their locations");
            return;
         }

         for (const auto& [name, value] : j.items())
         {
            if (name == "comment")
               continue;

            auto path = Path("targets", name);
            const TargetDef* target = FindTargetDef(name);
            if (!target)
            {
               Error(path, UnknownTarget(name));
               continue;
            }

            const json* jlocations = &value;
            auto locations_path = path;
            if (value.is_object())
            {
               // { "mode": "replace", "locations": [...] }
               CheckKeys(value, path, { "mode", "locations" });
               auto mode = GetString(value, "mode", path, false).value_or("prepend");
               if (mode == "replace")
                  m_File.replace_targets.insert(name);
               else if (mode != "prepend")
                  Error(Path(path, "mode"), "must be \"prepend\" or \"replace\"");

               if (!value.contains("locations"))
               {
                  Error(Path(path, "locations"), "missing");
                  continue;
               }
               jlocations = &value["locations"];
               locations_path = Path(path, "locations");
            }

            if (!jlocations->is_array())
            {
               Error(locations_path, "must be a list of locations");
               continue;
            }

            auto& locations = m_File.targets[name];
            for (size_t i = 0; i < jlocations->size(); i++)
            {
               auto location = ParseLocation((*jlocations)[i], Path(locations_path, i), *target);
               if (location)
                  locations.push_back(std::move(*location));
            }
         }
      }

      std::optional<MappingLocation> ParseLocation(const json& j, const std::string& path, const TargetDef& target)
      {
         if (!j.is_object())
         {
            Error(path, "must be an object");
            return std::nullopt;
         }

         CheckKeys(j, path, { "property", "attribute", "type_attribute", "classification", "field", "on", "value_types", "unit", "parse", "list_index", "map" });

         size_t errors = m_Errors.size();

         MappingLocation location;
         location.table = m_File.name;
         location.origin = PathToString(m_File.path) + " " + path;

         int nKinds = (int)j.contains("property") + (int)j.contains("attribute") + (int)j.contains("type_attribute") + (int)j.contains("classification");
         if (nKinds != 1)
         {
            Error(path, "must have exactly one of \"property\", \"attribute\", \"type_attribute\", or \"classification\"");
            return std::nullopt;
         }

         if (j.contains("property"))
         {
            location.kind = MappingLocation::Kind::Property;
            const auto& jproperty = j["property"];
            auto property_path = Path(path, "property");
            if (!jproperty.is_object())
            {
               Error(property_path, "must be an object with \"pset\" and \"name\"");
            }
            else
            {
               CheckKeys(jproperty, property_path, { "pset", "name" });
               location.pset = GetString(jproperty, "pset", property_path, true).value_or("");
               location.name = GetString(jproperty, "name", property_path, true).value_or("");
            }
            GetOwner(j, "on", path, location.on);
         }
         else if (j.contains("attribute") || j.contains("type_attribute"))
         {
            bool bType = j.contains("type_attribute");
            location.kind = bType ? MappingLocation::Kind::TypeAttribute : MappingLocation::Kind::Attribute;
            location.name = GetString(j, bType ? "type_attribute" : "attribute", path, true).value_or("");
         }
         else
         {
            location.kind = MappingLocation::Kind::Classification;
            const auto& jclassification = j["classification"];
            auto classification_path = Path(path, "classification");
            if (!jclassification.is_object())
            {
               Error(classification_path, "must be an object (it may be empty)");
            }
            else
            {
               CheckKeys(jclassification, classification_path, { "system", "identification" });
               location.classification_system = GetString(jclassification, "system", classification_path, false).value_or("");
               location.classification_identification = GetString(jclassification, "identification", classification_path, false).value_or("");
            }

            auto field = GetString(j, "field", path, false).value_or("Name");
            if (field == "Name")
               location.classification_field_is_name = true;
            else if (field == "Identification")
               location.classification_field_is_name = false;
            else
               Error(Path(path, "field"), "must be \"Name\" or \"Identification\"");
         }

         if (j.contains("on") && location.kind != MappingLocation::Kind::Property)
            Error(Path(path, "on"), "only a property location can have \"on\"");

         if (j.contains("field") && location.kind != MappingLocation::Kind::Classification)
            Error(Path(path, "field"), "only a classification location can have \"field\"");

         if (j.contains("value_types"))
         {
            const auto& jtypes = j["value_types"];
            if (!jtypes.is_array())
            {
               Error(Path(path, "value_types"), "must be a list of IFC value types");
            }
            else
            {
               for (size_t i = 0; i < jtypes.size(); i++)
               {
                  auto declaration = jtypes[i].is_string() ? find_declaration(jtypes[i].get<std::string>()) : nullptr;
                  if (!declaration || declaration->as_entity())
                     Error(Path(Path(path, "value_types"), i), "isn't an IFC value type");
                  else
                     location.value_types.push_back(upper(jtypes[i].get<std::string>()));
               }
            }
         }

         if (auto unit_name = GetString(j, "unit", path, false))
         {
            location.unit = FindTableUnit(*unit_name);
            if (!HasUnit(target.kind))
            {
               Error(Path(path, "unit"), "\"" + std::string(target.name) + "\" is a " + std::string(GetValueKindName(target.kind)) + " and has no unit");
               location.unit = nullptr;
            }
            else if (!location.unit)
            {
               std::vector<std::string> names;
               for (const auto* name : { "Pa", "kPa", "MPa", "psi", "ksi", "mm", "cm", "m", "in", "ft", "rad", "deg", "N", "kN", "lbf", "kip" })
               {
                  if (FindTableUnit(name)->kind == target.kind)
                     names.emplace_back(name);
               }
               Error(Path(path, "unit"), "unknown unit \"" + *unit_name + "\". Units for a " + std::string(GetValueKindName(target.kind)) + ": " + join(names));
            }
            else if (location.unit->kind != target.kind)
            {
               Error(Path(path, "unit"), "\"" + *unit_name + "\" is a " + std::string(GetValueKindName(location.unit->kind)) + " unit, but \"" + std::string(target.name) + "\" is a " + std::string(GetValueKindName(target.kind)));
               location.unit = nullptr;
            }
         }

         if (j.contains("parse"))
         {
            const auto& jparse = j["parse"];
            auto parse_path = Path(path, "parse");
            if (jparse.is_string() && jparse.get<std::string>() == "number")
            {
               location.parse = MappingLocation::Parse::Number;
               if (!IsNumeric(target.kind))
                  Error(parse_path, "\"number\" is only for targets that are numbers");
            }
            else if (jparse.is_string() && jparse.get<std::string>() == "feet_inches")
            {
               location.parse = MappingLocation::Parse::FeetInches;
               if (target.kind != ValueKind::Length)
                  Error(parse_path, "\"feet_inches\" is only for length targets");
            }
            else if (jparse.is_object())
            {
               CheckKeys(jparse, parse_path, { "regex", "group" });
               location.parse = MappingLocation::Parse::Regex;
               if (auto pattern = GetString(jparse, "regex", parse_path, true))
               {
                  try
                  {
                     location.regex = std::regex(*pattern, std::regex::ECMAScript | std::regex::icase);
                  }
                  catch (const std::regex_error& e)
                  {
                     Error(Path(parse_path, "regex"), std::string("isn't a valid regular expression: ") + e.what());
                  }
               }

               if (jparse.contains("group"))
               {
                  if (!jparse["group"].is_number_unsigned())
                     Error(Path(parse_path, "group"), "must be a number, 0 or more");
                  else
                     location.regex_group = jparse["group"].get<size_t>();
               }
            }
            else
            {
               Error(parse_path, "must be \"number\", \"feet_inches\", or { \"regex\": \"...\", \"group\": 1 }");
            }
         }

         if (j.contains("list_index"))
         {
            if (!j["list_index"].is_number_unsigned())
               Error(Path(path, "list_index"), "must be a number, 0 or more");
            else
               location.list_index = j["list_index"].get<size_t>();
         }

         if (j.contains("map"))
         {
            const auto& jmap = j["map"];
            auto map_path = Path(path, "map");
            if (!jmap.is_object() || jmap.empty())
            {
               Error(map_path, "must be an object of model values and target values");
            }
            else if (target.kind != ValueKind::Boolean && target.kind != ValueKind::Text && target.kind != ValueKind::TextList)
            {
               Error(map_path, "only boolean and text targets can have a map");
            }
            else
            {
               for (const auto& [key, value] : jmap.items())
               {
                  if (target.kind == ValueKind::Boolean && value.is_boolean())
                     location.map.emplace_back(lower(key), value.get<bool>());
                  else if (target.kind != ValueKind::Boolean && value.is_string())
                     location.map.emplace_back(lower(key), value.get<std::string>());
                  else
                     Error(Path(map_path, key), target.kind == ValueKind::Boolean ? "must be true or false" : "must be a text value");
               }
            }
         }

         if (errors != m_Errors.size())
            return std::nullopt;

         return location;
      }
   };

   MappingTableFile load_file(const std::filesystem::path& path, MappingTableSource source, const LoadContext& context)
   {
      std::error_code ec;
      if (!std::filesystem::exists(path, ec))
         fail(path, source, context, "The file was not found.", true);

      FILE* fp = nullptr;
      errno_t err = _wfopen_s(&fp, path.c_str(), L"rb");
      if (err != 0 || fp == nullptr)
      {
         char buffer[256];
         strerror_s(buffer, err);
         fail(path, source, context, std::string("The file can't be read: ") + buffer, true);
      }

      std::string content;
      char buffer[4096];
      size_t count;
      while ((count = fread(buffer, 1, sizeof(buffer), fp)) != 0)
         content.append(buffer, count);
      fclose(fp);

      json j;
      try
      {
         j = json::parse(content);
      }
      catch (const json::parse_error& e)
      {
         fail(path, source, context, std::string("The file isn't valid JSON: ") + e.what(), true);
      }

      MappingTableFile file;
      file.path = path;
      file.source = source;

      CTableParser parser(file);
      parser.Parse(j);
      const auto& errors = parser.GetErrors();
      if (!errors.empty())
      {
         std::ostringstream os;
         os << "The table has " << errors.size() << (errors.size() == 1 ? " problem:" : " problems:");
         for (const auto& error : errors)
            os << std::endl << "   " << error;
         fail(path, source, context, os.str(), false);
      }

      return file;
   }

   bool same_file(const std::filesystem::path& a, const std::filesystem::path& b)
   {
      std::error_code ec;
      return std::filesystem::equivalent(a, b, ec);
   }
}

const TableUnit* FindTableUnit(std::string_view name)
{
   using namespace WBFL::Units;
   static const TableUnit units[] = {
      { "Pa", ValueKind::Stress, [](Float64 v) {return ConvertToSysUnits(v, Measure::Pa); } },
      { "kPa", ValueKind::Stress, [](Float64 v) {return ConvertToSysUnits(v, Measure::kPa); } },
      { "MPa", ValueKind::Stress, [](Float64 v) {return ConvertToSysUnits(v, Measure::MPa); } },
      { "psi", ValueKind::Stress, [](Float64 v) {return ConvertToSysUnits(v, Measure::PSI); } },
      { "ksi", ValueKind::Stress, [](Float64 v) {return ConvertToSysUnits(v, Measure::KSI); } },
      { "mm", ValueKind::Length, [](Float64 v) {return ConvertToSysUnits(v, Measure::Millimeter); } },
      { "cm", ValueKind::Length, [](Float64 v) {return ConvertToSysUnits(v, Measure::Centimeter); } },
      { "m", ValueKind::Length, [](Float64 v) {return ConvertToSysUnits(v, Measure::Meter); } },
      { "in", ValueKind::Length, [](Float64 v) {return ConvertToSysUnits(v, Measure::Inch); } },
      { "ft", ValueKind::Length, [](Float64 v) {return ConvertToSysUnits(v, Measure::Feet); } },
      { "rad", ValueKind::Angle, [](Float64 v) {return ConvertToSysUnits(v, Measure::Radian); } },
      { "deg", ValueKind::Angle, [](Float64 v) {return ConvertToSysUnits(v, Measure::Degree); } },
      { "N", ValueKind::Force, [](Float64 v) {return ConvertToSysUnits(v, Measure::Newton); } },
      { "kN", ValueKind::Force, [](Float64 v) {return ConvertToSysUnits(v, Measure::Kilonewton); } },
      { "lbf", ValueKind::Force, [](Float64 v) {return ConvertToSysUnits(v, Measure::Pound); } },
      { "kip", ValueKind::Force, [](Float64 v) {return ConvertToSysUnits(v, Measure::Kip); } },
   };

   auto found = std::find_if(std::begin(units), std::end(units), [name](const auto& unit) {return unit.name == name; });
   return found == std::end(units) ? nullptr : found;
}

bool IsExportQuantityType(const std::string& type)
{
   return type == "IFCQUANTITYLENGTH" || type == "IFCQUANTITYAREA" || type == "IFCQUANTITYVOLUME" || type == "IFCQUANTITYWEIGHT" || type == "IFCQUANTITYCOUNT";
}

bool IsExportValueType(const std::string& type)
{
   // the value types CIfcPropertyWriter can create (see CreateIfcValue in IfcPropertyWriter.h)
   static const std::set<std::string> types{
      "IFCLABEL", "IFCTEXT", "IFCIDENTIFIER", "IFCBOOLEAN", "IFCINTEGER", "IFCREAL", "IFCCOUNTMEASURE",
      "IFCPRESSUREMEASURE", "IFCLENGTHMEASURE", "IFCPOSITIVELENGTHMEASURE", "IFCNONNEGATIVELENGTHMEASURE",
      "IFCPLANEANGLEMEASURE", "IFCPOSITIVEPLANEANGLEMEASURE", "IFCRATIOMEASURE", "IFCPOSITIVERATIOMEASURE",
      "IFCAREAMEASURE", "IFCVOLUMEMEASURE", "IFCMASSMEASURE", "IFCFORCEMEASURE",
   };
   return types.find(type) != types.end();
}

std::string MappingLocation::Describe() const
{
   std::ostringstream os;
   switch (kind)
   {
   case Kind::Property:
      os << pset << "." << name;
      if (on == PropertyOwner::Type)
         os << " (type)";
      else if (on == PropertyOwner::Material)
         os << " (material)";
      break;
   case Kind::Attribute:
      os << "attribute " << name;
      break;
   case Kind::TypeAttribute:
      os << "type attribute " << name;
      break;
   case Kind::Classification:
      os << "classification";
      if (!classification_system.empty())
         os << " " << classification_system;
      if (!classification_identification.empty())
         os << " " << classification_identification;
      os << (classification_field_is_name ? " Name" : " Identification");
      break;
   }

   if (list_index)
      os << "[" << *list_index << "]";

   if (unit)
      os << " [" << unit->name << "]"; // brackets, because property names often have the unit in parentheses

   return os.str();
}

std::string ElementSelector::Describe() const
{
   if (!any_of.empty())
   {
      std::vector<std::string> alternatives;
      for (const auto& alternative : any_of)
         alternatives.push_back(alternative.Describe());
      return join(alternatives, " or ");
   }

   std::ostringstream os;
   os << entity;
   if (!predefined_type.empty())
      os << "." << predefined_type;
   for (const auto& [name, value] : attributes)
      os << " with " << name << " = \"" << value << "\"";
   if (!classification.empty())
      os << " classified " << classification;
   return os.str();
}

std::string PathToString(const std::filesystem::path& path)
{
   try
   {
      return path.string();
   }
   catch (...)
   {
      auto text = path.u8string();
      return std::string(text.begin(), text.end());
   }
}

std::filesystem::path CIfcMappingTable::GetStandardTablePath()
{
   wchar_t module_path[MAX_PATH];
   DWORD length = ::GetModuleFileNameW((HMODULE)&__ImageBase, module_path, MAX_PATH);
   std::filesystem::path path(std::wstring(module_path, length));
   return path.parent_path() / L"MappingTables" / L"Standard.json";
}

std::unique_ptr<CIfcMappingTable> CIfcMappingTable::Load(const std::filesystem::path& selected_path, MappingTableSource source)
{
   auto table = std::make_unique<CIfcMappingTable>();

   LoadContext context;
   std::filesystem::path path = selected_path;
   if (path.empty())
   {
      path = GetStandardTablePath();
      source = MappingTableSource::InstalledStandard;
   }
   context.root_source = source;

   while (true)
   {
      auto file = load_file(path, source, context);
      auto extends = file.extends;
      table->m_Files.push_back(std::move(file));

      if (extends.empty())
         break;

      std::filesystem::path next = (extends == "standard") ? GetStandardTablePath() : path.parent_path() / path_from_utf8(extends);
      context.chain.push_back(path);
      source = (extends == "standard") ? MappingTableSource::InstalledStandard : MappingTableSource::Extends;

      for (const auto& previous : context.chain)
      {
         if (same_file(previous, next))
         {
            std::vector<std::string> files;
            for (const auto& p : context.chain)
               files.push_back(PathToString(p));
            files.push_back(PathToString(next));
            fail(next, source, context, "The tables extend each other in a cycle: " + join(files, " -> "), false);
         }
      }

      path = next;
   }

   table->Merge();
   return table;
}

void CIfcMappingTable::Merge()
{
   for (const auto& target : GetTargetDefs())
   {
      std::vector<MappingLocation> locations;
      for (const auto& file : m_Files)
      {
         auto found = file.targets.find(target.name);
         if (found != file.targets.end())
         {
            locations.insert(locations.end(), found->second.begin(), found->second.end());
            if (file.replace_targets.find(target.name) != file.replace_targets.end())
               break; // this table replaces the locations of the tables it extends
         }
         else
         {
            // no explicit locations, so the properties bound to the target are the locations, in declaration order
            for (size_t i = 0; i < file.property_sets.size(); i++)
            {
               const auto& pset = file.property_sets[i];
               if (pset.quantities)
                  continue; // the reader doesn't read quantities yet
               for (size_t k = 0; k < pset.properties.size(); k++)
               {
                  const auto& property = pset.properties[k];
                  if (property.target != &target || !property.import)
                     continue;

                  MappingLocation location;
                  location.kind = MappingLocation::Kind::Property;
                  location.pset = pset.name;
                  location.name = property.name;
                  location.on = pset.attach;
                  location.table = file.name;
                  location.origin = PathToString(file.path) + " property_sets[" + std::to_string(i) + "].properties[" + std::to_string(k) + "]";
                  locations.push_back(location);
               }
            }
         }
      }
      m_Locations.emplace(std::string(target.name), std::move(locations));
   }

   // the selected table's selectors win
   for (const auto& file : m_Files)
   {
      for (const auto& [role, selector] : file.elements)
         m_Selectors.emplace(role, selector);
   }

   // export property sets: from the base table up, an extending table's property set replaces or removes the one
   // with the same name, element role, and attach, in place. Others are added at the end
   auto same = [](const PropertySetDeclaration* a, const PropertySetDeclaration& b) {return a->name == b.name && a->applies_to == b.applies_to && a->attach == b.attach; };
   for (auto file = m_Files.rbegin(); file != m_Files.rend(); file++)
   {
      for (const auto& pset : file->property_sets)
      {
         auto found = std::find_if(m_PropertySets.begin(), m_PropertySets.end(), [&](const auto* p) {return same(p, pset); });
         if (pset.remove)
         {
            if (found != m_PropertySets.end())
               m_PropertySets.erase(found);
         }
         else if (found != m_PropertySets.end())
         {
            *found = &pset;
         }
         else
         {
            m_PropertySets.push_back(&pset);
         }
      }

      for (const auto& system : file->classification_systems)
      {
         auto found = std::find_if(m_ClassificationSystems.begin(), m_ClassificationSystems.end(), [&](const auto* s) {return s->name == system.name; });
         if (found != m_ClassificationSystems.end())
            *found = &system;
         else
            m_ClassificationSystems.push_back(&system);
      }

      for (const auto& classification : file->classifications)
      {
         auto found = std::find_if(m_Classifications.begin(), m_Classifications.end(), [&](const auto* c) {return c->applies_to == classification.applies_to && c->identification == classification.identification; });
         if (classification.remove)
         {
            if (found != m_Classifications.end())
               m_Classifications.erase(found);
         }
         else if (found != m_Classifications.end())
         {
            *found = &classification;
         }
         else
         {
            m_Classifications.push_back(&classification);
         }
      }
   }
}

std::vector<const ClassificationDeclaration*> CIfcMappingTable::GetClassifications(ElementKind role) const
{
   std::vector<const ClassificationDeclaration*> classifications;
   std::copy_if(m_Classifications.begin(), m_Classifications.end(), std::back_inserter(classifications), [role](const auto* c) {return c->applies_to == role; });
   return classifications;
}

std::vector<const PropertySetDeclaration*> CIfcMappingTable::GetPropertySets(ElementKind role, PropertyOwner attach) const
{
   std::vector<const PropertySetDeclaration*> psets;
   std::copy_if(m_PropertySets.begin(), m_PropertySets.end(), std::back_inserter(psets), [role, attach](const auto* pset) {return pset->applies_to == role && pset->attach == attach; });
   return psets;
}

const std::vector<MappingLocation>& CIfcMappingTable::GetLocations(const TargetDef& target) const
{
   static const std::vector<MappingLocation> none;
   auto found = m_Locations.find(target.name);
   return found == m_Locations.end() ? none : found->second;
}

const ElementSelector* CIfcMappingTable::GetSelector(ElementKind role) const
{
   auto found = m_Selectors.find(role);
   return found == m_Selectors.end() ? nullptr : &(found->second);
}

void CIfcMappingTable::LogFiles() const
{
   for (const auto& file : m_Files)
   {
      std::ostringstream os;
      os << "IFC mapping table \"" << file.name << "\" (version " << file.version << "): " << PathToString(file.path) << ", chosen by " << source_description(file.source);
      WBFL::System::Logger::Info(os.str());
   }
}
