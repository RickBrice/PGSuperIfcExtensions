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
#include "IfcTableIds.h"
#include "IfcTableFormat.h"
#include "IdsBuilder.h"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>

using json = nlohmann::ordered_json;
using namespace ids_builder;

namespace
{
   // The writer records what a general IDS can't say in the instructions of a facet, so a table generated
   // back from the IDS is the same: "PGSuper: target=girder.fc; attach=type; ..."
   const std::string INSTRUCTIONS_TAG = "PGSuper:";

   std::string upper(std::string s)
   {
      std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {return (char)std::toupper(c); });
      return s;
   }

   bool same_text(const std::string& a, const std::string& b)
   {
      return upper(a) == upper(b);
   }

   std::string instructions(const std::vector<std::pair<std::string, std::string>>& items)
   {
      std::string text;
      for (const auto& [key, value] : items)
      {
         if (value.empty())
            continue;
         text += (text.empty() ? INSTRUCTIONS_TAG + " " : std::string("; ")) + key + "=" + value;
      }
      return text;
   }

   // the key=value items of PGSuper instructions, or nothing if the instructions aren't PGSuper's
   std::map<std::string, std::string> parse_instructions(const std::optional<std::string>& text)
   {
      std::map<std::string, std::string> items;
      if (!text || !text->starts_with(INSTRUCTIONS_TAG))
         return items;

      std::stringstream ss(text->substr(INSTRUCTIONS_TAG.size()));
      std::string item;
      while (std::getline(ss, item, ';'))
      {
         auto eq = item.find('=');
         if (eq == std::string::npos)
            continue;
         auto trim = [](std::string s) { s.erase(0, s.find_first_not_of(' ')); s.erase(s.find_last_not_of(' ') + 1); return s; };
         items[trim(item.substr(0, eq))] = trim(item.substr(eq + 1));
      }
      return items;
   }

   // quantities are checked as the measure they hold, and generated back from it
   const std::vector<std::pair<std::string, std::string>>& quantity_measures()
   {
      static const std::vector<std::pair<std::string, std::string>> map{
         { "IFCQUANTITYLENGTH", "IFCLENGTHMEASURE" },
         { "IFCQUANTITYAREA", "IFCAREAMEASURE" },
         { "IFCQUANTITYVOLUME", "IFCVOLUMEMEASURE" },
         { "IFCQUANTITYWEIGHT", "IFCMASSMEASURE" },
         { "IFCQUANTITYCOUNT", "IFCCOUNTMEASURE" },
      };
      return map;
   }

   std::string measure_of_quantity(const std::string& quantity_type)
   {
      for (const auto& [q, m] : quantity_measures())
         if (q == quantity_type) return m;
      return quantity_type;
   }

   std::string quantity_of_measure(const std::string& measure)
   {
      for (const auto& [q, m] : quantity_measures())
         if (m == measure) return q;
      return "";
   }

   std::string value_text(const TargetValue& value)
   {
      if (std::holds_alternative<std::string>(value)) return std::get<std::string>(value);
      if (std::holds_alternative<bool>(value)) return std::get<bool>(value) ? "true" : "false";
      if (std::holds_alternative<Int64>(value)) return std::to_string(std::get<Int64>(value));
      std::ostringstream os;
      os << std::defaultfloat << std::setprecision(10) << std::get<Float64>(value);
      return os.str();
   }

   std::string condition_name(PropertySetDeclaration::Condition condition)
   {
      switch (condition)
      {
      case PropertySetDeclaration::Condition::Quantities: return "quantities";
      case PropertySetDeclaration::Condition::Always: return "always";
      default: return "classify";
      }
   }
}

void WriteTableAsIds(const CIfcMappingTable& table, std::ostream& os, std::vector<std::string>& notes)
{
   try
   {
      XercesGuard xercesGuard;

      const auto& files = table.GetFiles();
      IDS::info info(files.front().name + " (PGSuper IFC mapping table)");
      info.description("What PGSuper exports with the IFC mapping table \"" + files.front().name +
         "\": the property sets, quantity sets, and classifications of each element role. No design values.");
      info.date(Today());

      IDS::specificationsType specifications;
      for (const auto& [role, selector] : table.GetSelectors())
      {
         std::string role_name(GetElementRoleName(role));

         // property sets on occurrences and types (IDS checks both); material properties can't be checked through the element
         std::vector<const PropertySetDeclaration*> psets = table.GetPropertySets(role, PropertyOwner::Occurrence);
         for (const auto* pset : table.GetPropertySets(role, PropertyOwner::Type))
            psets.push_back(pset);
         for (const auto* pset : table.GetPropertySets(role, PropertyOwner::Material))
            notes.push_back("Material property set " + pset->name + " of the " + role_name + " isn't in the IDS: IDS property facets check an element's own and type property sets, not its material's");

         auto classifications = table.GetClassifications(role);
         if (psets.empty() && classifications.empty())
            continue;

         // one specification per alternative of the selector
         std::vector<const ElementSelector*> alternatives;
         if (selector.any_of.empty())
            alternatives.push_back(&selector);
         for (const auto& alternative : selector.any_of)
            alternatives.push_back(&alternative);

         for (size_t i = 0; i < alternatives.size(); i++)
         {
            const auto& a = *alternatives[i];
            IDS::applicabilityType applicability = MakeApplicability(upper(a.entity).c_str(), a.predefined_type.empty() ? nullptr : a.predefined_type.c_str(), "0");
            for (const auto& [name, value] : a.attributes)
               AddAttr(applicability, name, SimpleValue(value));
            if (!a.classification.empty())
               notes.push_back("The classification in the selector of the " + role_name + " isn't in the IDS applicability");

            IDS::requirements requirements;
            for (const auto* pset : psets)
            {
               for (const auto& property : pset->properties)
               {
                  std::vector<std::pair<std::string, std::string>> items;
                  if (property.target)
                     items.emplace_back("target", std::string(property.target->name));
                  if (property.target && !property.import)
                     items.emplace_back("import", "false");
                  if (pset->attach == PropertyOwner::Type)
                     items.emplace_back("attach", "type");
                  if (pset->condition != PropertySetDeclaration::Condition::Classify)
                     items.emplace_back("condition", condition_name(pset->condition));
                  if (pset->quantities)
                  {
                     items.emplace_back("set", "quantities");
                     items.emplace_back("method", pset->method);
                  }
                  items.emplace_back("pset_uri", pset->uri);
                  items.emplace_back("enumeration", property.enumeration_name);
                  if (property.value && !property.enumeration_values.empty())
                  {
                     // the IDS value is the constant, so the enumeration's values go in the instructions
                     std::string values;
                     for (const auto& v : property.enumeration_values)
                        values += (values.empty() ? "" : "|") + v;
                     items.emplace_back("enumeration_values", values);
                  }

                  std::optional<IDS::idsValue> value;
                  if (property.value)
                     value = SimpleValue(value_text(*property.value));
                  else if (!property.enumeration_values.empty())
                     value = EnumerationValue(property.enumeration_values);

                  Card card = (property.target || property.value) ? Card::Required : Card::Optional;
                  std::string type = pset->quantities ? measure_of_quantity(property.type) : property.type;
                  requirements.property().push_back(PropertyReq(pset->name, property.name, type.c_str(), value, card, instructions(items), property.uri));
               }
            }

            for (const auto* classification : classifications)
            {
               auto facet = ClassificationReq(classification->system, classification->identification, classification->location);
               auto text = instructions({ { "name", classification->name } });
               if (!text.empty())
                  facet.instructions(text);
               requirements.classification().push_back(facet);
            }

            std::string name = role_name + (alternatives.size() == 1 ? std::string() : " (" + std::to_string(i + 1) + ")");
            auto spec = MakeSpec(std::move(applicability), name, name, "PGSuper element role " + role_name + ": " + a.Describe(), std::move(requirements));
            spec.instructions(instructions({ { "role", role_name } }));
            specifications.specification().push_back(spec);
         }
      }

      IDS::ids document(info, specifications);
      Write(os, document);
   }
   catch (const xml_schema::exception& e)
   {
      std::ostringstream msg;
      msg << "The IDS can't be written: " << e;
      throw std::runtime_error(msg.str());
   }
}

namespace
{
   struct BindingProperty
   {
      std::string pset;
      std::string name;
      std::string target;
      json location; // unit and parse for the import (a location of the targets section)
   };

   struct BindingSpecification
   {
      std::string specification; // name or identifier in the IDS
      std::string role;
      std::vector<BindingProperty> properties;
   };

   std::vector<BindingSpecification> load_binding(const std::filesystem::path& path)
   {
      std::vector<BindingSpecification> binding;
      if (path.empty())
         return binding;

      std::ifstream file(path);
      if (!file)
         throw std::runtime_error("The binding file can't be read: " + PathToString(path));

      json j;
      try
      {
         j = json::parse(file);
      }
      catch (const json::parse_error& e)
      {
         throw std::runtime_error("The binding file " + PathToString(path) + " isn't valid JSON: " + e.what());
      }

      std::vector<std::string> errors;
      auto error = [&errors](const std::string& where, const std::string& what) { errors.push_back(where + ": " + what); };

      if (j.value("format", "") != "PGSuperIfcBinding" || j.value("version", 0) != 1)
         error("(top level)", "\"format\" must be \"PGSuperIfcBinding\" and \"version\" 1");

      if (j.contains("specifications") && j["specifications"].is_array())
      {
         for (size_t i = 0; i < j["specifications"].size(); i++)
         {
            const auto& js = j["specifications"][i];
            std::string where = "specifications[" + std::to_string(i) + "]";
            BindingSpecification spec;
            spec.specification = js.value("specification", "");
            spec.role = js.value("role", "");
            ElementKind kind;
            if (spec.specification.empty())
               error(where, "\"specification\" (the IDS specification name or identifier) is missing");
            if (!spec.role.empty() && !GetElementKind(spec.role, kind))
               error(where + ".role", "unknown element role \"" + spec.role + "\"");

            if (js.contains("properties"))
            {
               for (size_t k = 0; k < js["properties"].size(); k++)
               {
                  const auto& jp = js["properties"][k];
                  std::string pwhere = where + ".properties[" + std::to_string(k) + "]";
                  BindingProperty property;
                  property.pset = jp.value("pset", "");
                  property.name = jp.value("name", "");
                  property.target = jp.value("target", "");
                  if (property.pset.empty() || property.name.empty())
                     error(pwhere, "\"pset\" and \"name\" are required");
                  if (!property.target.empty() && property.target != "?" && !FindTargetDef(property.target))
                     error(pwhere + ".target", "unknown target \"" + property.target + "\"");
                  for (const char* key : { "unit", "parse" })
                  {
                     if (jp.contains(key))
                        property.location[key] = jp[key];
                  }
                  spec.properties.push_back(property);
               }
            }
            binding.push_back(spec);
         }
      }

      if (!errors.empty())
      {
         std::string message = "The binding file " + PathToString(path) + " has problems:";
         for (const auto& e : errors)
            message += "\n   " + e;
         throw std::runtime_error(message);
      }
      return binding;
   }

   std::optional<std::string> simple(const IDS::idsValue& value)
   {
      if (value.simpleValue())
         return std::string(*value.simpleValue());
      return std::nullopt;
   }

   // the values of an <xs:restriction> with <xs:enumeration> facets, or nothing if it's another kind of restriction
   std::vector<std::string> enumeration_values(const IDS::idsValue& value)
   {
      std::vector<std::string> values;
      if (!value.any().present())
         return values;

      const xercesc::DOMElement& restriction = value.any().get();
      for (auto* node = restriction.getFirstChild(); node; node = node->getNextSibling())
      {
         if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE)
            continue;
         auto* element = static_cast<const xercesc::DOMElement*>(node);
         char* local = xercesc::XMLString::transcode(element->getLocalName());
         std::string local_name(local ? local : "");
         xercesc::XMLString::release(&local);
         if (local_name != "enumeration")
            return {};
         char* v = xercesc::XMLString::transcode(element->getAttribute(XStr("value")));
         values.emplace_back(v ? v : "");
         xercesc::XMLString::release(&v);
      }
      return values;
   }

   json constant(const std::string& text, const std::string& type)
   {
      if (type == "IFCLABEL" || type == "IFCTEXT" || type == "IFCIDENTIFIER")
         return text;
      if (type == "IFCBOOLEAN")
         return upper(text) == "TRUE";
      try
      {
         if (type == "IFCINTEGER" || type == "IFCCOUNTMEASURE")
            return std::stoll(text);
         return std::stod(text);
      }
      catch (...)
      {
         return text;
      }
   }

   std::string role_of_selector(const CIfcMappingTable& standard, const std::string& entity, const std::string& predefined, const std::vector<std::pair<std::string, std::string>>& attributes)
   {
      for (const auto& [role, s] : standard.GetSelectors())
      {
         if (!s.any_of.empty() || !same_text(s.entity, entity) || !same_text(s.predefined_type, predefined) || s.attributes.size() != attributes.size())
            continue;
         bool bSame = true;
         for (size_t i = 0; i < attributes.size(); i++)
            bSame = bSame && s.attributes[i].first == attributes[i].first && same_text(s.attributes[i].second, attributes[i].second);
         if (bSame)
            return std::string(GetElementRoleName(role));
      }
      return "";
   }

   // the target of the standard table's property with the same property set and name
   std::string standard_target(const CIfcMappingTable& standard, ElementKind role, const std::string& pset, const std::string& name)
   {
      for (auto attach : { PropertyOwner::Occurrence, PropertyOwner::Type })
      {
         if (auto exported = standard.FindExportedProperty(role, attach, pset, name); exported && exported->property->target)
            return std::string(exported->property->target->name);
      }
      return "";
   }
}

namespace
{
   // a declared property as a table would have it
   json property_json(const PropertyDeclaration& property, bool bQuantity)
   {
      json p;
      p["name"] = property.name;
      p["type"] = IfcSpelling(property.type); // IfcReal, IfcQuantityLength, ... as tables are written
      if (property.target)
         p["target"] = std::string(property.target->name);
      if (property.value)
         std::visit([&p](const auto& v) { p["value"] = v; }, *property.value);
      if (property.target && !property.import)
         p["import"] = false;
      if (!property.uri.empty())
         p["uri"] = property.uri;
      if (!property.enumeration_name.empty())
         p["enumeration"] = { { "name", property.enumeration_name }, { "values", property.enumeration_values } };
      return p;
   }

   // Reduces a generated set to what differs from the standard table: nothing if the standard table exports the same set with
   // the same properties, the standard set with the new properties added if the IDS has more, or the set itself if it's new
   std::optional<json> delta_set(const json& set, bool bQuantity, const CIfcMappingTable& standard, std::vector<std::string>& report)
   {
      ElementKind role;
      GetElementKind(set["applies_to"].get<std::string>(), role);
      PropertyOwner attach = (set.value("attach", "occurrence") == "type") ? PropertyOwner::Type : PropertyOwner::Occurrence;

      const PropertySetDeclaration* standard_set = nullptr;
      for (const auto* pset : standard.GetPropertySets(role, attach))
      {
         if (pset->name == set["name"].get<std::string>() && pset->quantities == bQuantity)
            standard_set = pset;
      }
      if (!standard_set)
         return set; // new

      const char* items_key = bQuantity ? "quantities" : "properties";
      json merged = json::array();
      for (const auto& property : standard_set->properties)
         merged.push_back(property_json(property, bQuantity));

      bool bNew = false;
      for (const auto& p : set[items_key])
      {
         auto name = p["name"].get<std::string>();
         auto found = std::find_if(standard_set->properties.begin(), standard_set->properties.end(), [&name](const auto& sp) {return sp.name == name; });
         if (found == standard_set->properties.end())
         {
            merged.push_back(p);
            bNew = true;
         }
         else if (!bQuantity && upper(found->type) != upper(p["type"].get<std::string>()))
         {
            report.push_back("The IDS expects " + standard_set->name + "." + name + " of the " + set["applies_to"].get<std::string>() + " as " + p["type"].get<std::string>() +
               ", but the standard table exports it as " + found->type + ". The standard table's type is kept.");
         }
      }

      if (!bNew)
         return std::nullopt; // the standard table exports it already

      json result = set;
      if (standard_set->condition != PropertySetDeclaration::Condition::Classify)
         result["condition"] = condition_name(standard_set->condition);
      if (!standard_set->uri.empty())
         result["uri"] = standard_set->uri;
      if (bQuantity && !standard_set->method.empty())
         result["method"] = standard_set->method;
      result[items_key] = merged;
      return result;
   }
}

IdsToTableResult GenerateTableFromIds(const std::filesystem::path& ids_path, const std::filesystem::path& binding_path, const CIfcMappingTable& standard, const std::string& table_name, bool bExtendStandard)
{
   IdsToTableResult result;
   auto binding = load_binding(binding_path);

   XercesGuard xercesGuard;
   std::unique_ptr<IDS::ids> document;
   {
      std::ifstream file(ids_path, std::ios::binary);
      if (!file)
         throw std::runtime_error("The IDS can't be read: " + PathToString(ids_path));
      try
      {
         document = IDS::ids_(file, xml_schema::flags::dont_validate | xml_schema::flags::dont_initialize);
      }
      catch (const xml_schema::exception& e)
      {
         std::ostringstream msg;
         msg << "The IDS " << PathToString(ids_path) << " can't be read: " << e;
         throw std::runtime_error(msg.str());
      }
   }

   json table;
   table["comment"] = "Generated from the IDS " + PathToString(ids_path.filename()) + (binding_path.empty() ? std::string() : " with the binding file " + PathToString(binding_path.filename())) +
      ". Assign targets and element roles in the binding file and generate the table again, rather than editing this table.";
   table["format"] = "PGSuperIfcMapping";
   table["version"] = 1;
   table["name"] = table_name.empty() ? std::string(document->info().title()) : table_name;
   if (bExtendStandard)
      table["extends"] = "standard";
   json elements = json::object();
   json property_sets = json::array();
   json quantity_sets = json::array();
   json classification_systems = json::array();
   json classifications = json::array();
   json targets = json::object();
   std::set<std::string> systems;

   auto& report = result.report;

   for (const auto& spec : document->specifications().specification())
   {
      std::string spec_name(spec.name());
      std::string spec_identifier = spec.identifier() ? std::string(*spec.identifier()) : std::string();
      const auto& applicability = spec.applicability();
      std::string entity = applicability.entity() ? simple(applicability.entity()->name()).value_or("") : std::string();
      std::string predefined = (applicability.entity() && applicability.entity()->predefinedType()) ? simple(*applicability.entity()->predefinedType()).value_or("") : std::string();
      std::vector<std::pair<std::string, std::string>> attributes;
      for (const auto& a : applicability.attribute())
      {
         auto name = simple(a.name());
         auto value = a.value() ? simple(*a.value()) : std::nullopt;
         if (name && value)
            attributes.emplace_back(*name, *value);
      }

      // the binding file's entry for the specification
      const BindingSpecification* bound_spec = nullptr;
      for (const auto& b : binding)
      {
         if (b.specification == spec_name || (!spec_identifier.empty() && b.specification == spec_identifier))
            bound_spec = &b;
      }

      // the element role: the writer's instructions, the binding file, or the standard table's selector
      std::string role_name = parse_instructions(spec.instructions() ? std::optional<std::string>(std::string(*spec.instructions())) : std::nullopt)["role"];
      if (role_name.empty() && bound_spec)
         role_name = bound_spec->role;
      if (role_name.empty())
         role_name = role_of_selector(standard, entity, predefined, attributes);

      ElementKind role;
      if (role_name.empty() || !GetElementKind(role_name, role))
      {
         report.push_back("Specification \"" + spec_name + "\" (" + entity + (predefined.empty() ? "" : "." + predefined) + ") has no element role. " +
            "To use it, add it to the binding file: { \"specification\": \"" + spec_name + "\", \"role\": \"?\" }");
         continue;
      }

      // the selector, when it isn't the standard table's
      json selector;
      selector["entity"] = IfcSpelling(entity); // an IDS has IFCSLAB, a table IfcSlab
      if (!predefined.empty())
         selector["predefined_type"] = predefined;
      if (!attributes.empty())
      {
         selector["attributes"] = json::object();
         for (const auto& [name, value] : attributes)
            selector["attributes"][name] = value;
      }
      if (!bExtendStandard || role_of_selector(standard, entity, predefined, attributes) != role_name)
         elements[role_name] = selector;

      if (!spec.requirements())
         continue;

      const auto& requirements = *spec.requirements();
      std::map<std::string, size_t> set_index; // "quantities|attach|name" -> index in property_sets or quantity_sets
      for (const auto& p : requirements.property())
      {
         auto pset = simple(p.propertySet());
         auto name = simple(p.baseName());
         if (!pset || !name)
         {
            report.push_back("Specification \"" + spec_name + "\": a property facet with a pattern for its property set or name isn't supported, and is left out");
            continue;
         }

         std::string cardinality = p.cardinality();
         if (cardinality == "prohibited")
         {
            report.push_back("Specification \"" + spec_name + "\": " + *pset + "." + *name + " is prohibited; a mapping table can't express that, so it's left out");
            continue;
         }

         auto items = parse_instructions(p.instructions() ? std::optional<std::string>(std::string(*p.instructions())) : std::nullopt);
         const bool bPGSuperInstructions = !items.empty(); // before items[] adds keys
         std::string type = p.dataType() ? upper(std::string(*p.dataType())) : std::string("IFCLABEL");

         // a quantity: the writer says so, or a quantity set holds a measure
         std::string quantity_type = quantity_of_measure(type);
         bool bQuantity = items["set"] == "quantities" || (items["set"].empty() && pset->starts_with("Qto_") && !quantity_type.empty());

         // the target: the writer's instructions, the binding file, or the standard table's property with the same property set and name
         std::string target = items["target"];
         const BindingProperty* bound_property = nullptr;
         if (bound_spec)
         {
            for (const auto& bp : bound_spec->properties)
            {
               if (bp.pset == *pset && bp.name == *name)
                  bound_property = &bp;
            }
         }
         if (target.empty() && bound_property && bound_property->target != "?")
            target = bound_property->target;
         if (target.empty())
            target = standard_target(standard, role, *pset, *name);

         json property;
         property["name"] = *name;
         property["type"] = IfcSpelling(bQuantity ? quantity_type : type);

         bool bValue = false;
         if (p.value())
         {
            if (auto text = simple(*p.value()))
            {
               if (target.empty())
               {
                  property["value"] = constant(*text, type);
                  bValue = true;
               }

               // the values of the enumeration a constant belongs to (see WriteTableAsIds)
               if (!items["enumeration_values"].empty())
               {
                  std::vector<std::string> values;
                  std::stringstream ss(items["enumeration_values"]);
                  std::string v;
                  while (std::getline(ss, v, '|'))
                     values.push_back(v);
                  property["enumeration"] = { { "name", items["enumeration"].empty() ? "PEnum_" + *name : items["enumeration"] }, { "values", values } };
               }
            }
            else
            {
               auto values = enumeration_values(*p.value());
               if (values.empty())
                  report.push_back("Specification \"" + spec_name + "\": the value restriction of " + *pset + "." + *name + " isn't supported, and is left out (the property is kept)");
               else
                  property["enumeration"] = { { "name", items["enumeration"].empty() ? "PEnum_" + *name : items["enumeration"] }, { "values", values } };
            }
         }

         if (!target.empty())
         {
            property["target"] = target;
            if (items["import"] == "false")
               property["import"] = false;
            result.bound++;
         }
         else if (!bValue)
         {
            result.unbound++;
            if (cardinality == "required")
            {
               report.push_back("Specification \"" + spec_name + "\" (" + role_name + "): " + *pset + "." + *name + " is required but has no target, so it's exported without a value. " +
                  "To bind it, add it to the binding file: { \"pset\": \"" + *pset + "\", \"name\": \"" + *name + "\", \"target\": \"?\" }");
            }
         }

         if (p.uri() && !bQuantity)
            property["uri"] = std::string(*p.uri());
         if (p.instructions() && !bPGSuperInstructions)
            property["comment"] = std::string(*p.instructions()); // the IDS author's instructions

         // the binding's unit and parse for the import
         if (bound_property && !bound_property->location.empty() && !target.empty())
         {
            json location = bound_property->location;
            location["property"] = { { "pset", *pset }, { "name", *name } };
            targets[target].push_back(location);
         }

         // the property set (or quantity set) the property belongs to
         std::string attach = items["attach"].empty() ? "occurrence" : items["attach"];
         std::string key = std::string(bQuantity ? "q|" : "p|") + attach + "|" + *pset;
         json& sets = bQuantity ? quantity_sets : property_sets;
         auto found = set_index.find(key);
         if (found == set_index.end())
         {
            json set;
            set["name"] = *pset;
            set["applies_to"] = role_name;
            if (!bQuantity && attach != "occurrence")
               set["attach"] = attach;
            if (!items["condition"].empty())
               set["condition"] = items["condition"];
            if (!items["pset_uri"].empty())
               set["uri"] = items["pset_uri"];
            if (bQuantity && !items["method"].empty())
               set["method"] = items["method"];
            set[bQuantity ? "quantities" : "properties"] = json::array();
            sets.push_back(set);
            found = set_index.emplace(key, sets.size() - 1).first;
         }
         sets[found->second][bQuantity ? "quantities" : "properties"].push_back(property);
      }

      for (const auto& c : requirements.classification())
      {
         auto system = simple(c.system());
         auto identification = c.value() ? simple(*c.value()) : std::nullopt;
         if (!system || !identification)
         {
            report.push_back("Specification \"" + spec_name + "\": a classification facet without a system and a value isn't supported, and is left out");
            continue;
         }

         // the standard table's classifications are there already
         if (bExtendStandard)
         {
            auto existing = standard.GetClassifications(role);
            if (std::any_of(existing.begin(), existing.end(), [&](const auto* ec) {return ec->system == *system && ec->identification == *identification; }))
               continue;
         }

         json classification;
         classification["applies_to"] = role_name;
         classification["system"] = *system;
         classification["identification"] = *identification;
         auto items = parse_instructions(c.instructions() ? std::optional<std::string>(std::string(*c.instructions())) : std::nullopt);
         if (!items["name"].empty())
            classification["name"] = items["name"];
         if (c.uri())
            classification["location"] = std::string(*c.uri());
         classifications.push_back(classification);

         bool bKnown = std::any_of(standard.GetClassificationSystems().begin(), standard.GetClassificationSystems().end(), [&](const auto* s) {return s->name == *system; });
         if (!bKnown && systems.insert(*system).second)
            classification_systems.push_back({ { "name", *system } });
      }

      if (!requirements.material().empty() || !requirements.attribute().empty() || !requirements.partOf().empty() || !requirements.entity().empty())
         report.push_back("Specification \"" + spec_name + "\": material, attribute, entity, and partOf requirements aren't part of a mapping table, and are left out");
   }

   // an extending table has only what differs from the standard table
   if (bExtendStandard)
   {
      for (auto [sets, bQuantity] : { std::make_pair(&property_sets, false), std::make_pair(&quantity_sets, true) })
      {
         json reduced = json::array();
         for (const auto& set : *sets)
         {
            if (auto delta = delta_set(set, bQuantity, standard, report))
               reduced.push_back(*delta);
         }
         *sets = reduced;
      }
   }

   if (!elements.empty()) table["elements"] = elements;
   if (!property_sets.empty()) table["property_sets"] = property_sets;
   if (!quantity_sets.empty()) table["quantity_sets"] = quantity_sets;
   if (!classification_systems.empty()) table["classification_systems"] = classification_systems;
   if (!classifications.empty()) table["classifications"] = classifications;
   if (!targets.empty()) table["targets"] = targets;

   result.table_json = FormatMappingTable(table); // the style of Standard.json
   return result;
}
