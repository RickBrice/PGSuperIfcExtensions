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
#include "IfcTargetReader.h"
#include "IfcTargetHints.h"
#include "Utilities.h"
#include "USBridge_Classifications.h"

#include <charconv>

namespace
{
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

   std::string trim(const std::string& s)
   {
      auto begin = s.find_first_not_of(" \t\r\n");
      if (begin == std::string::npos)
         return "";
      auto end = s.find_last_not_of(" \t\r\n");
      return s.substr(begin, end - begin + 1);
   }

   // A value as found in the model
   struct RawValue
   {
      std::string ifc_type; // upper case IFC type of the value (e.g. IFCPRESSUREMEASURE). Empty for attribute values
      std::variant<std::monostate, int64_t, double, bool, std::string> value;
      IfcSchema::IfcUnit unit; // the property's own unit

      std::string Text() const
      {
         std::ostringstream os;
         std::visit([&os](const auto& v)
            {
               using T = std::decay_t<decltype(v)>;
               if constexpr (std::is_same_v<T, bool>)
                  os << (v ? "true" : "false");
               else if constexpr (!std::is_same_v<T, std::monostate>)
                  os << v;
            }, value);
         return os.str();
      }

      bool IsNumber() const { return std::holds_alternative<int64_t>(value) || std::holds_alternative<double>(value); }
      double Number() const { return std::holds_alternative<int64_t>(value) ? (double)std::get<int64_t>(value) : std::get<double>(value); }
   };

   // Converts an attribute value (of an IfcValue's wrapped value, or of an entity attribute) to a RawValue
   bool to_raw(const ifcopenshell::attribute_value& attribute, RawValue& raw)
   {
      if (attribute.isNull())
         return false;

      switch (attribute.type())
      {
      case ifcopenshell::Argument_INT:
         raw.value = (int64_t)attribute;
         return true;
      case ifcopenshell::Argument_DOUBLE:
         raw.value = (double)attribute;
         return true;
      case ifcopenshell::Argument_BOOL:
         raw.value = (bool)attribute;
         return true;
      case ifcopenshell::Argument_LOGICAL:
      {
         boost::logic::tribool value = attribute;
         if (boost::logic::indeterminate(value))
            return false;
         raw.value = (bool)value;
         return true;
      }
      case ifcopenshell::Argument_STRING:
         raw.value = (std::string)attribute;
         return true;
      case ifcopenshell::Argument_ENUMERATION:
         raw.value = std::string(((ifcopenshell::enumeration_reference)attribute).value());
         return true;
      default:
         return false;
      }
   }

   bool to_raw(IfcSchema::IfcValue value, IfcSchema::IfcUnit unit, RawValue& raw)
   {
      if (!value)
         return false;

      raw.ifc_type = upper(value.declaration().name());
      raw.unit = unit;
      return to_raw(value.get_attribute_value(0), raw);
   }

   // The value of a property. List and enumerated values give the element at list_index (default: the first)
   bool property_value(IfcSchema::IfcProperty property, const MappingLocation& location, RawValue& raw, std::string& problem)
   {
      if (auto single = property.as<IfcSchema::IfcPropertySingleValue>())
         return to_raw(single.NominalValue(), single.Unit(), raw);

      std::optional<std::vector<IfcSchema::IfcValue>> values;
      IfcSchema::IfcUnit unit;
      if (auto list = property.as<IfcSchema::IfcPropertyListValue>())
      {
         values = list.ListValues();
         unit = list.Unit();
      }
      else if (auto enumerated = property.as<IfcSchema::IfcPropertyEnumeratedValue>())
      {
         values = enumerated.EnumerationValues();
      }
      else
      {
         problem = "is a " + property.declaration().name() + ", which can't be read";
         return false;
      }

      if (!values || values->empty())
         return false;

      size_t index = location.list_index.value_or(0);
      if (values->size() <= index)
      {
         problem = "has " + std::to_string(values->size()) + " values, so there is no value " + std::to_string(index);
         return false;
      }

      return to_raw((*values)[index], unit, raw);
   }

   IfcSchema::IfcProperty find_property(const std::vector<IfcSchema::IfcProperty>& properties, const std::string& name)
   {
      for (auto& property : properties)
      {
         if (property.Name() == name)
            return property;
      }
      return {};
   }

   // The property of a property location
   IfcSchema::IfcProperty find_property(IfcSchema::IfcObject object, const MappingLocation& location)
   {
      if (location.on == PropertyOwner::Material)
      {
         auto material = GetMaterial<IfcSchema>(object);
         if (!material)
            return {};

         for (auto& material_properties : material.HasProperties())
         {
            if (material_properties.Name() == location.pset)
            {
               if (auto property = find_property(material_properties.Properties(), location.name))
                  return property;
            }
         }
         return {};
      }

      // Property sets of the object override the property sets of its type (IFC 4.3 5.1.3.6 IfcObject)
      if (location.on == PropertyOwner::Occurrence)
      {
         for (auto& rel : object.IsDefinedBy())
         {
            auto pset = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySet>();
            if (pset && pset.Name() == location.pset)
            {
               if (auto property = find_property(pset.HasProperties(), location.name))
                  return property;
            }
         }
      }

      for (auto& rel : object.IsTypedBy())
      {
         auto psets = rel.RelatingType().HasPropertySets();
         if (!psets)
            continue;

         for (auto& definition : *psets)
         {
            auto pset = definition.as<IfcSchema::IfcPropertySet>();
            if (pset && pset.Name() == location.pset)
            {
               if (auto property = find_property(pset.HasProperties(), location.name))
                  return property;
            }
         }
      }

      return {};
   }

   bool attribute_value(express::entity entity, const std::string& name, RawValue& raw)
   {
      if (!entity)
         return false;

      try
      {
         return to_raw(entity.get(name), raw);
      }
      catch (...)
      {
         return false; // the entity doesn't have this attribute
      }
   }

   // The values found at a location
   std::vector<RawValue> find_values(IfcSchema::IfcObject object, const MappingLocation& location, std::string& problem)
   {
      std::vector<RawValue> values;
      RawValue raw;
      switch (location.kind)
      {
      case MappingLocation::Kind::Property:
         if (auto property = find_property(object, location); property && property_value(property, location, raw, problem))
            values.push_back(raw);
         break;

      case MappingLocation::Kind::Attribute:
         if (attribute_value(object.as<express::entity>(), location.name, raw))
            values.push_back(raw);
         break;

      case MappingLocation::Kind::TypeAttribute:
         for (auto& rel : object.IsTypedBy())
         {
            if (attribute_value(rel.RelatingType().as<express::entity>(), location.name, raw))
            {
               values.push_back(raw);
               break;
            }
         }
         break;

      case MappingLocation::Kind::Classification:
         for (auto& rel : object.HasAssociations())
         {
            auto rel_classification = rel.as<IfcSchema::IfcRelAssociatesClassification>();
            auto reference = rel_classification ? rel_classification.RelatingClassification().as<IfcSchema::IfcClassificationReference>() : IfcSchema::IfcClassificationReference{};
            if (!reference)
               continue;

            if (!location.classification_identification.empty() && reference.Identification().value_or("") != location.classification_identification)
               continue;

            if (!location.classification_system.empty())
            {
               auto classification = reference.ReferencedSource().as<IfcSchema::IfcClassification>();
               if (!classification || classification.Name() != location.classification_system)
                  continue;
            }

            auto text = location.classification_field_is_name ? reference.Name() : reference.Identification();
            if (text)
            {
               RawValue value;
               value.value = *text;
               values.push_back(value);
            }
         }
         break;
      }
      return values;
   }

   bool parse_number(const std::string& text, double& value)
   {
      auto s = trim(text);
      if (s.empty())
         return false;
      auto result = std::from_chars(s.data(), s.data() + s.size(), value);
      return result.ec == std::errc() && result.ptr == s.data() + s.size();
   }

   // Parses a length such as 8'-6", 8' 6 1/2", 8.5", or 6'. Returns false if the text isn't a length in feet and inches.
   // bMarked is false if the text is a number without ' or " (its unit is then the location unit)
   bool parse_feet_inches(const std::string& text, double& inches, bool& bMarked)
   {
      static const std::regex pattern(R"(^\s*(?:(\d+(?:\.\d+)?)\s*')?\s*-?\s*(?:(\d+(?:\.\d+)?)(?:\s+(\d+)\s*/\s*(\d+))?\s*")?\s*$)");
      std::smatch match;
      if (!std::regex_match(text, match, pattern) || (!match[1].matched && !match[2].matched))
      {
         bMarked = false;
         return false;
      }

      bMarked = true;
      inches = 0;
      if (match[1].matched)
         inches += 12.0 * std::stod(match[1].str());
      if (match[2].matched)
         inches += std::stod(match[2].str());
      if (match[3].matched && std::stod(match[4].str()) != 0)
         inches += std::stod(match[3].str()) / std::stod(match[4].str());
      return true;
   }

   struct MeasureInfo
   {
      ValueKind kind;
      IfcSchema::IfcUnitEnum::Value unit_type;
   };

   // IFC measure types that hold a number of a value kind with a unit
   std::optional<MeasureInfo> get_measure_info(const std::string& ifc_type)
   {
      using U = IfcSchema::IfcUnitEnum;
      static const std::map<std::string, MeasureInfo> measures{
         { "IFCPRESSUREMEASURE", { ValueKind::Stress, U::IfcUnit_PRESSUREUNIT } },
         { "IFCMODULUSOFELASTICITYMEASURE", { ValueKind::Stress, U::IfcUnit_PRESSUREUNIT } },
         { "IFCLENGTHMEASURE", { ValueKind::Length, U::IfcUnit_LENGTHUNIT } },
         { "IFCPOSITIVELENGTHMEASURE", { ValueKind::Length, U::IfcUnit_LENGTHUNIT } },
         { "IFCNONNEGATIVELENGTHMEASURE", { ValueKind::Length, U::IfcUnit_LENGTHUNIT } },
         { "IFCPLANEANGLEMEASURE", { ValueKind::Angle, U::IfcUnit_PLANEANGLEUNIT } },
         { "IFCPOSITIVEPLANEANGLEMEASURE", { ValueKind::Angle, U::IfcUnit_PLANEANGLEUNIT } },
         { "IFCFORCEMEASURE", { ValueKind::Force, U::IfcUnit_FORCEUNIT } },
      };
      auto found = measures.find(ifc_type);
      return found == measures.end() ? std::nullopt : std::make_optional(found->second);
   }

   // Converts value, in the IFC unit (or project unit) for unit_type, to system units
   Float64 convert_measure(const CIfcImportUnits& units, ValueKind kind, IfcSchema::IfcUnitEnum::Value unit_type, Float64 value, IfcSchema::IfcUnit unit)
   {
      Float64 factor = units.GetConversionFactor(unit, unit_type);
      switch (kind)
      {
      case ValueKind::Stress: return WBFL::Units::ConvertToSysUnits(value, WBFL::Units::Pressure(factor, _T("IFC")));
      case ValueKind::Length: return WBFL::Units::ConvertToSysUnits(value, WBFL::Units::Length(factor, _T("IFC")));
      case ValueKind::Angle: return WBFL::Units::ConvertToSysUnits(value, WBFL::Units::Angle(factor, _T("IFC")));
      case ValueKind::Force: return WBFL::Units::ConvertToSysUnits(value, WBFL::Units::Force(factor, _T("IFC")));
      default: ASSERT(false); return value;
      }
   }

   bool is_text_type(const std::string& ifc_type)
   {
      return ifc_type.empty() || ifc_type == "IFCLABEL" || ifc_type == "IFCTEXT" || ifc_type == "IFCIDENTIFIER";
   }

   // IFC types that hold a plain number (no unit): the location's unit applies
   bool is_plain_number_type(const std::string& ifc_type)
   {
      return ifc_type.empty() || ifc_type == "IFCREAL" || ifc_type == "IFCINTEGER" || ifc_type == "IFCNUMERICMEASURE" ||
         ifc_type == "IFCCOUNTMEASURE" || ifc_type == "IFCPOSITIVEINTEGER";
   }

   // Converts a value found in the model to a target value. Returns false, with the reason in problem, if the value can't be used
   bool convert(const CIfcImportUnits& units, const RawValue& raw, const MappingLocation& location, const TargetDef& target, TargetValue& result, std::string& problem)
   {
      if (!location.value_types.empty() && !raw.ifc_type.empty() &&
         std::find(location.value_types.begin(), location.value_types.end(), raw.ifc_type) == location.value_types.end())
      {
         problem = "is a " + raw.ifc_type + ", which isn't one of the location's value types";
         return false;
      }

      std::string text = raw.Text();
      if (location.parse == MappingLocation::Parse::Regex)
      {
         std::smatch match;
         if (!std::regex_search(text, match, location.regex) || match.size() <= location.regex_group || !match[location.regex_group].matched)
         {
            problem = "doesn't match the location's regular expression";
            return false;
         }
         text = match[location.regex_group].str();
      }

      if (!location.map.empty())
      {
         auto key = lower(trim(text));
         auto found = std::find_if(location.map.begin(), location.map.end(), [&key](const auto& item) {return item.first == key; });
         if (found == location.map.end())
         {
            problem = "isn't one of the values in the location's map";
            return false;
         }
         result = found->second;
         return true;
      }

      switch (target.kind)
      {
      case ValueKind::Boolean:
         if (!std::holds_alternative<bool>(raw.value))
         {
            problem = "isn't a boolean. Add a \"map\" to the location to say which values mean true and false";
            return false;
         }
         result = std::get<bool>(raw.value);
         return true;

      case ValueKind::Text:
      case ValueKind::TextList:
         text = trim(text);
         if (text.empty())
         {
            problem = "is empty";
            return false;
         }
         result = text;
         return true;

      default:
         break;
      }

      // numeric targets
      ASSERT(IsNumeric(target.kind));
      double number = 0;
      bool bUnitApplied = false; // true if number is in system units

      auto measure = get_measure_info(raw.ifc_type);
      if (measure && location.parse != MappingLocation::Parse::Regex)
      {
         if (measure->kind != target.kind)
         {
            problem = "is a " + raw.ifc_type + ", which isn't a " + std::string(GetValueKindName(target.kind));
            return false;
         }
         // a measure with a unit: its own unit or the project unit. It always wins over the location's unit
         number = convert_measure(units, target.kind, measure->unit_type, raw.Number(), raw.unit);
         bUnitApplied = true;
      }
      else if (raw.IsNumber() && location.parse != MappingLocation::Parse::Regex)
      {
         if (!is_plain_number_type(raw.ifc_type) && HasUnit(target.kind))
         {
            problem = "is a " + raw.ifc_type + ", which isn't a " + std::string(GetValueKindName(target.kind));
            return false;
         }
         number = raw.Number();
      }
      else if (std::holds_alternative<std::string>(raw.value) || location.parse == MappingLocation::Parse::Regex)
      {
         if (!is_text_type(raw.ifc_type) && location.parse != MappingLocation::Parse::Regex)
         {
            problem = "is a " + raw.ifc_type + ", which isn't a " + std::string(GetValueKindName(target.kind));
            return false;
         }

         if (location.parse == MappingLocation::Parse::FeetInches)
         {
            double inches;
            bool bMarked;
            if (parse_feet_inches(text, inches, bMarked))
            {
               number = WBFL::Units::ConvertToSysUnits(inches, WBFL::Units::Measure::Inch);
               bUnitApplied = true;
            }
            else if (!parse_number(text, number))
            {
               problem = "isn't a length in feet and inches (e.g. 8'-6\" or 8.5\")";
               return false;
            }
         }
         else if (!parse_number(text, number))
         {
            problem = "isn't a number. If the number is part of the text, use \"parse\": { \"regex\": ... } to pick it out";
            return false;
         }
      }
      else
      {
         problem = "isn't a number";
         return false;
      }

      if (HasUnit(target.kind) && !bUnitApplied)
      {
         if (!location.unit)
         {
            problem = "has no unit. Add \"unit\" to the location (e.g. \"ksi\")";
            return false;
         }
         number = location.unit->to_system_units(number);
      }

      if (target.kind == ValueKind::Count)
      {
         if (fabs(number - std::round(number)) > 1.0e-9)
         {
            problem = "isn't a whole number";
            return false;
         }
         result = (Int64)std::llround(number);
      }
      else
      {
         result = (Float64)number;
      }
      return true;
   }

   std::string object_name(IfcSchema::IfcObject object)
   {
      std::ostringstream os;
      os << object.declaration().name() << " #" << object.id();
      if (object.Name())
         os << " \"" << *object.Name() << "\"";
      return os.str();
   }

   std::string predefined_type(IfcSchema::IfcObject object)
   {
      // the predefined type is from the type object if there is one (see GetPredefinedType)
      RawValue raw;
      for (auto& rel : object.IsTypedBy())
      {
         if (attribute_value(rel.RelatingType().as<express::entity>(), "PredefinedType", raw))
            return upper(raw.Text());
      }

      if (attribute_value(object.as<express::entity>(), "PredefinedType", raw))
         return upper(raw.Text());

      return "";
   }
}

std::string TargetReading::Source() const
{
   std::ostringstream os;
   os << "\"" << raw << "\" from " << location->Describe() << " in mapping table \"" << location->table << "\"";
   return os.str();
}

CIfcTargetReader::CIfcTargetReader(const CIfcMappingTable& table, const CIfcImportUnits& units) :
   m_Table(table),
   m_Units(units)
{
}

const TargetDef& CIfcTargetReader::GetTarget(std::string_view target) const
{
   auto* def = FindTargetDef(target);
   ASSERT(def); // targets are named in code, so this is a programming error
   if (!def)
      throw std::invalid_argument("unknown target " + std::string(target));
   return *def;
}

std::optional<TargetReading> CIfcTargetReader::Read(std::string_view target, IfcSchema::IfcObject object) const
{
   std::vector<TargetReading> readings;
   Read(GetTarget(target), object, false, readings);
   return readings.empty() ? std::nullopt : std::make_optional(readings.front());
}

std::vector<TargetReading> CIfcTargetReader::ReadAll(std::string_view target, IfcSchema::IfcObject object) const
{
   std::vector<TargetReading> readings;
   Read(GetTarget(target), object, true, readings);
   return readings;
}

void CIfcTargetReader::Read(const TargetDef& target, IfcSchema::IfcObject object, bool bAll, std::vector<TargetReading>& readings) const
{
   if (!object)
      return;

   for (const auto& location : m_Table.GetLocations(target))
   {
      std::string problem;
      auto values = find_values(object, location, problem);
      if (!problem.empty())
         LogProblem(location.Describe() + " of " + object_name(object) + " " + problem + ". It isn't used for " + std::string(target.description) + " (" + location.origin + ")");

      for (const auto& raw : values)
      {
         TargetValue value;
         if (!convert(m_Units, raw, location, target, value, problem))
         {
            LogProblem("The value \"" + raw.Text() + "\" of " + location.Describe() + " of " + object_name(object) + " " + problem + ". It isn't used for " + std::string(target.description) + " (" + location.origin + ")");
            continue;
         }

         bool bDuplicate = std::any_of(readings.begin(), readings.end(), [&value](const auto& reading) {return reading.value == value; });
         if (!bDuplicate)
            readings.push_back({ value, &location, raw.Text() });

         if (!bAll)
            return;
      }
   }
}

void CIfcTargetReader::ReportNotFound(std::string_view target, IfcSchema::IfcObject object, const std::string& element_name) const
{
   const auto& def = GetTarget(target);
   const auto& locations = m_Table.GetLocations(def);

   if (1 < ++m_NotFound[std::string(target)])
      return; // only the first element is logged

   std::ostringstream os;
   os << "The " << def.description << " of " << element_name << " wasn't found. ";
   if (locations.empty())
   {
      os << "The mapping table has no locations for " << def.name << ".";
   }
   else
   {
      os << "Locations tried: ";
      for (const auto& location : locations)
         os << location.Describe() << (&location == &locations.back() ? "." : "; ");
   }
   WBFL::System::Logger::Info(os.str());

   auto hints = CIfcTargetHints::Find(def, object);
   if (!hints.empty())
   {
      std::ostringstream hint;
      hint << "Hint: properties of " << element_name << " that may hold the " << def.description << ": ";
      for (const auto& h : hints)
         hint << h << (&h == &hints.back() ? ". " : "; ");
      hint << "They aren't used. To use one, add it to the locations of " << def.name << " in the mapping table.";
      WBFL::System::Logger::Info(hint.str());
   }
}

void CIfcTargetReader::LogNotFoundSummary() const
{
   for (const auto& [target, count] : m_NotFound)
   {
      if (count < 2)
         continue;

      std::ostringstream os;
      os << "The " << GetTarget(target).description << " also wasn't found for " << (count - 1) << " more " << (count == 2 ? "element." : "elements.");
      WBFL::System::Logger::Info(os.str());
   }
}

bool CIfcTargetReader::HasSelector(ElementKind role) const
{
   return m_Table.GetSelector(role) != nullptr;
}

bool CIfcTargetReader::Matches(ElementKind role, IfcSchema::IfcObject object) const
{
   auto* selector = m_Table.GetSelector(role);
   return selector && object && Matches(*selector, object);
}

bool CIfcTargetReader::Matches(const ElementSelector& selector, IfcSchema::IfcObject object) const
{
   if (!selector.any_of.empty())
      return std::any_of(selector.any_of.begin(), selector.any_of.end(), [this, &object](const auto& alternative) {return Matches(alternative, object); });

   if (!object.declaration().is(selector.entity))
      return false;

   if (!selector.predefined_type.empty() && predefined_type(object) != selector.predefined_type)
      return false;

   for (const auto& [name, value] : selector.attributes)
   {
      RawValue raw;
      if (!attribute_value(object.as<express::entity>(), name, raw) || lower(raw.Text()) != lower(value))
         return false;
   }

   if (!selector.classification.empty() && !HasClassification<IfcSchema>(object, selector.classification))
      return false;

   return true;
}

std::vector<IfcSchema::IfcObject> CIfcTargetReader::Select(ElementKind role, ifcopenshell::file& file) const
{
   std::vector<IfcSchema::IfcObject> objects;
   auto* selector = m_Table.GetSelector(role);
   if (!selector)
      return objects;

   // the entities to look through
   std::vector<std::string> entities;
   std::function<void(const ElementSelector&)> collect = [&](const ElementSelector& s)
      {
         if (s.any_of.empty())
         {
            if (std::find(entities.begin(), entities.end(), s.entity) == entities.end())
               entities.push_back(s.entity);
         }
         for (const auto& alternative : s.any_of)
            collect(alternative);
      };
   collect(*selector);

   std::set<uint32_t> ids;
   for (const auto& entity : entities)
   {
      for (auto& instance : file.instances_by_type(entity))
      {
         auto object = instance.as<IfcSchema::IfcObject>();
         if (object && Matches(*selector, object) && ids.insert(object.id()).second)
            objects.push_back(object);
      }
   }

   std::sort(objects.begin(), objects.end(), [](const auto& a, const auto& b) {return a.id() < b.id(); });
   return objects;
}

void CIfcTargetReader::LogProblem(const std::string& problem) const
{
   if (m_Problems.insert(problem).second)
      WBFL::System::Logger::Warning(problem);
}
