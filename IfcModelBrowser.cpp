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
#include "IfcModelBrowser.h"
#include "IfcTargetReader.h"
#include "IfcTargetHints.h"
#include "Utilities.h"

#include <iomanip>
#include <sstream>

namespace
{
   // A value (the wrapped value of an IfcValue) as text
   std::string attribute_text(const ifcopenshell::attribute_value& attribute)
   {
      if (attribute.isNull())
         return "";

      std::ostringstream os;
      switch (attribute.type())
      {
      case ifcopenshell::Argument_INT:
         os << (int64_t)attribute;
         break;
      case ifcopenshell::Argument_DOUBLE:
         os << std::setprecision(10) << (double)attribute;
         break;
      case ifcopenshell::Argument_BOOL:
         os << ((bool)attribute ? "true" : "false");
         break;
      case ifcopenshell::Argument_LOGICAL:
      {
         boost::logic::tribool value = attribute;
         os << (boost::logic::indeterminate(value) ? "unknown" : ((bool)value ? "true" : "false"));
         break;
      }
      case ifcopenshell::Argument_STRING:
         os << "\"" << (std::string)attribute << "\"";
         break;
      case ifcopenshell::Argument_ENUMERATION:
         os << ((ifcopenshell::enumeration_reference)attribute).value();
         break;
      default:
         os << "(a value the editor can't show)";
         break;
      }
      return os.str();
   }

   std::string value_text(IfcSchema::IfcValue value)
   {
      return value ? attribute_text(value.get_attribute_value(0)) : std::string();
   }

   std::string value_type(IfcSchema::IfcValue value)
   {
      std::string type = value ? value.declaration().name() : std::string();
      std::transform(type.begin(), type.end(), type.begin(), ::toupper);
      return type;
   }

   // A property as a ModelProperty (value and type)
   void describe(IfcSchema::IfcProperty property, ModelProperty& mp)
   {
      mp.name = property.Name();
      if (auto single = property.as<IfcSchema::IfcPropertySingleValue>())
      {
         mp.value = value_text(single.NominalValue());
         mp.ifc_type = value_type(single.NominalValue());
         if (single.Unit())
            mp.ifc_type += " (with a unit)";
         return;
      }

      std::optional<std::vector<IfcSchema::IfcValue>> values;
      std::string kind;
      if (auto list = property.as<IfcSchema::IfcPropertyListValue>())
      {
         values = list.ListValues();
         kind = "list";
      }
      else if (auto enumerated = property.as<IfcSchema::IfcPropertyEnumeratedValue>())
      {
         values = enumerated.EnumerationValues();
         kind = "enumeration";
      }
      else
      {
         mp.ifc_type = property.declaration().name();
         mp.value = "(can't be read)";
         return;
      }

      std::string text;
      if (values)
      {
         for (size_t i = 0; i < values->size(); i++)
            text += (i == 0 ? "" : ", ") + value_text((*values)[i]);
      }
      mp.value = "[" + text + "]";
      mp.ifc_type = kind + (values && !values->empty() ? " of " + value_type(values->front()) : std::string());
   }

   // The locations of a target as the importer tries them, each once (a table can declare the same property for two roles)
   std::string locations_text(const CIfcMappingTable& table, const TargetDef& target)
   {
      std::vector<std::string> seen;
      std::string text;
      for (const auto& location : table.GetLocations(target))
      {
         auto description = location.Describe();
         if (std::find(seen.begin(), seen.end(), description) != seen.end())
            continue;
         seen.push_back(description);
         text += " " + description + ";";
      }
      return text;
   }

   std::string element_label(IfcSchema::IfcObject object)
   {
      std::ostringstream os;
      os << object.declaration().name() << " #" << object.id();
      if (auto name = object.Name(); name && !name->empty())
         os << " '" << *name << "'";
      return os.str();
   }

   std::string reading_text(const TargetDef& target, const TargetReading& reading)
   {
      return FormatTargetValue(target, reading.value) + ": " + reading.Source();
   }

   // Captures what the reader logs (e.g. values that can't be used) while it exists
   class CLogCapture
   {
   public:
      CLogCapture() { m_pOld = WBFL::System::Logger::SetOutput(&m_Stream); }
      ~CLogCapture() { WBFL::System::Logger::SetOutput(m_pOld); }
      std::string Text() const { return m_Stream.str(); }
   private:
      std::ostringstream m_Stream;
      std::ostream* m_pOld = nullptr;
   };
}

std::string FormatTargetValue(const TargetDef& target, const TargetValue& value)
{
   using namespace WBFL::Units;
   std::ostringstream os;
   if (auto number = std::get_if<Float64>(&value))
   {
      os << std::setprecision(6);
      switch (target.kind)
      {
      case ValueKind::Stress: os << ConvertFromSysUnits(*number, Measure::KSI) << " ksi (" << ConvertFromSysUnits(*number, Measure::MPa) << " MPa)"; break;
      case ValueKind::Length: os << ConvertFromSysUnits(*number, Measure::Inch) << " in (" << *number << " m)"; break;
      case ValueKind::Angle: os << ConvertFromSysUnits(*number, Measure::Degree) << " deg (" << *number << " rad)"; break;
      case ValueKind::Force: os << ConvertFromSysUnits(*number, Measure::Kip) << " kip (" << ConvertFromSysUnits(*number, Measure::Kilonewton) << " kN)"; break;
      case ValueKind::Area: os << ConvertFromSysUnits(*number, Measure::Feet2) << " ft2 (" << *number << " m2)"; break;
      case ValueKind::Mass: os << ConvertFromSysUnits(*number, Measure::PoundMass) << " lbm (" << *number << " kg)"; break;
      default: os << *number; break;
      }
   }
   else if (auto integer = std::get_if<Int64>(&value))
   {
      os << *integer;
   }
   else if (auto flag = std::get_if<bool>(&value))
   {
      os << (*flag ? "true" : "false");
   }
   else if (auto text = std::get_if<std::string>(&value))
   {
      os << "\"" << *text << "\"";
   }
   return os.str();
}

std::unique_ptr<CIfcModel> CIfcModel::Open(const std::filesystem::path& path)
{
   std::error_code ec;
   if (!std::filesystem::exists(path, ec))
      throw std::runtime_error("The model wasn't found: " + PathToString(path));

   auto extension = path.extension().wstring();
   std::transform(extension.begin(), extension.end(), extension.begin(), ::towlower);
   if (extension == L".zip")
      throw std::runtime_error("The model is in a .zip file. Extract the .ifc file from it and open that.");

   std::unique_ptr<CIfcModel> model(new CIfcModel);
   model->m_Path = path;
   model->m_pFile = std::make_unique<ifcopenshell::file>(path.string()); // IfcOpenShell takes a narrow (ANSI) path, as the importer passes it
   if (!model->m_pFile->good())
      throw std::runtime_error("The model can't be read as an IFC file: " + PathToString(path));

   std::string schema = model->m_pFile->schema()->name();
   if (schema != "IFC4X3_ADD2")
      throw std::runtime_error("The model is " + schema + ". The IFC extension reads IFC 4.3 (IFC4X3_ADD2) models.");

   CLogCapture log; // the units log the units they assume
   model->m_Units.Init(*model->m_pFile);
   return model;
}

CIfcModel::~CIfcModel() = default;

std::vector<ModelElement> CIfcModel::GetElements(const CIfcMappingTable& table, ElementKind role)
{
   std::vector<ModelElement> elements;
   CIfcTargetReader reader(table, m_Units);
   if (!reader.HasSelector(role))
      return elements;

   for (auto& object : reader.Select(role, *m_pFile))
      elements.push_back({ (int)object.id(), element_label(object) });
   return elements;
}

IfcSchema::IfcObject CIfcModel::GetObject(int id)
{
   try
   {
      return m_pFile->instance_by_id(id).as<IfcSchema::IfcObject>();
   }
   catch (...)
   {
      return {};
   }
}

std::vector<ModelProperty> CIfcModel::GetProperties(int id)
{
   std::vector<ModelProperty> properties;
   auto object = GetObject(id);
   if (!object)
      return properties;

   auto add_pset = [&properties](IfcSchema::IfcPropertySet pset, PropertyOwner owner)
   {
      for (auto& property : pset.HasProperties())
      {
         ModelProperty mp;
         mp.owner = owner;
         mp.pset = pset.Name().value_or("");
         describe(property, mp);
         properties.push_back(mp);
      }
   };

   for (auto& rel : object.IsDefinedBy())
   {
      if (auto pset = rel.RelatingPropertyDefinition().as<IfcSchema::IfcPropertySet>())
         add_pset(pset, PropertyOwner::Occurrence);
   }

   for (auto& rel : object.IsTypedBy())
   {
      if (auto psets = rel.RelatingType().HasPropertySets())
      {
         for (auto& definition : *psets)
         {
            if (auto pset = definition.as<IfcSchema::IfcPropertySet>())
               add_pset(pset, PropertyOwner::Type);
         }
      }
   }

   if (auto material = GetMaterial<IfcSchema>(object))
   {
      for (auto& material_properties : material.HasProperties())
      {
         for (auto& property : material_properties.Properties())
         {
            ModelProperty mp;
            mp.owner = PropertyOwner::Material;
            mp.pset = material_properties.Name().value_or("");
            describe(property, mp);
            properties.push_back(mp);
         }
      }
   }
   return properties;
}

std::string CIfcModel::DescribeReading(const CIfcMappingTable& table, std::string_view target_name, int id)
{
   const TargetDef* target = FindTargetDef(target_name);
   auto object = GetObject(id);
   if (!target || !object)
      return "";

   CLogCapture log;
   CIfcTargetReader reader(table, m_Units);
   std::ostringstream os;
   if (target->kind == ValueKind::TextList)
   {
      auto readings = reader.ReadAll(target->name, object);
      if (readings.empty())
         os << "The table doesn't find it.";
      for (const auto& reading : readings)
         os << "The table reads " << reading_text(*target, reading) << std::endl;
   }
   else if (auto reading = reader.Read(target->name, object))
   {
      os << "The table reads " << reading_text(*target, *reading);
   }
   else
   {
      os << "The table doesn't find it.";
   }

   auto locations = table.GetLocations(*target);
   if (locations.empty())
   {
      os << std::endl << "The table has no locations for it.";
   }
   else
   {
      os << std::endl << "Locations, in order:" << locations_text(table, *target);
   }

   auto hints = CIfcTargetHints::Find(*target, object);
   if (!hints.empty())
   {
      os << std::endl << "Hints (properties that look like they may hold it):";
      for (const auto& hint : hints)
         os << " " << hint << ";";
   }

   auto problems = log.Text();
   if (!problems.empty())
      os << std::endl << problems;
   return os.str();
}

std::string CIfcModel::TryTable(const CIfcMappingTable& table)
{
   CLogCapture log;
   CIfcTargetReader reader(table, m_Units);
   std::ostringstream os;
   os << "Model: " << PathToString(m_Path) << std::endl;
   os << "Table: \"" << table.GetFiles().front().name << "\"";
   for (size_t i = 1; i < table.GetFiles().size(); i++)
      os << ", which extends \"" << table.GetFiles()[i].name << "\"";
   os << std::endl;

   std::vector<std::string> not_selected;
   for (int kind = (int)ElementKind::Project; kind <= (int)ElementKind::Barrier; kind++)
   {
      ElementKind role = (ElementKind)kind;
      std::string role_name(GetElementRoleName(role));

      std::vector<const TargetDef*> targets;
      for (const auto& target : GetTargetDefs())
      {
         if (target.element == GetTargetElement(role) && !table.GetLocations(target).empty())
            targets.push_back(&target);
      }
      if (targets.empty() || !reader.HasSelector(role))
         continue;

      auto objects = reader.Select(role, *m_pFile);
      if (objects.empty())
      {
         not_selected.push_back(role_name + " (" + table.GetSelector(role)->Describe() + ")");
         continue;
      }

      os << std::endl << role_name << ": " << objects.size() << (objects.size() == 1 ? " element" : " elements") << ", e.g. " << element_label(objects.front()) << std::endl;
      for (const auto* target : targets)
      {
         size_t found = 0;
         std::optional<TargetReading> first;
         for (auto& object : objects)
         {
            if (auto reading = reader.Read(target->name, object))
            {
               if (!first)
                  first = reading;
               found++;
            }
         }

         os << "   " << target->name << ": ";
         if (first)
         {
            os << "found for " << found << " of " << objects.size() << ". First: " << reading_text(*target, *first) << std::endl;
         }
         else
         {
            os << "not found. Tried:" << locations_text(table, *target);
            auto hints = CIfcTargetHints::Find(*target, objects.front());
            if (!hints.empty())
            {
               os << " Hints:";
               for (const auto& hint : hints)
                  os << " " << hint << ";";
            }
            os << std::endl;
         }
      }
   }

   if (!not_selected.empty())
   {
      os << std::endl << "Element roles with targets but no elements in the model:" << std::endl;
      for (const auto& role : not_selected)
         os << "   " << role << std::endl;
   }

   auto problems = log.Text();
   if (!problems.empty())
      os << std::endl << "Problems found while reading:" << std::endl << problems;
   return os.str();
}
