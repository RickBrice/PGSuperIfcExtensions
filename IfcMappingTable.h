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
CLASS
   CIfcMappingTable

   A mapping table binds targets (PGSuper data items, see IfcTargets.h) to
   IFC locations, and says which IFC elements play which element role.
   Tables are JSON files. A table can extend another table (e.g. an agency
   table extends the standard table); the extending table's import locations
   are tried first, and its element roles replace the base table's.

   See devdocs/MappingTablesDesign.md for the format.
*****************************************************************************/

#include "IfcTargets.h"

#include <filesystem>
#include <map>
#include <optional>
#include <regex>
#include <stdexcept>

// How the path of a table file was chosen
enum class MappingTableSource
{
   CommandLine,          // /IfcMapping=<file>
   ConfigurationSetting, // the mapping table setting (Options > IFC Mapping Table, registry)
   InstalledStandard,    // the standard table installed with the extension
   Extends               // named by the "extends" of another table
};

// How a table file was chosen, for logs and messages (e.g. "the command line (/IfcMapping)")
std::string MappingTableSourceDescription(MappingTableSource source);

// A unit a table can state for a model value that has no unit of its own (e.g. "ksi")
struct TableUnit
{
   std::string_view name;
   ValueKind kind;
   Float64(*to_system_units)(Float64 value);
};

// The unit with the given name, or nullptr
const TableUnit* FindTableUnit(std::string_view name);

// Where the properties of a location are looked for
enum class PropertyOwner
{
   Occurrence, // property sets of the element, then of its type object (as IfcObject.IsDefinedBy/IsTypedBy)
   Type,       // property sets of the type object only
   Material    // material properties of the element's material
};

// A place in an IFC model where a target value may be found
struct MappingLocation
{
   enum class Kind { Property, Attribute, TypeAttribute, Classification };
   enum class Parse { None, Number, FeetInches, Regex };

   Kind kind = Kind::Property;
   std::string pset;           // Property: property set name
   std::string name;           // Property: property name; Attribute and TypeAttribute: attribute name
   PropertyOwner on = PropertyOwner::Occurrence;
   std::string classification_system;          // Classification: name of the IfcClassification (empty: any)
   std::string classification_identification;  // Classification: identification of the reference (empty: any)
   bool classification_field_is_name = true;   // Classification: read Name (true) or Identification (false)
   std::vector<std::string> value_types;       // accepted IFC value types, upper case (empty: any that fits the target)
   const TableUnit* unit = nullptr;            // unit of a model value that has no unit of its own
   Parse parse = Parse::None;                  // text -> value
   std::regex regex;                           // Parse::Regex
   size_t regex_group = 1;
   std::optional<size_t> list_index;           // element of a list or enumerated value
   std::vector<std::pair<std::string, TargetValue>> map; // model value (lower case) -> target value

   std::string table;   // name of the table that defined the location
   std::string origin;  // file and JSON path, for messages

   // e.g. "IaDOT_PPCB.5_Final Concrete Strength, Fc [ksi]"
   std::string Describe() const;
};

// Which IFC elements play an element role
struct ElementSelector
{
   std::string entity;         // IFC entity, subtypes included (e.g. "IfcSlab")
   std::string predefined_type; // upper case, from the type object if there is one (empty: any)
   std::vector<std::pair<std::string, std::string>> attributes; // attribute name, value (exact, ignoring case)
   std::string classification; // identification of a classification reference (empty: any)
   std::vector<ElementSelector> any_of; // alternatives. When not empty, the other members are empty

   std::string table; // name of the table that defined the selector
   std::string Describe() const;
};

// A property of a property set declared for export. When it has a target, it's also an import location (unless import is false).
// A property with neither a target nor a value is written without a value
struct PropertyDeclaration
{
   std::string name;
   std::string type; // IFC value type, upper case (e.g. "IFCPRESSUREMEASURE")
   const TargetDef* target = nullptr;
   std::optional<TargetValue> value; // a constant, written as it is
   bool import = true;
   std::string uri; // the property's Specification: a bSDD property URI, an explicit URI, or empty
   std::string enumeration_name;
   std::vector<std::string> enumeration_values;
};

struct PropertySetDeclaration
{
   // the export option that includes the property set
   enum class Condition
   {
      Classify,   // usBridge classification (CIfcExportOptions::classify), the default
      Quantities, // quantities (CIfcExportOptions::include_quantities)
      Always
   };

   std::string name;
   ElementKind applies_to = ElementKind::Bridge;
   PropertyOwner attach = PropertyOwner::Occurrence;
   Condition condition = Condition::Classify;
   std::string uri; // the property set's Description: a bSDD class URI, an explicit URI, or empty
   bool quantities = false; // an IfcElementQuantity: the properties are quantities (type IfcQuantityLength, ...)
   std::string method; // quantity set MethodOfMeasurement (e.g. "BaseQuantities")
   bool shared = false;
   bool remove = false; // an extending table removes the base table's property set with the same name, element role, and attach
   std::vector<PropertyDeclaration> properties;
};

// True if the IFC value type can be written by the exporter (e.g. "IFCLABEL", "IFCPRESSUREMEASURE")
bool IsExportValueType(const std::string& type);

// True if the quantity type can be written by the exporter (e.g. "IFCQUANTITYLENGTH")
bool IsExportQuantityType(const std::string& type);

// A classification system (IfcClassification), associated with the project when it's exported
struct ClassificationSystemDeclaration
{
   std::string name;
   std::string source;
   std::string edition;
   std::string edition_date;
   std::string specification; // a URI, or empty
};

// A classification reference of an element role
struct ClassificationDeclaration
{
   ElementKind applies_to = ElementKind::Bridge;
   PropertySetDeclaration::Condition condition = PropertySetDeclaration::Condition::Classify;
   std::string system; // name of the classification system
   std::string identification;
   std::string name;
   std::string location; // a URI, or empty
   bool remove = false; // an extending table removes the base table's classification with the same element role and identification
};

// One table file
struct MappingTableFile
{
   std::filesystem::path path;
   MappingTableSource source = MappingTableSource::InstalledStandard;
   std::string name;
   int version = 0;
   std::string extends;

   std::map<std::string, std::vector<MappingLocation>, std::less<>> targets; // explicit import locations
   std::set<std::string, std::less<>> replace_targets; // targets whose locations replace those of the base table
   std::map<ElementKind, ElementSelector> elements;
   std::vector<PropertySetDeclaration> property_sets; // and quantity sets
   std::vector<ClassificationSystemDeclaration> classification_systems;
   std::vector<ClassificationDeclaration> classifications;
};

class CIfcMappingTableException : public std::runtime_error
{
public:
   using std::runtime_error::runtime_error;
};

class CIfcMappingTable
{
public:
   // Loads a table and the tables it extends. An empty path loads the installed standard table.
   // Throws CIfcMappingTableException with a message for users and support staff if a table can't be found, read, or used.
   static std::unique_ptr<CIfcMappingTable> Load(const std::filesystem::path& path, MappingTableSource source);

   // Loads the table an import or export uses, first found wins: command_line_file (/IfcMapping=, if not empty),
   // the mapping table setting, the installed standard table. Throws CIfcMappingTableException like Load
   static std::unique_ptr<CIfcMappingTable> LoadActive(const std::filesystem::path& command_line_file);

   // The standard table installed with the extension (MappingTables\Standard.json next to the extension DLL)
   static std::filesystem::path GetStandardTablePath();

   // The mapping table setting (Options > IFC Mapping Table), stored in the registry for all BridgeLink applications.
   // Empty means the installed standard table
   static std::filesystem::path GetTableSetting();
   static void SetTableSetting(const std::filesystem::path& path);

   // Import locations of a target, in the order they are tried
   const std::vector<MappingLocation>& GetLocations(const TargetDef& target) const;

   // The selector for an element role, or nullptr if no table defines it
   const ElementSelector* GetSelector(ElementKind role) const;

   // The selectors of all element roles the tables define
   const std::map<ElementKind, ElementSelector>& GetSelectors() const { return m_Selectors; }

   // The property sets exported for an element role and owner, in table order. An extending table's property set
   // replaces the base table's property set with the same name, element role, and attach, or removes it
   std::vector<const PropertySetDeclaration*> GetPropertySets(ElementKind role, PropertyOwner attach) const;

   // The classification systems. An extending table's system replaces the base table's system with the same name
   const std::vector<const ClassificationSystemDeclaration*>& GetClassificationSystems() const { return m_ClassificationSystems; }

   // The classification references of an element role, in table order
   std::vector<const ClassificationDeclaration*> GetClassifications(ElementKind role) const;

   // A property (or quantity) the table exports, and its property set (or quantity set)
   struct ExportedProperty
   {
      const PropertySetDeclaration* pset;
      const PropertyDeclaration* property;
   };

   // The properties and quantities a target is exported to for an element role and owner, in table order
   // (e.g. for the design-value IDS, which checks the values where the export put them)
   std::vector<ExportedProperty> GetExportedProperties(ElementKind role, PropertyOwner attach, std::string_view target) const;

   // A property the table exports for an element role and owner, by property set and property name, or nullopt
   std::optional<ExportedProperty> FindExportedProperty(ElementKind role, PropertyOwner attach, std::string_view pset, std::string_view property) const;

   // The table files, the selected table first
   const std::vector<MappingTableFile>& GetFiles() const { return m_Files; }

   // Writes the table files (path, how it was chosen, name, version) to the log
   void LogFiles() const;

private:
   std::vector<MappingTableFile> m_Files; // the selected table first, then the tables it extends
   std::map<std::string, std::vector<MappingLocation>, std::less<>> m_Locations; // merged import locations
   std::map<ElementKind, ElementSelector> m_Selectors; // merged selectors
   std::vector<const PropertySetDeclaration*> m_PropertySets; // merged export property sets, in table order
   std::vector<const ClassificationSystemDeclaration*> m_ClassificationSystems; // merged
   std::vector<const ClassificationDeclaration*> m_Classifications; // merged, in table order

   void Merge();
};

// Path as text, for messages
std::string PathToString(const std::filesystem::path& path);
