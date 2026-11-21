# Mapping Tables: Design

A mapping table says where PGSuper data lives in an IFC model: which IFC elements play which role, where the import reads each value, and which property sets, quantity sets, and classifications the export writes. The import, the export, the design-value IDS, and the general IDS commands all use the same table, so they agree on where data is.

Designed 2026-09-29 (design note M0, steps M0-M7 of "Mapping tables and IDS" in [IfcImportPlan.md](IfcImportPlan.md)) and built the same week; this note describes the design as it is. Decisions A1-A5 in [IfcImportPlan.md](IfcImportPlan.md) and G1-G7 below apply throughout. The user documentation of the table and its editor is in `Documentation/MappingTableEditor`.

## How it works

1. **Targets.** PGSuper has a fixed list of data items, defined in code (`IfcTargets.cpp`): `girder.fc`, `girder.fci`, `girder.type_names`, `bearing.fixed_x`, `bridge.number_of_spans`, and so on (120 targets). Each has a value kind (stress, length, text, count, ...), the element it belongs to, and, if it's exported, a getter. Tables can't add targets; they only say where targets live in IFC. [Targets.md](Targets.md) lists every target and the PGSuper data behind it.
2. **Tables.** A mapping table is a JSON file (`.ifcmap`, A3). It says:
   - which IFC elements play each element role (girder = `IfcBeam.BEAM`, deck = `IfcSlab.FLOOR`, ...)
   - which property sets, quantity sets, and classifications are written for each role, and which targets their properties hold
   - where the import looks for each target, in order
3. **Import.** For each target, the import tries the table's locations in order and uses the first value found. It converts the value to PGSuper's system units (SI) and logs where it came from.
4. **Export.** For each element role, the export writes the property sets the table declares. A property bound to a target gets the target's value; a constant gets its value; a property with neither is written with no value (it declares the property).
5. **Layering.** The standard table (`MappingTables/Standard.ifcmap`) describes TPF/usBridge, which is what PGSuper exports. An agency table extends it: on import its locations are tried first; on export it can add, replace, or remove property sets.
6. **Selection.** One configured table (File > Configure BridgeLink). An import or export can choose another table, a project remembers its table, and `/IfcMapping=` overrides them (A1, A2, M6).
7. **IDS.** The design-value IDS takes its property locations from the export's table (M4); a table can be written as a general IDS, and generated from one (M5).

Tables contain only names a person wrote (G3). The import doesn't guess where data is: when a value isn't found, it suggests likely properties in the log, to help whoever writes the table.

## Decisions (2026-09-29)

- **G1 Spatial containment:** selectors cover entity, predefined type, classification, and attribute values. Containment rules, such as girders being in the superstructure, stay in code as a prescribed structure (see "The exporter's IFC conventions"). Revisit if models with other structures turn up.
- **G2 One value, several properties:** PGSuper doesn't carry form-stripping, lifting, and release strengths separately, so the export writes f'ci to all three properties. On import, a target is bound only to the property most closely related to it: `girder.fci` is read from `ReleaseStrength`. The other bindings are export only (`"import": false`).
- **G3 No name guessing for data; guesses become hints:** imported data comes only from exact locations in the table. Guessing can find the wrong value without telling anyone, e.g. a "Form Type" property taken as the girder type.
  - Locations use exact property set and property names. The table format has no wildcards or name patterns.
  - The standard table binds only TPF/usBridge locations. Agency tables bind agency names, such as Iowa `IaDOT_PPCB."2_Type"`.
  - When a target isn't found in any location, the import lists properties that look like they might hold it (`CIfcTargetHints`), and never uses them. Example: `Girder type not found. Properties that may hold it: IaDOT_PPCB.2_Type = "BTB45". Add one to the mapping table.` The hint rules: text properties ending in "Type" or containing "Shape", `ObjectType`, and non-usBridge classification names (girder type); properties ending in "Fixity" (bearing fixity); "haunch" in `Name`/`ObjectType` (haunch elements).
- **G4 `NumberOfSpans` conflict:** when `bridge.number_of_spans` doesn't match the modeled substructure (abutment and pier parts − 1), the import logs an error and uses the model.
- **G5 Table files:** the standard table is installed as a file, not embedded. Missing, unreadable, or invalid tables must give users and support staff enough information to fix the problem (see "Loading tables and errors").
- **G6 Units:**
  - PGSuper works in SI internally, and every value is converted to SI as soon as it is read.
  - A table names a unit only for a model value that has no unit of its own: a plain number, or text. Iowa writes f'c as the plain number 5.0, meaning ksi.
  - Unit names come from a short, fixed list mapped to WBFL units.
  - Tables contain no other numbers. Export rounding (e.g. `RoadwayWidth` to 0.1 ft in display units) is part of the target (`display_round`), in code.
- **G7 Export extensions (2026-10-02):** tables still can't add targets. Another plug-in can, through the model builder interface ([ExtensionsDesign.md](ExtensionsDesign.md)): its targets and element roles have its prefix (e.g. `targetApp.`), are export only, and exist only during its own exports. Its table says so with `"requires": ["<extension id>"]`, so any other use of the table refuses it with a clear message.

## Targets

```cpp
struct TargetDef                   // IfcTargets.h
{
   std::string_view name;          // "girder.fci"
   ElementKind element;            // the element role it belongs to
   ValueKind kind;                 // sets the unit type and the IFC value types that can hold it
   std::string_view description;   // for messages and the editor
   ExportUnit display_unit;        // export: the PGSuper display unit of a number
   Float64 display_round;          // export: rounding increment in display units (0: none)
   bool needs_analysis;            // export: left out (no value) when the options leave out analysis results
   ExportGetter get;               // export: the value for an element; empty for import-only targets
   std::string_view provider;      // the export extension that registered it; empty for this project's targets
};
```

- **Value kinds:** `Stress`, `Length`, `Angle`, `Force`, `Area`, `Mass`, `Volume`, `Ratio`, `Count`, `Boolean`, `Text`, `TextList` (every location's value is collected, e.g. candidate girder type names, `CIfcTargetReader::ReadAll`), and, for export extensions, `Time`, `LinearForce`, `MassDensity`, `Real`, `Moment`.
- **Values:** `TargetValue` is a `Float64` (SI), an `Int64`, a `bool`, or a `std::string`.
- **Getters** take an `ExportContext` (the broker, the export options, and the keys of the element being written: segment, pier, strand, bar, bearing, barrier side, diaphragm, ...) and return an `ExportValue`: a value, a list (written as `IfcPropertyListValue`), `NoValue()` (the property is written without a value), or `Absent()` (the property is left out, e.g. `VehicularLiveLoad` without a design live load, or a bar shape dimension the shape doesn't have). Export options that change a value (e.g. `include_camber`) are applied by the getters.
- **There are no setters.** Import code asks for a target by name where it has the element (`CIfcImporter::GetTargetReader().Read("girder.fci", beam)`) and sets the PGSuper data itself. Setters never read template values ("keep import code independent of template values").

## Table format

A table has an IFC side and a PGSuper side:
- `elements`: selectors, i.e. which IFC entities play each element role.
- `property_sets`, `quantity_sets`, `classification_systems`, `classifications`: the export (A5). Any property can be bound to a target.
- `targets`: the import order of each target: an ordered list of locations. A location can be a declared property, or import only: a vendor property set, a stated unit, a text parser.

A target not listed in `targets` is imported from the properties bound to it (except those with `"import": false`), in declaration order.

### Standard table (excerpt)

```json
{
  "format": "PGSuperIfcMapping",
  "version": 1,
  "name": "Standard (TPF/usBridge)",
  "bsdd_dictionaries": {"usBridge": "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/"},
  "bsdd_dictionary": "usBridge",
  "elements": {
    "girder": {"entity": "IfcBeam", "predefined_type": "BEAM"},
    "deck": {"entity": "IfcSlab", "predefined_type": "FLOOR"},
    "bearing": {"entity": "IfcBearing"},
    "tendon_bundle": {"entity": "IfcElementAssembly", "predefined_type": "TENDON_BUNDLE"}
  },
  "property_sets": [
    {
      "name": "Pset_PrecastConcreteElementGeneral", "applies_to": "girder",
      "properties": [
        {"name": "TypeDesignation", "type": "IfcLabel", "target": "girder.family_name"},
        {"name": "CornerChamfer", "type": "IfcPositiveLengthMeasure"},
        {"name": "FormStrippingStrength", "type": "IfcPressureMeasure", "target": "girder.fci", "import": false},
        {"name": "ReleaseStrength", "type": "IfcPressureMeasure", "target": "girder.fci"}
      ]
    },
    {
      "name": "Pset_MaterialConcrete", "applies_to": "girder", "attach": "material", "condition": "always",
      "properties": [
        {"name": "CompressiveStrength", "type": "IfcPressureMeasure", "target": "girder.fc"},
        {"name": "MaxAggregateSize", "type": "IfcPositiveLengthMeasure", "target": "girder.max_aggregate_size"}
      ]
    },
    {
      "name": "usBrPset_PrecastConcreteBeam", "applies_to": "girder",
      "properties": [
        {"name": "ShapeName", "type": "IfcLabel", "target": "girder.type_names", "uri": "bsdd"},
        {"name": "MinimumTimeToDeckPlacement", "type": "IfcTimeMeasure", "uri": "bsdd:MinimumTimetoDeckPlacement"}
      ]
    }
  ],
  "quantity_sets": [
    {
      "name": "Qto_BeamBaseQuantities", "applies_to": "girder", "condition": "quantities", "method": "BaseQuantities",
      "quantities": [{"name": "Length", "type": "IfcQuantityLength", "target": "girder.length"}]
    }
  ],
  "classification_systems": [
    {"name": "usBridge", "source": "usBridge", "edition": "1", "edition_date": "2026-04-20", "specification": "bsdd"}
  ],
  "classifications": [
    {"applies_to": "girder", "system": "usBridge", "identification": "usBridge_GirderPrecastConcrete", "name": "GirderPrecastConcrete", "location": "bsdd"}
  ],
  "targets": {
    "girder.type_names": [
      {"type_attribute": "Name"},
      {"property": {"pset": "usBrPset_PrecastConcreteBeam", "name": "ShapeName"}}
    ],
    "bearing.fixed_x": [
      {"property": {"pset": "Pset_BearingCommon", "name": "DisplacementAccommodated"}, "list_index": 0, "map": {"true": false, "false": true}}
    ]
  }
}
```

### Agency table (Iowa excerpt)

```json
{
  "format": "PGSuperIfcMapping",
  "version": 1,
  "name": "Iowa DOT",
  "extends": "standard",
  "elements": {
    "haunch": {"entity": "IfcSlab", "predefined_type": "USERDEFINED", "attributes": {"ObjectType": "Haunch"}},
    "tendon": {"entity": "IfcTendon", "attributes": {"Name": "Prestressing Strand"}}
  },
  "targets": {
    "girder.fc": [{"property": {"pset": "IaDOT_PPCB", "name": "5_Final Concrete Strength, Fc"}, "unit": "ksi"}],
    "girder.type_names": [{"property": {"pset": "IaDOT_PPCB", "name": "2_Type"}}],
    "bearing.fixed_x": [{"property": {"pset": "Bearing Item - ItemTypes", "name": "3_Fixity"}, "map": {"Fixed": true, "Expansion": false}}],
    "tendon.tension_force": [{"property": {"pset": "IaDOT_Strand", "name": "07-Initial Prestress"}, "parse": {"regex": "([0-9]+(\\.[0-9]+)?)", "group": 1}, "unit": "kip"}],
    "deck.gross_depth": [{"property": {"pset": "IaDOT_Deck", "name": "3_Minimum Thickness"}, "parse": "feet_inches", "unit": "in"}]
  }
}
```

The agency tables are in `Tests/ImportValidation/mappings`: Iowa, PennDOT, and TxDOT (from the names in their example models), and `TestExtension.ifcmap` for the test export extension.

### General rules

- **Unknown keys are errors.** A misspelled key (e.g. `"unti"`) would otherwise be ignored without notice. Any object may have a `"comment"`.
- **Element roles:** fixed by PGSuper (`ElementKind`, `GetElementRoleName`; export extensions add their own, G7).
  - A selector has `entity` (subtypes included), and optionally `predefined_type` (from the type object if there is one; as in IDS, a USERDEFINED element's user-defined type matches too: its type object's `ElementType`, else its `ObjectType`, e.g. `"TENDON_BUNDLE"`), `attributes` (`{ "ObjectType": "Haunch" }`, exact, ignoring case), and `classification` (a reference identification).
  - `any_of` lists alternative selectors.
  - The loader checks entity and attribute names against the IFC schema.
  - `applies_to` (property sets, quantity sets, classifications) can be a list of roles; the loader makes one declaration per role.
- **bSDD dictionaries:** `"bsdd_dictionaries"` (top level) names the bSDD dictionaries the table's URIs use, `{ "usBridge": "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/" }`, each URI ending with `/` (class and property names are appended after `class/` and `prop/`). `"bsdd_dictionary"` names the default. A table has the dictionaries of the tables it extends too, and can add one, replace one with the same name, or name another default; the standard table has `usBridge`, its default, so an agency table's `"bsdd"` URIs are usBridge's. The forms: `"bsdd"` (the default dictionary and the declaration's own name), `"bsdd:Name"` (another name), `"bsdd@dictionary"`, and `"bsdd@dictionary:Name"`. A classification system's `specification` is `"bsdd"` or `"bsdd@dictionary"`; a classification's `location` is its class. A dictionary that isn't declared, or `"bsdd"` without a default, is a table error. The URIs are resolved once the `extends` chain is loaded (`resolve_bsdd_uris`), so no dictionary is in code. Several dictionaries since 2026-11-20 (the user: more than one bSDD dictionary may be used).
- **Target lists:** a target's locations are a list, or `{ "mode": "replace", "locations": [...] }` to replace the base table's locations instead of going before them.
- **Export extensions:** `"requires": ["<extension id>"]` (G7).
- **Not supported:** quantity locations for the import, and table-driven element `Name`/`ObjectType` (they stay in code; see "Export").

### Location reference

| Field | Meaning |
|---|---|
| `property: {pset, name}` | a single, list, or enumerated property, by exact property set and property name |
| `on` | `occurrence` (the default: the occurrence's property sets, then its type's), `type`, or `material` |
| `attribute` / `type_attribute` | an attribute of the occurrence, or of its type object |
| `classification: {system?, identification?}` + `field` | a classification reference and the field to read (`Identification` or `Name`) |
| `value_types` | accepted IFC value types (default: any compatible with the target kind) |
| `unit` | the unit of a model value that has no unit of its own: a plain number or text (G6). A unit in the model always wins. Messages show it in brackets, e.g. `_PS Concrete Beams.Concrete Strength at 28 days (ksi) [ksi]` |
| `parse` | text → value: `number`, `feet_inches` (`8'-6"`, `8.5"`), or `{ "regex": "...", "group": 1 }` |
| `list_index` | an element of a list or enumerated value |
| `map` | maps a model value to the target value. Text matches exactly, ignoring case |

### Units (G6)

The import converts a value to SI (system units) as soon as it's read:

| Model value | Example | Unit used |
|---|---|---|
| measure with its own unit | PGSuper export: `ReleaseStrength` with a ksi unit on the property | the property's unit (`CIfcImportUnits`) |
| measure without its own unit | a pressure measure with no unit on the property | the IFC project unit for that measure |
| plain number or text | Iowa `5_Final Concrete Strength, Fc` = 5.0, PennDOT "8" | the location's `unit`. It's an error if the target needs a unit and the location doesn't give one |

Unit names (`GetTableUnits`):

| Kind | Names |
|---|---|
| Stress | `Pa`, `kPa`, `MPa`, `psi`, `ksi` |
| Length | `mm`, `cm`, `m`, `in`, `ft` |
| Angle | `rad`, `deg` |
| Force | `N`, `kN`, `lbf`, `kip` |

Each name maps to a `WBFL::Units::Measure`. The loader rejects unknown names, and units of the wrong kind for the target (e.g. `in` for `girder.fc`).

### Property declaration reference (export)

| Field | Meaning |
|---|---|
| `name`, `type` | property name and IFC value type (`IfcPressureMeasure`, `IfcLabel`, ...; `IsExportValueType`), or quantity type (`IfcQuantityLength`, ...) |
| `target` | bound target. Without one, the property is written with no value |
| `value` | a constant, written for every element (e.g. `SurfaceFinish = "Raked"`). A property has a target or a value, not both |
| `import` | `false` makes the binding export only (G2) |
| `uri` | `"bsdd"`, `"bsdd:Name"`, `"bsdd@dictionary"`, `"bsdd@dictionary:Name"` (see "bSDD dictionaries"), or an explicit URI, written as the description |
| `enumeration` | name and values of an `IfcPropertyEnumeration`. The value is written as an enumerated value |

Property set fields:
- `name`, `applies_to` (element roles), `attach` (`occurrence`, `type`, or `material`), `uri` (`"bsdd"`: the dictionary, `class/`, and the property set name)
- `condition`: the export option that includes the set: `"classify"` (the default, `options.classify`), `"quantities"` (`options.include_quantities`), or `"always"`
- `shared`: `true` makes one instance related to every element of the role. Only allowed when no property depends on the element; the loader checks this.
- `remove`: `true` leaves the base table's set with the same name, role, and attach out of the export
- quantity sets: `method` (e.g. `BaseQuantities`)

## Layering

- `extends` names the base table: `"standard"` is the installed standard table. A file path is also allowed, relative to the extending table.
- **Selectors:** an agency table adds element roles or replaces a base role's selector.
- **Import:** an agency table's `targets.<t>` list goes **before** the base's list. With `"mode": "replace"` it replaces the base list.
- **Export:** a property set in the agency table replaces the base set with the same name, `applies_to`, and `attach`, in place; `"remove": true` drops it; new sets are added at the end. Classification systems replace by name; classifications follow the property set rules.
- **bSDD dictionaries:** the selected table's dictionaries are added to (and replace by name) those of the tables it extends; each table's URIs are resolved with its own and its base tables' dictionaries.

## Loading tables and errors (G5)

**Where tables come from**, first found wins (`CIfcMappingTable::LoadActive`):
1. `/IfcMapping=<file>` on the command line.
2. The table chosen for one import or export: the IFC Import window, or the Mapping Table tab of the IFC export options. The export tab defaults to the project's table (below). A command line export without `/IfcMapping` uses the project's table.
3. The configured table: the IFC Extension page of File > Configure BridgeLink (registry `IfcExtensions\MappingTable`, A1, M6).
4. The installed standard table: `MappingTables\Standard.ifcmap` in the extension's install folder (next to the DLL).

An `extends` chain is followed from there. `"standard"` always means file 4.

**The project's table.** The project (`CIfcExtensionAgent`, `IIfcProjectMappingTable`) saves the table of the interactive IFC import that made it and the table of its last interactive IFC export, as a path or `standard`. Exports default to the last export table, else the import table, else the configured table. A changed export table marks the project modified. Command line imports and exports don't change it.

**Recent tables.** The tables chosen in these windows and used by interactive imports and exports are kept in one list (registry `RecentMappingTable1..8`), shown by table name.

There is no automatic choice of a table from the IFC model (A2).

**Logged at the start of every import and export:** each table file used (full path, how it was chosen, table name, and version) and the extension version, so support staff can tell from a log which tables produced a result.

**Errors.** Each message says what failed, the full path, how the path was chosen (`MappingTableSource`), and what to do:

| Problem | Message |
|---|---|
| file not found | `File: D:\Tables\IowaDOT.ifcmap` / `The file was chosen by the configured mapping table (File > Configure BridgeLink).` / `The file was not found.` / `What to do: Correct the table, or choose a different table with File > Configure BridgeLink (choose the standard table to go back to it).` The "chosen by" and "What to do" text follows how the table was chosen. |
| standard table missing | the path, and that the IFC extension installation is incomplete |
| can't read the file | the path plus the system error text (access denied, locked, ...) |
| invalid JSON | the path plus line, column, and the parser's message |
| not a mapping table, or an unsupported version | the format and version found, and the version this IFC extension reads |
| `extends` not found, or a cycle | the chain of files |
| validation | every problem, each with its JSON path (e.g. `targets.girder.fc[0].unit`): an unknown target or element role, a unit or IFC type that doesn't fit the target, an unknown unit name, a bad regex, a shared set with element-specific values, a bSDD dictionary that isn't declared, a table that requires an export extension outside that extension's export |

**What happens then:**
- A table error stops the import or export. There is no silent fallback to the standard table, because the results would differ without notice.
- Interactive: a message box shows the error. Headless (`/IfcImport`, `/IfcExport`): the error goes to the log and the command fails.
- The windows that choose a table (`CMappingTableChooser`: the Configure BridgeLink page, the IFC Import window, the export Mapping Table tab) load and validate a table when the user picks it, and show the same messages. A table that can't be used can't be chosen.

## Import

```cpp
class CIfcTargetReader     // IfcTargetReader.h
{
public:
   // the first value found in the target's locations, converted to system units
   std::optional<TargetReading> Read(std::string_view target, IfcSchema::IfcObject object) const;
   // every value found, in location order (TextList targets)
   std::vector<TargetReading> ReadAll(std::string_view target, IfcSchema::IfcObject object) const;
   // the elements that play a role (Select, Matches), and the elements other structures own (ExcludeElements)
   // the log for a target that wasn't found, with hints (ReportNotFound)
};
```

- `CIfcImporter::GetTargetReader()` gives import code the reader of the current import, as `GetUnits()` does for units. The table is loaded before the model is read, so a table problem stops the import before any work is done.
- A `TargetReading` has the value, the location it was found at (and the table that declared it), and the raw value, for the log: "f'ci = 6.8 ksi from IaDOT_PPCB.6_Concrete Release Strength, Fci (Iowa DOT)".
- A value that is found but can't be parsed or converted is logged with its location, and the next location is tried.
- **Hints (G3):** when no location gives a value, `CIfcTargetHints` looks for likely properties on the element and logs them as suggestions, with their values. Hints never set data. A target without a hint rule logs "not found". Only the first element is logged in full; the others are counted and summarized at the end of the import ("... also wasn't found for 19 more elements").
- Geometry wins over properties for layout (D1): a property such as `deck.gross_depth` is a cross check of the measured value, used when the value can't be measured ([ImportGeometry.md](ImportGeometry.md)).

## Export

The exporter decides which elements play which role, and when their property sets are written; the table decides what's written (A5). Element creation, geometry, the spatial structure, and the order elements are created in stay in code ("The exporter's IFC conventions").

- **The session** (`CIfcExportSession`, `IfcPropertyWriter.h/.cpp`) loads the table for one export, with the export profile's hooks ([ExportProfilesPlan.md](ExportProfilesPlan.md)), and keeps what's shared within the export: enumerations, classification systems and references, bar and tendon types, and type names.
- **Writers** (`IfcPropertyWriter.h`): `WritePropertySets` (the role's occurrence property sets and quantity sets, each created once and related to all the objects of the call), `CreateTypePropertySets` (`IfcTypeObject.HasPropertySets`), `CreateSharedPropertySets` (`"shared": true`, related to many elements through `RebarRelationshipBatch`), `WriteMaterialProperties` (`IfcMaterialProperties`, when the element that uses the material creates it), `WriteClassificationSystems`, and `Classify`.
- **Types:** strands and bars with different type property values get different types (`GetTypePropertySignature`), since elements of one type must share their type's values.
- **Units:** each numeric target names its PGSuper display unit (span length, deflection, stress, angle, small and big area, mass, volume, component dimension), and optionally a rounding increment in display units. With `display_units_for_properties` and US display units, the value is converted and rounded, and the property gets its own unit: an `IfcConversionBasedUnit` with a name IDS checkers know ("ksi", "inch", ...), since they convert by the unit's name (`GetPGSuperUnit`, `Units.h`). Otherwise the value is in system units with no property unit, and the project units (SI, `AddUsedProjectUnits`) apply.
- **Constants, enumerations, and no value:** as declared. Each `IfcPropertyEnumeration` is created once and reused.
- **Round trips:** a copy of the source model gets the table's property sets, quantity sets, material properties, and classifications on the elements PGSuper models ([RoundTripDesign.md](RoundTripDesign.md)).
- **What stays in code:** element `Name` and `ObjectType` (tied to IDS pinning, `BeamLabels.h`, and to the import's designation parsing); the bar count quantity of grouped bars (`Qto_ReinforcingElementBaseQuantities.Count`); `Pset_Stationing` on the referents and `Pset_LinearReferencingMethod` on the alignment (`Referents.h`: stationing is part of the alignment structure, as on the import side).

### Checking the export

- `Tests/ImportValidation/compare_export.py` reduces an IFC file to a canonical listing, one line per property: property sets on occurrences, types, and materials, quantity sets, classifications, and `ObjectType`. Elements are keyed by entity, `Name`, and their order among elements with the same entity and name, because GlobalIds change on every export; the property sets of an element are sorted, properties keep their order.
- `export_baseline/<run>.txt.gz` holds the listings of the round-trip exports. `run_validation.py` compares every round-trip export with its baseline ("export same as baseline", or the number of differing lines, with the details in `results/<run>.export-diff.txt`). `--update-export-baseline` records a new baseline after an intended change; every baseline change is reviewed.

## Design-value IDS (M4)

`CIdsExporter` (`IdsExporter.cpp`) loads the export's mapping table (`CIdsExportOptions::mapping_file`) and takes from it the property set, property, and data type of every requirement, and the classification references.
- **Design values by target:** e.g. `girder.fci` gives `ReleaseStrength` (the location the import reads, used in minimal mode) and, with the beam properties, `FormStrippingStrength` and `LiftingStrength`. Quantities are checked as the measure they hold (`IfcQuantityLength` -> `IFCLENGTHMEASURE`). A target the table doesn't export gives no requirement.
- **Constants:** every property with a constant value the table declares for the girder, the girder type, and the strand material, with the table's value; numbers with the IDS tolerance.
- **Pinning:** each specification applies to one element, by `Name` (the default) or `GlobalId` (round trips).
- **Command line:** `/IfcIds=<file.ids>` writes it with an `/IfcExport`, with the default IDS options.
- **Validation:** `run_validation.py` writes the IDS with every round-trip export, compares it with `export_baseline/<run>.ids.txt` (`compare_ids.py`: facets sorted within each specification, the date left out), and checks the export with it ("design-value IDS: 0 failures").

The design-value IDS is model specific: it isn't a general IDS and never a source for mapping tables.

## The exporter's IFC conventions

What a mapping table can't change: how the model is built. These are the exporter's conventions, shared by every export profile (the standard export and every export extension's; [ExportProfilesPlan.md](ExportProfilesPlan.md)), not usBridge's, although usBridge's IDS and IDM shaped several of them. Make one a profile option only when a standard needs something different (decided 2026-11-20).

- **Spatial structure:** `IfcProject` > `IfcSite` > `IfcBridge` (`GIRDER`), which aggregates the `SUPERSTRUCTURE` and `SUBSTRUCTURE` `IfcBridgePart`s. The superstructure aggregates the `DECK` part; the substructure a `PIER` or `ABUTMENT` part for each pier (which aggregates its `FOUNDATION` part) and a `USERDEFINED` part for each PGSplice erection tower. Referents position the bridge, the piers, and the foundations (`IfcRelPositions`).
- **Containment:** the girder segments (`IfcBeam`), the spliced girders' `IfcElementAssembly` `GIRDER`s, and the pier and intermediate diaphragms (`IfcBeam` `DIAPHRAGM`) are in the superstructure; the deck slab and the barriers in the deck part; the bearings and the bearing seat survey points in their pier. A girder-only export has the girder in the project.
- **Types:** each girder has its `IfcBeamType`. A strand's and a bar's data (property sets, material, classification) are on its type (`IfcTendonType`, `IfcReinforcingBarType`), and elements with different type values get different types ("<name> (2)", ...): the types come from the strand or bar path and the type property values (`GetTendonType`, `GetBarType`; see "Strand and bar data on their types" in `CLAUDE.md`).
- **Assemblies:** a segment's strands are an `IfcElementAssembly` with the shared `IfcElementAssemblyType` whose `ElementType` is `TENDON_BUNDLE`; its bars are an `IfcElementAssembly` `REINFORCEMENT_UNIT`. Both are aggregated by the girder.
- **Materials:** concrete materials by strength (`GetConcreteMaterial`; each barrier's from its railing system's strength), a material for each bar type and grade, and one for each strand material, each with the table's material property sets.
- **Not modeled:** the construction sequence (usBridge's IDM keeps it out of the model, `ConstructionSequence.h`).

The table decides the rest: which elements play which element role in an import (selectors), and every property set, quantity set, classification, and constant value of the export.

## Tables and general IDS (M5)

A general IDS says what an export must contain, with no design values. It is not the design-value IDS, which pins values to single elements of one model and is not an input here.

### Table -> general IDS
`/IfcTableToIds=<out.ids> <template.pgt> [/IfcMapping=<table.ifcmap>]` writes the table (the configured table if none is given) as a general IDS, with `WriteTableAsIdsFile` (`IfcTableIds.cpp`):
- **Specifications:** one for each element role that has requirements. It applies to the role's elements through the role's selector (entity, predefined type, and attribute values).
- **Requirements:** properties and quantities with a target are required; constants are required with their value; properties with no value are optional; classifications are required.
- **PGSuper instructions:** a general IDS can't say everything a table does, so the writer puts the rest in the facet's `instructions`, e.g. `PGSuper: target=girder.fci; attach=type; condition=quantities; set=quantities; method=...; pset_uri=...; enumeration=PEnum_...; enumeration_values=A|B`. Specifications get `PGSuper: role=girder`; classifications get `name=`. Other tools ignore the instructions; the generator reads them back.
- **Not written:** material property sets. IDS property facets check an element's own property sets and its type's, not its material's. The log lists each one left out.
- An export extension's table is written the same way with its targets (`IIfcModelBuilder::WriteTableIds`, or `/IfcTableToIds` with `/IfcTestExtension` for the test extension).

### General IDS -> table
`/IfcIdsToTable=<in.ids> <template.pgt> /IfcTable=<out.ifcmap> [/IfcBinding=<binding.json>] [/IfcTableExtends=standard|none] [/IfcLog=<file.log>]` generates a table with `GenerateTableFromIds`:
- **Element role of a specification:** the `role=` instruction, else the binding file, else the standard table's role with the same selector. A specification with no role is reported, with a binding entry to fill in.
- **Target of a property facet:** the `target=` instruction, else the binding file (which can also give the unit, e.g. ksi for a unitless number), else the standard table's property with the same property set and name (e.g. an agency's `ReleaseStrength` gets `girder.fci`). A facet with no target is reported, with a binding entry; it's still declared, so the export writes the property without a value.
- **Classifications:** each one the IDS requires is declared, and its system is added.
- **Targets:** a facet bound with a unit gets a `targets` entry, so the import reads it too.
- **Extending (the default):** the table extends the standard table and has only what differs from it: element roles with another selector, new property sets and quantity sets, standard sets with more properties (the standard set with the new properties added), and new classifications. A property the IDS types differently from the standard table is reported, and the standard table's type is kept.
- **Standalone (`/IfcTableExtends=none`):** the table has everything the IDS has. The round trip check uses it.
- **Checking the result:** the generated table is loaded and validated like any other table, and the log says "The generated table is valid." or gives the table's error.
- **Style:** the table is written like `Standard.ifcmap` (`FormatMappingTable`), with entity and type names in the schema's spelling (`IfcSlab`, `IfcReal`) rather than the IDS's upper case (`IfcSpelling`).

The agency IDS is never changed (A4). Manual assignments live in the binding file, so regenerating the table after the IDS changes keeps them:
```json
{
  "format": "PGSuperIfcBinding", "version": 1,
  "specifications": [
    { "specification": "Prestressed girders", "role": "girder",
      "properties": [ { "pset": "IaDOT_PPCB", "name": "5_Final Concrete Strength, Fc", "target": "girder.fc", "unit": "ksi" } ] }
  ]
}
```
`specification` matches the IDS specification's name or identifier.

### Validation
- **Round trip:** `check_table_ids.py` (the "Standard table <-> IDS" row of `run_validation.py`) writes the standard table as a general IDS, generates a standalone table from it, and compares what the two tables declare for every role: the same 361 declarations (material property sets not compared). Without PGSuper instructions, the roles and targets all come back through the standard table.
- **Binding:** tested by hand with an agency-style IDS and binding file (a girder property set with a bound strength and an unbound fabricator, `ReleaseStrength`, a haunch specification bound by role, an agency classification, and a specification with no role): the table extends the standard table, reports the unbound facet and the specification with no role, leaves out `ReleaseStrength` (the standard table already exports it), and is valid.
- **Not tested yet:** a real agency IDS ([MasterPlan.md](MasterPlan.md)).

## The mapping table editor (M7)

A BridgeLink application plugin in the IFC extension DLL (`IfcExtensions.Manifest.BridgeLink`, `CATID_BridgeLinkAppPlugin`), so it uses the table code as it is and installs with the extension. It edits the table's JSON as an ordered tree, so key order, comments, and what it doesn't show are kept, and writes files in the style of `Standard.ifcmap` (`FormatMappingTable`; `/IfcFormatTable` does the same from the command line). Table > Validate loads the table with `CIfcMappingTable` from memory, as an import or export would. Its user documentation is `Documentation/MappingTableEditor` (Doxygen; the element roles and targets page is generated from `IfcTargets.cpp` by `make_reference.py`; `BuildDocumentation.bat` builds it and copies it to BridgeLink, in the documentation version from `WBFL_VERSION_STRING` in WBFLVersion.h).

**MFC module states** (fixed 2026-10-15): the extension is a regular MFC DLL, and each module state has its own window handle map. The frame and views are created by EAF in the application's module state, the panes in the DLL's (where their dialog templates are), so views are used in the application's state (`AfxGetAppModuleState`) and panes in the DLL's; a pane gets its parent with `CWnd::FromHandle` in the DLL's state, because `CWnd::CreateDlgIndirect` validates the parent in the current state.

**Not published yet** ([MasterPlan.md](MasterPlan.md)): it goes to the online documentation root (BridgeLink's default `https://www.wsdot.wa.gov/eesc/bridge/software/Documentation/`, or `Settings\DocumentationRoot`) under `IfcMappingTableEditor/<major.minor>/`, with a copy of the map named `IfcMappingTableEditor.dm.html` (EAF retries with `.html`); until then Help says the topic isn't found. Add `Docs\IfcMappingTableEditor\<version>\` to `UpdateBinFiles.bat` and the installer for local documentation.

## Code layout

| File | Contents |
|---|---|
| `IfcTargets.h/.cpp` | `ExportContext`, `TargetDef`, the target list, the element role names (`ElementKind`, `ValueKind`, `TargetValue`, `ExportValue` are in the extension SDK, `Include\IfcExport\IfcExtensionTypes.h`) |
| `IfcMappingTable.h/.cpp` | finding table files, JSON load (nlohmann-json, from vcpkg), the `extends` merge, bSDD dictionaries, validation, error messages |
| `IfcTargetReader.h/.cpp` | the import engine |
| `IfcTargetHints.h/.cpp` | hint rules for targets that weren't found (G3) |
| `IfcPropertyWriter.h/.cpp` | the export session and the property set, classification, and type writers |
| `ExportProfile.h/.cpp` | `CExportProfile`: the standard export or an export extension's ([ExportProfilesPlan.md](ExportProfilesPlan.md)) |
| `IfcExtensionRegistry.h/.cpp`, `Include\IfcExport` | export extensions' targets and roles, and the SDK ([ExtensionsDesign.md](ExtensionsDesign.md)) |
| `TestExportExtension.h/.cpp` | the test export extension (`/IfcTestExtension`, `check_extension.py`) |
| `Classifications.h` | the classification helpers (`AssociateClassification`, `HasClassification`) |
| `IdsBuilder.h` | IDS document helpers (facets, restrictions, writing) shared by the design-value IDS and M5 |
| `IdsExporter.h/.cpp` | the design-value IDS (M4) |
| `IfcTableIds.h/.cpp` | table -> general IDS writer and general IDS -> table generator (M5) |
| `IfcTableFormat.h/.cpp` | `FormatMappingTable`: writes a table in the style of `Standard.ifcmap` (the editor, `/IfcFormatTable`) |
| `MappingTableChooser.h/.cpp` | `CMappingTableChooser`: the controls that choose a table (default table, recent tables, browse, validation) (M6) |
| `MappingTableWindows.h/.cpp` | the Configure BridgeLink page and its callback, the IFC Import window, and the export Mapping Table tab (M6); the configured table itself is `CIfcMappingTable::GetTableSetting`/`SetTableSetting` |
| `IfcProjectMappingTable.h/.cpp` | `IIfcProjectMappingTable`, the tables a project remembers, and the export default (M6) |
| `MappingEditorDoc.h/.cpp` | mapping table editor (M7): the BridgeLink application plugin, doc template, and document (the table's JSON) |
| `MappingEditorViews.h/.cpp` | editor frame, section tree, and the view that hosts the panes |
| `MappingEditorPanes.h/.cpp` | editor panes (table, selector, target locations) and the location dialog |
| `MappingEditorExportPanes.h/.cpp` | editor export sections: property and quantity sets, classification systems, classifications, and the base table's entries |
| `MappingEditorUtil.h` | helpers shared by the editor panes |
| `MappingEditorModelPane.h/.cpp` | editor Model item: browse a model, see what the table reads, add a property as a location |
| `IfcModelBrowser.h/.cpp` | `CIfcModel`: an IFC model opened by the editor, its elements and properties, and "try the table on the model" (with `CIfcTargetReader` and `CIfcTargetHints`) |
| `IfcRoundTrip.h/.cpp`, `IfcRoundTripMerge.h/.cpp` | the round trip registry and the merge into a copy of the source model ([RoundTripDesign.md](RoundTripDesign.md)) |
| `IfcExtensions.Manifest.BridgeLink` | registers the editor with BridgeLink (copied next to the DLL by the post-build step) |
| `MappingTables/Standard.ifcmap` | the standard table, installed next to the DLL (post-build copy for development builds) |
| `Tests/ImportValidation/mappings/*.ifcmap` | the Iowa, PennDOT, TxDOT, and test extension tables; `models.json` names a model's `mapping` |
