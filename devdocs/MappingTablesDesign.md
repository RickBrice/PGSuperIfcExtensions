# Mapping Tables: Design Note (M0)

Drafted 2026-09-29. This is step M0 of "Mapping tables and IDS" in [IfcImportPlan.md](IfcImportPlan.md). It covers the element registry, the target vocabulary, the JSON table format, and how the import and export engines use a table. Decisions A1–A5 in the plan and G1–G6 below apply throughout.

Status: **reviewed 2026-09-29.** The questions of the first draft are answered (G1–G6).

## How it works

1. **Targets.** PGSuper has a fixed list of data items, defined in code: `girder.fc`, `girder.fci`, `girder.type_names`, `bearing.fixed_x`, `bridge.number_of_spans`, and so on. For each one, the code knows:
   - its value kind (stress, length, text, count, ...)
   - how to get it from PGSuper (export)
   - where to put it in `CBridgeDescription2` (import)

   Tables can't add targets. They only say where targets live in IFC.
2. **Tables.** A mapping table is a JSON file (A3). It says:
   - which IFC elements play which role (girder = `IfcBeam.BEAM`, deck = `IfcSlab.FLOOR`, ...)
   - which property sets are written for each role, and which targets their properties hold
   - where the importer looks for each target, in order
3. **Import.** For each target, the importer tries the table's locations in order and uses the first value found. It converts the value to PGSuper's internal SI units and logs where the value came from.
4. **Export.** For each element role, the exporter writes the property sets the table declares. A property tied to a target gets the target's value. A property with no target is written with no value, as some are today.
5. **Layering.** The standard table describes TPF/usBridge, which is what PGSuper exports today. An agency table extends it:
   - On import, the agency's locations are tried first.
   - On export, the agency table can add, replace, or remove property sets.
6. **Selection.** The active table comes from the setting on the configuration wizard page, or from `/IfcMapping=` (A1, A2).
7. **IDS.** The design-value IDS takes its property locations from the same table (M4).

Until the mapping table editor exists (requirement R3 in the plan, step M7), tables are written by hand.

Tables contain only names the user specifies (G3). The importer doesn't guess where data is. When a value isn't found, it only suggests likely properties in the log, to help whoever writes the table.

## Decisions (2026-09-29)

- **G1 Spatial containment:** selectors cover entity, predefined type, classification, and attribute or property values. Containment rules, such as girders being in the superstructure, stay in code as a prescribed structure. Revisit if models with other structures turn up.
- **G2 One value, several properties:** PGSuper doesn't carry form-stripping, lifting, and release strengths separately, so the export writes f'ci to all three properties. On import, a target is bound only to the property most closely related to it: `girder.fci` is read from `ReleaseStrength`. The other bindings are export only (`"import": false`).
- **G3 No name guessing for data; guesses become hints (option D):** imported data comes only from exact locations in the table. Guessing can find the wrong value without telling anyone, e.g. a "Form Type" property taken as the girder type.
  - Locations use exact property set and property names. The table format has no wildcards or name patterns.
  - The standard table binds only TPF/usBridge locations. Agency tables bind agency names, such as Iowa `IaDOT_PPCB."2_Type"`.
  - Today's guesses are kept, but only for **hints in the import log**. When a target isn't found in any table location, the importer lists properties that look like they might hold it, and never uses them. Example: `Girder type not found. Properties that may hold it: IaDOT_PPCB.2_Type = "BTB45". Add one to the mapping table.`
  - The guesses used for hints are: text properties ending in "Type" or containing "Shape", `ObjectType`, and non-usBridge classification names (girder type); properties ending in "Fixity" (bearing fixity); "haunch" in `Name`/`ObjectType` (haunch elements).
  - The hints go into the structured import report and BCF later (R2), so the table author and the modeler both see them.
- **G4 `NumberOfSpans` conflict:** when `usBrPset_BridgeGeometry.NumberOfSpans` doesn't match the modeled substructure (ABUTMENT + PIER bridge parts − 1), the importer logs an error and uses the model. Today it throws (`IfcBridgeImporter.cpp:74`). The error goes into the structured import report and BCF when those exist (R2).
- **G5 Table files:** the standard table is installed as a file, not embedded. Missing, unreadable, or invalid tables must give users and support staff enough information to fix the problem (see "Loading tables and errors").
- **G6 Units:**
  - PGSuper works in SI internally, and every value is converted to SI as soon as it is read.
  - A table names a unit only for a model value that has no unit of its own: a plain number, or text. Iowa writes f'c as the plain number 5.0, meaning ksi.
  - Unit names come from a short, fixed list mapped to WBFL units (option (a)).
  - Tables contain no other numbers. Export rounding (e.g. `RoadwayWidth` to 0.1 ft in display units) stays in code.

## Where the data is today

### Import
Every pset and property the importer reads today, and what happens to it:

| Target | Where it's read | Code | With tables |
|---|---|---|---|
| `bridge.number_of_spans` | `IfcBridge` `usBrPset_BridgeGeometry.NumberOfSpans` (IfcInteger); cross check only | `IfcBridgeImporter.cpp:69` | standard table; a conflict is logged as an error (G4) |
| `girder.fci` | `IfcBeam` `Pset_PrecastConcreteElementGeneral.ReleaseStrength` | `IfcBridgeImporter.cpp:300` | standard table |
| `girder.fc` | beam material `Pset_MaterialConcrete.CompressiveStrength` | `IfcBridgeImporter.cpp:313` | standard table |
| `girder.type_names` | `IfcBeamType.Name`, `ObjectType`, text properties named "...Type" or "...Shape...", non-usBridge classification names | `GetGirderTypeNames`, `IfcBridgeImporter.cpp:730` | standard table: `IfcBeamType.Name`, `usBrPset_PrecastConcreteBeam.ShapeName`. Agency names go in agency tables; the guesses become log hints (G3) |
| `bearing.fixed_x`, `bearing.fixed_y` | `Pset_BearingCommon.DisplacementAccommodated` (boolean list), else any text property named "...Fixity" | `get_bearing_fixity`, `Piers.cpp:307` | standard table: `DisplacementAccommodated`. Iowa and PennDOT "Fixity" go in agency tables; "...Fixity" becomes a log hint |
| `girder.assembly_place`, `girder.casting_method` | `Pset_ConcreteElementGeneral.AssemblyPlace`/`CastingMethod`, or classification `GirderPrestressedConcrete` (`HasValidGirdersByTPF`, not called today) | `IfcBridgeImporter.cpp:386` | standard table |
| `girder.designation` | `IfcBeam` `Pset_PrecastConcreteElementGeneral.DesignLocationNumber`, else parsed from `Name` | `get_girder_key`, `Utilities.h` | standard table; the `Name` parsing stays in code. (Found while implementing M1; missed in the first draft) |
| haunch selector | `IfcSlab`/`IfcBuildingElementPart` with "haunch" in `ObjectType` or `Name` | `DeckSlab.cpp:269` | agency tables (Iowa: `IfcSlab.USERDEFINED` with ObjectType "Haunch"); "haunch" in a name becomes a log hint. Haunches in the deck solid, or parts aggregated under the deck, stay in code |

These stay in code:
- `Pset_Stationing.Station`/`IncomingStation` on referents (`Piers.cpp:83`, `IfcAlignmentImporter.cpp:633`). IFC 4.3 defines stationing this way; it isn't agency data.
- Girder designations parsed from `Name` in `get_girder_layout` (checked against location; location wins).

### Export
- **Property sets:** about 30 `Create_*` functions in `PropertySets.h`. `Materials.h` adds `Pset_MaterialConcrete` (`GetConcreteMaterial`) and `Pset_MaterialSteel`. `QuantitySets.h` has four `Qto_*`, and `USBridge_Classifications.h` has 17 `Classify_usBridge_*`. They're attached in `IfcExporter.cpp`: beams at 3239–3258, slab 2478–2485, piers 2724, bridge parts 3363–3392, barriers 3456–3464, bridge 3685–3694, project 3893–3914, rebar and strands 370–1560.
- **Features a table must be able to express:**
  - Properties with no value (`IfcValue{}`) that are written only to declare the pset's structure, e.g. `CornerChamfer`, `PieceMark`.
  - bSDD URIs, as property descriptions (`BSDD_PROPERTY(...)`) and pset descriptions (e.g. `usBrPset_Roadway`).
  - Enumerated values (`usBrPEnum_ElementStatus`, `usBrPEnum_SubstructureType`) and list values.
  - Material psets (`IfcMaterialProperties`), type psets (`HasPropertySets`), and occurrence psets (`IfcRelDefinesByProperties`).
  - Psets with no element-specific values, created once and related to many elements (rebar and tendons, `IfcExporter.cpp:405-416, 658-675`). This matters for file size.
  - One PGSuper value written to several properties (G2).
- **Handled in code, not in the table:**
  - export option gates (`CamberAtMidspan` is empty unless `options.include_camber`)
  - units: display units with an `IfcConversionBasedUnit` per property when `display_units_for_properties` and the unit mode is US, otherwise SI with no property unit
  - rounding (G6)
- **IDS:** `IdsExporter.cpp:543-754` repeats about 45 of these pset/property/data-type triples.

## Element registry

One identity per PGSuper item that becomes an IFC element. Tables use it (`applies_to`), and so will GlobalId persistence (R1).

**Status (M1):** the import engine doesn't need the registry, so it's built with R1 and the export engine (M3). M1 has only `ElementKind`, in `IfcTargets.h`, with a `Haunch` role added.

```cpp
enum class ElementKind { Project, Site, Bridge, BridgePart, Pier, Foundation, Alignment, Referent,
                         Girder, Deck, Bearing, Barrier };

struct ElementId
{
   ElementKind kind;
   BridgePartRole role = BridgePartRole::None; // Superstructure, Substructure, Deck (BridgePart only)
   PierIDType pier = INVALID_ID;              // Pier, Foundation, Referent, Bearing
   SegmentIDType segment = INVALID_ID;        // Girder (one IfcBeam per segment)
   pgsTypes::PierFaceType face = pgsTypes::Back;         // Bearing
   IndexType index = 0;                       // Bearing number at a girder end, Barrier side (0 left, 1 right)
   auto operator<=>(const ElementId&) const = default;
};
```

- **Keys are PGSuper IDs, not indices.** Inserting or removing a span doesn't shift other elements (R1).
- **Registry:** `CIfcElementRegistry` holds `ElementId ↔ IFC instance` for the current import or export, and `ElementId → GlobalId` for persistence.
  - **Import:** filled as elements are identified: bridge parts by predefined type, girders by `get_girder_layout` (beam id → segment ID once the bridge description exists), bearings by `locate_bearings`, deck by `FindDeck`.
  - **Export:** filled as elements are created. The GlobalId provider (R1) takes an `ElementId`.
- **Outside the registry:** materials, psets, relationships, types, tendons, and rebar. They get new GlobalIds on every export (B2). Material psets are reached through their element ("girder material").
- **Persistence** (R1 work, not M1): `CIfcExtensionAgent::Save/Load` gets an `IfcElementRegistry` unit with the `ElementId → GlobalId` map. The agent unit version goes to 2.0; loading version 1.0 gives an empty map.

## Targets

A target is a PGSuper data item with a fixed name, the element it belongs to, a value kind, and a getter and/or setter. Targets are defined in code (`IfcTargets.cpp`). Tables refer to them by name.

```cpp
enum class ValueKind { Stress, Length, Angle, Ratio, Count, Boolean, Text, TextList };

using TargetValue = std::variant<Float64, Int64, bool, std::string, std::vector<std::string>>; // Float64 in SI (system units)

struct TargetDef
{
   std::string_view name;        // "girder.fci"
   ElementKind element;
   ValueKind kind;               // sets the WBFL unit type and the default IFC measure type
   bool collect_all = false;     // import gathers every location's value instead of the first (TextList)
   std::function<std::optional<TargetValue>(const ExportContext&, const ElementId&)> get;  // empty: import only
   std::function<void(ImportContext&, const ElementId&, const TargetValue&)> set;          // empty: export only
};
```

- The getter returns `nullopt` for "no value". The property is still written with no value, as the export does today.
- The getter applies export options (`include_camber`, ...), reading `CIfcExportOptions` from the context. Tables don't know about export options.
- Setters write into `CBridgeDescription2` through the `ImportContext`. They never read template values ("keep import code independent of template values").

**Status (M1):**
- `TargetDef` has only the name, element, value kind, and a description for messages. Getters and setters come with the export engine (M3), when there are enough targets to justify them. Import code sets the values it reads.
- `TargetValue` has no list alternative. `ReadAll` returns one reading per value instead.
- `ValueKind` also has `Force`.

### Vocabulary for M1–M2 (import)

| Target | Element | Kind | Import setter | Export getter |
|---|---|---|---|---|
| `bridge.number_of_spans` | Bridge | Count | cross check against the modeled substructure; on a mismatch, log an error and use the model (G4) | `IBridge::GetSpanCount` |
| `girder.fc` | Girder | Stress | `Segment.Material.Concrete.Fc` | `IMaterials::GetSegmentFc28` |
| `girder.fci` | Girder | Stress | `Segment.Material.Concrete.Fci` | `GetSegmentFc` at release |
| `girder.designation` | Girder | Text | girder key cross check in `get_girder_layout` (location wins) | the girder label |
| `girder.type_names` | Girder | TextList (collect all) | candidate names for `GetGirderLibraryEntry` | `GetGirderName` |
| `girder.assembly_place`, `girder.casting_method` | Girder | Text | precast check | "FACTORY", "PRECAST" |
| `bearing.fixed_x`, `bearing.fixed_y` | Bearing | Boolean | `CBearingData2::FixedX/FixedY` | from support fixity (Phase 2 boundary conditions) |
| `deck.gross_depth` | Deck | Length | cross check against the measured depth (D1: geometry wins, a conflict over 1/8 in is logged); used when the depth can't be measured. Agency tables only: the standard location `Qto_SlabBaseQuantities.Depth` is exported as a 0.0 placeholder today (fix in M3) | `GetGrossSlabDepth` |

M3 adds the export-only targets behind today's psets: `girder.fc_lifting`, `girder.fc_hauling`, `girder.jacking_stress`, `girder.camber_at_release`, `girder.camber_after_losses`, `girder.screed_camber`, `girder.camber_ratio`, `girder.batter`, `girder.span`, `girder.slope`, `girder.roll`, `girder.bunk_point`, `girder.family_name`, `girder.design_location`, `bridge.length`, `bridge.roadway_width`, `bridge.start_station`, `bridge.end_station`, `bridge.max_skew`, `pier.substructure_type`, `material.max_aggregate_size`, strand and rebar material targets, and so on. The full list comes from walking `PropertySets.h` in M3.

## Table format (JSON, A3)

A table has an IFC side and a PGSuper side:
- `elements`: selectors, i.e. which IFC entities play each element role.
- `property_sets`, `quantity_sets`, `classifications`, `attributes`: the export structure (A5). Any property can be bound to a target.
- `targets`: import order per target: an ordered list of locations. A location can refer to a declared property, or be import only: a vendor pset, a stated unit, a text parser.

A target not listed in `targets` is imported from its bound properties (except those with `"import": false`), in declaration order.

### Standard table (excerpt)

```json
{
  "format": "PGSuperIfcMapping",
  "version": 1,
  "name": "Standard (TPF/usBridge)",
  "elements": {
    "girder": { "entity": "IfcBeam", "predefined_type": "BEAM" },
    "deck": { "entity": "IfcSlab", "predefined_type": "FLOOR" },
    "bearing": { "entity": "IfcBearing" }
  },
  "property_sets": [
    {
      "name": "Pset_PrecastConcreteElementGeneral",
      "applies_to": "girder",
      "attach": "occurrence",
      "properties": [
        { "name": "TypeDesignation", "type": "IFCLABEL", "target": "girder.family_name" },
        { "name": "CornerChamfer", "type": "IFCPOSITIVELENGTHMEASURE" },
        { "name": "FormStrippingStrength", "type": "IFCPRESSUREMEASURE", "target": "girder.fci", "import": false },
        { "name": "LiftingStrength", "type": "IFCPRESSUREMEASURE", "target": "girder.fci", "import": false },
        { "name": "ReleaseStrength", "type": "IFCPRESSUREMEASURE", "target": "girder.fci" },
        { "name": "CamberAtMidspan", "type": "IFCRATIOMEASURE", "target": "girder.camber_ratio" }
      ]
    },
    {
      "name": "Pset_MaterialConcrete",
      "applies_to": "girder",
      "attach": "material",
      "properties": [
        { "name": "CompressiveStrength", "type": "IFCPRESSUREMEASURE", "target": "girder.fc" },
        { "name": "MaxAggregateSize", "type": "IFCPOSITIVELENGTHMEASURE", "target": "material.max_aggregate_size" }
      ]
    },
    {
      "name": "usBrPset_PrecastConcreteBeam",
      "applies_to": "girder",
      "properties": [
        { "name": "ShapeName", "type": "IFCLABEL", "target": "girder.type_names", "uri": "bsdd" }
      ]
    },
    {
      "name": "usBrPset_BridgeGeometry",
      "applies_to": "bridge",
      "properties": [
        { "name": "RoadwayWidth", "type": "IFCPOSITIVELENGTHMEASURE", "target": "bridge.roadway_width", "uri": "bsdd" },
        { "name": "NumberOfSpans", "type": "IFCINTEGER", "target": "bridge.number_of_spans", "uri": "bsdd" }
      ]
    }
  ],
  "classifications": [
    { "applies_to": "girder", "system": "usBridge", "identification": "usBridge_GirderPrecastConcrete" }
  ],
  "attributes": [
    { "applies_to": "girder", "attribute": "Name", "value": "{girder_label}" }
  ],
  "targets": {
    "girder.type_names": [
      { "type_attribute": "Name" },
      { "property": { "pset": "usBrPset_PrecastConcreteBeam", "name": "ShapeName" } }
    ],
    "bearing.fixed_x": [
      { "property": { "pset": "Pset_BearingCommon", "name": "DisplacementAccommodated" }, "list_index": 0,
        "map": { "true": false, "false": true } }
    ],
    "bearing.fixed_y": [
      { "property": { "pset": "Pset_BearingCommon", "name": "DisplacementAccommodated" }, "list_index": 1,
        "map": { "true": false, "false": true } }
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
    "haunch": { "entity": "IfcSlab", "predefined_type": "USERDEFINED", "attributes": { "ObjectType": "Haunch" } }
  },
  "targets": {
    "girder.fc": [ { "property": { "pset": "IaDOT_PPCB", "name": "5_Final Concrete Strength, Fc" }, "unit": "ksi" } ],
    "girder.fci": [ { "property": { "pset": "IaDOT_PPCB", "name": "6_Concrete Release Strength, Fci" }, "unit": "ksi" } ],
    "girder.type_names": [ { "property": { "pset": "IaDOT_PPCB", "name": "2_Type" } } ],
    "bearing.fixed_x": [ { "property": { "pset": "Bearing Item - ItemTypes", "name": "3_Fixity" },
                           "map": { "Fixed": true, "Expansion": false } } ],
    "bearing.fixed_y": [ { "property": { "pset": "Bearing Item - ItemTypes", "name": "3_Fixity" },
                           "map": { "Fixed": true, "Expansion": false } } ],
    "deck.gross_depth": [ { "property": { "pset": "IaDOT_Deck", "name": "3_Minimum Thickness" },
                            "parse": "feet_inches", "unit": "in" } ]
  }
}
```

PennDOT is the same with:
- `_PS Concrete Beams."Concrete Strength at 28 days (ksi)"` (text, `"parse": "number"`, `"unit": "ksi"`)
- `_PS Concrete Beams.Type` and `"Beam Designation"`
- `_Bearing Assembly.Fixity`
- `_Deck."Minimum Thickness (in)"`

The values in `map` (e.g. every value `3_Fixity` actually takes) are taken from the models when the tables are written in M1–M2.

### General rules

- **Unknown keys are errors.** A misspelled key (e.g. `"unti"`) would otherwise be ignored without notice. Any object may have a `"comment"`.
- **Element roles:**
  - A selector has `entity` (subtypes included), and optionally `predefined_type` (from the type object if there is one), `attributes` (`{ "ObjectType": "Haunch" }`, exact, ignoring case), and `classification` (a reference identification).
  - `any_of` lists alternative selectors.
  - The loader checks entity and attribute names against the IFC schema.
  - Roles the importer uses in M1: `girder` (together with superstructure containment, G1), `deck`, `haunch`, `bearing`.
- **Target lists:** a target's locations are a list, or `{ "mode": "replace", "locations": [...] }` to replace the base table's locations instead of going before them.
- **Not in M1:** quantity locations (`quantity: {qto, name}`), and the `quantity_sets`, `classifications`, and `attributes` export sections. They come with the export engine (M3).

### Location reference

| Field | Meaning |
|---|---|
| `property: {pset, name}` | a single, list, or enumerated property, by exact pset and property name |
| `on` | `occurrence` (default: occurrence psets first, then type psets, as `GetPropertySet` does today), `type`, or `material` |
| `attribute` / `type_attribute` | an attribute of the occurrence, or of its type object |
| `classification: {system?, identification?}` + `field` | a classification reference and the field to read (`Identification` or `Name`) |
| `quantity: {qto, name}` | an `IfcElementQuantity` quantity |
| `value_types` | accepted IFC value types (default: any compatible with the target kind) |
| `unit` | the unit of a model value that has no unit of its own: a plain number or text (G6). A unit in the model always wins. Messages show it in brackets, e.g. `_PS Concrete Beams.Concrete Strength at 28 days (ksi) [ksi]` |
| `parse` | text → value: `number`, `feet_inches` (`8'-6"`, `8.5"`), or `{ "regex": "...", "group": 1 }` |
| `list_index` | an element of a list value |
| `map` | maps a model value to the target value. Text matches exactly, ignoring case |
| `import` | `false` on a declared property that the importer must not read (default `true`) |

### Units (G6)

The importer converts a value to SI (system units) as soon as it's read:

| Model value | Example | Unit used |
|---|---|---|
| measure with its own unit | PGSuper export: `ReleaseStrength` with a ksi unit on the property | the property's unit (`CIfcImportUnits`, as today) |
| measure without its own unit | a pressure measure with no unit on the property | the IFC project unit for that measure (as today) |
| plain number or text | Iowa `5_Final Concrete Strength, Fc` = 5.0, PennDOT "8" | the location's `unit`. It's an error if the target needs a unit and the location doesn't give one |

Unit names:

| Kind | Names |
|---|---|
| Stress | `Pa`, `kPa`, `MPa`, `psi`, `ksi` |
| Length | `mm`, `cm`, `m`, `in`, `ft` |
| Angle | `rad`, `deg` |
| Force | `N`, `kN`, `lbf`, `kip` |

Each name maps to a `WBFL::Units::Measure`. The loader rejects unknown names, and units of the wrong kind for the target (e.g. `in` for `girder.fc`). The list grows with the targets.

### Property declaration reference (export)

| Field | Meaning |
|---|---|
| `name`, `type` | property name and IFC value type (`IFCPRESSUREMEASURE`, `IFCLABEL`, ...) |
| `target` | bound target. Without one, the property is written with no value |
| `import` | `false` makes the binding export only (G2) |
| `uri` | `"bsdd"` (built from the pset and property name, as `BSDD_PROPERTY` does) or an explicit URI, written as the description |
| `enumeration` | name and values of an `IfcPropertyEnumeration`. The value is written as an enumerated value |

Pset fields:
- `name`, `applies_to` (element role), `attach` (`occurrence`, `type`, or `material`), `uri`
- `shared`: `true` makes one instance related to all elements of the role. Only allowed when no property depends on the element; the loader checks this.

## Layering

- `extends` names the base table: `"standard"` is the installed standard table. A file path is also allowed, relative to the extending table.
- **Selectors:** an agency table adds element roles or replaces a base role.
- **Import:** an agency table's `targets.<t>` list goes **before** the base's list. With `"mode": "replace"` it replaces the base list.
- **Export:** a pset in the agency table replaces the base pset with the same name and `applies_to`. `"remove": true` drops it. New psets are added. Classifications and attributes follow the same rules.

## Loading tables and errors (G5)

**Where tables come from**, first found wins:
1. `/IfcMapping=<file>` on the command line.
2. The configuration wizard setting (registry, A1). Added in M6; until then, only the command line.
3. The installed standard table: `MappingTables\Standard.json` in the extension's install folder (next to the DLL).

An `extends` chain is followed from there. `"standard"` always means file 3.

**Logged at the start of every import and export:**
- each table file used: full path, how the path was chosen ("command line", "configuration setting", "installed standard table"), table name, and version
- the extension version

Support staff can then tell from a log which tables produced a result.

**Errors.** Each message says what failed, the full path, how the path was chosen, and what to do:

| Problem | Example message |
|---|---|
| file not found | `IFC mapping table not found: D:\Tables\IowaDOT.json (from the configuration setting). Choose a different file in the BridgeLink configuration, or remove the setting to use the standard table.` |
| standard table missing | `The standard IFC mapping table is missing: C:\...\MappingTables\Standard.json. The IFC extension installation is incomplete; reinstall it.` |
| can't read the file | the path plus the system error text (access denied, locked, ...) |
| invalid JSON | the path plus line, column, and the parser's message |
| not a mapping table, or an unsupported version | `... has format "X" version 3; this version of the IFC extension reads PGSuperIfcMapping version 1.` |
| `extends` not found, or a cycle | the chain of files |
| validation | every problem, each with its JSON path (e.g. `targets.girder.fc[0].unit`): an unknown target or element role, a unit or IFC type that doesn't fit the target, an unknown unit name, a bad regex, a shared pset with element-specific values |

**What happens then:**
- A table error stops the import or export. There is no silent fallback to the standard table, because the results would differ without notice.
- Interactive: a message box shows the error and the log path. Headless (`/IfcImport`, `/IfcExport`): the error goes to the log and the command fails.
- Findings that aren't errors (e.g. a declared property with no target) are logged as warnings and don't stop anything.
- The configuration page (M6) loads and validates a file when the user picks it, and shows the same messages.

## Import engine (M1)

```cpp
class CIfcMappingTable;   // loaded, merged, validated table
class CIfcTargetReader
{
public:
   CIfcTargetReader(const CIfcMappingTable& table, const CIfcImportUnits& units);
   // first value found in the target's location list, converted to system units
   std::optional<TargetReading> Read(std::string_view target, IfcSchema::IfcObject object) const;
   // every value found, in location order (collect_all targets)
   std::vector<TargetReading> ReadAll(std::string_view target, IfcSchema::IfcObject object) const;
};

struct TargetReading
{
   TargetValue value;
   const MappingLocation* location; // where it was found, and in which table
   std::string raw;                 // "6.8" as found, for the report
};
```

As built in M1, the reader also:
- selects the elements that play a role (`Select`, `Matches`)
- reports targets that weren't found (`ReportNotFound`). Only the first element is logged in full, with hints; the others are counted and summarized at the end of the import ("... also wasn't found for 19 more elements"), so the log stays readable.

`CIfcImporter::GetTargetReader()` gives import code the reader for the current import, as `GetUnits()` does for units. The table is loaded before the model is read, so a table problem stops the import before any work is done.

- Import code calls `Read("girder.fci", beam)` in place of the hard-coded `GetMeasureProperty(..., "Pset_PrecastConcreteElementGeneral", "ReleaseStrength")`. Where the value is set stays in import code for M1. Setters in `TargetDef` come in when there are enough targets to justify it.
- A value that is found but can't be parsed or converted is logged with its location, and the next location is tried.
- **Hints (G3):** when no location gives a value, `CIfcTargetHints` looks for likely properties on the element and logs them as suggestions, with their values. Hints never set data. Each target can have a hint rule in code (girder type, bearing fixity, haunch elements to start with); targets without one just log "not found".
- Each reading goes to the import log now ("f'ci = 6.8 ksi from IaDOT_PPCB.6_Concrete Release Strength, Fci (Iowa DOT)"), and to the structured import report and BCF later (R2, Phase 0.5).
- Performance: M1 walks `IsDefinedBy` for each lookup, as the importer did before. Indexing psets once per object waits until there are enough targets for it to matter.

**M1 acceptance:**
- The standard table reproduces today's reads of standard locations.
- The guessing in `GetGirderTypeNames`, the "...Fixity" fallback in `get_bearing_fixity`, and the haunch name check move to `CIfcTargetHints` and no longer set data (G3).
- Without the guesses, the Iowa and PennDOT girder type names and bearing fixity would be lost. So M1 includes first Iowa and PennDOT tables with the girder type, fixity, and haunch locations, and `run_validation.py` imports those models with them. With those tables, the results are unchanged for all models (0 mismatches where there are 0 today).
- Without those tables, the Iowa and PennDOT import logs show hints naming `IaDOT_PPCB.2_Type`, `3_Fixity`, `_PS Concrete Beams.Type`, and so on. `run_validation.py` checks this with one extra run of each model using only the standard table.
- The import log names the table file and the source of each value.
- A `NumberOfSpans` mismatch is logged as an error and the model is used (G4).
- Table loading errors were checked by hand in M1 with broken tables: a missing file, invalid JSON, a table with 7 problems (misspelled key, unknown target, unknown unit, unit on a boolean, bad map value, unknown entity, unknown attribute), an unsupported version, and an `extends` cycle. They become automated tests in Phase 8.

## Export engine (M3, outline)

- For each registry element and each pset whose `applies_to` matches: evaluate the bound targets' getters, convert to display or system units (the `display_units_for_properties` rule), and write.
- Shared psets are created once and related to every matching element in one `IfcRelDefinesByProperties`.
- **Check before removing hard-coded definitions:** the standard table must reproduce today's export exactly. A script diffs psets (names, values, units, descriptions) between the old and new exports for all `Tests/ImportValidation` projects, in both property-unit modes. The round trip must still have 0 mismatches.
- **IDS (M4):** `IdsExporter` asks the table for the pset, property, and data type of each design value. Pinning and values stay as they are.

## Code layout

| File | Contents |
|---|---|
| `IfcElementRegistry.h/.cpp` | `ElementId`, `CIfcElementRegistry` (with R1 and M3) |
| `IfcTargets.h/.cpp` | `ElementKind`, `ValueKind`, `TargetValue`, `TargetDef`, the target list |
| `IfcMappingTable.h/.cpp` | finding table files, JSON load (nlohmann-json, already in vcpkg), the `extends` merge, validation, error messages |
| `IfcTargetReader.h/.cpp` | import engine |
| `IfcTargetHints.h/.cpp` | hint rules for targets that weren't found (G3) |
| `MappingTables/Standard.json` | standard table, installed next to the DLL (post-build copy for development builds) |
| `Tests/ImportValidation/mappings/Iowa.json`, `PennDOT.json` | agency tables (first versions in M1, completed in M2); `models.json` gets a `mapping` entry per model |
| `Tests/ImportValidation/models.json` | `PennDOT-StandardTable` and `Iowa-StandardTable` runs with `expect_log`: text the log must contain (the hints) |
