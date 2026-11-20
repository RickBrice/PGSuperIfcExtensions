# Round Trip: Design (R1)

The requirement is R1 in [IfcImportPlan.md](IfcImportPlan.md): a project imported from an IFC model can be exported as that model with PGSuper's design in it, keeping the elements' `GlobalId`s and losing nothing in the model. Designed 2026-09-30; the first version ("properties plus missing elements", the user, 2026-10-26) is built. The user documentation is `Documentation/MappingTableEditor/src/round_trip.md`.

## What a round trip guarantees

1. **The model itself isn't changed.** The export is a copy of it (B1), so IfcOpenShell's IfcDiff (`src/ifcdiff/ifcdiff.py`) shows what PGSuper changed.
2. **Element GlobalIds don't change.** Every element of the model keeps its GlobalId.
3. **Nothing is lost.** Elements PGSuper models keep every property set they came with; PGSuper updates only the properties the mapping table declares. Everything else in the model stays as it is: elements PGSuper doesn't model (wingwalls, piles, approach slabs), schedules, tasks, cost items, documents, groups, and their relationships.
4. **PGSuper's design is in the result:** the values the table binds to targets, and the elements PGSuper designs that the model doesn't have (see "The merge").

## How it works

**Export, then merge into a copy of the source model.** The exporter builds PGSuper's model in memory as for any export. The merge (`IfcRoundTripMerge.h/.cpp`, run by `CIfcExporter::BuildModel` when `CIfcExportOptions::source_model` is set) then applies it to a copy of the source model, element by element, and writes the result. A plain export skips the merge, so the export and IDS baselines keep checking the exporter.

### The registry

`CRoundTripRegistry` (`IfcRoundTrip.h/.cpp`) maps keys of PGSuper's elements to the GlobalIds of the source model's elements. The keys use PGSuper's IDs, not indices, so editing spans or girders doesn't shift them (`roundtrip_key`, e.g. `girder:<segment ID>`, `bearings:<segment ID>:start`, `diaphragms:<pier ID>`, `intermediate_diaphragms:<start pier ID>`, `tendons:<segment ID>`, `bars:<segment ID>`).

- **The import records** (`CIfcImporter::GetRoundTripRecorder`, `CRoundTripRecorder`) the project, site, bridge, bridge parts, piers and abutments, girders, deck, exterior barriers, bearings (by girder end), pier diaphragms, intermediate diaphragms (by span), and which girders have strands or bars in the model.
- **The project saves** the registry with the source model's path, size, and SHA-256 (`IIfcRoundTrip`, `CIfcExtensionAgent` unit version 3).
- **The export names** each element it writes by its key (`IfcExportIdentities`, `CIfcExportOptions::export_identities`), and the merge pairs them with the registry's GlobalIds. PGSuper's own model keeps new GlobalIds on every export; only the copy has the model's.

### The merge

- **Elements in both:** the model's element gets the table's property sets, quantity sets, material properties, and classifications:
  - properties bound to a target get PGSuper's value, converted to the model property's unit (or the model's project unit). If the model's property has another kind of value (text where PGSuper writes a pressure), PGSuper's replaces it.
  - constants and properties with no value are added only where the model doesn't have the property; the model's value is kept.
  - a classification of the table is added where the element doesn't have it.
  - a property set shared with other elements is copied for the element before it changes; a type's property set is updated on the type unless elements of that type need different values (then on the element, which overrides it); a shared material whose elements need different values is copied, one copy per set of values (Iowa's span 2 girders, 10 ksi in a 5 ksi material).
  - the model's names, geometry, placements, types, and everything else are kept.
- **Elements PGSuper adds:** the strands and bars of girders without them in the model, bearings at girder ends without them, pier diaphragms at piers without them, and intermediate diaphragms in spans without them. A girder with strands or bars in the model keeps the model's; diaphragms the model has at a pier get PGSuper's diaphragm values; a span with its own intermediate diaphragms keeps them as they are (logged: they aren't matched one by one yet).
  - They're copied with IfcOpenShell's `add_entity` (which converts `IfcLengthMeasure` attributes to the model's length unit; property and quantity values are converted by the merge).
  - Their placements are made relative to the model's element they belong to (or the world) wherever they pointed outside the added elements (PGSuper's girder placement, linear placements along PGSuper's alignment), from the computed world matrices; the export's world coordinates are the model's.
  - Their relationships are rebuilt with only the added elements; representation contexts and classification systems are the model's.
  - PGSuper's tendon bundles and reinforcement cages are matched with the model girder's assemblies of the same role (by their strands or bars, else by the table's `tendon_bundle` and `reinforcement_cage` selectors), and updated like other matched elements, so they get the classification and property sets even when PGSuper has nothing to add. PGSuper's strands and bars go into a matched assembly that's empty, so it keeps its GlobalId and type.
- **The log:** every change goes to the export log (`<file>.export.log` beside the IFC file from the export dialog): each property changed (old and new value), each property and classification added, each element added, and anything left alone and why.

### The design-value IDS of a round trip

`CIdsExportOptions::round_trip`: the girders are pinned by GlobalId (the copy keeps the model's names); their specifications don't require PGSuper's material (the model's is kept); the girder type and concrete material specifications are left out; the strand and bar specifications are left out when a girder keeps the model's strands or bars.

### User interface and command line

- **Export options:** the IFC tab's Round Trip group, "Update a copy of the IFC model this project was imported from", is on when the model is where it was and unchanged since the import; "..." locates a moved or changed model, which is used as it is. While it's on, the tab's other options are disabled: a round trip uses their defaults. The export won't replace the model itself.
- **Command line:** `/IfcSourceModel` (the project's model, which must be unchanged) or `/IfcSourceModel=<file>`. An export extension's export can't be a round trip.

## Decisions (2026-10-26)

1. **The source model** is found by the path stored in the project, checked by size and SHA-256. A copy isn't stored with the project: IFC models are often tens to hundreds of MB, and the modeler's file stays the reference.
2. **Elements deleted in PGSuper** are to be removed, with what only they use and their membership in shared relationships, each removal logged. Not built yet (they're kept).
3. **Geometry and placement:** PGSuper is to be authoritative for the geometry of the elements it models, replacing the representation and placement of the existing instance and keeping its name, description, tag, and `ObjectType`. Not built yet (the model's geometry is kept).
4. **Declared properties:** targets win; constants and properties with no value never overwrite the model's values (built).
5. **Types and materials:** the model's are kept and only their declared properties updated, with copies where the elements sharing them need different values (built).

## Not done yet

- Geometry and placement of the model's elements (decision 3).
- Elements deleted in PGSuper (decision 2).
- Elements added in PGSuper that the model doesn't have, such as girders, spans, decks, barriers, and piers (logged, not added). Once added they'd be recorded, so later exports repeat their GlobalIds.
- Haunches: a model's separate haunch elements (Iowa's) are kept, but the deck's representation isn't changed to leave them out.
- The alignment and the georeferencing: the model's are kept. To replace them when the user changed them in PGSuper, compare with a fingerprint recorded at import.
- Stable GlobalIds for plain exports (a PGSuper project with no source model exporting the same GlobalIds every time).
- PGSplice: closure joints, tendons, ducts, and anchors as keyed elements, and the bearings of an import into PGSplice ([PGSplicePlan.md](PGSplicePlan.md), S2 and S5).
- A fixture with schedules and tasks for the retention checks (a script that adds an `IfcWorkSchedule` with `IfcTask`s assigned to the girders and the deck, a cost item, and a document reference to an exported model).

## Validation

`Tests/ImportValidation/check_roundtrip.py` (run by `run_validation.py`; see its README):
- **IFC models (PennDOT, Iowa, TxDOT):** each is exported as a copy after its import and checked: the model file unchanged, every GlobalId kept, no entity type with fewer instances, and the copy passes its design-value IDS.
- **Added reinforcement:** the PGSuper model's export, stripped of its strands and bars, gets them back exactly where an export of the project puts them.

## Risks

- **Copying between files:** adding an entity from PGSuper's model copies what it references, and representations can share items (e.g. a profile). Each shared item is copied once.
- **Removing entities** (decision 2) would leave dangling references if anything else points to them: check inverse references before removing anything.
- **Performance:** large models (hundreds of MB, or many bars). The merge works per element with GlobalId lookups, not whole-file scans (Iowa's 65 MB model takes about 30 s).
- **Models PGSuper can't fully read** (e.g. a girder matched to the wrong span): the registry records what the import decided, so a wrong match is updated with the wrong girder's data. The import log reports its matching, and the merge log lists each updated element with its PGSuper identity.
