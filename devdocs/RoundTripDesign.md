# Round Trip: Design Note and Plan (R1)

Draft, 2026-09-30. The requirement is R1 in [IfcImportPlan.md](IfcImportPlan.md): importing an IFC model into PGSuper and exporting it again keeps the elements' `GlobalId`s, and loses nothing in the model. This note proposes how, with recommended answers to the open questions (marked **Recommended**; they need your review before implementation).

## What must be true after a round trip

1. **Element GlobalIds don't change.** Every element PGSuper models (bridge, bridge parts, alignment, referents, girders, deck, bearings, barriers) has the GlobalId it had in the source model.
2. **Nothing is lost.**
   - Elements PGSuper models keep every property set they came with. PGSuper updates only the properties the mapping table declares; a declared property set that's already on the element is updated in place, keeping its other properties.
   - Everything else in the model stays as it is: elements PGSuper doesn't model (wingwalls, piles, footings, approach slabs), schedules, tasks, cost items, documents, approvals, groups, systems, and their relationships to the elements.
3. **PGSuper's changes are in the result.** Values PGSuper owns (the declared properties with targets), the geometry and placement of the elements it models, elements added in PGSuper, and elements deleted in PGSuper.
4. **The source model isn't changed.** The result is a new file (B1).
5. **A PGSuper project with no source model** exports the same element GlobalIds every time.

## Approach

**Export, then merge into the source model.** The exporter stays as it is: it builds PGSuper's model in memory, as today. A new merge step then applies that model to a copy of the source model, element by element, and writes the result. A plain export (no source model) skips the merge.

Why not make the exporter write into the source model directly: the exporter creates entities in about 100 places, all assuming a new file. Merging keeps the exporter, the mapping tables, and the validation baselines unchanged, and puts every round-trip rule in one module that can be tested on its own.

**An element registry ties the two models together.** It maps a PGSuper element identity to a GlobalId and is saved in the PGSuper project:
- **Import** records the GlobalId of every IFC element that becomes PGSuper data.
- **Export** gives each element the registry's GlobalId (a new one if there's none) and records new ones, so the next export repeats them (requirement 5).
- **Merge** then matches the two models by GlobalId: an element in both is updated in place, an element only in PGSuper's model is added, and an element only in the source is kept, unless the registry says it was PGSuper's and PGSuper deleted it.

## Components

### 1. Element registry (`IfcElementRegistry.h/.cpp`, planned since M0)
- **Element identity** (`ElementId`): element role (the mapping tables' `ElementKind`) plus PGSuper's stable IDs, not indices, so inserting or removing a span or a girder doesn't shift the others:

  | Role | Key |
  |---|---|
  | project, site, bridge | (one each) |
  | superstructure, substructure, deck part | (one each) |
  | pier, abutment (bridge part), positioning referent | pier ID |
  | alignment | (one; its components by index) |
  | girder | segment ID |
  | girder assembly (spliced) | girder ID |
  | closure joint | closure ID |
  | deck | (one) |
  | bearing | pier ID, face, girder ID, bearing index |
  | barrier | side |

  Strands, tendons, and reinforcing bars have no stable identity in PGSuper (they're regenerated from the girder's data); they aren't in the registry (see "Strands and reinforcement").
- **Persisted** in the project by `CIfcExtensionAgent` (`IAgentPersist`): a new unit in the agent's data (agent version 1.0 → 2.0; version 1.0 projects load with an empty registry).
- Also persisted: the **source model** (path, size, time, SHA-256; see question 1) and a fingerprint of the imported alignment and georeferencing (see "Alignment and georeferencing").

### 2. Import records the identities
The project importer knows which IFC element became which PGSuper item (e.g. `get_girder_layout` assigns each `IfcBeam` to a span and a girder). After the bridge is built, a final pass records `ElementId → GlobalId` for every element it used, and the source model. Elements the import reads but PGSuper doesn't model one to one (e.g. Iowa's haunch slabs, measured into the deck) are recorded with their role, so the merge knows they were read (see "Haunches").

### 3. Export uses the registry
One GlobalId provider: the exporter asks it for an element's GlobalId instead of calling `ifcopenshell::global_id()`, for the elements in the table above (about 25 of the ~100 call sites; property sets, relationships, types, and materials keep new GlobalIds, B2). After the export, new identities are recorded in the registry. This marks the project modified, which is correct: the next export must repeat them.

The design-value IDS already pins girders by GlobalId as an option; with stable GlobalIds that becomes a sound default for round trips.

### 4. Merge (`IfcRoundTripMerge.h/.cpp`)
Input: PGSuper's model (in memory), a copy of the source model, the registry, and the mapping table. For each element of PGSuper's model:
- **In both** (same GlobalId): update the source instance in place:
  - geometry and placement replaced (question 3);
  - the declared properties updated in the element's, its type's, and its material's property sets (question 4); property sets the table declares that the element doesn't have are added; nothing else on the element is touched;
  - classifications the table declares are added if missing (existing ones are kept).
- **Only in PGSuper's model**: added, with its relationships (spatial containment in the corresponding source bridge part, aggregation, type, material, property sets, classifications).
- **Only in the source, recorded in the registry, not in PGSuper's model**: PGSuper deleted it (question 2).
- **Only in the source, not in the registry**: kept as it is.

Everything that isn't an element PGSuper models is kept: the merge only ever touches the instances above and their representations, placements, and declared property values.

**Units:** the source model's units are kept. PGSuper's model is built with the source's project units (the exporter's unit setup takes them from the source model), so copied geometry needs no conversion. Property values are written with their own units, as the export does today.

**Log:** every change goes to the export log: elements updated, added, and removed, properties changed (old and new value), relationships changed, and every decision the merge couldn't make cleanly.

### 5. User interface and command line
- **Export options:** when the project has a source model, the IFC export offers "Update the source model (keeps everything in it)" (the default) or "Export PGSuper's model only". The source model's path is shown.
- **Source model missing or changed** (question 1): interactive, the user locates it or chooses a plain export; from the command line the export fails with a message, unless `/IfcSourceModel=<file>` gives it (or `/IfcNoSourceModel` asks for a plain export).

## Recommended answers to the open questions

1. **Where the source model comes from at export.** **Recommended:** store its path, size, time, and SHA-256 in the project. At export, if the file is there and unchanged, merge into it. If it's missing or changed, ask the user to locate it (or choose a plain export); from the command line, fail with a message, or use `/IfcSourceModel=`. Don't store a copy in the project: IFC models are often tens to hundreds of MB, and the modeler's file stays the reference.
2. **Elements deleted in PGSuper.** **Recommended:** remove the element, what only it uses (its representation, placement, and property sets not shared with other elements), and its membership in shared relationships (a task's assignment, a group, an aggregation); delete a relationship left with no members. Log every removal, including each task, group, or other entity it was removed from, so the modeler can review it. Keeping deleted girders in the model would leave a model that contradicts the design.
3. **Geometry and placement.** **Recommended:** PGSuper is authoritative for the geometry of the elements it models: replace the representation and placement of the existing instance (keeping the instance, its GlobalId, and its relationships), whatever representation the source used. The replaced representation items are removed unless something else uses them. `Name`, `Description`, `Tag`, and `ObjectType` are kept from the source (the modeler's naming; selectors may depend on `ObjectType`); new elements get PGSuper's names.
4. **Declared properties.** **Recommended:**
   - Properties with a target: PGSuper's value always wins (it's PGSuper's data), including `"import": false` properties. When the source property has a unit or a value type of the same kind (e.g. a pressure in ksi), keep them and convert PGSuper's value; when the kind differs (e.g. a text "8" where PGSuper writes a pressure), replace the property and log it.
   - Constants (e.g. `Pset_BeamCommon.Status = "NEW"`): keep the source's value when the property exists (the source knows better, e.g. `EXISTING`); write the constant only when it's missing.
   - Placeholders (no value): never overwrite a source value.
5. **Types and materials.** **Recommended:** keep the source's type objects and materials for existing elements (they often carry agency data), and update only their declared properties. When a source material or type is shared by elements whose PGSuper values now differ (e.g. two girders sharing a material, one of which PGSuper changed to 10 ksi), give the differing elements their own copy, and log it. New elements use PGSuper's types and materials.
6. **Validation.** **Recommended:** a round-trip check in `Tests/ImportValidation` (see "Validation").

## Special cases

### Haunches
PGSuper exports the deck and its haunches as one `IfcSlab`. A source model can model haunches as separate elements (Iowa has 21). **Recommended:** when the source has separate haunch elements, keep them and replace the deck's geometry with the deck without haunches (a new exporter option for the deck representation); log that the haunch elements weren't updated. Updating the haunch elements' geometry from PGSuper's haunch depths is a later refinement.

### Strands and reinforcement
The import doesn't read strands or reinforcement yet (Phase 6). **Recommended:** in a round trip, PGSuper's strands and bars are written only for girders whose source element has none; a girder that has tendons or bars in the source keeps them, and the log says PGSuper's weren't written. Once the import reads them, they can be matched like other elements.

### Alignment and georeferencing
PGSuper's alignment and georeferencing come from the source model. **Recommended:** keep the source's (a civil model's alignment is richer than PGSuper's) unless the user changed them in PGSuper: compare with the fingerprint recorded at import, and replace them only when they differ (logged).

### Elements PGSuper added
A span or girder added in PGSuper has no GlobalId in the registry: it gets a new one, is added to the source model in the corresponding bridge part, and is recorded, so later exports repeat it.

## Validation

A new round-trip check in `run_validation.py`, for each IFC model (PennDOT, Iowa, and a fixture with schedules and tasks):
1. Import the model, export the project with the source model (`/IfcExport` merges), and compare:
   - **GlobalIds** of the elements PGSuper models: all the same.
   - **Kept elements' data:** every property set and property the table doesn't declare has exactly the source's value; declared properties have PGSuper's values.
   - **Everything else:** entity counts by type are the same for every entity type PGSuper doesn't create; the relationships to non-PGSuper entities (e.g. `IfcRelAssignsToProcess` of a task) have the same members.
   - The report lists every difference; the check passes when there are none besides PGSuper's declared values and geometry.
2. For PGSuper models: export twice and compare element GlobalIds (identical).
3. Edits: a scripted edit (delete a girder, add a span) between import and export, with the expected removals and additions.

**Fixture with schedules and tasks:** a small script (IfcOpenShell's Python, already in `D:\IfcOpenShell`) adds an `IfcWorkSchedule` with `IfcTask`s assigned to the girders and the deck, a cost item, a document reference, and an unrelated property set to an exported PGSuper model. It's the model the retention checks run on.

## Milestones

| | Milestone | Done when |
|---|---|---|
| RT0 | This note, reviewed; the open questions decided | |
| RT1 | Element registry, persisted in the project; import records identities and the source model | A saved and reopened imported project has a GlobalId for every imported element |
| RT2 | Export uses the registry | PGSuper models export the same element GlobalIds twice; IFC → PGSuper → IFC (plain export) gives the source's GlobalIds to the elements PGSuper models |
| RT3 | Merge: update in place, add, keep everything else; units; log | The round-trip check passes for PennDOT and Iowa (no losses) |
| RT4 | Deleted and added elements; haunches; strands and reinforcement; alignment and georeferencing | The edit checks pass; Iowa's haunches are kept |
| RT5 | Export options, locating a missing source model, `/IfcSourceModel=`, `/IfcNoSourceModel` | |
| RT6 | Schedules and tasks fixture, and its retention checks | The fixture's schedules, tasks, and their assignments survive a round trip |

Plain exports (no source model) stay the same through RT0–RT6, so the existing export and IDS baselines keep checking the exporter.

## Risks
- **IfcOpenShell copying between files:** adding an entity from PGSuper's model to the source model copies what it references; representations can share items (e.g. a profile). The merge must copy each shared item once (a map from PGSuper's instances to the copies), as `RebarRelationshipBatch` already does for relationships.
- **Removing entities** leaves dangling references if anything else points to them; the merge checks inverse references before removing anything.
- **Performance:** large models (hundreds of MB, or many bars). The merge works per element with GlobalId lookups (`instance_by_guid`), not whole-file scans.
- **Source models PGSuper can't fully read** (e.g. a girder PGSuper matched to the wrong span): the registry records what the import decided, so a wrong match is updated with the wrong girder's data. The import log already reports its matching; the merge log lists each updated element with its PGSuper identity so it can be checked.
