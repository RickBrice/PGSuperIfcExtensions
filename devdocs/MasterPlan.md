# Master Plan

All open work on the IFC import and export, in one list. Each item is here once, with its status and priority; the detailed plans and the design notes it links to have the rest. When an item is done, take it off this list and record what was done where its design is described. Organized 2026-11-20.

**Priorities** (the user's order, 2026-10-15; also in `CLAUDE.md`):
1. PGSuper import and export, with maximum usBridge IDS coverage.
2. Girder mapping and modeling with the PGSuper girder library.
3. Getting information from geometry when properties lack it.
4. PGSplice support equivalent to PGSuper's.

**Status:** **Next** (the next to work on), **Planned** (a plan exists), **Not started**, **Parked** (on hold, the user's choice or waiting for something), **Decision** (waiting on an answer), **Outside** (in another repository, or for someone else).

What's done is summarized in [IfcImportPlan.md](IfcImportPlan.md), "Status"; how things work is in the design notes ([Introduction.md](Introduction.md)).

## Large plans

| Plan | Status | Priority | Where |
|---|---|---|---|
| PGSplice import and export: S1 and S3 done; S0, S2, S4, S5 to do (below) | Planned | 4 | [PGSplicePlan.md](PGSplicePlan.md) |
| Stirrups: bar shapes and layouts as plug-in input, analysis through `CShearData2`, design, and a general bar export | Planned (open questions) | 1, 3 | [StirrupPlan.md](StirrupPlan.md) |
| 3D view of the bridge model | Exploratory (open questions) | | [3DViewPlan.md](3DViewPlan.md) |
| R2: import issues as BCF for the modeler | Not started (C1-C5 decided) | 3 | [IfcImportPlan.md](IfcImportPlan.md), R2 |
| Unit tests in Visual Studio (Phase 8): a static library of the testable core, a test project, and a first batch of tests (tables, bSDD dictionaries, the formatter, table <-> IDS, the target reader, bar geometry, extension registration) | Planned (one of the last steps) | | [UnitTestPlan.md](UnitTestPlan.md) |

## Import

| Item | Status | Priority | Details |
|---|---|---|---|
| End blocks and tapers from the geometry (sections near the ends; end block length, transition, width) | Next | 2, 3 | [GirderImportPlan.md](GirderImportPlan.md), G1 |
| Superelevation events (`IfcReferent` SUPERELEVATIONEVENT, `Pset_Superelevation`): export them, import them with the deck as fallback and check | Not started | 1, 3 | [ImportGeometry.md](ImportGeometry.md), roadway section |
| The girder orientation (roll): measured, not the factory's first orientation | Not started | 3 | [ImportGeometry.md](ImportGeometry.md), Section 9 |
| Strands: other representations (indexed poly curves); a model with debonding; round trips by girder family | Not started | 2, 3 | [GirderImportPlan.md](GirderImportPlan.md), G2 |
| Girder bars off the CL girder without a mirror bar (PGSuper rows are centered) | Not started | 3 | [GirderImportPlan.md](GirderImportPlan.md), G3 |
| Intermediate diaphragms: import them (fit the library's diaphragm rules to the model's), precast internal diaphragms of multi-web girders, tell them from pier diaphragms without `DiaphragmType` | Not started | 3 | [IfcImportPlan.md](IfcImportPlan.md), Phase 7 |
| Pier diaphragm dimensions (PGSuper computes them); PennDOT's end diaphragms are `IfcWall`s its table doesn't select | Not started | 3 | [IfcImportPlan.md](IfcImportPlan.md), Phase 2 |
| Integral connections: how to tell them from continuous ones (exported as continuous) | Decision | 1 | [IfcImportPlan.md](IfcImportPlan.md), Phase 2 |
| Sidewalks and interior barriers; barriers that vary along the bridge; barrier concrete | Not started | 3 | [IfcImportPlan.md](IfcImportPlan.md), Phase 7 |
| Stay-in-place deck panels; deck concrete f'c | Not started | 3 | [IfcImportPlan.md](IfcImportPlan.md), Phase 3 |
| The slab offset tolerance as an import option (D2) | Not started | | [IfcImportPlan.md](IfcImportPlan.md), Decisions |
| Pier caps and columns: the physical pier model | Not started (low) | | [IfcImportPlan.md](IfcImportPlan.md), Phase 7 |
| Alignments with spirals: export and evaluate with IfcOpenShell | Not started | | [IfcImportPlan.md](IfcImportPlan.md), Phase 0 |
| One reusable mesh → top surface → boundary → segments helper | Not started | | [IfcImportPlan.md](IfcImportPlan.md), Phase 1 |
| Girder fitting speed (about 2.5 s per distinct section): screen factories by topology; threads later | Parked | 2 | Details, below |

## Export

| Item | Status | Priority | Details |
|---|---|---|---|
| F1, girders as constructed or with camber (the haunch gap; the PennDOT haunch about 0.2 in thinner than the model's, and its overhang 0.1 to 0.25 in thinner within an inch of the edges, not traced) | Decision | 1 | [IfcImportPlan.md](IfcImportPlan.md), F1 |
| The deck reinforcement (bars from the lump sums, or a property set), so round trips check its import | Not started | 1 | [IfcImportPlan.md](IfcImportPlan.md), Phase 3 |
| Stirrups for any girder, with their real legs; a separate splitting zone; extended bars drawn projecting | Planned | 1 | [StirrupPlan.md](StirrupPlan.md); [GirderImportPlan.md](GirderImportPlan.md), G3 |
| File size: share repeated value entities (geometry, properties), optionally `.ifcZIP` | Parked (the user, 2026-11-19) | | Details, below |

## Round trips

| Item | Status | Priority | Details |
|---|---|---|---|
| Geometry and placement of the model's elements (decision 3) | Not started | 1 | [RoundTripDesign.md](RoundTripDesign.md) |
| Elements deleted in PGSuper (decision 2); elements added (girders, spans, decks, barriers, piers) | Not started | 1 | [RoundTripDesign.md](RoundTripDesign.md) |
| Haunches modeled separately (the deck's representation without them) | Not started | | [RoundTripDesign.md](RoundTripDesign.md) |
| The alignment and georeferencing when changed in PGSuper | Not started | | [RoundTripDesign.md](RoundTripDesign.md) |
| Stable GlobalIds for plain exports | Not started | | [RoundTripDesign.md](RoundTripDesign.md) |
| A fixture with schedules and tasks, for the retention checks | Not started | | [RoundTripDesign.md](RoundTripDesign.md), Validation |

## Mapping tables and IDS

| Item | Status | Priority | Details |
|---|---|---|---|
| A test with a real agency IDS (generate a table, bind it, import with it) | Not started | 1 | [MappingTablesDesign.md](MappingTablesDesign.md), M5 |
| The mapping table editor's documentation: publish it online, and add it to `UpdateBinFiles.bat` and the installer | Not started | | [MappingTablesDesign.md](MappingTablesDesign.md), M7 |

## Export extensions

| Item | Status | Priority | Details |
|---|---|---|---|
| The design-value IDS for extension targets | Not started | | [ExtensionsDesign.md](ExtensionsDesign.md) |
| Extension hooks for PGSplice bridges | Not started | 4 | [ExtensionsDesign.md](ExtensionsDesign.md) |
| General PGSuper data that extensions register as their own targets, into this project's catalog; a declarative `IServices::CreateElement` for the elements they build directly; then no IfcOpenShell types in the hooks | Not started | | [ExtensionsDesign.md](ExtensionsDesign.md) |
| Release API version 1 (it's provisional until then) | Not started | | [ExtensionsDesign.md](ExtensionsDesign.md), ABI |

## PGSplice

The steps are in [PGSplicePlan.md](PGSplicePlan.md).

| Item | Status | Priority |
|---|---|---|
| S0: the three test models, `.spl` models in `run_validation.py` (with PGSplice's configuration switched and restored), `pgs_extract.py` for PGSplice input | Planned | 4 |
| S2: closure joints (type, classification, material, reinforcement, solid placed between the segment ends), segment classification, post-tensioning, erection tower height and supported segments, the construction sequence, the single-girder export, the design-value IDS, the round trip keys | Planned | 4 |
| S3 leftovers: tendons, the timeline, closure joint concrete, slab offsets per segment, a temporary support's end distance, skewed temporary supports | Planned | 4 |
| S4: post-tensioning import and the timeline | Planned | 4 |
| S5: the data importer and round trips under PGSplice | Planned | 4 |
| Later: an `IfcWorkSchedule` export option; where the duct properties PGSuper doesn't have come from | Not started | 4 |

## Validation and tests

| Item | Status | Priority | Details |
|---|---|---|---|
| The E2 validation template (deliberately unusual values, fewer "default matches") | Not started | | [IfcImportPlan.md](IfcImportPlan.md), Decisions, E2 |
| Expected values for PennDOT and Iowa from their property sets (E1) | Not started | | [IfcImportPlan.md](IfcImportPlan.md), Phase 0 |
| More test models: PGSuper regression projects by girder factory, and imports under other configurations | Not started | 2 | Details, below |
| Document the scoring systems (validation scores, name similarity, section comparison and fitting, usBridge counts, workbook categories) | Not started | | Details, below |
| Unit tests (Phase 8) | Planned | | [UnitTestPlan.md](UnitTestPlan.md) |

## Outside this repository

| Item | Who | Details |
|---|---|---|
| IfcOpenShell drops the deck slab solid of the `PGSuper_Skew_Straight_AlongPier` export (the OpenCASCADE loft fails silently); needs a debug build | IfcOpenShell | Details, below |
| `loft.cpp:176-180` on the `rab_infrastructure` branch prints the loft sections to `std::wcout` (leftover debug output) | IfcOpenShell | |
| The GenericBridgeTools test program crashes in its older section cut tests (`TestPrecastGirderBridge` next) | WBFL | |
| `GirderLibraryEntry` strand coordinate functions return the count, not the index their comments say | PGSuper (the user) | Details, below |
| Extended strands with individual strands (`sdtDirectStrandInput`): `GetExtendedStrandCount`, `IsExtendedStrand`, `GetNumExtendedStrands` may look at the wrong rows | PGSuper | [GirderImportPlan.md](GirderImportPlan.md), G2 |
| usBridge: `BeamLongitudinalSlope` is an `IfcNormalisedRatioMeasure` (0 to 1), so a girder sloping down can't be written; the strand and bar type specifications contradict the "Reinforcing Bar" and "Tendon" ones; the PGSplice shortcomings (temporary supports, construction sequence, closure joint pours, duct and anchor properties) | usBridge authors | `CLAUDE.md`; [PGSplicePlan.md](PGSplicePlan.md), "usBridge shortcomings" |
| GRF005 checks only the unit ratio of `IfcMapConversion.Scale` (buildingSMART/ifc-gherkin-rules#522) | buildingSMART | `CLAUDE.md` |
| Bonsai's IDS auditor matches a tendon type with a material specification (left visible on purpose) | Bonsai | `CLAUDE.md` |

## Details of the smaller items

Items without a design note of their own.

- **Export file size: repeated value entities** (found 2026-11-19, not now, the user): about half of the `PGSuper_Import_Model` export with rebar (3.1 MB of 6.1 MB) is entities that exactly repeat one already in the file, mostly from the bar types (one per bar shape and length) and the confinement bars at every stirrup. The export dialog leaves the bars out by default, but `/IfcExport` and round trips include them. A value cache in `CIfcExportSession` (as it has for enumerations) would return an existing entity with the same values:
  1. Geometry (safe for every reader, about 1.3 MB): `IfcDirection`, `IfcCartesianPoint`, and `IfcCartesianTransformationOperator3D`, mostly from `MapBarType` (`Rebar.h`); the "up" and axis directions of placements and girder cross sections; hierarchy_helper's `addLocalPlacement`/`addPlacement3d`; identical `IfcMappedItem`s between girders.
  2. Properties and quantities (about 1.5 MB): `IfcPropertySingleValue`, `IfcPropertyEnumeratedValue`, and `IfcQuantityCount` from `IfcPropertyWriter.h`. Sharing them between property sets is valid IFC, but the round-trip merge edits property values in place, so it would have to copy a shared property before editing it; editors such as Bonsai would change every set that shares an edited property. Sharing whole identical property sets between bar types is the alternative.
  3. Optionally `.ifcZIP` output (IFC text compresses to about a tenth).
- **Girder fitting speed:** about 2.5 s per distinct section in Release, mostly screening about 50 factory and start candidates. Next: measure per factory, then screen factories by topology before fitting (voids, open top, number of webs). Threads stay an option (each worker would need its own COM initialization and sections, the factories would need checking for thread safety, and `QuietChecks` changes a process-wide CRT setting).
- **More test models from PGSuper projects:** round-trip a few of the 730 projects in `PGSuper/RegressionTest` for each girder factory (I-beam 293, box beam 90, voided slab 85, deck bulb tee 60, U-beam 82, double tee 29, ribbed 13, decked slab 9), and import some under another configuration (e.g. a TxDOT project under WSDOT; agency libraries are in `PGSuper/Configurations`), so entries are fitted and created rather than matched.
- **Document the scoring systems** (asked 2026-10-15): the validation scores of `compare_bridge.py` and their tolerances; the girder name similarity (`girder_name_similarity`, 0.5 minimum); the section comparison and ranking (`CompareSections`, `SectionDifference`, `SectionTolerances`) and the fitting objective (`FitDimensions`); the usBridge IDS counts and the workbook comparison categories (`check_standard.py`).
- **IfcOpenShell drops the main deck slab solid** of the `PGSuper_Skew_Straight_AlongPier` export without a message, so its deck imports as a haunch (the run's mismatches). The deck sections are valid (constant point counts, no self intersections, inside the deck edges); the failure is in the OpenCASCADE loft conversion (`kernels/opencascade/loft.cpp`), which returns false silently. Needs an IfcOpenShell debug build. Also, `loft.cpp:176-180` on the `rab_infrastructure` branch prints the first two loft sections to `std::wcout` on every loft (leftover debug output).
- **PGSuper documentation error** (for the user to fix in PGSuper): in `Include/psgLib/GirderLibraryEntry.h`, `AddStraightStrandCoordinates`, `AddHarpedStrandCoordinates`, and `AddTemporaryStrandCoordinates` say they return "the strand grid index for the new strand position", and `AddStrandToPermStrandGrid` says it returns "the position in the permanent strand grid", but `GirderLibraryEntry.cpp` returns the number of positions, one more than the index. `GenerateStrandGrid` subtracts one. Either the comments or the return values should change (callers such as TxDOT TOGA may depend on the current values).
