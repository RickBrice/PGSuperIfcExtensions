# IFC Import and Export Plan

Drafted 2026-09-23 as the import improvement plan; it now tracks the IFC import and export as a whole. How the import reads geometry is in [ImportGeometry.md](ImportGeometry.md), the mapping tables in [MappingTablesDesign.md](MappingTablesDesign.md), and the targets in [Targets.md](Targets.md).

## Status (2026-11-20)

Everything needed to resume is in the repositories; working rules, build and test notes are in `CLAUDE.md`. All repositories are on `develop`.

**Priorities (2026-10-15):** (1) PGSuper import and export with maximum usBridge IDS coverage, (2) girder mapping and modeling with the PGSuper library (Phase 5), (3) getting information from geometry when properties lack it (Phases 4-6), (4) PGSplice support equivalent to PGSuper's ([PGSplicePlan.md](PGSplicePlan.md)).

**Done:**
- **Phase 0, test setup:** headless import and export, and `Tests/ImportValidation` (round trips of PGSuper projects, imports of the PennDOT, Iowa, and TxDOT models, the standard-table runs, export and IDS baselines, the checks against the usBridge IDS files and workbook, schema and unit checks).
- **Phase 1, correctness:** units, girder layout, spans, spacing.
- **Phase 2, superstructure geometry:** pier skew, girder end distance and bearing offset, uniform spacing (normal to the alignment or along the CL piers), bearing data, and boundary conditions (continuity from the pier diaphragms, hinges and rollers from the support conditions). Integral connections are exported as continuous for now.
- **Phase 3, deck:** cast-in-place deck dimensions and slab offset; deck reinforcement as lump sums; the deck and haunches export as one `IfcSlab`, and the import measures haunches however a model has them; decks of several slabs (joints across or along the bridge).
- **Phase 4, roadway from the deck:** the roadway section (cross slopes and superelevation transitions), and the profile when the model has no vertical alignment.
- **Phase 5, girder section** ([GirderImportPlan.md](GirderImportPlan.md), G1): the section from the model's swept profile or a cut through its solid, matching library entries by geometry, and new entries with fitted dimensions and a generated strand grid, named "<type> (IFC)".
- **Phase 6, strands** (G2): from their centerlines, on the library grid when every strand is on it, else individual strands; material, jacking, debonding, extended strands.
- **Girder reinforcement** (G3): longitudinal bars, stirrups, confinement bars, and splitting reinforcement, by geometry.
- **Phase 7, exterior traffic barriers** (matched to library entries by section, or new entries made from it), and the export of PGSuper's intermediate diaphragms (2026-11-19).
- **Several structures in one bridge, support lines from the geometry** (2026-11-04, the TxDOT US 77 model): the structure on the alignment is imported (`select_structure`), and supports are found where girder ends and bearings meet when the bridge parts don't say.
- **R1 round trip, first version** (2026-10-26): an export can be a copy of the IFC model the project was imported from, with PGSuper's design in it ([RoundTripDesign.md](RoundTripDesign.md)).
- **Mapping tables M0-M7** (2026-09-29/30): import and export engines, the standard table and the agency tables, the general IDS commands, choosing the table, and the mapping table editor ([MappingTablesDesign.md](MappingTablesDesign.md)). bSDD dictionaries are table data (several, 2026-11-20).
- **usBridge conformance:** strand and bar data on their types, bearing and project data, `BeamLongitudinalSlope` as a magnitude, `SacrificialWearingThickness`, the diaphragms' quantities, barrier concrete, and property units IDS checkers read right (2026-11-19). The usBridge IDS fails 10 of 599 specifications for the PGSuper exports (9 for the skew models); see the `usBridge-ids.txt` results for which, and "Strand and bar data on their types" in `CLAUDE.md` for the inconsistency that's left visible.
- **Export extensions and export profiles** (2026-11-19/20): another plug-in can add to the export ([ExtensionsDesign.md](ExtensionsDesign.md)); every export is a profile of one build, tested with a test extension ([ExportProfilesPlan.md](ExportProfilesPlan.md)).
- **PGSplice:** registration (S1), erection towers in the export, and the import of spliced girders (S3), 2026-11-10/11 ([PGSplicePlan.md](PGSplicePlan.md)).
- **Import log:** the Import Results window can save the log, and View > IFC Import Log shows the project's import log again until BridgeLink closes (`IIfcImportLog`; not saved with the project).
- **Other repositories:** WBFL `ISectionCutTool::CreateSlabShapeWithoutHaunches`/`CreateHaunchShape` (PGSuper `IShapes::GetSlabShapeWithoutHaunches`/`GetHaunchShape`), skewed section cuts in COGO and GenericBridgeTools, `UBeam`/`UBeam2` with vertical webs, and the AutoSave status bar fix; the PGSuper Test and TxDOT agent command line fix (F2); IfcOpenShell's false GEO 242 warning (592ea02f7, #9590).

**Open work:** every open item, with its status and priority, is in [MasterPlan.md](MasterPlan.md).

**How to validate:** Release build of BridgeLink with this extension (`%ARPDIR%\BridgeLink\RegFreeCOM\x64\Release`), then `python Tests/ImportValidation/run_validation.py` (see its README). Last run 2026-11-20 (`WSDOT:WSDOT` configuration): the PGSuper round trips 460 to 478 of 623 values match, with no mismatches except the `PGSuper-ExteriorGirders` run's 2 (skew models: 273, 301, and 302 of 619; `PGSuper-Skew`'s 8 mismatches are bearing offsets without bearings, by design); PennDOT 78 of 115, Iowa 250 of 293, and TxDOT 1673 of 2602 values set by the import; every export the same as its baseline, valid, with valid property units and its own design-value IDS passing; the round trips into the source models keep them and pass their design-value IDS; the export extension API checks pass; the standard table round trips through a general IDS (361 declarations). Under the Regression configuration, the PGSuper round trips have 13 harped strand mismatches, by design (WSDOT's bundled harped strands don't fit a generated grid, G2). Use `--config-file` with a copy of the template that has unusual values to show a value is set by the import rather than a template default.

## Requirements

### R1: Clean round trip
Importing an IFC model into PGSuper and exporting it again must keep the elements' `GlobalId`s and lose nothing in the model: downstream tools (model federation, change tracking, BCF issues, IDS pinned by GlobalId) identify elements by GlobalId. A round trip updates a copy of the source model rather than creating a new one: the elements and property sets that exist keep their instances (and GlobalIds), PGSuper replaces what it owns (the values of the declared properties; later, the geometry of the elements it models), and adds what's new. The source model isn't changed (B1). Built as described in [RoundTripDesign.md](RoundTripDesign.md), which lists what's still to do.

### R2: Import issues as BCF for the modeler (optional)
When the import finds missing data or has to make an assumption, the user can have it write a BCF file (BIM Collaboration Format) with one topic per issue. The modeler opens it in their authoring or review tool, and each topic selects the affected elements by IFC GlobalId and zooms to them. The issues go back to the people who can fix the model, instead of staying in the PGSuper import log.

**Issues the import detects** (today they're log lines): no alignment in the model; a pier station not given; `NumberOfSpans` missing or in conflict with the model; a girder designation that disagrees with its location; a girder type not in the library (an entry is created); f'c or f'ci not found; a span with fewer than 2 girders; a deck dimension that can't be measured; a deck that doesn't reach the bearings; a unit not defined in the model; geometry that couldn't be processed; strands and bars left out, and why.

**Design:**
- **A structured import report first.** Replace the ad-hoc `Logger` calls with an issue collector. Each issue has a category (missing data, assumption, conflict, geometry), a severity (error, warning, info), a title and description, the affected IFC elements (GlobalIds), the PGSuper item, and the value used and where it came from (model, geometry, derived, default). The collector feeds the results dialog, the import log, a JSON report for `Tests/ImportValidation`, and the BCF writer.
- **BCF writer:** one topic per issue kind, with all affected elements selected in its viewpoint, rather than one topic per element (21 girders missing f'ci is one issue, not 21); the description lists each element and the value used; topic type and priority from the severity; labels `PGSuper import` plus the category; one viewpoint per topic (the components select the elements by `IfcGuid`, plus a camera fitted to the elements' bounding box); no snapshots at first.
- **Options:** `CIfcImportOptions::bcf_file` (empty = no BCF, the default; C5) and the BCF version (2.1 or 3.0; C1); an import options checkbox, with the path defaulting to `<model>.bcf`; command line `/IfcBcf=<file.bcf>`.
- **Mapping tables:** a "missing data" topic can say where the table (and the agency IDS it came from) expects the value, e.g. "f'ci expected in `IaDOT_PPCB`.`6_Concrete Release Strength, Fci`".
- **Validation:** `run_validation.py` writes the BCF for each IFC model and checks it: zip structure, XML against the BCF schemas, one topic per expected issue kind, and every referenced GlobalId present in the model.

**Decided:** C1-C5 (see Decisions). Test with Blender/Bonsai (C4).

### R3: Mapping table editor
Done (2026-09-30, M7): a BridgeLink application, like the PGSuper Library Editor, that creates and edits tables, picks locations from a model, tries a table on a model, and generates tables from IDS files and writes them as IDS. See [MappingTablesDesign.md](MappingTablesDesign.md), "The mapping table editor (M7)", and the user documentation in `Documentation/MappingTableEditor`. Decided 2026-09-29: a BridgeLink application (R3.1), files only, no web or catalog servers (R3.2), MFC (R3.3).

## Decisions (2026-09-24)
Answers to the open questions (`IFC_Import_Open_Questions.docx`). The IDs are used in the rest of the design notes.

**Mapping tables and IDS**
- **A1 Location:** the configured mapping table is chosen on the IFC Extension page of File > Configure BridgeLink (registered by the mapping table editor plugin, which lives for the whole application) and stored in the registry (`HKCU\...\BridgeLink\IfcExtensions\MappingTable`, shared by PGSuper and PGSplice). The file can be anywhere. `/IfcMapping=<file>` overrides it on the command line. (First planned as a page of the configuration wizard, then an Options menu command; changed 2026-10-15.)
- **A2 Selection:** the configured table, the table chosen for one import or export, or the project's table (M6). No automatic detection.
- **A3 Format:** JSON.
- **A4 IDS binding:** a separate binding (overlay) file. Agency IDS files are not modified.
- **A5 Exported content:** property sets, quantity sets, and classifications are table driven. Geometry, the spatial structure, and element names and `ObjectType` stay in code.

**Clean round trip (R1)**
- **B1 Unmodeled elements:** carry them through. Keep the source model and merge the PGSuper export into a copy of it. The original IFC stays unmodified.
- **B2 Secondary objects:** only element GlobalIds need to be stable; there are no derived GlobalIds. A round trip keeps everything in the source model, so property sets, relationships, types, and materials already in it keep their instances and GlobalIds; new ones get new GlobalIds (revised 2026-09-30).

**BCF (R2)**
- **C1 Version:** both 2.1 and 3.0, selectable.
- **C2 Topics:** one topic per issue kind, selecting all affected elements.
- **C3 Zip:** a vcpkg library (libzip or minizip-ng). PGSuper's MakePgz zips with a bundled copy of Lucian Wischik's `zip.cpp`/`unzip.cpp`; aim for one zip solution across the projects.
- **C4 Tools:** probably Blender/Bonsai. Test the BCF there.
- **C5 Default:** off. Turned on in the import options.

**Import behavior**
- **D1 Conflicts:** geometry wins for layout (stations, spacing, dimensions), property sets win for materials, and every conflict is reported.
- **D2 Slab offset tolerance:** 1/4 in, to become an import option.
- **D3 Girders not in the library:** command line: create the entry and log it. Interactive: create the entry; the user picks among factories that fit equally well.

**Validation**
- **E1 Drawings:** not available. Expected values for PennDOT and Iowa come from the model property sets only.
- **E2 Validation template:** create a template with deliberately unusual values in PGSuper and add it to `Tests/ImportValidation`.

**Exporter and PGSuper**
- **F1 Deck haunch vs girder camber:** left for now. The question is whether girders are modeled as constructed (no camber) or in the expected final condition (with camber); it may become an export option.
- **F2 TxDOT command line:** done 2026-09-24 (PGSuper `develop`, 1c45952285): the Test agent (`/TxA...`) and the TxDOT agent (`/TxTOGA`) claim only their own commands.
- **F3 Exporter placeholders:** fixed (the stress unit name, the barrier properties).

**Other decisions** (in `CLAUDE.md`, "Design decisions to keep"): the import template, the design-value IDS vs a general IDS, library reinforcement as a seed value, confinement bars, intermediate diaphragms, rebar off by default in the export dialog, girder families, the usBridge workbook vs the IDS, the tendon type and material names, and the GRF005 warning.

## Test models

| | PGSuper exports | PennDOT Rearick Rd | Iowa US59 | TxDOT US 77 |
|---|---|---|---|---|
| Source | the round-trip projects in `Tests/ImportValidation/models` | `BIM for Bridges\Precast Concrete Girder Modeling\PennDOT Rearick Rd\323010451_STR1.ifc` | `...\Iowa ADCMS\Iowa ADCMS Pilot_US59 over IA92_2025_07_09.ifc` | `...\TxDOT` (zipped in `models`) |
| Authoring tool | PGSuper | OpenBridge (Quadri) | unknown, 65 MB | unknown |
| Alignment | yes, with `IfcReferent` stationing on each pier | none | none (has `IfcMapConversion`) | yes |
| Spans / girders | 3 spans, 4/5/4 girders (the main model) | 1 span, 5 girders, integral abutments | 3 spans × 7 girders, 3-column piers with caps | 27 spans; two structures side by side, each deck in 9 slabs |
| Beam geometry | `IfcSectionedSolidHorizontal` | `IfcFacetedBrep` | mapped `IfcTriangulatedFaceSet` | |
| Beam type name | `IfcBeamType` "WF66G" | `_PS Concrete Beams.Type` = "PA Bulb-Tee 33/31.25" | `IaDOT_PPCB.2_Type` = "BTB45" and classification | in its table (`mappings/TxDOT.ifcmap`) |
| f'c / f'ci | standard property sets, per-property units | custom property set, text values ("8", "6.8") | custom property set, unitless numbers in ksi | in its table |
| Deck | 1 FLOOR slab with the haunches | 1 FLOOR + approach slabs | 1 FLOOR + 21 haunch slabs (`IfcSlab.USERDEFINED` "Haunch") | 9 slabs per structure, between joints |

## Phases

How each step reads the geometry is in [ImportGeometry.md](ImportGeometry.md); the girder section, strands, and girder bars in [GirderImportPlan.md](GirderImportPlan.md). What's left in each phase:

### Phase 0: Test setup
- **Headless import:** `BridgeLink.exe /IfcImport=<model.ifc> <template.pgt|template.spt> [/IfcOut=<project.pgs|project.spl>] [/IfcLog=<import.log>] [/IfcMapping=<table.ifcmap>]`. `CIfcExtensionAgent` implements `IEAFProcessCommandLine`; it imports into the template-created project, saves it, and EAF exits. `/IfcImport` comes first, for PGSuper builds older than F2.
- **Headless export:** `BridgeLink.exe /IfcExport=<model.ifc> <project.pgs> [/IfcPropertyUnits=Display|System] [/IfcLog=<export.log>] [/IfcMapping=<table.ifcmap>] [/IfcIds=<design.ids>] [/IfcSourceModel[=<model.ifc>]]`, with the default `CIfcExportOptions`.
- **Validation tooling:** `Tests/ImportValidation` (see its README).

### Phase 1: Correctness of what we already import
Done: one geometry pass over the girders (`get_girder_layout`, the single source of girder keys), girders assigned to spans by pier stations and ordered by offset, model designations cross checked (the location wins), units from the property, else the project unit (`CIfcImportUnits`).

### Phase 2: Superstructure geometry
Done: pier skew (from the positioning referent's `RefDirection`, the CL bearing lines, or the girder ends), end distance and bearing offset from the bearings, uniform or general spacing, bearing data (`CBearingData2`), and boundary conditions:
- **Boundary conditions** (decided by the user, 2026-10-23): the export writes PGSuper's pier diaphragms as `IfcBeam.DIAPHRAGM` (usBridge's "Concrete Diaphragm", in the superstructure, not extended to the support) and the pier's boundary condition as `usBrPset_BearingCommon.SupportCondition` on the bearings and `PGSuperPset_Diaphragm.SupportCondition` on the diaphragms. `Pset_BearingCommon` is the elastomeric bearing's fixity (`CBearingData2::FixedX/FixedY`), not the pier's boundary condition. The import makes a pier continuous after deck when a diaphragm holds the girder ends of both spans, else uses the bearings' support condition, else the diaphragms', else a hinge at the first abutment and rollers elsewhere, each logged with its evidence.

### Phase 3: Deck
Done: cast-in-place deck dimensions (gross depth, edge depths, overhang taper, haunch shape and fillet), slab offset at the CL bearings, the deck and haunches in one `IfcSlab` in the export (cut along skewed lines on skewed bridges), deck reinforcement as lump sums, decks of several slabs.

### Phase 4: Alignment and profile when the model has none
Done: the roadway section and the profile from the deck, a straight alignment from the deck when the user asks. Superelevation events aren't written or read yet: `IfcReferent` SUPERELEVATIONEVENT with `Pset_Superelevation` (`IfcAlignmentCant` is for railways), with the deck as the fallback and check, as for the profile, and deck station property sets (Iowa) for the reference station.

### Phase 5: Girder section and material mapping
Done (G1). Geometry decides; names are clues.

### Phase 6: Prestressing
Done (G2).

### Phase 7: Secondary elements and loads
- **Exterior traffic barriers:** done; sidewalks and interior barriers are logged, not imported.
- **Intermediate diaphragms:** exported (2026-11-19); round trips add them to a model without them. The import recognizes them and logs them, but doesn't import them: PGSuper's intermediate diaphragms come from the girder library entries' diaphragm layout rules, so importing them means fitting rules to the model's diaphragms.

### Phase 8: Visual Studio unit tests (one of the last steps)
The checks that don't need a running BridgeLink move to Visual Studio unit tests, with a static library of the testable core: see [UnitTestPlan.md](UnitTestPlan.md).

## Working conventions

The working rules are in `CLAUDE.md`. In short: commit code and docs separately from `Tests/ImportValidation/results`, and don't push (the user pushes after the PGSuper regression tests); the importer starts from `IfcImportTemplate.pgt` on purpose (template-only workarounds go with it, e.g. `CIfcBridgeImporter::RemoveUnusedTemplateEntries`); geometry for the export comes from WBFL GenericBridge tools through PGSuper interfaces, and new interface methods go at the end of the interface; entities outside the Alignment-based View (AbV), such as `IfcRelConnectsWithRealizingElements`, are out of scope; IfcOpenShell's upstream branch is `v0.9.0`, and new IfcOpenShell files and commits need its AI-assistance note (its AGENTS.md).
