# IFC → PGSuper Model Import Improvement Plan

Drafted 2026-09-23. Last completed work before this plan: estimating deck slab edges from geometry.

## Resume here (status 2026-09-25)

Work was paused here to resume in a new session, possibly on a different computer. Everything below is in the repositories; nothing depends on local notes.

**Repositories** (all on `develop`, committed and pushed as of 2026-09-25):
- PGSuperIfcExtensions: importer/exporter work through bearing data (see Phases 0-3 and 2 below).
- WBFL: `ISectionCutTool::CreateSlabShapeWithoutHaunches`/`CreateHaunchShape`, `GetExteriorGirderPoint` (uses the exterior girder whose line the cut crosses), the `CreateSlabShape` refactor (`CreateSlabTop`), and the fillet point elevation fix. PGSuper regression tests passed (minor numerical differences only).
- PGSuper: `IShapes::GetSlabShapeWithoutHaunches`/`GetHaunchShape` (added at the end of the interface), and the Test/TxDOT agent command line fix (F2).
- IfcOpenShell (`F:\IfcOpenShell`, branch `rab_infrastructure`): the false GEO 242 "IfcCurveSegment belongs to multiple IfcCompositeCurve instances" warning is fixed (592ea02f7, #9590), with the example `src/examples/IfcCompositeCurveSegments.cpp`. `v0.9.0` is the upstream branch that matters (compare and target PRs against it). The IfcOpenShell build must be current, or the fix and example don't take effect.

**Decisions waiting on answers:**
- Phase 2, boundary conditions and connections (see Phase 2): (1) is `Pset_BearingCommon` in the AbV? (2) should the exporter emit PGSuper's pier diaphragms as `IfcBeam .DIAPHRAGM.`? (3) default when there is no fixity information? (4) is "after deck" the right default for continuity timing?
- Phase 5, girder library entries (see Phase 5): (1) a new PGSuper interface to create girder library entries (creation is UI-only today)? (2) clone the closest entry of the family as the base for data the IFC doesn't define? (3) naming of created entries? (4) prefer an existing entry that matches geometrically over creating one, even when names differ?
- F1 (girders as constructed vs. with camber) is still open; the haunch bottom follows PGSuper's haunch depth, so a gap to the cambered girder is expected for now.

**Parked / to do:**
- IfcOpenShell drops the main deck slab solid of the `PGSuper_Skew_Straight_AlongPier` export without a message (see Phase 2). Needs an IfcOpenShell debug build to find the failing OpenCASCADE loft step.
- `kernels/opencascade/loft.cpp:176-180` on IfcOpenShell `rab_infrastructure` prints every loft's sections to `std::wcout` (leftover debug output).
- Phase 8: Visual Studio unit tests for the checks now run by hand (one of the last steps).
- E2 validation template with unusual values; expected values for PennDOT and Iowa.

**Next:** Phase 5 (girder library matching and entry creation) once its questions are answered; Phase 2 boundary conditions once those are answered.

**How to validate:** Release build of BridgeLink with this extension (`%ARPDIR%\BridgeLink\RegFreeCOM\x64\Release`), then `python Tests/ImportValidation/run_validation.py` (see its README). Models are in `Tests/ImportValidation/models.json`: the PGSuper round trip (0 mismatches), the three skew models (`PGSuper-Skew`, `-Bearings` with 0 mismatches, `-AlongPier` with the parked deck issue), PennDOT, and Iowa. Use `--config-file` with a copy of the template that has unusual values to show a value is set by the import rather than a template default.

**Working conventions:**
- Commit code and docs separately from `Tests/ImportValidation/results` (code first, then an "Updates IFC import validation results ..." commit). Don't push; the user pushes after PGSuper regression tests. Don't commit WBFL or PGSuper changes until the user says so.
- The importer starts from `IfcImportTemplate.pgt` on purpose (fault tolerant, handles data drift); building `CBridgeDescription2` from scratch is the long-term goal. Keep import code independent of template values.
- PGSuper configuration may be switched for tests with `/Configuration`; record the current one and switch back (normally `Regression:Regression`).
- PGSuper system units are SI (kg, m, Pa), in memory and in `.pgs` files.
- Geometry for the export comes from WBFL GenericBridge tools, exposed through PGSuper interfaces that don't reveal the implementation. New interface methods go at the end of the interface.
- `IfcRelConnectsWithRealizingElements` and other entities outside the Alignment-based View (AbV) are out of scope.
- New files in IfcOpenShell need the "generated with the assistance of an AI coding tool" note, and commits that change code need it in the message body (see its AGENTS.md).

## Requirements

### R1: Clean round trip - GlobalIds don't change
Importing an IFC model into PGSuper and exporting it again must give the elements the same `GlobalId` they had in the imported model. Downstream tools (model federation, change tracking, BCF issues, IDS pinned by GlobalId) identify elements by GlobalId, so a round trip through PGSuper must not look like new elements.

**What it takes:**
- **Store GlobalIds at import.** Record the GlobalId of every IFC element that becomes PGSuper data, in the PGSuper project (`CIfcExtensionAgent` already persists project data through `IAgentPersist`). Key it by PGSuper's stable IDs, not indices, so edits that insert or remove spans or girders don't shift the association:
  - `IfcProject`, `IfcSite`, `IfcBridge`, `IfcBridgePart`s (superstructure, substructure, deck, per-pier parts): by role and pier ID
  - `IfcAlignment` and its components
  - `IfcBeam` girders: by segment ID
  - `IfcSlab` deck
  - `IfcBearing`: pier ID + face + girder ID
  - `IfcWall` barriers: by side
  - `IfcReferent`s, and positioning referents for piers: pier ID
- **Reuse them at export.** The exporter creates GlobalIds in about 100 places (`ifcopenshell::global_id()`: 50 in `IfcExporter.cpp`, 26 in `PropertySets.h`, 13 in `Referents.h`, the rest in `QuantitySets.h`, `RebarRelationshipBatch.h`, `USBridge_Classifications.h`, `Materials.h`). Route them all through one GlobalId provider. Wherever the exporter creates an element with `ifcopenshell::global_id()`, it first looks up a stored GlobalId for the corresponding PGSuper item. Elements without a stored GlobalId get a new one. Stored GlobalIds whose PGSuper item was deleted are dropped.
- **Secondary objects** (property sets, relationships, types, materials) get new GlobalIds on each export. Only element GlobalIds need to be stable (B2).
- **Unmodeled elements** (wingwalls, piles, footings, approach slabs, rebar from other tools) are carried through: keep the source model and merge the PGSuper export into it (B1). The user can also keep the original IFC unmodified for safe keeping.
- **Mapping tables** (see "Mapping tables and IDS"): the target/element identity used for GlobalIds should be the same element identity the tables use, so both are built on one element registry.

**Validation:** add an IFC → PGSuper → IFC check to `Tests/ImportValidation`. For each IFC model, import it, export the imported project (`/IfcExport`), and compare the GlobalIds of corresponding elements (matched by entity type and role/position). The report lists changed, missing and new GlobalIds, including carried-through elements. For PGSuper models, export twice and check that the element GlobalIds are identical between the two exports.

### R2: Import issues as BCF for the modeler (optional)
When the import finds missing data or has to make an assumption, the user can have it write a BCF file (BIM Collaboration Format) with one topic per issue. The modeler opens it in their authoring or review tool, and each topic selects the affected elements by IFC GlobalId and zooms to them. The issues go back to the people who can fix the model, instead of staying in the PGSuper import log.

**Issues the importer already detects** (today they are only log lines):

| Issue | Affected elements | What the import did |
|---|---|---|
| No alignment in the model | `IfcBridge` / deck | Derived a straight alignment from the deck |
| Pier station not given (`IfcReferent` / `Pset_Stationing`) | `IfcBridgePart` pier | Estimated it from the bearing geometry, or assumed 100 ft spacing |
| `NumberOfSpans` missing | `IfcBridge` | Derived it from the piers |
| Girder designation disagrees with its location | `IfcBeam` | Used the location |
| Girder type not in the library | `IfcBeam` / `IfcBeamType` | Used the most similar or a default library girder |
| f'c / f'ci not found | `IfcBeam` / material | Kept the default |
| Span with fewer than 2 girders | span girders | Kept the default girder count |
| Deck dimension not measurable (gross/edge depth, taper) | `IfcSlab` | Kept the default |
| Haunch not part of the deck solid | `IfcSlab` | Kept the default haunch shape and fillet |
| Deck doesn't reach the bearings | `IfcSlab`, `IfcBeam` | Extrapolated the slab offset |
| Unit not defined in the model | property owner | Assumed SI |
| Geometry that couldn't be processed | the element | Skipped |

**Design:**
- **Structured import report first** (Phase 0.5). Replace ad-hoc `Logger::Info` calls with an issue collector. Each issue has:
  - a category (missing data, assumption, conflict, geometry)
  - a severity (error, warning, info)
  - a title and description
  - the affected IFC elements (GlobalIds)
  - the PGSuper item
  - the value used and where it came from (model, geometry, derived, default)
  
  The collector feeds the results dialog, the import log, a JSON report for `Tests/ImportValidation`, and the BCF writer.
- **BCF writer:**
  - one topic per issue kind, with all affected elements selected in its viewpoint, rather than one topic per element (21 girders missing f'ci is one issue, not 21)
  - the description lists each element and the value used
  - topic type and priority come from the severity
  - labels: `PGSuper import` plus the category
  - one viewpoint per topic: the components select the elements by `IfcGuid`, plus an orthogonal or perspective camera fitted to the elements' bounding box (the importer already has their meshes)
  - no snapshots at first
- **Options:**
  - `CIfcImportOptions::bcf_file` (empty = no BCF, the default; C5) and the BCF version (2.1 or 3.0; C1)
  - an import options checkbox, with the path defaulting to `<model>.bcf`
  - command line `/IfcBcf=<file.bcf>`
- **Mapping tables:** once tables exist, a "missing data" topic can say where the active table (and the agency IDS it came from) expects the value, e.g. "f'ci expected in `IaDOT_PPCB`.`6_Concrete Release Strength, Fci`".
- **Validation:** `run_validation.py` writes the BCF for each IFC model and checks it: zip structure, XML against the BCF schemas, one topic per expected issue kind, and every referenced GlobalId present in the model.

**Decided:** C1–C5 (see Decisions). Test with Blender/Bonsai (C4).

## Decisions (2026-09-24)
Answers to the open questions (`IFC_Import_Open_Questions.docx`). The IDs are used in the rest of this plan.

**Mapping tables and IDS**
- **A1 Location:** The BridgeLink configuration system (catalog servers) does not extend to extension agents. The IFC extension adds a page to the BridgeLink configuration wizard where the user picks the mapping file. The file can be anywhere, and the setting is stored in the registry. `/IfcMapping=<file>` overrides it on the command line.
- **A2 Selection:** The mapping table comes from that registry setting (or `/IfcMapping=`). No automatic detection.
- **A3 Format:** JSON.
- **A4 IDS binding:** A separate binding (overlay) file. Agency IDS files are not modified.
- **A5 Exported content:** Property sets, quantity sets, classifications, and element names/ObjectType are table driven. Geometry and the spatial structure stay in code.

**Clean round trip (R1)**
- **B1 Unmodeled elements:** Carry them through. Keep the source model and merge the PGSuper export into it. The user must also be able to keep the original IFC unmodified for safe keeping.
- **B2 Secondary objects:** Only element GlobalIds need to be stable. No derived GlobalIds for property sets, relationships, types, or materials.

**BCF (R2)**
- **C1 Version:** Both 2.1 and 3.0, selectable.
- **C2 Topics:** One topic per issue kind, selecting all affected elements.
- **C3 Zip:** A vcpkg library (libzip or minizip-ng). PGSuper's MakePgz already zips with a bundled copy of Lucian Wischik's `zip.cpp`/`unzip.cpp` (repackaged zlib). Aim for one standard, consistent zip solution across the projects.
- **C4 Tools:** Probably Blender/Bonsai. Test the BCF there.
- **C5 Default:** Off. Turned on in the import options.

**Import behavior**
- **D1 Conflicts:** Geometry wins for layout (stations, spacing, dimensions), property sets win for materials, and every conflict is reported.
- **D2 Slab offset tolerance:** 1/4 in. Make it a configurable value: an import option at first, user interface later.
- **D3 Girders not in the library:** Interactive: prompt to pick a similar entry or create one. Command line: create the entry and log it.

**Validation**
- **E1 Drawings:** Not available. Expected values for PennDOT and Iowa come from the model property sets only.
- **E2 Validation template:** Yes. Create a template with deliberately unusual values in PGSuper and add it to `Tests/ImportValidation`.

**Exporter and PGSuper**
- **F1 Deck haunch vs girder camber:** Leave it for now. The underlying question is whether girders should be modeled as constructed (no camber) or in the expected final condition (with camber). It will be discussed with another engineer and may become an export option.
- **F2 TxDOT command line:** Change the TxDOT agent (PGSuper repository) so it only claims its own flags. **Done 2026-09-24** (PGSuper `develop`, 1c45952285): it was the Test agent (`/TxA...`) that rejected file-first command lines; the TxDOT agent (`/TxTOGA`) had the same flaw. Both now claim only their own commands.
- **F3 Exporter placeholders:** Fix the stress unit name now. The barrier properties later.

## Test models

| | PGSuper export | PennDOT Rearick Rd | Iowa US59 |
|---|---|---|---|
| Path | `Desktop\PGSuper_Import_Model.ifc` | `BIM for Bridges\Precast Concrete Girder Modeling\PennDOT Rearick Rd\323010451_STR1.ifc` | `BIM for Bridges\Precast Concrete Girder Modeling\Iowa ADCMS\Iowa ADCMS Pilot_US59 over IA92_2025_07_09.ifc` |
| Authoring tool | PGSuper 9.1 | OpenBridge (Quadri) | unknown, 65 MB |
| Alignment | yes, plus `IfcReferent` stationing on each pier | **none** | **none** (has `IfcMapConversion`) |
| Spans / girders | 3 spans, 4/5/4 girders | 1 span, 5 girders, integral abutments | 3 spans × 7 girders, 3-column piers with caps |
| Beam geometry | `IfcSectionedSolidHorizontal` | `IfcFacetedBrep` | mapped `IfcTriangulatedFaceSet` |
| Beam type name | `IfcBeamType` "WF66G" | pset `_PS Concrete Beams.Type` = "PA Bulb-Tee 33/31.25" | pset `IaDOT_PPCB.2_Type` = "BTB45" and classification |
| Beam names | "Span 1, Girder 1" | "Beam2" (no span) | "Beam 7 Span 3" |
| f'c / f'ci | standard psets, per-property units (ksi) | custom pset, string values ("8", "6.8") | custom pset, unitless numbers in ksi |
| Deck | 1 FLOOR slab | 1 FLOOR + approach slabs | 1 FLOOR + 21 haunch slabs (`IfcSlab.USERDEFINED` "Haunch") |
| Other | 369 strand `IfcTendon`, 2 barrier `IfcWall` | parapets, railings, diaphragm | 661 STRAND tendons, diaphragms, pier caps/columns |

## Baseline (2026-09-23)

The IfcOpenShell 0.9 crash is fixed:
- The `ifcopenshell_geometry_mapping_ifc4x3_add2.dll` plugin was not deployed, so the geometry iterator threw.
- `ImportFromIFC` now catches `std::exception`.

The Iowa model imports end to end. What its log shows:

- **Works:** girder keys are parsed from names ("Beam 7 Span 3"); the span count is derived from the number of piers; pier stations are estimated from bearing geometry; the deck edge step runs without errors.
- **Nothing mapped for girder type or concrete strengths.** There's no `IfcBeamType`, so every beam gets "Unknown" and falls back to the first library I-beam entry (`Unknown_Girder_Type`). f'c/f'ci are in `IaDOT_PPCB` and the girder type is in `IaDOT_PPCB.2_Type`/classification ("BTB45"), but none of these are read.
- **`GetBeamTypeCount()` returns 0 when there are no `IfcBeamType`s**, so `UseSameGirderForEntireBridge` becomes false and the library lookup runs once per beam. Also, the same-girder branch uses `*beams.begin()`, which is *all* `IfcBeam`s (pier caps and diaphragms included), not `prestressed_beams`.
- **Wrong or noisy log messages:**
  - "Pset_PrecastConcreteElementGeneral not found" is logged 3× per beam, because the girder key is computed three times.
  - `GetGirderLibraryEntry` logs each message twice.
  - The f'c "not found" message names the wrong pset (it should be `Pset_MaterialConcrete`).
  - The pier-station message says "IfcBridgePart.PIER" for abutments.
  - The derived alignment is never mentioned.
- **Slow geometry.** Beam geometry takes ~1.5 s per beam (21 beams ≈ 35 s) and bearings ~1 s per pier. Each product's geometry is run through the iterator separately, and the deck is processed twice (once to derive the alignment, once for the edges). Caching meshes by product id and running the iterator once over all needed products would help.
- **Not yet checked:** the imported PGSuper model's layout (stations, spacing, deck edges) against the drawings. That's Phase 0's job.

## Known bugs and fragile spots

1. ~~`GetProperty` ignores `IfcPropertySingleValue.Unit`~~ **Fixed 2026-09-23:** `CIfcImportUnits` (`IfcImportUnits.h/.cpp`) reduces any IFC unit (SI with prefix, conversion based, derived) to a factor to fundamental SI units and converts with `WBFL::Units::ConvertToSysUnits`. Values without their own unit use the project unit. f'c, f'ci, pier and alignment stations, and referent distance-along are converted. On the exporter side, `AddUsedProjectUnits` declares SI project units for unit types used by values that have no unit of their own.
2. ~~Span count in `get_beam_spacing`~~ **Fixed:** girders are assigned to spans by pier stations (see below).
3. ~~Assumed girder keys computed three ways~~ **Fixed:** `get_girder_layout` is the single source of girder keys.
4. ~~Girder start/end by distance from the origin~~ **Fixed:** chosen by station.
5. ~~Spacing measured as 3D chord, stored as AlongItem/AtPierLine~~ **Fixed:** measured normal to the alignment at the pier line (NormalToItem/AtPierLine).
6. ~~Crash with no ABUTMENT parts~~ **Fixed.**
7. Girder orientation is hard-coded to the first supported option.
8. The mesh → decompose → top surface → boundary → segments pipeline is copied three times.

## Phases

### Phase 0: Test setup
- **Headless import (done 2026-09-23):** `BridgeLink.exe /IfcImport=<model.ifc> <template.pgt> [/IfcOut=<project.pgs>] [/IfcLog=<import.log>]`.
  - `CIfcExtensionAgent` implements `IEAFProcessCommandLine`. It imports into the template-created project, saves it, and EAF exits.
  - `/IfcImport` had to come first because PGSuper's Test agent claimed any command line whose first parameter was not a flag. Fixed in PGSuper (F2); putting the flag first still works with older PGSuper builds.
  - `CIfcImportOptions::interactive = false` suppresses the results dialog and the bridge picker.
  - First round-trip results on the PGSuper model: pier stations match; f'ci comes in as ≈4 Pa (the units bug); deck edges come in as 199 raw points vs 1 in the source; the WF66G/WF100G girder types are not in the `IfcImportTemplate.pgt` library.
- **Headless export:** `BridgeLink.exe /IfcExport=<model.ifc> <project.pgs> [/IfcPropertyUnits=Display|System] [/IfcLog=<export.log>]` exports with the default `CIfcExportOptions`. The round trip (.pgs → IFC → .pgs) runs without the UI in both property-unit modes.
- **Validation tooling (done 2026-09-23):** `Tests/ImportValidation/` has `pgs_extract.py`, `compare_bridge.py`, `run_validation.py` and `models.json`; see its README.
  - Round trip of the PGSuper model: 66 of 168 values set by the importer match under the Regression configuration, and 79 under WSDOT (Local), where the girder types match.
  - Still open: slab offset (template default), and spacing written along the pier instead of normal to it.
  - **Deck edges fixed 2026-09-23:** the left and right edges are kept as separate station/offset profiles, each simplified (Douglas-Peucker, 1/8 in), with parallel transitions where the offset doesn't change. The PGSuper, PennDOT and Iowa decks each come in as one parallel point (was 199 for PGSuper). This replaces `condense()` and its symmetric-deck fallback, and `LeftEdge` gets the correct sign when both edges are on the same side of the alignment.
  - Follow-ups:
    - A validation template with deliberately unusual values, created in PGSuper (E2), so fewer values are ambiguous "default matches" (e.g. the 15 ft deck edges and the end distances in `IfcImportTemplate.pgt` coincide with the PGSuper model).
    - Evaluate alignments with spirals: export the alignment to IFC and evaluate it with IfcOpenShell.
    - `expected/PennDOT.json` and `expected/Iowa.json` from the model property sets (no drawings; E1).
- Python ifcopenshell "expected values" scripts for the PennDOT and Iowa models.
- Structured import report: parameter | value | source (pset/geometry/assumed/default) | confidence. This issue collector is also the source for the BCF file (R2).

### Phase 1: Correctness of what we already import
- Fix bugs 1–6, plus the baseline findings: the beam-type count/`beams.begin()` bug, compute each girder key only once, and fix the log messages.
- Single geometry pass: iterate once over beams, bearings and the slab, and cache meshes by product id. Consider `num_threads > 1`, but the progress stream isn't thread-safe. Honour per-property units, then project units. Add PRESSUREUNIT, MASSDENSITY and FORCE.
- **Done 2026-09-23:** `get_girder_layout` (BeamSpacing.cpp) runs one geometry pass over the girders, assigns each to the span containing its mid-point, orders girders left to right by offset, and cross-checks model designations (a mismatch is logged, the location wins). Girder counts, spacing and girder properties all use it. `UseSameGirderForEntireBridge` is set when every girder resolves to the same library entry (was: count of IfcBeamType). Round trip: 85 of 168 match, 0 mismatches, spacing matches exactly. PennDOT spacing is 6.50 ft (was 6.00 ft; the girder axes are 1.98 m apart).
- One place that owns the girder layout: get each beam's CL end points from geometry → station/offset → span by pier stations → order by offset → cross-check against parsed names.
- Pull the geometry pipeline out into one reusable helper.

### Phase 2: Superstructure geometry
- **Done 2026-09-24 (pier skew):** `set_pier_orientation` (Piers.cpp) takes the CL pier direction from the first of: an explicit `RefDirection` of the pier's positioning `IfcReferent` (the exporter now writes the local x-axis along the CL pier, z up), the CL bearing lines, the line through the girder ends, or normal to the alignment (logged as assumed). The other sources are cross checked (0.5 deg). The skew is from the left normal of the alignment (PGSuper's NORMAL pier line), positive counterclockwise, and is written as "dd.ddd L|R". End distance and bearing offset are then measured normal to that CL pier. Without bearings, the template end distance is kept and the bearing offset locates the girder ends (girder ends match the model for the design-then-update workflow); the implied CL bearings are used for the slab offset.
  - Test models: `PGSuper_Skew_Straight.pgs` (no bearing dimensions, so no `IfcBearing` export) and `PGSuper_Skew_Straight_Bearings.pgs`: straight N 30 W alignment, 4 spans, orientations 20 L, NORMAL, 15 R, 30 R, N 35 E. `pgs_extract` compares `pier[i].skew` (angle) instead of the orientation text, so a bearing and the equivalent skew compare equal.
  - Results: all 5 skews in both models. With bearings, all 20 end distances and bearing offsets match. Without bearings, the girder ends match on all 10 faces. (The no-bearings model first had Pier 2 back with an end distance larger than the bearing offset, a girder end past the CL pier, which PGSuper corrected silently; it was changed to match the bearings model.)
- **Done 2026-09-24 (end distance and bearing offset):** `locate_bearings` (Piers.cpp) finds the CL bearing at each girder end as the average of the `IfcBearing` centers (bounding box) near the CL girder at that end, projected onto the CL girder. `set_end_distance_and_bearing_offset` fits the CL bearing line through the CL bearings on each pier face; its direction also gives the pier line, so skewed piers work before skew is imported. End distance (FromBearingNormalToPier) and bearing offset (NormalToPier) are the medians over the girders, and any spread over 1/4 in is logged. An abutment's unused face gets the same values. The slab offset is now measured at these CL bearings.
  - Results: the PGSuper round trip matches all 16 values; with a scratch template (12 in / 24 in), all 16 are recovered (20.5 in / 32.5 in), so they aren't default matches. PennDOT: 12.0 in end distance, 0 bearing offset. Iowa: 8.0 in end distance, 0 bearing offset at the abutments, 12.0 in at the piers. (PennDOT and Iowa pier stations are estimated from the bearings, so the abutment bearing offset is 0 by construction.)
- Pier/abutment orientation (skew) from bearing centroids or girder end segments → `CPierData2::SetOrientation`.
- **Done 2026-09-24 (uniform spacing):** when every span has girders and all spacings (normal to the alignment at the pier lines) are within 1/16 in, the bridge gets `sbsUniform` with the average spacing (NormalToItem, AtPierLine); otherwise spacing is general, per pier face. Skew with bearings model: 0 mismatches. PennDOT 6.500 ft, Iowa 6.694 ft uniform. A bridge that is uniform only along skewed pier lines imports as general (still correct geometry).
- **Done 2026-09-24 (uniform spacing along the CL piers):** spacing is also measured along each CL pier (`get_beam_spacing` with the pier directions, which are now found before the spacing). Uniform normal to the alignment is tried first, then uniform along the CL piers (AlongItem). Test model `PGSuper_Skew_Straight_AlongPier.pgs` (the bearings model with the bridge spacing measured along the CL pier, so the girders fan between skewed piers): uniform 6.000 ft along the CL piers is recovered. This model also found a bug in WBFL `GetExteriorGirderPoint` (it used the first exterior girder, span 1, for the whole bridge), fixed.
  - **To do (parked):** IfcOpenShell drops the main deck slab solid of the AlongPier export without a message and keeps only the haunches, so its deck imports as a haunch (64 mismatches). The deck sections have constant point counts, no self intersections, and stay inside the deck edges. IfcConvert fails even for a single pair of consecutive sections; identical or symmetric soffits work, and moving the soffit points of the working bearings model still works. The failure is in the OpenCASCADE loft conversion (`kernels/opencascade/loft.cpp`), which returns false silently; a debug build is needed to see which step fails.
  - **To do:** `loft.cpp:176-180` on the IfcOpenShell `rab_infrastructure` branch prints the first two loft sections to `std::wcout` on every loft (leftover debug output).
- Girder end distance and bearing offset from beam end, bearing centroid and pier line.
- **Done 2026-09-24 (bearing data):** each `IfcBearing` belongs to the girder end whose CL is nearest to it. `set_bearing_data` (Piers.cpp) sets `CBearingData2` (basic definition) from the first bearing at each girder end: length along and width across the girder, height, round when there are 8 or more distinct plan points at the same distance from the center, and the number and spacing of the bearings. Fixity from `Pset_BearingCommon.DisplacementAccommodated`, or from a text property named `...Fixity` ("Fixed", "Expansion"/"Free"; PennDOT, Iowa) until the mapping tables exist. The bearing is set for the bridge, each pier face, or each girder, depending on how much it varies. Detailed bearing data (elastomer layers) isn't in the IFC and is not imported. Results: the PGSuper round trip and the skew bearings model match all bearing values; Iowa 9.0 x 28.5 x 1.0 in at the piers and 3.0 x 28.5 x 3.0 in at the abutments (matching its `Bearing Item - ItemTypes` properties), fixed; PennDOT 12.0 x 31.5 x 0.75 in, fixed.
- **Boundary conditions and connections (proposed, pending decisions):** `IfcRelConnectsWithRealizingElements` is not in the AbV, and properties should not be required of modelers, so the boundary condition at each support is inferred from the modeled elements and every decision is reported with its evidence (log, and BCF per R2):
  - Evidence per support: (A) deck continuity (separate slabs / DECK_SEGMENT parts ending at the support, a gap in the deck, an expansion joint element); (B) a diaphragm (IfcBeam .DIAPHRAGM., or a wall/member typed or named diaphragm, e.g. PennDOT `Abut#_EndDiaphragm`) enclosing the girder ends of both spans; (C) the diaphragm or girder ends touching or inside the substructure (pier cap, abutment wall, footing) with no bearing gap; (D) bearings present/absent and their fixity; (E) explicit hints (properties, classifications, vendor text like PennDOT `End - Integral`) through the mapping tables, used when present and cross checked.
  - Rules: deck discontinuous → simple (hinge/roller from fixity); else a diaphragm across the support holding both girder ends → integral if it touches the substructure or there are no bearings, otherwise continuous; else simple spans under a continuous deck (weak, reported). Abutments: embedded/touching → integral, else hinge/roller. Before/after deck: no reliable geometric signal, default after deck (logged).
  - PGSuper round trip: the exporter would emit PGSuper's pier diaphragms as `IfcBeam .DIAPHRAGM.` (one across continuous/integral piers, one per span at hinge/roller piers; extended to the support for integral) and bearing fixity; the deck stays one slab (it is continuous in PGSuper even over simple spans).
  - Real models: Iowa has pier and abutment `IfcBeam .DIAPHRAGM.`, pier caps, one deck, all bearings "Fixed"; PennDOT has abutment end diaphragms as `IfcWall`, one deck, "Fixed" bearings.
  - Open questions: (1) is `Pset_BearingCommon` in the AbV? (2) should the exporter emit the diaphragms as above? (3) default when there is no fixity information (e.g. hinge at the first abutment, rollers elsewhere)? (4) is "after deck" the right default?

### Phase 3: Deck
- **Done 2026-09-24 (cast-in-place decks):** `analyze_deck_section` (DeckSlab.cpp) samples the deck and girder geometry along lines normal to the alignment at the quarter points of each span and takes the median.
  - **Gross depth:** deck thickness 2 in outside the interior flange tips (clear of fillets and crown points).
  - **Overhang edge depth:** thickness just inside each deck edge.
  - **Overhang taper:** the overhang soffit is fitted and extrapolated to the exterior flange tip. It's classified per `SectionCutTool`: about equal to the edge depth → none; meets the haunch bottom (or the girder top when the haunch isn't in the deck solid) → to top of top flange; deeper → to bottom of top flange.
  - **Haunch shape and fillet:** only when the deck solid contains the haunches. The width of the transition from gross depth to the haunch bottom gives square or filleted, and the Fillet.
  - **Slab offset ("A"):** top of deck minus top of girder above each CL bearing, which is located with the girder end distance. Set for the bridge, per bearing line, or per segment depending on the spread. The tolerance was 1/16 in; per D2 it becomes 1/4 in and an import option (user interface later). Where the deck doesn't reach the bearing (integral abutments: PennDOT, Iowa), it's extrapolated linearly from the first 3 ft where it does.
  - **Edge depth:** where the soffit within 3 in of the edge is straight (fit residual ≤ 1/32 in), the fitted soffit is extrapolated to the edge; otherwise (chamfers, drips, e.g. PennDOT) the depth 1/4 in inside the edge is used. Deck edge stations and offsets are rounded to 1 µm to remove geometry noise.
  - **Fixed 2026-09-24 (WBFL `SectionCutTool::CreateSlabShape`):** the fillet end points (1 and 6) took the deck elevation above the flange tip instead of their own offset, so on a crowned deck the bay soffit wasn't parallel to the deck top (7.5 ± 0.015 in with a 2% crown and 0.75 in fillet). Deck depths are vertical in the section plane.
  - Results: the PGSuper round trip recovers every value exactly (gross 7.500 in, edges 7.000 in, taper, filleted with 0.75 in, A = 9.500 in). PennDOT: 8.00 in gross, 10.02 in edges, square haunch, A about 10.0 in per segment. Iowa: 8.47 in gross, about 10 in edges, A 9.0 to 10.7 in per segment.
- **Done 2026-09-24 (haunch modeling):** the deck and haunches are cast together, so they are one `IfcSlab` in the export, and the importer measures all of the deck concrete however it is modeled.
  - Importer: `get_deck_concrete_mesh` (DeckSlab.cpp) adds separate haunch elements (parts aggregated under the deck, and `IfcSlab`/`IfcBuildingElementPart` with "haunch" in ObjectType or Name; Iowa has 21) to the deck mesh for the section measurements. Haunches in the deck solid (PennDOT) or as other representation items (PGSuper export) are already in the deck mesh. The slab offset is only measured where the concrete is at least the gross depth, so a haunch beyond the end of the slab isn't taken as the deck.
  - Exporter: the deck representation has the slab without haunches, swept along the alignment, plus a haunch solid over each mating surface of each segment, swept along the segment and clipped to the slab stations. The geometry comes from new WBFL `ISectionCutTool::CreateSlabShapeWithoutHaunches`/`CreateHaunchShape` (the slab shape has the same points at every station, with exterior girder lines extended at skewed ends), exposed as PGSuper `IShapes::GetSlabShapeWithoutHaunches`/`GetHaunchShape`. The haunch bottom is PGSuper's haunch depth (the gap to the cambered girder is left for F1). The skewed-deck limitation (no haunches) is gone.
  - Results: PGSuper round trip unchanged (0 mismatches); the skew models now recover gross depth, edge depths, taper, and haunch shape.
- Still to do:
  - stay-in-place deck panels
  - the CL bearing location depends on end distances, which are still template defaults (Phase 2)
  - deck concrete f'c
- **Exporter finding:** at mid-span of the PGSuper export, the deck solid's haunch bottom is about 1.6 in above the exported girder top (the deck haunch is about 0.4 in thick, A = 10.0 in measured from the girder). The deck uses PGSuper's computed haunch, while the girder geometry includes camber. At the bearings they agree. Left as is for now (F1): whether girders are modeled as constructed or in the final condition is still to be discussed, and may become an export option.

### Phase 4: Alignment/profile when the model has none
- Derive the vertical profile and crown slope from the deck top surface.
- Use deck station psets (Iowa) for the reference station.

### Phase 5: Girder section and material mapping
- **Current behavior (interim, 2026-09-23):** `GetGirderLibraryEntry` matches names only: `IfcBeamType.Name`, `ObjectType`, "...Type"/"...Shape..." text properties and non-usBridge classification names, exact (ignoring case and punctuation) or most similar (longest common substring >= 0.5), otherwise the first I-Beam. Each girder type is matched and logged once.
- **What the models carry (2026-09-25):**
  - PGSuper export: `IfcBeamType` named for the library entry (WF100G, WF66G) and `usBrPset_PrecastConcreteBeam.ShapeName`; each girder is an `IfcSectionedSolidHorizontal` of `IfcArbitraryClosedProfileDef` polygons (the exact PGSuper section, including end blocks).
  - PennDOT: no `IfcBeamType`; `IfcFacetedBrep`; `_PS Concrete Beams` properties (`Type`, `Beam Designation`, f'c, f'ci).
  - Iowa: no `IfcBeamType`; `IfcMappedItem` of swept `IfcArbitraryClosedProfileDef` sections; `IaDOT_PPCB.2_Type` (e.g. BTB45) and classification "Beam, PPC, BTB45".
  - PGSuper: about 15 girder families (I-beam, bulb tee, NU, U, box, voided slab, decked bulb tees, double tee, multi-web, ...), each a beam factory with named dimensions (I-beam: C1, D1-D6, H, T1, T2, W1-W5). Library entries are created only through the UI today; `LibraryFw` has `NewEntry`/`CloneEntry`/`AddEntry` and `GirderLibraryEntry` has `SetBeamFactory`/`SetDimension(name, value, bAdjustStrands)`, but there is no programmatic path that fills in a complete, valid entry.
- **Steps:**
  1. **Section from the model:** a cross section polygon (outer boundary and voids) for each girder in girder section coordinates (top center, plumb): from the swept profiles when the representation has them (PGSuper export, Iowa), otherwise by cutting the girder mesh normal to its axis (PennDOT). Cut at several stations so end blocks and tapers are seen, and use the prismatic middle portion for the section.
  2. **Match an existing library entry by geometry, not just name:** compare the section with the section of each library entry in the configuration (from its beam factory and dimensions): depth, flange widths and thicknesses, area, moment of inertia, centroid, within tolerances. A name match is confirmed by geometry; a geometric match with a different name is used and logged. Replaces the "first I-Beam" fallback. Conflicts reported (D1).
  3. **Classify the family** when nothing matches: from the section topology (number of webs, voids, open top, top flange width relative to spacing for decked sections) and the families available in the configuration.
  4. **Fit the family dimensions:** initial values from measured features (depth, flange widths, thicknesses, web width, fillets/chamfers), refined by minimizing the difference between the measured polygon and the factory's section for the dimensions (area of the symmetric difference). Accept only if area, I, and centroid are within tolerance; otherwise report that the section doesn't fit the family.
  5. **Create the library entry (D3):** a new programmatic path in PGSuper (a new interface, since creation is UI-only today): clone the closest existing entry of the family (keeps strand grids, shear, debonding, and other data the IFC doesn't define), set the fitted dimensions with `bAdjustStrands`, name it from the IFC type name (unique if it collides), add it to the project library, validate it, and log it. Command line: create and log. Interactive: prompt to pick one of the geometric candidates or create the entry (later step).
  6. **Validation:** the PGSuper round trip under the Regression configuration (WF girders not in the library) must create entries whose dimensions match the source project's entries; `pgs_extract` compares the girder entry dimensions. PennDOT and Iowa: fitted dimensions compared with their standard shapes (PA bulb tee, Iowa BTB45) where they are known.
- **Open questions:** (1) a new PGSuper interface to create girder library entries (e.g. clone an entry and set dimensions)? (2) clone the closest entry of the family as the base for data the IFC doesn't define? (3) name for created entries (IFC type name, with a suffix when it collides)? (4) prefer an existing entry that matches geometrically over creating one, even when the names differ?
- **Configuration for exact-match testing:** `BridgeLink.exe /App=PGSuper /Configuration="WSDOT (Local)":"WSDOT (Local)"`. Record the current `CatalogServer2`/`Publisher2` (`HKCU\...\PGSuper\Options`) first and switch back with `/Configuration` when done; the dev machine is normally `"Regression":"Regression"`.
- Orientation from the measured beam roll/plumb.
- Material lookup chain: standard psets → configurable vendor pset aliases (with string/unit-in-name parsing) → default. Vendor psets and girder name aliases are handled by mapping tables (see "Mapping tables and IDS").

### Phase 6: Prestressing
- Strand `IfcTendon` → section cuts → `CStrandData` direct fill or strand rows. Harped strands and harp points from the tendon path, debonding from tendon length. Strand material, Pjack from `InitialTension`.

### Phase 7: Secondary elements and loads
- Barriers → traffic barrier library match or weight per length. Diaphragms → diaphragm loads. Pier caps and columns → physical pier model (low priority).

### Phase 8: Visual Studio unit tests (one of the last steps)
- Turn the checks run during this work into Visual Studio unit tests (Microsoft C++ unit test framework), so they run from Test Explorer and in CI instead of by hand:
  - the round trips and model imports in `Tests/ImportValidation` (`run_validation.py`), with the expected values and tolerances from `compare_bridge.py`
  - the geometry checks done by hand: deck section outlines (constant point count, no self intersection, valid for IfcOpenShell), haunch profiles, girder layout and spacing, pier skew from referents/bearings/girder ends, end distance and bearing offset
  - exports that IfcOpenShell must be able to convert (every representation item produces geometry)
  - the new WBFL `ISectionCutTool::CreateSlabShapeWithoutHaunches`/`CreateHaunchShape` and PGSuper `IShapes` methods, in their own test projects

## Mapping tables and IDS (import and export)

**Problem.** Both directions are hard coded. The importer reads a fixed set of standard psets. The exporter writes about 165 hard-coded property definitions (`PropertySets.h`, `QuantitySets.h`, `Materials.h`, `USBridge_Classifications.h`). The IDS exporter (`IdsExporter.cpp`) repeats the same pset and property names a third time. Models from other agencies put the same data elsewhere:

| Target | Standard (TPF/usBridge) | Iowa | PennDOT |
|---|---|---|---|
| Girder f'c | `Pset_MaterialConcrete.CompressiveStrength` (material, per-property unit) | `IaDOT_PPCB."5_Final Concrete Strength, Fc"` = 5.0 (no unit) | `_PS Concrete Beams."Concrete Strength at 28 days (ksi)"` = "8" (text) |
| Girder f'ci | `Pset_PrecastConcreteElementGeneral.ReleaseStrength` | `IaDOT_PPCB."6_Concrete Release Strength, Fci"` | `"Concrete Strength at Strand Release (ksi)"` = "6.8" |
| Girder type | `IfcBeamType.Name` | `IaDOT_PPCB."2_Type"` | `_PS Concrete Beams.Type` |
| Deck thickness | geometry | `IaDOT_Deck."3_Minimum Thickness"` = `8.5"` | `_Deck."Minimum Thickness (in)"` = "8" |

**Goal.** One declarative mapping between PGSuper data and IFC locations that drives import, export and IDS generation, and that can be generated from an agency's IDS. Agencies control what the exporter writes and what the importer reads without code changes.

**Two kinds of IDS - keep them apart:**
- **General-purpose (agency) IDS:** what a regular IDS is. Its specifications apply to classes of elements (every precast girder, every deck) and require properties to exist with a given data type, sometimes with a range or enumeration. It isn't tied to one model. This is the kind a mapping table is **generated from**, and the kind a table can be **written out as**.
- **PGSuper design-value IDS (`IdsExporter`):** a special IDS for one exported model. Each specification is pinned to one element (by `Name` or `GlobalId`) and requires the **exact design values** PGSuper computed for it (f'c, f'ci, camber, strand data, ...). It is used to check that a model built downstream (e.g. in a model authoring tool) still carries the design values. It is **not** a source for mapping tables. It only takes its property locations from the table, so it checks the values where the export put them.

### Concepts
- **Targets** (in code): a fixed vocabulary of PGSuper data items, e.g. `girder.fc`, `girder.fci`, `girder.type`, `girder.camber_at_release`, `deck.gross_depth`, `bearing.fixity`, `bridge.number_of_spans`. Each target has:
  - the element it belongs to (bridge, pier, girder, deck, bearing, material, ...)
  - a kind (stress, length, text, count, ...), which sets its WBFL unit type and IFC measure type
  - a getter (for export) and/or setter (for import) bound to `CBridgeDescription2` or the PGSuper interfaces
  
  Geometry, the spatial structure, and anything computed from geometry stay in code. Tables only address targets.
- **Element selectors**: which IFC entities play which role (girder = `IfcBeam.BEAM` in the superstructure, deck = `IfcSlab.FLOOR`, haunch = `IfcSlab.USERDEFINED` "Haunch", ...). These are expressed the way IDS applicability is: entity, predefined type, classification, attributes.
- **Locations**: where a target lives in IFC. That's a pset property (on the occurrence, type, or material), an attribute (`Name`, `ObjectType`), a classification reference, or a quantity. Each location has the IFC data type (`IFCPRESSUREMEASURE`, `IFCLABEL`, ...) and a unit policy:
  - IFC unit if present
  - a stated unit for unitless numbers or text (Iowa, PennDOT)
  - on export: display or system units per `display_units_for_properties`
  
  Text values may need a parser (number, feet-inches `8'-6"`, a pattern).
- **Mapping table**: selectors plus an ordered list of locations for each target, in a JSON file (A3). Tables are layered: an agency table on top of the standard table. The active table is set on an IFC page in the BridgeLink configuration wizard and stored in the registry, and `/IfcMapping=` overrides it (A1, A2).

### Bidirectional
- **Import:** for each target, try the locations in order and use the first value found. Record which table and location supplied it (feeds the structured import report).
- **Export:** write each target to its **primary** location, i.e. the first location that is marked exportable. Also write the pset and quantity-set structure the table declares, and the classifications. Locations that only make sense for reading (name parsing, patterns, text parsers) are marked import-only.
- **Design-value IDS:** `IdsExporter` takes the location of each design value (pset, property, data type) from the same table, instead of hard-coding it. The per-element pinning and the exact values stay model-specific. The IFC export, the design-value IDS and the importer then always agree on where data lives.
- The standard table must reproduce today's export exactly before any hard-coded definitions are removed. Check this by diffing exported IFC files (property sets, values, units) and with the round trip.

### Generating tables from IDS
A general-purpose (agency) IDS already states the IFC side: applicability (entity, predefined type, classification) and required properties (pset, name, data type). It does not say which PGSuper value a property holds. Proposed generation:
1. Each IDS `<specification>` applicability becomes an element selector.
2. Each `<property>` / `<attribute>` / `<classification>` requirement becomes a location.
3. The location is bound to a target by, in order:
   - a binding (overlay) file that assigns IDS facets to targets. Agency IDS files are not modified (A4)
   - the same pset/property as a location in the standard table
   - otherwise the facet is listed as unbound for a person to assign
4. The data type comes from the IDS. Units: IDS values are in SI, and the stated unit for unitless or text values comes from a human-edited overlay.
5. Output: a table file that can be reviewed and edited. Regenerating keeps the manual assignments (merge, not overwrite).

The reverse also works: a table can be written out as a general-purpose IDS (applicability by element class, required properties and data types, no design values), so an agency can publish what it expects from PGSuper models. Writing the standard table out and generating a table back from that IDS is also the test for the generator.

### Steps
- **M0:** design note, element registry (element identities shared with R1 GlobalId persistence), and target vocabulary for what import and export use today (girder f'c, f'ci, type, camber, deck, bearings, bridge geometry pset).
- **M1:** import engine plus a standard table that reproduces today's import. The validation results must not change. Add the sources to the import log.
- **M2:** Iowa and PennDOT tables (f'c, f'ci, girder type, deck thickness). The validation inventories should show these values move from template default to set by the import.
- **M3:** export engine for property sets, quantity sets, classifications, and element names/ObjectType (A5). The standard table reproduces today's export (IFC diff plus round trip). Then remove the hard-coded definitions.
- **M4:** the design-value IDS exporter takes its property locations from the table (values and per-element pinning unchanged).
- **M5:** general IDS → table generator (instructions tag, standard-table matching, unbound report) and table → general IDS writer. Test by writing the standard table out as a general IDS and generating it back, and with an agency IDS. The design-value IDS is not an input.
- **M6:** the IFC page in the BridgeLink configuration wizard (mapping file in the registry) and `/IfcMapping=` (A1, A2).
- **Validation:** `Tests/ImportValidation` gets export checks as well:
  - run ifctester with the design-value IDS against the exported model (the exported values match the design)
  - run ifctester with a general IDS written from the table (the exported structure matches the table)
  - compare exported psets against the table

### Decided
A1–A5 (see Decisions).

## Open decisions
All decided on 2026-09-24 (see Decisions). New open questions go here.
