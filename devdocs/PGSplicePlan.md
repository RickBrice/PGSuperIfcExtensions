# PGSplice Import and Export Plan

Drafted 2026-10-15, revised 2026-11-10 (export first, then import; closure joints and post-tensioning made concrete from the usBridge workbook and IDS). Reviewed with the user 2026-11-10; no open questions.

## Goal

PGSuper is a subset of PGSplice. Both use `CBridgeDescription2`. What PGSplice adds that matters for IFC:

| | PGSuper | PGSplice |
|---|---|---|
| Segments per girder | 1 | 1 or more, joined by **closure joints** |
| Girder group | one span | one or more spans |
| Supports | piers, abutments | also temporary supports (erection towers, strongbacks) |
| Prestressing | pretensioned strands | also **post-tensioning**: girder tendons (span the girder, through closure joints) and segment tendons (one segment), each in a duct |
| Construction | assumed sequence | a timeline (segment erection, closure casting, tendon stressing, ...) |

One code path written for the PGSplice model, with PGSuper as the special case (one segment, one span per group, no post-tensioning). Everything that works for PGSuper today must keep working: all current validation runs are the regression test for this work. The new elements get roles in the mapping table, so their classes, property sets, and classifications are table data, not C++.

## Where things stand (checked 2026-11-10)

- **Not registered for PGSplice.** `IfcExtensions.Manifest.PGSuper` lists only the four PGSuper categories, so PGSplice has no IFC menu items, no command line, no georeferencing, and no mapping table setting. `CIfcExtensionAgent::RegisterUIExtensions` goes through `IExtendPGSuperUI` only.
- **Export has a PGSplice branch that has never run** (`CreateBridge`, `IsPGSpliceDocument`): one `IfcElementAssembly.GIRDER` per girder, an `IfcBeam` per segment with a "Spliced_Girder_Type" `GIRDER_SEGMENT` type, prismatic/parabolic/double-parabolic segment solids, and each closure joint as an `IfcBeam` (ObjectType `CLOSUREJOINT`, own solid from `CreateClosureJointRepresentation`, `closure_joint` role property sets). Segment GlobalIds are already kept for round trips (`segment_global_ids`). Missing: post-tensioning, temporary supports, closure joint material, reinforcement, and type; all-segment output in the single-girder export (`CreateGirder` writes segment 0, `IfcExporter.cpp:3131`).
- **Import is PGSuper only:** span-indexed layout, `GetSegment(0)` (`IfcBridgeImporter.cpp:375, 962`, `GirderBars.cpp:231`, `Strands.cpp:893`), no closure joints, temporary supports, or tendons. `CPGSuperProjectImporter` only.
- **Tests:** `run_validation.py` and `pgs_extract.py`/`compare_bridge.py` know `.pgs` only. PGSuper has 41 PGSplice regression models (`PGSuper/RegressionTest/Current/PGSplice`).

## What usBridge asks for (workbook 2.0 and `usBrSpec_*.ids`)

**Segments and closure joints.** Both are `IfcBeam.GIRDER_SEGMENT`, told apart by classification: "GirderSegmentConcrete" for a precast segment, "GirderClosurePour" for a closure joint. A one-segment girder stays `IfcBeam.BEAM` "GirderPrecastConcrete" (the workbook notes that covers pretensioned and post-tensioned non-spliced girders, so a PGSplice girder with one segment and a tendon is still a `BEAM`). Today's closure joint (ObjectType `CLOSUREJOINT`, no type) becomes an `IfcBeam` with its own "Closure Joint Type" (`GIRDER_SEGMENT`), classified GirderClosurePour, with the cast-in-place concrete material and `Pset_ConcreteElementGeneral`. The segment and closure joint beams are aggregated by the girder's `IfcElementAssembly.GIRDER` in order along the girder.

**Post-tensioning.** The workbook has a complete set:

| Element | IFC | Key properties (required unless noted) |
|---|---|---|
| Tendon assembly (one per duct) | `IfcElementAssembly.USERDEFINED` `TENDON_ASSEMBLY`, "PrestressingTendonAssembly" | |
| Tendon | `IfcTendon` with an `IfcTendonType` `STRAND` | `usBrPset_ACI_TendonMaterial` (Specification, TendonGrade), `usBrPset_ACI_TendonType.TendonType` = `DUCTED-BONDED`, `usBrPset_ACI_Tendon.StrandOrBarCount`; optional effective forces, TendonMark |
| Duct | `IfcTendonConduit` `DUCT` | `usBrPset_ACI_TendonConduitType`: DuctType, DuctFlexibility (RIGID/FLEXIBLE), DuctShape (ROUND/FLAT/OVAL), DuctSpecification(+Version), Narrow dimensions at and away from ribs; optional friction and wobble coefficients |
| Anchors | `IfcTendonAnchor` `TENSIONING_END` / `FIXED_END` | `usBrPset_ACI_TendonAnchorType`: EncapsulationType, DefaultJackingForce; optional SeatingDistance. `usBrPset_ACI_TendonAnchor` (optional): JackingForce, CalculatedElongation, StressingIndex |

Mapped to PGSuper's data (`CPTData`/`CSegmentPTData`, `CDuctData`, the duct library entry, `GetDuctCenterline`):
- **One `IfcTendon` per duct** (the user's decision, 2026-11-10), StrandOrBarCount = number of strands; the strands in a duct aren't modeled individually. PGSuper doesn't place strands inside a duct, so modeling each strand would be invented geometry. The tendon body is a swept disk along the tendon centerline (the duct centerline offset by the eccentricity of the strands in the duct), with the strand bundle's equivalent diameter. The ducts are swept hollow disks from the duct library entry's diameters.
- **Anchors:** PGSuper's jacking end (start, end, both) gives `TENSIONING_END` at each jacked end and `FIXED_END` otherwise. Pj (calculated or input) goes to JackingForce; anchor set to SeatingDistance; the stressing event's order in the timeline to StressingIndex.
- **Friction** (`μ`, wobble `k`) is PGSuper's loss criteria, per project, not per duct; write it on the conduit type as optional values.
- **Placement in the spatial tree:** a girder tendon's assembly is aggregated by the girder's `IfcElementAssembly`, because it runs through segments and closure joints; a segment tendon's assembly nests in its segment `IfcBeam` the way pretensioned strands do.
- **Missing from PGSuper:** DuctType, DuctSpecification(+Version), the dimensions at ribs, EncapsulationType. Duct library entries have a nominal and inside diameter only. For now (the user, 2026-11-10) these properties are exported empty (the property is written with no value), so the usBridge specifications that require a value fail visibly. Where their values come from is to be addressed later (Future considerations; usBridge shortcoming 4).

**Temporary supports** (the user's decision, 2026-11-10). usBridge has no class for them (usBridge shortcoming 1).
- **Erection towers** are `IfcBridgePart` with PredefinedType `USERDEFINED` and ObjectType `ERECTION_TOWER`, and usBrPset_Common Status "Temporary". They're placed and stationed like piers (an `IfcReferent` with `Pset_Stationing`), and the segment ends they carry rest on them. A new table role `erection_tower` holds the class, so it can change without C++. The `pier` selector (PredefinedType `PIER`) doesn't match them, so today's pier import doesn't mistake them for piers.
- **Strongbacks** aren't modeled. A closure joint between permanent supports that isn't over an erection tower is taken to be on a strongback, both on export (nothing written) and on import (a strongback temporary support at the closure joint).

**Construction sequence** (the user's decision, 2026-11-10). The exporter declares the Alignment-based View (AbV), and `usBrSpec_AbV Exclusions.ids` excludes `IfcWorkSchedule`, `IfcWorkPlan`, `IfcRelAssignsToProduct`/`Process`/`Control`, `IfcRelSequence`, and `IfcSchedulingTime`, so a schedule can't be written (an `IfcTask` alone can't refer to elements). usBridge covers the sequence only in part:
- `usBrPset_ACI_TendonAnchor.StressingIndex` (optional): the order of stressing. Written on each tensioning anchor from the order of the stressing events.
- `usBrPset_ConcretePour.PourID`/`PourSequenceNumber` (optional, described for `IfcSlab`): written on the deck from PGSplice's deck casting regions, and on each closure joint from its casting event (the user, 2026-11-10; the IDS accepts the property set on any element, though the workbook describes it only for `IfcSlab`, see usBridge shortcomings).
- `ConstructionStageNumber` (on `IfcBridge`) is staged construction (e.g. phased widening), and `Status` is a state; neither is a sequence.
- The rest of the sequence (erection on the segments, casting on the closure joints, tower installation and removal, and the event times PGSplice's time-step analysis needs) goes in a PGSplice property set defined in the mapping table: data PGSuper has and the standard lacks, stored as table data like the other targets, so another standard can map it elsewhere.
- The import builds the timeline from whichever of these the model has and fills the gaps, the times in particular, with the default sequence and intervals (S4).
- An `IfcWorkSchedule` export may be added later as an option, for files that don't declare AbV (see Future considerations).

## Phases

Export first (S0-S2), then import (S3-S5). Each phase ends with the PGSuper validation runs unchanged.

### S0. Test models and tools
PGSplice analyzes; it doesn't design, and its time-step analysis is slow. These tests check the IFC implementation only (what's written, and what an import reads back), never analysis results. The overnight PGSuper regression tests catch analysis regressions.
- **Models:** copies of existing PGSuper regression models (`PGSuper/RegressionTest/Current/PGSplice`), chosen for one IFC feature each and cut down (fewest girders and spans that still show the feature) where it helps:
  - `Segment_Tendons.spl`: segment tendons.
  - `Multi_stage_partial_length_pt.spl` or `TwoSpanTubGirders.spl`: girder tendons through closure joints.
  - `CantileverSplicedGirder.spl`: variable depth segments, temporary supports.

  The user approved these three (2026-11-10). A purpose-built model is added only when no regression model has a feature.
- **No analysis in a test run:** under PGSplice the export writes input values only. Targets that need analysis or design (cambers, effective tendon forces, required strengths, the design values of the design-value IDS) are left out for PGSplice documents, so an export never starts the time-step analysis. Check this with the run time of the first export.
- `run_validation.py`: accept `"spl"` models (round trip through PGSplice). Before the first PGSplice model it reads PGSplice's configuration (`PGSplice\Options`, `CatalogServer`/`Publisher`), switches to `"Regression":"Regression"`, and switches back to what it read when done, as it does for PGSuper (the user, 2026-11-10).
- `pgs_extract.py`/`compare_bridge.py`: input data only: groups, segments per girder (length, variation), closure joints (length, concrete, event), temporary supports (station, type), ducts and tendons (type, centerline points, strands, Pj, jacking end, event). Keyed by group, girder, segment, and duct index so PGSuper reports don't change.

### S1. Registration
**Done 2026-11-10** (except the user's check of the menus):
- The manifest registers the exporter and the extension agent in `CATID_PGSpliceDataExporter` and `CATID_PGSpliceExtensionAgent` too. The data importer stays PGSuper-only until S5 (its import assumes one segment per span); `CATID_PGSpliceProjectImporter` comes in S4.
- `CIfcExtensionAgent::GetExtendUI` uses `IExtendPGSpliceUI` or `IExtendPGSuperUI` by the document type.
- **No analysis in an export:** a PGSplice export left the cambers out (`TargetDef::needs_analysis`, `CIfcExportOptions::include_analysis_values`) and still ran the time-step analysis (6 minutes for `Segment_Tendons.spl`). The cause was `IBridge::GetBearingElevationDetails(..., bIgnoreUnrecoverableDeformations = false)`: with haunch depth input its finished elevations need a full analysis. A PGSplice export now asks for the geometry-only elevations; PGSuper exports are unchanged. The export takes 12 s. If an analysis turns out to be unavoidable for some value, the test models get a coarse timeline (few intervals); the analysis's accuracy doesn't matter here (the user, 2026-11-10).
- `/IfcExport` on `Segment_Tendons.spl` writes a `GIRDER` assembly per girder with Segment 1, Closure Joint 1, and Segment 2, closure joints with their solids, and the pretensioned strands. File > Export > IFC in PGSplice is for the user to check.
- **Configuration:** PGSplice has its own configuration (`HKCU\...\PGSplice\Options`, `CatalogServer`/`Publisher`, not PGSuper's `CatalogServer2`/`Publisher2`). The regression models need the Regression configuration (with another library, loading stops at a Library Entry Conflict dialog): `BridgeLink.exe /App=PGSplice /Configuration="Regression":"Regression"`, with a `Regression` catalog server for PGSplice (`FILSYS|%ARPDIR%\PGSuper\RegressionTest\PGSuperRegTestLibrary.lbr|%ARPDIR%\PGSuper\Configurations\PGSplice\WSDOT`).

### S2. Export
1. **Run what's there** on the S0 models, fix what breaks, and check with the usBridge IDS (`check_ids.py`).
2. **Closure joints:** a closure joint type (`GIRDER_SEGMENT`), "GirderClosurePour", cast-in-place concrete material (`Pset_ConcreteElementGeneral`, StrengthClass from the closure joint's f'c), and the closure joint's longitudinal bars and stirrups (the commented-out block in `CreateBridge`) through the same code as segment reinforcement, all in this step (the user, 2026-11-10).
3. **Segments:** classify multi-segment girders' segments "GirderSegmentConcrete"; one-segment girders stay `BEAM`. Review the targets for single-segment assumptions (`girder.span`, cambers) and decide which values go on the girder assembly.
4. **Post-tensioning** (new `PostTensioning.cpp/.h`, the pattern of `Strands.cpp`): tendon assembly, duct, tendon, anchors as above; new `ElementKind`s and table roles `tendon_assembly`, `pt_tendon`, `tendon_conduit`, `tendon_anchor`; targets for strand count, Pj, jacking end, duct diameters, anchor set, stressing event. `Standard.ifcmap` gets the property sets; `Targets.md` documents them.
5. **Erection towers** as `IfcBridgePart` (above), with their station, tower height if PGSuper has it, and the segments they support.
6. **Construction sequence:** StressingIndex, PourID/PourSequenceNumber, and the PGSplice sequence property set (above; after the anchors and towers it's written on), as targets in `IfcTargets.cpp` with their names in `Standard.ifcmap`.
7. **Single-girder export** (`CreateGirder`): all segments, closure joints, tendons.
8. **Design-value IDS:** for PGSplice, input values only (no design or analysis values): specs for segments, closure joints, and tendons.
9. Round trip R1 (`/IfcSourceModel`): segments already have GlobalIds; add closure joints, tendons, ducts, and anchors to the keyed elements.

### S3. Import core: groups, segments, closure joints
Refactor `CIfcBridgeImporter` from span-indexed to group → girder → segment:
- A girder is an `IfcElementAssembly.GIRDER` aggregating segment and closure joint beams (told apart by the table's `girder_segment` and `closure_joint` selectors), or a lone `IfcBeam` (one segment).
- **Groups** from girder continuity: a girder that runs over a pier (through a closure joint or a continuous segment) puts that pier inside its group.
- **Closure joints** from geometry, not names: the gap between consecutive segments along the girder line. Its center sets the closure joint location; the pier or temporary support under it sets where the segment ends are supported. Closure joint concrete from its material and property sets.
- **Segments:** end positions along the girder, and variation (prismatic, linear, parabolic) from the solid, using the girder section work in `ImportGeometry.md`.
- **Temporary supports:** an `erection_tower` part gives an erection tower at its station. A closure joint between permanent piers with no erection tower under it gives a strongback.
- Replace every `GetSegment(0)`; the pretensioning, bar, and stirrup imports run per segment.
- **Done when** the PGSuper runs are unchanged and the S0 models import their group, segment, and closure joint layout (PGSplice importer from S4 needed to see it end to end; until then, check with the data structures in a test build).

### S4. PGSplice project importer and post-tensioning
- `CPGSpliceProjectImporter` in `CATID_PGSpliceProjectImporter` with `IfcImportTemplate.spt` (template decision in `CLAUDE.md`).
- **Post-tensioning import:**
  - A tendon is girder or segment by what it lies in (the girder assembly, or one segment's solid), with the geometry as the check: a tendon that crosses a closure joint is a girder tendon.
  - Duct geometry: fit PGSuper's duct types to the tendon or duct centerline: linear (points), parabolic (low points, inflection points), offset. A parabolic fit within tolerance, else linear with enough points.
  - Strands from StrandOrBarCount, else from the tendon's area and the strand size; Pj from JackingForce or DefaultJackingForce; jacking end from the anchor types.
  - Duct size: match a duct library entry by diameter, or create one (as girder library matching does, Phase 5 of `IfcImportPlan.md`).
- **Timeline:** from the sequence properties (StressingIndex, PourSequenceNumber, the PGSplice property set) where the model has them; the gaps from a default: erect segments (on temporary supports), cast closure joints, stress girder tendons in StressingIndex order and remove the temporary supports in the same event, cast deck, open to traffic, with default intervals between events (the user's choice, 2026-11-10). An `IfcWorkSchedule`, if a later option writes one, is read first.
- A model with spliced girders imported by PGSuper says it needs PGSplice instead of importing wrongly; a PGSuper-shaped model imports into either.

### S5. Data importer and round trip under PGSplice
Import into an open PGSplice project and the R1 round trip ([RoundTripDesign.md](RoundTripDesign.md)).

## Likely PGSuper API needs
- A valid default timeline for a bridge description, so the importer needn't rebuild PGSplice's timeline rules.
- Creating segment layouts and temporary supports from positions, with the UI's checks.
- The tendon centerline (duct centerline plus strand eccentricity in the duct) through a documented call.
- Which values need analysis results, so a PGSplice export can leave them out without starting the analysis.

## Open questions
1. ~~S0 models~~ Decided 2026-11-10: the three regression models in S0, cut down where it helps.
2. ~~Tendons~~ Decided 2026-11-10: one `IfcTendon` per duct with StrandOrBarCount; the strands in a duct aren't modeled individually.
3. ~~Temporary supports~~ Decided 2026-11-10: erection towers as `IfcBridgePart` (`USERDEFINED`, `ERECTION_TOWER`); strongbacks not modeled, assumed at closure joints between permanent supports that aren't on an erection tower.
4. ~~Duct properties PGSuper doesn't have~~ Decided 2026-11-10: exported empty for now; their source is a later question (Future considerations).
5. ~~Default timeline~~ Decided 2026-11-10: the S4 sequence, with the temporary supports removed in the tendon stressing event. Closure joints get `usBrPset_ConcretePour` (PourID, PourSequenceNumber).
6. ~~Closure joint reinforcement~~ Decided 2026-11-10: with the rest of the closure joint work (S2 step 2).

## usBridge shortcomings (to raise with the usBridge authors)
Found while planning PGSplice support (2026-11-10). The plan works around each one; the workaround is noted.
1. **Temporary supports:** there's no class or classification for an erection tower or a strongback, only the "Temporary" Status. Workaround: erection towers as `IfcBridgePart` `USERDEFINED` `ERECTION_TOWER` (an `erection_tower` table role, so it can change without C++); strongbacks not modeled, assumed at closure joints between permanent supports that aren't on an erection tower. If usBridge defines them, strongbacks may need explicit elements (position, the segments they connect, removal event).
2. **Construction sequence:** the erection sequence is left to referenced documents (the workbook's `IfcDocumentReference` example "Proposed Erection Sequence"). In the model there are only StressingIndex (on the optional `usBrPset_ACI_TendonAnchor`) and PourID/PourSequenceNumber; nothing for segment erection, temporary support installation and removal, or the times between events. The Alignment-based View excludes `IfcWorkSchedule` and the assignment relationships, so the sequence can't be carried as a schedule either. Workaround: a PGSplice sequence property set in the mapping table.
3. **Concrete pours on closure joints:** `usBrPset_ConcretePour` is described only for `IfcSlab`, but a girder closure joint (`IfcBeam.GIRDER_SEGMENT`, GirderClosurePour) is a cast-in-place pour with a sequence too. Workaround: the property set is written on closure joints (the IDS applies to any element with it, so it validates).
4. **Required duct and anchor properties a design program doesn't have:** DuctType, DuctSpecification and DuctSpecificationVersion, NarrowDimensionAtRibs, and EncapsulationType are required, but they're usually the post-tensioning supplier's choice, not known at design. Workaround: exported with no value, so the specifications fail visibly.

## Future considerations
- **`IfcWorkSchedule` option (the user, 2026-11-10):** a later export option may write the timeline as an `IfcWorkSchedule` with `IfcTask`s assigned to the segments, closure joints, tendons, and towers, with their times. It can't be used in a file that declares AbV, so the option needs a different view definition (or a separate file) and its own validation.
- **Duct properties PGSuper doesn't have (2026-11-10):** where their values come from is to be decided: mapping table defaults, project data stored by the plugin (as georeferencing is), or new PGSuper duct library data (usBridge shortcoming 4).
