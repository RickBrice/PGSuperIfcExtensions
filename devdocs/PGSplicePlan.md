# PGSplice Import and Export Plan

Drafted 2026-10-15. Waiting for review.

## Goal

The extension was written for PGSuper. PGSplice is the general case, and PGSuper is a subset of it:

| | PGSuper | PGSplice |
|---|---|---|
| Segments per girder | 1 | 1 or more, with closure joints between them |
| Girder group | exactly one span | one or more spans |
| Supports | piers and abutments | also temporary supports (segment ends between piers) |
| Prestressing | pretensioning | pretensioning and post-tensioning (girder tendons, segment tendons, ducts) |
| Construction | assumed sequence | a timeline (time-step analysis) |

Both use the same bridge model (`CBridgeDescription2`). The aim is one code path written for the PGSplice model, where PGSuper is the special case: one segment per girder, groups one span long, no post-tensioning. It should not be two parallel code paths.

Following the direction in `CLAUDE.md`, the new element kinds (segments, closure joints, temporary supports, ducts, post-tensioning tendons) get roles in the mapping table, so their IFC classes, property sets and classifications are table data, not C++.

## Where things stand (2026-10-15)

**Registration (the first blocker).** `IfcExtensions.Manifest.PGSuper` registers the exporter, data importer, project importer, and extension agent only in the PGSuper categories. PGSplice loads none of them:
- There's no IFC export or import in PGSplice's menus.
- The extension agent isn't loaded, so there's no command line (`/IfcExport` on an `.spl` stops at an "invalid command line" message box), no georeferencing (`IGeoreferencing`), and the project doesn't remember its mapping table.
- `CIfcExtensionAgent::RegisterUIExtensions` goes through `IExtendPGSuperUI` only, so the alignment editing callback would need `IExtendPGSpliceUI` under PGSplice.

**Export: mostly general already.**
- `CreateBridge` has a PGSplice branch (`IsPGSpliceDocument`). It writes one `IfcElementAssembly.GIRDER` per girder with an `IfcBeam` for each segment, a shared "Precast Girder Type" with `GIRDER_SEGMENT`, and closure joints as `IfcBeam` with ObjectType `CLOSUREJOINT`, each with geometry (`CreateClosureJointRepresentation`).
- Segment geometry handles prismatic, parabolic and double-parabolic segments (`CreatePrecastSegmentRepresentation`).
- `IdsExporter` loops over segments and uses `GIRDER_SEGMENT` for spliced girders.
- **Missing:**
  - post-tensioning: ducts, girder and segment tendons, anchorages. Nothing is exported.
  - temporary supports
  - the timeline: `ConstructionSequence.h` only does PGSuper's assumed sequence, and its calls are commented out.
  - segment-level values for closure joints in the table (only a `closure_joint` role exists)
- The `CreateGirder` path (export of one girder, `options.girderKey`) writes only segment 0.
- Never run on a PGSplice project, because it can't be reached (see Registration).

**Import: PGSuper only.**
- `CIfcBridgeImporter::Import` assumes each span is a girder group: `layout` is indexed by span, and girder counts and spacing are set per span group.
- Each girder is one segment (`GetSegment(0)`).
- There's no notion of a girder as an assembly of segments, of closure joints, of temporary supports, or of post-tensioning.
- `CPGSuperProjectImporter` is registered only for PGSuper and starts from `IfcImportTemplate.pgt`.

**Test tools: PGSuper only.**
- `run_validation.py` takes `"pgs"` models.
- `pgs_extract.py` reads one `PrecastSegment` per girder and groups as written (fine), but no segments, closure joints, temporary supports, tendons, or timeline.
- PGSuper has 41 PGSplice regression models (`PGSuper/RegressionTest/Current/PGSplice/*.spl`) to choose round-trip models from.

**Run time.** PGSplice uses time-step analysis, which is slow. Some export targets need analysis results (camber, screed camber). Test models must be small, or analysis-dependent values must be optional for these runs.

## Phases

### S0. Test setup
- Choose 2-3 small PGSplice models covering:
  - several segments per girder, with closure joints over piers
  - a group over several spans, with temporary supports
  - girder and segment tendons
  - a variable-depth segment

  Candidates: `Segment_Tendons.spl`, `TwoSpanTubGirders.spl`, `Multi_stage_partial_length_pt.spl`, `CantileverSplicedGirder.spl`. Check their run times first.
- `run_validation.py`: accept `"spl"` models (round trip with the PGSplice template), and keep the PGSuper runs unchanged.
- `pgs_extract.py`/`compare_bridge.py`: extract and compare the following. Matched by girder and segment index, so PGSuper reports don't change.
  - groups (start and end pier)
  - segments per girder (length, end positions, variation)
  - closure joints
  - temporary supports (station, type)
  - ducts and tendons (count, path type, strands, jacking)
  - the timeline events, summarized
- A way to keep runs short: an export option to leave out analysis-dependent values, or models whose analysis is fast.

### S1. Registration and host integration
- Add the PGSplice categories to `IfcExtensions.Manifest.PGSuper`:
  - `CATID_PGSpliceDataExporter`
  - `CATID_PGSpliceDataImporter`
  - `CATID_PGSpliceExtensionAgent`
  - later, `CATID_PGSpliceProjectImporter` (S4)
- Extension agent: use `IExtendPGSpliceUI` when the document is PGSplice.
- Check that every interface the agent and exporter use exists in PGSplice.
- Acceptance: `/IfcExport` works on the S0 models, and File > Export > IFC appears in PGSplice. PGSuper results are unchanged.

### S2. Export completeness
- Run the export on the S0 models and fix what breaks. Compare the export with `check_ids.py` against the usBridge IDS. The workbook has `IfcBeam.GIRDER_SEGMENT` "GirderSegmentConcrete" and "GirderClosurePour".
- Table roles and classifications:
  - `girder_segment`: segments of a multi-segment girder. PGSuper girders stay `girder` (`IfcBeam.BEAM`, GirderPrecastConcrete).
  - `closure_joint`: classification and property sets.
  - `temporary_support`
- **Temporary supports:** element class and representation, to decide with the usBridge workbook.
- **Post-tensioning:**
  - ducts as `IfcTendonConduit` with their path geometry
  - girder and segment tendons as `IfcTendon` with `POSTTENSIONED` predefined type and their strands
  - anchorages (`IfcTendonAnchor`) if wanted
  - targets for duct size, number of strands, jacking force or stress, friction and anchor set, and the jacking end
- **Targets:**
  - Review the 62 targets for single-segment assumptions, e.g. `girder.span` and the camber targets for a segment that isn't the whole span.
  - Decide what group-level values go on the `IfcElementAssembly` girder.
- `CreateGirder` (single girder export): all segments and closure joints.
- Timeline (optional, later): an `IfcWorkSchedule` from PGSplice's events. This replaces the commented-out `ConstructionSequence.h` calls for both applications.

### S3. Import core: generalize to groups and segments
Refactor `CIfcBridgeImporter` from a span-indexed layout to group → girder → segment:
- **Girders:**
  - An `IfcElementAssembly.GIRDER` aggregating `IfcBeam` segments is one girder.
  - A lone `IfcBeam` is a one-segment girder.
  - The table selectors decide which elements are segments, girders and closure joints.
- **Groups:** a group runs from the pier where its girders start to the pier where they end. A girder that runs continuously over a pier means that pier is inside the group.
- **Segments:**
  - segment ends along the girder
  - closure joints between segments (at a pier or a temporary support)
  - segment variation from the geometry or properties
- **Temporary supports:** wherever segment ends fall between permanent piers.
- **Pier connections:** continuous segment, closure joint at the pier, or the group ends there.
- **Acceptance:** all current validation runs are unchanged (PGSuper models are the one-segment, one-span-per-group case), and the S0 PGSplice models import their group and segment layout.

### S4. PGSplice Project Importer
- A new component `CPGSpliceProjectImporter` in `CATID_PGSpliceProjectImporter` (pattern: `IEPluginExample/PGSpliceProjectImporter`), with a template `IfcImportTemplate.spt` (see the template decision in `CLAUDE.md`).
- The shared import core (S3) fills in the bridge. The project importer chooses the template, the document type, and the steps only PGSplice has.
- **Timeline:** PGSplice needs one before it can analyze. The importer builds it from the model's `IfcWorkSchedule` if there is one, or else a default sequence: erect segments, cast closure joints, stress tendons, cast deck, open to traffic.
- **Post-tensioning:** ducts, tendons and jacking from the IFC.
- **Which application:**
  - An IFC model with spliced girders imported with the PGSuper project importer is reported as needing PGSplice, rather than being imported wrongly.
  - A PGSuper-shaped model can be imported by either application.

### S5. Data importer and round trip under PGSplice
The data importer (import into an open project) and the R1 round trip ([RoundTripDesign.md](RoundTripDesign.md)) for PGSplice projects, after S3 and S4.

## Likely PGSuper API needs
These match the standards direction in `CLAUDE.md`: what can be done through a documented API doesn't need code in the extension.
- **Timeline:** creating a default timeline for a bridge description, or documented rules for a valid one, so the importer doesn't have to rebuild PGSplice's timeline logic.
- **Temporary supports and segments:** creating segment layouts and temporary supports from positions, with the checks the UI does.
- **Girder library entries for spliced families** (the spliced I-beam, U-beam and NU-beam factories). This is the same question as Phase 5 in [IfcImportPlan.md](IfcImportPlan.md).
- **Analysis-free values:** a documented way to tell which values need analysis results, so an export can leave them out quickly.

## Open questions
1. Which PGSplice models to use for S0 (and their run times)?
2. Post-tensioning in IFC: ducts as `IfcTendonConduit` and tendons as `IfcTendon` `POSTTENSIONED`? Model each strand, or one tendon per duct?
3. Temporary supports: which IFC class and usBridge classification? Exported as elements, or only as positions?
4. Default timeline for imported PGSplice projects: the sequence in S4, or from an agency template?
5. Should a PGSuper-shaped IFC model (one segment per girder) be importable into PGSplice too? This plan assumes yes.
6. Order: S1-S2 (export) before S3-S4 (import), as listed, or the PGSplice project importer sooner?
