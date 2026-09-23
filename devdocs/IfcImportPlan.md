# IFC → PGSuper Model Import Improvement Plan

Drafted 2026-09-23. Last completed work before this plan: estimating deck slab edges from geometry.

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

1. `GetProperty` (`Properties.h`) ignores `IfcPropertySingleValue.Unit`, so the PGSuper export's `ReleaseStrength = 4.0 ksi` is imported as f'ci = 4 Pa (same for f'c). `InitUnits` handles only length and angle.
2. `get_beam_spacing`: `nSpans = get_pier_count(file) + 1` should be `- 1` (pier count includes abutments).
3. Assumed girder keys are computed differently in `Import()`, `SetGirderProperties()` and `get_beam_spacing()`. Beam ordering comes from iteration order, not position.
4. The girder start/end is decided by distance from the origin (`BeamSpacing.cpp`). It should use alignment station.
5. Spacing is the 3D distance between the centres of girder ends, but it's stored as `AlongItem`/`AtPierLine`. It should be normal to the alignment (or along the pier line) at the CL bearing.
6. `piers.insert(abutments.front())` crashes if there are no ABUTMENT parts.
7. Girder orientation is hard-coded to the first supported option.
8. The mesh → decompose → top surface → boundary → segments pipeline is copied three times.

## Phases

### Phase 0: Test setup
- Round-trip test: PGSuper export → import → compare `CBridgeDescription2` field by field.
- Python ifcopenshell "expected values" scripts for the PennDOT and Iowa models.
- Structured import report: parameter | value | source (pset/geometry/assumed/default) | confidence.

### Phase 1: Correctness of what we already import
- Fix bugs 1–6, plus the baseline findings: the beam-type count/`beams.begin()` bug, compute each girder key only once, and fix the log messages.
- Single geometry pass: iterate once over beams, bearings and the slab, and cache meshes by product id. Consider `num_threads > 1`, but the progress stream isn't thread-safe. Honour per-property units, then project units. Add PRESSUREUNIT, MASSDENSITY and FORCE.
- One place that owns the girder layout: get each beam's CL end points from geometry → station/offset → span by pier stations → order by offset → cross-check against parsed names.
- Pull the geometry pipeline out into one reusable helper.

### Phase 2: Superstructure geometry
- Pier/abutment orientation (skew) from bearing centroids or girder end segments → `CPierData2::SetOrientation`.
- Spacing measured correctly at the CL bearing, normal to the alignment. Detect uniform spacing (`sbsUniform`).
- Girder end distance and bearing offset from beam end, bearing centroid and pier line.
- `CBearingData2` from the bearing bounding box and fixity psets.
- Integral abutments → `SetBoundaryConditionType`.

### Phase 3: Deck
- Gross depth and overhang edge depth from cross-sections of the slab mesh. Fall back to psets.
- Deck type (CIP/SIP) from form-type psets.
- Haunch/slab offset from deck soffit minus girder top (Iowa: use the haunch slabs).
- Deck concrete f'c. Replace the stopgap in `condense()`.

### Phase 4: Alignment/profile when the model has none
- Derive the vertical profile and crown slope from the deck top surface.
- Use deck station psets (Iowa) for the reference station.

### Phase 5: Girder section and material mapping
- Section match order: `IfcBeamType.Name` → vendor pset/classification value (Iowa `IaDOT_PPCB.2_Type` / classification "Beam, PPC, BTB45"; PennDOT `_PS Concrete Beams.Type`) → user alias table → cross-section shape matching → prompt the user. Replace the "first I-Beam" fallback. Because Iowa has no `IfcBeamType` at all, this is the most visible gap in the baseline and a candidate for moving earlier.
- Orientation from the measured beam roll/plumb.
- Material lookup chain: standard psets → configurable vendor pset aliases (with string/unit-in-name parsing) → default.

### Phase 6: Prestressing
- Strand `IfcTendon` → section cuts → `CStrandData` direct fill or strand rows. Harped strands and harp points from the tendon path, debonding from tendon length. Strand material, Pjack from `InitialTension`.

### Phase 7: Secondary elements and loads
- Barriers → traffic barrier library match or weight per length. Diaphragms → diaphragm loads. Pier caps and columns → physical pier model (low priority).

## Open decisions
1. Where the vendor-pset and girder-name alias tables live (import options UI, JSON file, or PGSuper library). Leaning JSON.
2. When geometry and psets disagree: suggest geometry wins for layout, psets win for materials, and every conflict is flagged in the report.
