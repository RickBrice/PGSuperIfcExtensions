# Girder Import Plan: Section, Strands, Transverse Reinforcement

Drafted 2026-10-15 and reviewed the same day (see "Decisions"). **Status (2026-10-16):** G1 steps 1-3 are done (the section from the model, matching by geometry, fitting dimensions and cloning the closest entry); step 4, a new entry with a generated strand grid, is next. G2 started 2026-10-20 (see G2). G3 is not started. Priorities 2 and 3 in `CLAUDE.md` (girder mapping with the PGSuper library, information from geometry). It builds on Phase 5 (girder section) and Phase 6 (prestressing) of [IfcImportPlan.md](IfcImportPlan.md).

Before G1, the importer set the girder type by name only (`GetGirderLibraryEntry`, falling back to the first I-beam). The strands, stirrups, confinement, and splitting bars still come from the import template unchanged (G2 and G3), so an imported girder has the template's reinforcement.

## What the test models have (2026-10-15)

| | PGSuper export | Iowa US59 | PennDOT |
|---|---|---|---|
| Girder geometry | `IfcSectionedSolidHorizontal` of the exact PGSuper sections | `IfcMappedItem` of swept sections | `IfcFacetedBrep` |
| Girder type name | `IfcBeamType` = library entry (WF66G, WF100G) | `IaDOT_PPCB.2_Type` (BTB45) | `_PS Concrete Beams.Type` |
| Strands | 369 `IfcTendon`, `IfcSweptDiskSolid` along a polyline (start, harp points, end), types with debonding (`usBrPset_TendonDebonding`) | 661 `IfcTendon` STRAND in "Strands_Beam n_Span m" assemblies, mapped `IfcExtrudedAreaSolid` (straight, apparently scaled by the mapping), `IaDOT_Strand` properties as text (`0.5"`, `270`, `Low Relaxation`); 81 more are lifting loops (NOTDEFINED). No harping or debonding words in the file | none |
| Girder bars | 2354 `IfcReinforcingBar` with types: G2 stirrups, G3 top bars, G9/G10 confinement, longitudinal; `usBrPset_ACI_ReinforcingBar.BarUse` (STIRRUP, TIE, ...) | 6080 bars in "Rebar_Beam n_Span m" assemblies, mapped extrusions, ProStructures marks and ACI bend shapes (`4d1` IA56, `4e1` S10, `4h1` S11, `6b3` 20, `4c1` straight), lengths as text (`2:3`) | none |

## PGSuper data the import fills

- **Section:** the girder library entry of each girder (`GirderLibraryEntry`, a beam factory with named dimensions). Entries are created with the existing library API (see "Library entry creation").
- **Strands:** `CStrandData`. With `pgsTypes::sdtDirectStrandInput` (available in PGSuper and PGSplice) each `CStrandRow` is one strand: its type (straight, harped, temporary), its horizontal position `m_Z`, its elevation at the start, the two harp points, and the end (`m_Y`, from the top or bottom face), its debonding at each end, and extension. The harp points are on `CStrandData` (`SetHarpPoints`). Strand material per type, and jacking.
  The other definition types fill the library entry's strand grid (`sdtTotal`, `sdtStraightHarped`, `sdtDirectSelection`).
- **Transverse reinforcement:** `CShearData2`: shear zones (length, bar size, spacing, number of legs, horizontal interface bars), the splitting zone (bar size, spacing, length, number of bars), the confinement zone (bar size, spacing, length), symmetric or not.

## G1. Girder section

Phase 5 of [IfcImportPlan.md](IfcImportPlan.md) is the plan; in short:
1. Get the section polygon of each girder from the model: the swept profiles (PGSuper export, Iowa), or a cut through the solid (PennDOT). Several cuts to see end blocks and tapers; the prismatic middle portion is the section.
2. Match a library entry by geometry (depth, flange widths and thicknesses, area, I, centroid), confirming a name match. A geometric match with another name is used and logged.
3. If nothing matches: classify the family and fit its dimensions (minimize the difference between the polygon and the factory's section).
4. Create the library entry: clone the closest entry of the family and set the dimensions, or start from the factory's defaults. The existing PGSuper library API is enough (see "Library entry creation").

Validation: under the Regression configuration the WF girders aren't in the library, so the round trip must create entries that match the source project's.

## G2. Strands

1. **Find the strands of each girder:** the table's tendon role (selector), through the girder's assemblies (PGSuper export: tendon bundle aggregated in the beam; Iowa: "Strands_Beam n_Span m" assemblies, related to the beam by name or position). Lifting loops and other tendons are left out by the selector (e.g. Iowa's NOTDEFINED).
2. **The centerline path of each strand in girder coordinates:** from the swept disk directrix (polyline or indexed curve), or the axis of an extruded solid, through the mapping and placements. Generic, like the bar paths of the export (`BarShape.h`), so the same code reads any representation that has an axis.
3. **Classify:**
   - harped: the path has interior break points; the harp points are where they are (all harped strands share them in PGSuper: use the common points, report outliers)
   - temporary: the table's `tendon.status` location (`Pset_ElementComponentCommon.Status` = TEMPORARY), else position in the top flange (a hint, not data, per G3 of the design note)
   - debonded: the table's debond locations (`usBrPset_TendonDebonding`), else none. Strands shorter than the girder, or pieces of a strand (Iowa's short extrusions?), are a question to settle with the Iowa model
4. **Material and jacking:** diameter and grade from the type or properties (`IfcTendonType.NominalDiameter`, `Pset_MaterialSteel`, Iowa's text `0.5"` and `270` through the table) to the PGSuper strand pool; jacking from `IfcTendon.TensionForce`/`PreStress` or `InitialTension`.
5. **Fill `CStrandData`:** with the library grid when every strand is on it (`sdtDirectSelection`), else with `sdtDirectStrandInput` (decision 2).

Validation: `pgs_extract.py` extracts the strands (count per type, positions at the ends and harp points, debonding, harp points, material, jacking) and the round trip compares them. Iowa: an inventory and the counts from the model.

**Status (2026-10-20):**
- **Validation:** the strands are compared physically rather than from the `.pgs` data (`Tests/ImportValidation/ifc_strands.py`): the round trip exports the imported project again and compares its strands with the first export's (per girder: counts by type, material, prestress, tension force, centerlines within 1/8 in, debonding, harp points). PGSuper's strand definitions (grid fills, rows, individual strands, harped offsets, extended strands) don't need re-implementing in Python.
- **Done (`Strands.cpp`), steps 1-4 and the individual-strand half of step 5:** strands by the tendon role, assigned to girders through their assemblies (else by geometry); centerlines from swept disk polylines in girder coordinates; harped by a bend, temporary by `tendon.status`; material from the diameter and `tendon.fpu`/`tendon.grade`/`tendon.corrosion_treatment`; jacking from `PreStress` (else `TensionForce`); debonding from `tendon.debond_start`/`end`; filled as individual strands (`sdtDirectStrandInput`) with common harp points.
- **Grid fill (step 5, 2026-10-20):** when every strand is at a grid position of the girder's library entry, the grid is filled (`sdtDirectSelection`), else individual strands (the reason is logged). A pair position holds two strands at ±X (both at the center where its X is 0) or one strand at the center. The coincident case is the WSDOT library (the user, 2026-10-20): WSDOT bundles harped strands at the harp points, 12 per bundle, and its library models every harped position at X = 0 there. This is unique to WSDOT; other libraries, e.g. TxDOT, keep a 2 in grid at the harp points, so don't treat bundling as general (and generated grids, step 4, use the 2 in grid). The fill is found at the girder start first, where the strands of a pair are apart in both kinds of library; harped strands may have a common vertical offset at each end and harp point (found by trying each grid position for the lowest strand), which the entry must allow. The offsets are input as PGSuper measures them from the strands: top strand to top of girder at the ends (`hsoTOP2TOP`), bottom strand to bottom of girder at the harp points (`hsoBOTTOM2BOTTOM`). Debonded strands must be at debondable straight positions, both strands of a pair alike.
- **PGSuper round trip:** every girder is filled into its entry's grid (including the fitted clone WF100G) and every strand value matches (311 of 473 values, 0 mismatches).
- **Expected results by configuration (the user, 2026-10-20):** the PGSuper round-trip models (`PGSuper_Import_Model` and similar) use WF66G and WF100G, which are defined in the WSDOT configuration. Two tests: (1) with the WSDOT:WSDOT configuration (or WSDOT (Local)), the girders match the library entries exactly and the strands fill their grids; (2) without it (e.g. Regression), the import creates project-specific library entries with default strand grids at 2 in spacing (G1 step 4), and the imported strands are input as individual strands (user-defined strands), not in a grid. Today, without the WSDOT configuration, the girders match the template's own entry "Unknown_Girder_Type" (a WF66G under another name) and WF100G is a clone of it with WSDOT's grid, so test (2) fills WSDOT grids; matching the template's entry also goes against keeping the import independent of template values. Fix with G1 step 4.
- **Iowa strands (2026-10-20):** the model draws each strand as a chain of circle extrusions (straight pieces) and revolved solids (bends) in a mapped item, with bent projections beyond the girder ends, and draws its girders and strands along the profile. The importer chains the straight pieces through the bends, clips the centerline to the girder's plan footprint (densified, so strands whose vertices are all outside the girder are found), logs the projecting strands (extended strands aren't imported yet), skips tendons that don't run the girder's length, and measures depths from the girder's top sampled along its length (so strands parallel to a curved girder are straight). The model's predefined types are unreliable (span 1's prestressing strands are mostly NOTDEFINED and the lifting loops are STRAND), so the Iowa table selects strands by name ("Prestressing Strand") and maps `tendon.diameter` (the model draws 0.6" strands as 5/8" circles), `tendon.grade`, and `tendon.tension_force` ("07-Initial Prestress", text in kips). Result with the Iowa table: every girder has the model's strand count (span 2: 38 straight and 12 harped with harp points at 42.65 and 63.65 ft; spans 1 and 3: 12 straight), input as individual strands (no library entry has their grid). Question for the model: in Beam 1 Span 1 six bottom strands vary 1.5-3 in in depth along the girder while its other strands don't, so they're imported as harped. With the standard table, span 1 gets one strand per girder (the NOTDEFINED strands aren't selected) and no jacking force.
- **Extended strands (the user, 2026-10-20):** PGSuper stores only whether a strand is extended at each end (`CStrandRow::m_bIsExtendedStrand`, the grid's extended strand lists), not the shape of the extension. The import can set the flag from a projection, but Iowa's bent detail can't be kept, and an export won't reproduce it.
  - **Done (2026-10-20):** a straight strand that projects beyond a girder end (more than 1 in outside the girder's plan footprint) is extended at that end, unless it's debonded there; for a grid fill both strands of a pair must agree. PGSuper extends only straight strands (`AddStraightStrand` takes the flags), so projecting harped and temporary strands are logged and only their part in the girder is imported. Iowa: 6 extended straight strands at each end of every girder (span 2's 6 projecting harped strands per girder can't be extended; Beam 1 Span 1's projecting strands are its six depth-varying ones, imported as harped).
- **Harp points in grid mode (2026-10-20):** a grid fill's harp points come from the library entry (location, reference, measure, minimum distance), so the import computes them as WBFL's `CStrandModelBase::GetHarpPointLocation` does (end distances from the girder's bearings in the model) and fills the grid only if they're within 1/8 in of the model's; otherwise individual strands, which carry the model's harp points. With WSDOT (Local) the WF entries' harp points are the model's; the generated entries' default harp points (26.2 and 71.6 ft for the 97.8 ft WF66G girders) aren't.
- **Next:** other strand representations (indexed poly curves); a model with debonding (none of the round-trip models has any); round trips by girder family (box beams, voided slabs, U-beams).

## G3. Stirrups, confinement, and splitting bars

1. **Find the bars of each girder** (the cage assemblies, as for strands).
2. **The role of each bar:** from the table where the model says it (e.g. `usBrPset_ACI_ReinforcingBar.BarUse` STIRRUP/TIE in the PGSuper export, bar marks or bend shapes for an agency such as Iowa's IA56), else from geometry: vertical legs through the web are stirrups (legs = vertical legs crossing mid-depth), bars enclosing the bottom flange near the ends are confinement, vertical bars within h/4 of the ends are splitting reinforcement. Each geometric guess is logged.
3. **Zones:** bars of a role are ordered along the girder and grouped where size, legs, and spacing stay the same; zones are measured from the girder ends. Symmetric zones when both ends match.
4. **Fill `CShearData2`:** shear zones, the splitting zone (or `bUsePrimaryForSplitting` when the stirrups are the splitting bars), the confinement zone. Horizontal interface bars: the stirrup legs extending into the deck.

The bar paths are read the same way as strands (step 2 of G2), from the bar types (PGSuper export) or the mapped extrusions (Iowa).

Validation: `pgs_extract.py` extracts the shear data and the round trip compares it. The PGSuper export's stirrups are the WSDOT G2/G9/G10 details made up by the exporter (see `CreateStirrups`), so the round trip checks the zones, not the bar shapes.

## Order

G1 first: the strand elevations are measured from the girder's top or bottom and the zones need the girder depth. G2 next (it carries the path reading that G3 reuses), then G3. The round-trip models exercise all three; Iowa exercises the agency paths (mapped extrusions, text properties, marks).

## Decisions (2026-10-15)

1. **Library entries:** a PGSuper interface to create girder library entries is OK. Its details are worked out before coding (see "Library entry creation" below).
2. **Strands: use the library's strand grid when possible.** Two use cases:
   - **(a) The girder has strands: the import is for analysis.** The strands are imported as they are (on the grid when every strand is on it, else as individual strands); the import doesn't change them.
   - **(b) The girder has no strands: the import is for design.** The priority is a library entry that already exists (by geometry, confirming the name), then cloning the closest entry and adjusting its dimensions, then a new entry with a strand grid generated from the section perimeter. The strands are left for design in PGSuper.
3. **Bar roles in mapping tables:** tabled. PGSuper's rebar model is simple, so a table section for bar roles may not be worth it.
4. **Longitudinal girder bars:** imported too (`CLongitudinalRebarData`).

## Library entry creation

**What PGSuper already has (2026-10-15).** Everything needed is public today:
- `ILibrary::GetGirderLibrary()` gives the project's girder library for changing; the library framework has `CloneEntry`, `AddEntry`, `NewEntry`.
- `GirderLibraryEntry`: `SetBeamFactory`, `SetDimension(name, value, bAdjustStrands)`, `ClearAllStrands`, `AddStraightStrandCoordinates`/`AddHarpedStrandCoordinates`/`AddTemporaryStrandCoordinates`, harping point and adjustment limits, `ValidateData`.
- `psgLib::BeamFamilyManager` lists the girder families and their factories. A `BeamFactory` builds the section for a set of dimensions (`CreateGirderSection`), and lists, validates, and defaults its dimensions (`GetDimensionNames`, `ValidateDimensions`, `GetDefaultDimensions`).

**The steps for a girder** (case (b), and case (a) when the strands need a girder):
1. **Measure** the section from the model (G1 step 1): outline polygon and voids, depth, flange widths and thicknesses, area, I, centroid.
2. **Existing entry:** compare with every girder entry of the project library (its factory's section for its dimensions). An entry within tolerance is used. A name match (IFC type name, or the table's girder type locations) is confirmed by geometry; a name match that doesn't fit is reported, and the geometric match is used.
3. **Clone and adjust:** the closest entry of the matching family is cloned and its dimensions fitted (`SetDimension` with strand adjustment); the result must validate (`ValidateDimensions`, `ValidateData`) and fit within tolerance.
4. **New entry:** a family is chosen from the section topology, its dimensions fitted starting from the factory defaults, and a strand grid generated from the section perimeter (rules below).
5. **Strands:** see decision 6 below.

**Decided (2026-10-15):**
1. **Where the code lives:** in the extension, with the existing PGSuper API (option i).
2. **Matching tolerances:** 1/8 in for depth, flange widths, and thicknesses; 1% for area and I. The same for every family for now, but kept as data so they can differ by family later.
3. **Fitting:** all of the family's dimensions. Girders of the same family can differ in more than depth (an Iowa and a WSDOT girder can both be WF girders with very different flanges).
4. **Names:** the IFC type name with " (IFC)" (e.g. "BTB45 (IFC)"; the suffix added by the user, 2026-10-20, so a project opened later under a configuration with an entry of the type's name doesn't stop at a library conflict), with " (2)", and so on, when the name is already used by a different girder.
5. **Generated strand grid:** 2 in spacing with 2 in cover to the strand center. Straight strands in the bottom flange, harped strands in the web, temporary strands in the top flange. No strands are debondable (harping is assumed, not debonding; harping and debonding could become a configuration option later).
6. **The user's intent isn't known** when a project is created from an IFC file, so:
   - **The girder has strands:** if they're all on the grid of the girder's library entry, a grid fill; otherwise individual strands.
   - **The girder has no strands:** the library entry gets a strand grid if it doesn't have one.
   - **A cloned entry** keeps the clone's grid. **A new entry** gets a generated grid (2 in spacing) even when the girder's strands are imported as individual strands: designers will move to the grid input, and it needs a grid.
7. **Where entries go:** the project's own library (saved with the project), never the master library.

**Revised (2026-10-20, the user):** step 3 (clone and adjust) is dropped. When no library entry matches the section exactly, the import creates a new project entry from the beam factory with its dimensions fitted to the section and a generated 2 in strand grid (step 4). Cloning only made sense to keep the clone's grid, and a generated grid makes the result independent of whichever entry happened to be closest. The strands of a girder on a new entry fill its generated grid when they fit it, else they're input as individual strands (user-defined strands), like any other girder's (revised again by the user the same day). WSDOT's harped strands (±1 in at the ends, bundled at X = 0 at the harp points) don't fit a generated 2 in grid, so without the WSDOT configuration they're individual strands, which have one offset per strand: the PGSuper round trip under Regression has the harped strands 1 in off at the ends (13 mismatches), by design. This replaces decision 6's "a cloned entry keeps the clone's grid".

**Guidance (2026-10-15):**
- **No configuration is assumed.** The import must work with any PGSuper configuration and publisher (WSDOT, TxDOT, California, ...): it compares with whatever girder entries and beam factories the project's library has, and has no agency's girder names or shapes in code. WSDOT (Local) is only the test library (the Regression configuration has two girders).
- **Names are clues, not data.** A girder name often hints at the height (e.g. "WF100G", "NU2000"), but not always. The height in a name may help rank candidates; it never replaces the section.
- **I-beams first, nothing I-beam only.** PGSuper also has slabs, voided slabs, box beams, U-beams, decked and double tees, and others. The work is tested on I-beams (the most common) but must not prevent the other families: the section cut, comparison, and fitting are shape-generic (voids, composite sections, any factory's dimensions). Shape-specific rules (the generated strand grid of step 4, the bar roles of G3) are per family with a generic fallback (e.g. straight strands on a grid in the bottom of any section, clear of voids), so families can be added without changing the rest. To check when other families come up: the girder layout finds girder ends from the top surface outline (U-beams have two top flanges).
- **I-beams are typically WF or NU types.** NU girders have rounded corners at the web-flange junctions (fillets); the section comparison sees them as width changes near the flanges, and the family's fillet dimension is fitted like the others.

## How dimensions are fitted (G1 steps 3-4, 2026-10-15)

Every girder factory makes its section the same way (WBFL GeomModel, e.g. `PrecastBeam` for the WF factory, `NUBeam` for the NU factories; GenericBridge wraps them): the left half of the outline is a chain of corners (and fillet arcs) from the top center to the bottom center, whose positions are sums of the dimensions, mirrored about the vertical axis. The fitter (`OutlineDimensionFitter` in `GirderSection.cpp`) uses that, without knowing any family:
- **Start:** the candidate's dimensions with the height solved. A factory with a library entry is fitted from both the entry's dimensions and the factory's defaults, and the better fit is kept (2026-10-16): a fit settles near its start, so the WSDOT Precast I-Beam started from a 58 in entry missed the 31 in PennDOT fit that its own defaults (and the BridgeSight wrapper of the same factory) found. A fit from either start clones the entry. The dimensions that change the height the most are tried in turn (the height of an NU girder is the sum of its depths: the web's can reach a height the flanges' can't).
- **Symmetry:** only the left halves of the outlines are compared, so an edge is never paired with the mirror image of an edge.
- **Two least-squares fits, the better kept:**
  - corners first: the corners of the model's half outline (Douglas-Peucker at the length tolerance) paired in order with the factory's corners, then the outline distance. For polygonal sections (WF), the corners fix the dimensions.
  - the outline distance alone, for curved sections (NU fillets), whose corners are many short edges that pair poorly.
- **Outline distance:** points spaced evenly along each half outline measured to the nearest edge of the other, both ways. Comparing widths depth by depth (the first fitter) holds a fit where a flange edge lines up with the wrong edge (e.g. a 34 in wide "web" standing in for a bottom bulb).
- **Zero dimensions:** many factory defaults have zero lengths (a chamfer, a flange's horizontal part). Least-squares steps stop at zero instead of being rejected (otherwise the fit never moves). A factory's outline leaves out the corners of a zero dimension, and one zero dimension can hide another (with no horizontal offset, the lower part of a WF web is just more web, so its depth changes nothing), so the dimensions to fit are found with every zero length at 1 in, and the corner fit starts there.
- **Sanity:** a factory may accept dimensions whose section crosses itself or has non-finite coordinates (one hung an import); such sections are rejected. In debug builds WBFL checks some of these in the shape code (e.g. `Polygon::AddPoint` asserts finite points; a U-beam with vertical webs gets non-finite points), so while a trial section is made the checks report to the debugger output instead of stopping (`QuietChecks` in `GirderSection.cpp`).
- **Constant depth (2026-10-16):** the section is fitted at one place, so a dimension that only shapes the far end (a tapered I-beam's `D7_End`) keeps its start value. `MakeConstantDepth` solves such a dimension so `BeamFactory::GetBeamHeight` is the same at both ends without changing the fitted section, or rejects the fit. (`BeamFactory::IsPrismatic` isn't the test: it means constant section properties, so end blocks make a box beam non-prismatic, and the bulb tee factory always says no.)
- **Acceptance** (`FittedSectionComparison`): the tolerances of decision 2. Widths may differ over 1/2 in of depth in total (local details), or over any depth where the outlines are within 1/8 in of each other: on a shallow taper a small vertical offset is a large width difference (the PennDOT sections are about 0.09 in lower on one side, from the crown).
- **Same girder, no new fit (2026-10-16):** a girder whose section is the same as one already matched (`SameGirderSection`: outlines within 1/8 in, area and I within 1%, centroid within 1/8 in, which allows for sections cut from geometry, e.g. PennDOT's girders sit 0.1 in lower on one side) uses that girder's library entry without fitting. The girders of a span, and often of the bridge, are one shape: PennDOT's five girders get one entry (was three, "(2)" and "(3)") and the import takes half the time.
- **Girder type groups:** when girders have different entries, each span's girders are regrouped so adjacent girders with the same entry form a group (PGSuper asserts that a group's girders share an entry). The round trip `PGSuper_ExteriorGirders.pgs` (the main model with a 6 in wider top flange on span 1's exterior girders, entry "WF66G Exterior") tests it.
- **Progress:** fitting takes seconds per girder (longer in debug builds), so the progress window names each screening and full fit.
- **Choice:** an accepted fit of a library entry (cloned) is preferred to an accepted fit of a factory's defaults (a new entry), since the clone keeps the entry's other data. Fits from other factories as good as the chosen one (accepted, the same kind, outlines within the length tolerance of each other) are ties: a factory made for girders of constant depth (`HasConstantDepthDefaults`) is preferred to one made for varying depth (e.g. the Tapered I-Beam ties with the Precast I-Beam only by being fitted with a constant depth); then the user picks the factory (WBFL MFCTools `AfxChoose`, asked once per set of factories in an import), and an import without a user (the command line) takes the one published by WSDOT (`BeamFactory::GetPublisher`), else the best fit. The log names the ties and how one was chosen.

Results (WSDOT (Local) library): the three PennDOT girders fit the WF factory (W32BTG cloned) with the dimensions read from the outline (D1 4.52, D2 1.28, D3 2, D4 3.5, D5 3, D6 9, T1 = T2 = 8, W2 15, W3 2, W4 3.5, W5 9, C1 0.77 in). The Iowa girders fit the NU factories (W1 34, W2 30, T 6.51, R1 8.2, R2 8.4 in; outlines within 0.08 in); the library has no NU entry, so a new entry is needed (step 4). Several factories make the same shape (the Nebraska NU Girder (WSDOT), NU Beam and CalTrans I-Beam (BridgeSight) factories all use `NUBeam`): a command line import takes the Nebraska NU Girder.

## Open questions

None. The draft's four questions were answered on 2026-10-15:
1. **Library entries:** allowed; the existing library API turned out to be enough (Decision 1, "Library entry creation").
2. **Strand definition:** the library's grid when every strand is on it, otherwise individual strands (`sdtDirectStrandInput`) (Decision 2, and decision 6 of "Library entry creation").
3. **Bar roles in mapping tables:** tabled (Decision 3).
4. **Longitudinal girder bars:** imported (`CLongitudinalRebarData`, Decision 4).
