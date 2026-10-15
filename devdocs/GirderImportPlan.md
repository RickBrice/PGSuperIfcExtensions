# Girder Import Plan: Section, Strands, Transverse Reinforcement

Drafted 2026-10-15. Waiting for review. Priorities 2 and 3 in `CLAUDE.md` (girder mapping with the PGSuper library, information from geometry). It builds on Phase 5 (girder section) and Phase 6 (prestressing) of [IfcImportPlan.md](IfcImportPlan.md).

Today the importer sets the girder type by name only (`GetGirderLibraryEntry`, falling back to the first I-beam) and leaves the strands, stirrups, confinement, and splitting bars of the import template unchanged. An imported girder has the template's reinforcement.

## What the test models have (2026-10-15)

| | PGSuper export | Iowa US59 | PennDOT |
|---|---|---|---|
| Girder geometry | `IfcSectionedSolidHorizontal` of the exact PGSuper sections | `IfcMappedItem` of swept sections | `IfcFacetedBrep` |
| Girder type name | `IfcBeamType` = library entry (WF66G, WF100G) | `IaDOT_PPCB.2_Type` (BTB45) | `_PS Concrete Beams.Type` |
| Strands | 369 `IfcTendon`, `IfcSweptDiskSolid` along a polyline (start, harp points, end), types with debonding (`usBrPset_TendonDebonding`) | 661 `IfcTendon` STRAND in "Strands_Beam n_Span m" assemblies, mapped `IfcExtrudedAreaSolid` (straight, apparently scaled by the mapping), `IaDOT_Strand` properties as text (`0.5"`, `270`, `Low Relaxation`); 81 more are lifting loops (NOTDEFINED). No harping or debonding words in the file | none |
| Girder bars | 2354 `IfcReinforcingBar` with types: G2 stirrups, G3 top bars, G9/G10 confinement, longitudinal; `usBrPset_ACI_ReinforcingBar.BarUse` (STIRRUP, TIE, ...) | 6080 bars in "Rebar_Beam n_Span m" assemblies, mapped extrusions, ProStructures marks and ACI bend shapes (`4d1` IA56, `4e1` S10, `4h1` S11, `6b3` 20, `4c1` straight), lengths as text (`2:3`) | none |

## PGSuper data the import fills

- **Section:** the girder library entry of each girder (`GirderLibraryEntry`, a beam factory with named dimensions). Creating entries has no programmatic path today (Phase 5).
- **Strands:** `CStrandData`. With `pgsTypes::sdtDirectStrandInput` (available in PGSuper and PGSplice) each `CStrandRow` is one strand: its type (straight, harped, temporary), its horizontal position `m_Z`, its elevation at the start, the two harp points, and the end (`m_Y`, from the top or bottom face), its debonding at each end, and extension. The harp points are on `CStrandData` (`SetHarpPoints`). Strand material per type, and jacking.
  The other definition types fill the library entry's strand grid (`sdtTotal`, `sdtStraightHarped`, `sdtDirectSelection`).
- **Transverse reinforcement:** `CShearData2`: shear zones (length, bar size, spacing, number of legs, horizontal interface bars), the splitting zone (bar size, spacing, length, number of bars), the confinement zone (bar size, spacing, length), symmetric or not.

## G1. Girder section

Phase 5 of [IfcImportPlan.md](IfcImportPlan.md) is the plan; in short:
1. Get the section polygon of each girder from the model: the swept profiles (PGSuper export, Iowa), or a cut through the solid (PennDOT). Several cuts to see end blocks and tapers; the prismatic middle portion is the section.
2. Match a library entry by geometry (depth, flange widths and thicknesses, area, I, centroid), confirming a name match. A geometric match with another name is used and logged.
3. If nothing matches: classify the family and fit its dimensions (minimize the difference between the polygon and the factory's section).
4. Create the library entry (clone the closest entry of the family, set the dimensions). This needs a new PGSuper interface (open question 1).

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

## Library entry creation (details to agree before coding)

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
4. **Names:** the IFC type name (e.g. "BTB45"), with " (2)", and so on, when the name is already used by a different girder.
5. **Generated strand grid:** 2 in spacing with 2 in cover to the strand center. Straight strands in the bottom flange, harped strands in the web, temporary strands in the top flange. No strands are debondable (harping is assumed, not debonding; harping and debonding could become a configuration option later).
6. **The user's intent isn't known** when a project is created from an IFC file, so:
   - **The girder has strands:** if they're all on the grid of the girder's library entry, a grid fill; otherwise individual strands.
   - **The girder has no strands:** the library entry gets a strand grid if it doesn't have one.
   - **A cloned entry** keeps the clone's grid. **A new entry** gets a generated grid (2 in spacing) even when the girder's strands are imported as individual strands: designers will move to the grid input, and it needs a grid.
7. **Where entries go:** the project's own library (saved with the project), never the master library.

## Open questions

1. **Library entries:** may we add a PGSuper interface to create a girder library entry (clone an entry, set its dimensions)? (Phase 5 questions 1-4: cloning the closest entry, naming, preferring a geometric match over creating one.)
2. **Strand definition:** always `sdtDirectStrandInput` (exact positions, but the strands can't be adjusted with the library's strand grid in PGSuper), or the library grid when every strand is on it (familiar input for PGSuper users)?
3. **Bar roles:** should the agency bar roles (Iowa's marks and bend shapes) go in the agency mapping tables, as the girder type does? This would mean a table section for bar roles.
4. **Longitudinal girder bars** (the PGSuper export's longitudinal rows, Iowa's `5b1`, `6b3`): import them too (`CLongitudinalRebarData`)?
