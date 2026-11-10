# Stirrup Plan: Bars, Layouts, and General IFC Export

Draft, 2026-11-10. Assessment and plan only; nothing is implemented. Goals (the user): (a) stirrup input that matches actual bars and layouts, (b) bar geometry for IFC, (c) the bars used in PGSuper's analysis, (d) PGSuper's shear design mapped to the new input, (e) a general IFC export instead of the WSDOT-specific one.

## Decisions (the user, 2026-11-10)

- **A plug-in, not PGSuper core**, as far as possible. Core changes are limited to extension points: a plug-in can add pages to the girder, segment, and closure joint dialogs and hide the default Transverse Reinforcement page, and the same at the girder library level, in both PGSuper and the PGSLibraryEditor.
- **Bars PGSuper doesn't analyze** (e.g. the WSDOT G3 top flange bars, end block bars): the plug-in gives a UI to input and describe them, and exports them to IFC. The importer should capture them from a model too, so a plain export (not only a round trip) writes them.
- **Orientation:** normal to the girder, parallel to the skew, or splayed, set per spacing zone.
- **Design (d):** first version below; other options are noted for later.

## Assessment

### PGSuper's model

`CShearData2` (one per segment and per closure joint) is zones of counts: primary zones with a length, spacing, bar size, `nVertBars` (legs, a `Float64`), `nHorzInterfaceBars`, and `ConfinementBarSize`; separate horizontal interface zones; an additional splitting zone and an additional confinement zone; `bAreZonesSymmetrical`; one bar material. It doesn't say what a bar looks like or how bars are grouped ("2 legs" is one U bar or two L bars) and has no orientation. The girder library entry has a `CShearData2` too, a seed value.

The analysis needs only scalars, all through `IStirrupGeometry` (BridgeAgent): `GetVertStirrupAvs` and `GetAlpha` (always 90 degrees) for shear capacity, the horizontal interface Avs and counts, `GetSplittingAv`, `GetConfinementAv`, and the bar diameter for the detailing checks. The shear design works on a copy (`STIRRUPCONFIG::SHEARZONEDATA` in `GDRCONFIG`), searching bar size, legs, and spacing; its result is a `CShearData2` on `pgsSegmentDesignArtifact`, copied into the girder when the design is accepted. Other readers: the reports, the girder elevation views, the WSDOT, TxDOT, and KDOT agents, and this plug-in's import (`Stirrups.cpp`).

### The exporter

`CreateStirrups` (`IfcExporter.cpp`) writes WSDOT details for WF girders only: one G2 hairpin per stirrup place whatever the number of legs, its legs always into the deck (`H1 = Hg + A + 3"`, from the WSDOT standard drawing); G3 top flange bars, which aren't PGSuper data; G9 confinement bars with a fixed 9.125 in bulb height (PGSuper doesn't expose it); G10 bars; and the WSDOT splay at skewed ends (`fn_bar_angle`, a hard-coded 10 ft transition). What's already general: a bar is its centerline path (`BarPath`/`BarDefinition`, `BarShape.h`), and its IFC type, length, mass, and usBridge properties follow from it (`GetBarType`, `Rebar.h`).

### The round trip

A girder whose source model has bars keeps all of them, and PGSuper's aren't written (`RoundTripDesign.md`). So bars PGSuper doesn't analyze survive a round trip today, but so do the model's stirrups after the user changes them in PGSuper (logged, not updated). A plain export loses the bars PGSuper doesn't model.

### What's there for plug-ins

- `IExtendUI` (`IFace/ExtendUI.h`): `IEditGirderCallback`, `IEditSegmentCallback`, `IEditClosureJointCallback`, and `IEditSplicedGirderCallback` add a page to the girder dialogs, with `OnOK` returning a transaction (undo/redo), and `ICopyGirderPropertiesCallback` takes part in Copy Girder Properties. A plug-in can't hide or replace a default page.
- `IAgentPersist`: the plug-in stores its data in the project (as the georeferencing is).
- The library has no plug-ins at the entry level. `LibraryMgrPluginExample` is a document plug-in for the library editor (menus and views), not an extension of a library entry or its dialog (`GirderMainSheet`, `psgLib`).

## Approach: the plug-in owns the bars, `CShearData2` is computed from them

The plug-in's stirrup layout is the input; on every change it writes the `CShearData2` that gives the same analysis results. PGSuper's analysis, design, reports, and other agents keep reading `CShearData2` unchanged, so (c) needs no analysis changes in the first version.

### The plug-in's model

1. **Bar shape:** a centerline path of lines and arcs, optionally with its ACI/CRSI bend shape code and dimensions. A shape's dimensions are expressions of named anchors of the girder section (top of girder, web faces, top of the bottom flange, cover, haunch depth "A", top of deck, bar diameter, bend radius), so one shape fits every girder of a family and follows variable depth and end blocks. Anchors come from the beam factory's named dimensions (the library entry's dimensions, read by name; no core change) and, for the generic fallback, from the section's shape (`IShapes`). Each leg has a role from the shape's geometry: a shear leg (crosses the web at mid-depth), into the deck (horizontal interface), around the bottom flange (confinement), or none.
2. **Bar set:** the bars at one stirrup place: shapes with a size, a position across the girder, and mirroring (e.g. one G2 hairpin, a G9, and a G10; two mirrored L bars). This ends the "legs" ambiguity.
3. **Zones:** the familiar table (length, spacing, symmetric or not), each zone with a bar set and an orientation (normal, parallel to the skew, or splayed from one end's angle to the other's across the zone).
4. **Detailing bars:** bars PGSuper doesn't analyze, placed with the stirrups (at every place of chosen zones), at stations, or along the girder. They're exported and never counted in `CShearData2`.
5. **Catalogs as data:** agency details (WSDOT G1 to G10, TxDOT, Iowa, ...) are catalog data, not code, and a girder library entry has its default bar sets and zones (a seed value, as the library's reinforcement is now). Each girder family has a generic fallback from the section: a hairpin per web (I-beams one, U-beams and boxes two), a confinement bar following the bottom flange outline offset by the cover. No agency's details are assumed.

### Computing `CShearData2` (c)

Per zone: the bar size of the shear legs and `nVertBars` from their area (fractional legs when a place has bars of different sizes: `nVertBars` is a `Float64`); `nHorzInterfaceBars` from the legs into the deck; `ConfinementBarSize` from the bar enclosing the bottom strands (the import's geometric test, `Stirrups.cpp`); `bUsePrimaryForSplitting` (the stirrups near the ends are the splitting bars, as in the import); `bAreZonesSymmetrical` and the material from the layout. Spacing is along the girder for every orientation (to confirm against LRFD for splayed bars). Not expressible without core changes, so later: inclined stirrups (`GetAlpha`), anchorage of the legs into the deck, checks from the real geometry (bend radius, cover, strand clearance, LRFD 5.10.8 anchorage).

### Keeping the two consistent

`CShearData2` also changes outside the plug-in: an accepted design, Copy Girder Properties, a change of girder type, the default page when it isn't hidden, a file saved by a PGSuper without the plug-in (the broker keeps the plug-in's data unchanged then; see "Core changes", item 5). The plug-in records a fingerprint of the `CShearData2` it wrote; when the bridge description events show a different one, it reconciles (below) and reports it.

### Design (d), first version

Design stays on bar size, legs, and spacing (`ShearDesignTool` unchanged). When an accepted design changes a girder's `CShearData2`, the plug-in assigns each zone a bar set from the girder's catalog that gives its legs, legs into the deck, and confinement bar, keeping the zone's orientation; with none, the generic bar set of the family, reported. Detailing bars stay. **Later options:** the design searching the catalog's bar sets (ordered by area) instead of size and leg pairs, which needs a design extension point in core; designing the orientation or the detailing bars.

### IFC export (e)

One general loop: each zone, each stirrup place, each bar of the bar set, its shape at the section there, oriented by the zone, placed through `GetBarType`; then the detailing bars. Bar use, position, shape code, and mark come from the shape data. Projects without the plug-in's data (and older projects) export the generic family bars from their `CShearData2`, so every family gets stirrups. `CreateStirrups`' WSDOT code goes: the WSDOT details become a catalog.

### Import and round trips

The import keeps the model's actual bars: each stirrup place's bars matched to catalog shapes, else a shape made from the path, and the zones from the places (as now); the bars that aren't stirrups or confinement become detailing bars, so a plain export writes them. The round trip can then write PGSuper's bars in place of the model's stirrups (pairing them through the registry) while the detailing bars keep their GlobalIds, closing the stale-stirrup case and the open items of G3 in `GirderImportPlan.md` (other leg counts and splitting bars that don't round-trip).

## Core changes (PGSuper and psgLib)

1. **Hide a default page:** a girder, segment, or closure joint dialog callback says which default pages it replaces (the Transverse Reinforcement page; its horizontal interface zones go with it). An ID per default page, or a new callback interface so existing plug-ins aren't affected.
2. **Library entry extension data:** a girder library entry keeps plug-in data by plug-in ID, saved in the library and project files, kept as is when the plug-in isn't loaded (shared master libraries, a PGSLibraryEditor without it), copied with the entry, and part of the entry comparison (library conflicts on opening a project or with the import template: the plug-in compares its data, else byte for byte).
3. **Library entry dialog extension:** a component category that `psgLib` enumerates, so the girder entry dialog has the plug-in's pages (and hides the default one) in both PGSuper and the PGSLibraryEditor. The library part of the plug-in must not need the broker (there's none in the library editor), so it's a separate component in the DLL.
4. **Seeding:** where PGSuper copies the library's transverse reinforcement into a girder (a new girder, a change of girder type), the plug-in copies its library data too: an event or callback at those places, unless the bridge description events are enough (to check).
5. To check: whether a hint tells an accepted design from other girder changes.

**A project saved without the plug-in (the user, 2026-11-10):** no core change is needed. The WBFL::EAF broker keeps the data unit of an agent that isn't loaded as an XML blob and writes it back unchanged when the project is saved. The project may have been edited in the meantime, so when the plug-in loads its data it must check it against PGSuper's (the girders, segments, and closure joints that exist, their `CShearData2`, girder types, and sections) and reconcile: data for elements that are gone is dropped, elements without data get the generic layout from their `CShearData2`, and a changed `CShearData2` is reconciled as described in "Keeping the two consistent". Each step is reported.

## Order

1. **Export first, no core change:** move the bar geometry into a shared module, shapes on section anchors, the WSDOT details as a catalog, the generic family fallback, all from the current `CShearData2`. Removes the WSDOT code; every family gets stirrups; the export baselines show the change.
2. **Core extension points** (1 to 3 above), uncommitted until the user says so, with the regression tests.
3. **The plug-in's girder pages:** zones, bar sets, orientation, detailing bars, with a section view of the bars; persistence; `CShearData2` computed from them; reconciling.
4. **The library pages and seeding.**
5. **Design mapping** (first version).
6. **Import of the actual bars** and detailing bars; round trip with PGSuper's stirrups replacing the model's.
7. **PGSplice:** segments and closure joints (their callbacks exist); other girder families.

## Open questions

- Where the catalog lives: in the library entry's extension data only, or also as a shared catalog file (agency details reused across entries)?
- The bar mark: data in the catalog, or generated?
- How a girder's layout is copied between girders and girder types when their sections differ (anchors make shapes portable; zones may not fit).
