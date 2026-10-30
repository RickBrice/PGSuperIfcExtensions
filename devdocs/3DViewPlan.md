# 3D View Plan

Drafted 2026-10-30 as an exploratory plan (the user: no programming yet). A 3D view of the bridge model in PGSuper and PGSplice.

## What exists (2026-10-30)

- PGSuper's bridge views (plan, section, girder elevation: `PGSuperAppPlugin`, e.g. `BridgePlanView`, `BridgeSectionView`) are 2D, drawn with WBFL's display object framework (DManip). PGSuper, WBFL, and BridgeLink don't use OpenGL or OpenCASCADE.
- OpenCASCADE comes with this plugin only, through IfcOpenShell's geometry kernel (`ifcopenshell_geometry_kernel_opencascade`, the `TK*` libraries in `PGSuperIfcExtensions.vcxproj`). The visualization libraries (`TKV3d`, `TKService`) are linked; the OpenGL driver (`TKOpenGl`) isn't.
- The IFC exporter already turns PGSuper data into geometry: girder sections from `IShapes::GetSegmentShape` at points of interest with the skew cut angle (`CreateSectionProfile`, `GirderSegment_SectionedSolidHorizontal`), the deck and haunch, barriers, piers, bearings, strands, and bars. The import meshes the model with the IfcOpenShell iterator (see [ImportGeometry.md](ImportGeometry.md)).

## The two approaches

1. **Native:** PGSuper's own data (broker interfaces, `WBFL::GenericBridge`) turned into geometry and shown by a viewer independent of this plugin.
2. **IFC:** the project exported to IFC (or an imported IFC model) meshed by IfcOpenShell with OpenCASCADE and shown.

| | Native | IFC |
|---|---|---|
| Updates with edits | Yes: it listens to the broker like the other views | A snapshot: re-export after each change (geometry takes ~1.5 s per beam today, [IfcImportPlan.md](IfcImportPlan.md)) |
| Selection | Maps to girder, segment, and pier keys, so it can open PGSuper's editors | Maps to IFC entities; back to PGSuper items needs the planned GlobalId mapping |
| Stages, PGSplice, results | The model at any interval (erection, closure joints, tendons, temporary supports); colors by stress or deflection | What the exporter writes: no intervals or results |
| Dependencies | None on this plugin | Requires this plugin |
| Geometry code | New, but thin: PGSuper's interfaces do the hard parts (sections at POIs, skew, strand and tendon geometry) | Exists, validated by the round trips |
| Robustness | Can be triangulated directly: consecutive section polygons joined into triangle strips, no solid modeling, nothing to fail | Inherits OpenCASCADE failures (the deck slab of the `PGSuper_Skew_Straight_AlongPier` export is dropped without a message) |
| Other models (Iowa, PennDOT) | No | Yes: useful for checking imports |

## Recommendation

The native approach, packaged as its own BridgeLink extension, with a viewer that takes triangle meshes tagged with object ids, so this plugin can later feed it IFC models as a second source.

- A PGSuper view must update with edits, select into PGSuper's editors, follow the construction stages, and work without optional extensions. The IFC route shows a snapshot of a lossy export.
- The native route needs no OpenCASCADE for geometry: a segment is a series of section polygons (the same `IShapes` sections the exporter writes), joined into strips. OpenCASCADE is only a candidate renderer.
- The IFC route's real advantage, viewing other people's models, is added on top of a mesh viewer (step 5 below) instead of being the foundation.
- Duplication with the exporter stays small since both read the same PGSuper interfaces. Later, the geometry step can become a shared intermediate model (sections placed along girder lines) that the viewer and the IFC writer both use.

## Plan

1. **Spikes before committing:**
   - OpenCASCADE's viewer (`V3d`/`AIS` with triangulations, `TKOpenGl`) in a BridgeLink MDI child frame: OpenGL context, DPI, docking.
   - Installer size and license: OpenCASCADE is LGPL 2.1 with an exception, probably fine for dynamic linking, but WSDOT decides what the installer ships.
   - Speed with the round-trip model's ~2,400 bars and ~400 strands (instancing).
   - A small Direct3D or OpenGL mesh renderer for comparison, in case OpenCASCADE's viewer is too heavy to ship.
2. **Architecture:** a separate extension DLL with a 3D view, so OpenCASCADE (if chosen) stays out of PGSuper's core. Three layers:
   - **Geometry from PGSuper:** girder sections at points of interest, the deck and haunch, barriers, piers, bearings, strands and tendons, bars.
   - **Triangulation:** meshes from the section polygons.
   - **Viewer:** meshes with an object id, a category (girder, deck, strand, ...), and the stage they exist in.
3. **Increments:**
   1. Girders only, in a fixed stage.
   2. Deck, haunch, and barriers.
   3. Piers and bearings.
   4. Strands, tendons, and bars, each switched on and off.
   5. Selection opens PGSuper's editors.
   6. An interval slider for the construction stages (PGSplice).
   7. Coloring by analysis results.
4. **In this plugin, later:** an IFC source that feeds the IfcOpenShell iterator's meshes (already computed by the import) to the same viewer, to look at imported and third-party models and to compare a model with its re-export.

## Open questions

- Will WSDOT ship OpenCASCADE with PGSuper?
- Is the view for PGSuper and PGSplice only, or a BridgeLink service other applications could use?
