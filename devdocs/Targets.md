# Mapping table targets and the PGSuper data behind them

A **target** is a PGSuper data item with a fixed, readable name such as `girder.fci` or `bearing.fixed_x`. Mapping tables, the Mapping Table Editor and the user documentation (`Documentation/MappingTableEditor/src/targets.md`) use only these names. This note is for developers: it says which PGSuper data each name stands for, and where the C++ connects them.

## Where the connection is made in the code

**One list, in code.** Every target is a `TargetDef` in `GetTargetDefs()` (`IfcTargets.cpp`; the struct is in `IfcTargets.h`):

| Field | Meaning |
|---|---|
| `name` | The name the tables use, e.g. `"girder.fci"` |
| `element` | The element role it belongs to (`ElementKind`). An abutment uses the pier targets (`GetTargetElement`) |
| `kind` | `ValueKind`: Stress, Length, Angle, Count, Boolean, Text, TextList, and so on. It decides the unit conversion and the IFC measure type |
| `description` | Shown in the editor and in log messages |
| `display_unit`, `display_round` | Export only: the display unit and rounding used when properties are written in display units |
| `get` | Export only: a lambda that reads the value from PGSuper, given an `ExportContext` (the broker, export options, and the segment, pier, strand, bar, bearing, or barrier being written). Empty for an import-only target |

A table can't add targets. `CIfcMappingTable` resolves each name it reads with `FindTargetDef` (`IfcMappingTable.cpp`, when properties and locations are loaded), so an unknown name is a table error. The editor (`MappingEditorPanes.cpp`, `MappingEditorModelPane.cpp`), the model browser (`IfcModelBrowser.cpp`) and the IDS writer (`IfcTableIds.cpp`) look names up the same way.

**Export: the getter is the mapping.** The exporter (`IfcExporter.cpp`) fills an `ExportContext` for each element it writes. For example, for a bearing it sets `context.bearing` to the pier's `CBearingData2` and `context.bearing_elevation` to the `IBridge::GetBearingElevationDetails` result. It then calls `WritePropertySets` (`IfcPropertyWriter.h`). For each property the table declares for the element's role, the writer calls `property.target->get(context)` (`IfcPropertyWriter.h`, and `GetTypePropertySignature` in `IfcPropertyWriter.cpp`). That converts to display units if asked (`ConvertToDisplayUnits`) and writes the IFC property. The getter returns `ExportValue`: a value, a list, `NoValue()` (write the property empty) or `Absent()` (leave it out). So for export, the code of the getter in `IfcTargets.cpp` is the definitive answer to "what PGSuper data is this target?".

**Import: the call site is the mapping.** There are no setters in `TargetDef`. Import code asks for a target by name:

```cpp
const auto& reader = CIfcImporter::GetTargetReader();
if (auto reading = reader.Read("girder.fci", beam))     // or ReadAll for a TextList target
   segment.Material.Concrete.Fci = std::get<Float64>(reading->value);  // already in system units (SI)
else
   reader.ReportNotFound("girder.fci", beam, beam_name);  // logs the locations tried and hints
```

`CIfcTargetReader` (`IfcTargetReader.h/.cpp`) tries the target's locations from the table in order. It converts the first value found to PGSuper system units (with the property's own unit, else the IFC project unit, else the unit the table states) and maps text values through the location's `map`. It returns a `TargetReading` with the value and its source, for the log. The code that called it decides what to do with the value: set PGSuper input data, cross-check geometry, or pick a library entry. So for import, find the target's name in the import code (the call sites are listed below).

The import writes **input data**, mostly members of `CBridgeDescription2` and its parts, which are then set with `IBridgeDescription::SetBridgeDescription`. The export usually reads **results or computed values** through broker interfaces such as `IMaterials`, `IBridge` and `ICamber`. So a target's export and import may touch different PGSuper objects for the same quantity. `girder.fci` is an example: the export reads `IMaterials::GetSegmentFc` at release, and the import sets the segment's input `Fci`.

### Adding a target

1. Add a `TargetDef` to `GetTargetDefs()` (`IfcTargets.cpp`), with a getter if it's exported.
2. If it's imported, call `reader.Read`/`ReadAll` with its name where the importer has the element, and set the PGSuper data.
3. Give it locations or properties in `MappingTables/Standard.ifcmap` (and in the agency tables in `Tests/ImportValidation/mappings` if they need it).
4. Add a row to the tables below.

## The targets

All values are in PGSuper system units (SI). "Export" is what the `get` lambda in `IfcTargets.cpp` reads; "Import" is what the import does with the value. A dash means the target isn't used in that direction.

### Project and bridge

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `project.job_number` | Text | `IProjectProperties::GetJobNumber` (empty: no value) | – |
| `project.export_date` | Text | today's date, YYYY-MM-DD (not PGSuper data) | – |
| `bridge.number_of_spans` | Count | `IBridge::GetSpanCount` | Cross check only (`CIfcBridgeImporter::Import`, `IfcBridgeImporter.cpp`): the number of spans comes from the modeled piers. A different value is logged as an error |
| `bridge.length` | Length | `IBridge::GetLength` | – |
| `bridge.roadway_width` | Length | `IBridge::GetCurbToCurbWidth(0.0)` (at the start of the bridge) | – |
| `bridge.start_station` | Length | `IBridge::GetPierStation(0)` | – |
| `bridge.end_station` | Length | `IBridge::GetPierStation(nSpans)` | – |
| `bridge.max_skew` | Angle | the largest `IBridge::GetPierSkew`, last pier excluded | – |
| `bridge.name` | Text | `IProjectProperties::GetBridgeName` | – |
| `bridge.id` | Text | `IProjectProperties::GetBridgeID` | – |
| `bridge.design_live_load` | Text | the first of `ILiveLoads::GetLiveLoadNames(pgsTypes::lltDesign)` (absent if there is no design live load) | – |

### Piers and abutments

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `pier.station` | Length | `IBridge::GetPierStation(pier)` | – (pier stations come from `Pset_Stationing` on the pier's `IfcReferent`, `Piers.cpp`; that is IFC structure, not a table target) |
| `pier.skew` | Angle | `IBridge::GetPierSkew(pier)` | – (the skew is measured from referents, bearings, and girder ends, `set_pier_orientation` in `Piers.cpp`) |

### Girders (one `IfcBeam` per segment)

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `girder.designation` | Text | `SEGMENT_LABEL(segmentKey)` with numeric labels, e.g. "Span 1, Girder 2" | `get_girder_key` (`BeamSpacing.cpp`): parsed into a `CGirderKey` (group and girder index) for the girder layout. Otherwise the beam's Name is tried, then the layout's geometry |
| `girder.type_names` | TextList | `CSplicedGirderData::GetGirderName()` (the girder library entry name, from `IBridgeDescription::GetGirder`) | `GetGirderTypeNames` (`IfcBridgeImporter.cpp`): every name found, generic words dropped. They are clues for `GetGirderLibraryEntry`, which picks or creates the `GirderLibraryEntry` set with `CBridgeDescription2::SetGirderLibraryEntry` or `CSplicedGirderData::SetGirderLibraryEntry`. The section geometry decides, and the names break ties and name new entries |
| `girder.family_name` | Text | `CBridgeDescription2::GetGirderFamilyName()` | – (the family comes from the chosen library entry) |
| `girder.fc` | Stress | `IMaterials::GetSegmentFc28(segmentKey)` | `SetGirderProperties` (`IfcBridgeImporter.cpp`): `CPrecastSegmentData::Material.Concrete.Fc` (segment 0 of each girder) |
| `girder.fci` | Stress | `IMaterials::GetSegmentFc(segmentKey, IIntervals::GetPrestressReleaseInterval(segmentKey))` | `SetGirderProperties`: `CPrecastSegmentData::Material.Concrete.Fci` |
| `girder.strength_class` | Text | `GetSegmentFc28` formatted in the display stress unit | – |
| `girder.max_aggregate_size` | Length | `IMaterials::GetSegmentMaxAggrSize` | – |
| `girder.assembly_place` | Text | – (written as a fixed value by the table) | `IsValidBridge` (`IfcBridgeImporter.cpp`): with `girder.casting_method`, checks that the beams are factory-made precast girders ("FACTORY", "PRECAST") |
| `girder.casting_method` | Text | – (fixed value in the table) | See `girder.assembly_place` |
| `girder.span` | Length | `IBridge::GetSegmentSpanLength` | – |
| `girder.slope` | Angle | `atan(IBridge::GetSegmentSlope)` | – |
| `girder.roll` | Angle | `atan(IGirder::GetOrientation)` | – |
| `girder.batter` | Angle | `atan(IBridge::GetSegmentSlope)` when the export battered the girder ends (`CIfcExportOptions::batter_ends`), else 0 | – |
| `girder.bunk_point` | Length | the larger of `CPrecastSegmentData::HandlingData.LeadingSupportPoint` and `TrailingSupportPoint` | – |
| `girder.jacking_stress` | Stress | `IStrandGeometry::GetJackingStress(segmentKey, pgsTypes::Permanent)` | – |
| `girder.camber_ratio` | Ratio | at the released mid-span POI: `IProductForces::GetDeflection` for pretension plus girder self-weight at release, plus `IGirder::GetPrecamber`, divided by `IBridge::GetSegmentPlanLength`. No value unless the export includes camber | – |
| `girder.camber_at_release` | Length | `ICamber::GetInitialCamber` at the released mid-span POI | – |
| `girder.camber_after_losses` | Length | `ICamber::GetExcessCamber(poi, CreepTime::Max)` | – |
| `girder.screed_camber` | Length | `ICamber::GetScreedCamber(poi, CreepTime::Max)` | – |
| `girder.length` | Length | `IBridge::GetSegmentPlanLength` | – |
| `girder.cross_section_area` | Area | `ISectionProperties::GetAg` at release, released mid-span POI | – |
| `girder.outer_surface_area` | Area | `GetSegmentPlanLength` × `ISectionProperties::GetPerimeter` | – |
| `girder.gross_weight` | Mass | `ISectionProperties::GetSegmentWeight` / g | – |

### Deck

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `deck.gross_depth` | Length | – | `ImportSlab` (`IfcBridgeImporter.cpp`): cross check of the depth measured from the deck geometry, which wins (a difference over 1/8 in is logged). It sets `CDeckDescription2::GrossDepth` only when the depth can't be measured |
| `deck.depth` | Length | `CDeckDescription2::GrossDepth` (input) | – |
| `deck.fc` | Stress | `IMaterials::GetDeckFc28` | – |
| `deck.strength_class` | Text | `GetDeckFc28` formatted in the display stress unit | – |
| `deck.length` | Length | `IBridge::GetLength` | – |
| `deck.perimeter` | Length | perimeter of `IBridge::GetSlabPerimeter` (50 points along each edge) | – |
| `deck.gross_area` | Area | area of the same outline | – |
| `deck.max_aggregate_size` | Length | `IMaterials::GetDeckMaxAggrSize` | – |

### Traffic barriers

PGSuper assumes the barriers are the same concrete as the deck.

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `barrier.type` | Text | `TrafficBarrierEntry::GetName()` of `CRailingSystem::GetExteriorRailing()` (left or right railing system of `CBridgeDescription2`) | `import_traffic_barriers` (`TrafficBarriers.cpp`): a clue. The named entry is used if its section matches the barrier's geometry, and it names a new `TrafficBarrierEntry` when none matches. The entry is set with `CRailingSystem::SetExteriorRailing` |
| `barrier.height` | Length | highest point of `TrafficBarrierEntry::GetBarrierPoints` (above the top of the deck) | – |
| `barrier.fc` | Stress | `IMaterials::GetDeckFc28` | – |
| `barrier.strength_class` | Text | the same, formatted | – |
| `barrier.max_aggregate_size` | Length | `IMaterials::GetDeckMaxAggrSize` | – |

### Strands (`IfcTendon`) and their material

The context carries the `WBFL::Materials::PsStrand`, the strand type, and the strand index.

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `tendon.fy` | Stress | `PsStrand::GetYieldStrength` | – |
| `tendon.fpu` | Stress | `PsStrand::GetUltimateStrength` | `get_strand_material` (`Strands.cpp`): the strand grade (Grade 250 below 260 ksi, Grade 300 from 285 ksi, else Grade 270), before `tendon.grade` |
| `tendon.grade` | Text | "ASTM A416 Grade " + `PsStrand::GetGrade(grade, US units)` | `get_strand_material`: the grade when there's no `tendon.fpu` ("250" or "300" in the text, else Grade 270), and stress relieved when the text says so (else low relaxation) |
| `tendon.diameter` | Length | `PsStrand::GetNominalDiameter`; the import prefers it to the geometry, since models may draw strands at another size (Iowa draws 0.6" strands as 5/8" circles) | `get_strand_material`: the strand size, before the `IfcTendonType.NominalDiameter` and the drawn diameter; the nearest strand in the pool within 0.5 mm |
| `tendon.tension_force` | Force | the strand's jacking force, `IStrandGeometry::GetPjack` / strand count of its type (as `IfcTendon.TensionForce`); the import uses it when the model has neither `PreStress` nor `TensionForce` | `get_jacking_force` (`Strands.cpp`): a strand's jacking force when it has neither `PreStress` nor `TensionForce` |
| `tendon.corrosion_treatment` | Text | "NONE" or "EPOXYCOATED" from `PsStrand::GetCoating` | `get_strand_material`: "EPOXYCOATED" picks a smooth epoxy coated strand (IFC doesn't say smooth or grit impregnated) |
| `tendon.status` | Text | "TEMPORARY" for `pgsTypes::Temporary` strands, else "NEW" (`Pset_ElementComponentCommon.Status` on the strand and `usBrPset_Common.Status` on its type, so temporary strands get their own type) | `import_strands` (`Strands.cpp`): "TEMPORARY" makes a temporary strand |
| `tendon.debond_start`, `tendon.debond_end` | Length | `IStrandGeometry::IsStrandDebonded` start and end lengths (with usBridge classification only) | `import_strands`: the debonded lengths at the start and end |
| `tendon.debond_midspan_to_start`, `tendon.debond_midspan_to_end` | Length | temporary strands only: half the plan length less 10 ft (assumes the WSDOT detail, bonded for 10 ft at each end) | – |

### Reinforcing bars (`IfcReinforcingBar`) and their material

The context carries the `WBFL::Materials::Rebar` and a `BarDefinition` (`BarShape.h`, made from PGSuper's longitudinal and transverse reinforcement by `GetBarType` in `Rebar.h`).

The girder bar imports (`GirderBars.h`, for `LongitudinalRebar.h` and `Stirrups.h`) read `rebar.bar_size` (the number in the text, e.g. "#5", "5"), and `rebar.grade` and `rebar.fy` for the girder's bar materials (`CLongitudinalRebarData::BarType`, `BarGrade`; `CShearData2::ShearBarType`, `ShearBarGrade`): ASTM A615 unless the grade text says A706 or A1035, and the grade from `rebar.fy`, else the number in the grade text.

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `rebar.fy`, `rebar.fu` | Stress | `Rebar::GetYieldStrength`, `GetUltimateStrength` | `rebar.fy` only, `get_bar_material` (`GirderBars.cpp`): the bar grade, nearest of 40, 60, 75, 80, 100, and 120 ksi |
| `rebar.elongation` | Ratio | `Rebar::GetElongation` | – |
| `rebar.grade` | Text | `WBFL::LRFD::RebarPool::GetMaterialName(type, grade)` | `get_bar_material`: the bar type (A615 unless the text says A706 or A1035), and the grade from the number in the text when there's no `rebar.fy` |
| `rebar.specification`, `rebar.specification_edition` | Text | `GetRebarSpecification`, `GetRebarSpecificationEdition` (`SteelSpecifications.h`) | – |
| `rebar.bar_mark` | Text | `BarDefinition::mark` | – |
| `rebar.bar_size` | Text | `RebarPool::GetBarSize(rebar size)`, e.g. "#5" | `get_bar_size` (`GirderBars.cpp`): the bar size from the number in the text, before the nominal diameter of the bar type or the bar and the drawn diameter |
| `rebar.bar_mass` | Mass | `BarDefinition::GetMass()` (centerline length × nominal area × density) | – |
| `rebar.bar_element`, `rebar.bar_use`, `rebar.bar_position` | Text | `BarDefinition::element`, `use`, `position` (no value when empty) | – |
| `rebar.cover_top`, `cover_side`, `cover_bottom`, `cover_end` | Length | `BarDefinition::cover_*` (absent when the bar has none) | – |

### Bearings

The context carries the `CBearingData2` of the girder at the pier face (`CPierData2::GetBearingData`) and its `BearingElevationDetails` (`IBridge::GetBearingElevationDetails`).

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `bearing.fixed_x` | Boolean | `CBearingData2::FixedX` (fixed along the girder) | `get_bearing_fixity` → `set_bearing_data` (`Piers.cpp`): sets `CBearingData2::FixedX`, applied with `CBridgeDescription2::SetBearingData` or `CPierData2::SetBearingData`. When only one of `fixed_x`/`fixed_y` is found (e.g. a single "Fixity" property), it's used for both |
| `bearing.fixed_y` | Boolean | `CBearingData2::FixedY` (fixed across the girder) | the same, `CBearingData2::FixedY` |
| `bearing.displacement_accommodated` | Boolean list | `{ !FixedX, !FixedY, false }` (along, across, normal to the girder) | – (the property is `"import": false`; the standard table reads its items 0 and 1 into `fixed_x`/`fixed_y` with `list_index` locations whose `map` inverts the value: displacement accommodated means not fixed) |
| `bearing.rotation_accommodated` | Boolean list | `{ true, true, true }` (elastomeric bearing) | – |
| `bearing.height` | Length | `CBearingData2::Height` | – (bearing dimensions are measured from the bearing geometry, `set_bearing_data`) |
| `bearing.length`, `bearing.width` | Length | `CBearingData2::Length`, `Width` (absent for a round bearing) | – |
| `bearing.diameter` | Length | `CBearingData2::Length` of a round bearing (`Shape == bsRound`) | – |
| `bearing.plate_count` | Count | `CBearingData2::NumIntLayers + 1` | – |
| `bearing.plate_thickness` | Length list | `ShimThickness`, one per plate | – |
| `bearing.cover_thickness` | Length | `CBearingData2::CoverThickness` | – |
| `bearing.layer_thickness` | Length | `CBearingData2::ElastomerThickness` | – |
| `bearing.top_elevation` | Length | `BearingElevationDetails::TopBrgElevation` | – |
| `bearing.bottom_elevation` | Length | `BearingElevationDetails::BrgSeatElevation` | – |
| `bearing.girder_slope` | Ratio | `BearingElevationDetails::GirderGrade` | – |
| `bearing.support_condition` | Text | the pier's boundary condition (`CPierData2::GetBoundaryConditionType`), not the elastomeric bearing's fixity: "Fixed" for a hinge, "Expansion" for a roller, "Fixed" for continuous and integral connections (PGSuper's analysis pins the girders there) | `set_boundary_conditions` (`Piers.cpp`): hinge or roller for the pier when no diaphragm makes it continuous ("Fixed", "Hinge", "Pinned" → hinge; "Expansion", "Roller", "Free", "Sliding" → roller), set with `CPierData2::SetBoundaryConditionType` |

### Pier diaphragms (`IfcBeam.DIAPHRAGM`)

The context carries the pier index. PGSuper uses the deck concrete for the diaphragms. The exporter writes the diaphragms from `IBridge::GetPierDiaphragmSize` (`CreatePierDiaphragms`); the import reads continuity from their geometry (a diaphragm that the girder ends of both spans are in makes the pier continuous after the deck), see Phase 2 in [IfcImportPlan.md](IfcImportPlan.md).

| Target | Kind | Export: PGSuper data | Import |
|---|---|---|---|
| `diaphragm.support_condition` | Text | the pier's boundary condition: "Fixed" (hinge, and an integral abutment), "Expansion" (roller), or "Continuous" (continuous and integral connections; integral ones are exported as continuous for now) | `set_boundary_conditions`, after the diaphragm geometry and the bearings' `bearing.support_condition`: hinge, roller, or continuous after the deck |
| `diaphragm.fc` | Stress | `IMaterials::GetDeckFc28` | – |
| `diaphragm.strength_class` | Text | the same, formatted | – |
| `diaphragm.max_aggregate_size` | Length | `IMaterials::GetDeckMaxAggrSize` | – |

## Data imported without targets

Most of the import comes from geometry and IFC structure, not from table targets, and so is not in this list: the alignment and profile, pier stations and skews, girder layout and spacing, end distances and bearing offsets, bearing dimensions, continuity from the pier diaphragms, deck edges, depth, haunches and slab offset, roadway section, and the barrier sections. These are described in [IfcImportPlan.md](IfcImportPlan.md).
