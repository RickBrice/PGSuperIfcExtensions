Element Roles and Targets {#reference}
==========================

This list is made from the IFC extension's source code when the documentation is built, so it matches this version of the extension.

Element Roles
-------------

An element role is a kind of IFC element that PGSuper reads or writes. A table's selectors say which IFC elements play each role (see @ref element_roles).

| Element role |
|---|
| `project` |
| `site` |
| `bridge` |
| `superstructure` |
| `substructure` |
| `deck_part` |
| `pier` |
| `abutment` |
| `foundation` |
| `alignment` |
| `referent` |
| `girder` |
| `girder_assembly` |
| `closure_joint` |
| `tendon` |
| `tendon_bundle` |
| `rebar` |
| `reinforcement_cage` |
| `deck` |
| `haunch` |
| `bearing` |
| `barrier` |

Targets
-------

A target is a PGSuper value with a fixed name. A table can't add targets; it says where each target is found in IFC (see @ref targets) and which exported properties hold it (see @ref property_sets).

| Target | Element role | Kind | Description |
|---|---|---|---|
| `bridge.number_of_spans` | `bridge` | count | number of spans |
| `bridge.length` | `bridge` | length | bridge length |
| `bridge.roadway_width` | `bridge` | length | roadway width at the start of the bridge |
| `bridge.start_station` | `bridge` | length | station of the start of the bridge |
| `bridge.end_station` | `bridge` | length | station of the end of the bridge |
| `bridge.max_skew` | `bridge` | angle | largest pier skew |
| `bridge.design_live_load` | `bridge` | text | design vehicular live load |
| `pier.station` | `pier` | length | station of the CL pier |
| `pier.skew` | `pier` | angle | pier skew |
| `girder.designation` | `girder` | text | girder designation (e.g. \ |
| `girder.type_names` | `girder` | list of text | girder type name |
| `girder.family_name` | `girder` | text | girder family (e.g. I-Beam) |
| `girder.fc` | `girder` | stress | girder concrete strength, f'c |
| `girder.fci` | `girder` | stress | girder concrete strength at release, f'ci |
| `girder.strength_class` | `girder` | text | girder concrete strength class |
| `girder.max_aggregate_size` | `girder` | length | girder concrete maximum aggregate size |
| `girder.assembly_place` | `girder` | text | girder assembly place (e.g. FACTORY) |
| `girder.casting_method` | `girder` | text | girder casting method (e.g. PRECAST) |
| `girder.span` | `girder` | length | girder span length |
| `girder.slope` | `girder` | angle | girder slope angle |
| `girder.roll` | `girder` | angle | girder roll angle |
| `girder.batter` | `girder` | angle | girder end batter |
| `girder.bunk_point` | `girder` | length | girder hauling support location |
| `girder.jacking_stress` | `girder` | stress | jacking stress of the permanent strands |
| `girder.camber_ratio` | `girder` | ratio | camber at release divided by the girder length |
| `girder.camber_at_release` | `girder` | length | camber at release |
| `girder.camber_after_losses` | `girder` | length | excess camber |
| `girder.screed_camber` | `girder` | length | screed camber |
| `girder.length` | `girder` | length | girder plan length |
| `girder.cross_section_area` | `girder` | Area | girder cross section area at mid-span |
| `girder.outer_surface_area` | `girder` | Area | girder outer surface area (length times perimeter) |
| `girder.gross_weight` | `girder` | Mass | girder mass (weight and mass are the same in IfcQuantityWeight) |
| `deck.gross_depth` | `deck` | length | deck gross depth |
| `deck.fc` | `deck` | stress | deck concrete strength, f'c |
| `deck.strength_class` | `deck` | text | deck concrete strength class |
| `deck.length` | `deck` | length | deck length (bridge length) |
| `deck.depth` | `deck` | length | deck gross depth (input) |
| `deck.perimeter` | `deck` | length | perimeter of the deck plan outline |
| `deck.gross_area` | `deck` | Area | area of the deck plan outline |
| `deck.max_aggregate_size` | `deck` | length | deck concrete maximum aggregate size |
| `barrier.fc` | `barrier` | stress | barrier concrete strength, f'c |
| `barrier.strength_class` | `barrier` | text | barrier concrete strength class |
| `barrier.type` | `barrier` | text | exterior barrier (traffic barrier library entry) |
| `barrier.height` | `barrier` | length | exterior barrier height above the top of the deck (from its shape) |
| `barrier.max_aggregate_size` | `barrier` | length | barrier concrete maximum aggregate size |
| `tendon.fy` | `tendon` | stress | strand yield strength |
| `tendon.fpu` | `tendon` | stress | strand ultimate strength |
| `tendon.grade` | `tendon` | text | strand grade (e.g. ASTM A416 Grade 270) |
| `tendon.corrosion_treatment` | `tendon` | text | strand corrosion treatment (NONE or EPOXYCOATED) |
| `tendon.status` | `tendon` | text | strand status (TEMPORARY for temporary strands) |
| `tendon.debond_start` | `tendon` | length | debond length at the start of a debonded strand |
| `tendon.debond_end` | `tendon` | length | debond length at the end of a debonded strand |
| `tendon.debond_midspan_to_start` | `tendon` | length | debonded length from mid-span toward the start of a temporary strand |
| `tendon.debond_midspan_to_end` | `tendon` | length | debonded length from mid-span toward the end of a temporary strand |
| `rebar.fy` | `rebar` | stress | reinforcing bar yield strength |
| `rebar.fu` | `rebar` | stress | reinforcing bar ultimate strength |
| `rebar.elongation` | `rebar` | ratio | reinforcing bar elongation at fracture |
| `rebar.grade` | `rebar` | text | reinforcing bar material (e.g. AASHTO M31 (A615) - Grade 60) |
| `rebar.specification` | `rebar` | text | reinforcing bar specification (e.g. ASTM A615 (AASHTO M31)) |
| `rebar.specification_edition` | `rebar` | text | edition of the reinforcing bar specification |
| `bearing.fixed_x` | `bearing` | true/false | bearing fixed along the girder |
| `bearing.fixed_y` | `bearing` | true/false | bearing fixed across the girder |
