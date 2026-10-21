# PGSuper-Display: round trip

- Imported: `results/PGSuper-Display.pgs`
- Expected: `models/PGSuper_Import_Model.pgs`

Alignment, profile, and deck edges are compared at at[0] = 20+00.00, at[1] = 20+50.00, at[2] = 21+00.00, at[3] = 21+65.00, at[4] = 22+30.00, at[5] = 22+80.00, at[6] = 23+30.00

**358 of 525 match** - mismatch: 13, not imported: 0, missing: 0, not comparable: 0, default match: 154, extra: 0

A *default match* equals the expected value only because the template has the same value.

## mismatch (13)

| Value | Expected | Imported | Note |
|---|---|---|---|
| `group[0].girder[0].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[0].girder[1].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[0].girder[2].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[0].girder[3].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[1].girder[0].strands.harped.paths` | 9 strands | 9 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[1].girder[1].strands.harped.paths` | 9 strands | 9 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[1].girder[2].strands.harped.paths` | 9 strands | 9 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[1].girder[3].strands.harped.paths` | 9 strands | 9 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[1].girder[4].strands.harped.paths` | 9 strands | 9 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[2].girder[0].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[2].girder[1].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[2].girder[2].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |
| `group[2].girder[3].strands.harped.paths` | 8 strands | 8 strands | 8 expected strands without a match within 0.118 in, 8 imported strands left over, largest distance to the closest strand 1.000 in |

## default match (154)

| Value | Expected | Imported | Note |
|---|---|---|---|
| `bridge.girder_family` | I-Beam | I-Beam | template default |
| `bridge.girder_orientation` | 0 | 0 | template default |
| `bridge.slab_offset_type` | 0 | 0 | template default |
| `bridge.fillet` | 0.750 in | 0.750 in | template default |
| `pier[0].connection_type` | 1 | 1 | template default |
| `pier[0].back.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[0].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[0].back.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[0].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[0].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[0].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[0].ahead.bearing.count` | 1 | 1 | template default |
| `pier[0].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[0].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[0].ahead.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[0].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[0].ahead.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[0].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[1].connection_type` | 6 | 6 | template default |
| `pier[1].back.bearing.shape` | 0 | 0 | template default |
| `pier[1].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[1].back.bearing.count` | 1 | 1 | template default |
| `pier[1].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[1].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[1].back.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[1].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[1].back.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[1].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[1].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[1].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[1].ahead.bearing.count` | 1 | 1 | template default |
| `pier[1].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[1].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[1].ahead.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[1].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[1].ahead.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[1].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[2].connection_type` | 6 | 6 | template default |
| `pier[2].back.bearing.shape` | 0 | 0 | template default |
| `pier[2].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[2].back.bearing.count` | 1 | 1 | template default |
| `pier[2].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[2].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[2].back.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[2].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[2].back.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[2].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[2].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[2].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[2].ahead.bearing.count` | 1 | 1 | template default |
| `pier[2].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[2].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[2].ahead.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[2].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[2].ahead.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[2].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[3].connection_type` | 6 | 6 | template default |
| `pier[3].back.bearing.shape` | 0 | 0 | template default |
| `pier[3].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[3].back.bearing.count` | 1 | 1 | template default |
| `pier[3].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[3].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[3].back.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[3].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[3].back.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[3].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[3].ahead.end_distance` | 20.500 in | 20.500 in | template default |
| `pier[3].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[3].ahead.bearing_offset` | 32.500 in | 32.500 in | template default |
| `pier[3].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `group[0].piers` | 0-1 | 0-1 | template default |
| `group[1].girder_count` | 5 | 5 | template default |
| `deck.type` | 0 | 0 | template default |
| `deck.gross_depth` | 7.500 in | 7.500 in | template default |
| `deck.left_overhang_edge_depth` | 7.000 in | 7.000 in | template default |
| `deck.right_overhang_edge_depth` | 7.000 in | 7.000 in | template default |
| `deck.left_overhang_taper` | 1 | 1 | template default |
| `deck.right_overhang_taper` | 1 | 1 | template default |
| `deck.haunch_shape` | 1 | 1 | template default |
| `deck.fc` | 4.000 ksi | 4.000 ksi | template default |
| `deck.edge_point_count` | 1 | 1 | template default |
| `railing.left.exterior` | 42&dq; Single Slope | 42&dq; Single Slope | template default |
| `railing.right.exterior` | 42&dq; Single Slope | 42&dq; Single Slope | template default |
| `pier[0].skew` | 0.0000 deg | 0.0000 deg | template default |
| `pier[1].skew` | 0.0000 deg | 0.0000 deg | template default |
| `pier[2].skew` | 0.0000 deg | 0.0000 deg | template default |
| `pier[3].skew` | 0.0000 deg | 0.0000 deg | template default |
| `at[0].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[0].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[0].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[0].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `at[1].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[1].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[1].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[1].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `at[2].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[2].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[2].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[2].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `at[3].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[3].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[3].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[3].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `at[4].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[4].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[4].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[4].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `at[5].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[5].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[5].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[5].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `at[6].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[6].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[6].roadway.at_left_edge` | -0.300 ft | -0.300 ft | template default |
| `at[6].roadway.at_right_edge` | -0.300 ft | -0.300 ft | template default |
| `group[0].girder[0].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[0].girder[0].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[0].girder[0].strands.temporary.count` | 0 | 0 | template default |
| `group[0].girder[1].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[0].girder[1].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[0].girder[1].strands.temporary.count` | 0 | 0 | template default |
| `group[0].girder[2].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[0].girder[2].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[0].girder[2].strands.temporary.count` | 0 | 0 | template default |
| `group[0].girder[3].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[0].girder[3].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[0].girder[3].strands.temporary.count` | 0 | 0 | template default |
| `group[1].girder[0].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[1].girder[0].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[1].girder[0].strands.temporary.debonded.count` | 0 | 0 | template default |
| `group[1].girder[1].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[1].girder[1].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[1].girder[1].strands.temporary.debonded.count` | 0 | 0 | template default |
| `group[1].girder[2].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[1].girder[2].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[1].girder[2].strands.temporary.debonded.count` | 0 | 0 | template default |
| `group[1].girder[3].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[1].girder[3].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[1].girder[3].strands.temporary.debonded.count` | 0 | 0 | template default |
| `group[1].girder[4].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[1].girder[4].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[1].girder[4].strands.temporary.debonded.count` | 0 | 0 | template default |
| `group[2].girder[0].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[2].girder[0].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[2].girder[0].strands.temporary.count` | 0 | 0 | template default |
| `group[2].girder[1].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[2].girder[1].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[2].girder[1].strands.temporary.count` | 0 | 0 | template default |
| `group[2].girder[2].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[2].girder[2].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[2].girder[2].strands.temporary.count` | 0 | 0 | template default |
| `group[2].girder[3].strands.straight.debonded.count` | 0 | 0 | template default |
| `group[2].girder[3].strands.harped.debonded.count` | 0 | 0 | template default |
| `group[2].girder[3].strands.temporary.count` | 0 | 0 | template default |

## match (358)

| Value | Expected | Imported | Note |
|---|---|---|---|
| `bridge.same_girder_for_entire_bridge` | 0 | 0 |  |
| `bridge.girder_spacing_type` | 1 | 1 |  |
| `bridge.slab_offset` | 9.500 in | 9.500 in |  |
| `bridge.pier_count` | 4 | 4 |  |
| `pier[0].station` | 20+00.00 | 20+00.00 |  |
| `pier[0].ahead.bearing.width` | 47.000 in | 47.000 in |  |
| `pier[0].ahead.bearing.height` | 1.725 in | 1.725 in |  |
| `pier[0].ahead.spacing[0]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[0].ahead.spacing[1]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[0].ahead.spacing[2]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[1].station` | 21+00.00 | 21+00.00 |  |
| `pier[1].back.bearing.width` | 47.000 in | 47.000 in |  |
| `pier[1].back.bearing.height` | 1.724 in | 1.725 in |  |
| `pier[1].back.spacing[0]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[1].back.spacing[1]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[1].back.spacing[2]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[1].ahead.bearing.width` | 47.000 in | 47.000 in |  |
| `pier[1].ahead.bearing.height` | 1.724 in | 1.725 in |  |
| `pier[1].ahead.spacing[0]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[1].ahead.spacing[1]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[1].ahead.spacing[2]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[1].ahead.spacing[3]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[2].station` | 22+30.00 | 22+30.00 |  |
| `pier[2].back.bearing.width` | 47.000 in | 47.000 in |  |
| `pier[2].back.bearing.height` | 1.724 in | 1.725 in |  |
| `pier[2].back.spacing[0]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[2].back.spacing[1]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[2].back.spacing[2]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[2].back.spacing[3]` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) |  |
| `pier[2].ahead.bearing.width` | 47.000 in | 47.000 in |  |
| `pier[2].ahead.bearing.height` | 1.725 in | 1.725 in |  |
| `pier[2].ahead.spacing[0]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[2].ahead.spacing[1]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[2].ahead.spacing[2]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[3].station` | 23+30.00 | 23+30.00 |  |
| `pier[3].back.bearing.width` | 47.000 in | 47.000 in |  |
| `pier[3].back.bearing.height` | 1.724 in | 1.725 in |  |
| `pier[3].back.spacing[0]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[3].back.spacing[1]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `pier[3].back.spacing[2]` | 7.000 ft (type 0, location 0) | 7.000 ft (type 0, location 0) |  |
| `bridge.group_count` | 3 | 3 |  |
| `group[0].girder_count` | 4 | 4 |  |
| `group[0].girder[0].type` | WF66G | WF66G (IFC) |  |
| `group[0].girder[0].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[0].girder[0].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[0].girder[0].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[0].girder[0].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[0].girder[1].type` | WF66G | WF66G (IFC) |  |
| `group[0].girder[1].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[0].girder[1].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[0].girder[1].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[0].girder[1].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[0].girder[2].type` | WF66G | WF66G (IFC) |  |
| `group[0].girder[2].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[0].girder[2].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[0].girder[2].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[0].girder[2].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[0].girder[3].type` | WF66G | WF66G (IFC) |  |
| `group[0].girder[3].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[0].girder[3].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[0].girder[3].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[0].girder[3].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[1].piers` | 1-2 | 1-2 |  |
| `group[1].girder[0].type` | WF100G | WF100G (IFC) |  |
| `group[1].girder[0].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[1].girder[0].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[1].girder[0].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[1].girder[0].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[1].girder[1].type` | WF100G | WF100G (IFC) |  |
| `group[1].girder[1].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[1].girder[1].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[1].girder[1].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[1].girder[1].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[1].girder[2].type` | WF100G | WF100G (IFC) |  |
| `group[1].girder[2].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[1].girder[2].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[1].girder[2].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[1].girder[2].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[1].girder[3].type` | WF100G | WF100G (IFC) |  |
| `group[1].girder[3].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[1].girder[3].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[1].girder[3].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[1].girder[3].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[1].girder[4].type` | WF100G | WF100G (IFC) |  |
| `group[1].girder[4].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[1].girder[4].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[1].girder[4].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[1].girder[4].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[2].piers` | 2-3 | 2-3 |  |
| `group[2].girder_count` | 4 | 4 |  |
| `group[2].girder[0].type` | WF66G | WF66G (IFC) |  |
| `group[2].girder[0].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[2].girder[0].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[2].girder[0].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[2].girder[0].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[2].girder[1].type` | WF66G | WF66G (IFC) |  |
| `group[2].girder[1].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[2].girder[1].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[2].girder[1].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[2].girder[1].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[2].girder[2].type` | WF66G | WF66G (IFC) |  |
| `group[2].girder[2].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[2].girder[2].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[2].girder[2].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[2].girder[2].fci` | 4.000 ksi | 4.000 ksi |  |
| `group[2].girder[3].type` | WF66G | WF66G (IFC) |  |
| `group[2].girder[3].slab_offset_start` | 9.500 in | 9.500 in |  |
| `group[2].girder[3].slab_offset_end` | 9.500 in | 9.500 in |  |
| `group[2].girder[3].fc` | 5.000 ksi | 5.000 ksi |  |
| `group[2].girder[3].fci` | 4.000 ksi | 4.000 ksi |  |
| `at[0].alignment.point` | (50000.000, 50000.000) ft | (50000.000, 50000.000) ft |  |
| `at[0].alignment.direction` | 18.0031 deg | 18.0031 deg |  |
| `at[0].profile.elevation` | 104.830 ft | 104.830 ft |  |
| `at[0].profile.grade` | 0.6600 % | 0.6600 % |  |
| `at[1].alignment.point` | (50047.474, 50015.691) ft | (50047.474, 50015.691) ft |  |
| `at[1].alignment.direction` | 18.5761 deg | 18.5761 deg |  |
| `at[1].profile.elevation` | 105.117 ft | 105.117 ft |  |
| `at[1].profile.grade` | 0.4900 % | 0.4900 % |  |
| `at[2].alignment.point` | (50094.789, 50031.856) ft | (50094.789, 50031.856) ft |  |
| `at[2].alignment.direction` | 19.1490 deg | 19.1490 deg |  |
| `at[2].profile.elevation` | 105.320 ft | 105.320 ft |  |
| `at[2].profile.grade` | 0.3200 % | 0.3200 % |  |
| `at[3].alignment.point` | (50156.052, 50053.576) ft | (50156.052, 50053.576) ft |  |
| `at[3].alignment.direction` | 19.8939 deg | 19.8939 deg |  |
| `at[3].profile.elevation` | 105.456 ft | 105.456 ft |  |
| `at[3].profile.grade` | 0.0990 % | 0.0990 % |  |
| `at[4].alignment.point` | (50217.027, 50076.091) ft | (50217.027, 50076.091) ft |  |
| `at[4].alignment.direction` | 20.6387 deg | 20.6387 deg |  |
| `at[4].profile.elevation` | 105.449 ft | 105.449 ft |  |
| `at[4].profile.grade` | -0.1220 % | -0.1220 % |  |
| `at[5].alignment.point` | (50263.730, 50093.948) ft | (50263.730, 50093.948) ft |  |
| `at[5].alignment.direction` | 21.2117 deg | 21.2117 deg |  |
| `at[5].profile.elevation` | 105.345 ft | 105.345 ft |  |
| `at[5].profile.grade` | -0.2920 % | -0.2920 % |  |
| `at[6].alignment.point` | (50310.251, 50112.272) ft | (50310.251, 50112.272) ft |  |
| `at[6].alignment.direction` | 21.7846 deg | 21.7846 deg |  |
| `at[6].profile.elevation` | 105.157 ft | 105.157 ft |  |
| `at[6].profile.grade` | -0.4620 % | -0.4620 % |  |
| `group[0].girder[0].strands.count` | 28 | 28 |  |
| `group[0].girder[0].strands.straight.count` | 20 | 20 |  |
| `group[0].girder[0].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[0].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[0].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[0].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[0].girder[0].strands.harped.count` | 8 | 8 |  |
| `group[0].girder[0].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[0].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[0].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[0].strands.harp_points` | 39.114, 58.675 ft | 39.114, 58.675 ft |  |
| `group[0].girder[1].strands.count` | 28 | 28 |  |
| `group[0].girder[1].strands.straight.count` | 20 | 20 |  |
| `group[0].girder[1].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[1].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[1].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[1].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[0].girder[1].strands.harped.count` | 8 | 8 |  |
| `group[0].girder[1].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[1].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[1].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[1].strands.harp_points` | 39.173, 58.757 ft | 39.173, 58.757 ft |  |
| `group[0].girder[2].strands.count` | 28 | 28 |  |
| `group[0].girder[2].strands.straight.count` | 20 | 20 |  |
| `group[0].girder[2].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[2].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[2].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[2].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[0].girder[2].strands.harped.count` | 8 | 8 |  |
| `group[0].girder[2].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[2].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[2].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[2].strands.harp_points` | 39.229, 58.842 ft | 39.229, 58.842 ft |  |
| `group[0].girder[3].strands.count` | 28 | 28 |  |
| `group[0].girder[3].strands.straight.count` | 20 | 20 |  |
| `group[0].girder[3].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[3].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[3].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[3].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[0].girder[3].strands.harped.count` | 8 | 8 |  |
| `group[0].girder[3].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[0].girder[3].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[0].girder[3].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[0].girder[3].strands.harp_points` | 39.285, 58.927 ft | 39.285, 58.927 ft |  |
| `group[1].girder[0].strands.count` | 29 | 29 |  |
| `group[1].girder[0].strands.straight.count` | 18 | 18 |  |
| `group[1].girder[0].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[0].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[0].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[0].strands.straight.paths` | 18 strands | 18 strands |  |
| `group[1].girder[0].strands.harped.count` | 9 | 9 |  |
| `group[1].girder[0].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[0].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[0].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[0].strands.temporary.count` | 2 | 2 |  |
| `group[1].girder[0].strands.temporary.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[0].strands.temporary.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[0].strands.temporary.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[0].strands.temporary.paths` | 2 strands | 2 strands |  |
| `group[1].girder[0].strands.harp_points` | 51.073, 76.611 ft | 51.073, 76.611 ft |  |
| `group[1].girder[1].strands.count` | 29 | 29 |  |
| `group[1].girder[1].strands.straight.count` | 18 | 18 |  |
| `group[1].girder[1].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[1].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[1].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[1].strands.straight.paths` | 18 strands | 18 strands |  |
| `group[1].girder[1].strands.harped.count` | 9 | 9 |  |
| `group[1].girder[1].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[1].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[1].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[1].strands.temporary.count` | 2 | 2 |  |
| `group[1].girder[1].strands.temporary.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[1].strands.temporary.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[1].strands.temporary.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[1].strands.temporary.paths` | 2 strands | 2 strands |  |
| `group[1].girder[1].strands.harp_points` | 51.135, 76.703 ft | 51.135, 76.703 ft |  |
| `group[1].girder[2].strands.count` | 29 | 29 |  |
| `group[1].girder[2].strands.straight.count` | 18 | 18 |  |
| `group[1].girder[2].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[2].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[2].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[2].strands.straight.paths` | 18 strands | 18 strands |  |
| `group[1].girder[2].strands.harped.count` | 9 | 9 |  |
| `group[1].girder[2].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[2].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[2].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[2].strands.temporary.count` | 2 | 2 |  |
| `group[1].girder[2].strands.temporary.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[2].strands.temporary.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[2].strands.temporary.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[2].strands.temporary.paths` | 2 strands | 2 strands |  |
| `group[1].girder[2].strands.harp_points` | 51.198, 76.798 ft | 51.198, 76.798 ft |  |
| `group[1].girder[3].strands.count` | 29 | 29 |  |
| `group[1].girder[3].strands.straight.count` | 18 | 18 |  |
| `group[1].girder[3].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[3].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[3].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[3].strands.straight.paths` | 18 strands | 18 strands |  |
| `group[1].girder[3].strands.harped.count` | 9 | 9 |  |
| `group[1].girder[3].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[3].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[3].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[3].strands.temporary.count` | 2 | 2 |  |
| `group[1].girder[3].strands.temporary.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[3].strands.temporary.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[3].strands.temporary.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[3].strands.temporary.paths` | 2 strands | 2 strands |  |
| `group[1].girder[3].strands.harp_points` | 51.260, 76.890 ft | 51.260, 76.890 ft |  |
| `group[1].girder[4].strands.count` | 29 | 29 |  |
| `group[1].girder[4].strands.straight.count` | 18 | 18 |  |
| `group[1].girder[4].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[4].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[4].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[4].strands.straight.paths` | 18 strands | 18 strands |  |
| `group[1].girder[4].strands.harped.count` | 9 | 9 |  |
| `group[1].girder[4].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[4].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[4].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[4].strands.temporary.count` | 2 | 2 |  |
| `group[1].girder[4].strands.temporary.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[1].girder[4].strands.temporary.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[1].girder[4].strands.temporary.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[1].girder[4].strands.temporary.paths` | 2 strands | 2 strands |  |
| `group[1].girder[4].strands.harp_points` | 51.322, 76.985 ft | 51.322, 76.985 ft |  |
| `group[2].girder[0].strands.count` | 28 | 28 |  |
| `group[2].girder[0].strands.straight.count` | 20 | 20 |  |
| `group[2].girder[0].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[0].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[0].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[0].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[2].girder[0].strands.harped.count` | 8 | 8 |  |
| `group[2].girder[0].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[0].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[0].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[0].strands.harp_points` | 39.114, 58.675 ft | 39.114, 58.675 ft |  |
| `group[2].girder[1].strands.count` | 28 | 28 |  |
| `group[2].girder[1].strands.straight.count` | 20 | 20 |  |
| `group[2].girder[1].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[1].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[1].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[1].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[2].girder[1].strands.harped.count` | 8 | 8 |  |
| `group[2].girder[1].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[1].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[1].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[1].strands.harp_points` | 39.170, 58.757 ft | 39.170, 58.757 ft |  |
| `group[2].girder[2].strands.count` | 28 | 28 |  |
| `group[2].girder[2].strands.straight.count` | 20 | 20 |  |
| `group[2].girder[2].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[2].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[2].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[2].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[2].girder[2].strands.harped.count` | 8 | 8 |  |
| `group[2].girder[2].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[2].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[2].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[2].strands.harp_points` | 39.229, 58.842 ft | 39.229, 58.842 ft |  |
| `group[2].girder[3].strands.count` | 28 | 28 |  |
| `group[2].girder[3].strands.straight.count` | 20 | 20 |  |
| `group[2].girder[3].strands.straight.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[3].strands.straight.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[3].strands.straight.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[3].strands.straight.paths` | 20 strands | 20 strands |  |
| `group[2].girder[3].strands.harped.count` | 8 | 8 |  |
| `group[2].girder[3].strands.harped.material` | Grade 270 Low Relaxation 0.6" | Grade 270 Low Relaxation 0.6" |  |
| `group[2].girder[3].strands.harped.prestress` | 202.500 ksi | 202.500 ksi |  |
| `group[2].girder[3].strands.harped.tension_force` | 43.943 kip | 43.943 kip |  |
| `group[2].girder[3].strands.harp_points` | 39.285, 58.924 ft | 39.285, 58.924 ft |  |
| `group[0].girder[0].longitudinal_bars.count` | 6 | 6 |  |
| `group[0].girder[0].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[0].girder[0].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[0].girder[0].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[0].girder[1].longitudinal_bars.count` | 6 | 6 |  |
| `group[0].girder[1].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[0].girder[1].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[0].girder[1].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[0].girder[2].longitudinal_bars.count` | 6 | 6 |  |
| `group[0].girder[2].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[0].girder[2].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[0].girder[2].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[0].girder[3].longitudinal_bars.count` | 6 | 6 |  |
| `group[0].girder[3].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[0].girder[3].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[0].girder[3].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[1].girder[0].longitudinal_bars.count` | 6 | 6 |  |
| `group[1].girder[0].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[1].girder[0].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[1].girder[0].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[1].girder[1].longitudinal_bars.count` | 6 | 6 |  |
| `group[1].girder[1].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[1].girder[1].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[1].girder[1].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[1].girder[2].longitudinal_bars.count` | 6 | 6 |  |
| `group[1].girder[2].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[1].girder[2].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[1].girder[2].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[1].girder[3].longitudinal_bars.count` | 6 | 6 |  |
| `group[1].girder[3].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[1].girder[3].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[1].girder[3].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[1].girder[4].longitudinal_bars.count` | 6 | 6 |  |
| `group[1].girder[4].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[1].girder[4].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[1].girder[4].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[2].girder[0].longitudinal_bars.count` | 6 | 6 |  |
| `group[2].girder[0].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[2].girder[0].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[2].girder[0].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[2].girder[1].longitudinal_bars.count` | 6 | 6 |  |
| `group[2].girder[1].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[2].girder[1].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[2].girder[1].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[2].girder[2].longitudinal_bars.count` | 6 | 6 |  |
| `group[2].girder[2].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[2].girder[2].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[2].girder[2].longitudinal_bars.paths` | 6 bars | 6 bars |  |
| `group[2].girder[3].longitudinal_bars.count` | 6 | 6 |  |
| `group[2].girder[3].longitudinal_bars.sizes` | 6 #5 | 6 #5 |  |
| `group[2].girder[3].longitudinal_bars.material` | AASHTO M31 (A615) - Grade 60 | AASHTO M31 (A615) - Grade 60 |  |
| `group[2].girder[3].longitudinal_bars.paths` | 6 bars | 6 bars |  |
