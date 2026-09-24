# PGSuper-Skew-Bearings-System: round trip

- Imported: `results/PGSuper-Skew-Bearings-System.pgs`
- Expected: `models/PGSuper_Skew_Straight_Bearings.pgs`

Alignment, profile, and deck edges are compared at at[0] = 0+00.00, at[1] = 0+55.00, at[2] = 1+10.00, at[3] = 1+65.00, at[4] = 2+20.00, at[5] = 2+75.00, at[6] = 3+30.00, at[7] = 3+85.00, at[8] = 4+40.00

**58 of 294 match** - mismatch: 0, not imported: 24, missing: 0, not comparable: 0, default match: 212, extra: 0

A *default match* equals the expected value only because the template has the same value.

## not imported (24)

| Value | Expected | Imported | Note |
|---|---|---|---|
| `pier[1].connection_type` | 2 | 6 | template default |
| `pier[2].connection_type` | 2 | 6 | template default |
| `pier[3].connection_type` | 2 | 6 | template default |
| `pier[4].connection_type` | 1 | 6 | template default |
| `group[0].girder[0].type` | WF42G | Unknown_Girder_Type | template default |
| `group[0].girder[1].type` | WF42G | Unknown_Girder_Type | template default |
| `group[0].girder[2].type` | WF42G | Unknown_Girder_Type | template default |
| `group[0].girder[3].type` | WF42G | Unknown_Girder_Type | template default |
| `group[0].girder[4].type` | WF42G | Unknown_Girder_Type | template default |
| `group[1].girder[0].type` | WF42G | Unknown_Girder_Type | template default |
| `group[1].girder[1].type` | WF42G | Unknown_Girder_Type | template default |
| `group[1].girder[2].type` | WF42G | Unknown_Girder_Type | template default |
| `group[1].girder[3].type` | WF42G | Unknown_Girder_Type | template default |
| `group[1].girder[4].type` | WF42G | Unknown_Girder_Type | template default |
| `group[2].girder[0].type` | WF42G | Unknown_Girder_Type | template default |
| `group[2].girder[1].type` | WF42G | Unknown_Girder_Type | template default |
| `group[2].girder[2].type` | WF42G | Unknown_Girder_Type | template default |
| `group[2].girder[3].type` | WF42G | Unknown_Girder_Type | template default |
| `group[2].girder[4].type` | WF42G | Unknown_Girder_Type | template default |
| `group[3].girder[0].type` | WF42G | Unknown_Girder_Type | template default |
| `group[3].girder[1].type` | WF42G | Unknown_Girder_Type | template default |
| `group[3].girder[2].type` | WF42G | Unknown_Girder_Type | template default |
| `group[3].girder[3].type` | WF42G | Unknown_Girder_Type | template default |
| `group[3].girder[4].type` | WF42G | Unknown_Girder_Type | template default |

## default match (212)

| Value | Expected | Imported | Note |
|---|---|---|---|
| `bridge.girder_family` | I-Beam | I-Beam | template default |
| `bridge.girder_orientation` | 0 | 0 | template default |
| `bridge.same_girder_for_entire_bridge` | -1 | -1 | template default |
| `bridge.girder_spacing_type` | 0 | 0 | template default |
| `bridge.girder_spacing` | 6.000 ft (type 0, location 0) | 6.000 ft (type 0, location 0) | template default |
| `bridge.slab_offset_type` | 0 | 0 | template default |
| `bridge.slab_offset` | 11.000 in | 11.000 in | template default |
| `bridge.fillet` | 0.750 in | 0.719 in | template default |
| `pier[0].station` | 0+00.00 | 0+00.00 | template default |
| `pier[0].connection_type` | 1 | 1 | template default |
| `pier[0].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[0].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[0].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[0].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[0].ahead.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[0].ahead.bearing.count` | 1 | 1 | template default |
| `pier[0].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[0].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[0].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[0].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[1].back.bearing.shape` | 0 | 0 | template default |
| `pier[1].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[1].back.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[1].back.bearing.count` | 1 | 1 | template default |
| `pier[1].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[1].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[1].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[1].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[1].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[1].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[1].ahead.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[1].ahead.bearing.count` | 1 | 1 | template default |
| `pier[1].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[1].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[1].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[1].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[2].back.bearing.shape` | 0 | 0 | template default |
| `pier[2].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[2].back.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[2].back.bearing.count` | 1 | 1 | template default |
| `pier[2].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[2].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[2].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[2].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[2].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[2].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[2].ahead.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[2].ahead.bearing.count` | 1 | 1 | template default |
| `pier[2].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[2].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[2].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[2].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[3].back.bearing.shape` | 0 | 0 | template default |
| `pier[3].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[3].back.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[3].back.bearing.count` | 1 | 1 | template default |
| `pier[3].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[3].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[3].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[3].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[3].ahead.bearing.shape` | 0 | 0 | template default |
| `pier[3].ahead.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[3].ahead.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[3].ahead.bearing.count` | 1 | 1 | template default |
| `pier[3].ahead.bearing.fixed_x` | 0 | 0 | template default |
| `pier[3].ahead.bearing.fixed_y` | 0 | 0 | template default |
| `pier[3].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[3].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[4].back.bearing.shape` | 0 | 0 | template default |
| `pier[4].back.bearing.length` | 12.000 in | 12.000 in | template default |
| `pier[4].back.bearing.width` | 12.000 in | 12.000 in | template default |
| `pier[4].back.bearing.count` | 1 | 1 | template default |
| `pier[4].back.bearing.fixed_x` | 0 | 0 | template default |
| `pier[4].back.bearing.fixed_y` | 0 | 0 | template default |
| `pier[4].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[4].back.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `pier[4].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier | template default |
| `pier[4].ahead.bearing_offset_measure` | NormalToPier | NormalToPier | template default |
| `group[0].piers` | 0-1 | 0-1 | template default |
| `group[0].girder_count` | 5 | 5 | template default |
| `group[0].girder[0].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[0].girder[0].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[0].girder[0].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[0].girder[0].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[0].girder[1].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[0].girder[1].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[0].girder[1].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[0].girder[1].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[0].girder[2].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[0].girder[2].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[0].girder[2].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[0].girder[2].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[0].girder[3].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[0].girder[3].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[0].girder[3].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[0].girder[3].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[0].girder[4].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[0].girder[4].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[0].girder[4].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[0].girder[4].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[1].girder_count` | 5 | 5 | template default |
| `group[1].girder[0].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[1].girder[0].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[1].girder[0].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[1].girder[0].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[1].girder[1].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[1].girder[1].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[1].girder[1].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[1].girder[1].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[1].girder[2].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[1].girder[2].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[1].girder[2].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[1].girder[2].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[1].girder[3].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[1].girder[3].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[1].girder[3].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[1].girder[3].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[1].girder[4].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[1].girder[4].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[1].girder[4].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[1].girder[4].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[2].girder_count` | 5 | 5 | template default |
| `group[2].girder[0].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[2].girder[0].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[2].girder[0].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[2].girder[0].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[2].girder[1].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[2].girder[1].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[2].girder[1].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[2].girder[1].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[2].girder[2].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[2].girder[2].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[2].girder[2].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[2].girder[2].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[2].girder[3].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[2].girder[3].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[2].girder[3].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[2].girder[3].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[2].girder[4].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[2].girder[4].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[2].girder[4].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[2].girder[4].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[3].girder_count` | 5 | 5 | template default |
| `group[3].girder[0].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[3].girder[0].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[3].girder[0].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[3].girder[0].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[3].girder[1].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[3].girder[1].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[3].girder[1].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[3].girder[1].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[3].girder[2].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[3].girder[2].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[3].girder[2].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[3].girder[2].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[3].girder[3].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[3].girder[3].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[3].girder[3].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[3].girder[3].fci` | 5.000 ksi | 5.000 ksi | template default |
| `group[3].girder[4].slab_offset_start` | 11.000 in | 11.000 in | template default |
| `group[3].girder[4].slab_offset_end` | 11.000 in | 11.000 in | template default |
| `group[3].girder[4].fc` | 6.000 ksi | 6.000 ksi | template default |
| `group[3].girder[4].fci` | 5.000 ksi | 5.000 ksi | template default |
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
| `pier[1].skew` | 0.0000 deg | 0.0000 deg | template default |
| `at[0].alignment.point` | (0.000, 0.000) ft | (0.000, 0.000) ft | template default |
| `at[0].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[0].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[0].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[0].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[1].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[1].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[1].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[1].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[2].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[2].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[2].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[2].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[3].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[3].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[3].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[3].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[4].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[4].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[4].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[4].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[5].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[5].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[5].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[5].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[6].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[6].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[6].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[6].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[7].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[7].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[7].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[7].deck.right_edge` | 15.000 ft | 15.000 ft | template default |
| `at[8].profile.elevation` | 0.000 ft | 0.000 ft | template default |
| `at[8].profile.grade` | 0.0000 % | 0.0000 % | template default |
| `at[8].deck.left_edge` | 15.000 ft | 15.000 ft | template default |
| `at[8].deck.right_edge` | 15.000 ft | 15.000 ft | template default |

## match (58)

| Value | Expected | Imported | Note |
|---|---|---|---|
| `bridge.pier_count` | 5 | 5 |  |
| `pier[0].back.end_distance` | 15.000 in | 15.000 in |  |
| `pier[0].back.bearing_offset` | 20.000 in | 20.000 in |  |
| `pier[0].ahead.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[0].ahead.end_distance` | 15.000 in | 15.000 in |  |
| `pier[0].ahead.bearing_offset` | 20.000 in | 20.000 in |  |
| `pier[1].station` | 1+10.00 | 1+10.00 |  |
| `pier[1].back.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[1].back.end_distance` | 15.000 in | 15.000 in |  |
| `pier[1].back.bearing_offset` | 15.000 in | 15.000 in |  |
| `pier[1].ahead.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[1].ahead.end_distance` | 18.000 in | 18.000 in |  |
| `pier[1].ahead.bearing_offset` | 30.000 in | 30.000 in |  |
| `pier[2].station` | 2+20.00 | 2+20.00 |  |
| `pier[2].back.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[2].back.end_distance` | 18.000 in | 18.000 in |  |
| `pier[2].back.bearing_offset` | 30.000 in | 30.000 in |  |
| `pier[2].ahead.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[2].ahead.end_distance` | 10.000 in | 10.000 in |  |
| `pier[2].ahead.bearing_offset` | 15.000 in | 15.000 in |  |
| `pier[3].station` | 3+30.00 | 3+30.00 |  |
| `pier[3].back.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[3].back.end_distance` | 6.000 in | 6.000 in |  |
| `pier[3].back.bearing_offset` | 6.000 in | 6.000 in |  |
| `pier[3].ahead.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[3].ahead.end_distance` | 6.000 in | 6.000 in |  |
| `pier[3].ahead.bearing_offset` | 6.000 in | 6.000 in |  |
| `pier[4].station` | 4+40.00 | 4+40.00 |  |
| `pier[4].back.bearing.height` | 2.000 in | 2.000 in |  |
| `pier[4].back.end_distance` | 25.000 in | 25.000 in |  |
| `pier[4].back.bearing_offset` | 18.000 in | 18.000 in |  |
| `pier[4].ahead.end_distance` | 25.000 in | 25.000 in |  |
| `pier[4].ahead.bearing_offset` | 18.000 in | 18.000 in |  |
| `bridge.group_count` | 4 | 4 |  |
| `group[1].piers` | 1-2 | 1-2 |  |
| `group[2].piers` | 2-3 | 2-3 |  |
| `group[3].piers` | 3-4 | 3-4 |  |
| `pier[0].skew` | 20.0000 deg | 20.0000 deg |  |
| `pier[2].skew` | -15.0000 deg | -15.0000 deg |  |
| `pier[3].skew` | -30.0000 deg | -30.0000 deg |  |
| `pier[4].skew` | 25.0000 deg | 25.0000 deg |  |
| `at[0].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[1].alignment.point` | (-27.500, 47.631) ft | (-27.500, 47.631) ft |  |
| `at[1].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[2].alignment.point` | (-55.000, 95.263) ft | (-55.000, 95.263) ft |  |
| `at[2].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[3].alignment.point` | (-82.500, 142.894) ft | (-82.500, 142.894) ft |  |
| `at[3].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[4].alignment.point` | (-110.000, 190.526) ft | (-110.000, 190.526) ft |  |
| `at[4].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[5].alignment.point` | (-137.500, 238.157) ft | (-137.500, 238.157) ft |  |
| `at[5].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[6].alignment.point` | (-165.000, 285.788) ft | (-165.000, 285.788) ft |  |
| `at[6].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[7].alignment.point` | (-192.500, 333.420) ft | (-192.500, 333.420) ft |  |
| `at[7].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
| `at[8].alignment.point` | (-220.000, 381.051) ft | (-220.000, 381.051) ft |  |
| `at[8].alignment.direction` | 120.0000 deg | 120.0000 deg |  |
