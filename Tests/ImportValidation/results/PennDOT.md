# PennDOT: IFC import

- Imported: `results/PennDOT.pgs`
- No expected values - this shows which values the import set

Alignment, profile, and deck edges are evaluated at at[0] = -0+06.00, at[1] = 0+24.00, at[2] = 0+54.00

**49 of 102 values set by the import**, 53 are template defaults

## template default (53)

| Value | Template | Imported |
|---|---|---|
| `bridge.girder_family` | I-Beam | I-Beam |
| `bridge.girder_orientation` | 0 | 0 |
| `bridge.same_girder_for_entire_bridge` | -1 | -1 |
| `bridge.girder_spacing_type` | 0 | 0 |
| `bridge.fillet` | 0.750 in | 0.750 in |
| `bridge.pier_count` | 2 | 2 |
| `pier[0].connection_type` | 1 | 1 |
| `pier[0].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier |
| `pier[0].back.bearing_offset_measure` | NormalToPier | NormalToPier |
| `pier[0].ahead.bearing.shape` | 0 | 0 |
| `pier[0].ahead.bearing.length` | 12.000 in | 12.000 in |
| `pier[0].ahead.bearing.count` | 1 | 1 |
| `pier[0].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier |
| `pier[0].ahead.bearing_offset_measure` | NormalToPier | NormalToPier |
| `pier[1].connection_type` | 1 | 6 |
| `pier[1].back.bearing.shape` | 0 | 0 |
| `pier[1].back.bearing.length` | 12.000 in | 12.000 in |
| `pier[1].back.bearing.count` | 1 | 1 |
| `pier[1].back.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier |
| `pier[1].back.bearing_offset_measure` | NormalToPier | NormalToPier |
| `pier[1].ahead.end_distance_measure` | FromBearingNormalToPier | FromBearingNormalToPier |
| `pier[1].ahead.bearing_offset_measure` | NormalToPier | NormalToPier |
| `bridge.group_count` | 1 | 1 |
| `group[0].piers` | 0-1 | 0-1 |
| `group[0].girder_count` | 5 | 5 |
| `group[0].girder[0].type` | Unknown_Girder_Type | Unknown_Girder_Type |
| `group[0].girder[0].fc` | 6.000 ksi | 6.000 ksi |
| `group[0].girder[0].fci` | 5.000 ksi | 5.000 ksi |
| `group[0].girder[1].type` | Unknown_Girder_Type | Unknown_Girder_Type |
| `group[0].girder[1].fc` | 6.000 ksi | 6.000 ksi |
| `group[0].girder[1].fci` | 5.000 ksi | 5.000 ksi |
| `group[0].girder[2].type` | Unknown_Girder_Type | Unknown_Girder_Type |
| `group[0].girder[2].fc` | 6.000 ksi | 6.000 ksi |
| `group[0].girder[2].fci` | 5.000 ksi | 5.000 ksi |
| `group[0].girder[3].type` | Unknown_Girder_Type | Unknown_Girder_Type |
| `group[0].girder[3].fc` | 6.000 ksi | 6.000 ksi |
| `group[0].girder[3].fci` | 5.000 ksi | 5.000 ksi |
| `group[0].girder[4].type` | Unknown_Girder_Type | Unknown_Girder_Type |
| `group[0].girder[4].fc` | 6.000 ksi | 6.000 ksi |
| `group[0].girder[4].fci` | 5.000 ksi | 5.000 ksi |
| `deck.type` | 0 | 0 |
| `deck.fc` | 4.000 ksi | 4.000 ksi |
| `deck.edge_point_count` | 1 | 1 |
| `railing.left.exterior` | 42&dq; Single Slope | 42&dq; Single Slope |
| `railing.right.exterior` | 42&dq; Single Slope | 42&dq; Single Slope |
| `pier[0].skew` | 0.0000 deg | 0.0000 deg |
| `pier[1].skew` | 0.0000 deg | 0.0000 deg |
| `at[0].profile.elevation` | 0.000 ft | 0.000 ft |
| `at[0].profile.grade` | 0.0000 % | 0.0000 % |
| `at[1].profile.elevation` | 0.000 ft | 0.000 ft |
| `at[1].profile.grade` | 0.0000 % | 0.0000 % |
| `at[2].profile.elevation` | 0.000 ft | 0.000 ft |
| `at[2].profile.grade` | 0.0000 % | 0.0000 % |

## set by the import (49)

| Value | Template | Imported |
|---|---|---|
| `bridge.girder_spacing` | 6.000 ft (type 0, location 0) | 6.500 ft (type 0, location 0) |
| `bridge.slab_offset_type` | 0 | 2 |
| `bridge.slab_offset` | 11.000 in |  |
| `pier[0].station` | 0+00.00 | -0+06.00 |
| `pier[0].back.end_distance` | 20.500 in | 12.000 in |
| `pier[0].back.bearing_offset` | 32.500 in | -0.000 in |
| `pier[0].ahead.bearing.width` | 12.000 in | 31.500 in |
| `pier[0].ahead.bearing.height` | 0.000 in | 0.750 in |
| `pier[0].ahead.bearing.fixed_x` | 0 | -1 |
| `pier[0].ahead.bearing.fixed_y` | 0 | -1 |
| `pier[0].ahead.end_distance` | 20.500 in | 12.000 in |
| `pier[0].ahead.bearing_offset` | 32.500 in | -0.000 in |
| `pier[1].station` | 0+00.00 | 0+54.00 |
| `pier[1].back.bearing.width` | 12.000 in | 31.500 in |
| `pier[1].back.bearing.height` | 0.000 in | 0.750 in |
| `pier[1].back.bearing.fixed_x` | 0 | -1 |
| `pier[1].back.bearing.fixed_y` | 0 | -1 |
| `pier[1].back.end_distance` | 20.500 in | 12.000 in |
| `pier[1].back.bearing_offset` | 32.500 in | -0.000 in |
| `pier[1].ahead.end_distance` | 20.500 in | 12.000 in |
| `pier[1].ahead.bearing_offset` | 32.500 in | -0.000 in |
| `group[0].girder[0].slab_offset_start` | 11.000 in | 10.067 in |
| `group[0].girder[0].slab_offset_end` | 11.000 in | 10.172 in |
| `group[0].girder[1].slab_offset_start` | 11.000 in | 9.983 in |
| `group[0].girder[1].slab_offset_end` | 11.000 in | 10.098 in |
| `group[0].girder[2].slab_offset_start` | 11.000 in | 10.036 in |
| `group[0].girder[2].slab_offset_end` | 11.000 in | 10.133 in |
| `group[0].girder[3].slab_offset_start` | 11.000 in | 9.983 in |
| `group[0].girder[3].slab_offset_end` | 11.000 in | 10.098 in |
| `group[0].girder[4].slab_offset_start` | 11.000 in | 10.067 in |
| `group[0].girder[4].slab_offset_end` | 11.000 in | 10.172 in |
| `deck.gross_depth` | 7.500 in | 8.000 in |
| `deck.left_overhang_edge_depth` | 7.000 in | 10.019 in |
| `deck.right_overhang_edge_depth` | 7.000 in | 10.020 in |
| `deck.left_overhang_taper` | 1 | 2 |
| `deck.right_overhang_taper` | 1 | 2 |
| `deck.haunch_shape` | 1 | 0 |
| `at[0].alignment.point` | (-6.000, 0.000) ft | (1526519.628, 459886.167) ft |
| `at[0].alignment.direction` | 0.0000 deg | -34.1725 deg |
| `at[0].deck.left_edge` | 15.000 ft | 15.219 ft |
| `at[0].deck.right_edge` | 15.000 ft | 15.219 ft |
| `at[1].alignment.point` | (24.000, 0.000) ft | (1526544.449, 459869.316) ft |
| `at[1].alignment.direction` | 0.0000 deg | -34.1725 deg |
| `at[1].deck.left_edge` | 15.000 ft | 15.219 ft |
| `at[1].deck.right_edge` | 15.000 ft | 15.219 ft |
| `at[2].alignment.point` | (54.000, 0.000) ft | (1526569.269, 459852.465) ft |
| `at[2].alignment.direction` | 0.0000 deg | -34.1725 deg |
| `at[2].deck.left_edge` | 15.000 ft | 15.219 ft |
| `at[2].deck.right_edge` | 15.000 ft | 15.219 ft |
