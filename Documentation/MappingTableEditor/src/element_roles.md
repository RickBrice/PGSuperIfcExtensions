Element Roles {#element_roles}
=============

An element role is a kind of IFC element that PGSuper reads or writes: the girders, the deck, the haunches, the bearings, and so on. The table's **selector** for a role says which IFC elements play it. For example, the standard table selects girders as `IfcBeam` elements with the predefined type `BEAM` in the superstructure, and the deck as the `IfcSlab` with the predefined type `FLOOR`.

The list of element roles is fixed by PGSuper (see @ref reference). A table can't add roles; it says which elements play them.

Select **Element roles** in the tree to see the roles. Roles in bold have a selector in this table; the others use the base table's selector. Select a role to see or change its selector.

The selector
------------

| Field | Meaning |
|---|---|
| Use the base table's selector | The role uses the selector of the table this table extends. |
| Define the selector in this table | The role uses the selector below, in place of the base table's. |
| Entity | The IFC entity, for example `IfcBeam`. Subtypes are included: `IfcSlab` also selects `IfcSlabStandardCase` elements. |
| Predefined type | The element's predefined type, for example `BEAM`. If the element has a type object, its predefined type is used. Leave it empty for any predefined type. |
| Classification | The identification of a classification reference the element must have. Leave it empty for any classification. |
| Attributes | Attributes the element must have, one per line, as `Name=Value`. The values match exactly, ignoring case. For example, `ObjectType=Haunch` selects the `IfcSlab` elements whose `ObjectType` is Haunch. |

An element plays the role if it matches every field that is given.

**Alternatives.** A selector can list alternatives (`any_of` in the file): an element plays the role if it matches any of them. The editor shows alternatives; they are edited in the file.

Entity and attribute names are checked against the IFC schema when the table is validated (Table > Validate).

Example
-------

The Iowa DOT models each haunch as an `IfcSlab` with the `ObjectType` Haunch. The Iowa table selects them with:

| Field | Value |
|---|---|
| Entity | `IfcSlab` |
| Attributes | `ObjectType=Haunch` |

PGSuper measures these slabs with the deck, because the haunches are cast with the deck.

Strands and reinforcing bars
----------------------------

The import reads the girders' strands (the `tendon` role) and reinforcing bars (the `rebar` role) from their geometry: the centerline of each strand or bar, from its body representation. An element with a body item for each bar (bars of one type written as one element) is several bars. A strand or bar belongs to the girder it's part of through its assemblies, else to the girder whose geometry contains it.

- **Strands:** the strand positions, harping, debonding, and extended strands. The targets of the `tendon` role give the material, the jacking force, the temporary strands, and the debond lengths.
- **Longitudinal bars:** bars along the girder, as rows referenced to the top or bottom of the girder.
- **Stirrups:** bars with vertical legs, as shear zones. A leg that extends above the girder engages the deck. Bars around the bottom strands at the stirrups are the zones' confinement bars. The stirrups near the girder ends are the splitting reinforcement; end ties that don't reach the deck are left out.
- **Deck bars:** bars inside the deck slab. PGSuper's deck reinforcement is simplified, so the longitudinal bars are imported as lump sums (bar area per width of deck): the bars that run the length of the bridge are each mat's primary reinforcement, and the shorter bars that cross a pier are its supplemental reinforcement, with their cutoffs. The cover is the depth of the bars' centroid. Transverse bars aren't PGSuper data.

The bar size is read from `rebar.bar_size`, else from the bar's nominal diameter or its geometry. The bar material is read from `rebar.grade` and `rebar.fy`. A girder gets only the strands and bars that are in the model. A deck without longitudinal bars in the model keeps the existing deck reinforcement.

The import log describes what was read for each girder, and why any strand or bar was left out.
