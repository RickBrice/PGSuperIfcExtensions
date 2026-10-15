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
