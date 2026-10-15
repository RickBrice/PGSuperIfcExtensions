Targets {#targets}
=======

A target is a PGSuper value with a fixed name, such as `girder.fci` (the girder concrete strength at release) or `bearing.fixed_x`. Each target belongs to an element role and has a kind of value: a stress, a length, a count, text, and so on. The list of targets is fixed by PGSuper (see @ref reference). A table can't add targets; it says where each target is found in an IFC model.

Select **Targets** in the tree to see the targets, grouped by element role. Targets in bold have places to look in this table.

Where the import looks
----------------------

Each target has an ordered list of **locations**, places in the model where the value may be: a property, an attribute, or a classification. For each element, the import tries the locations in order. **The first location that has a value is used**, and the import log says which location it was.

Select a target to see its locations:

- **Locations in this table, in the order they're tried:** add, edit, remove, and reorder this table's locations. See @ref locations.
- **Base table locations:** the locations of the table this table extends, in the order they are tried.
- **Replace the base table's locations:** when checked, only this table's locations are tried. Otherwise this table's locations are tried first, then the base table's.

A target also reads the exported properties that hold it (see @ref property_sets), unless the property says the import doesn't read it. When a target has no locations of its own, those properties are where the import looks.

When a value isn't found
------------------------

The import doesn't guess where a value is: it looks only in the table's locations. When a target isn't found, the import log says so, and suggests properties of the element that look like they hold the value. The suggestions help whoever writes the table; use @ref model to see the same suggestions without running an import.

Values that differ from the geometry
------------------------------------

Some targets are measured from the geometry as well, for example the deck's gross depth. The geometry wins, and the import log reports when the property and the geometry differ.
