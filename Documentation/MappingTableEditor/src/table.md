The Table {#table}
=========

Select **Table** at the top of the tree to see the table's name, the table it extends, and its comment.

| Field | Meaning |
|---|---|
| Name | The name of the table. It is written to the import and export logs, so whoever reads a log can tell which table was used. |
| Extends | The table this table is based on. Choose **standard** for the standard table installed with the IFC extension, or a table file. A file path is relative to this table's file. Leave it empty for a table that stands on its own. |
| Comment | Anything you want to say about the table, such as who maintains it and what it is for. |

Extending a table
-----------------

A table that extends another table has only what is different from it:

- **Element roles:** a selector in this table replaces the base table's selector for that role. Roles without a selector in this table use the base table's.
- **Targets:** the places in this table are tried before the base table's places, unless this table replaces the base table's places. See @ref targets.
- **Property sets, quantity sets, and classifications:** an entry in this table replaces the base table's entry with the same name for the same element role. An entry can also leave the base table's entry out of the export. New entries are added. See @ref property_sets.

The standard table describes the TPF/usBridge property sets PGSuper exports. Extending it keeps an agency table small, and an agency table keeps up with new versions of the standard table without being changed.

Files
-----

Tables are JSON files with the extension `.ifcmap`. A table with the extension `.json` can be used and opened too. New tables are saved as `.ifcmap`.

| Command | What it does |
|---|---|
| File > New | A new table that extends the standard table. |
| File > Open | Opens a table file. |
| File > Save, File > Save As | Writes the table. |

The editor keeps the order of the entries in the file, comments, and anything in the file it doesn't show. It writes the file in the same style as the standard table, so a table edited in the editor and a table edited by hand look the same.

Validating the table
--------------------

Table > Validate checks the table the same way the import and the export do. The table is valid if it can be used. Otherwise each problem is listed with its place in the file, for example `targets.girder.fci[0].unit`, and what is wrong:

- a target or element role that doesn't exist
- an IFC entity, attribute, or value type that isn't in the IFC schema
- a unit that is unknown, or doesn't fit the target (a length unit for a strength)
- a regular expression that isn't valid
- a shared property set with a property whose value differs from element to element
- a base table that can't be found, or tables that extend each other in a circle

A table that isn't valid stops an import or export. PGSuper doesn't fall back to the standard table, because the results would differ without notice.
