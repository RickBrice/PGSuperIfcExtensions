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

bSDD dictionaries {#bsdd_dictionaries}
-----------------

Property sets, properties, classification systems, and classifications can have URIs in a bSDD dictionary (the buildingSMART Data Dictionary). A table names its dictionaries at the top of the file:

    "bsdd_dictionaries": {"usBridge": "https://identifier.buildingsmart.org/uri/usTransportation/usBridge/1/"},
    "bsdd_dictionary": "usBridge",

`bsdd_dictionaries` gives each dictionary a short name and its URI, ending with `/`. `bsdd_dictionary` is the default. A URI then names a dictionary instead of being written out:

| URI | Means |
|---|---|
| `bsdd` | The default dictionary and the entry's own name, for example `.../usBridge/1/prop/BarMark` for the property `BarMark`. |
| `bsdd:Other` | The default dictionary and another name. |
| `bsdd@name` | The dictionary with that short name, and the entry's own name. |
| `bsdd@name:Other` | That dictionary and another name. |

The standard table has the usBridge dictionary, its default, so a table that extends it writes usBridge URIs with `bsdd`. A table has the dictionaries of the table it extends, and can add its own, replace one by giving its short name another URI, or name another default. For example, a table that adds property sets from another dictionary declares that dictionary, `"bsdd_dictionaries": {"agency": "https://identifier.buildingsmart.org/uri/myAgency/bridges/2/"}`, and writes their URIs as `bsdd@agency`.

The editor doesn't show the dictionaries; they're edited in the file. When the editor copies an entry from the base table, it writes its URIs in the short form for the table's dictionaries.

Tables for an export extension
------------------------------

Another plug-in can add to PGSuper's IFC export with an export extension, for example an exchange with another program. Its table extends the standard table and names the extension in `"requires"`; it can use the extension's own element roles and targets. Such a table works only with the export command of the plug-in that provides the extension. The import, the editor's validation, and PGSuper's own export report the table as needing the extension.

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
- a URI that names a bSDD dictionary the table and the tables it extends don't have, or `bsdd` without a default dictionary
- a table that needs an export extension (`"requires"`), outside that extension's export command

A table that isn't valid stops an import or export. PGSuper doesn't fall back to the standard table, because the results would differ without notice.
