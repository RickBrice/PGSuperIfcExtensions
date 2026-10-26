Choosing the Table PGSuper Uses {#choosing_table}
===============================

PGSuper and PGSplice use one configured mapping table for IFC imports and exports, and for the design-value IDS. An import or an export can use another table, and a project remembers its table.

The configured table
--------------------

File > Configure BridgeLink has an **IFC Extension** page that chooses the configured table:

- **Standard table:** the table installed with the IFC extension. It describes the TPF/usBridge property sets.
- **Another table:** a table file. Choose a recently used table from the list, by its name, or choose a file with the **...** button.

The page loads and validates the table when you choose it, and shows the table and the tables it extends, or what is wrong with it. A table that can't be used can't be configured.

The configured table is saved for the user and applies to PGSuper and PGSplice.

Importing
---------

After you choose the IFC model to import (File > Import > IFC Model Importer), the **IFC Import** window chooses the table for this import:

- **The configured table:** the window shows its name.
- **Another table:** the standard table, a recently used table, or a file chosen with the **...** button.

Choosing another table here doesn't change the configured table. The window validates the table, as the configuration page does, and won't import with a table that can't be used.

The new project remembers the table the import used.

Exporting
---------

The **Mapping Table** tab of the export options (File > Export > Bridge Model to IFC) chooses the table for the IFC model and the design-value IDS. Its default is the first of these:

1. The project's last IFC export table.
2. The table of the IFC import that made the project.
3. The configured table.

Choose **Another table** to use a different table for this export: the standard table, a recently used table, or a file. The project remembers the table of its last export, so its next export uses that table. Save the project to keep it.

Recently used tables
--------------------

The tables chosen in any of these windows, and the tables imports and exports used, are kept in one list of recently used tables, shown by the table name. The list is saved for the user.

Which table is used
-------------------

The first of these is used:

1. `/IfcMapping=<table.ifcmap>` on the command line.
2. The table chosen in the IFC Import window or on the Mapping Table tab of the export options. For an export from the command line without `/IfcMapping`, the project's table.
3. The configured table.
4. The standard table installed with the IFC extension.

The tables the chosen table extends are then loaded with it. **standard** always means the installed standard table.

PGSuper doesn't choose a table from the contents of an IFC model.

The start of every import and export log lists each table file that was used, how it was chosen, its name, and the version of the IFC extension. Support staff can tell from a log which tables produced a result.

When a table can't be used
--------------------------

A table that can't be found, can't be read, or isn't valid stops the import or export, with a message that says what is wrong, where the file is, how it was chosen, and what to do. PGSuper doesn't fall back to the standard table, because the results would differ without notice. Correct the table (Table > Validate in the editor lists the problems), or choose another table.

A project that remembers a table that has since been moved or deleted shows the problem on the Mapping Table tab; choose another table there.

Command line
------------

| Command | What it does |
|---|---|
| `/IfcMapping=<table.ifcmap>` | Uses this table for an `/IfcImport` or `/IfcExport`, in place of the configured table or the project's table. |
| `/IfcSourceModel[=<model.ifc>]` | With `/IfcExport`, writes a copy of the model the project was imported from (or of the given model) with PGSuper's design in it. The table decides which values the copy gets. See @ref round_trip. |
| `/IfcIdsToTable=<agency.ids>` | Makes a table from an IDS (`/IfcTable=` names the table to write). See @ref ids. |
| `/IfcTableToIds=<table.ids>` | Writes the configured table, or the `/IfcMapping=` table, as a general IDS. See @ref ids. |
| `/IfcFormatTable=<table.ifcmap>` | Rewrites a table file in the standard table's style, as the editor does. |

Imports and exports from the command line don't change the table a project remembers, or the list of recently used tables.
