Choosing the Table PGSuper Uses {#choosing_table}
===============================

PGSuper and PGSplice use one mapping table for every IFC import and export, and for the design-value IDS.

Options > IFC Mapping Table
---------------------------

In PGSuper or PGSplice, Options > IFC Mapping Table chooses the table:

- **Standard table:** the table installed with the IFC extension. It describes the TPF/usBridge property sets.
- **Agency table:** a table file. The window loads and validates the table when you choose it, and shows the table and the tables it extends, or what is wrong with it. A table that can't be used isn't saved as the setting.

The setting is saved for the user and applies to PGSuper and PGSplice.

Which table is used
-------------------

The first of these is used:

1. `/IfcMapping=<table.json>` on the command line.
2. The table chosen with Options > IFC Mapping Table.
3. The standard table installed with the IFC extension.

The tables the chosen table extends are then loaded with it. **standard** always means the installed standard table.

The start of every import and export log lists each table file that was used, how it was chosen, its name, and the version of the IFC extension. Support staff can tell from a log which tables produced a result.

When a table can't be used
--------------------------

A table that can't be found, can't be read, or isn't valid stops the import or export, with a message that says what is wrong, where the file is, how it was chosen, and what to do. PGSuper doesn't fall back to the standard table, because the results would differ without notice. Correct the table (Table > Validate in the editor lists the problems), or choose another table with Options > IFC Mapping Table.

Command line
------------

| Command | What it does |
|---|---|
| `/IfcMapping=<table.json>` | Uses this table for an `/IfcImport` or `/IfcExport`, in place of the setting. |
| `/IfcIdsToTable=<agency.ids>` | Makes a table from an IDS (`/IfcTable=` names the table to write). See @ref ids. |
| `/IfcTableToIds=<table.ids>` | Writes the active table, or the `/IfcMapping=` table, as a general IDS. See @ref ids. |
| `/IfcFormatTable=<table.json>` | Rewrites a table file in the standard table's style, as the editor does. |
