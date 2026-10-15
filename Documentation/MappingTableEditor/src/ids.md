Tables and IDS {#ids}
==============

An IDS (Information Delivery Specification) says what an IFC model must contain, for example which property sets each kind of element must have. Many agencies publish an IDS for their bridge models. The editor can make a mapping table from an agency's IDS, and write a table out as an IDS.

Generating a table from an IDS
------------------------------

Table > Generate from IDS makes a table from a general IDS. The generated table replaces the table in the editor as a new, untitled table; you are asked first whether to save the changes to the table in the editor.

| Field | Meaning |
|---|---|
| IDS | The IDS file (`.ids`). |
| Binding file | Optional. A JSON file that says which element role each IDS specification is for, and which target each property holds, where the IDS doesn't say. |
| Extend the standard table | The generated table extends the standard table and has only what differs from it. Otherwise the table stands on its own. |

When the table is generated, a report lists how many properties got a target and how many didn't, and what the generator found. Properties without a target are written with no value. Check the generated table, give targets to the properties that hold PGSuper values (see @ref property_sets), and validate it (Table > Validate) before using it.

The same can be done from the command line:

    BridgeLink.exe /IfcIdsToTable=<agency.ids> <template.pgt> /IfcTable=<table.json> [/IfcBinding=<binding.json>] [/IfcTableExtends=standard|none] [/IfcLog=<file.log>]

Writing a table as an IDS
-------------------------

Table > Write as General IDS writes the table as a general IDS (`.ids`): one specification for each element role that has requirements, applying to the role's elements through its selector. Properties and quantities bound to a target are required, constants are required with their value, properties with no value are optional, and classifications are required. Material property sets aren't written, because an IDS checks an element's own and its type's property sets, not its material's. Use it to check exported models with an IDS checker, or to tell modelers what PGSuper writes and reads.

An IDS can't say everything a table says. What it can't say (targets, where a set is attached, conditions, and so on) is written in the instructions of each facet, so a table generated from the IDS gets it back.

The same can be done from the command line. It writes the active table, or the table given with `/IfcMapping=`:

    BridgeLink.exe /IfcTableToIds=<table.ids> <template.pgt> [/IfcMapping=<table.json>]

Design-value IDS
----------------

The IFC export can also write a design-value IDS: an IDS with the values PGSuper designed, for checking a model against the design. That IDS takes its property locations from the active table. It is written by the export, not by the editor.
