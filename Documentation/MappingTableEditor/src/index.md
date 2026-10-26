IFC Mapping Table Editor {#mainpage}
========================

The IFC Mapping Table Editor creates and changes IFC mapping tables. A mapping table says where PGSuper's IFC import reads each value, and which property sets, quantity sets, and classifications PGSuper's IFC export writes.

Bridge owners and their consultants put their data in different places in an IFC model. One agency puts the release strength of a girder in a property set of its own, another in a text property with the units in the property name. A mapping table describes those places, so the import finds the values without anyone changing PGSuper. The same table tells the export how to write the model the way the agency wants it.

What is in a mapping table
--------------------------

| Part | What it says | See |
|---|---|---|
| Element roles | Which IFC elements are the girders, the deck, the bearings, and so on | @ref element_roles |
| Targets | Where the import looks for each PGSuper value, in the order the places are tried | @ref targets, @ref locations |
| Property sets and quantity sets | The property sets and quantity sets the export writes for each element role, and which PGSuper value each property holds | @ref property_sets |
| Classifications | The classification systems and classification references the export writes | @ref classifications |

A table usually **extends** another table. The IFC extension installs a **standard table**, which describes the TPF/usBridge property sets that PGSuper exports. An agency table extends the standard table and contains only what is different: its own places to look for values, and the property sets it adds, changes, or leaves out. See @ref table.

Starting the editor
-------------------

The editor is a BridgeLink application installed with the IFC extension.

- **A new table:** File > New, then choose **IFC Mapping Table**. A new table extends the standard table.
- **An existing table:** File > Open, and choose the table (`.ifcmap`).
- **A table from an IDS:** Table > Generate from IDS makes a table from an agency's IDS (Information Delivery Specification). See @ref ids.

The editor window
-----------------

The tree on the left lists the parts of the table:

- **Table:** the table's name, the table it extends, and a comment.
- **Element roles:** the selector for each element role.
- **Targets:** the targets, grouped by element role, with the places the import looks for each one.
- **Property sets**, **Quantity sets**, **Classification systems**, and **Classifications:** what the export writes, grouped by element role.
- **Model:** an IFC model, opened to look at its elements and properties and to try the table on it. See @ref model.

Select an item to see and change it on the right. **Items in bold are defined in this table.** The others come from the table it extends; they are shown so you can see the whole result, and most of them can be changed in this table.

Changes are made to the table in the editor. File > Save writes the table. Table > Validate checks the table the same way the import and the export do, and lists every problem with where it is in the file.

Using a table
-------------

PGSuper and PGSplice use the mapping table configured with File > Configure BridgeLink. An IFC import or export can choose another table, and a project remembers its table. See @ref choosing_table.

A project imported from an IFC model can be exported as a copy of that model with PGSuper's design in it. The table decides which values the copy gets. See @ref round_trip.

Topics
------

- @subpage table
- @subpage element_roles
- @subpage targets
- @subpage locations
- @subpage property_sets
- @subpage classifications
- @subpage ids
- @subpage model
- @subpage choosing_table
- @subpage round_trip
- @subpage reference
