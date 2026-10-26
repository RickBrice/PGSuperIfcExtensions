Round trips {#round_trip}
===========

A project imported from an IFC model can be exported as a copy of that model with PGSuper's design in it. The model itself isn't changed. The copy keeps everything in the model, and every element keeps its GlobalId, so the modeler can see exactly what PGSuper changed (for example with IfcOpenShell's IfcDiff) and bring the design into their model.

Exporting a copy of the model
-----------------------------

1. Import the IFC model (File > Import > IFC Model Importer). The project remembers the model: where it is, and which of its elements became which girders, piers, bearings, and so on. Save the project to keep it.
2. Design in PGSuper.
3. Export (File > Export > Bridge Model to IFC). At the top of the **IFC** tab, under **Round Trip**, **Update a copy of the IFC model this project was imported from** is checked when the model is where it was and hasn't changed since the import. The box below it shows the model. While it's checked, the tab's other options are disabled: they're for exporting PGSuper's model, and a round trip uses their defaults. The Mapping Table and IDS tabs apply.
4. Choose the file to write. The default name is the model's name followed by the project's name. The export won't replace the model itself.

If the model has moved or changed, the box says so and the option is off. Use **...** to locate it; the model you choose is used as it is.

Clear the option to export PGSuper's model instead, as for a project that wasn't imported from IFC.

What the copy gets from PGSuper
-------------------------------

**On the model's elements** (the project, site, bridge, bridge parts, piers and abutments, girders, deck, exterior barriers, bearings, and pier diaphragms that the import used), the property sets, quantity sets, material properties, and classifications of the mapping table the export uses (see @ref choosing_table and @ref property_sets):

| Property in the table | In the copy |
|---|---|
| Bound to a target | PGSuper's value. It replaces the model's value, in the unit of the model's property (or the model's project unit). If the model has the property with another kind of value (for example text where PGSuper writes a pressure), PGSuper's property replaces it. |
| A constant | Added where the model doesn't have the property. The model's value is kept. |
| No value | Added where the model doesn't have the property. The model's value is kept. |

A classification of the table is added where the element doesn't have it.

The model's names, geometry, placements, types, and everything else on the element are kept. Where an element shares a property set, a type, or a material with other elements and needs other values, the element gets its own copy: for example, girders whose concrete strength differs from the strength of the material they share get one copy of the material with their strength.

**Elements PGSuper designs that the model doesn't have:**

- The strands and reinforcing bars of a girder that has none in the model. A girder that has strands or bars in the model keeps the model's.
- The bearings at a girder end that has none in the model.
- The pier diaphragms at a pier that has none in the model. Diaphragms the model has at a pier get PGSuper's diaphragm values instead.

They are placed in the model's coordinates and units, related to the model's elements (for example a girder's strands are part of the model's girder), and use the model's representation contexts and classification systems.

What the copy doesn't get yet
-----------------------------

- Changes to the geometry and placement of the model's elements. The model's geometry is kept, even where PGSuper's differs.
- Elements added in PGSuper that the model doesn't have, such as a girder or a span. They're listed in the export log.
- Elements removed in PGSuper. The model's are kept.
- The alignment and the georeferencing. The model's are kept.

What changed
------------

The export says how many elements were updated and added, and writes the details to a log beside the IFC file (`<file>.export.log`): each property changed (with the old and new values), each property and classification added, each element added, and anything the export left alone and why.

To compare the copy with the model, IfcOpenShell's IfcDiff lists the elements added, removed, and changed, and what changed on each:

    python ifcdiff.py <model.ifc> <copy.ifc> -o diff.json -r "attributes type property container aggregate classification"

Add `geometry` to the list to compare shapes too.

The design-value IDS
--------------------

When the IDS tab asks for a design-value IDS, it's written for the copy: the girders are identified by their GlobalIds (the copy keeps the model's girder names), and it checks the values PGSuper wrote. It doesn't require PGSuper's girder materials and types, which the copy doesn't have, nor PGSuper's strands and bars when a girder keeps the model's. See @ref ids.

Command line
------------

    BridgeLink.exe /IfcExport=<copy.ifc> <project.pgs> /IfcSourceModel [/IfcIds=<design.ids>]

`/IfcSourceModel` updates a copy of the model the project remembers, which must be where it was and unchanged since the import. `/IfcSourceModel=<model.ifc>` uses the given model instead. The export log lists every change.
