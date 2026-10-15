Working with a Model {#model}
====================

The editor can open an IFC model to look at its elements and properties, to add locations by picking properties, and to try the table on the model. The model is only read; it isn't changed.

Opening a model
---------------

Table > Open Model opens an IFC 4.3 model (`.ifc`). A zipped model (`.zip`) must be extracted first. The **Model** item at the bottom of the tree shows the model's file name.

Looking at the elements
-----------------------

Select **Model** in the tree.

| Field | Meaning |
|---|---|
| Element role | The elements of this role, selected with the table being edited. A role with no elements may have a selector that doesn't fit the model; see @ref element_roles. |
| Element | One of those elements. |
| Properties of the element | The element's own properties, its type's, and its material's, with their values. |
| Target | A target of the role. The pane shows what the table reads for it from the element: the value in US and SI units and the location it came from, every location that was tried, and the properties that look like they hold the value. |

Adding a location from the model
--------------------------------

Choose the target, then double-click a property in the list, or select it and click **Add as Location**. The Location window opens with the property filled in. Add a unit or a parse rule if the value needs one (see @ref locations), and click OK. The location is added to the target, and the pane shows what the table now reads.

Trying the table on the model
-----------------------------

Table > Try on Model reads every target from every element of each element role, as an import would, without importing. The report lists:

- for each target, how many elements it was found for, the first value, and the locations that were tried
- the properties that look like they hold targets that weren't found
- the element roles with no elements in the model
- the problems found while reading, such as a text that couldn't be read as a number

The table must be valid to try it (Table > Validate).

Use the report to see what an import would get before running one, and to find what is still missing from the table.
