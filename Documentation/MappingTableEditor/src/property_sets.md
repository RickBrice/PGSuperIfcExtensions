Property Sets and Quantity Sets {#property_sets}
===============================

The **Property sets** and **Quantity sets** items list what the export writes for each element role. Entries in bold are in this table. The other entries come from the table it extends.

Select an element role under Property sets or Quantity sets, then **Add** to add a set for that role.

Entries from the base table
---------------------------

Select an entry from the base table to see it. Two things can be done with it:

- **Change in This Table:** copies the entry into this table, where it can be changed. The copy replaces the base table's entry in the export.
- **Leave Out of the Export:** adds an entry to this table that leaves the base table's set out of the export. It is shown as "(left out)". Delete that entry to write the set again.

A set in this table
-------------------

| Field | Meaning |
|---|---|
| Name | The property set name, for example `Pset_ConcreteElementGeneral`. A set with the same name and element role as a base table's set replaces it. |
| Applies to | The element roles the set is written for. |
| Attach to | Where the set is written: on the element (`occurrence`), on its type object (`type`), or on its material (`material`). |
| Condition | The export option that includes the set: `classify` (written when the export classifies elements, the default), `quantities` (written when the export includes quantities), or `always`. |
| URI | The set's URI, written as its description, for example its bSDD URI. |
| Shared | Write one instance of the set, related to every element of the role. Only allowed when no property's value depends on the element. |
| Comment | Anything you want to say about the set. |
| Properties | The set's properties, in the order they are written. Add, edit, remove, and reorder them. |

**Delete from This Table** removes the set from this table. A base table's set with the same name is written again.

The Property window {#property_dialog}
-------------------

| Field | Meaning |
|---|---|
| Name | The property name. |
| Type | The IFC value type, for example `IfcPressureMeasure` or `IfcLabel`. |
| No value | The property is written with no value. It declares the property, for example one an agency requires but PGSuper doesn't have. |
| Bound to a target | The property holds a PGSuper value. The export writes the target's value, in PGSuper's units. See @ref reference for the targets. |
| An import reads the target from it | When checked, the import also reads the target from this property (after the target's own locations). Uncheck it for a property that only the export uses. |
| A constant | The property is written with this value for every element. |
| URI | The property's URI, written as its description. `bsdd` builds the bSDD URI from the set and property names. |
| Enumeration name, Enumeration values | Makes the property an enumerated property. The values are one per line. |
| Comment | Anything you want to say about the property. |

Quantity sets have the same fields. Their properties are quantities (length, area, volume, and so on).
