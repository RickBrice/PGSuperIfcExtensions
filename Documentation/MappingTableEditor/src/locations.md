Locations {#locations}
=========

A location is one place in the model where the import looks for a target's value. Add or edit a location on the target's page (see @ref targets) to open the Location window.

Where the value is
------------------

| Choice | Fields | Meaning |
|---|---|---|
| Property | Property set, Property or attribute, Property sets of | A property, by the exact names of its property set and property. Single, list, and enumerated properties can be read. |
| Attribute | Property or attribute | An attribute of the element, for example `ObjectType` or `Tag`. |
| Type attribute | Property or attribute | An attribute of the element's type object, for example the type's `Name`. |
| Classification | Classification system, Identification, Read | A classification reference of the element. Give the system and the identification to pick one reference, or leave them empty for any. **Read** says what to read from it: its Identification or its Name. |

**Property sets of** says whose property sets are read:

| Choice | Meaning |
|---|---|
| the element, then its type | The element's own property sets, then its type's (the default). An element's own property set overrides its type's, as IFC intends. |
| the element's type only | Only the type object's property sets. |
| the element's material | The property sets of the element's material. |

Reading the value
-----------------

| Field | Meaning |
|---|---|
| Unit of a plain number | The unit of a value that has no unit of its own: a plain number, or a number in text. Give it for a value like `5.0` meaning 5.0 ksi. A unit in the model always wins. A target that needs a unit can't read a plain number without one. |
| Parse text as | How to read a number from text: **number**, **feet and inches** (such as `8'-6"` or `8.5"`), or **regex** (a regular expression). **(no parsing)** reads the value as it is. |
| Regex, Group | The regular expression and the group that holds the value, for text such as `Fc=8.0 ksi`. |
| List index | Which item of a list value to read. Leave it empty for the whole value. |
| Value types | The IFC value types that are accepted, for example `IfcPressureMeasure`. Leave it empty for any value type that fits the target. |
| Map | Model values and the target values they stand for, one per line: `model value = target value`. For example, `Fixed = true` for a fixity text read into `bearing.fixed_x`. Text matches exactly, ignoring case. |
| Comment | Anything you want to say about the location. |

Units
-----

Values are converted to PGSuper's units as soon as they are read. A measure value with a unit of its own, or with the model's project unit, needs no unit in the location.

| Kind of target | Units |
|---|---|
| Stress | `Pa`, `kPa`, `MPa`, `psi`, `ksi` |
| Length | `mm`, `cm`, `m`, `in`, `ft` |
| Angle | `rad`, `deg` |
| Force | `N`, `kN`, `lbf`, `kip` |

A unit that doesn't fit the target, such as `in` for a concrete strength, is reported when the table is validated.

Examples
--------

| Model | Target | Location |
|---|---|---|
| Iowa DOT: `IaDOT_PPCB` . `5_Final Concrete Strength, Fc` = 5.0 | `girder.fc` | Property; Unit of a plain number: `ksi` |
| PennDOT: `_PS Concrete Beams` . `Concrete Strength at 28 days (ksi)` = "8" (text) | `girder.fc` | Property; Parse text as: number; Unit of a plain number: `ksi` |
| Bearings with a text property `Fixity` = "Fixed" | `bearing.fixed_x` | Property; Map: `Fixed = true`, `Expansion = false` |

The easiest way to add a location is to pick the property from a model. See @ref model.
