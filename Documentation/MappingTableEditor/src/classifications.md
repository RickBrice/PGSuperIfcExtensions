Classifications {#classifications}
===============

The export can classify each element, for example with the usBridge classification. The table lists the classification systems and the classification references written for each element role.

Classification systems
----------------------

Select **Classification systems** and **Add** to add a system, or select a system to change it.

| Field | Meaning |
|---|---|
| Name | The system's name, for example `usBridge`. Classifications refer to the system by name. |
| Source | Who publishes the system. |
| Edition, Edition date | The system's edition and its date. |
| Specification | Where the system is described, for example a URL. |
| Comment | Anything you want to say about the system. |

Classifications
---------------

Classifications are listed by element role. Select an element role and **Add** to add a classification for it, or select a classification to change it. Classifications from the base table can be changed in this table or left out of the export, the same as property sets (see @ref property_sets).

| Field | Meaning |
|---|---|
| Applies to | The element roles the classification is written for. |
| Condition | The export option that includes it: `classify` (the default), `quantities`, or `always`. |
| System | The classification system. |
| Identification | The classification's identification, for example `usBridge_GirderPrecastConcrete`. |
| Name | The classification's name. |
| Location | Where the classification is described, for example its bSDD URI. |
| Comment | Anything you want to say about the classification. |

**Delete from This Table** removes the entry from this table.

Classifications in the import
-----------------------------

A classification reference can also be a place where the import finds a value, and an element role's selector can require a classification. See @ref locations and @ref element_roles.
