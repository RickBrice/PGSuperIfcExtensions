# Unit Test Plan: Visual Studio Tests for the Testable Core

Drafted 2026-11-20. Phase 8 of [IfcImportPlan.md](IfcImportPlan.md): move the checks that don't need a running BridgeLink into Visual Studio unit tests (the Microsoft C++ unit test framework, as WBFL's `CogoUnitTests`), so they run from Test Explorer and in CI, and are fast enough to run on every build. Not started; one of the last steps.

## What moves, and what stays in Python

**Moves** (C++ logic that needs neither a broker nor a PGSuper project):
- Mapping table loading and validation, the `extends` merge, bSDD dictionaries, `requires`.
- The table formatter, and the table <-> general IDS round trip.
- The import's target reader on small IFC models.
- Bar geometry (`BarShape.h`), the extension registry's rules, the geometry helpers.

**Stays in `Tests/ImportValidation`** (needs BridgeLink running a PGSuper project, or Python tools):
- The round trips and model imports and their scores (`compare_bridge.py`), the strand and bar comparisons.
- The export and IDS baselines, the design-value IDS checks, the round trips into source models (`check_roundtrip.py`).
- The usBridge IDS checks (ifctester), the workbook comparison, schema and unit checks.
- The export extension hooks test (`check_extension.py`).

## Why a test project can't use the DLL today

A test project calls the code it tests by linking to it. `PGSuperIfcExtensions.dll` exports only EAF's component factory (`create_class_object`, `BOOST_DLL_ALIAS` in `WBFL/Include/EAF/ComponentModule.h`). It creates the six EAF components of the DLL's object map: the data importer, the exporter, the PGSuper and PGSplice project importers, the extension agent, and the mapping table editor plugin. Those need a running BridgeLink. The classes worth unit testing are internal: no `__declspec(dllexport)`, no import library.

## The approach: a static library of the testable core

Move the testable core into a static library, `PGSuperIfcCore` (a name to settle), linked by both the DLL and a new test project, `PGSuperIfcExtensionsUnitTests`. Nothing changes for BridgeLink: the DLL links the same code.

The alternatives, and why not:
- **Link the DLL build's `.obj` files into the test project:** no restructuring, but fragile. The precompiled header has to match, and the objects drag in MFC, EAF, and the rest of the DLL.
- **Compile the sources into the test project too:** simple, but two projects then have to keep the same file lists and settings.

**What goes in the library** (to confirm when it starts):

| Files | Contents |
|---|---|
| `IfcMappingTable.cpp`, `IfcTargets.cpp`, `IfcExtensionRegistry.cpp` | tables, targets, element roles, extension registration |
| `IfcTargetReader.cpp`, `IfcTargetHints.cpp`, `IfcImportUnits.cpp` | the import's reading of targets |
| `IfcTableFormat.cpp`, `IfcTableIds.cpp` (with `IdsBuilder.h`) | formatting tables, table <-> general IDS |
| `BarShape.cpp` | bar paths, lengths, masses |

**What the library must not depend on:**
- **The application:** `CIfcMappingTable::GetTableSetting`, `SetTableSetting`, and the recent table list use `EAFGetApp()` (the registry). Move them to a file that stays in the DLL.
- **The DLL's location:** the installed standard table's path (`GetStandardTablePath`) comes from `__ImageBase`. Tests load tables from a given path or text (`CIfcMappingTable::Load`).

The library uses MFC and the WBFL, PGSuper, and IfcOpenShell headers as the DLL does (`IfcExtensionSdk.props`). The target getters reach PGSuper through the broker, but tests that don't call them don't need one.

## The tests

**First batch** (self-contained, and today checked only by hand or indirectly):

| Area | Tests |
|---|---|
| Table loading and validation | each error and its message, with its JSON path: an unknown key, target, or element role; a unit that's unknown or doesn't fit the target; a bad regex; a shared set with an element-specific property; an `extends` that isn't found or makes a cycle; invalid JSON; an unsupported format or version; a table that `"requires"` an extension, without it. The `extends` merge: locations prepended or replaced, property sets replaced or removed in place, selectors |
| bSDD dictionaries | the four URI forms, a dictionary inherited from the base table, one replaced by name, another default; an undeclared dictionary and `"bsdd"` without a default are errors |
| The table formatter | `Standard.ifcmap` and the agency tables format to themselves byte for byte (`FormatMappingTable`) |
| Table <-> general IDS | the standard table, written as a general IDS and generated back, declares the same (now `check_table_ids.py`); the binding file assigns roles and targets |
| The target reader | on small IFC models built in memory with IfcOpenShell: locations tried in order, the property's unit, else the project unit, else the location's; `parse` (number, feet-inches, regex), `map`, `list_index`, `on` (occurrence, type, material), `value_types`; a USERDEFINED selector matches the type object's `ElementType`; the hints |
| Bar geometry | `BarPath` lengths with lines and arcs, `Scaled`, `GetSignature` to 0.01 mm; `GetBarMass` against the ASTM A615 nominal masses |
| Extension registration | prefix, duplicates, an unknown role, a target without a getter, import locations for an extension target; the API version and build fingerprint refusals (with a fake extension, as the test extension does through BridgeLink) |

**Later:**
- The geometry helpers on synthetic meshes: cuts (`CutMesh`, `CutMeshSegments`), Douglas-Peucker simplification, the outline segmentation of [ImportGeometry.md](ImportGeometry.md).
- The checks done by hand: deck section outlines (constant point count, no self intersection, valid for IfcOpenShell), haunch profiles, girder layout and spacing, pier skew, end distance and bearing offset.
- Section fitting (`OutlineDimensionFitter`) on known sections. It needs PGSuper's beam factories loaded, so it's the hardest.
- The WBFL `ISectionCutTool::CreateSlabShapeWithoutHaunches`/`CreateHaunchShape` and PGSuper `IShapes` methods, in WBFL's and PGSuper's own test projects.

## Steps

1. Split the library out of the DLL project: the files above, with the application-dependent parts left in the DLL. The DLL and the validation results must be unchanged.
2. Add the test project to `PGSuperIfcExtensions.sln` and `BridgeLinkCompleteEx.sln`, after the library.
3. Write the first batch, one area at a time.
4. Run the tests from Test Explorer, and from the command line (`vstest.console.exe`, as for WBFL's tests), so a build script can run them.
5. Decide which hand checks and validation rows the tests replace; `check_table_ids.py` is the first candidate.

## Open questions

- The library's name, and whether it takes more of the import (the girder section code) once the first batch works.
- Where CI runs them, and with which IfcOpenShell build.
