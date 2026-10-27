# PGSuper IFC Extensions
IFC import and export for PGSuper/PGSplice

This project develops IFC features for [BridgeLink](https://github.com/wsdot/bridgelink) [PGSuper and PGSplice](https://github.com/wsdot/pgsuper): importing IFC bridge models into PGSuper, exporting PGSuper projects to IFC (IFC 4.3, with the usBridge property sets and classifications), mapping tables that say where PGSuper data lives in a model, round trips that update the model a project was imported from, and IDS files for checking models.

## Building

The project is built with Visual Studio 2026 (x64). It needs:
- **IfcOpenShell v0.9.0**, built from source: see [devdocs/Dependencies.md](devdocs/Dependencies.md), which points to IfcOpenShell's own build instructions.
- **WBFL, PGSuper, and BridgeLink**, in the folder named by the `ARPDIR` environment variable.
- **PROJ**, from vcpkg (below).

## PROJ for georeferencing

The [PROJ](https://proj.org/) library (a C library for georeferencing coordinate transformations) is used by PGSuperIfcExtensions.

This project uses [vcpkg](https://vcpkg.io) (in manifest mode) to pull in
PROJ. Dependencies are declared in `vcpkg.json` and restored automatically
when you build - you just need vcpkg installed and registered once.

### Quick Setup

1. Edit `setup-dev-environment.ps1`, specifically `$vcpkgDir = "C:\dev\vcpkg"`, to choose where to install vcpkg. (Don't commit your change to the script, so we can all start from the same place.)
2. Run `setup-dev-environment.ps1` (right-click -> Run with PowerShell). The script clones Microsoft's vcpkg code, runs its bootstrapper, and integrates vcpkg into Visual Studio. This is a one-time operation.
3. Open `PGSuperIfcExtensions.sln` in Visual Studio 2026 and build.

That's it. The first build will be slower than usual - vcpkg downloads and
compiles PROJ the first time only. After that, builds are normal speed.

### Troubleshooting

- **Build fails looking for proj.h / can't find PROJ**: confirm you ran
  `vcpkg integrate install` (step 2 above) - this is a one-time,
  per-machine step and won't happen automatically just from cloning.
- **Missing DLL when running the .exe outside Visual Studio**: `proj.dll`
  should be copied next to the .exe automatically as part of the build.
  If it's missing, check the project's post-build steps ran (rebuild
  rather than incremental build).

## Documentation

- **Developer documentation:** the design notes and plans are in [devdocs](devdocs/Introduction.md). Run `doxygen Doxyfile.dox` in `devdocs` to build them as HTML (`devdocs/API/html`).
- **User documentation:** in `Documentation` (e.g. the Mapping Table Editor).
- **Validation:** `Tests/ImportValidation` (see its README) runs the import and export round trips and the checks against the usBridge IDS files.

## IDS

The export can write an IDS with the project's design values (concrete strengths, strand forces, and so on), so a model made later in a model authoring tool can be checked against the design. That IDS is specific to the model. General IDS files, for checking models before they're imported, are written from mapping tables (and mapping tables from them): see [devdocs/MappingTablesDesign.md](devdocs/MappingTablesDesign.md).
