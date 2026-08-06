# PGSuper IFC Extensions
Experimental features for PGSuper/PGSplice

This project is a sandbox for developing experimental IFC features for [BridgeLink](github.com/wsdot/bridgelink) [PGSuper and PGSplice](github.com/wsdot/pgsuper). 

## PROJ for georeferencing

The [PROJ](https://proj.org/) library (a C library for general georeferencing coordinate transformations) is used by PGSuperIfcExtensions.

This project uses [vcpkg](https://vcpkg.io) (in manifest mode) to pull in
PROJ. Dependencies are declared in `vcpkg.json` and restored automatically
when you build - you just need vcpkg installed and registered once.

### Quick Setup

1. Edit the `setup-dev-enviornment.ps1`, specifically `$vcpkgDir = "C:\dev\vcpkg"` to direct where you want to install vcpkg. (Don't commit your change to the script so we can all start from the same place).
2. Run `setup-dev-environment.ps1` (right-click -> Run with PowerShell). This scripts clones Microsoft's vcpkg code, runs it bootstrapper, and then integrates vcpkg into VisualStudio. This is a one-time operation.
3. Open the PGSuperIfcExtensions `.sln` in Visual Studio 2022 and build.

That's it. The first build will be slower than usual - vcpkg downloads and
compiles PROJ the first time only. After that, builds are normal speed.

### Troubleshooting

- **Build fails looking for proj.h / can't find PROJ**: confirm you ran
  `vcpkg integrate install` (step 2/3 above) - this is a one-time,
  per-machine step and won't happen automatically just from cloning.
- **Missing DLL when running the .exe outside Visual Studio**: `proj.dll`
  should be copied next to the .exe automatically as part of the build.
  If it's missing, check the project's post-build steps ran (rebuild
  rather than incremental build).


## Things to remember for the IDS

* Beams must be typed with IfcBeamType
* Beams must have Pset_PrecastConcreteELementGeneral property set with
    * ReleaseStrength
    * (final strength?)
    * DesignLocationNumber and its expected to be "Span i, Girder j" so use a reg expression
* Beams must be classified with usBridge_GirderPrestressedConcrete from US Bridge data dictionary on bSDD
* Beams must have Pset_ConcreteElementGeneral with
    * AssemblyPlace = FACTORY
    * CastingMethod = PRECAST


We want to have an input IDS to validate models before importing them.

We also want to output an IDS that contains things like the design values for concrete strength, strand forces, camber, etc so the proper design results can be evaluated against the final model which will be created in a model authoring tool

