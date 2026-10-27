# Dependencies

PGSuperIfcExtensions is built with Visual Studio 2026 (platform toolset v145, x64) and uses IfcOpenShell v0.9.0, built from source. PROJ comes from vcpkg (see the project's README.md).

## IfcOpenShell

Build IfcOpenShell and its dependencies by following IfcOpenShell's own instructions, `win/readme.md` in the IfcOpenShell repository (`%IFCOPENSHELL_DIR%\win\readme.md`). In short, from the x64 Native Tools Command Prompt for VS 2026, in `%IFCOPENSHELL_DIR%\win`:

~~~
python build-deps.py vs2026-x64 Debug
python build-deps.py vs2026-x64 Release
python run-cmake.py vs2026-x64
~~~

then build IfcOpenShell (`_build-vs2026-x64\IfcOpenShell.sln`, or `build-ifcopenshell.py` and `install-ifcopenshell.py`) in Debug and Release. Debug and Release builds of the dependencies coexist. The readme has the options, the requirements, and how to troubleshoot a build.

Use the `v0.9.0` branch (or a branch made from it). Clone with `--recursive` so the submodules come too.

### What this project expects

- **Environment variables:** `IFCOPENSHELL_DIR`, the root of the IfcOpenShell repository, and `BOOSTDIR`, the Boost that IfcOpenShell's `build-deps` installs (`%IFCOPENSHELL_DIR%\_deps\boost_<version>`). The project links Boost iostreams from `%BOOSTDIR%\stage\vs2026-x64\lib`, which `build-deps` builds. Don't mix it with another Boost installation: mixed Boost versions cause runtime problems.
- **Folders** (made by the IfcOpenShell scripts with the `vs2026-x64` generator): `_deps-vs2026-x64-installed` (the dependencies), `_build-vs2026-x64` (the build), and `_installed-vs2026-x64` (the installed libraries). `PGSuperIfcExtensions.vcxproj` refers to these names, so a different generator means changing the project.
- **Safe instances:** IfcOpenShell's CMake option `IFOPSH_SAFE_INSTANCE` (on by default since the 0.9.1 plugin ABI v2 change) changes the exported signatures of every entity class. `PGSuperIfcExtensions.vcxproj` defines `IFOPSH_SAFE_INSTANCE` to match. If IfcOpenShell is configured with `IFOPSH_SAFE_INSTANCE=OFF`, remove the define from the project, or the link fails with unresolved `instance_data` constructors.
- **Debug libraries:** on Windows, IfcOpenShell's Debug libraries and DLLs have a `_d` suffix (`ifcopenshell_parse_d.lib`, ...).
- **Python:** the IfcOpenShell build includes a Python (`%IFCOPENSHELL_DIR%\_deps\python.<version>`) with ifcopenshell; the import validation scripts in `Tests/ImportValidation` run with it (see its README).

### Git safe directories

The dependency scripts clone repositories and check out branches. If git refuses because a directory is owned by another user, mark the IfcOpenShell directories as safe, e.g.

~~~
git config --global --add safe.directory %IFCOPENSHELL_DIR%
~~~

and the same for each directory git names (under `_deps` and `src`), or trust every directory with `git config --global --add safe.directory *`.
