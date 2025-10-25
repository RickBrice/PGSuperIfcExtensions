The PGSuperIfcExtension is based on IfcOpenShell v0.8.0. Before building this PGSuper plug-in, the IfcOpenShell dependency needs to be installed from source and compiled.

**NOTE** If you have a version of the boost libraries already installed, consider using the version installed with IFCOS instead. Mixing versions of boost libraries leads to runtime issues.

Designate the repository directories as `safe` for git. This is a one-time setup step. The build scripts will fail if git cannot clone the repositories and check out branches.
Option 1 - Make git trust any directory
~~~
git config --global --add safe.directory *
~~~

This may work too, but I haven't tested it
~~~
git config --global --add safe.directory F:/IfcOpenShell/*
~~~

Option 2 - Make git trust the following directories (more directories may be added to this list as IfcOpenShell evolves with time)
~~~
git config --global --add safe.directory F:/IfcOpenShell
git config --global --add safe.directory F:/IfcOpenShell/_deps/mpfr
git config --global --add safe.directory F:/IfcOpenShell/_deps/mpir
git config --global --add safe.directory F:/IfcOpenShell/_deps/OpenCOLLADA
git config --global --add safe.directory F:/IfcOpenShell/_deps/occt_git
git config --global --add safe.directory F:/IfcOpenShell/_deps/cgal
git config --global --add safe.directory F:/IfcOpenshell/_deps-vs2022-x64-installed/Eigen
git config --global --add safe.directory F:/IfcOpenshell/src/ifcconvert/cityjson
~~~


Get the sources
Checkout the v0.8.0 branch and updated submodules (The cityjson submodule gets missed during the oringal clone)
~~~
git clone --recursive https://github.com/IfcOpenShell/IfcOpenShell.git
~~~

Uninstall Python and use the version that comes with IfcOpenShell. If you want to use a different version of Python, you'll need to dig into the documentation and figure it out. 

Open the `x64 Native Tools Command Prompt for VS 2022` window.
~~~
Start > Visual Studio 2022 > x64 Native Tools Command Prompt for VS 2022
~~~

This is an optional step. One of my computers has a 24 core processor. The default is to use all 24 cores, but this ends up with out of memory errors. I find that using 12 cores works well.
~~~
set IFCOS_NUM_BUILD_PROCS=12
~~~

Build the IfcOpenShell dependencies by running the following commands.
This will take a long time. 
I have experienced the scripts failing the first time I run them. Just run them a second time.

First build the Debug version of the dependencies.

~~~
cd F:\IfcOpenShell\win
build-deps.cmd vs2022-x64 Debug
~~~

Python has been installed. You need to manually uninstall it. Go to Control Panel and remove Python.

Rename the IfcOpenShell folder to IfcOpenShell_Debug. You will need to close the command prompt, and maybe even reboot.

Now we need to start over for the Release dependencies.

~~~
Start > Visual Studio 2022 > x64 Native Tools Command Prompt for VS 2022
git clone --recursive https://github.com/IfcOpenShell/IfcOpenShell.git
set IFCOS_NUM_BUILD_PROCS=12
cd F:\IfcOpenShell\win
build-deps.cmd vs2022-x64 Release
run-cmake.cmd vs2022-x64
~~~

Next, merge the debug libraries into the release folder structure.
1. Go to the folder IfcOpenShell_Debug\_deps-vs2022-x64-installed\rocksdb\lib. Make of copy of rocksdb.lib and rename it to rocksdbd.lib. Move rocksdbd.lib to IfcOpenShell\_deps-vs2022-x64-installed\rocksdb\lib.
2. Go to the folder IfcOpenShell_Debug\_deps-vs2022-x64-installed\zstd\lib. Make of copy of zstd_static.lib and rename it to zstd_staticd.lib. Move zstd_staticd.lib to IfcOpenShell\_deps-vs2022-x64-installed\zstd\lib.
3. Go to the folder IfcOpenShell_Debug\_deps-vs2022-x64-installed\HDF5-1_13_1-win64\lib. Copy all of the lib files (their names should end with _D.lib) to IfcOpenShell\_deps-vs2022-x64-installed\HDF5-1_13_1-win64\lib
4. Go to the folder IfcOpenShell_Debug\_deps-vs2022-x64-installed\OpenCOLLADA\lib\opencollada. Copy all of the lib files (their names should end with d.lib) to IfcOpenShell\_deps-vs2022-x64-installed\OpenCOLLADA\lib\opencollada

5. Go to the folder IfcOpenShell\_deps\boost_1_86_0 and run b2.exe. This will compile all of the boost libraries need by BridgeLink.


Now we are ready to build IfcOpenShell. The solution file is IfcOpenShell\_build-vs2022-x64\IfcOpenShell.sln

Select the Debug configuration and build only the IfcParse project.

Select the Release configuration and build the INSTALL project, found in the folder CMakePredefinedTargets.