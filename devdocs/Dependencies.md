The PGSuperIfcExtension is based on IfcOpenShell v0.8.0. Before building this PGSuper plug-in, the IfcOpenShell dependency needs to be installed from source and compiled.

Designate the repository directories as `safe` for git. This is a one-time setup step. The build scripts will fail if git cannot clone the repositories and check out branches.
Option 1 - Make git trust any directory
~~~
git config --global --add safe.directory *
~~~

This may work too, but I haven't tested it
~~~
git config --global --add safe.directory F:/IfcOpenShell/*
~~~

Option 2 - Make git trust the following directories
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
cd IfcOpenShell
git checkout v0.8.0
git submodule update --init --recursive
~~~

The build system uses an old version of Python. We want to use the latest version of python, make sure you have it installed. I've installed it in F:\Python\Python312 (for version 3.12).

Set the following environment variable in the Visual Studio command window to prevent IfcOpenShell from installing python.

Open the `x64 Native Tools Command Prompt for VS 2022` window.
~~~
Start > Visual Studio 2022 > x64 Native Tools Command Prompt for VS 2022
set IFCOS_INSTALL_PYTHON=FALSE
~~~

Build the IfcOpenShell dependencies by running the following commands.
This will take a long time. 
When building the Debug dependencies, there will be 3 assert windows that you need to press the Ignore button.
~~~
cd F:\IfcOpenshell\win
build-deps.cmd vs2022-x64 Debug
build-deps.cmd vs2022-x64 Release
~~~

Next run batch file for cmake to create the visual studio solution file.

If you did not have IfcOpenShell install Python, you need to tell it what Python version you
have installed and where it is located.

Use the following commands
~~~
echo PY_VER_MAJOR_MINOR=312>> BuildDepsCache-x64.txt
echo PYTHONHOME=F:\Python\Python312>> BuildDepsCache-x64.txt
~~~

Ok, ready to run cmake
~~~
run-cmake.bat vs2022-x64
~~~

Fire up Visual Studio and open the solution file. The solution file is in `F:\IfcOpenShell\build-vs2022-x64\IfcOpenShell.sln`

The HDF5 libraries need to be changed to the debugging version for multiple projects. The list include IfcHouse, IfcAdvancedHouse, _ifcopenshell_wrapper, IfcGeomServer, IfcConvert. An _D needs to be appended to the HDF file name in the linker settings.

Now build the Debug and Release configurations.

Also see https://github.com/IfcOpenShell/IfcOpenShell/issues/4584

---

## Obsolete
Below are notes for things that are now obsolete. Since the v0.8.0 branch of IfcOpenShell is the bleeding edge, this information is retained here in case we need it in the future.

### Boost 1.78 is now the default library

Edit the win\build-deps.cmd and win\run-cmake.bat files. In these files search for BOOST_VERSION and change the version from 1.74.0 to 1.78.0.


### The cgal kernel is now compiling with VS2022.

The geometry_kernel_cgal and geometry_kernel_cgal_simple projects are not compatible with VS2022. Change the toolset to VS2019
~~~
Select geometry_kernel_cgal and geometry_kernel_cgal_simple
Right click, select Properties
Change configurations to "All Configurations"
Select Configuration Properties > General
Change Platform Toolset to "Visual Studio 2019 (v142)"
~~~

### City Json project are now correct

Before building, the cityjson_converter project is missing an include path. Right-click on the cityjson_convter project and select `Properties > C/C++ > General > Additional Include Directories` and add 
~~~
F:\IfcOpenshell\_deps-vs2022-x64-installed\json
~~~
