The PGSuperIfcExtension is based on IfcOpenShell v0.8.0. Before building this PGSuper plug-in, the IfcOpenShell dependency needs to be installed from source and compiled.

The basic instructions are at https://blenderbim.org/docs-python/ifcopenshell/installation.html.

Make the following adjustments after cloning the repository.

Checkout the v0.8.0 branch and updated submodules (The cityjson submodule gets missed during the oringal clone)
~~~
$ git clone --recursive https://github.com/IfcOpenShell/IfcOpenshell.git
$ cd IfcOpenshell
$ git checkout v0.8.0
$ git submodule update --init --recursive
~~~

Open the `x64 Native Tools Compand Prompt for VS 2022` window. Build the dependencies and run cmake.
~~~
Start > Visual Studio 2022 > x64 Native Tools Command Prompt for VS 2022
$ F:
$ cd F:\IfcOpenshell\win
$ build-deps.cmd vs2022-x64 Debug
$ build-deps.cmd vs2022-x64 Release
$ run-cmake.bat vs2022-x64
~~~

Open the solution file in Visual Studio:

~~~
$ ..\build-vs2022-x64\IfcOpenShell.sln
~~~

Before building, the cityjson_converter project is missing an include path. Right-click on the cityjson_convter project and select `Properties > C/C++ > General > Additional Include Directories` and add 
~~~
F:\IfcOpenshell\_deps-vs2022-x64-installed\json
~~~


The HDF5 libraries need to be changed to the debugging version for multiple projects. The list include IfcHouse, IfcAdvancedHouse, _ifcopenshell_wrapper, IfcGeomServer, IfcConvert. An _D needs to be appended to the HDF file name in the linker settings.

Now build the Debug and Release configurations.

---

## Obsolete

Boost 1.78 is now the default library

Edit the win\build-deps.cmd and win\run-cmake.bat files. In these files search for BOOST_VERSION and change the version from 1.74.0 to 1.78.0.


The cgal kernel is now compiling with VS2022.

The geometry_kernel_cgal and geometry_kernel_cgal_simple projects are not compatible with VS2022. Change the toolset to VS2019
~~~
Select geometry_kernel_cgal and geometry_kernel_cgal_simple
Right click, select Properties
Change configurations to "All Configurations"
Select Configuration Properties > General
Change Platform Toolset to "Visual Studio 2019 (v142)"
~~~
