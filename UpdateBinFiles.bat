REM - Script to prepare for Release
REM - Stages everything the IFC extensions need at run time in bin\Extensions\Experimental\x64.
REM - BridgeLink's BuildRelease.bat copies bin to \ARP_RELEASE\PGSuper\<version>\bin, where the
REM - installer (BridgeLink\Installer\ExperimentalExtensions.wxs) picks the files up.
REM - Keep this list and ExperimentalExtensions.wxs in step.

SET BINTARGET=bin\Extensions\Experimental\x64
SET REGFREECOM=%ARPDIR%\BridgeLink\RegFreeCOM
SET SOURCE=%REGFREECOM%\x64\Release

if not exist %BINTARGET% mkdir %BINTARGET%
if not exist %BINTARGET%\MappingTables mkdir %BINTARGET%\MappingTables
if not exist %BINTARGET%\proj mkdir %BINTARGET%\proj

REM - Experimental Extensions (manifests register the COM components, no registry entries)
xcopy /Y /d %SOURCE%\PGSuperIfcExtensions.dll           %BINTARGET%\
xcopy /Y /d %SOURCE%\IfcExtensions.Manifest.PGSuper     %BINTARGET%\
xcopy /Y /d %SOURCE%\IfcExtensions.Manifest.BridgeLink  %BINTARGET%\
xcopy /Y /d %SOURCE%\IfcImportTemplate.pgt              %BINTARGET%\
xcopy /Y /d %SOURCE%\IfcImportTemplate.spt              %BINTARGET%\

REM - Mapping tables (the standard table is found next to the DLL, MappingTables\Standard.ifcmap)
xcopy /Y /d %SOURCE%\MappingTables\*.ifcmap             %BINTARGET%\MappingTables\

REM - IfcOpenShell (OpenCASCADE is linked statically)
for %%f in (ifcopenshell_plugin.dll ifcopenshell_parse.dll ifcopenshell_parse_schema_ifc4x3_add2.dll ifcopenshell_geometry.dll ifcopenshell_geometry_kernel_opencascade.dll ifcopenshell_geometry_mapping_ifc4x3_add2.dll) do xcopy /Y /d %SOURCE%\%%f %BINTARGET%\

REM - PROJ (georeferencing) and the DLLs it needs, from vcpkg
for %%f in (proj_9.dll sqlite3.dll tiff.dll libcurl.dll z.dll jpeg62.dll liblzma.dll) do xcopy /Y /d %SOURCE%\%%f %BINTARGET%\

REM - PROJ data: proj.db is the EPSG database (read from <app folder>\proj\proj.db); the vcpkg CMake files aren't needed
for %%f in (proj.db proj.ini CH GL27 ITRF2000 ITRF2008 ITRF2014 ITRF2020 nad.lst nad27 nad83 other.extra world deformation_model.schema.json projjson.schema.json triangulation.schema.json copyright) do xcopy /Y /d %SOURCE%\proj\%%f %BINTARGET%\proj\

REM - Report anything missing, so the installer build doesn't fail later
SET MISSING=0
for %%f in (PGSuperIfcExtensions.dll IfcExtensions.Manifest.PGSuper IfcExtensions.Manifest.BridgeLink IfcImportTemplate.pgt IfcImportTemplate.spt MappingTables\Standard.ifcmap ifcopenshell_plugin.dll ifcopenshell_parse.dll ifcopenshell_parse_schema_ifc4x3_add2.dll ifcopenshell_geometry.dll ifcopenshell_geometry_kernel_opencascade.dll ifcopenshell_geometry_mapping_ifc4x3_add2.dll proj_9.dll sqlite3.dll tiff.dll libcurl.dll z.dll jpeg62.dll liblzma.dll proj\proj.db) do if not exist %BINTARGET%\%%f (echo MISSING: %BINTARGET%\%%f & SET MISSING=1)
if "%MISSING%"=="1" echo One or more IFC extension files are missing. Build PGSuperIfcExtensions (Release, x64) first.
