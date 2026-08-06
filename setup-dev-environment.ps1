<#
.SYNOPSIS
    One-time developer machine setup for this repo's C++ dependencies (vcpkg).

.DESCRIPTION
    Clones and bootstraps vcpkg if it isn't already present, then registers
    it with Visual Studio via 'vcpkg integrate install'. Safe to re-run.

.USAGE
    Right-click this file -> Run with PowerShell
    (or from a PowerShell prompt:  .\setup-dev-environment.ps1)
#>

$ErrorActionPreference = "Stop"
$vcpkgDir = "C:\dev\vcpkg"

Write-Host "=== Dev environment setup ===" -ForegroundColor Cyan

# 1. Clone vcpkg if it isn't already there
if (Test-Path "$vcpkgDir\.git") {
    Write-Host "vcpkg already present at $vcpkgDir - skipping clone."
} else {
    Write-Host "Cloning vcpkg to $vcpkgDir ..."
    git clone https://github.com/microsoft/vcpkg $vcpkgDir
}

# 2. Bootstrap it if the executable isn't built yet
if (Test-Path "$vcpkgDir\vcpkg.exe") {
    Write-Host "vcpkg.exe already built - skipping bootstrap."
} else {
    Write-Host "Bootstrapping vcpkg ..."
    & "$vcpkgDir\bootstrap-vcpkg.bat"
}

# 3. Register vcpkg with Visual Studio (system/user-wide, one time)
Write-Host "Running vcpkg integrate install ..."
& "$vcpkgDir\vcpkg.exe" integrate install

Write-Host ""
Write-Host "=== Setup complete ===" -ForegroundColor Green
Write-Host "Next: open the .sln in Visual Studio 2022 and build."
Write-Host "The first build will take longer than usual - vcpkg needs to"
Write-Host "download and compile PROJ before your code can build against it."
