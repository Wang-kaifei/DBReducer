# Windows 上编译 DBReducer.exe（CMake，配置见 DBReducer/CMakeLists.txt）
#
# 用法:  .\build.ps1 -OutDir <编译目录> [-Toolchain <mingw64\bin>]
# 产出:  <编译目录>\bin\DBReducer.exe（全静态，不依赖任何 dll；v0.4.0 起不再有 libUtils.dll）
#
# 编译器用 msys2 的 mingw64 g++（与 v0.3.x 发布版同一套，12.2）。

param(
    [Parameter(Mandatory = $true)][string]$OutDir,
    [string]$Toolchain = "C:\msys64\mingw64\bin",
    [int]$Jobs = 8
)
$ErrorActionPreference = "Stop"
$gpp = Join-Path $Toolchain "g++.exe"
$make = Join-Path $Toolchain "mingw32-make.exe"
if (-not (Test-Path $make)) { $make = Join-Path $Toolchain "make.exe" }
foreach ($f in $gpp, $make) { if (-not (Test-Path $f)) { throw "找不到 $f" } }
if (-not (Get-Command cmake -ErrorAction SilentlyContinue)) { throw "找不到 cmake" }
$env:PATH = "$Toolchain;$env:PATH"      # cc1plus 等需要同目录的 dll

if (-not [System.IO.Path]::IsPathRooted($OutDir)) { $OutDir = Join-Path (Get-Location) $OutDir }
$OutDir = [System.IO.Path]::GetFullPath($OutDir)

# PowerShell 5.1 会把原生程序写到 stderr 的内容当成错误；这里只按退出码判断成败
$ErrorActionPreference = "Continue"
$fwd = { param($p) $p -replace '\\', '/' }      # CMake 里反斜杠会被当成转义
& cmake -S (Join-Path $PSScriptRoot "DBReducer") -B $OutDir -G "MinGW Makefiles" `
    "-DCMAKE_MAKE_PROGRAM=$(& $fwd $make)" "-DCMAKE_C_COMPILER=$(& $fwd (Join-Path $Toolchain 'gcc.exe'))" `
    "-DCMAKE_CXX_COMPILER=$(& $fwd $gpp)" 2>&1 | ForEach-Object { "$_" }
if ($LASTEXITCODE -ne 0) { throw "cmake 配置失败" }
& cmake --build $OutDir -j $Jobs 2>&1 | ForEach-Object { "$_" }
if ($LASTEXITCODE -ne 0) { throw "编译失败" }
Get-ChildItem (Join-Path $OutDir "bin") -Filter *.exe | ForEach-Object { "OK  $($_.FullName)" }
