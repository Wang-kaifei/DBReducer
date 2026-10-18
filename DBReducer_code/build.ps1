# 用 msys2 的 MinGW-w64 g++ 直接编译 DBReducer（不依赖 cmake/make）
# 编译参数与 DBReducer/CMakeLists.txt、utils/CMakeLists.txt 保持一致；
# 和 CMake 一样直接链接 libUtils.dll 本体（不是导入库）。
#
# 用法:
#   .\build.ps1 -OutDir ..\dist\updated
#       从源码编 libUtils.dll + DBReducer.exe，运行时 dll 取自工具链
#   .\build.ps1 -OutDir ..\dist\baseline -UtilsDll <某目录\libUtils.dll>
#       不编 utils，链接给定的 libUtils.dll；该 dll 及其同目录的运行时 dll 一起拷到产出目录。
#       用发布版目录下的 libUtils.dll 时，产出的 exe 与发布版只差 PE 时间戳和校验和（2026-10-06 实测）。
#
# 产出: OutDir\DBReducer.exe + libUtils.dll + 运行时 dll，可整个目录拷走使用。

param(
    [Parameter(Mandatory = $true)][string]$OutDir,
    [string]$UtilsDll = "",
    [string]$Toolchain = "C:\msys64\mingw64\bin"
)
$ErrorActionPreference = "Stop"

$src = $PSScriptRoot
$gpp = Join-Path $Toolchain "g++.exe"
if (-not (Test-Path $gpp)) { throw "找不到 $gpp" }
$env:PATH = "$Toolchain;$env:PATH"      # cc1plus 等需要同目录的 dll

if (-not [System.IO.Path]::IsPathRooted($OutDir)) { $OutDir = Join-Path (Get-Location) $OutDir }
$OutDir = [System.IO.Path]::GetFullPath($OutDir)
$obj = Join-Path $OutDir "obj"
New-Item -ItemType Directory -Force $obj | Out-Null

# 版本信息写进 exe（状态文件的 tool_build 字段）
$hash = (& git -C $src describe --always --dirty 2>$null)
if (-not $hash) { $hash = "nogit" }
$build = "{0} {1}" -f (Get-Date -Format "yyyy-MM-dd"), $hash

function Run($argv) {
    # PowerShell 5.1 会把原生程序写到 stderr 的警告包装成错误记录；这里只按退出码判断成败
    $old = $ErrorActionPreference; $ErrorActionPreference = "Continue"
    $out = & $gpp @argv 2>&1 | ForEach-Object { "$_" }
    $code = $LASTEXITCODE; $ErrorActionPreference = $old
    if ($code -ne 0) { $out | Write-Host; throw "g++ 失败: $($argv -join ' ')" }
    $warn = @($out | Where-Object { $_ -match 'warning:' }).Count
    if ($warn -gt 0) { Write-Host ("  {0}: {1} 条警告" -f (Split-Path $argv[-1] -Leaf), $warn) }
}

$dll = Join-Path $OutDir "libUtils.dll"
if ($UtilsDll) {
    # ---- 链接现成的 libUtils.dll ----
    $UtilsDll = (Resolve-Path $UtilsDll).Path
    $runtimeDir = Split-Path $UtilsDll
    if ($UtilsDll -ne $dll) { Copy-Item $UtilsDll $dll -Force }
} else {
    # ---- utils -> libUtils.dll ----
    $runtimeDir = $Toolchain
    $utilFlags = @("-O3", "-Wall", "-std=gnu++14")
    $uo = @()
    foreach ($f in "IOUtil.cpp", "StringUtil.cpp") {
        $o = Join-Path $obj ($f -replace '\.cpp$', '.o')
        Run ($utilFlags + @("-c", (Join-Path $src "utils\$f"), "-o", $o))
        $uo += $o
    }
    Run (@("-O3", "-Wall", "-shared", "-o", $dll) + $uo)
}

# ---- DBReducer.exe ----
# 版本信息走生成的头文件：PowerShell 5.1 传给原生程序的 -D"带引号的值" 会被吞掉引号
Set-Content -Path (Join-Path $obj "dbr_build_info.h") -Value "#define DBR_BUILD `"$build`"" -Encoding ASCII
$dbrFlags = @("-fopenmp", "-O3", "-Wall", "-std=gnu++20", "-I$obj")
# 先按 CMakeLists.txt 里的顺序（链接顺序影响二进制），新增的源文件排在后面
$cmakeOrder = "main.cpp", "Rebuild.cpp", "Param.cpp", "Rebuild-spec.cpp"
$all = Get-ChildItem (Join-Path $src "DBReducer") -Filter *.cpp
$dbrSrc = @($cmakeOrder | ForEach-Object { $n = $_; $all | Where-Object Name -eq $n }) +
          @($all | Where-Object { $cmakeOrder -notcontains $_.Name } | Sort-Object Name)
$do = @()
foreach ($f in $dbrSrc) {
    $o = Join-Path $obj ($f.BaseName + ".o")
    Run ($dbrFlags + @("-c", $f.FullName, "-o", $o))
    $do += $o
}
$exe = Join-Path $OutDir "DBReducer.exe"
# 末尾这串系统库是 CMake 在 MinGW 下默认附加的，影响导入表顺序；去掉它产出就和发布版不一致
$cmakeStdLibs = "-lkernel32", "-luser32", "-lgdi32", "-lwinspool", "-lshell32", "-lole32", "-loleaut32", "-luuid", "-lcomdlg32", "-ladvapi32"
Run (@("-fopenmp", "-O3", "-Wall", "-fopenmp") + $do + @("-o", $exe, $dll) + $cmakeStdLibs)

# ---- 运行时 dll ----
foreach ($r in "libgcc_s_seh-1.dll", "libstdc++-6.dll", "libwinpthread-1.dll", "libgomp-1.dll") {
    Copy-Item (Join-Path $runtimeDir $r) $OutDir -Force
}

"OK  $exe"
"    build = $build"
