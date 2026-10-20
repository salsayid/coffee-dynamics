param(
    [Parameter(Mandatory = $true)][string]$Script,
    [string]$Function = "main()",
    [int]$TimeoutMin = 120,
    [switch]$Gui
)

# ansa64.bat mangles any argument containing "|", so set its environment here and call the exe directly
$ansaDir = "C:\Users\aalsayid\AppData\Local\Apps\BETA_CAE_Systems\ansa_v24.1.6\"
$shared = "C:\Users\aalsayid\AppData\Local\Apps\BETA_CAE_Systems\shared_v24.1.6\"
$env:CAE_TOOLS_VERSION = "24.1.6"
$env:ANSA_EXEC_DIR = $ansaDir
$env:ANSA_EXEC_PATH = $ansaDir
$env:BETA_SHARED_DIR = $shared
$env:FFMPEG_EXT_DIR = $shared
$env:QTDIR = "${shared}win64"
$env:QT_PLUGIN_PATH = "${shared}win64\plugins"
$env:QTWEBENGINE_DISABLE_SANDBOX = "1"
$env:CT_LIB_PATH = "${shared}lib\CT_2022_SP8"
$env:JT_OPEN_PATH = "${shared}win64"
$env:RHINO_LIB_PATH = "${shared}win64\cadexe"
$env:U3DPDF_LIB_PATH = "${shared}win64"
$env:ASSETS_DIR = "${shared}assets"
$env:PRE_BETA_LD_LIBRARY_PATH = $env:PATH
$env:PYTHONHOME = "${shared}python\win64\"
$env:ANSA_HOME = "${ansaDir}config\"
$env:ANSA_CARDS_HELP = "${ansaDir}config\ANSA_CARDS_HELP\"
$env:ANSA_HELP = "${ansaDir}config\helpfiles\"
$env:ANSA_DOCS = "${ansaDir}docs\"
$env:HDF5_DISABLE_VERSION_CHECK = "2"
if (-not $env:HOME) { $env:HOME = "$env:HOMEDRIVE$env:HOMEPATH" }

$exe = "${shared}win64\ansa_win64.exe"
$scriptPath = (Resolve-Path $Script).Path
$logDir = Join-Path (Split-Path $PSScriptRoot) "logs"
New-Item -ItemType Directory -Force $logDir | Out-Null
$name = [IO.Path]::GetFileNameWithoutExtension($scriptPath)
$log = Join-Path $logDir "ansa_$name.log"
$err = Join-Path $logDir "ansa_$name.err"

$argList = @("-nolauncher")
if (-not $Gui) { $argList += "-nogui" }
$argList += @("-execscript", "`"$scriptPath|$Function`"")

$p = Start-Process -FilePath $exe -ArgumentList $argList -PassThru -NoNewWindow `
    -RedirectStandardOutput $log -RedirectStandardError $err
if (-not $p.WaitForExit($TimeoutMin * 60 * 1000)) {
    Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
    Write-Output "TIMEOUT after $TimeoutMin min; log: $log"
    Get-Content $log -Tail 25
    exit 1
}
Write-Output "exit=$($p.ExitCode) log=$log"
Get-Content $log -Tail 25
Get-Content $err -Tail 10
