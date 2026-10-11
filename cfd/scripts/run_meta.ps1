param(
    [Parameter(Mandatory = $true)][string]$Session,
    [int]$TimeoutMin = 60,
    # batch (-b) mode on this machine draws only outlines; GUI mode renders filled contours. Session must end with "options quit".
    [switch]$Gui
)

$meta = "C:\Users\aalsayid\AppData\Local\Apps\BETA_CAE_Systems\meta_post_v24.1.6\meta_post64.bat"
$ses = (Resolve-Path $Session).Path
$logDir = Join-Path (Split-Path $PSScriptRoot) "logs"
New-Item -ItemType Directory -Force $logDir | Out-Null
$name = [IO.Path]::GetFileNameWithoutExtension($ses)
$log = Join-Path $logDir "meta_$name.log"
$err = Join-Path $logDir "meta_$name.err"
$work = Split-Path $ses

$args_ = @("-nolauncher", "-noses")
if (-not $Gui) { $args_ += "-b" } else { $args_ += @("-geom", "1600x1100") }
$args_ += @("-s", "`"$ses`"")
$p = Start-Process -FilePath $meta -ArgumentList $args_ `
    -WorkingDirectory $work -PassThru -NoNewWindow -RedirectStandardOutput $log -RedirectStandardError $err
if (-not $p.WaitForExit($TimeoutMin * 60 * 1000)) {
    Get-CimInstance Win32_Process | Where-Object { $_.Name -match 'meta_post' -and $_.CommandLine -match [regex]::Escape($ses) } |
        ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
    Stop-Process -Id $p.Id -Force -ErrorAction SilentlyContinue
    Write-Output "TIMEOUT after $TimeoutMin min; log: $log"
    Get-Content $log -Tail 30
    exit 1
}
Write-Output "exit=$($p.ExitCode) log=$log"
Get-Content $log -Tail 30
Get-Content $err -Tail 10
