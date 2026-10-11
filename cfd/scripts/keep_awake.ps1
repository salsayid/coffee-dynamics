param([int]$Hours = 72)

# Holds an ES_SYSTEM_REQUIRED request so Windows does not idle-sleep; released when this process exits.
Add-Type -Namespace Win32 -Name Power -MemberDefinition @"
[DllImport("kernel32.dll")] public static extern uint SetThreadExecutionState(uint esFlags);
"@
$ES_CONTINUOUS = [uint32]"0x80000000"
$ES_SYSTEM_REQUIRED = [uint32]"0x00000001"
$end = (Get-Date).AddHours($Hours)
while ((Get-Date) -lt $end) {
    [void][Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS -bor $ES_SYSTEM_REQUIRED)
    Start-Sleep -Seconds 60
}
[void][Win32.Power]::SetThreadExecutionState($ES_CONTINUOUS)
