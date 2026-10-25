# 把 BizGold 服务端注册为 Windows 计划任务常驻运行（不依赖打开的终端窗口；登录即可用）。
# 用法：powershell -ExecutionPolicy Bypass -File scripts\service.ps1 -Install   # 注册并启动
#       powershell -ExecutionPolicy Bypass -File scripts\service.ps1 -Status    # 查看状态与健康
#       powershell -ExecutionPolicy Bypass -File scripts\service.ps1 -Restart   # 重启（更新代码并 pnpm build 后用）
#       powershell -ExecutionPolicy Bypass -File scripts\service.ps1 -Uninstall # 停止并删除任务
# 密钥请放在 %USERPROFILE%\.bizgold\settings.json（设置页保存即可）或系统环境变量，不要写进本脚本。
param(
  [switch]$Install,
  [switch]$Uninstall,
  [switch]$Restart,
  [switch]$Status,
  [string]$TaskName = 'BizGoldServer',
  [string]$LogFile = (Join-Path $env:LOCALAPPDATA 'BizGold\server.log')
)
$ErrorActionPreference = 'Stop'
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$root = Split-Path -Parent $PSScriptRoot

function Get-Health {
  try { $h = Invoke-RestMethod -Uri 'http://127.0.0.1:5178/api/health' -TimeoutSec 10; "健康：ok=$($h.ok) 密钥=$($h.apiKeySource) PDF=$($h.pdfConverter) 运行中=$($h.activeRuns)" }
  catch { "健康：未响应（$($_.Exception.Message)）" }
}

function Install-Task {
  if (-not (Test-Path "$root\apps\server\dist\index.js")) { throw '先执行 pnpm build（或 scripts\setup.ps1）再注册服务' }
  New-Item -ItemType Directory -Force -Path (Split-Path -Parent $LogFile) | Out-Null
  $cmd = "& { Set-Location '$root'; `$env:PYTHONUTF8='1'; node apps\server\dist\index.js *>> '$LogFile' }"
  $action = New-ScheduledTaskAction -Execute 'powershell.exe' -Argument "-NoProfile -WindowStyle Hidden -ExecutionPolicy Bypass -Command $cmd" -WorkingDirectory $root
  $trigger = New-ScheduledTaskTrigger -AtLogOn -User $env:USERNAME
  $settings = New-ScheduledTaskSettingsSet -ExecutionTimeLimit (New-TimeSpan -Days 30) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries -StartWhenAvailable -RestartCount 3 -RestartInterval (New-TimeSpan -Minutes 1)
  Register-ScheduledTask -TaskName $TaskName -Action $action -Trigger $trigger -Settings $settings -Force | Out-Null
  Start-ScheduledTask -TaskName $TaskName
  Start-Sleep -Seconds 6
  Write-Output "已注册并启动计划任务 $TaskName（日志：$LogFile）"
  Write-Output (Get-Health)
}

function Stop-Task {
  if (Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue) {
    Stop-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    # Stop-ScheduledTask 只结束 powershell 宿主，node 可能仍占着 5178 端口
    Get-NetTCPConnection -LocalPort 5178 -State Listen -ErrorAction SilentlyContinue | ForEach-Object { Stop-Process -Id $_.OwningProcess -Force -ErrorAction SilentlyContinue }
    Start-Sleep -Seconds 2
  }
}

if ($Install) { Install-Task }
elseif ($Restart) { Stop-Task; Start-ScheduledTask -TaskName $TaskName; Start-Sleep -Seconds 6; Write-Output (Get-Health) }
elseif ($Uninstall) { Stop-Task; Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false -ErrorAction SilentlyContinue; Write-Output "已删除计划任务 $TaskName" }
else {
  $t = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
  if ($t) { Write-Output "计划任务 $TaskName：$($t.State)" } else { Write-Output "计划任务 $TaskName 未注册（用 -Install 注册）" }
  Write-Output (Get-Health)
}
