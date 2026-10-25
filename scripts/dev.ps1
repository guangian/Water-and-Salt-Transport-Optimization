# 开发模式：同时启动后端（5178）与前端（5173）
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$env:PYTHONUTF8 = '1'
pnpm dev
