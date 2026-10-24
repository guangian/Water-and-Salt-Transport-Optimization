# 生产模式：构建前端并由后端单端口托管（http://localhost:5178）
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
$env:PYTHONUTF8 = '1'
if (-not (Test-Path "$root\apps\web\dist\index.html")) { pnpm build }
pnpm start
