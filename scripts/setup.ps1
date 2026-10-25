# 金策 BizGold 一键环境准备（Windows PowerShell）
# 用法：在项目根目录执行  powershell -ExecutionPolicy Bypass -File scripts\setup.ps1
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
Set-Location $root
Write-Host "== 金策 BizGold 环境准备 ==" -ForegroundColor Cyan

function Need($cmd, $hint) {
  $c = Get-Command $cmd -ErrorAction SilentlyContinue
  if (-not $c) { Write-Host "缺少 $cmd ：$hint" -ForegroundColor Red; return $false }
  Write-Host "✓ $cmd => $($c.Source)"; return $true
}
$ok = $true
$ok = (Need node "请安装 Node.js 20+ https://nodejs.org") -and $ok
$ok = (Need pnpm "npm i -g pnpm 或 corepack enable") -and $ok
$ok = (Need python "请安装 Python 3.10+（Anaconda 亦可）") -and $ok
if (-not (Test-Path 'C:\Program Files\Git\bin\bash.exe')) { Write-Host "提示：未找到 Git Bash（C:\Program Files\Git\bin\bash.exe），智能体将改用 PowerShell 工具；建议安装 Git for Windows。" -ForegroundColor Yellow } else { Write-Host "✓ Git Bash" }
if (-not $ok) { exit 1 }

Write-Host "`n[1/3] pnpm install（含 Claude Agent SDK 平台二进制，首次约 300 MB）..." -ForegroundColor Cyan
pnpm install
if ($LASTEXITCODE -ne 0) { Write-Host "pnpm install 失败，可尝试：pnpm install --registry https://registry.npmmirror.com" -ForegroundColor Red; exit 1 }

Write-Host "`n[2/3] Python 依赖..." -ForegroundColor Cyan
$pip = "-i https://pypi.tuna.tsinghua.edu.cn/simple"
python -m pip install -r agent\toolkit\requirements.txt $pip.Split(' ')
if ($LASTEXITCODE -ne 0) { Write-Host "pip 安装失败，请检查网络或手动执行 python -m pip install -r agent\toolkit\requirements.txt" -ForegroundColor Red; exit 1 }

Write-Host "`n[3/3] 工具包自检..." -ForegroundColor Cyan
$env:PYTHONPATH = "$root\agent\toolkit"; $env:PYTHONUTF8 = '1'
python -m bizkit.fonts | Out-String | Write-Host

Write-Host "`n完成。下一步：" -ForegroundColor Green
Write-Host "  1) 设置密钥（任选其一）：在系统环境变量设置 BIZGOLD_API_KEY / BIZGOLD_API_BASE_URL，或启动后在界面『设置』页填写（保存在 %USERPROFILE%\.bizgold\settings.json）。"
Write-Host "  2) 开发模式：pnpm dev   → 浏览器打开 http://localhost:5173"
Write-Host "  3) 生产模式：pnpm build; pnpm start → http://localhost:5178"
