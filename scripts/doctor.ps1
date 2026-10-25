# 一键体检（Windows）：检查运行 BizGold 所需的一切，缺什么直接给修复命令。
# 用法：powershell -ExecutionPolicy Bypass -File scripts\doctor.ps1 [-BaseUrl https://llm.goaichat.top]
param([string]$BaseUrl = $env:BIZGOLD_API_BASE_URL)
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$root = Split-Path -Parent $PSScriptRoot
$fail = 0
function Ok($m) { Write-Host "  [OK]   $m" -ForegroundColor Green }
function Warn($m) { Write-Host "  [WARN] $m" -ForegroundColor Yellow }
function Bad($m, $fix) { Write-Host "  [FAIL] $m" -ForegroundColor Red; if ($fix) { Write-Host "         修复：$fix" -ForegroundColor Gray }; $script:fail++ }
function Ver($cmd, $arg) { try { $o = & $cmd $arg 2>$null | Where-Object { $_ -match '\d' } | Select-Object -First 1; if ($o) { $o.ToString().Trim() } } catch { $null } }
function Http($url) { if (Get-Command curl.exe -ErrorAction SilentlyContinue) { $o = & curl.exe -sS --noproxy '*' --max-time 8 $url 2>$null; if ($LASTEXITCODE -eq 0 -and $o) { return ($o | ConvertFrom-Json) } }; return (Invoke-RestMethod -Uri $url -TimeoutSec 8) }

Write-Host "== 运行环境 =="
$node = Ver node '--version'
if ($node -and [int]($node.TrimStart('v').Split('.')[0]) -ge 20) { Ok "Node $node" } else { Bad "Node 20+ 未安装（当前：$node）" "https://nodejs.org 安装 LTS 版" }
$pnpm = Ver pnpm '--version'
if ($pnpm) { Ok "pnpm $($pnpm -replace 'Version\s*','')" } else { Bad "pnpm 未安装" "npm i -g pnpm" }
$py = Ver python '--version'
if ($py -match '3\.(1[0-9]|[2-9][0-9])') { Ok "$py" } else { Bad "Python 3.10+ 未找到（当前：$py）" "https://www.python.org/downloads/ 安装并勾选 Add to PATH；或设置环境变量 BIZGOLD_PYTHON 指向解释器" }
$gitBash = @("$env:ProgramFiles\Git\bin\bash.exe", "${env:ProgramFiles(x86)}\Git\bin\bash.exe", "$env:LOCALAPPDATA\Programs\Git\bin\bash.exe") | Where-Object { Test-Path $_ } | Select-Object -First 1
if ($gitBash) { Ok "Git Bash：$gitBash" } else { Bad "Git for Windows 未安装（智能体的 Bash 工具依赖它）" "https://git-scm.com/download/win" }

Write-Host "== Office / PDF 转换 =="
$word = Get-ItemProperty 'HKLM:\SOFTWARE\Classes\Word.Application\CurVer' -ErrorAction SilentlyContinue
$wps = Get-Command wps -ErrorAction SilentlyContinue
$soffice = Get-Command soffice -ErrorAction SilentlyContinue
if ($word) { Ok "Microsoft Word 可用（docx → PDF 走 Word）" } elseif ($wps) { Ok "WPS 可用" } elseif ($soffice) { Ok "LibreOffice 可用" } else { Bad "没有 Word / WPS / LibreOffice，无法生成 PDF" "安装任一即可；LibreOffice 免费：https://www.libreoffice.org/" }
Add-Type -AssemblyName System.Drawing -ErrorAction SilentlyContinue
$fonts = try { (New-Object System.Drawing.Text.InstalledFontCollection).Families.Name } catch { @() }
if ($fonts -and ($fonts -contains '微软雅黑' -or $fonts -contains 'Microsoft YaHei' -or $fonts -contains '思源黑体 CN' -or $fonts -contains 'Noto Sans CJK SC')) { Ok "中文字体就绪" } else { Warn "未检测到微软雅黑 / 思源黑体，图表与 PPT 中文可能显示为方块" }

Write-Host "== 项目与依赖 =="
if (Test-Path "$root\node_modules\.pnpm") { Ok "node_modules 已安装" } else { Bad "依赖未安装" "pnpm install" }
if ((Test-Path "$root\apps\web\dist\index.html") -and (Test-Path "$root\apps\server\dist\index.js")) { Ok "前后端已构建" } else { Bad "尚未构建" "pnpm build" }
$pyMods = & python -c "import docx, pptx, openpyxl, matplotlib, pypdf, pdfplumber; print('ok')" 2>$null
if ($pyMods -eq 'ok') { Ok "Python 工具包依赖齐全" } else { Bad "Python 依赖缺失" "python -m pip install -r agent\toolkit\requirements.txt" }

Write-Host "== 密钥与网关 =="
$settings = Join-Path $env:USERPROFILE '.bizgold\settings.json'
$keyFrom = if ($env:BIZGOLD_API_KEY) { '环境变量' } elseif ((Test-Path $settings) -and ((Get-Content $settings -Raw) -match '"apiKey"\s*:\s*"sk-')) { "设置文件 $settings" } else { $null }
if ($keyFrom) { Ok "API 密钥来源：$keyFrom" } else { Bad "未配置 API 密钥" "设置环境变量 BIZGOLD_API_BASE_URL / BIZGOLD_API_KEY，或启动后在「设置」页填写保存" }
if (-not $BaseUrl -and (Test-Path $settings)) { try { $BaseUrl = (Get-Content $settings -Raw | ConvertFrom-Json).api.baseUrl } catch {} }
if ($BaseUrl) {
  try { $r = Invoke-WebRequest -Uri "$BaseUrl/v1/models" -TimeoutSec 15 -UseBasicParsing -ErrorAction Stop; Ok "网关可达：$BaseUrl（HTTP $($r.StatusCode)）" }
  catch { if ($_.Exception.Response -and [int]$_.Exception.Response.StatusCode -in 401,403) { Ok "网关可达：$BaseUrl（需要密钥，正常）" } else { Bad "网关不可达：$BaseUrl（$($_.Exception.Message)）" "检查网络 / 代理，或在设置里切换备用地址" } }
} else { Warn "未知网关地址，跳过连通性测试" }

Write-Host "== 服务 =="
try { $h = Http 'http://127.0.0.1:5178/api/health'; Ok "服务端在线：密钥=$($h.apiKeySource) PDF=$($h.pdfConverter) 运行中=$($h.activeRuns)" }
catch { Warn "服务端未运行（pnpm start，或 scripts\service.ps1 -Install 注册常驻）" }

Write-Host ""
if ($fail -eq 0) { Write-Host "体检通过：可以打开 http://localhost:5178 开始使用。" -ForegroundColor Green } else { Write-Host "有 $fail 项需要处理（见上方 FAIL 行）。" -ForegroundColor Red; exit 1 }
