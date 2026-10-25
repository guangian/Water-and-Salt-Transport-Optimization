# 端到端实跑（Windows）：创建示例项目（马大姐电解质软糖 · 品牌策划）→ 上传资料 → 启动运行 → 轮询进度直到完成。
# 用法：powershell -ExecutionPolicy Bypass -File scripts\e2e_run.ps1 [-BaseUrl http://127.0.0.1:5178/api] [-ReviewRounds 2] [-PageTarget 55] [-TimeoutMin 240] [-Smoke]
param(
  [string]$BaseUrl = 'http://127.0.0.1:5178/api',
  [int]$ReviewRounds = 2,
  [int]$PageTarget = 55,
  [int]$TimeoutMin = 240,
  [switch]$Smoke
)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $PSScriptRoot
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8

Write-Output "== 健康检查 =="
$h = Invoke-RestMethod -Uri "$BaseUrl/health" -TimeoutSec 60
Write-Output ("ok=$($h.ok) python=$($h.python.version) pdf=$($h.pdfConverter) gitBash=$($h.gitBash) key=$($h.apiKeySource)")

$brief = @"
河北康贝尔食品有限公司是“马大姐”品牌所有者，国内老牌糖果企业（1996 年创立，酥糖、牛轧糖等传统品类知名度高，渠道以商超、流通批发为主）。2025 年推出“电解质软糖”新品（含钠、钾、镁等电解质，面向运动补水、夏季户外、加班熬夜等场景），希望借此实现品牌年轻化、进入功能性零食赛道，并打开线上与新零售渠道。命题要求：围绕电解质软糖进行品牌策划，覆盖全国市场，重点解决三个问题：① 老牌糖果品牌如何在 Z 世代人群中建立“功能性零食”新认知；② 电解质软糖相对电解质饮料、运动饮料与竞品软糖（如 BUFFX、Minayo 等）的差异化定位与价格带；③ 以有限预算（300 万元级）在 12 个月内完成从 0 到 1 的声量与销量目标。企业现有资产：全国商超与流通渠道、保定生产基地、“马大姐”品牌知名度（在 35 岁以上人群中较高）、少量抖音 / 天猫旗舰店运营经验。约束：不可过多陈述企业现有成就；需给出品牌元素（名称延展、视觉符号、口号、包装）、传播执行、预算与 KPI、品牌资产保护与风险管理。
"@
$body = @{
  name = '电解一夏，糖力觉醒——马大姐电解质软糖品牌策划书'
  teamName = '糖力实验室品牌策划工作室'
  track = 'brand'; group = '本科组'
  subject = @{ type = 'brand'; name = '马大姐电解质软糖'; company = '河北康贝尔食品有限公司'; isAssignedCase = $true
               caseTitle = '2026 届品牌策划赛道命题：马大姐电解质软糖项目策划（目标市场：全国）'; brief = $brief; region = '河北保定'; industry = '休闲食品 / 功能性糖果' }
  team = @(
    @{ role = '总策划 / 队长'; major = '市场营销'; strength = '品牌战略与调研' },
    @{ role = '创意与设计'; major = '视觉传达'; strength = '品牌视觉与短视频' },
    @{ role = '传播与渠道'; major = '新媒体运营'; strength = '抖音 / 小红书投放' },
    @{ role = '财务与数据'; major = '会计学'; strength = '预算与数据分析' })
  params = @{ researchMode = 'synthetic'; pageTarget = $PageTarget; pitchMinutes = 10; reviewRounds = $ReviewRounds; targetScore = 85; language = 'zh-CN' }
} | ConvertTo-Json -Depth 6

Write-Output "== 创建项目 =="
$p = Invoke-RestMethod -Uri "$BaseUrl/projects" -Method Post -ContentType 'application/json; charset=utf-8' -Body ([System.Text.Encoding]::UTF8.GetBytes($body)) -TimeoutSec 60
$projectId = $p.id
Write-Output "项目 id: $projectId"

Write-Output "== 上传资料 =="
$file = Join-Path $root 'docs\examples\马大姐_企业资料.txt'
& curl.exe -sS -X POST "$BaseUrl/projects/$projectId/files" -F "file=@$file" | Out-String | Write-Output

Write-Output "== 启动运行（smoke=$Smoke）=="
$runBody = if ($Smoke) { '{"mode":"run","smoke":true}' } else { '{"mode":"run"}' }
Invoke-RestMethod -Uri "$BaseUrl/projects/$projectId/run" -Method Post -ContentType 'application/json' -Body $runBody -TimeoutSec 60 | ConvertTo-Json -Compress | Write-Output

$start = Get-Date
while ($true) {
  Start-Sleep -Seconds 60
  $pr = Invoke-RestMethod -Uri "$BaseUrl/projects/$projectId" -TimeoutSec 60
  $stages = ($pr.stages.PSObject.Properties | ForEach-Object { "$($_.Name):$($_.Value.status.Substring(0, [Math]::Min(4, $_.Value.status.Length)))" }) -join '|'
  $u = $pr.usage
  $elapsed = [int]((Get-Date) - $start).TotalMinutes
  Write-Output "[$elapsed min] $($pr.status) | $stages | turns=$($u.turns),in=$($u.inputTokens),out=$($u.outputTokens),score=$($pr.review.latestScore)"
  if ($pr.status -in @('completed', 'failed', 'stopped')) { break }
  if ($elapsed -ge $TimeoutMin) { Write-Output "已超过 $TimeoutMin 分钟，脚本停止轮询（运行不会被中断，可在工作台继续查看或手动停止）"; break }
}
Write-Output "== 最终状态 =="
$pr = Invoke-RestMethod -Uri "$BaseUrl/projects/$projectId" -TimeoutSec 60
$pr | Select-Object status, review, artifacts | ConvertTo-Json -Depth 6
Write-Output "工作区：$($h.workspaceDir)\$projectId"
