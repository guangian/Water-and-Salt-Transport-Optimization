#!/usr/bin/env bash
# 端到端实跑：通过 REST API 创建示例项目（马大姐电解质软糖 · 品牌策划）、上传资料、启动运行、跟踪进度直至完成。
# 用法：BASE=http://localhost:5178/api bash scripts/e2e_run.sh [--smoke] [--track brand] [--timeout-min 240]
set -euo pipefail
BASE=${BASE:-http://localhost:5178/api}
SMOKE=false; TRACK=brand; TIMEOUT_MIN=240; REVIEW_ROUNDS=${REVIEW_ROUNDS:-2}; PAGE_TARGET=${PAGE_TARGET:-55}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --smoke) SMOKE=true; shift;;
    --track) TRACK="$2"; shift 2;;
    --timeout-min) TIMEOUT_MIN="$2"; shift 2;;
    *) echo "未知参数 $1"; exit 1;;
  esac
done
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
LOG_DIR=${LOG_DIR:-/tmp/bizgold_e2e}; mkdir -p "$LOG_DIR"

echo "== 健康检查 =="; curl -sS "$BASE/health" | head -c 600; echo
echo "== 创建项目 =="
if [[ "$TRACK" == "brand" ]]; then
read -r -d '' BODY <<'JSON' || true
{
  "name": "电解一夏，糖力觉醒——马大姐电解质软糖品牌策划书",
  "teamName": "糖力实验室品牌策划工作室",
  "track": "brand",
  "group": "本科组",
  "subject": {
    "type": "brand", "name": "马大姐电解质软糖", "company": "河北康贝尔食品有限公司", "isAssignedCase": true,
    "caseTitle": "2026 届品牌策划赛道命题：马大姐电解质软糖项目策划（目标市场：全国）",
    "brief": "河北康贝尔食品有限公司是“马大姐”品牌所有者，国内老牌糖果企业（1996 年创立，酥糖、牛轧糖等传统品类知名度高，渠道以商超、流通批发为主）。2025 年推出“电解质软糖”新品（含钠、钾、镁等电解质，面向运动补水、夏季户外、加班熬夜等场景），希望借此实现品牌年轻化、进入功能性零食赛道，并打开线上与新零售渠道。命题要求：围绕电解质软糖进行品牌策划，覆盖全国市场，重点解决三个问题：① 老牌糖果品牌如何在 Z 世代人群中建立“功能性零食”新认知；② 电解质软糖相对电解质饮料、运动饮料与竞品软糖（如 BUFFX、Minayo 等）的差异化定位与价格带；③ 以有限预算（300 万元级）在 12 个月内完成从 0 到 1 的声量与销量目标。企业现有资产：全国商超与流通渠道、保定生产基地、“马大姐”品牌知名度（在 35 岁以上人群中较高）、少量抖音 / 天猫旗舰店运营经验。约束：不可过多陈述企业现有成就；需给出品牌元素（名称延展、视觉符号、口号、包装）、传播执行、预算与 KPI、品牌资产保护与风险管理。",
    "region": "河北保定", "industry": "休闲食品 / 功能性糖果"
  },
  "team": [
    { "role": "总策划 / 队长", "major": "市场营销", "strength": "品牌战略与调研" },
    { "role": "创意与设计", "major": "视觉传达", "strength": "品牌视觉与短视频" },
    { "role": "传播与渠道", "major": "新媒体运营", "strength": "抖音 / 小红书投放" },
    { "role": "财务与数据", "major": "会计学", "strength": "预算与数据分析" }
  ],
  "params": { "researchMode": "synthetic", "pageTarget": PAGE_TARGET_PLACEHOLDER, "pitchMinutes": 10, "reviewRounds": REVIEW_ROUNDS_PLACEHOLDER, "targetScore": 85, "language": "zh-CN" }
}
JSON
else
  echo "暂只内置 brand 示例"; exit 1
fi
BODY=${BODY/PAGE_TARGET_PLACEHOLDER/$PAGE_TARGET}; BODY=${BODY/REVIEW_ROUNDS_PLACEHOLDER/$REVIEW_ROUNDS}
PROJECT=$(curl -sS -X POST "$BASE/projects" -H 'content-type: application/json' -d "$BODY")
PID=$(echo "$PROJECT" | python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])')
echo "项目 id: $PID"
echo "== 上传资料 =="
curl -sS -X POST "$BASE/projects/$PID/files" -F "file=@$ROOT/docs/examples/马大姐_企业资料.txt" | head -c 400; echo
echo "== 启动运行（smoke=$SMOKE）=="
if $SMOKE; then RUNBODY='{"mode":"run","smoke":true}'; else RUNBODY='{"mode":"run"}'; fi
curl -sS -X POST "$BASE/projects/$PID/run" -H 'content-type: application/json' -d "$RUNBODY"; echo
# 事件流落盘
( curl -sS -N "$BASE/projects/$PID/events" > "$LOG_DIR/$PID.events.log" 2>&1 ) &
EV_PID=$!
START=$(date +%s)
while true; do
  sleep 60
  P=$(curl -sS "$BASE/projects/$PID")
  STATUS=$(echo "$P" | python3 -c 'import json,sys; p=json.load(sys.stdin); print(p["status"])')
  STAGES=$(echo "$P" | python3 -c 'import json,sys; p=json.load(sys.stdin); print(" ".join(f"{k}:{v.get(\"status\",\"-\")[:4]}" for k,v in p.get("stages",{}).items()))')
  USAGE=$(echo "$P" | python3 -c 'import json,sys; p=json.load(sys.stdin); u=p.get("usage",{}); print(f"turns={u.get(\"turns\",0)} in={u.get(\"inputTokens\",0)} out={u.get(\"outputTokens\",0)}")')
  ELAPSED=$(( ($(date +%s) - START) / 60 ))
  echo "[$ELAPSED min] $STATUS | $STAGES | $USAGE"
  if [[ "$STATUS" == "completed" || "$STATUS" == "failed" || "$STATUS" == "stopped" ]]; then break; fi
  if (( ELAPSED >= TIMEOUT_MIN )); then echo "超时，停止"; curl -sS -X POST "$BASE/projects/$PID/stop"; break; fi
done
kill $EV_PID 2>/dev/null || true
echo "== 最终状态 =="; curl -sS "$BASE/projects/$PID" | python3 -c 'import json,sys; p=json.load(sys.stdin); print(json.dumps({k:p.get(k) for k in ("status","stages","usage","review","artifacts")}, ensure_ascii=False, indent=1))'
echo "== 产出树 =="; curl -sS "$BASE/projects/$PID/artifacts" | python3 -c '
import json,sys
def walk(n,d=0):
    print("  "*d + n["name"] + ("/" if n["type"]=="dir" else f"  ({n.get(\"size\",0)})"))
    for c in n.get("children",[]) or []: walk(c,d+1)
t=json.load(sys.stdin); walk(t) if isinstance(t,dict) else [walk(x) for x in t]' | head -120
echo "事件日志：$LOG_DIR/$PID.events.log"
