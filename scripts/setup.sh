#!/usr/bin/env bash
# 金策 BizGold 环境准备（Linux / macOS）
set -e
cd "$(dirname "$0")/.."
command -v node >/dev/null || { echo "缺少 node（20+）"; exit 1; }
command -v pnpm >/dev/null || { echo "缺少 pnpm（npm i -g pnpm）"; exit 1; }
PY=${BIZGOLD_PYTHON:-python3}
command -v "$PY" >/dev/null || { echo "缺少 python3"; exit 1; }
pnpm install
"$PY" -m pip install -r agent/toolkit/requirements.txt
PYTHONPATH=agent/toolkit "$PY" -m bizkit.fonts
echo "完成。pnpm dev → http://localhost:5173 ；或 pnpm build && pnpm start → http://localhost:5178"
