#!/usr/bin/env bash
# 一键体检（Linux / macOS）：检查运行 BizGold 所需的一切，缺什么直接给修复命令。
# 用法：bash scripts/doctor.sh [网关地址]
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BASE_URL="${1:-${BIZGOLD_API_BASE_URL:-}}"
FAIL=0
ok() { printf '  \033[32m[OK]\033[0m   %s\n' "$1"; }
warn() { printf '  \033[33m[WARN]\033[0m %s\n' "$1"; }
bad() { printf '  \033[31m[FAIL]\033[0m %s\n' "$1"; [ -n "${2:-}" ] && printf '         修复：%s\n' "$2"; FAIL=$((FAIL + 1)); }

echo "== 运行环境 =="
NODE_V="$(node --version 2>/dev/null || true)"
if [ -n "$NODE_V" ] && [ "${NODE_V#v}" != "$NODE_V" ] && [ "$(echo "${NODE_V#v}" | cut -d. -f1)" -ge 20 ]; then ok "Node $NODE_V"; else bad "Node 20+ 未安装（当前：${NODE_V:-无}）" "https://nodejs.org 或 nvm install 20"; fi
if command -v pnpm >/dev/null 2>&1; then ok "pnpm $(pnpm --version)"; else bad "pnpm 未安装" "npm i -g pnpm"; fi
PY="${BIZGOLD_PYTHON:-python3}"
PY_V="$($PY --version 2>/dev/null || true)"
if echo "$PY_V" | grep -Eq '3\.(1[0-9]|[2-9][0-9])'; then ok "$PY_V（$PY）"; else bad "Python 3.10+ 未找到（当前：${PY_V:-无}）" "安装 python3，或设置 BIZGOLD_PYTHON 指向解释器"; fi

echo "== Office / PDF 转换 =="
if command -v soffice >/dev/null 2>&1 || command -v libreoffice >/dev/null 2>&1; then ok "LibreOffice 可用（docx → PDF）"; else bad "没有 LibreOffice，无法生成 PDF" "Ubuntu: sudo apt install -y libreoffice-writer libreoffice-impress；macOS: brew install --cask libreoffice"; fi
if command -v fc-list >/dev/null 2>&1; then
  if fc-list :lang=zh family 2>/dev/null | grep -Eqi 'Noto Sans CJK|Source Han|WenQuanYi|PingFang|Hiragino|Microsoft YaHei|SimHei'; then ok "中文字体就绪"; else bad "缺中文字体，图表 / PPT 会显示方块" "Ubuntu: sudo apt install -y fonts-noto-cjk fonts-wqy-microhei"; fi
else warn "无 fc-list，跳过字体检查（macOS 自带苹方即可）"; fi

echo "== 项目与依赖 =="
[ -d "$ROOT/node_modules/.pnpm" ] && ok "node_modules 已安装" || bad "依赖未安装" "pnpm install"
if [ -f "$ROOT/apps/web/dist/index.html" ] && [ -f "$ROOT/apps/server/dist/index.js" ]; then ok "前后端已构建"; else bad "尚未构建" "pnpm build"; fi
if "$PY" -c 'import docx, pptx, openpyxl, matplotlib, pypdf, pdfplumber' 2>/dev/null; then ok "Python 工具包依赖齐全"; else bad "Python 依赖缺失" "$PY -m pip install -r agent/toolkit/requirements.txt"; fi

echo "== 密钥与网关 =="
SETTINGS="$HOME/.bizgold/settings.json"
if [ -n "${BIZGOLD_API_KEY:-}" ]; then ok "API 密钥来源：环境变量"
elif [ -f "$SETTINGS" ] && grep -q '"apiKey"[[:space:]]*:[[:space:]]*"sk-' "$SETTINGS"; then ok "API 密钥来源：$SETTINGS"
else bad "未配置 API 密钥" "export BIZGOLD_API_BASE_URL=… BIZGOLD_API_KEY=…，或启动后在「设置」页填写保存"; fi
if [ -z "$BASE_URL" ] && [ -f "$SETTINGS" ]; then BASE_URL="$("$PY" -c 'import json,sys;print(json.load(open(sys.argv[1])).get("api",{}).get("baseUrl",""))' "$SETTINGS" 2>/dev/null || true)"; fi
if [ -n "$BASE_URL" ]; then
  CODE="$(curl -s -o /dev/null -w '%{http_code}' --max-time 15 "$BASE_URL/v1/models" || echo 000)"
  case "$CODE" in 200|401|403) ok "网关可达：$BASE_URL（HTTP $CODE）";; *) bad "网关不可达：$BASE_URL（HTTP $CODE）" "检查网络 / 代理，或在设置里切换备用地址";; esac
else warn "未知网关地址，跳过连通性测试"; fi

echo "== 服务 =="
if H="$(curl -s --max-time 5 http://127.0.0.1:5178/api/health)"; then ok "服务端在线：$(echo "$H" | "$PY" -c 'import json,sys;d=json.load(sys.stdin);print("密钥=%s PDF=%s 运行中=%s"%(d.get("apiKeySource"),d.get("pdfConverter"),d.get("activeRuns")))' 2>/dev/null)"
else warn "服务端未运行（pnpm start）"; fi

echo
if [ "$FAIL" -eq 0 ]; then printf '\033[32m体检通过：可以打开 http://localhost:5178 开始使用。\033[0m\n'; else printf '\033[31m有 %s 项需要处理（见上方 FAIL 行）。\033[0m\n' "$FAIL"; exit 1; fi
