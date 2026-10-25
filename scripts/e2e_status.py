#!/usr/bin/env python3
"""读取 stdin 的 project.json，打印一行状态摘要（供 e2e_run.sh 轮询）。"""
import json, sys
p = json.load(sys.stdin)
st = "|".join(f"{k}:{v.get('status', '-')[:4]}" + (f"(r{v.get('round')})" if v.get('round') else '') for k, v in p.get("stages", {}).items())
u = p.get("usage", {})
rv = p.get("review", {}) or {}
print(p.get("status", "?"), st, f"turns={u.get('turns', 0)},in={u.get('inputTokens', 0)},out={u.get('outputTokens', 0)},score={rv.get('latestScore')}")
