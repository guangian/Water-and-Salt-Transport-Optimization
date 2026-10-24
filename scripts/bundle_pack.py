#!/usr/bin/env python3
"""离线分发打包：把项目源码（不含 node_modules / workspace / 大体积原始资料）打成 tar.gz，
再按 base64 文本分片，便于通过只能传输文本的通道（粘贴、IM、远程写文件工具）搬运到另一台机器。

用法：
  python scripts/bundle_pack.py --root . --out /tmp/bizgold_bundle [--chunk-mb 1.5] [--include-raw-dir docs/research/raw]
解包：目标机运行 scripts/bundle_unpack.py（见该文件说明）。
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
import os
import tarfile
import time
from pathlib import Path

EXCLUDE_DIRS = {"node_modules", "workspace", ".git", "dist", "__pycache__", ".venv", ".cache", ".claude-home", ".pytest_cache"}
EXCLUDE_FILES = {".env", "A2_DONE.txt", "B_DONE.txt", "C_DONE.txt"}
EXCLUDE_SUFFIX = {".pyc", ".log"}


def should_skip(rel: Path, include_raw: bool) -> bool:
    parts = rel.parts
    if any(p in EXCLUDE_DIRS for p in parts):
        return True
    if rel.name in EXCLUDE_FILES or rel.suffix in EXCLUDE_SUFFIX or rel.name.startswith("tmp_"):
        return True
    if not include_raw and len(parts) >= 3 and parts[0] == "docs" and parts[1] == "research" and parts[2] == "raw":
        return True
    # 工具包示例工作区里的渲染产物可由测试重新生成，不随包分发
    if "examples" in parts and (("产出" in parts) or ("08_提交包" in parts) or rel.suffix.lower() in {".png", ".pdf", ".docx", ".pptx", ".xlsx", ".zip"}):
        return True
    return False


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chunk-mb", type=float, default=1.5)
    ap.add_argument("--include-raw", action="store_true", help="包含 docs/research/raw 原始资料（体积大）")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    buf = io.BytesIO()
    count = 0
    with tarfile.open(fileobj=buf, mode="w:gz", format=tarfile.PAX_FORMAT) as tar:
        for dirpath, dirnames, filenames in os.walk(root):
            dirnames[:] = [d for d in dirnames if d not in EXCLUDE_DIRS]
            for fn in sorted(filenames):
                p = Path(dirpath) / fn
                rel = p.relative_to(root)
                if should_skip(rel, args.include_raw):
                    continue
                tar.add(p, arcname=str(rel).replace(os.sep, "/"))
                count += 1
    data = buf.getvalue()
    sha = hashlib.sha256(data).hexdigest()
    b64 = base64.b64encode(data).decode("ascii")
    chunk = int(args.chunk_mb * 1024 * 1024)
    parts = [b64[i:i + chunk] for i in range(0, len(b64), chunk)]
    for f in out.glob("bundle_part_*.txt"):
        f.unlink()
    for i, part in enumerate(parts, 1):
        (out / f"bundle_part_{i:03d}.txt").write_text(part, encoding="ascii")
    manifest = {"createdAt": time.strftime("%Y-%m-%dT%H:%M:%S"), "files": count, "tarGzBytes": len(data), "sha256": sha, "parts": len(parts), "chunkChars": chunk}
    (out / "bundle_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    main()
