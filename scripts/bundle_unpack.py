#!/usr/bin/env python3
"""离线分发解包：把 bundle_pack.py 生成的 bundle_part_*.txt + bundle_manifest.json 还原为项目目录。

用法（Windows PowerShell 示例）：
  python scripts\bundle_unpack.py --parts D:\tmp\bundle --target "D:\0、项目列表\16、商赛agent"
校验 sha256 后用 tar 解到 target（已存在的同名文件会被覆盖，其他文件保留）。
"""
from __future__ import annotations

import argparse
import base64
import hashlib
import io
import json
import tarfile
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parts", required=True, help="存放 bundle_part_*.txt 与 bundle_manifest.json 的目录")
    ap.add_argument("--target", required=True, help="解压目标目录")
    args = ap.parse_args()

    parts_dir = Path(args.parts)
    manifest = json.loads((parts_dir / "bundle_manifest.json").read_text(encoding="utf-8"))
    files = sorted(parts_dir.glob("bundle_part_*.txt"))
    if len(files) != manifest["parts"]:
        raise SystemExit(f"分片数量不符：期望 {manifest['parts']}，实际 {len(files)}")
    b64 = "".join(f.read_text(encoding="ascii").strip() for f in files)
    data = base64.b64decode(b64)
    sha = hashlib.sha256(data).hexdigest()
    if sha != manifest["sha256"]:
        raise SystemExit(f"校验失败：sha256 {sha} != {manifest['sha256']}")
    target = Path(args.target)
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        members = tar.getmembers()
        for m in members:
            # 防御：不允许越出目标目录
            if m.name.startswith("/") or ".." in Path(m.name).parts:
                raise SystemExit(f"非法路径：{m.name}")
        tar.extractall(target)
    print(json.dumps({"ok": True, "files": len(members), "target": str(target), "sha256": sha}, ensure_ascii=False))


if __name__ == "__main__":
    main()
