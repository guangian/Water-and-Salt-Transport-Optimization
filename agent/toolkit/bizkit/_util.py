"""公共小工具：UTF-8 标准输出、JSON 读写、路径解析、日志。所有 CLI 模块统一从这里引入。"""
from __future__ import annotations

import io
import json
import os
import sys
from pathlib import Path
from typing import Any, Union


def force_utf8_stdout() -> None:
    """Windows 控制台默认 GBK，打印中文会报 UnicodeEncodeError；统一改为 UTF-8。"""
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name)
        try:
            if hasattr(stream, "reconfigure"):
                stream.reconfigure(encoding="utf-8", errors="replace")
            else:
                setattr(sys, name, io.TextIOWrapper(stream.buffer, encoding="utf-8", errors="replace"))
        except Exception:
            pass


def read_json(path: Union[str, Path]) -> Any:
    with open(path, "r", encoding="utf-8-sig") as f:
        return json.load(f)


def write_json(path: Union[str, Path], data: Any) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def read_text(path: Union[str, Path]) -> str:
    with open(path, "r", encoding="utf-8-sig") as f:
        return f.read()


def write_text(path: Union[str, Path], text: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)


def resolve(base: Union[str, Path], p: Union[str, Path]) -> Path:
    """相对路径相对 base 解析；绝对路径原样返回。"""
    p = Path(p)
    return p if p.is_absolute() else (Path(base) / p).resolve()


def workspace_root(start: Union[str, Path, None] = None) -> Path:
    """向上查找含 `任务书.md` 或 `CLAUDE.md` 的目录作为工作区根；找不到用当前目录。"""
    cur = Path(start or os.getcwd()).resolve()
    for d in [cur] + list(cur.parents):
        if (d / "任务书.md").exists() or (d / "CLAUDE.md").exists():
            return d
    return cur


def find_workspace_root(start: Union[str, Path, None], probe: Union[str, Path, None] = None) -> Path:
    """定位工作区根，用于解析 charts.json / book.json 里“相对工作区根”的路径。

    与 `workspace_root` 的区别：① 额外把 `project.json` 当作根标记（示例工作区与单元测试没有任务书）；
    ② 可给一个 probe 相对路径，依次在候选根（标记目录 → start 的父目录 → start → 当前目录）里探测它是否存在，
    用“哪个根下能找到该文件”来消歧，避免 cwd 不是工作区时解析到错误位置。
    """
    start_dir = Path(start or os.getcwd()).resolve()
    if start_dir.is_file():
        start_dir = start_dir.parent
    marked: list[Path] = []
    for d in [start_dir] + list(start_dir.parents):
        if any((d / m).exists() for m in ("任务书.md", "CLAUDE.md", "project.json")):
            marked.append(d)
    candidates = marked + [start_dir.parent, start_dir, Path(os.getcwd()).resolve()]
    if probe is not None:
        for c in candidates:
            if (c / probe).exists():
                return c
    return candidates[0]


def log(msg: str) -> None:
    print(f"[bizkit] {msg}", flush=True)
