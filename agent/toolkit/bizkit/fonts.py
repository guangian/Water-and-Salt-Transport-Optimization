"""中文字体发现与映射（Windows / Linux / macOS 通用）。

用途：
- matplotlib：`apply_matplotlib()` 一次性设置中文字体与负号显示；
- python-docx / python-pptx：`doc_fonts()` 返回正文/标题/西文字体名（文件里只写字体名，由阅读端解析）。

Windows 优先：正文宋体（SimSun）、标题黑体（SimHei）、界面微软雅黑（Microsoft YaHei）；
Linux 退化：Noto Sans/Serif CJK SC → 文泉驿微米黑 → Droid Sans Fallback；
macOS：Songti SC / Heiti SC / PingFang SC。
"""
from __future__ import annotations

import os
import platform
from functools import lru_cache
from typing import Dict, List, Optional

# 候选顺序：(展示名, matplotlib 可识别的 family 名)
_SERIF_CANDIDATES = ["SimSun", "Songti SC", "Noto Serif CJK SC", "Source Han Serif SC", "AR PL UMing CN", "WenQuanYi Micro Hei", "Droid Sans Fallback"]
_SANS_CANDIDATES = ["Microsoft YaHei", "SimHei", "PingFang SC", "Heiti SC", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei", "Droid Sans Fallback"]
_LATIN_CANDIDATES = ["Times New Roman", "Liberation Serif", "DejaVu Serif"]


@lru_cache(maxsize=1)
def installed_font_names() -> List[str]:
    """返回 matplotlib 能找到的全部字体 family 名（去重）。"""
    try:
        from matplotlib import font_manager
        names = {f.name for f in font_manager.fontManager.ttflist}
        # Windows 下额外扫描系统字体目录，matplotlib 缓存可能滞后
        if platform.system() == "Windows":
            windir = os.environ.get("WINDIR", r"C:\Windows")
            for p in font_manager.findSystemFonts(fontpaths=[os.path.join(windir, "Fonts")]):
                try:
                    names.add(font_manager.FontProperties(fname=p).get_name())
                except Exception:
                    pass
        return sorted(names)
    except Exception:
        return []


def _first_available(cands: List[str], fallback: str) -> str:
    names = set(installed_font_names())
    for c in cands:
        if c in names:
            return c
    return fallback


@lru_cache(maxsize=1)
def pick_fonts() -> Dict[str, str]:
    """选择可用字体：serif（正文）、sans（标题/图表）、latin（西文）。"""
    serif = _first_available(_SERIF_CANDIDATES, _first_available(_SANS_CANDIDATES, "sans-serif"))
    sans = _first_available(_SANS_CANDIDATES, serif)
    latin = _first_available(_LATIN_CANDIDATES, "serif")
    return {"serif": serif, "sans": sans, "latin": latin}


def apply_matplotlib(prefer: str = "sans") -> str:
    """设置 matplotlib 全局中文字体，返回所用 family 名。"""
    import matplotlib
    fonts = pick_fonts()
    family = fonts["sans"] if prefer == "sans" else fonts["serif"]
    chain = [family] + [f for f in (_SANS_CANDIDATES + _SERIF_CANDIDATES) if f != family and f in set(installed_font_names())]
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = chain + ["DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    return family


def doc_fonts(body: Optional[str] = None, heading: Optional[str] = None, latin: Optional[str] = None) -> Dict[str, str]:
    """返回写入 Office 文件的字体名。

    Office 文件内只记录字体名，由打开它的机器解析；因此 Windows 目标机默认写“宋体/黑体/Times New Roman”，
    即使当前生成机（Linux）没有安装这些字体也没关系。可通过参数或环境变量覆盖：
    BIZKIT_BODY_FONT / BIZKIT_HEADING_FONT / BIZKIT_LATIN_FONT。
    """
    return {
        "body": body or os.environ.get("BIZKIT_BODY_FONT", "宋体"),
        "heading": heading or os.environ.get("BIZKIT_HEADING_FONT", "黑体"),
        "ui": os.environ.get("BIZKIT_UI_FONT", "微软雅黑"),
        "latin": latin or os.environ.get("BIZKIT_LATIN_FONT", "Times New Roman"),
    }


def report() -> Dict[str, object]:
    """环境自检用：当前机器可用的中文字体与选择结果。"""
    names = installed_font_names()
    cjk = [n for n in names if any(k in n for k in ("SimSun", "SimHei", "YaHei", "Songti", "Heiti", "PingFang", "CJK", "Han", "WenQuanYi", "Droid", "宋体", "黑体", "雅黑"))]
    return {"system": platform.system(), "picked": pick_fonts(), "cjkInstalled": cjk, "docFonts": doc_fonts()}


if __name__ == "__main__":
    import json
    print(json.dumps(report(), ensure_ascii=False, indent=2))
