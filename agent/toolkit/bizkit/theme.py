"""金策 BizGold 视觉体系：配色、字号、尺寸。所有渲染模块（图表 / Word / PPT / Excel）共用。

设计语言：海军蓝为主色（专业、可信），金色为强调色（“国金”），灰蓝做辅助；
图表序列色 6 个，亮度递进、黑白打印可分；语义色用于正/负/警示。
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

# ---- 颜色（十六进制，不带 #）----
NAVY = "14305C"        # 主色：深海军蓝
NAVY_DARK = "0B1F3F"
GOLD = "C9A227"        # 强调：金
GOLD_LIGHT = "E8D48B"
SLATE = "5B6B82"       # 辅助灰蓝
SLATE_LIGHT = "C7D0DD"
INK = "1F2937"         # 正文深灰
MUTED = "6B7280"       # 次要文字
LINE = "D9DEE7"        # 分隔线 / 表格线
PAPER = "FFFFFF"
PAPER_ALT = "F5F7FA"   # 斑马纹 / 卡片底
POS = "2E7D5B"         # 正向 / 增长
NEG = "B23B3B"         # 负向 / 风险
WARN = "D98E04"

# 图表序列色：主蓝 → 金 → 青 → 灰蓝 → 砖红 → 橄榄（同屏定性色建议 ≤4）
SERIES: List[str] = ["14305C", "C9A227", "2A9D8F", "5B6B82", "B5533C", "7A8B3A"]
SERIES_LIGHT: List[str] = ["4A6FA5", "E0C46B", "7FC8BE", "9AA7B9", "D98C7A", "AFC06E"]
# 顺序色（亮度单调）
SEQUENTIAL: List[str] = ["E8EEF6", "C3D2E6", "9DB4D6", "6F91BF", "456EA3", "14305C"]
# 发散色（对零对称）
DIVERGING: List[str] = ["B23B3B", "D98C7A", "F3E3DF", "F5F7FA", "C3D2E6", "6F91BF", "14305C"]


def hex2rgb(h: str) -> Tuple[int, int, int]:
    h = h.lstrip("#")
    return int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)


def mpl(h: str) -> str:
    """matplotlib 颜色字符串。"""
    return "#" + h.lstrip("#")


@dataclass
class Typo:
    """字号体系（pt）。Word 用正文小四 = 12pt；PPT 用大字号。"""
    doc_body: float = 12.0        # 小四
    doc_small: float = 10.5       # 五号（图注 / 来源 / 表内）
    doc_h1: float = 16.0          # 三号
    doc_h2: float = 14.0          # 四号
    doc_h3: float = 12.0          # 小四加粗
    doc_title: float = 26.0       # 封面标题
    ppt_title: float = 32.0
    ppt_subtitle: float = 20.0
    ppt_body: float = 20.0
    ppt_small: float = 14.0
    ppt_kpi: float = 44.0
    chart_label: float = 9.0
    chart_title: float = 11.0
    chart_tick: float = 8.5
    chart_source: float = 7.5


@dataclass
class Theme:
    name: str = "bizgold"
    primary: str = NAVY
    primary_dark: str = NAVY_DARK
    accent: str = GOLD
    accent_light: str = GOLD_LIGHT
    secondary: str = SLATE
    secondary_light: str = SLATE_LIGHT
    ink: str = INK
    muted: str = MUTED
    line: str = LINE
    paper: str = PAPER
    paper_alt: str = PAPER_ALT
    positive: str = POS
    negative: str = NEG
    warning: str = WARN
    series: List[str] = field(default_factory=lambda: list(SERIES))
    series_light: List[str] = field(default_factory=lambda: list(SERIES_LIGHT))
    sequential: List[str] = field(default_factory=lambda: list(SEQUENTIAL))
    diverging: List[str] = field(default_factory=lambda: list(DIVERGING))
    typo: Typo = field(default_factory=Typo)

    # 图表默认尺寸（毫米）：A4 版心宽约 160mm；PPT 16:9 画布 338.7×190.5mm
    chart_w_mm: float = 160.0
    chart_h_mm: float = 90.0
    dpi: int = 300


THEMES = {"bizgold": Theme()}


def get_theme(name: str = "bizgold") -> Theme:
    return THEMES.get(name, THEMES["bizgold"])


def mm2in(mm: float) -> float:
    return mm / 25.4
