"""charts.json → 300 dpi PNG 图表（matplotlib Agg），契约 §7.1。

用法：
    python -m bizkit.charts 04_图表/charts.json [--only fig_a,fig_b]

输入 charts.json：
    { "outputDir": "04_图表", "dpi": 300, "theme": "bizgold",
      "charts": [ { "id": "fig_x", "type": "bar", "title": "…", "data": {…},
                    "xLabel": "", "yLabel": "", "unit": "亿元", "valueLabels": true,
                    "source": "资料来源：…", "note": "…", "width": 160, "height": 90, "legend": "top" } ] }

支持 type：bar / grouped_bar / stacked_bar / hbar / line / area / pie_as_bar / donut / waterfall / funnel /
tornado / breakeven / cashflow_combo / scatter_matrix / positioning_map / gantt / radar / bmc。
每种 type 的 data 结构见本模块 `DATA_SPECS`（校验失败会给出中文提示）。一张图失败不影响其他图，最后汇总。

输出：outputDir/<id>.png 与 outputDir/index.json（id → 文件、标题、来源、类型、尺寸）。
标题不画在图内（由文档图注承担），资料来源画在图底部左下角。
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
import traceback
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patches  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

from . import fonts, theme as theme_mod  # noqa: E402
from ._util import find_workspace_root, force_utf8_stdout, log, read_json, write_json  # noqa: E402

Number = float

# 黑白打印可分：第 3 个序列起叠加填充图案 / 线型
HATCHES = ["", "", "//", "..", "xx", "\\\\"]
LINESTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1))]
MARKERS = ["o", "s", "^", "D", "v", "P"]
CIRCLED = "①②③④⑤⑥⑦⑧⑨⑩"

# 各 type 的 data 结构说明（校验报错时附带，便于模型自行修正）
DATA_SPECS: Dict[str, str] = {
    "bar": "{categories:[…], values:[…]} 或 {categories, series:[{name, values}]}（单序列）",
    "grouped_bar": "{categories:[…], series:[{name, values:[…]}, …]}（序列 ≤ 6）",
    "stacked_bar": "{categories:[…], series:[{name, values:[…]}, …]}（序列 ≤ 6）",
    "hbar": "{categories:[…], values:[…], sort?: 'desc'|'asc'|'none', highlight?:[类目名]}",
    "line": "{categories:[…], series:[{name, values:[…]}, …]}",
    "area": "{categories:[…], series:[{name, values:[…]}, …], stacked?: true}",
    "pie_as_bar": "{categories:[…], values:[…]}（自动按占比降序 + 累计占比线）",
    "donut": "{categories:[…≤4], values:[…], centerText?: '…'}（>4 类自动改为 pie_as_bar）",
    "waterfall": "{items:[{name, value, total?: true}, …]}（total 项从 0 画起，value 省略时取累计值）",
    "funnel": "{stages:[{name, value}, …]}（按顺序递减）",
    "tornado": "{base: 数值, items:[{name, low, high}, …]}",
    "breakeven": "{fixedCost, unitPrice, unitVarCost, maxUnits?, unitName?: '件'}",
    "cashflow_combo": "{categories:[…], values:[当期现金流…], cumulative?:[…]}",
    "scatter_matrix": "{xLabel, yLabel, points:[{name, x, y, size?, self?: true}], quadrants?:[右上, 左上, 左下, 右下], xMid?, yMid?}",
    "positioning_map": "{xAxis:{low, high}, yAxis:{low, high}, points:[{name, x, y, size?, self?}], xRange?:[a,b], yRange?:[a,b]}",
    "gantt": "{tasks:[{name, start:'2026-04', end:'2026-06', group?, done?}], milestones?:[{name, date}]}",
    "radar": "{axes:[…≥3], series:[{name, values:[…]}], max?: 5}",
    "bmc": "{keyPartners:[…], keyActivities, keyResources, valuePropositions, customerRelationships, channels, customerSegments, costStructure, revenueStreams}",
}

# 没给 width/height 时按图型取更合适的默认尺寸（毫米）
DEFAULT_SIZES: Dict[str, Tuple[float, float]] = {
    "radar": (120, 105),
    "donut": (120, 90),
    "bmc": (170, 115),
    "gantt": (170, 90),
    "funnel": (140, 90),
    "scatter_matrix": (150, 110),
    "positioning_map": (140, 120),
    "tornado": (160, 95),
}


class ChartDataError(ValueError):
    """data 结构不符合约定；消息面向模型/用户，中文可读。"""


@dataclass
class ChartSpec:
    id: str
    type: str
    title: str
    data: Dict[str, Any]
    x_label: str = ""
    y_label: str = ""
    unit: str = ""
    value_labels: bool = True
    source: str = ""
    note: str = ""
    width_mm: float = 160.0
    height_mm: float = 90.0
    legend: str = "top"
    decimals: Optional[int] = None
    raw: Dict[str, Any] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, d: Dict[str, Any], th: theme_mod.Theme) -> "ChartSpec":
        if not isinstance(d, dict):
            raise ChartDataError("图表条目必须是对象")
        cid = str(d.get("id") or "").strip()
        if not cid:
            raise ChartDataError("缺少 id")
        if not re.fullmatch(r"[\w\-\u4e00-\u9fff]+", cid):
            raise ChartDataError(f"id「{cid}」只能包含字母、数字、下划线、连字符")
        ctype = str(d.get("type") or "").strip()
        if ctype not in RENDERERS:
            raise ChartDataError(f"不支持的 type「{ctype}」，可选：{', '.join(RENDERERS)}")
        data = d.get("data")
        if not isinstance(data, dict):
            raise ChartDataError(f"data 必须是对象，结构：{DATA_SPECS[ctype]}")
        dw, dh = DEFAULT_SIZES.get(ctype, (th.chart_w_mm, th.chart_h_mm))
        w = _num(d.get("width"), dw, "width")
        h = _num(d.get("height"), dh, "height")
        if not (40 <= w <= 400 and 30 <= h <= 400):
            raise ChartDataError("width/height 应在 40–400 毫米之间")
        dec = d.get("decimals")
        if dec is not None and (not isinstance(dec, int) or dec < 0 or dec > 4):
            raise ChartDataError("decimals 应为 0–4 的整数")
        return cls(
            id=cid, type=ctype, title=str(d.get("title") or ""), data=data,
            x_label=str(d.get("xLabel") or ""), y_label=str(d.get("yLabel") or ""),
            unit=str(d.get("unit") or ""), value_labels=bool(d.get("valueLabels", True)),
            source=str(d.get("source") or ""), note=str(d.get("note") or ""),
            width_mm=w, height_mm=h, legend=str(d.get("legend") or "top"), decimals=dec, raw=d,
        )


# ----------------------------------------------------------------------------------------------
# 基础工具
# ----------------------------------------------------------------------------------------------

def _num(v: Any, default: float, name: str) -> float:
    if v is None or v == "":
        return default
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        raise ChartDataError(f"{name} 必须是数字")
    return float(v)


def _require(data: Dict[str, Any], key: str, ctype: str) -> Any:
    if key not in data or data[key] is None:
        raise ChartDataError(f"data 缺少「{key}」。结构：{DATA_SPECS[ctype]}")
    return data[key]


def _as_floats(values: Any, name: str, allow_none: bool = False) -> List[float]:
    if not isinstance(values, (list, tuple)) or len(values) == 0:
        raise ChartDataError(f"「{name}」必须是非空数组")
    out: List[float] = []
    for i, v in enumerate(values):
        if v is None and allow_none:
            out.append(math.nan)
            continue
        if isinstance(v, bool) or not isinstance(v, (int, float)) or (isinstance(v, float) and math.isnan(v)):
            raise ChartDataError(f"「{name}」第 {i + 1} 个值「{v}」不是数字")
        out.append(float(v))
    return out


def _as_strs(values: Any, name: str) -> List[str]:
    if not isinstance(values, (list, tuple)) or len(values) == 0:
        raise ChartDataError(f"「{name}」必须是非空数组")
    return [str(v) for v in values]


def _series(data: Dict[str, Any], ctype: str, max_series: int = 6, single: bool = False) -> Tuple[List[str], List[Tuple[str, List[float]]]]:
    """统一解析 {categories, values} / {categories, series:[{name, values}]}。"""
    cats = _as_strs(_require(data, "categories", ctype), "categories")
    if "series" in data and data["series"] is not None:
        raw = data["series"]
        if not isinstance(raw, list) or not raw:
            raise ChartDataError(f"「series」必须是非空数组。结构：{DATA_SPECS[ctype]}")
        series: List[Tuple[str, List[float]]] = []
        for i, s in enumerate(raw):
            if not isinstance(s, dict):
                raise ChartDataError(f"series 第 {i + 1} 项必须是 {{name, values}} 对象")
            vals = _as_floats(s.get("values"), f"series[{i}].values", allow_none=True)
            if len(vals) != len(cats):
                raise ChartDataError(f"序列「{s.get('name', i + 1)}」有 {len(vals)} 个值，与 categories 的 {len(cats)} 个不一致")
            series.append((str(s.get("name") or f"系列{i + 1}"), vals))
    elif "values" in data:
        vals = _as_floats(data["values"], "values", allow_none=True)
        if len(vals) != len(cats):
            raise ChartDataError(f"values 有 {len(vals)} 个值，与 categories 的 {len(cats)} 个不一致")
        series = [(str(data.get("name") or ""), vals)]
    else:
        raise ChartDataError(f"data 需要「values」或「series」。结构：{DATA_SPECS[ctype]}")
    if single and len(series) > 1:
        raise ChartDataError(f"{ctype} 只支持单序列，多序列请改用 grouped_bar / stacked_bar / line")
    if len(series) > max_series:
        raise ChartDataError(f"序列数 {len(series)} 超过同屏 {max_series} 个定性色上限，请合并为「其他」或拆成两张图")
    return cats, series


def fmt_num(v: float, decimals: Optional[int] = None, unit: str = "") -> str:
    """千分位数字格式；decimals 为 None 时按量级自动取位数。"""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if decimals is None:
        av = abs(v)
        if float(v).is_integer():
            decimals = 0
        elif av >= 100:
            decimals = 0
        elif av >= 10:
            decimals = 1
        else:
            decimals = 2 if av < 1 else 1
    s = f"{v:,.{decimals}f}"
    if decimals > 0 and "." in s:
        s = s.rstrip("0").rstrip(".") if abs(v) < 100 else s
    if unit in ("%", "％"):
        s += "%"
    return s


def _auto_decimals(values: Sequence[float], override: Optional[int]) -> Optional[int]:
    if override is not None:
        return override
    vals = [v for v in values if v is not None and not math.isnan(v)]
    if not vals:
        return 0
    if all(float(v).is_integer() for v in vals):
        return 0
    return None


def _tick_formatter(unit: str = "") -> FuncFormatter:
    suffix = "%" if unit in ("%", "％") else ""

    def _f(v: float, _pos: int) -> str:
        if abs(v) >= 1000 or float(v).is_integer():
            return f"{v:,.0f}{suffix}"
        return f"{v:,.1f}".rstrip("0").rstrip(".") + suffix

    return FuncFormatter(_f)


_CJK_RE = re.compile(r"[\u3000-\u303f\u4e00-\u9fff\uff00-\uffef\u3400-\u4dbf]")


def _text_width(s: str) -> float:
    """以“中文字宽”为 1 的近似宽度：CJK 1.0，西文/数字 0.55，空格 0.3。"""
    return sum(1.0 if _CJK_RE.match(ch) else (0.3 if ch == " " else 0.55) for ch in s)


def _wrap(text: str, width: int) -> str:
    """中英文混排换行：CJK 字符间可断，西文/数字词内不断（避免 2026 → 20|26）。"""
    if width <= 0 or _text_width(text) <= width:
        return text
    tokens = re.findall(r"[\u3000-\u303f\u4e00-\u9fff\uff00-\uffef\u3400-\u4dbf]|\s+|[^\s\u3000-\u303f\u4e00-\u9fff\uff00-\uffef\u3400-\u4dbf]+", text)
    lines: List[str] = []
    cur = ""
    for tok in tokens:
        if cur and _text_width(cur + tok) > width and tok.strip():
            lines.append(cur.rstrip())
            cur = tok.lstrip()
        else:
            cur += tok
    if cur.strip():
        lines.append(cur.rstrip())
    return "\n".join(lines) or text


def _wrap_lines(text: str, width: int) -> List[str]:
    return _wrap(text, width).split("\n")


def _axis_title(label: str, unit: str) -> str:
    if label and unit and unit not in label:
        return f"{label}（{unit}）"
    if label:
        return label
    return f"（{unit}）" if unit else ""


def _mpl(h: str) -> str:
    return theme_mod.mpl(h)


class Canvas:
    """一张图的生命周期：建图 → 作图 → 页脚（注 / 来源）→ 保存。统一视觉规范在这里实现。"""

    def __init__(self, spec: ChartSpec, th: theme_mod.Theme, dpi: int):
        self.spec, self.th, self.dpi = spec, th, dpi
        self.typo = th.typo
        self.fig = plt.figure(figsize=(theme_mod.mm2in(spec.width_mm), theme_mod.mm2in(spec.height_mm)), dpi=dpi)
        self.fig.patch.set_facecolor("white")
        self.footer_lines = self._footer_lines()

    def _footer_lines(self) -> List[str]:
        lines: List[str] = []
        if self.spec.note:
            n = self.spec.note.strip()
            lines.append(n if n.startswith(("注", "备注")) else f"注：{n}")
        if self.spec.source:
            s = self.spec.source.strip()
            lines.append(s if s.startswith(("资料来源", "数据来源", "来源", "图片来源")) else f"资料来源：{s}")
        return lines

    def footer_frac(self) -> float:
        """页脚占画布高度的比例：每行约 3.4mm + 边距。"""
        if not self.footer_lines:
            return 0.0
        return min(0.35, (len(self.footer_lines) * 3.4 + 2.0) / self.spec.height_mm)

    def axes(self, **kw: Any):
        ax = self.fig.add_subplot(111, **kw)
        self.style_axes(ax)
        return ax

    def style_axes(self, ax, grid: str = "y") -> None:
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(_mpl(self.th.secondary_light))
            ax.spines[side].set_linewidth(0.8)
        ax.tick_params(labelsize=self.typo.chart_tick, colors=_mpl(self.th.ink), length=3, width=0.6)
        ax.set_axisbelow(True)
        if grid in ("y", "both"):
            ax.yaxis.grid(True, color=_mpl(self.th.line), linewidth=0.6)
        if grid in ("x", "both"):
            ax.xaxis.grid(True, color=_mpl(self.th.line), linewidth=0.6)

    def legend(self, ax, n_items: int, handles=None, labels=None) -> None:
        pos = (self.spec.legend or "top").lower()
        if pos == "none" or n_items <= 1 and handles is None:
            return
        kw: Dict[str, Any] = dict(fontsize=self.typo.chart_tick, frameon=False, handlelength=1.6, columnspacing=1.2, handletextpad=0.5)
        if handles is not None:
            kw["handles"] = handles
        if labels is not None:
            kw["labels"] = labels
        if pos == "top":
            ax.legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=min(max(n_items, 1), 4), borderaxespad=0.2, **kw)
        elif pos == "bottom":
            ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=min(max(n_items, 1), 4), **kw)
        elif pos == "right":
            ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), **kw)
        else:
            ax.legend(loc="best", **kw)

    def save(self, path: Path) -> Tuple[int, int]:
        bottom = self.footer_frac()
        try:
            self.fig.tight_layout(rect=(0, bottom, 1, 1), pad=0.6)
        except Exception:
            self.fig.subplots_adjust(bottom=bottom + 0.12)
        y = 0.012
        for line in reversed(self.footer_lines):
            self.fig.text(0.012, y, line, fontsize=self.typo.chart_source, color=_mpl(self.th.muted), ha="left", va="bottom")
            y += 3.4 / self.spec.height_mm
        self.fig.savefig(path, dpi=self.dpi, facecolor="white")
        w, h = self.fig.get_size_inches()
        return int(round(w * self.dpi)), int(round(h * self.dpi))

    def close(self) -> None:
        plt.close(self.fig)


def _cat_wrap_width(n_cats: int, width_mm: float, pt: float) -> int:
    """类目标签的换行宽度：轴宽 / 类目数 ÷ 字宽。"""
    per_cat_mm = max(8.0, (width_mm * 0.82) / max(n_cats, 1))
    return max(3, int(per_cat_mm / (pt * 0.3528) * 0.95))


def _color_cycle(th: theme_mod.Theme, n: int) -> List[str]:
    return [_mpl(c) for c in th.series[:n]]


def _label_series_end(ax, xs: Sequence[float], ys_list: List[List[float]], names: List[str], colors: List[str], fontsize: float) -> None:
    """在折线 / 面积末端直接标注序列名，并在显示坐标里错开避免重叠。"""
    ends = []
    for ys, name, col in zip(ys_list, names, colors):
        for x, y in zip(reversed(list(xs)), reversed(ys)):
            if not math.isnan(y):
                ends.append([x, y, name, col])
                break
    if not ends:
        return
    ax.figure.canvas.draw()
    trans = ax.transData
    inv = trans.inverted()
    ends.sort(key=lambda e: e[1])
    disp = [trans.transform((e[0], e[1]))[1] for e in ends]
    min_gap = fontsize * ax.figure.dpi / 72 * 1.15
    for i in range(1, len(disp)):
        if disp[i] - disp[i - 1] < min_gap:
            disp[i] = disp[i - 1] + min_gap
    for e, dy in zip(ends, disp):
        y_new = inv.transform((0, dy))[1]
        moved = abs(y_new - e[1]) > 1e-9
        ax.annotate(e[2], xy=(e[0], e[1]), xytext=(8, dy - trans.transform((e[0], e[1]))[1]),
                    textcoords="offset pixels", fontsize=fontsize, color=e[3], va="center", ha="left", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", color=e[3], lw=0.5, alpha=0.6) if moved else None)


# ----------------------------------------------------------------------------------------------
# 各图型渲染器（签名统一：(canvas, spec) -> None）
# ----------------------------------------------------------------------------------------------

def _bars(cv: Canvas, spec: ChartSpec, mode: str) -> None:
    th, typo = cv.th, cv.typo
    cats, series = _series(spec.data, spec.type, single=(mode == "single"))
    ax = cv.axes()
    n, m = len(cats), len(series)
    x = list(range(n))
    colors = _color_cycle(th, m)
    dec = _auto_decimals([v for _, vals in series for v in vals], spec.decimals)
    highlight = set(str(h) for h in (spec.data.get("highlight") or []))
    label_fs = max(8.0, typo.chart_label)

    if mode == "grouped":
        width = 0.8 / m
        for i, (name, vals) in enumerate(series):
            xs = [xi - 0.4 + width * (i + 0.5) for xi in x]
            bars = ax.bar(xs, [0 if math.isnan(v) else v for v in vals], width=width * 0.92, color=colors[i], label=name,
                          hatch=HATCHES[i] if m >= 3 else "", edgecolor="white", linewidth=0.6)
            if spec.value_labels and n * m <= 40:
                ax.bar_label(bars, labels=[fmt_num(v, dec, spec.unit) for v in vals], padding=2, fontsize=label_fs, color=_mpl(th.ink))
    elif mode == "stacked":
        bottoms = [0.0] * n
        for i, (name, vals) in enumerate(series):
            clean = [0 if math.isnan(v) else v for v in vals]
            bars = ax.bar(x, clean, width=0.62, bottom=bottoms, color=colors[i], label=name,
                          hatch=HATCHES[i] if m >= 3 else "", edgecolor="white", linewidth=0.6)
            if spec.value_labels and n * m <= 40:
                total_max = max(sum(s[1][j] for s in series if not math.isnan(s[1][j])) for j in range(n)) or 1
                labels = [fmt_num(v, dec, spec.unit) if v / total_max >= 0.07 else "" for v in clean]
                ax.bar_label(bars, labels=labels, label_type="center", fontsize=label_fs, color="white")
            bottoms = [b + v for b, v in zip(bottoms, clean)]
        if spec.value_labels:
            for xi, tot in zip(x, bottoms):
                ax.text(xi, tot, fmt_num(tot, dec, spec.unit), ha="center", va="bottom", fontsize=label_fs, color=_mpl(th.ink), fontweight="bold")
        ax.set_ylim(0, max(bottoms) * 1.12 if max(bottoms) > 0 else 1)
    else:
        name, vals = series[0]
        clean = [0 if math.isnan(v) else v for v in vals]
        bar_colors = [_mpl(th.accent) if c in highlight else (_mpl(th.negative) if v < 0 else colors[0]) for c, v in zip(cats, clean)]
        bars = ax.bar(x, clean, width=0.6, color=bar_colors, label=name or None)
        if spec.value_labels:
            ax.bar_label(bars, labels=[fmt_num(v, dec, spec.unit) for v in vals], padding=2, fontsize=label_fs, color=_mpl(th.ink))

    ax.set_xticks(x)
    ax.set_xticklabels([_wrap(c, _cat_wrap_width(n, spec.width_mm, typo.chart_tick)) for c in cats])
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_ylabel(_axis_title(spec.y_label, spec.unit), fontsize=typo.chart_label)
    if spec.x_label:
        ax.set_xlabel(spec.x_label, fontsize=typo.chart_label)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(min(0, ymin), ymax * 1.08 if ymax > 0 else ymax)
    ax.margins(x=0.02 if n > 8 else 0.05)
    if m > 1 or (series[0][0] and mode == "single"):
        cv.legend(ax, m)


def render_bar(cv: Canvas, spec: ChartSpec) -> None:
    _bars(cv, spec, "single")


def render_grouped_bar(cv: Canvas, spec: ChartSpec) -> None:
    _bars(cv, spec, "grouped")


def render_stacked_bar(cv: Canvas, spec: ChartSpec) -> None:
    _bars(cv, spec, "stacked")


def render_hbar(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    cats, series = _series(spec.data, "hbar", single=True)
    vals = [0 if math.isnan(v) else v for v in series[0][1]]
    sort = str(spec.data.get("sort") or "none").lower()
    pairs = list(zip(cats, vals))
    if sort == "desc":
        pairs.sort(key=lambda p: p[1], reverse=True)
    elif sort == "asc":
        pairs.sort(key=lambda p: p[1])
    cats = [p[0] for p in pairs]
    vals = [p[1] for p in pairs]
    highlight = set(str(h) for h in (spec.data.get("highlight") or []))
    ax = cv.axes()
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color=_mpl(th.line), linewidth=0.6)
    y = list(range(len(cats)))
    colors = [_mpl(th.accent) if c in highlight else (_mpl(th.negative) if v < 0 else _mpl(th.primary)) for c, v in zip(cats, vals)]
    bars = ax.barh(y, vals, height=0.62, color=colors)
    dec = _auto_decimals(vals, spec.decimals)
    if spec.value_labels:
        ax.bar_label(bars, labels=[fmt_num(v, dec, spec.unit) for v in vals], padding=3, fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
    ax.set_yticks(y)
    ax.set_yticklabels([_wrap(c, 14) for c in cats])
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_xlabel(_axis_title(spec.x_label or spec.y_label, spec.unit), fontsize=typo.chart_label)
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(min(0, xmin), xmax * 1.12 if xmax > 0 else xmax)


def render_line(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    cats, series = _series(spec.data, "line")
    ax = cv.axes()
    n, m = len(cats), len(series)
    x = list(range(n))
    colors = _color_cycle(th, m)
    dec = _auto_decimals([v for _, vals in series for v in vals], spec.decimals)
    for i, (name, vals) in enumerate(series):
        ax.plot(x, vals, color=colors[i], linewidth=1.8, linestyle=LINESTYLES[i] if m >= 3 else "-",
                marker=MARKERS[i], markersize=3.8, label=name)
        if spec.value_labels and n <= 12 and m <= 2:
            for xi, v in zip(x, vals):
                if not math.isnan(v):
                    ax.annotate(fmt_num(v, dec, spec.unit), (xi, v), textcoords="offset points", xytext=(0, 6), ha="center",
                                fontsize=max(8.0, typo.chart_label), color=colors[i])
    ax.set_xticks(x)
    ax.set_xticklabels([_wrap(c, _cat_wrap_width(n, spec.width_mm, typo.chart_tick)) for c in cats])
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_ylabel(_axis_title(spec.y_label, spec.unit), fontsize=typo.chart_label)
    if spec.x_label:
        ax.set_xlabel(spec.x_label, fontsize=typo.chart_label)
    ax.margins(x=0.04, y=0.15)
    if bool(spec.data.get("fromZero", False)):
        ax.set_ylim(bottom=0)
    if m >= 3 and spec.data.get("endLabels", True):
        ax.margins(x=0.08)
        _label_series_end(ax, x, [vals for _, vals in series], [nm for nm, _ in series], colors, typo.chart_tick)
        if spec.legend == "top":
            cv.legend(ax, m)
    elif m > 1:
        cv.legend(ax, m)


def render_area(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    cats, series = _series(spec.data, "area")
    stacked = bool(spec.data.get("stacked", True))
    ax = cv.axes()
    n, m = len(cats), len(series)
    x = list(range(n))
    colors = _color_cycle(th, m)
    clean = [[0 if math.isnan(v) else v for v in vals] for _, vals in series]
    if stacked:
        polys = ax.stackplot(x, *clean, colors=colors, labels=[nm for nm, _ in series], alpha=0.9, edgecolor="white", linewidth=0.6)
        for i, poly in enumerate(polys):
            if m >= 3:
                poly.set_hatch(HATCHES[i])
        tops = [sum(c[j] for c in clean) for j in range(n)]
        if spec.value_labels and n <= 12:
            dec = _auto_decimals(tops, spec.decimals)
            for xi, t in zip(x, tops):
                ax.annotate(fmt_num(t, dec, spec.unit), (xi, t), textcoords="offset points", xytext=(0, 5), ha="center",
                            fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
        ax.set_ylim(0, max(tops) * 1.15 if max(tops) > 0 else 1)
    else:
        for i, (name, vals) in enumerate(series):
            ax.fill_between(x, clean[i], alpha=0.25, color=colors[i])
            ax.plot(x, clean[i], color=colors[i], linewidth=1.6, linestyle=LINESTYLES[i] if m >= 3 else "-", label=name)
        ax.set_ylim(bottom=0)
    ax.set_xticks(x)
    ax.set_xticklabels([_wrap(c, _cat_wrap_width(n, spec.width_mm, typo.chart_tick)) for c in cats])
    ax.set_xlim(0, n - 1)
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_ylabel(_axis_title(spec.y_label, spec.unit), fontsize=typo.chart_label)
    if spec.x_label:
        ax.set_xlabel(spec.x_label, fontsize=typo.chart_label)
    if m > 1:
        cv.legend(ax, m)


def render_pie_as_bar(cv: Canvas, spec: ChartSpec) -> None:
    """饼图替代：按占比降序的条形 + 累计占比线（帕累托）。"""
    th, typo = cv.th, cv.typo
    cats, series = _series(spec.data, "pie_as_bar", single=True, max_series=1)
    vals = [0 if math.isnan(v) else v for v in series[0][1]]
    if any(v < 0 for v in vals):
        raise ChartDataError("pie_as_bar 的 values 不能为负数（占比图）")
    total = sum(vals)
    if total <= 0:
        raise ChartDataError("pie_as_bar 的 values 合计必须大于 0")
    pairs = sorted(zip(cats, vals), key=lambda p: p[1], reverse=True)
    cats = [p[0] for p in pairs]
    vals = [p[1] for p in pairs]
    shares = [v / total * 100 for v in vals]
    cum = [sum(shares[: i + 1]) for i in range(len(shares))]
    ax = cv.axes()
    x = list(range(len(cats)))
    colors = [_mpl(th.primary)] * len(cats)
    highlight = set(str(h) for h in (spec.data.get("highlight") or []))
    colors = [_mpl(th.accent) if c in highlight else col for c, col in zip(cats, colors)]
    bars = ax.bar(x, vals, width=0.6, color=colors)
    dec = _auto_decimals(vals, spec.decimals)
    is_pct = spec.unit in ("%", "％")
    if spec.value_labels:
        # 数据本身已是百分比时不再重复标占比
        labels = [f"{s:.1f}%" if is_pct else f"{fmt_num(v, dec, '')}\n({s:.1f}%)" for v, s in zip(vals, shares)]
        ax.bar_label(bars, labels=labels, padding=2, fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
    ax.set_xticks(x)
    ax.set_xticklabels([_wrap(c, _cat_wrap_width(len(cats), spec.width_mm, typo.chart_tick)) for c in cats])
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_ylabel(_axis_title(spec.y_label, spec.unit), fontsize=typo.chart_label)
    ax.set_ylim(0, max(vals) * 1.3)
    ax2 = ax.twinx()
    ax2.plot(x, cum, color=_mpl(th.accent), marker="o", markersize=3.5, linewidth=1.6, label="累计占比")
    ax2.set_ylim(0, 105)
    ax2.set_ylabel("累计占比（%）", fontsize=typo.chart_label, color=_mpl(th.accent))
    ax2.tick_params(labelsize=typo.chart_tick, colors=_mpl(th.accent))
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_color(_mpl(th.accent))
    ax2.spines["left"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    for xi, c in zip(x, cum):
        ax2.annotate(f"{c:.0f}%", (xi, c), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=typo.chart_tick, color=_mpl(th.accent))
    bar_name = series[0][0] or spec.y_label or ("占比" if is_pct else "数值")
    handles = [patches.Patch(color=_mpl(th.primary), label=bar_name), plt.Line2D([0], [0], color=_mpl(th.accent), marker="o", label="累计占比")]
    cv.legend(ax, 2, handles=handles)


def render_donut(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    cats, series = _series(spec.data, "donut", single=True, max_series=1)
    if len(cats) > 4:
        spec.warnings.append(f"donut 有 {len(cats)} 类（>4），已自动改为 pie_as_bar")
        render_pie_as_bar(cv, spec)
        return
    vals = [0 if math.isnan(v) else v for v in series[0][1]]
    if any(v < 0 for v in vals) or sum(vals) <= 0:
        raise ChartDataError("donut 的 values 须为非负且合计 > 0")
    total = sum(vals)
    ax = cv.fig.add_subplot(111)
    ax.set_aspect("equal")
    colors = _color_cycle(th, len(cats))
    wedges, _ = ax.pie(vals, colors=colors, startangle=90, counterclock=False,
                       wedgeprops=dict(width=0.38, edgecolor="white", linewidth=1.5))
    for i, w in enumerate(wedges):
        if len(cats) >= 3:
            w.set_hatch(HATCHES[i])
    dec = _auto_decimals(vals, spec.decimals)
    is_pct = spec.unit in ("%", "％")
    for w, c, v in zip(wedges, cats, vals):
        ang = math.radians((w.theta2 + w.theta1) / 2)
        xx, yy = math.cos(ang), math.sin(ang)
        ha = "left" if xx >= 0 else "right"
        value_text = f"{v / total * 100:.1f}%" if is_pct else f"{fmt_num(v, dec, '')}{spec.unit}（{v / total * 100:.1f}%）"
        ax.annotate(f"{c}\n{value_text}",
                    xy=(xx * 0.98, yy * 0.98), xytext=(1.25 * xx, 1.15 * yy), ha=ha, va="center",
                    fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink),
                    arrowprops=dict(arrowstyle="-", color=_mpl(th.secondary_light), lw=0.7, connectionstyle=f"angle,angleA=0,angleB={math.degrees(ang):.0f}"))
    center = spec.data.get("centerText")
    if center is None:
        center = "" if is_pct else f"合计\n{fmt_num(total, dec, '')}{spec.unit}"
    ax.text(0, 0, str(center), ha="center", va="center", fontsize=typo.chart_title, color=_mpl(th.primary), fontweight="bold")
    ax.set_xlim(-1.9, 1.9)
    ax.set_ylim(-1.3, 1.3)
    ax.axis("off")


def render_waterfall(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    items = _require(spec.data, "items", "waterfall")
    if not isinstance(items, list) or len(items) < 2:
        raise ChartDataError(f"waterfall 的 items 至少 2 项。结构：{DATA_SPECS['waterfall']}")
    names: List[str] = []
    bottoms: List[float] = []
    heights: List[float] = []
    colors: List[str] = []
    labels: List[str] = []
    running = 0.0
    for i, it in enumerate(items):
        if not isinstance(it, dict) or "name" not in it:
            raise ChartDataError(f"items 第 {i + 1} 项需要 {{name, value}}")
        is_total = bool(it.get("total", False))
        if is_total:
            v = float(it["value"]) if isinstance(it.get("value"), (int, float)) else running
            names.append(str(it["name"]))
            bottoms.append(min(0.0, v))
            heights.append(abs(v))
            colors.append(_mpl(th.primary))
            labels.append(fmt_num(v, spec.decimals, spec.unit))
            running = v
        else:
            if not isinstance(it.get("value"), (int, float)):
                raise ChartDataError(f"items 第 {i + 1} 项「{it['name']}」的 value 必须是数字")
            v = float(it["value"])
            names.append(str(it["name"]))
            lo, hi = sorted((running, running + v))
            bottoms.append(lo)
            heights.append(hi - lo)
            colors.append(_mpl(th.positive) if v >= 0 else _mpl(th.negative))
            labels.append(("+" if v > 0 else "") + fmt_num(v, spec.decimals, spec.unit))
            running += v
    ax = cv.axes()
    x = list(range(len(names)))
    bars = ax.bar(x, heights, bottom=bottoms, width=0.6, color=colors)
    # 连接线：直观表达“上一柱终点 = 下一柱起点”
    tops = [b + h if c != _mpl(th.negative) else b for b, h, c in zip(bottoms, heights, colors)]
    level = 0.0
    for i in range(len(x)):
        if i > 0:
            ax.plot([x[i - 1] + 0.3, x[i] - 0.3], [level, level], color=_mpl(th.muted), linewidth=0.7, linestyle=":")
        level = tops[i] if colors[i] != _mpl(th.primary) else (bottoms[i] + heights[i] if bottoms[i] >= 0 else bottoms[i])
    if spec.value_labels:
        for xi, b, h, lab, col in zip(x, bottoms, heights, labels, colors):
            ax.text(xi, b + h, lab, ha="center", va="bottom", fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
    ax.axhline(0, color=_mpl(th.secondary), linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([_wrap(nm, _cat_wrap_width(len(names), spec.width_mm, typo.chart_tick)) for nm in names])
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_ylabel(_axis_title(spec.y_label, spec.unit), fontsize=typo.chart_label)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(min(0, ymin), ymax * 1.1 if ymax > 0 else ymax)
    handles = [patches.Patch(color=_mpl(th.positive), label="增加"), patches.Patch(color=_mpl(th.negative), label="减少"), patches.Patch(color=_mpl(th.primary), label="合计")]
    cv.legend(ax, 3, handles=handles)


def render_funnel(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    stages = _require(spec.data, "stages", "funnel")
    if not isinstance(stages, list) or len(stages) < 2:
        raise ChartDataError(f"funnel 的 stages 至少 2 层。结构：{DATA_SPECS['funnel']}")
    names, vals = [], []
    for i, s in enumerate(stages):
        if not isinstance(s, dict) or not isinstance(s.get("value"), (int, float)) or "name" not in s:
            raise ChartDataError(f"stages 第 {i + 1} 层需要 {{name, value(数字)}}")
        names.append(str(s["name"]))
        vals.append(float(s["value"]))
    if vals[0] <= 0 or any(v < 0 for v in vals):
        raise ChartDataError("funnel 各层数值须为非负且首层 > 0")
    ax = cv.fig.add_subplot(111)
    ax.axis("off")
    n = len(vals)
    seq = [_mpl(c) for c in reversed(th.sequential)][:n] if n <= len(th.sequential) else [_mpl(th.primary)] * n
    maxv = vals[0]
    dec = _auto_decimals(vals, spec.decimals)
    # 层间比例悬殊（如曝光→复购差 3 个量级）时线性宽度只剩细条，改用平方根缩放并在页脚注明
    sqrt_scale = min(vals) / maxv < 0.05
    if sqrt_scale:
        cv.footer_lines.insert(0, "注：各层宽度按数值平方根缩放，以保证可读；比例以标注数字为准")
    for i, (nm, v) in enumerate(zip(names, vals)):
        ratio = v / maxv
        w = max(math.sqrt(ratio) if sqrt_scale else ratio, 0.12)
        y = n - 1 - i
        ax.add_patch(patches.Rectangle((0.5 - w / 2, y + 0.08), w, 0.84, color=seq[i]))
        txt_color = "white" if i < n * 0.6 else _mpl(th.ink)
        ax.text(0.5, y + 0.5, f"{nm}  {fmt_num(v, dec, '')}{spec.unit if spec.unit not in ('%', '％') else '%'}", ha="center", va="center",
                fontsize=max(8.0, typo.chart_label), color=txt_color if w > 0.45 else _mpl(th.ink), fontweight="bold" if i == 0 else "normal")
        ax.text(1.03, y + 0.5, f"占首层 {v / maxv * 100:.1f}%", ha="left", va="center", fontsize=typo.chart_tick, color=_mpl(th.muted))
        if i > 0 and vals[i - 1] > 0:
            ax.text(-0.03, y + 0.98, f"转化率 {v / vals[i - 1] * 100:.1f}%", ha="right", va="center", fontsize=typo.chart_tick, color=_mpl(th.primary))
    ax.set_xlim(-0.3, 1.35)
    ax.set_ylim(0, n)


def render_tornado(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    base = _require(spec.data, "base", "tornado")
    if isinstance(base, bool) or not isinstance(base, (int, float)):
        raise ChartDataError("tornado 的 base 必须是数字（基准情形的结果值）")
    items = _require(spec.data, "items", "tornado")
    if not isinstance(items, list) or not items:
        raise ChartDataError(f"tornado 的 items 不能为空。结构：{DATA_SPECS['tornado']}")
    rows = []
    for i, it in enumerate(items):
        if not isinstance(it, dict) or not all(isinstance(it.get(k), (int, float)) for k in ("low", "high")) or "name" not in it:
            raise ChartDataError(f"items 第 {i + 1} 项需要 {{name, low(数字), high(数字)}}")
        rows.append((str(it["name"]), float(it["low"]), float(it["high"])))
    rows.sort(key=lambda r: abs(r[2] - r[1]))  # 最敏感的在最上面（barh 从下往上画）
    ax = cv.axes()
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color=_mpl(th.line), linewidth=0.6)
    y = list(range(len(rows)))
    dec = _auto_decimals([v for r in rows for v in r[1:]] + [float(base)], spec.decimals)
    for yi, (nm, lo, hi) in zip(y, rows):
        left, right = min(lo, hi), max(lo, hi)
        ax.barh(yi, left - base, left=base, height=0.6, color=_mpl(th.negative), label="不利情形" if yi == 0 else None)
        ax.barh(yi, right - base, left=base, height=0.6, color=_mpl(th.positive), label="有利情形" if yi == 0 else None)
        if spec.value_labels:
            ax.text(left, yi, fmt_num(left, dec, spec.unit) + " ", ha="right", va="center", fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
            ax.text(right, yi, " " + fmt_num(right, dec, spec.unit), ha="left", va="center", fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
    ax.axvline(base, color=_mpl(th.primary), linewidth=1.2)
    ax.text(base, len(rows) - 0.35, f"基准 {fmt_num(float(base), dec, spec.unit)}", ha="center", va="bottom", fontsize=typo.chart_tick, color=_mpl(th.primary))
    ax.set_yticks(y)
    ax.set_yticklabels([_wrap(r[0], 12) for r in rows])
    ax.xaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_xlabel(_axis_title(spec.x_label or spec.y_label, spec.unit), fontsize=typo.chart_label)
    lo_all = min(min(r[1], r[2]) for r in rows)
    hi_all = max(max(r[1], r[2]) for r in rows)
    pad = (hi_all - lo_all) * 0.18 or 1
    ax.set_xlim(lo_all - pad, hi_all + pad)
    ax.set_ylim(-0.6, len(rows) + 0.1)
    cv.legend(ax, 2)


def render_breakeven(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    d = spec.data
    f = _require(d, "fixedCost", "breakeven")
    p = _require(d, "unitPrice", "breakeven")
    v = _require(d, "unitVarCost", "breakeven")
    for k, val in (("fixedCost", f), ("unitPrice", p), ("unitVarCost", v)):
        if isinstance(val, bool) or not isinstance(val, (int, float)) or val < 0:
            raise ChartDataError(f"breakeven 的 {k} 必须是非负数字")
    if p <= v:
        raise ChartDataError("unitPrice 必须大于 unitVarCost，否则不存在盈亏平衡点")
    q_be = f / (p - v)
    max_units = d.get("maxUnits")
    if max_units is None:
        max_units = q_be * 2 if q_be > 0 else 100
    if not isinstance(max_units, (int, float)) or max_units <= 0:
        raise ChartDataError("maxUnits 必须是正数")
    unit_name = str(d.get("unitName") or "件")
    xs = [max_units * i / 100 for i in range(101)]
    rev = [p * q for q in xs]
    cost = [f + v * q for q in xs]
    ax = cv.axes()
    ax.plot(xs, rev, color=_mpl(th.primary), linewidth=1.8, label="销售收入")
    ax.plot(xs, cost, color=_mpl(th.negative), linewidth=1.8, linestyle="--", label="总成本")
    ax.plot(xs, [f] * len(xs), color=_mpl(th.secondary), linewidth=1.2, linestyle=":", label="固定成本")
    ax.fill_between(xs, rev, cost, where=[r >= c for r, c in zip(rev, cost)], color=_mpl(th.positive), alpha=0.12, interpolate=True)
    ax.fill_between(xs, rev, cost, where=[r < c for r, c in zip(rev, cost)], color=_mpl(th.negative), alpha=0.10, interpolate=True)
    be_amt = p * q_be
    ax.scatter([q_be], [be_amt], color=_mpl(th.accent), zorder=5, s=36, edgecolor=_mpl(th.primary))
    ax.annotate(f"盈亏平衡点\n{fmt_num(q_be, 0)} {unit_name} / {fmt_num(be_amt, None, '')}{spec.unit}",
                (q_be, be_amt), textcoords="offset points", xytext=(14, -28), fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink),
                arrowprops=dict(arrowstyle="-", color=_mpl(th.muted), lw=0.7))
    ax.text(xs[-1], rev[-1], " 盈利区", fontsize=typo.chart_tick, color=_mpl(th.positive), va="center", ha="left")
    ax.set_xlim(0, max_units)
    ax.set_ylim(0, max(rev[-1], cost[-1]) * 1.1)
    ax.xaxis.set_major_formatter(_tick_formatter(""))
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_xlabel(spec.x_label or f"销量（{unit_name}）", fontsize=typo.chart_label)
    ax.set_ylabel(_axis_title(spec.y_label or "金额", spec.unit), fontsize=typo.chart_label)
    ax.margins(x=0.08)
    cv.legend(ax, 3)


def render_cashflow_combo(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    cats = _as_strs(_require(spec.data, "categories", "cashflow_combo"), "categories")
    vals = _as_floats(_require(spec.data, "values", "cashflow_combo"), "values")
    if len(vals) != len(cats):
        raise ChartDataError(f"values 有 {len(vals)} 个值，与 categories 的 {len(cats)} 个不一致")
    if spec.data.get("cumulative") is not None:
        cum = _as_floats(spec.data["cumulative"], "cumulative")
        if len(cum) != len(cats):
            raise ChartDataError("cumulative 长度须与 categories 一致")
    else:
        cum = [sum(vals[: i + 1]) for i in range(len(vals))]
    ax = cv.axes()
    x = list(range(len(cats)))
    bars = ax.bar(x, vals, width=0.6, color=[_mpl(th.positive) if v >= 0 else _mpl(th.negative) for v in vals], label="当期现金流")
    ax.plot(x, cum, color=_mpl(th.primary), linewidth=1.8, marker="o", markersize=3.8, label="累计现金流")
    dec = _auto_decimals(vals + cum, spec.decimals)
    if spec.value_labels and len(cats) <= 16:
        ax.bar_label(bars, labels=[fmt_num(v, dec, spec.unit) for v in vals], padding=2, fontsize=max(8.0, typo.chart_label), color=_mpl(th.ink))
    i_min = min(range(len(cum)), key=lambda i: cum[i])
    ax.annotate(f"现金最低点 {fmt_num(cum[i_min], dec, spec.unit)}", (x[i_min], cum[i_min]), textcoords="offset points", xytext=(0, -16),
                ha="center", fontsize=max(8.0, typo.chart_label), color=_mpl(th.primary), fontweight="bold")
    for i in range(1, len(cum)):
        if cum[i - 1] < 0 <= cum[i]:
            # 标注放在穿越点左上方，避开柱顶的数据标签
            ax.annotate(f"累计转正：{cats[i]}", (x[i], cum[i]), textcoords="offset points", xytext=(-10, 26), ha="right",
                        fontsize=typo.chart_tick, color=_mpl(th.accent), fontweight="bold",
                        arrowprops=dict(arrowstyle="-", color=_mpl(th.accent), lw=0.6))
            break
    ax.axhline(0, color=_mpl(th.secondary), linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([_wrap(c, _cat_wrap_width(len(cats), spec.width_mm, typo.chart_tick)) for c in cats])
    ax.yaxis.set_major_formatter(_tick_formatter(spec.unit))
    ax.set_ylabel(_axis_title(spec.y_label, spec.unit), fontsize=typo.chart_label)
    if spec.x_label:
        ax.set_xlabel(spec.x_label, fontsize=typo.chart_label)
    ax.margins(y=0.2)
    cv.legend(ax, 2)


def _points(data: Dict[str, Any], ctype: str) -> List[Dict[str, Any]]:
    pts = _require(data, "points", ctype)
    if not isinstance(pts, list) or not pts:
        raise ChartDataError(f"points 不能为空。结构：{DATA_SPECS[ctype]}")
    out = []
    for i, p in enumerate(pts):
        if not isinstance(p, dict) or "name" not in p or not all(isinstance(p.get(k), (int, float)) for k in ("x", "y")):
            raise ChartDataError(f"points 第 {i + 1} 项需要 {{name, x(数字), y(数字)}}")
        size = p.get("size")
        if size is not None and (not isinstance(size, (int, float)) or size <= 0):
            raise ChartDataError(f"points「{p['name']}」的 size 必须是正数")
        out.append({"name": str(p["name"]), "x": float(p["x"]), "y": float(p["y"]), "size": float(size) if size else None, "self": bool(p.get("self", False))})
    return out


def _scatter_points(ax, cv: Canvas, pts: List[Dict[str, Any]]) -> None:
    th, typo = cv.th, cv.typo
    sizes = [p["size"] for p in pts if p["size"]]
    smax = max(sizes) if sizes else None
    for p in pts:
        s = 60 + 540 * (p["size"] / smax) if (smax and p["size"]) else 90
        color = _mpl(th.accent) if p["self"] else _mpl(th.primary)
        ax.scatter([p["x"]], [p["y"]], s=s, color=color, alpha=0.85 if p["self"] else 0.6, edgecolor="white" if not p["self"] else _mpl(th.primary), linewidth=1.0, zorder=4)
        ax.annotate(p["name"], (p["x"], p["y"]), textcoords="offset points", xytext=(0, 6 + math.sqrt(s) / 2), ha="center",
                    fontsize=max(8.0, typo.chart_label), color=_mpl(th.primary) if p["self"] else _mpl(th.ink), fontweight="bold" if p["self"] else "normal", zorder=5)


def render_scatter_matrix(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    pts = _points(spec.data, "scatter_matrix")
    quadrants = spec.data.get("quadrants")
    if quadrants is not None and (not isinstance(quadrants, list) or len(quadrants) != 4):
        raise ChartDataError("quadrants 须为 4 个标签：[右上, 左上, 左下, 右下]")
    ax = cv.axes()
    ax.yaxis.grid(False)
    xs = [p["x"] for p in pts]
    ys = [p["y"] for p in pts]
    xr = spec.data.get("xRange") or [min(xs), max(xs)]
    yr = spec.data.get("yRange") or [min(ys), max(ys)]
    xpad = (xr[1] - xr[0]) * 0.18 or 1
    ypad = (yr[1] - yr[0]) * 0.22 or 1
    x0, x1 = xr[0] - xpad, xr[1] + xpad
    y0, y1 = yr[0] - ypad, yr[1] + ypad
    xmid = float(spec.data.get("xMid", (xr[0] + xr[1]) / 2))
    ymid = float(spec.data.get("yMid", (yr[0] + yr[1]) / 2))
    ax.axvline(xmid, color=_mpl(th.secondary_light), linestyle="--", linewidth=0.9)
    ax.axhline(ymid, color=_mpl(th.secondary_light), linestyle="--", linewidth=0.9)
    _scatter_points(ax, cv, pts)
    if quadrants:
        corners = [(x1, y1, "right", "top"), (x0, y1, "left", "top"), (x0, y0, "left", "bottom"), (x1, y0, "right", "bottom")]
        for label, (cx, cy, ha, va) in zip(quadrants, corners):
            ax.text(cx, cy, f" {label} ", ha=ha, va=va, fontsize=typo.chart_tick, color=_mpl(th.muted), style="normal")
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel(str(spec.data.get("xLabel") or spec.x_label), fontsize=typo.chart_label)
    ax.set_ylabel(str(spec.data.get("yLabel") or spec.y_label), fontsize=typo.chart_label)
    ax.xaxis.set_major_formatter(_tick_formatter(""))
    ax.yaxis.set_major_formatter(_tick_formatter(""))
    if any(p["self"] for p in pts):
        handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=_mpl(th.accent), markersize=8, label="本品牌 / 本项目"),
                   plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=_mpl(th.primary), markersize=8, label="竞争者")]
        cv.legend(ax, 2, handles=handles)


def render_positioning_map(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    pts = _points(spec.data, "positioning_map")
    xa = _require(spec.data, "xAxis", "positioning_map")
    ya = _require(spec.data, "yAxis", "positioning_map")
    for nm, a in (("xAxis", xa), ("yAxis", ya)):
        if not isinstance(a, dict) or "low" not in a or "high" not in a:
            raise ChartDataError(f"{nm} 需要 {{low: '左/下极', high: '右/上极'}}")
    xr = spec.data.get("xRange") or [-1, 1]
    yr = spec.data.get("yRange") or [-1, 1]
    for p in pts:
        if not (xr[0] <= p["x"] <= xr[1] and yr[0] <= p["y"] <= yr[1]):
            raise ChartDataError(f"点「{p['name']}」坐标 ({p['x']}, {p['y']}) 超出范围 x{xr} y{yr}，请调整或给出 xRange/yRange")
    ax = cv.fig.add_subplot(111)
    ax.axis("off")
    xmid, ymid = (xr[0] + xr[1]) / 2, (yr[0] + yr[1]) / 2
    xpad, ypad = (xr[1] - xr[0]) * 0.22, (yr[1] - yr[0]) * 0.22
    ax.set_xlim(xr[0] - xpad, xr[1] + xpad)
    ax.set_ylim(yr[0] - ypad, yr[1] + ypad)
    ax.annotate("", xy=(xr[1] + xpad * 0.45, ymid), xytext=(xr[0] - xpad * 0.45, ymid), arrowprops=dict(arrowstyle="<->", color=_mpl(th.secondary), lw=1.0))
    ax.annotate("", xy=(xmid, yr[1] + ypad * 0.45), xytext=(xmid, yr[0] - ypad * 0.45), arrowprops=dict(arrowstyle="<->", color=_mpl(th.secondary), lw=1.0))
    lab_kw = dict(fontsize=typo.chart_label, color=_mpl(th.primary), fontweight="bold")
    ax.text(xr[0] - xpad * 0.5, ymid, _wrap(str(xa["low"]), 6), ha="right", va="center", **lab_kw)
    ax.text(xr[1] + xpad * 0.5, ymid, _wrap(str(xa["high"]), 6), ha="left", va="center", **lab_kw)
    ax.text(xmid, yr[1] + ypad * 0.55, str(ya["high"]), ha="center", va="bottom", **lab_kw)
    ax.text(xmid, yr[0] - ypad * 0.55, str(ya["low"]), ha="center", va="top", **lab_kw)
    ax.add_patch(patches.Rectangle((xr[0], yr[0]), xr[1] - xr[0], yr[1] - yr[0], fill=True, facecolor=_mpl(th.paper_alt), edgecolor="none", zorder=0))
    _scatter_points(ax, cv, pts)
    if any(p["self"] for p in pts):
        handles = [plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=_mpl(th.accent), markersize=8, label="本品牌"),
                   plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=_mpl(th.primary), markersize=8, label="竞争品牌")]
        ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 0), fontsize=typo.chart_tick, frameon=False)


def _parse_month(s: Any, name: str) -> Tuple[int, int, int]:
    """'2026-04' / '2026-04-15' / '2026Q2' → (年, 月, 日)。"""
    text = str(s).strip()
    m = re.fullmatch(r"(\d{4})[-/.年](\d{1,2})(?:[-/.月](\d{1,2})日?)?", text)
    if m:
        y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3) or 1)
        if not 1 <= mo <= 12 or not 1 <= d <= 31:
            raise ChartDataError(f"{name} 日期「{text}」月/日超出范围")
        return y, mo, d
    q = re.fullmatch(r"(\d{4})\s*[Qq季]\s*([1-4])", text) or re.fullmatch(r"(\d{4})\s*第?([1-4])季度?", text)
    if q:
        return int(q.group(1)), (int(q.group(2)) - 1) * 3 + 1, 1
    raise ChartDataError(f"{name} 日期「{text}」格式应为 YYYY-MM 或 YYYY-MM-DD")


def _month_pos(y: int, mo: int, d: int, y0: int, m0: int) -> float:
    return (y - y0) * 12 + (mo - m0) + (d - 1) / 30.0


def render_gantt(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    tasks = _require(spec.data, "tasks", "gantt")
    if not isinstance(tasks, list) or not tasks:
        raise ChartDataError(f"gantt 的 tasks 不能为空。结构：{DATA_SPECS['gantt']}")
    rows = []
    for i, t in enumerate(tasks):
        if not isinstance(t, dict) or "name" not in t or "start" not in t or "end" not in t:
            raise ChartDataError(f"tasks 第 {i + 1} 项需要 {{name, start, end}}")
        s = _parse_month(t["start"], f"tasks[{i + 1}].start")
        e = _parse_month(t["end"], f"tasks[{i + 1}].end")
        if s > e:
            raise ChartDataError(f"任务「{t['name']}」的 start 晚于 end")
        rows.append((str(t["name"]), s, e, str(t.get("group") or ""), bool(t.get("done", False))))
    milestones = []
    for i, ms in enumerate(spec.data.get("milestones") or []):
        if not isinstance(ms, dict) or "name" not in ms or "date" not in ms:
            raise ChartDataError(f"milestones 第 {i + 1} 项需要 {{name, date}}")
        milestones.append((str(ms["name"]), _parse_month(ms["date"], f"milestones[{i + 1}].date")))
    all_dates = [r[1] for r in rows] + [r[2] for r in rows] + [m[1] for m in milestones]
    y0, m0 = min(all_dates)[:2]
    y1, m1 = max(all_dates)[:2]
    span = (y1 - y0) * 12 + (m1 - m0) + 1
    groups = []
    for r in rows:
        if r[3] and r[3] not in groups:
            groups.append(r[3])
    if len(groups) > 6:
        raise ChartDataError(f"gantt 的 group 有 {len(groups)} 种，超过 6 色上限")
    gcolor = {g: _mpl(th.series[i]) for i, g in enumerate(groups)}
    ax = cv.axes()
    ax.yaxis.grid(False)
    ax.xaxis.grid(True, color=_mpl(th.line), linewidth=0.6)
    for i, (nm, s, e, g, done) in enumerate(rows):
        x_start = _month_pos(*s, y0, m0)
        x_end = _month_pos(*e, y0, m0) + (1 if len(str(tasks[i]["end"])) <= 7 else 0)  # 月粒度的 end 含整月
        x_end = max(x_end, x_start + 0.25)
        color = gcolor.get(g, _mpl(th.primary))
        ax.barh(i, x_end - x_start, left=x_start, height=0.6, color=color, alpha=0.95 if not done else 0.55, hatch="//" if done else "", edgecolor="white")
        period = f"{s[0]}.{s[1]:02d}–{e[0]}.{e[1]:02d}" + ("  ✓ 已完成" if done else "")
        ax.text(x_end + 0.12, i, period, va="center", ha="left", fontsize=typo.chart_tick - 0.5, color=_mpl(th.muted))
    # 里程碑放在最上面一行（y 轴已反转，负值在顶部），不与月份刻度打架
    for k, (nm, dt) in enumerate(milestones):
        xm = _month_pos(*dt, y0, m0)
        ax.plot([xm], [-0.75], marker="D", color=_mpl(th.accent), markersize=6, zorder=5)
        ax.text(xm, -1.05, nm, ha="center", va="bottom", fontsize=typo.chart_tick, color=_mpl(th.accent), fontweight="bold")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([_wrap(r[0], 11) for r in rows])
    ax.set_ylim(len(rows) - 0.4, -1.7 if milestones else -0.6)
    step = 1 if span <= 8 else (2 if span <= 16 else (3 if span <= 40 else 6))
    ticks = list(range(0, span + 1, step))
    ax.set_xticks(ticks)
    labels = []
    for t in ticks:
        yy, mm = y0 + (m0 - 1 + t) // 12, (m0 - 1 + t) % 12 + 1
        labels.append(f"{yy}-{mm:02d}")
    ax.set_xticklabels(labels, rotation=0 if len(ticks) <= 12 else 45, ha="center" if len(ticks) <= 12 else "right")
    ax.set_xlim(0, span + max(2.5, span * 0.3))
    for mon in range(span + 1):
        ax.axvline(mon, color=_mpl(th.line), linewidth=0.4, zorder=0)
    if groups:
        handles = [patches.Patch(color=gcolor[g], label=g) for g in groups]
        if any(r[4] for r in rows):
            handles.append(patches.Patch(facecolor="#BBBBBB", hatch="//", edgecolor="white", label="已完成"))
        cv.legend(ax, len(handles), handles=handles)


def render_radar(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    axes_names = _as_strs(_require(spec.data, "axes", "radar"), "axes")
    if len(axes_names) < 3:
        raise ChartDataError("radar 至少需要 3 个维度")
    raw = _require(spec.data, "series", "radar")
    if not isinstance(raw, list) or not raw:
        raise ChartDataError("radar 的 series 不能为空")
    if len(raw) > 4:
        raise ChartDataError(f"radar 序列 {len(raw)} 个过多（≤ 4 才可读），请拆图")
    series = []
    for i, s in enumerate(raw):
        if not isinstance(s, dict):
            raise ChartDataError(f"series 第 {i + 1} 项必须是 {{name, values}}")
        vals = _as_floats(s.get("values"), f"series[{i + 1}].values")
        if len(vals) != len(axes_names):
            raise ChartDataError(f"序列「{s.get('name', i + 1)}」有 {len(vals)} 个值，与 axes 的 {len(axes_names)} 个不一致")
        series.append((str(s.get("name") or f"系列{i + 1}"), vals))
    allv = [v for _, vals in series for v in vals]
    vmax = spec.data.get("max")
    if vmax is None:
        top = max(allv)
        vmax = 5 if top <= 5 else (10 if top <= 10 else (100 if top <= 100 else math.ceil(top / 10) * 10))
    if not isinstance(vmax, (int, float)) or vmax <= 0:
        raise ChartDataError("radar 的 max 必须是正数")
    n = len(axes_names)
    angles = [2 * math.pi * i / n for i in range(n)]
    ax = cv.fig.add_subplot(111, polar=True)
    ax.set_theta_offset(math.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, vmax)
    ticks = [vmax * k / 5 for k in range(1, 6)]
    ax.set_yticks(ticks)
    ax.set_yticklabels([fmt_num(t) for t in ticks], fontsize=typo.chart_tick - 1, color=_mpl(th.muted))
    ax.set_rlabel_position(90 / n + 90)
    ax.set_xticks(angles)
    ax.set_xticklabels([_wrap(a, 6) for a in axes_names], fontsize=typo.chart_label, color=_mpl(th.ink))
    ax.tick_params(axis="x", pad=8)
    ax.grid(color=_mpl(th.line), linewidth=0.7)
    ax.spines["polar"].set_color(_mpl(th.secondary_light))
    colors = _color_cycle(th, len(series))
    for i, (name, vals) in enumerate(series):
        data = vals + vals[:1]
        ang = angles + angles[:1]
        ax.plot(ang, data, color=colors[i], linewidth=1.8, linestyle=LINESTYLES[i] if len(series) >= 3 else "-", marker=MARKERS[i], markersize=3.5, label=name)
        ax.fill(ang, data, color=colors[i], alpha=0.12)
        if spec.value_labels and len(series) <= 2:
            for a, v in zip(angles, vals):
                ax.annotate(fmt_num(v, spec.decimals), (a, v), textcoords="offset points", xytext=(0, 5), ha="center", fontsize=max(8.0, typo.chart_tick), color=colors[i])
    if len(series) > 1 or spec.legend == "bottom":
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=min(len(series), 4), fontsize=typo.chart_tick, frameon=False)


BMC_BLOCKS = [
    ("keyPartners", "重要合作"), ("keyActivities", "关键业务"), ("keyResources", "核心资源"),
    ("valuePropositions", "价值主张"), ("customerRelationships", "客户关系"), ("channels", "渠道通路"),
    ("customerSegments", "客户细分"), ("costStructure", "成本结构"), ("revenueStreams", "收入来源"),
]


def render_bmc(cv: Canvas, spec: ChartSpec) -> None:
    th, typo = cv.th, cv.typo
    blocks: Dict[str, List[str]] = {}
    missing = [zh for key, zh in BMC_BLOCKS if key not in spec.data]
    if missing:
        raise ChartDataError(f"bmc 缺少九宫格字段：{', '.join(missing)}。结构：{DATA_SPECS['bmc']}")
    for key, zh in BMC_BLOCKS:
        v = spec.data[key]
        if isinstance(v, str):
            v = [v]
        if not isinstance(v, list):
            raise ChartDataError(f"bmc 的 {key}（{zh}）必须是字符串数组")
        blocks[key] = [str(x) for x in v]
    ax = cv.fig.add_subplot(111)
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    # 布局：上 4 行 ×5 列（中间两列各拆上下两格），下 2 行两格
    layout = {
        "keyPartners": (0, 2, 2, 4), "keyActivities": (2, 4, 2, 2), "keyResources": (2, 2, 2, 2),
        "valuePropositions": (4, 2, 2, 4), "customerRelationships": (6, 4, 2, 2), "channels": (6, 2, 2, 2),
        "customerSegments": (8, 2, 2, 4), "costStructure": (0, 0, 5, 2), "revenueStreams": (5, 0, 5, 2),
    }
    zh = dict(BMC_BLOCKS)
    fs = max(8.0, typo.chart_tick)
    # 画布实际占图面的比例（扣除页脚），用于把“数据坐标 1 单位”换算成毫米
    plot_h_mm = spec.height_mm * (1 - cv.footer_frac()) * 0.96
    unit_mm_x = spec.width_mm * 0.96 / 10
    unit_mm_y = plot_h_mm / 6
    line_units = (fs * 1.35 * 0.3528) / unit_mm_y
    for key, (x, y, w, h) in layout.items():
        face = _mpl(th.paper_alt) if key in ("valuePropositions",) else "white"
        ax.add_patch(patches.Rectangle((x, y), w, h, facecolor=face, edgecolor=_mpl(th.secondary_light), linewidth=0.9))
        ax.add_patch(patches.Rectangle((x, y + h - 0.42), w, 0.42, facecolor=_mpl(th.primary), edgecolor="none"))
        ax.text(x + 0.12, y + h - 0.21, zh[key], va="center", ha="left", fontsize=typo.chart_label, color="white", fontweight="bold")
        chars = max(4, int((w * unit_mm_x - 3.5) / (fs * 0.3528)))
        max_lines = max(2, int((h - 0.55) / line_units))
        lines: List[str] = []
        for item in blocks[key]:
            wrapped = _wrap_lines(item, chars - 2)
            lines.append("· " + wrapped[0])
            lines.extend("   " + l for l in wrapped[1:])
        if len(lines) > max_lines:
            lines = lines[: max_lines - 1] + ["  …"]
            spec.warnings.append(f"bmc「{zh[key]}」内容过多，已截断（最多约 {max_lines} 行，可增大 height 或精简文字）")
        ax.text(x + 0.12, y + h - 0.52, "\n".join(lines), va="top", ha="left", fontsize=fs, color=_mpl(th.ink), linespacing=1.35)


RENDERERS: Dict[str, Callable[[Canvas, ChartSpec], None]] = {
    "bar": render_bar, "grouped_bar": render_grouped_bar, "stacked_bar": render_stacked_bar, "hbar": render_hbar,
    "line": render_line, "area": render_area, "pie_as_bar": render_pie_as_bar, "donut": render_donut,
    "waterfall": render_waterfall, "funnel": render_funnel, "tornado": render_tornado, "breakeven": render_breakeven,
    "cashflow_combo": render_cashflow_combo, "scatter_matrix": render_scatter_matrix, "positioning_map": render_positioning_map,
    "gantt": render_gantt, "radar": render_radar, "bmc": render_bmc,
}


# ----------------------------------------------------------------------------------------------
# 批量渲染
# ----------------------------------------------------------------------------------------------

def render_one(spec: ChartSpec, out_dir: Path, th: theme_mod.Theme, dpi: int) -> Dict[str, Any]:
    cv = Canvas(spec, th, dpi)
    try:
        RENDERERS[spec.type](cv, spec)
        path = out_dir / f"{spec.id}.png"
        px = cv.save(path)
    finally:
        cv.close()
    return {
        "id": spec.id, "file": path.name, "title": spec.title, "source": spec.source, "type": spec.type,
        "widthMm": spec.width_mm, "heightMm": spec.height_mm, "px": list(px), "warnings": list(spec.warnings),
    }


def resolve_output_dir(cfg: Dict[str, Any], json_path: Path) -> Path:
    out = cfg.get("outputDir")
    if not out:
        return json_path.parent
    p = Path(out)
    if p.is_absolute():
        return p
    root = find_workspace_root(json_path.parent, probe=p)
    cand = root / p
    if cand.exists():
        return cand
    # charts.json 通常就放在 outputDir 里（04_图表/charts.json → outputDir "04_图表"）
    if json_path.parent.name == p.name:
        return json_path.parent
    return cand


def render_charts(cfg: Dict[str, Any], json_path: Path, only: Optional[Sequence[str]] = None) -> Dict[str, Any]:
    th = theme_mod.get_theme(str(cfg.get("theme") or "bizgold"))
    dpi = int(cfg.get("dpi") or th.dpi)
    out_dir = resolve_output_dir(cfg, json_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    fonts.apply_matplotlib("sans")
    plt.rcParams["hatch.linewidth"] = 0.6
    charts = cfg.get("charts")
    if not isinstance(charts, list):
        raise ChartDataError("charts.json 需要 charts 数组")
    wanted = set(only) if only else None
    index_path = out_dir / "index.json"
    index: Dict[str, Any] = {}
    if wanted and index_path.exists():
        try:
            index = read_json(index_path)
        except Exception:
            index = {}
    rendered: List[Dict[str, Any]] = []
    failed: List[Dict[str, Any]] = []
    warnings: List[str] = []
    for i, raw in enumerate(charts):
        cid = raw.get("id") if isinstance(raw, dict) else None
        if wanted and cid not in wanted:
            continue
        try:
            spec = ChartSpec.from_dict(raw, th)
            info = render_one(spec, out_dir, th, dpi)
            rendered.append(info)
            index[spec.id] = {k: info[k] for k in ("file", "title", "source", "type", "widthMm", "heightMm", "px")}
            for w in info["warnings"]:
                warnings.append(f"[{spec.id}] {w}")
            log(f"✓ {spec.id} ({spec.type}) {info['px'][0]}×{info['px'][1]}px")
        except ChartDataError as e:
            failed.append({"id": cid or f"#{i + 1}", "type": raw.get("type") if isinstance(raw, dict) else None, "error": str(e)})
            log(f"✗ {cid or i + 1}: {e}")
        except Exception as e:  # 渲染异常也不应拖垮整批
            failed.append({"id": cid or f"#{i + 1}", "type": raw.get("type") if isinstance(raw, dict) else None,
                           "error": f"渲染异常 {type(e).__name__}: {e}", "trace": traceback.format_exc(limit=3)})
            log(f"✗ {cid or i + 1}: {type(e).__name__}: {e}")
            plt.close("all")
    write_json(index_path, index)
    return {
        "ok": not failed, "outputDir": str(out_dir), "index": str(index_path), "count": len(rendered),
        "rendered": rendered, "failed": failed, "warnings": warnings, "generatedAt": date.today().isoformat(),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    force_utf8_stdout()
    ap = argparse.ArgumentParser(description="charts.json → PNG 图表")
    ap.add_argument("charts_json", help="charts.json 路径")
    ap.add_argument("--only", default="", help="只渲染这些 id（逗号分隔）")
    args = ap.parse_args(argv)
    json_path = Path(args.charts_json).resolve()
    try:
        cfg = read_json(json_path)
    except Exception as e:
        print(json.dumps({"ok": False, "error": f"无法读取 {json_path}: {e}"}, ensure_ascii=False))
        return 1
    only = [s.strip() for s in args.only.split(",") if s.strip()] or None
    try:
        result = render_charts(cfg, json_path, only)
    except ChartDataError as e:
        print(json.dumps({"ok": False, "error": str(e)}, ensure_ascii=False))
        return 1
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["ok"] else 2


if __name__ == "__main__":
    sys.exit(main())
