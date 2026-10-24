# -*- coding: utf-8 -*-
"""matplotlib 商业图表库（matplotlib 3.11.2 验证）。

- setup_chinese_font(): Windows/macOS/Linux 自动选择中文字体，缺失时重建字体缓存
- PALETTES: 商务蓝金 / 咨询风格 配色
- 11 类常用商业图函数，统一 300dpi 输出
"""
from __future__ import annotations

import os
import platform
import shutil
import subprocess
from typing import Sequence

import matplotlib

matplotlib.use("Agg")  # 服务器无显示环境
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

# --------------------------------------------------------------- 字体 ----

FONT_CANDIDATES = {
    "Windows": ["Microsoft YaHei", "SimHei", "DengXian", "SimSun", "KaiTi"],
    "Darwin": ["PingFang SC", "Hiragino Sans GB", "Heiti SC", "STHeiti", "Songti SC"],
    "Linux": ["Noto Sans CJK SC", "Source Han Sans SC", "Source Han Sans CN", "WenQuanYi Micro Hei",
              "WenQuanYi Zen Hei", "AR PL UMing CN", "Droid Sans Fallback"],
}


def _installed_font_names() -> set[str]:
    return {f.name for f in fm.fontManager.ttflist}


def _rebuild_font_cache():
    """apt 装完字体后 matplotlib 缓存(~/.cache/matplotlib/fontlist-*.json)仍是旧的，需重建。"""
    try:
        fm._load_fontmanager(try_read_cache=False)  # 私有 API，但 3.x 一直存在
    except Exception:
        cache_dir = matplotlib.get_cachedir()
        for f in os.listdir(cache_dir):
            if f.startswith("fontlist"):
                os.remove(os.path.join(cache_dir, f))


def setup_chinese_font(prefer: Sequence[str] | None = None, verbose=True) -> str:
    """返回最终选用的中文字体名，并写入 rcParams。找不到时抛 RuntimeError 并给出安装建议。"""
    system = platform.system()
    candidates = list(prefer or []) + FONT_CANDIDATES.get(system, []) + \
        [n for k in FONT_CANDIDATES.values() for n in k]
    names = _installed_font_names()
    chosen = next((c for c in candidates if c in names), None)
    if chosen is None:
        _rebuild_font_cache()
        names = _installed_font_names()
        chosen = next((c for c in candidates if c in names), None)
    if chosen is None:
        hint = {
            "Linux": "sudo apt-get install -y fonts-noto-cjk  （或 fonts-wqy-microhei），然后删除 ~/.cache/matplotlib/fontlist-*.json",
            "Windows": "系统应自带 微软雅黑/SimHei；若是精简版请安装中文语言包",
            "Darwin": "系统应自带 PingFang SC",
        }.get(system, "安装任意 CJK 字体")
        if system == "Linux" and shutil.which("fc-list"):
            out = subprocess.run(["fc-list", ":lang=zh", "family"], capture_output=True, text=True).stdout
            hint += f"\n  fc-list :lang=zh 当前可见: {sorted(set(out.strip().splitlines()))[:8]}"
        raise RuntimeError(f"未找到中文字体，候选: {candidates}\n建议: {hint}")
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False  # 否则负号显示成方框
    if verbose:
        print(f"[bizcharts] 中文字体: {chosen} ({system})")
    return chosen


# --------------------------------------------------------------- 配色 ----

PALETTES = {
    # 商务蓝金：主色藏青 + 点缀金色，适合商业计划书
    "blue_gold": {
        "primary": "#1F3864", "secondary": "#2E75B6", "accent": "#C9A227", "accent2": "#E2B93B",
        "positive": "#2E8B57", "negative": "#C0392B", "neutral": "#7F7F7F", "light": "#D6DCE5",
        "series": ["#1F3864", "#2E75B6", "#C9A227", "#8FAADC", "#7F7F7F", "#BDD7EE", "#E2B93B", "#404040"],
    },
    # 咨询公司风格：深海军蓝 + 电光蓝 + 浅青，少量强调色（参考 MBB 常见配色）
    "consulting": {
        "primary": "#051C2C", "secondary": "#2251FF", "accent": "#00A9F4", "accent2": "#AAE6F0",
        "positive": "#00A651", "negative": "#E5133B", "neutral": "#8C8C8C", "light": "#E6E6E6",
        "series": ["#051C2C", "#2251FF", "#00A9F4", "#AAE6F0", "#8C8C8C", "#C4C4C4", "#034B6F", "#1E6B52"],
    },
}
_P = PALETTES["blue_gold"]


def use_palette(name="blue_gold"):
    global _P
    _P = PALETTES[name]
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=_P["series"])
    return _P


def apply_style(base_font_size=10):
    """统一的商业图外观：去上/右边框、浅网格、标题左对齐。"""
    plt.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 300,
        "font.size": base_font_size,
        "axes.titlesize": base_font_size + 3,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.labelsize": base_font_size,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": "#E0E0E0",
        "grid.linewidth": 0.6,
        "axes.axisbelow": True,
        "legend.frameon": False,
        "legend.fontsize": base_font_size - 1,
        "xtick.labelsize": base_font_size - 1,
        "ytick.labelsize": base_font_size - 1,
    })
    use_palette()


def save_fig(fig, path, dpi=300, transparent=False):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white" if not transparent else "none",
                transparent=transparent)
    plt.close(fig)
    return path


def _source_note(ax, source: str | None):
    if source:
        ax.annotate(f"资料来源：{source}", xy=(0, -0.18), xycoords="axes fraction", fontsize=7.5,
                    color=_P["neutral"], ha="left", va="top")


def _fmt(v, unit=""):
    if abs(v) >= 1e8:
        return f"{v / 1e8:.1f}亿{unit}"
    if abs(v) >= 1e4:
        return f"{v / 1e4:.0f}万{unit}"
    return f"{v:,.0f}{unit}"


# ------------------------------------------------------- 1 市场规模堆叠柱 ----

def market_size_stacked_bar(years: Sequence, segments: dict[str, Sequence[float]], title="市场规模",
                            unit="亿元", cagr_label=True, source=None, figsize=(7, 4.2)):
    fig, ax = plt.subplots(figsize=figsize)
    bottom = np.zeros(len(years))
    x = np.arange(len(years))
    for i, (name, vals) in enumerate(segments.items()):
        vals = np.asarray(vals, dtype=float)
        ax.bar(x, vals, 0.6, bottom=bottom, label=name, color=_P["series"][i % len(_P["series"])])
        bottom += vals
    for xi, tot in zip(x, bottom):
        ax.text(xi, tot * 1.01, f"{tot:,.0f}", ha="center", va="bottom", fontsize=9, fontweight="bold")
    if cagr_label and len(years) > 1 and bottom[0] > 0:
        cagr = (bottom[-1] / bottom[0]) ** (1 / (len(years) - 1)) - 1
        ax.annotate(f"CAGR {cagr:.1%}", xy=(x[-1], bottom[-1]), xytext=(x[0] + 0.3, bottom[-1] * 1.08),
                    arrowprops=dict(arrowstyle="->", color=_P["accent"], lw=1.5), color=_P["accent"],
                    fontsize=10, fontweight="bold")
    ax.set_xticks(x, [str(y) for y in years])
    ax.set_ylabel(unit)
    ax.set_title(title)
    ax.set_ylim(0, bottom.max() * 1.2)
    ax.legend(ncol=min(4, len(segments)), loc="upper left", bbox_to_anchor=(0, -0.08))
    _source_note(ax, source)
    return fig


# ------------------------------------------------------------ 2 份额条形 ----

def share_barh(labels: Sequence[str], values: Sequence[float], title="市场份额", highlight: str | None = None,
               unit="%", source=None, figsize=(6.5, 4)):
    order = np.argsort(values)
    labels = [labels[i] for i in order]
    values = [values[i] for i in order]
    fig, ax = plt.subplots(figsize=figsize)
    colors = [_P["accent"] if l == highlight else _P["secondary"] for l in labels]
    bars = ax.barh(labels, values, color=colors, height=0.6)
    for b, v in zip(bars, values):
        ax.text(b.get_width() + max(values) * 0.01, b.get_y() + b.get_height() / 2, f"{v:.1f}{unit}",
                va="center", fontsize=9)
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    ax.set_xlim(0, max(values) * 1.18)
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# ------------------------------------------------------------ 3 增长折线 ----

def growth_line(x: Sequence, series: dict[str, Sequence[float]], title="增长趋势", unit="", annotate_last=True,
                yoy_bars: Sequence[float] | None = None, source=None, figsize=(7, 4)):
    fig, ax = plt.subplots(figsize=figsize)
    pos = np.arange(len(x))  # 用序号做横轴，年份只作刻度标签，便于控制留白
    for i, (name, ys) in enumerate(series.items()):
        c = _P["series"][i % len(_P["series"])]
        ax.plot(pos, ys, marker="o", lw=2.2, ms=5, color=c, label=name)
        if annotate_last:
            ax.annotate(f"{ys[-1]:,.0f}{unit}", (pos[-1], ys[-1]), textcoords="offset points", xytext=(6, 0),
                        va="center", fontsize=9, color=c, fontweight="bold")
    ax.set_xticks(pos, [str(v) for v in x])
    ax.set_xlim(-0.4, len(x) - 1 + (1.1 if annotate_last else 0.4))
    if yoy_bars is not None:
        ax2 = ax.twinx()
        ax2.bar(pos, yoy_bars, alpha=0.18, color=_P["accent"], width=0.5, label="同比增速")
        ax2.set_ylim(0, max(yoy_bars) * 3)
        ax2.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
        ax2.grid(False)
        ax2.spines["top"].set_visible(False)
        ax2.legend(loc="upper left", bbox_to_anchor=(0.4, -0.08))
    ax.set_title(title)
    ax.set_ylabel(unit)
    ax.legend(loc="upper left", bbox_to_anchor=(0, -0.08), ncol=3)
    _source_note(ax, source)
    return fig


# -------------------------------------------------------------- 4 瀑布图 ----

def waterfall(labels: Sequence[str], deltas: Sequence[float], title="利润瀑布", start_label=None,
              total_label="净利润", unit="万元", source=None, figsize=(7.5, 4.2)):
    """labels/deltas 的第一项是起始值(如 收入)，其余为增减项，自动追加合计柱。"""
    vals = list(deltas)
    cum = np.cumsum(vals)
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(vals) + 1)
    bottoms, heights, colors = [], [], []
    for i, v in enumerate(vals):
        if i == 0:
            bottoms.append(0)
            heights.append(v)
            colors.append(_P["primary"])
        else:
            bottoms.append(cum[i - 1] if v >= 0 else cum[i])
            heights.append(abs(v))
            colors.append(_P["positive"] if v >= 0 else _P["negative"])
    bottoms.append(0)
    heights.append(cum[-1])
    colors.append(_P["accent"])
    ax.bar(x, heights, 0.6, bottom=bottoms, color=colors, edgecolor="none")
    for i in range(len(vals)):  # 连接虚线
        ax.plot([x[i] + 0.3, x[i + 1] - 0.3], [cum[i], cum[i]], color=_P["neutral"], lw=0.8, ls="--")
    for i, (b, h) in enumerate(zip(bottoms, heights)):
        v = vals[i] if i < len(vals) else cum[-1]
        ax.text(x[i], b + h + max(heights) * 0.015, f"{v:+,.0f}" if 0 < i < len(vals) else f"{v:,.0f}",
                ha="center", va="bottom", fontsize=8.5)
    ax.set_xticks(x, list(labels) + [total_label], rotation=0)
    ax.set_ylabel(unit)
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# --------------------------------------------------------------- 5 漏斗 ----

def funnel(stages: Sequence[str], values: Sequence[float], title="转化漏斗", source=None, figsize=(6.5, 4.2)):
    fig, ax = plt.subplots(figsize=figsize)
    vmax = max(values)
    n = len(values)
    for i, (s, v) in enumerate(zip(stages, values)):
        w = max(v / vmax, 0.22)  # 最窄也保留可读宽度，真实比例写在右侧文字里
        y = n - 1 - i
        ax.barh(y, w, left=(1 - w) / 2, height=0.78, color=_P["series"][i % len(_P["series"])])
        ax.text(0.5, y, s, ha="center", va="center", color="white", fontsize=10, fontweight="bold")
        conv = f"（上一步转化率 {v / values[i - 1]:.1%}）" if i > 0 else ""
        ax.text(0.5 + w / 2 + 0.02, y, f"{_fmt(v)} {conv}", ha="left", va="center", fontsize=8.5, color="#404040")
    ax.set_xlim(0, 1.45)
    ax.set_ylim(-0.6, n - 0.4)
    ax.axis("off")
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# ------------------------------------------------------- 6 敏感性龙卷风 ----

def tornado(factors: Sequence[str], low: Sequence[float], high: Sequence[float], base: float,
            title="敏感性分析（NPV）", unit="万元", low_label="-10%", high_label="+10%", source=None,
            figsize=(7, 4.2)):
    spread = np.abs(np.asarray(high) - np.asarray(low))
    order = np.argsort(spread)
    fig, ax = plt.subplots(figsize=figsize)
    y = np.arange(len(factors))
    for yi, idx in zip(y, order):
        lo, hi = low[idx] - base, high[idx] - base
        ax.barh(yi, lo, color=_P["negative"], height=0.6, label=low_label if yi == 0 else None)
        ax.barh(yi, hi, color=_P["secondary"], height=0.6, label=high_label if yi == 0 else None)
        pad = spread.max() * 0.02
        for delta, val in ((lo, low[idx]), (hi, high[idx])):  # 标签贴在各自柱子的外端
            ax.text(delta - pad if delta < 0 else delta + pad, yi, f"{val:,.0f}",
                    ha="right" if delta < 0 else "left", va="center", fontsize=8.5)
    ax.axvline(0, color="black", lw=1)
    ax.set_yticks(y, [factors[i] for i in order])
    ax.set_xlabel(f"相对基准 {base:,.0f}{unit} 的变化")
    ax.grid(axis="y", visible=False)
    ax.grid(axis="x", visible=True)
    lim = spread.max() * 0.8
    ax.set_xlim(-lim, lim)
    ax.legend(loc="lower right")
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# ------------------------------------------------------------ 7 盈亏平衡 ----

def breakeven(fixed_cost: float, price: float, var_cost: float, max_units: int, title="盈亏平衡分析",
              unit="万元", qty_unit="件", source=None, figsize=(7, 4.2)):
    q = np.linspace(0, max_units, 200)
    rev = q * price
    cost = fixed_cost + q * var_cost
    be_q = fixed_cost / (price - var_cost)
    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(q, rev, color=_P["secondary"], lw=2.2, label="总收入")
    ax.plot(q, cost, color=_P["negative"], lw=2.2, label="总成本")
    ax.axhline(fixed_cost, color=_P["neutral"], lw=1.2, ls="--", label="固定成本")
    ax.fill_between(q, rev, cost, where=rev >= cost, color=_P["positive"], alpha=0.12, label="盈利区")
    ax.fill_between(q, rev, cost, where=rev < cost, color=_P["negative"], alpha=0.10, label="亏损区")
    ax.scatter([be_q], [be_q * price], color=_P["accent"], zorder=5, s=60)
    ax.annotate(f"盈亏平衡点\n{be_q:,.0f}{qty_unit} / {be_q * price:,.0f}{unit}", (be_q, be_q * price),
                textcoords="offset points", xytext=(15, -40), fontsize=9,
                arrowprops=dict(arrowstyle="->", color=_P["accent"]))
    ax.set_xlabel(f"销量（{qty_unit}）")
    ax.set_ylabel(unit)
    ax.legend(loc="upper left", ncol=3, bbox_to_anchor=(0, -0.14))
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# ---------------------------------------------------- 8 现金流柱线组合 ----

def cashflow_combo(periods: Sequence, inflow: Sequence[float], outflow: Sequence[float], title="现金流预测",
                   unit="万元", source=None, figsize=(7.5, 4.2)):
    inflow, outflow = np.asarray(inflow, float), np.asarray(outflow, float)
    net = inflow - outflow
    cum = np.cumsum(net)
    x = np.arange(len(periods))
    fig, ax = plt.subplots(figsize=figsize)
    ax.bar(x - 0.2, inflow, 0.38, color=_P["secondary"], label="现金流入")
    ax.bar(x + 0.2, -outflow, 0.38, color=_P["light"], edgecolor=_P["neutral"], label="现金流出")
    ax.plot(x, cum, color=_P["accent"], marker="o", lw=2.2, label="累计净现金流")
    for xi, c in zip(x, cum):
        ax.annotate(f"{c:,.0f}", (xi, c), textcoords="offset points", xytext=(0, 10), ha="center", fontsize=8.5,
                    color=_P["accent"], fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(x, [str(p) for p in periods])
    ax.set_ylabel(unit)
    ax.legend(loc="upper left", ncol=3, bbox_to_anchor=(0, -0.08))
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# ------------------------------------------------------ 9 竞争矩阵散点 ----

def competition_matrix(names: Sequence[str], x: Sequence[float], y: Sequence[float], size: Sequence[float] | None = None,
                       highlight: str | None = None, xlabel="价格定位", ylabel="产品/服务能力",
                       quadrant_labels=("高价低能", "高价高能", "低价低能", "低价高能"), title="竞争格局矩阵",
                       source=None, figsize=(6.5, 5.5)):
    fig, ax = plt.subplots(figsize=figsize)
    size = np.asarray(size if size is not None else [300] * len(names), float)
    size = 200 + 1200 * (size - size.min()) / (np.ptp(size) or 1)
    for n, xi, yi, s in zip(names, x, y, size):
        is_h = n == highlight
        ax.scatter(xi, yi, s=s, color=_P["accent"] if is_h else _P["secondary"], alpha=0.85 if is_h else 0.55,
                   edgecolor="white", linewidth=1.5, zorder=3)
        ax.annotate(n, (xi, yi), textcoords="offset points", xytext=(0, 14 if not is_h else 18), ha="center",
                    fontsize=9, fontweight="bold" if is_h else "normal")
    xm, ym = (min(x) + max(x)) / 2, (min(y) + max(y)) / 2
    ax.axvline(xm, color=_P["neutral"], lw=0.8, ls="--")
    ax.axhline(ym, color=_P["neutral"], lw=0.8, ls="--")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    pad_x, pad_y = (max(x) - min(x)) * 0.25, (max(y) - min(y)) * 0.25
    ax.set_xlim(min(x) - pad_x, max(x) + pad_x)
    ax.set_ylim(min(y) - pad_y, max(y) + pad_y)
    for (lx, ly, txt) in [(0.02, 0.97, quadrant_labels[3]), (0.98, 0.97, quadrant_labels[1]),
                          (0.02, 0.03, quadrant_labels[2]), (0.98, 0.03, quadrant_labels[0])]:
        ax.text(lx, ly, txt, transform=ax.transAxes, fontsize=8, color=_P["neutral"],
                ha="left" if lx < 0.5 else "right", va="top" if ly > 0.5 else "bottom")
    ax.grid(False)
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# --------------------------------------------------------------- 10 甘特图 ----

def gantt(tasks: Sequence[tuple[str, int, int, str]], total_periods=12, period_label="月", title="实施计划",
          phase_colors: dict[str, str] | None = None, source=None, figsize=(8, 4)):
    """tasks: [(任务名, 开始期(1-based), 持续期数, 阶段名)]"""
    phases = list(dict.fromkeys(t[3] for t in tasks))
    pc = phase_colors or {p: _P["series"][i % len(_P["series"])] for i, p in enumerate(phases)}
    fig, ax = plt.subplots(figsize=figsize)
    for i, (name, start, dur, phase) in enumerate(tasks):
        ax.barh(len(tasks) - 1 - i, dur, left=start - 1, height=0.55, color=pc[phase])
        ax.text(start - 1 + dur / 2, len(tasks) - 1 - i, f"{dur}{period_label}", ha="center", va="center",
                color="white", fontsize=8)
    ax.set_yticks(range(len(tasks)), [t[0] for t in reversed(tasks)])
    ax.set_xticks(range(total_periods + 1), [f"M{i}" if i else "" for i in range(total_periods + 1)])
    ax.set_xlim(0, total_periods)
    ax.grid(axis="x", visible=True)
    ax.grid(axis="y", visible=False)
    for p in phases:
        ax.bar(0, 0, color=pc[p], label=p)
    ax.legend(loc="upper left", ncol=len(phases), bbox_to_anchor=(0, -0.08))
    ax.set_title(title)
    _source_note(ax, source)
    return fig


# ------------------------------------------------- 11 商业模式画布九宫格 ----

CANVAS_LAYOUT = {  # (x, y, w, h) 以 10x6 网格计
    "重要伙伴": (0, 2, 2, 4), "关键业务": (2, 4, 2, 2), "核心资源": (2, 2, 2, 2),
    "价值主张": (4, 2, 2, 4), "客户关系": (6, 4, 2, 2), "渠道通路": (6, 2, 2, 2),
    "客户细分": (8, 2, 2, 4), "成本结构": (0, 0, 5, 2), "收入来源": (5, 0, 5, 2),
}


def business_model_canvas(blocks: dict[str, Sequence[str]], title="商业模式画布", figsize=(12, 7.2)):
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    for name, (x, y, w, h) in CANVAS_LAYOUT.items():
        head_color = _P["primary"] if name not in ("成本结构", "收入来源") else _P["accent"]
        ax.add_patch(FancyBboxPatch((x + 0.04, y + 0.04), w - 0.08, h - 0.08, boxstyle="round,pad=0,rounding_size=0.06",
                                    fc="white", ec=_P["light"], lw=1.2))
        ax.add_patch(FancyBboxPatch((x + 0.04, y + h - 0.44), w - 0.08, 0.4, boxstyle="round,pad=0,rounding_size=0.06",
                                    fc=head_color, ec="none"))
        ax.text(x + 0.15, y + h - 0.24, name, color="white", fontsize=10.5, fontweight="bold", va="center")
        items = blocks.get(name, [])
        for i, it in enumerate(items[:6]):
            ax.text(x + 0.18, y + h - 0.75 - i * 0.42, f"• {it}", fontsize=8.8, va="top", wrap=True)
    ax.set_title(title, fontsize=14, loc="left", fontweight="bold")
    return fig
