"""slides.json → 16:9 路演 PPT（python-pptx 原生图表 / 表格 / 演讲者备注 / 匿名清洗），契约 §7.3。

用法：
    python -m bizkit.pptx_render 06_路演/slides.json

输入 slides.json（`image` 相对 slides.json 所在目录；`output` 相对工作区根）：
    { "title": "…", "subtitle": "…", "aspect": "16:9", "theme": "bizgold", "minutes": 10, "footer": "…",
      "anonymize": { "forbidden": ["XX大学"], "whitelist": [] },
      "slides": [ { "layout": "cover | agenda | section | bullets | image | two_column | chart | table | kpi |
                               timeline | team | quote | closing | appendix", …, "note": "讲稿", "duration": 40 } ],
      "output": "产出/路演.pptx" }

输出：
    <output>                 pptx（空白 16:9 母版 + 绝对坐标排版，主题色 / 主题字体已写入 theme1.xml）
    <output>.preview.json    数组，每页 {index, layout, title, texts[], image, chart, table, kpis, steps, members, quote, note}，
                             前端用它拼 HTML 预览；image 为相对工作区根的路径
    stdout JSON              {ok, output, preview, pages, warnings[], anonymizeHits[], perPage[{index, chars, minFont}]}

设计规则（渲染前校验；超出先警告，再截断或缩字号，保证不溢出）：
    标题 32pt，过长缩到 ≥28pt，再长则换两行；正文 ≤120 字、要点 ≤6 条、每条 ≤40 字、正文字号 ≥18pt；
    表格 ≤7 列 ≤8 行、字号 ≥14pt；饼图 ≤4 类；kpi 2–4 张；timeline 3–6 步；其余小字（图注 / 页脚 / 标签）≥14pt。
    `anonymize.forbidden` 命中的词替换为「XX」并记入 anonymizeHits（含备注）。
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import os
import re
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from lxml import etree
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_LEGEND_POSITION, XL_MARKER_STYLE, XL_TICK_LABEL_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, MSO_AUTO_SIZE, PP_ALIGN
from pptx.opc.constants import RELATIONSHIP_TYPE as RT
from pptx.oxml import parse_xml
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

from . import fonts, theme as theme_mod
from ._util import find_workspace_root, force_utf8_stdout, log, read_json, write_json

# ---------------------------------------------------------------------------------------------- 常量
SLIDE_W_IN, SLIDE_H_IN = 13.333, 7.5          # 16:9
MARGIN_IN = 0.6                                # 标题 / 页脚左右边距
CONTENT_LEFT_IN = 0.7
CONTENT_W_IN = SLIDE_W_IN - 2 * CONTENT_LEFT_IN
CONTENT_TOP_IN = 1.7
CONTENT_BOTTOM_IN = 6.75
FOOTER_TOP_IN = 6.98

TITLE_PT, TITLE_MIN_PT = 32.0, 28.0
BODY_START_PT, BODY_MIN_PT = 24.0, 18.0
SMALL_PT = 14.0                                # 页脚 / 图注 / 表格下限
TABLE_PT, TABLE_MIN_PT = 18.0, 14.0
KPI_PT, KPI_MIN_PT = 44.0, 28.0

LIMITS = {
    "titleChars": 28, "bodyChars": 120, "bullets": 6, "bulletChars": 40, "imageBullets": 3,
    "tableCols": 7, "tableRows": 8, "pieCats": 4, "kpiMin": 2, "kpiMax": 4, "stepsMin": 3, "stepsMax": 6,
    "members": 6, "agendaItems": 8,
}

LAYOUTS = (
    "cover", "agenda", "section", "bullets", "image", "two_column", "chart", "table", "kpi", "timeline",
    "team", "quote", "closing", "appendix",
)

# 图表类型别名 → python-pptx 枚举（bar 按 Excel 口径为横条；column_* 为竖柱）
CHART_TYPES: Dict[str, XL_CHART_TYPE] = {
    "column_clustered": XL_CHART_TYPE.COLUMN_CLUSTERED, "column": XL_CHART_TYPE.COLUMN_CLUSTERED,
    "grouped_bar": XL_CHART_TYPE.COLUMN_CLUSTERED,
    "bar": XL_CHART_TYPE.BAR_CLUSTERED, "hbar": XL_CHART_TYPE.BAR_CLUSTERED,
    "stacked": XL_CHART_TYPE.COLUMN_STACKED, "stacked_bar": XL_CHART_TYPE.COLUMN_STACKED,
    "column_stacked": XL_CHART_TYPE.COLUMN_STACKED,
    "line": XL_CHART_TYPE.LINE_MARKERS,
    "pie": XL_CHART_TYPE.PIE, "donut": XL_CHART_TYPE.DOUGHNUT, "doughnut": XL_CHART_TYPE.DOUGHNUT,
}

_CJK_RE = re.compile(r"[\u2e80-\u9fff\uf900-\ufaff\uff00-\uffef\u3000-\u303f]")
_NUM_RE = re.compile(r"^[+\-−]?[\d,.]+\s*(%|‰|x|X|倍|万|亿|元|万元|亿元|人|份|个|天|月|年)?$")
A_NS = "http://schemas.openxmlformats.org/drawingml/2006/main"
C_NS = "http://schemas.openxmlformats.org/drawingml/2006/chart"


def I(v: float) -> int:  # noqa: E743 - 短名字让坐标表更易读
    return int(Inches(v))


# ---------------------------------------------------------------------------------------------- 文字度量
def em_width(text: str) -> float:
    """以字号为单位估算文本宽度：CJK 一个 em，西文字母 / 数字约半个 em。只用于排版预估，不要求精确。"""
    w = 0.0
    for ch in text:
        if _CJK_RE.match(ch):
            w += 1.0
        elif ch == " ":
            w += 0.3
        elif ch.isdigit():
            w += 0.58
        elif ch.isupper():
            w += 0.72
        elif ch.isalpha():
            w += 0.54
        else:
            w += 0.42
    return w


def count_chars(text: str) -> int:
    """「字数」口径：去掉空白后的字符数（中英文、数字、标点一视同仁，便于模型自行核对）。"""
    return len(re.sub(r"\s+", "", text or ""))


def lines_needed(text: str, size_pt: float, width_pt: float) -> int:
    if not text:
        return 1
    return max(1, math.ceil(em_width(text) * size_pt / max(width_pt, 1.0) - 1e-6))


def block_height(paras: Sequence[str], size_pt: float, width_pt: float, line_spacing: float = 1.2, gap_pt: float = 6.0) -> float:
    lines = sum(lines_needed(p, size_pt, width_pt) for p in paras)
    return lines * size_pt * line_spacing + max(len(paras) - 1, 0) * gap_pt


def fit_size(paras: Sequence[str], width_pt: float, height_pt: float, start: float, minimum: float,
             line_spacing: float = 1.2, gap_pt: float = 6.0) -> Tuple[float, bool]:
    """从 start 逐磅缩小到 minimum，返回 (字号, 是否放得下)。"""
    size = start
    while size > minimum:
        if block_height(paras, size, width_pt, line_spacing, gap_pt) <= height_pt:
            return size, True
        size -= 1
    return minimum, block_height(paras, minimum, width_pt, line_spacing, gap_pt) <= height_pt


def truncate_em(text: str, max_em: float) -> str:
    """按估算宽度截断并加省略号。"""
    if em_width(text) <= max_em:
        return text
    out, w = "", 0.0
    for ch in text:
        cw = em_width(ch)
        if w + cw > max_em - 1.0:
            break
        out += ch
        w += cw
    return out.rstrip(" ，。、；：,.;:") + "…"


def shrink_paras(paras: Sequence[str], width_pt: float, height_pt: float, size_pt: float,
                 line_spacing: float = 1.2, gap_pt: float = 6.0) -> List[str]:
    """字号已到下限仍放不下时，把可用行数分给各段并截断超出部分（每段至少 1 行）。"""
    n = len(paras)
    if n == 0:
        return []
    line_h = size_pt * line_spacing
    total_lines = max(n, int((height_pt - (n - 1) * gap_pt) // line_h))
    needed = [lines_needed(p, size_pt, width_pt) for p in paras]
    allowed = [1] * n
    remaining = total_lines - n
    changed = True
    while remaining > 0 and changed:          # 轮流给仍不够的段落加一行，保证分配均匀
        changed = False
        for i in range(n):
            if allowed[i] < needed[i] and remaining > 0:
                allowed[i] += 1
                remaining -= 1
                changed = True
    out = []
    for p, a, nd in zip(paras, allowed, needed):
        out.append(truncate_em(p, a * width_pt / size_pt) if nd > a else p)
    return out


def balanced_split(text: str, size_pt: float, width_pt: float) -> Optional[str]:
    """标题需要两行时，在最靠近中点的标点后断行（两半都要放得下一行）；没有合适断点返回 None，交给自然换行。"""
    best: Optional[Tuple[float, str, str]] = None
    for i, ch in enumerate(text[:-1]):
        if ch not in "：:，,；;、。—–！!？? ":
            continue
        a, b = text[: i + 1].rstrip(), text[i + 1:].lstrip()
        if not a or not b or em_width(a) * size_pt > width_pt or em_width(b) * size_pt > width_pt:
            continue
        score = abs(em_width(a) - em_width(b))
        if best is None or score < best[0]:
            best = (score, a, b)
    return f"{best[1]}\n{best[2]}" if best else None


def one_line_size(text: str, width_pt: float, start: float, minimum: float) -> float:
    """让单行文本放进宽度内的最大字号（夹在 [minimum, start]）。"""
    return float(max(minimum, min(start, math.floor(width_pt / max(em_width(text), 1e-6)))))


def _luminance(hex_color: str) -> float:
    r, g, b = theme_mod.hex2rgb(hex_color)
    return (0.299 * r + 0.587 * g + 0.114 * b) / 255.0


def _is_numeric(text: str) -> bool:
    return bool(_NUM_RE.match(text.strip())) if text and text.strip() else False


# ---------------------------------------------------------------------------------------------- 匿名
def anonymize_value(value: Any, forbidden: Sequence[str], whitelist: Sequence[str], page: int, path: str,
                    hits: List[Dict[str, Any]]) -> Any:
    """递归替换 dict / list / str 中的禁用词为「XX」，白名单词先用私用区占位符保护。"""
    if isinstance(value, str):
        terms = [t for t in forbidden if t]
        if not terms:
            return value
        protected: Dict[str, str] = {}
        text = value
        for k, w in enumerate(w for w in whitelist if w):
            token = f"\ue000{k}\ue001"
            protected[token] = w
            text = text.replace(w, token)
        for t in sorted(terms, key=len, reverse=True):
            c = text.count(t)
            if c:
                hits.append({"page": page, "field": path, "term": t, "count": c})
                text = text.replace(t, "XX")
        for token, w in protected.items():
            text = text.replace(token, w)
        return text
    if isinstance(value, list):
        return [anonymize_value(v, forbidden, whitelist, page, f"{path}[{i}]", hits) for i, v in enumerate(value)]
    if isinstance(value, dict):
        return {k: anonymize_value(v, forbidden, whitelist, page, f"{path}.{k}" if path else str(k), hits)
                for k, v in value.items()}
    return value


# ---------------------------------------------------------------------------------------------- 统计
@dataclass
class PageStat:
    index: int
    layout: str
    chars: int = 0
    min_font: Optional[float] = None

    def font(self, size: float) -> None:
        self.min_font = size if self.min_font is None else min(self.min_font, size)


@dataclass
class DeckStats:
    warnings: List[str] = field(default_factory=list)
    anonymize_hits: List[Dict[str, Any]] = field(default_factory=list)
    missing_images: List[str] = field(default_factory=list)
    pages: List[PageStat] = field(default_factory=list)
    preview: List[Dict[str, Any]] = field(default_factory=list)


# ---------------------------------------------------------------------------------------------- 渲染器
class DeckRenderer:
    def __init__(self, cfg: Dict[str, Any], json_path: Path):
        self.cfg = cfg
        self.json_path = json_path
        self.base_dir = json_path.parent                       # image 相对这里
        self.root = find_workspace_root(json_path.parent)      # output 相对工作区根
        self.th = theme_mod.get_theme(str(cfg.get("theme") or "bizgold"))
        df = fonts.doc_fonts(latin="Arial")                    # PPT 用界面字体微软雅黑 + Arial，而非 Word 的宋体 / Times
        self.ea = df["ui"]
        self.latin = df["latin"]
        self.footer_text = str(cfg.get("footer") or cfg.get("title") or "")
        anon = cfg.get("anonymize") or {}
        self.forbidden = [str(t) for t in (anon.get("forbidden") or [])]
        self.whitelist = [str(t) for t in (anon.get("whitelist") or [])]
        self.stats = DeckStats()
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = I(SLIDE_W_IN), I(SLIDE_H_IN)
        self._patch_theme()
        self.total = 0
        self.section_no = 0
        self.cur = PageStat(0, "")
        if str(cfg.get("aspect") or "16:9") != "16:9":
            self.warn(0, f"aspect「{cfg.get('aspect')}」暂不支持，按 16:9 渲染")

    # ---- 基础 ----
    def warn(self, page: int, msg: str) -> None:
        self.stats.warnings.append((f"第 {page} 页：" if page else "") + msg)

    def _patch_theme(self) -> None:
        """python-pptx 没有主题 API：直接改 theme1.xml 的配色与字体，让未显式设字体的文本（含图表）也落到微软雅黑 / Arial。"""
        part = self.prs.slide_master.part.part_related_by(RT.THEME)
        root = etree.fromstring(part.blob)
        ns = {"a": A_NS}
        th = self.th
        colors = {"dk2": th.primary, "lt2": th.paper_alt, "hlink": th.secondary, "folHlink": th.muted}
        for i, c in enumerate(th.series[:6], 1):
            colors[f"accent{i}"] = c
        for name, hexv in colors.items():
            el = root.find(f".//a:clrScheme/a:{name}", ns)
            if el is None:
                continue
            for child in list(el):
                el.remove(child)
            etree.SubElement(el, qn("a:srgbClr")).set("val", hexv)
        for tag in ("majorFont", "minorFont"):
            font_el = root.find(f".//a:fontScheme/a:{tag}", ns)
            if font_el is None:
                continue
            latin = font_el.find("a:latin", ns)
            if latin is not None:
                latin.set("typeface", self.latin)
            ea = font_el.find("a:ea", ns)
            if ea is not None:
                ea.set("typeface", self.ea)
            for f in font_el.findall("a:font", ns):
                if f.get("script") == "Hans":
                    font_el.remove(f)
            hans = etree.SubElement(font_el, qn("a:font"))
            hans.set("script", "Hans")
            hans.set("typeface", self.ea)
        part.blob = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)

    def _set_ea(self, rpr) -> None:
        """font.name 只写 a:latin；中文字形由 a:ea 决定。a:rPr 子元素有固定顺序，ea 必须紧跟 latin。"""
        for old in rpr.findall(qn("a:ea")):
            rpr.remove(old)
        ea = etree.Element(qn("a:ea"))
        ea.set("typeface", self.ea)
        latin = rpr.find(qn("a:latin"))
        if latin is not None:
            latin.addnext(ea)
            return
        for tag in ("a:cs", "a:sym", "a:hlinkClick", "a:hlinkMouseOver", "a:rtl", "a:extLst"):
            nxt = rpr.find(qn(tag))
            if nxt is not None:
                nxt.addprevious(ea)
                return
        rpr.append(ea)

    def _style_run(self, run, size: float, bold: bool = False, color: str = theme_mod.INK, italic: bool = False,
                   track: bool = True) -> None:
        f = run.font
        f.size = Pt(size)
        f.bold = bold
        if italic:
            f.italic = True
        f.color.rgb = RGBColor.from_string(color)
        f.name = self.latin
        self._set_ea(run._r.get_or_add_rPr())
        if track:
            self.cur.font(size)

    def _text(self, slide, x: float, y: float, w: float, h: float, text: Any, size: float = 18.0, bold: bool = False,
              color: str = theme_mod.INK, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP, line_spacing: float = 1.15,
              autofit: bool = False, track: bool = True, italic: bool = False, gap_pt: float = 0.0):
        """文本框。text 为 str 或 [(段落文本, {size,bold,color,align,italic}), …]。"""
        tb = slide.shapes.add_textbox(I(x), I(y), I(w), I(h))
        tf = tb.text_frame
        tf.word_wrap = True
        tf.vertical_anchor = anchor
        if autofit:  # 预估失误时让 PowerPoint / LibreOffice 再缩一次，兜底不溢出
            tf.auto_size = MSO_AUTO_SIZE.TEXT_TO_FIT_SHAPE
        paras = text if isinstance(text, list) else [(str(text), {})]
        for i, (t, ov) in enumerate(paras):
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.alignment = ov.get("align", align)
            p.line_spacing = line_spacing
            if gap_pt and i > 0:
                p.space_before = Pt(gap_pt)
            # 段内 "\n" 写成 a:br（软回车），标题的平衡断行靠它
            for k, piece in enumerate(t.split("\n")):
                if k:
                    p._p.add_br()
                r = p.add_run()
                r.text = piece
                self._style_run(r, ov.get("size", size), ov.get("bold", bold), ov.get("color", color), ov.get("italic", italic), track)
        return tb

    def _bullets(self, slide, x: float, y: float, w: float, h: float, items: Sequence[str], size: float,
                 color: str = theme_mod.INK, bullet_color: Optional[str] = None, gap_pt: float = 6.0,
                 line_spacing: float = 1.2):
        """带真正项目符号（a:buChar）的要点框；python-pptx 文本框默认没有项目符号。"""
        tb = slide.shapes.add_textbox(I(x), I(y), I(w), I(h))
        tf = tb.text_frame
        tf.word_wrap = True
        tf.auto_size = MSO_AUTO_SIZE.TEXT_TO_FIT_SHAPE
        for i, text in enumerate(items):
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.line_spacing = line_spacing
            if i > 0:
                p.space_before = Pt(gap_pt)
            ppr = p._p.get_or_add_pPr()
            ppr.set("marL", str(int(Pt(size * 1.3))))
            ppr.set("indent", str(-int(Pt(size * 1.3))))
            for tag in ("a:buNone", "a:buChar", "a:buClr", "a:buFont", "a:buSzPct"):
                for old in ppr.findall(qn(tag)):
                    ppr.remove(old)
            bu_clr = etree.SubElement(ppr, qn("a:buClr"))
            etree.SubElement(bu_clr, qn("a:srgbClr")).set("val", bullet_color or self.th.accent)
            etree.SubElement(ppr, qn("a:buSzPct")).set("val", "70000")
            etree.SubElement(ppr, qn("a:buFont")).set("typeface", "Arial")
            etree.SubElement(ppr, qn("a:buChar")).set("char", "■")
            r = p.add_run()
            r.text = text
            self._style_run(r, size, False, color)
        return tb

    def _rect(self, slide, x: float, y: float, w: float, h: float, fill: str, shape=MSO_SHAPE.RECTANGLE,
              line: Optional[str] = None, line_pt: float = 1.0):
        shp = slide.shapes.add_shape(shape, I(x), I(y), I(w), I(h))
        shp.fill.solid()
        shp.fill.fore_color.rgb = RGBColor.from_string(fill)
        if line:
            shp.line.color.rgb = RGBColor.from_string(line)
            shp.line.width = Pt(line_pt)
        else:
            shp.line.fill.background()
        shp.shadow.inherit = False
        if shape == MSO_SHAPE.ROUNDED_RECTANGLE:
            shp.adjustments[0] = 0.06   # 默认圆角过大，显得不够商务
        shp.text_frame.text = ""
        return shp

    def _badge(self, slide, x: float, y: float, d: float, text: str, fill: str, size: float = 16.0,
               color: str = theme_mod.PAPER) -> None:
        c = self._rect(slide, x, y, d, d, fill, MSO_SHAPE.OVAL)
        tf = c.text_frame
        tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        r = p.add_run()
        r.text = text
        self._style_run(r, size, True, color)

    def _notes(self, slide, note: str, duration: Any) -> str:
        text = (note or "").strip()
        if duration:
            text = (text + "\n" if text else "") + f"（建议时长：{duration} 秒）"
        if not text:
            text = "（本页无讲稿）"
        tf = slide.notes_slide.notes_text_frame
        tf.text = text
        return text

    def _blank(self):
        return self.prs.slides.add_slide(self.prs.slide_layouts[6])   # 6 = Blank；自建页用绝对坐标

    def _rel_ws(self, p: Path) -> str:
        try:
            return Path(os.path.relpath(p, self.root)).as_posix()
        except ValueError:          # Windows 跨盘符
            return p.as_posix()

    # ---- 页框架 ----
    def _title_block(self, slide, title: str, x: float, y: float, w: float, h: float, color: str,
                     start: float = TITLE_PT, minimum: float = TITLE_MIN_PT, anchor=MSO_ANCHOR.MIDDLE,
                     align=PP_ALIGN.LEFT, prefer_two_lines: bool = False) -> float:
        """标题：默认 32pt，放不下一行缩到 28pt，再长换两行（优先在标点处平衡断行；超过两行截断并警告）。

        深底大标题（封面 / 章节 / 致谢）prefer_two_lines=True：放不下时先试大字号两行，再考虑缩字号，视觉更有力。
        """
        width_pt = w * 72 - 16
        em = em_width(title)
        size = float(start)
        text = title
        if em * size > width_pt:
            split: Optional[str] = None
            if prefer_two_lines:
                # 从大到小找第一个能在标点处平衡断成两行的字号（不低于 minimum）
                s = float(start)
                while s >= minimum and split is None:
                    split = balanced_split(title, s, width_pt)
                    if split is None:
                        s -= 2
                if split:
                    size = s
            if split is None:
                size = one_line_size(title, width_pt, start, minimum)
                if em * size > width_pt:                       # 下限字号仍放不下 → 两行
                    if prefer_two_lines:
                        size = max(minimum, min(float(start), math.floor(2 * width_pt / em)))
                    split = balanced_split(title, size, width_pt)
                    if em > 2 * width_pt / size:
                        # 自然换行会在行尾浪费一点空间，留 8% 余量以免省略号掉到第三行
                        text = truncate_em(title, 2 * width_pt / size * 0.92)
                        self.warn(self.cur.index, f"标题超过两行，已截断为「{text}」")
                        split = None
            if split:
                text = split
        self._text(slide, x, y, w, h, text, size=size, bold=True, color=color, anchor=anchor, align=align, line_spacing=1.1)
        return size

    def _frame(self, slide, title: str, appendix: bool = False) -> None:
        th = self.th
        self._rect(slide, 0, 0, SLIDE_W_IN, 0.1, th.primary)
        title_w = SLIDE_W_IN - 2 * MARGIN_IN - (1.7 if appendix else 0)
        self._title_block(slide, title, MARGIN_IN, 0.3, title_w, 1.05, th.primary)
        self._rect(slide, MARGIN_IN, 1.42, 1.3, 0.05, th.accent)
        if appendix:
            tag = self._rect(slide, SLIDE_W_IN - MARGIN_IN - 1.5, 0.3, 1.5, 0.42, th.accent, MSO_SHAPE.ROUNDED_RECTANGLE)
            tf = tag.text_frame
            tf.margin_top = tf.margin_bottom = 0
            tf.vertical_anchor = MSO_ANCHOR.MIDDLE
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.CENTER
            r = p.add_run()
            r.text = "备用页"
            self._style_run(r, SMALL_PT, True, th.paper, track=False)
        self._footer(slide, dark=False)

    def _footer(self, slide, dark: bool) -> None:
        th = self.th
        color = th.secondary_light if dark else th.muted
        if not dark:
            self._rect(slide, MARGIN_IN, 6.9, SLIDE_W_IN - 2 * MARGIN_IN, 0.012, th.line)
        self._text(slide, MARGIN_IN, FOOTER_TOP_IN, 8.0, 0.35, self.footer_text, size=SMALL_PT, color=color,
                   anchor=MSO_ANCHOR.MIDDLE, track=False)
        self._text(slide, SLIDE_W_IN - MARGIN_IN - 2.0, FOOTER_TOP_IN, 2.0, 0.35, f"{self.cur.index} / {self.total}",
                   size=SMALL_PT, color=color, align=PP_ALIGN.RIGHT, anchor=MSO_ANCHOR.MIDDLE, track=False)

    # ---- 正文块（校验 + 缩放 + 截断）----
    def _fit_block(self, paras: Sequence[str], w_in: float, h_in: float, start: float, minimum: float,
                   line_spacing: float = 1.2, gap_pt: float = 6.0, what: str = "正文") -> Tuple[List[str], float]:
        paras = [p for p in paras]
        width_pt = w_in * 72 - 16 - (start * 1.3 if what == "要点" else 0)
        height_pt = h_in * 72 - 8
        size, ok = fit_size(paras, width_pt, height_pt, start, minimum, line_spacing, gap_pt)
        if not ok:
            paras = shrink_paras(paras, width_pt, height_pt, size, line_spacing, gap_pt)
            self.warn(self.cur.index, f"{what}放不下，字号已缩至 {size:g}pt 并截断超出部分（请精简文字）")
        elif size < start and size < 20 and what in ("要点", "正文"):
            self.warn(self.cur.index, f"{what}字号缩至 {size:g}pt（建议精简文字到 20pt 以上）")
        return paras, size

    def _body_chars(self, texts: Sequence[str], check: bool = True) -> int:
        n = sum(count_chars(t) for t in texts)
        self.cur.chars += n
        if check and n > LIMITS["bodyChars"]:
            self.warn(self.cur.index, f"正文 {n} 字 > {LIMITS['bodyChars']}")
        return n

    def _norm_bullets(self, raw: Any, limit: int, what: str = "要点") -> List[str]:
        items = [str(b).strip() for b in (raw or []) if str(b).strip()] if isinstance(raw, list) else ([str(raw)] if raw else [])
        if len(items) > limit:
            self.warn(self.cur.index, f"{what} {len(items)} 条 > {limit}，已保留前 {limit} 条")
            items = items[:limit]
        for k, b in enumerate(items, 1):
            if count_chars(b) > LIMITS["bulletChars"]:
                self.warn(self.cur.index, f"第 {k} 条{what} {count_chars(b)} 字 > {LIMITS['bulletChars']}")
        return items

    # ---- 图片 ----
    def _resolve_image(self, rel: str) -> Optional[Path]:
        p = Path(rel)
        cands = [p] if p.is_absolute() else [self.base_dir / p, self.root / p]
        for c in cands:
            if c.exists():
                return c.resolve()
        return None

    def _picture(self, slide, rel: str, x: float, y: float, w: float, h: float) -> Tuple[Optional[Path], Tuple[float, float, float, float]]:
        """按比例放进框内并居中；缺图画占位框并警告。返回 (绝对路径, 实际几何 (x,y,w,h) 英寸)。"""
        path = self._resolve_image(rel) if rel else None
        if path is None:
            self.stats.missing_images.append(rel)
            self.warn(self.cur.index, f"找不到图片：{rel}")
            box = self._rect(slide, x, y, w, h, self.th.paper_alt, MSO_SHAPE.RECTANGLE, line=self.th.line)
            tf = box.text_frame
            tf.vertical_anchor = MSO_ANCHOR.MIDDLE
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.CENTER
            r = p.add_run()
            r.text = f"【缺图：{rel or '未指定 image'}】"
            self._style_run(r, 18, False, self.th.negative, track=False)
            return None, (x, y, w, h)
        pic = slide.shapes.add_picture(str(path), I(x), I(y))
        scale = min(I(w) / pic.width, I(h) / pic.height)
        pw, ph = int(pic.width * scale), int(pic.height * scale)
        pic.width, pic.height = pw, ph
        pic.left = I(x) + (I(w) - pw) // 2
        pic.top = I(y) + (I(h) - ph) // 2
        return path, (pic.left / 914400, pic.top / 914400, pw / 914400, ph / 914400)

    # ---- 原生图表 ----
    def _norm_chart(self, raw: Any) -> Optional[Dict[str, Any]]:
        if not isinstance(raw, dict):
            self.warn(self.cur.index, "chart 字段缺失或不是对象")
            return None
        ctype = str(raw.get("type") or "column_clustered").strip().lower()
        if ctype not in CHART_TYPES:
            self.warn(self.cur.index, f"不支持的图表类型「{ctype}」，已改为 column_clustered（可选：{', '.join(sorted(set(CHART_TYPES)))}）")
            ctype = "column_clustered"
        cats = [str(c) for c in (raw.get("categories") or [])]
        series_in = raw.get("series") or []
        if not cats or not isinstance(series_in, list) or not series_in:
            self.warn(self.cur.index, "chart 需要 categories[] 与 series[{name, values}]")
            return None
        series: List[Dict[str, Any]] = []
        for k, s in enumerate(series_in):
            if not isinstance(s, dict):
                continue
            vals = list(s.get("values") or [])
            if len(vals) != len(cats):
                self.warn(self.cur.index, f"系列「{s.get('name', k + 1)}」有 {len(vals)} 个值，类目有 {len(cats)} 个，已对齐")
                vals = (vals + [None] * len(cats))[:len(cats)]
            clean = []
            for v in vals:
                try:
                    clean.append(None if v is None or v == "" else float(v))
                except (TypeError, ValueError):
                    clean.append(None)
            series.append({"name": str(s.get("name") or f"系列{k + 1}"), "values": clean})
        if not series:
            self.warn(self.cur.index, "chart.series 为空")
            return None
        if ctype in ("pie", "donut", "doughnut"):
            if len(cats) > LIMITS["pieCats"]:
                self.warn(self.cur.index, f"饼图有 {len(cats)} 类 > {LIMITS['pieCats']}，已改为 column_clustered")
                ctype = "column_clustered"
            elif len(series) > 1:
                self.warn(self.cur.index, "饼图只取第一个系列")
                series = series[:1]
        decimals = raw.get("decimals")
        if decimals is None:
            decimals = 0
            for s in series:
                for v in s["values"]:
                    if v is not None and abs(v - round(v)) > 1e-9:
                        decimals = max(decimals, min(2, len(f"{v:.2f}".rstrip("0").split(".")[-1])))
        return {"type": ctype, "categories": cats, "series": series, "unit": str(raw.get("unit") or ""), "decimals": int(decimals)}

    def _chart_ea_fonts(self, chart) -> None:
        """chart.font 只产生 a:latin；轴标签 / 图例 / 数据标签里的中文还要补 a:ea。"""
        for rpr in chart._chartSpace.iter(qn("a:defRPr"), qn("a:rPr")):
            self._set_ea(rpr)

    @staticmethod
    def _chart_no_frame(chart) -> None:
        """去掉图表外框与底色（python-pptx 的 Chart 没有 format 属性，c:spPr 需手工插在 c:chart 之后）。"""
        cs = chart._chartSpace
        for old in cs.findall(qn("c:spPr")):
            cs.remove(old)
        sppr = parse_xml(
            f'<c:spPr xmlns:c="{C_NS}" xmlns:a="{A_NS}"><a:noFill/><a:ln><a:noFill/></a:ln></c:spPr>'
        )
        cs.chart.addnext(sppr)

    def _add_chart(self, slide, spec: Dict[str, Any], x: float, y: float, w: float, h: float):
        th = self.th
        ctype = spec["type"]
        xl_type = CHART_TYPES[ctype]
        is_hbar = xl_type == XL_CHART_TYPE.BAR_CLUSTERED
        cd = CategoryChartData()
        # 横条图自下而上画类目：把数据倒序喂入，第一类自然出现在最上面，且不用动坐标轴交叉设置
        # （catAx.crosses / reverse_order 在 LibreOffice 与 PowerPoint 表现不一致）
        order = slice(None, None, -1) if is_hbar else slice(None)
        cd.categories = spec["categories"][order]
        for s in spec["series"]:
            cd.add_series(s["name"], s["values"][order])
        chart = slide.shapes.add_chart(xl_type, I(x), I(y), I(w), I(h), cd).chart
        chart.font.size = Pt(SMALL_PT)
        chart.font.name = self.latin
        chart.font.color.rgb = RGBColor.from_string(th.ink)
        chart.has_title = False          # 单系列时 PowerPoint 会自动把系列名当图表标题，页面标题已承担此职
        self._chart_no_frame(chart)
        self.cur.font(SMALL_PT)
        nf = "#,##0" + ("." + "0" * spec["decimals"] if spec["decimals"] else "")
        plot = chart.plots[0]
        is_pie = xl_type in (XL_CHART_TYPE.PIE, XL_CHART_TYPE.DOUGHNUT)
        is_line = xl_type == XL_CHART_TYPE.LINE_MARKERS
        is_stacked = xl_type == XL_CHART_TYPE.COLUMN_STACKED
        if is_pie:
            plot.vary_by_categories = True
            for i, pt in enumerate(plot.series[0].points):
                pt.format.fill.solid()
                pt.format.fill.fore_color.rgb = RGBColor.from_string(th.series[i % len(th.series)])
            plot.has_data_labels = True
            dl = plot.data_labels
            dl.show_percentage = True
            dl.show_value = False
            dl.show_category_name = False
            dl.number_format = "0.0%"
            dl.number_format_is_linked = False
            dl.position = XL_LABEL_POSITION.OUTSIDE_END
            dl.font.size = Pt(SMALL_PT)
            dl.font.color.rgb = RGBColor.from_string(th.ink)
            chart.has_legend = True
        else:
            for i, s in enumerate(plot.series):
                c = RGBColor.from_string(th.series[i % len(th.series)])
                if is_line:
                    s.format.line.color.rgb = c
                    s.format.line.width = Pt(2.5)
                    s.smooth = False
                    s.marker.style = XL_MARKER_STYLE.CIRCLE
                    s.marker.size = 8
                    s.marker.format.fill.solid()
                    s.marker.format.fill.fore_color.rgb = c
                    s.marker.format.line.color.rgb = c
                else:
                    s.format.fill.solid()
                    s.format.fill.fore_color.rgb = c
                    s.invert_if_negative = False
            if not is_line:
                plot.gap_width = 60 if is_stacked else 80
                if is_stacked:
                    plot.overlap = 100
            if is_stacked:
                # 堆叠段落里的标签：深色段用白字、浅色段用深字，系列级 dLbls 会覆盖图级设置，所以逐系列写全
                for i, s in enumerate(plot.series):
                    sdl = s.data_labels
                    sdl.show_value = True
                    sdl.number_format = nf
                    sdl.number_format_is_linked = False
                    sdl.position = XL_LABEL_POSITION.CENTER
                    sdl.font.size = Pt(SMALL_PT)
                    dark = _luminance(th.series[i % len(th.series)]) < 0.55
                    sdl.font.color.rgb = RGBColor.from_string(th.paper if dark else th.ink)
            else:
                plot.has_data_labels = True
                dl = plot.data_labels
                dl.show_value = True
                dl.number_format = nf
                dl.number_format_is_linked = False
                dl.position = XL_LABEL_POSITION.ABOVE if is_line else XL_LABEL_POSITION.OUTSIDE_END
                dl.font.size = Pt(SMALL_PT)
                dl.font.color.rgb = RGBColor.from_string(th.ink)
            chart.has_legend = len(spec["series"]) > 1
            va, ca = chart.value_axis, chart.category_axis
            va.has_major_gridlines = True
            va.major_gridlines.format.line.color.rgb = RGBColor.from_string(th.line)
            va.format.line.fill.background()
            va.tick_labels.font.size = Pt(SMALL_PT)
            va.tick_labels.font.color.rgb = RGBColor.from_string(th.muted)
            va.tick_labels.number_format = nf
            va.tick_labels.number_format_is_linked = False
            ca.tick_labels.font.size = Pt(SMALL_PT)
            ca.tick_labels.font.color.rgb = RGBColor.from_string(th.ink)
            ca.format.line.color.rgb = RGBColor.from_string(th.secondary_light)
            ca.has_major_gridlines = False
            # 类目标签固定贴在绘图区低端，否则有负值时会压在零线上与负数标签重叠
            ca.tick_label_position = XL_TICK_LABEL_POSITION.LOW
            if spec["unit"]:
                va.has_title = True
                va.axis_title.text_frame.text = spec["unit"]
                r = va.axis_title.text_frame.paragraphs[0].runs[0]
                self._style_run(r, SMALL_PT, False, th.muted)
        if chart.has_legend:
            chart.legend.position = XL_LEGEND_POSITION.BOTTOM
            chart.legend.include_in_layout = False
            chart.legend.font.size = Pt(SMALL_PT)
            chart.legend.font.color.rgb = RGBColor.from_string(th.ink)
        self._chart_ea_fonts(chart)
        return chart

    # ---- 表格 ----
    def _norm_table(self, slide_cfg: Dict[str, Any]) -> Optional[Tuple[List[str], List[List[str]]]]:
        cols = [str(c) for c in (slide_cfg.get("columns") or [])]
        rows_in = slide_cfg.get("rows") or []
        if not cols or not isinstance(rows_in, list):
            self.warn(self.cur.index, "table 需要 columns[] 与 rows[][]")
            return None
        n_cols_in = len(cols)
        if len(cols) > LIMITS["tableCols"]:
            self.warn(self.cur.index, f"表格 {len(cols)} 列 > {LIMITS['tableCols']}，已截去多余列")
            cols = cols[:LIMITS["tableCols"]]
        if len(rows_in) > LIMITS["tableRows"]:
            self.warn(self.cur.index, f"表格 {len(rows_in)} 行 > {LIMITS['tableRows']}，已截去多余行")
            rows_in = rows_in[:LIMITS["tableRows"]]
        rows: List[List[str]] = []
        ragged = 0
        for r in rows_in:
            cells = [("" if c is None else str(c)) for c in (r if isinstance(r, list) else [r])]
            if len(cells) != n_cols_in:
                ragged += 1
            rows.append((cells + [""] * len(cols))[:len(cols)])
        if ragged:
            self.warn(self.cur.index, f"表格有 {ragged} 行的格数与列数 {n_cols_in} 不一致，已补空 / 截断")
        return cols, rows

    def _add_table(self, slide, cols: List[str], rows: List[List[str]], x: float, y: float, w: float, h: float) -> None:
        th = self.th
        n_rows, n_cols = len(rows) + 1, len(cols)
        # 自动列宽：按各列最长内容的估算宽度分配，夹在 [0.9in, 4.5in]，再按比例缩放到总宽
        col_em = [max(2.0, max(em_width(c) for c in [cols[j]] + [r[j] for r in rows])) for j in range(n_cols)]
        weights = [min(max(e, 4.0), 22.0) for e in col_em]
        widths = [w * wt / sum(weights) for wt in weights]
        size = TABLE_PT
        while True:
            row_h = []
            for r in [cols] + rows:
                lines = max(lines_needed(c, size, cw * 72 - 14) for c, cw in zip(r, widths))
                row_h.append(max(0.45, lines * size * 1.2 / 72 + 0.16))
            if sum(row_h) <= h or size <= TABLE_MIN_PT:
                break
            size -= 1
        if sum(row_h) > h:
            self.warn(self.cur.index, f"表格内容过多，{TABLE_MIN_PT:g}pt 仍可能溢出（请删行或精简单元格文字）")
            scale = h / sum(row_h)
            row_h = [rh * scale for rh in row_h]
        else:
            # 行数少时把行高拉开（每行最多 0.75in），避免表格缩在页面上方留一大片空白
            target = min(h, 0.75 * n_rows)
            if sum(row_h) < target:
                scale = target / sum(row_h)
                row_h = [rh * scale for rh in row_h]
        gt = slide.shapes.add_table(n_rows, n_cols, I(x), I(y), I(w), I(sum(row_h)))
        table = gt.table
        tbl_pr = gt._element.graphic.graphicData.tbl.tblPr
        tbl_pr.set("firstRow", "1")
        tbl_pr.set("bandRow", "0")      # 关掉内置条纹，否则自定义底色被盖住
        for j, cw in enumerate(widths):
            table.columns[j].width = I(cw)
        for i, rh in enumerate(row_h):
            table.rows[i].height = I(rh)
        numeric_col = [all(_is_numeric(r[j]) for r in rows if r[j].strip()) and any(r[j].strip() for r in rows) for j in range(n_cols)]
        for i, r in enumerate([cols] + rows):
            for j, val in enumerate(r):
                cell = table.cell(i, j)
                cell.margin_left = cell.margin_right = I(0.08)
                cell.margin_top = cell.margin_bottom = I(0.04)
                cell.vertical_anchor = MSO_ANCHOR.MIDDLE
                cell.fill.solid()
                header = i == 0
                cell.fill.fore_color.rgb = RGBColor.from_string(th.primary if header else (th.paper_alt if i % 2 == 0 else th.paper))
                tf = cell.text_frame
                tf.word_wrap = True
                p = tf.paragraphs[0]
                if header:
                    p.alignment = PP_ALIGN.CENTER
                elif j == 0:
                    p.alignment = PP_ALIGN.LEFT
                else:
                    p.alignment = PP_ALIGN.RIGHT if numeric_col[j] else PP_ALIGN.CENTER
                run = p.add_run()
                run.text = val
                self._style_run(run, size, header or j == 0, th.paper if header else th.ink)

    # ---- 卡片 ----
    def _kpi_card(self, slide, x: float, y: float, w: float, h: float, value: str, label: str, accent: str) -> None:
        th = self.th
        self._rect(slide, x, y, w, h, th.paper_alt, MSO_SHAPE.ROUNDED_RECTANGLE)
        self._rect(slide, x + 0.25, y, w - 0.5, 0.07, accent)
        # 大数字必须一行放下：按框内净宽直接算字号（粗体略宽，留 8% 余量），夹在 [28, 44]，仍放不下才截断
        inner_w = w - 0.3
        avail_pt = (inner_w * 72 - 16) * 0.92
        vsize = one_line_size(value, avail_pt, KPI_PT, KPI_MIN_PT)
        if em_width(value) * vsize > avail_pt:
            value = truncate_em(value, avail_pt / vsize)
            self.warn(self.cur.index, f"数字「{value}」过长，已截断")
        label_h = h - 1.45
        labels, lsize = self._fit_block([label], inner_w, label_h, 18, 14, 1.15, 0, "标签")
        self._text(slide, x + 0.15, y + 0.35, inner_w, 1.1, value, size=vsize, bold=True, color=th.primary,
                   align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE, line_spacing=1.0)
        self._text(slide, x + 0.15, y + 1.4, inner_w, label_h, labels[0], size=lsize, color=th.ink,
                   align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.TOP)

    # ================================================================================ 页型
    def _render_slide(self, raw: Dict[str, Any]) -> None:
        idx = len(self.prs.slides) + 1
        layout = str(raw.get("layout") or "").strip()
        self.cur = PageStat(idx, layout)
        self.stats.pages.append(self.cur)
        hits: List[Dict[str, Any]] = []
        s = anonymize_value(copy.deepcopy(raw), self.forbidden, self.whitelist, idx, "", hits)
        self.stats.anonymize_hits.extend(hits)
        if layout not in LAYOUTS:
            self.warn(idx, f"未知 layout「{layout}」，按 bullets 渲染（可选：{', '.join(LAYOUTS)}）")
            layout = "bullets"
            self.cur.layout = layout
        title = str(s.get("title") or "").strip()
        if not title and layout not in ("cover", "closing"):
            self.warn(idx, f"{layout} 页缺少 title")
        if title and count_chars(title) > LIMITS["titleChars"]:
            self.warn(idx, f"标题 {count_chars(title)} 字 > {LIMITS['titleChars']}（标题应是一句结论）")
        slide = self._blank()
        preview: Dict[str, Any] = {
            "index": idx, "layout": layout, "title": title, "texts": [], "image": None, "chart": None, "table": None,
            "kpis": None, "steps": None, "members": None, "quote": None, "note": "",
        }
        builder = getattr(self, f"_slide_{layout}")
        builder(slide, s, title, preview)
        preview["note"] = self._notes(slide, str(s.get("note") or ""), s.get("duration"))
        self.stats.preview.append(preview)

    # ---- 深底页：cover / section / closing ----
    def _dark_bg(self, slide) -> None:
        self._rect(slide, 0, 0, SLIDE_W_IN, SLIDE_H_IN, self.th.primary)
        # 右下角一块更深的色带 + 金线，避免整页纯色过于单调
        self._rect(slide, 9.2, 0, SLIDE_W_IN - 9.2, SLIDE_H_IN, self.th.primary_dark)
        self._rect(slide, 9.2, 0, 0.05, SLIDE_H_IN, self.th.accent)

    def _slide_cover(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._dark_bg(slide)
        deck_title = str(self.cfg.get("title") or "")
        deck_sub = str(self.cfg.get("subtitle") or "")
        subtitle = str(s.get("subtitle") or "")
        if deck_title:
            self._text(slide, 0.9, 0.55, 8.0, 0.5, deck_title, size=16, color=th.accent_light, anchor=MSO_ANCHOR.MIDDLE)
        self._rect(slide, 0.9, 2.35, 0.08, 1.9, th.accent)
        self._title_block(slide, title or deck_title, 1.2, 2.2, 7.7, 1.9, th.paper, start=44, minimum=32,
                          anchor=MSO_ANCHOR.MIDDLE, prefer_two_lines=True)
        if subtitle:
            subs, ssize = self._fit_block([subtitle], 7.7, 0.9, 24, 18, 1.15, 0, "副标题")
            self._text(slide, 1.2, 4.25, 7.7, 0.9, subs[0], size=ssize, color=th.secondary_light, anchor=MSO_ANCHOR.TOP)
        self._rect(slide, 0.9, 5.55, 7.9, 0.03, th.accent)
        if deck_sub:
            self._text(slide, 0.9, 5.7, 7.9, 0.6, deck_sub, size=18, color=th.paper, anchor=MSO_ANCHOR.MIDDLE)
        minutes = self.cfg.get("minutes")
        if minutes:
            self._text(slide, 9.6, 5.7, 3.2, 0.6, f"路演 {minutes} 分钟", size=SMALL_PT, color=th.secondary_light,
                       align=PP_ALIGN.RIGHT, anchor=MSO_ANCHOR.MIDDLE, track=False)
        self._body_chars([subtitle], check=False)
        pv["texts"] = [t for t in (subtitle, deck_sub) if t]
        self._footer(slide, dark=True)

    def _slide_section(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self.section_no += 1
        self._dark_bg(slide)
        kicker = str(s.get("kicker") or "")
        self._text(slide, 0.9, 1.4, 4.0, 1.4, f"{self.section_no:02d}", size=72, bold=True, color=th.accent, anchor=MSO_ANCHOR.MIDDLE)
        self._rect(slide, 0.9, 2.95, 1.6, 0.06, th.accent)
        self._title_block(slide, title, 0.9, 3.1, 8.0, 1.5, th.paper, start=40, minimum=30, prefer_two_lines=True)
        if kicker:
            ks, ksize = self._fit_block([kicker], 8.0, 1.0, 22, 18, 1.15, 0, "副标题")
            self._text(slide, 0.9, 4.65, 8.0, 1.0, ks[0], size=ksize, color=th.accent_light)
        self._body_chars([kicker])
        pv["texts"] = [kicker] if kicker else []
        self._footer(slide, dark=True)

    def _slide_closing(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._dark_bg(slide)
        subtitle = str(s.get("subtitle") or "")
        contact = str(s.get("contact") or "")
        self._title_block(slide, title or "谢谢", 0.9, 2.3, 7.9, 1.5, th.paper, start=44, minimum=32, prefer_two_lines=True)
        self._rect(slide, 0.9, 3.95, 1.6, 0.06, th.accent)
        texts = []
        if subtitle:
            ss, ssize = self._fit_block([subtitle], 7.9, 0.9, 22, 18, 1.15, 0, "副标题")
            self._text(slide, 0.9, 4.15, 7.9, 0.9, ss[0], size=ssize, color=th.secondary_light)
            texts.append(subtitle)
        if contact:
            self._text(slide, 0.9, 5.2, 7.9, 0.6, contact, size=16, color=th.accent_light)
            texts.append(contact)
        deck_title = str(self.cfg.get("title") or "")
        if deck_title:
            self._text(slide, 0.9, 0.55, 8.0, 0.5, deck_title, size=16, color=th.accent_light, anchor=MSO_ANCHOR.MIDDLE)
        self._body_chars(texts, check=False)
        pv["texts"] = texts
        self._footer(slide, dark=True)

    # ---- agenda ----
    def _slide_agenda(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._frame(slide, title or "目录")
        items = self._norm_bullets(s.get("items"), LIMITS["agendaItems"], "目录项")
        self._body_chars(items)
        # 右侧深色面板放演示标题，左侧编号列表；> 4 项分两列
        panel_x = 9.3
        self._rect(slide, panel_x, CONTENT_TOP_IN, SLIDE_W_IN - MARGIN_IN - panel_x, CONTENT_BOTTOM_IN - CONTENT_TOP_IN, th.primary, MSO_SHAPE.ROUNDED_RECTANGLE)
        self._rect(slide, panel_x + 0.35, 2.3, 1.0, 0.05, th.accent)
        deck_title = str(self.cfg.get("title") or "")
        deck_sub = str(self.cfg.get("subtitle") or "")
        pw = SLIDE_W_IN - MARGIN_IN - panel_x - 0.7
        dsize = one_line_size(deck_title, pw * 72 - 16, 24, 18)      # 尽量一行，放不下再自然换行
        dts, dsize = self._fit_block([deck_title], pw, 2.2, dsize, 18, 1.15, 0, "副标题")
        self._text(slide, panel_x + 0.35, 2.5, pw, 2.2, dts[0], size=dsize, bold=True, color=th.paper)
        if deck_sub:
            self._text(slide, panel_x + 0.35, 4.9, pw, 1.4, deck_sub, size=16, color=th.accent_light)
        cols = 2 if len(items) > 4 else 1
        rows = math.ceil(len(items) / cols) if items else 0
        avail_h = CONTENT_BOTTOM_IN - CONTENT_TOP_IN - 0.2
        row_h = min(1.0, avail_h / max(rows, 1))
        col_w = (panel_x - CONTENT_LEFT_IN - 0.3) / cols
        y0 = CONTENT_TOP_IN + 0.1 + (avail_h - row_h * rows) / 2
        text_w = col_w - 0.8
        for i, it in enumerate(items):
            c, r = divmod(i, rows) if rows else (0, 0)
            x = CONTENT_LEFT_IN + c * col_w
            y = y0 + r * row_h
            self._badge(slide, x, y + (row_h - 0.55) / 2, 0.55, str(i + 1), th.accent, 16)
            its, isize = self._fit_block([it], text_w, row_h, 22, 18, 1.1, 0, "目录项")
            self._text(slide, x + 0.7, y, text_w, row_h, its[0], size=isize, color=th.ink, anchor=MSO_ANCHOR.MIDDLE)
        pv["texts"] = items

    # ---- bullets ----
    def _slide_bullets(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        self._frame(slide, title)
        self._bullets_body(slide, s, pv, CONTENT_TOP_IN + 0.1)

    def _bullets_body(self, slide, s: Dict[str, Any], pv: Dict[str, Any], top: float) -> None:
        th = self.th
        items = self._norm_bullets(s.get("bullets"), LIMITS["bullets"])
        highlight = str(s.get("highlight") or "").strip()
        label = str(s.get("highlightLabel") or "").strip()
        self._body_chars(items + ([highlight, label] if highlight else []))
        h = CONTENT_BOTTOM_IN - top
        if highlight:
            card_w = 3.9
            text_w = CONTENT_W_IN - card_w - 0.4
            cx = CONTENT_LEFT_IN + text_w + 0.4
            self._kpi_card(slide, cx, top + 0.25, card_w, 2.7, highlight, label, th.accent)
        else:
            text_w = CONTENT_W_IN
        if items:
            items, size = self._fit_block(items, text_w, h, BODY_START_PT, BODY_MIN_PT, 1.2, 8.0, "要点")
            self._bullets(slide, CONTENT_LEFT_IN, top, text_w, h, items, size, gap_pt=8.0)
        pv["texts"] = list(items)
        if highlight:
            pv["kpis"] = [{"value": highlight, "label": label}]

    # ---- image ----
    def _slide_image(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        self._frame(slide, title)
        self._image_body(slide, s, pv, CONTENT_TOP_IN)

    def _image_body(self, slide, s: Dict[str, Any], pv: Dict[str, Any], top: float) -> None:
        th = self.th
        items = self._norm_bullets(s.get("bullets"), LIMITS["imageBullets"], "读图要点")
        caption = str(s.get("caption") or "").strip()
        self._body_chars(items + [caption])
        h = CONTENT_BOTTOM_IN - top
        cap_h = 0.45 if caption else 0.0
        if items:
            img_w = 7.3
            img_box = (CONTENT_LEFT_IN, top + 0.05, img_w, h - cap_h - 0.1)
            text_x = CONTENT_LEFT_IN + img_w + 0.4
            text_w = CONTENT_W_IN - img_w - 0.4
            items, size = self._fit_block(items, text_w, h - 0.3, 22, BODY_MIN_PT, 1.2, 10.0, "要点")
            self._bullets(slide, text_x, top + 0.2, text_w, h - 0.3, items, size, gap_pt=10.0)
        else:
            img_box = (CONTENT_LEFT_IN, top + 0.05, CONTENT_W_IN, h - cap_h - 0.1)
        path, geom = self._picture(slide, str(s.get("image") or ""), *img_box)
        if caption:
            cx, cy, cw, ch = geom
            self._text(slide, cx, cy + ch + 0.05, cw, 0.4, caption, size=SMALL_PT, color=th.muted, align=PP_ALIGN.CENTER,
                       anchor=MSO_ANCHOR.TOP)
        pv["image"] = self._rel_ws(path) if path else None
        pv["texts"] = list(items) + ([caption] if caption else [])

    # ---- two_column ----
    def _slide_two_column(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._frame(slide, title)
        texts: List[str] = []
        gap = 0.3
        col_w = (CONTENT_W_IN - gap) / 2
        top = CONTENT_TOP_IN
        h = CONTENT_BOTTOM_IN - top
        for k, side in enumerate(("left", "right")):
            col = s.get(side) or {}
            if not isinstance(col, dict):
                col = {"bullets": col}
            heading = str(col.get("heading") or "").strip()
            items = self._norm_bullets(col.get("bullets"), LIMITS["bullets"], f"{'左' if k == 0 else '右'}栏要点")
            x = CONTENT_LEFT_IN + k * (col_w + gap)
            self._rect(slide, x, top, col_w, h, th.paper_alt, MSO_SHAPE.ROUNDED_RECTANGLE)
            accent = th.primary if k == 0 else th.accent
            self._rect(slide, x, top, 0.08, h, accent)
            if heading:
                hs, hsize = self._fit_block([heading], col_w - 0.5, 0.6, 22, 18, 1.1, 0, "栏标题")
                self._text(slide, x + 0.3, top + 0.2, col_w - 0.5, 0.6, hs[0], size=hsize, bold=True, color=th.primary, anchor=MSO_ANCHOR.MIDDLE)
                self._rect(slide, x + 0.4, top + 0.85, 0.9, 0.04, accent)
            body_top = top + 1.05
            if items:
                items, size = self._fit_block(items, col_w - 0.5, h - 1.25, 22, BODY_MIN_PT, 1.2, 8.0, "要点")
                self._bullets(slide, x + 0.3, body_top, col_w - 0.5, h - 1.25, items, size, gap_pt=8.0)
            texts.extend(([heading] if heading else []) + items)
        self._body_chars(texts)
        pv["texts"] = texts

    # ---- chart ----
    def _slide_chart(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        self._frame(slide, title)
        self._chart_body(slide, s, pv, CONTENT_TOP_IN)

    def _chart_body(self, slide, s: Dict[str, Any], pv: Dict[str, Any], top: float) -> None:
        th = self.th
        spec = self._norm_chart(s.get("chart"))
        insight = str(s.get("insight") or "").strip()
        self._body_chars([insight])
        h = CONTENT_BOTTOM_IN - top
        ins_h = 0.75 if insight else 0.0
        if spec:
            self._add_chart(slide, spec, CONTENT_LEFT_IN, top, CONTENT_W_IN, h - ins_h - 0.1)
        else:
            self._picture(slide, "", CONTENT_LEFT_IN, top, CONTENT_W_IN, h - ins_h - 0.1)
        if insight:
            y = CONTENT_BOTTOM_IN - ins_h
            self._rect(slide, CONTENT_LEFT_IN, y, CONTENT_W_IN, ins_h, th.paper_alt, MSO_SHAPE.ROUNDED_RECTANGLE)
            self._rect(slide, CONTENT_LEFT_IN, y, 0.08, ins_h, th.accent)
            ins, isize = self._fit_block([insight], CONTENT_W_IN - 0.5, ins_h, 18, 16, 1.1, 0, "结论")
            self._text(slide, CONTENT_LEFT_IN + 0.3, y, CONTENT_W_IN - 0.5, ins_h, ins[0], size=isize, bold=True,
                       color=th.primary, anchor=MSO_ANCHOR.MIDDLE)
        pv["chart"] = (dict(spec, insight=insight) if spec else None)
        pv["texts"] = [insight] if insight else []

    # ---- table ----
    def _slide_table(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        self._frame(slide, title)
        self._table_body(slide, s, pv, CONTENT_TOP_IN)

    def _table_body(self, slide, s: Dict[str, Any], pv: Dict[str, Any], top: float) -> None:
        norm = self._norm_table(s)
        if not norm:
            self._picture(slide, "", CONTENT_LEFT_IN, top, CONTENT_W_IN, CONTENT_BOTTOM_IN - top)
            return
        cols, rows = norm
        self._body_chars(cols + [c for r in rows for c in r], check=False)   # 表格受行列数约束，不受 120 字约束
        self._add_table(slide, cols, rows, CONTENT_LEFT_IN, top + 0.05, CONTENT_W_IN, CONTENT_BOTTOM_IN - top - 0.1)
        pv["table"] = {"columns": cols, "rows": rows}

    # ---- kpi ----
    def _slide_kpi(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._frame(slide, title)
        kpis_in = s.get("kpis") or []
        kpis = [{"value": str(k.get("value") or ""), "label": str(k.get("label") or "")} for k in kpis_in if isinstance(k, dict)]
        if len(kpis) < LIMITS["kpiMin"]:
            self.warn(self.cur.index, f"kpi 仅 {len(kpis)} 张（建议 {LIMITS['kpiMin']}–{LIMITS['kpiMax']} 张）")
        if len(kpis) > LIMITS["kpiMax"]:
            self.warn(self.cur.index, f"kpi {len(kpis)} 张 > {LIMITS['kpiMax']}，已保留前 {LIMITS['kpiMax']} 张")
            kpis = kpis[:LIMITS["kpiMax"]]
        self._body_chars([k["value"] + k["label"] for k in kpis])
        n = max(len(kpis), 1)
        gap = 0.35
        card_w = min(3.6, (CONTENT_W_IN - gap * (n - 1)) / n)
        card_h = 2.9
        x0 = CONTENT_LEFT_IN + (CONTENT_W_IN - (card_w * n + gap * (n - 1))) / 2
        y = CONTENT_TOP_IN + (CONTENT_BOTTOM_IN - CONTENT_TOP_IN - card_h) / 2
        for i, k in enumerate(kpis):
            accent = th.accent if i % 2 else th.primary
            self._kpi_card(slide, x0 + i * (card_w + gap), y, card_w, card_h, k["value"], k["label"], accent)
        pv["kpis"] = kpis
        pv["texts"] = [k["label"] for k in kpis]

    # ---- timeline ----
    def _slide_timeline(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._frame(slide, title)
        steps_in = s.get("steps") or []
        steps = [{"when": str(st.get("when") or ""), "what": str(st.get("what") or "")} for st in steps_in if isinstance(st, dict)]
        if len(steps) < LIMITS["stepsMin"]:
            self.warn(self.cur.index, f"timeline 仅 {len(steps)} 步（建议 {LIMITS['stepsMin']}–{LIMITS['stepsMax']} 步）")
        if len(steps) > LIMITS["stepsMax"]:
            self.warn(self.cur.index, f"timeline {len(steps)} 步 > {LIMITS['stepsMax']}，已保留前 {LIMITS['stepsMax']} 步")
            steps = steps[:LIMITS["stepsMax"]]
        self._body_chars([st["when"] + st["what"] for st in steps])
        n = len(steps)
        if n == 0:
            return
        line_y = 3.75
        x0, x1 = CONTENT_LEFT_IN + 0.9, SLIDE_W_IN - CONTENT_LEFT_IN - 0.9
        spacing = (x1 - x0) / (n - 1) if n > 1 else 0
        self._rect(slide, x0 - 0.3, line_y - 0.025, (x1 - x0) + 0.6, 0.05, th.secondary_light)
        cell_w = min(spacing if n > 1 else 4.0, 3.4)
        d = 0.62
        for i, st in enumerate(steps):
            cx = x0 + i * spacing if n > 1 else (x0 + x1) / 2
            self._badge(slide, cx - d / 2, line_y - d / 2, d, str(i + 1), th.accent if i % 2 else th.primary, 18)
            whens, wsize = self._fit_block([st["when"]], cell_w, 0.9, 20, 16, 1.1, 0, "时间")
            self._text(slide, cx - cell_w / 2, line_y - d / 2 - 1.0, cell_w, 0.9, whens[0], size=wsize, bold=True,
                       color=th.primary, align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.BOTTOM)
            whats, hsize = self._fit_block([st["what"]], cell_w, 2.2, 18, 16, 1.15, 0, "步骤")
            self._text(slide, cx - cell_w / 2, line_y + d / 2 + 0.12, cell_w, 2.2, whats[0], size=hsize, color=th.ink,
                       align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.TOP)
        pv["steps"] = steps
        pv["texts"] = [f"{st['when']} {st['what']}".strip() for st in steps]

    # ---- team ----
    def _slide_team(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._frame(slide, title)
        members_in = s.get("members") or []
        members = [{"role": str(m.get("role") or ""), "duty": str(m.get("duty") or "")} for m in members_in if isinstance(m, dict)]
        if len(members) > LIMITS["members"]:
            self.warn(self.cur.index, f"team {len(members)} 人 > {LIMITS['members']}，已保留前 {LIMITS['members']} 人")
            members = members[:LIMITS["members"]]
        self._body_chars([m["role"] + m["duty"] for m in members])
        n = len(members)
        if n == 0:
            return
        per_row = n if n <= 3 else math.ceil(n / 2)
        rows = math.ceil(n / per_row)
        gap = 0.3
        card_w = min(3.7, (CONTENT_W_IN - gap * (per_row - 1)) / per_row)
        avail_h = CONTENT_BOTTOM_IN - CONTENT_TOP_IN
        card_h = (avail_h - gap * (rows - 1)) / rows
        x0 = CONTENT_LEFT_IN + (CONTENT_W_IN - (card_w * per_row + gap * (per_row - 1))) / 2
        for i, m in enumerate(members):
            r, c = divmod(i, per_row)
            x = x0 + c * (card_w + gap)
            y = CONTENT_TOP_IN + r * (card_h + gap)
            self._rect(slide, x, y, card_w, card_h, th.paper_alt, MSO_SHAPE.ROUNDED_RECTANGLE)
            compact = rows > 1
            d = 0.7 if compact else 0.9
            self._badge(slide, x + 0.3, y + 0.25, d, (m["role"] or "？")[0], th.primary if i % 2 == 0 else th.accent, 20 if compact else 24)
            roles, rsize = self._fit_block([m["role"]], card_w - d - 0.6, d, 20, 16, 1.1, 0, "角色")
            self._text(slide, x + d + 0.45, y + 0.25, card_w - d - 0.6, d, roles[0], size=rsize, bold=True, color=th.primary, anchor=MSO_ANCHOR.MIDDLE)
            duty_top = y + d + 0.45
            duties, dsize = self._fit_block([m["duty"]], card_w - 0.5, card_h - (d + 0.65), 18, 16, 1.2, 0, "分工")
            self._text(slide, x + 0.25, duty_top, card_w - 0.5, card_h - (d + 0.65), duties[0], size=dsize, color=th.ink)
        pv["members"] = members
        pv["texts"] = [f"{m['role']}：{m['duty']}" for m in members]

    # ---- quote ----
    def _slide_quote(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        th = self.th
        self._frame(slide, title)
        quote = str(s.get("quote") or "").strip()
        by = str(s.get("by") or "").strip()
        if not quote:
            self.warn(self.cur.index, "quote 页缺少 quote")
        self._body_chars([quote, by])
        box_x, box_w = 1.6, SLIDE_W_IN - 3.2
        self._text(slide, box_x - 0.5, 1.75, 1.2, 1.3, "“", size=96, bold=True, color=th.accent, line_spacing=1.0, track=False)
        qs, qsize = self._fit_block([quote], box_w - 0.6, 2.4, 36, 24, 1.2, 0, "引言")
        self._text(slide, box_x + 0.3, 2.45, box_w - 0.6, 2.4, qs[0], size=qsize, bold=True, color=th.primary,
                   align=PP_ALIGN.CENTER, anchor=MSO_ANCHOR.MIDDLE)
        self._text(slide, box_x + box_w - 0.7, 4.0, 1.2, 1.3, "”", size=96, bold=True, color=th.accent, line_spacing=1.0,
                   align=PP_ALIGN.RIGHT, track=False)
        self._rect(slide, SLIDE_W_IN / 2 - 0.8, 5.05, 1.6, 0.05, th.accent)
        if by:
            self._text(slide, box_x, 5.25, box_w, 0.6, f"—— {by}", size=18, color=th.muted, align=PP_ALIGN.CENTER)
        pv["quote"] = {"text": quote, "by": by}
        pv["texts"] = [quote] + ([by] if by else [])

    # ---- appendix ----
    def _slide_appendix(self, slide, s: Dict[str, Any], title: str, pv: Dict[str, Any]) -> None:
        self._frame(slide, title, appendix=True)
        top = CONTENT_TOP_IN
        if s.get("image"):
            self._image_body(slide, s, pv, top)
        elif s.get("columns") or s.get("rows"):
            self._table_body(slide, s, pv, top)
        elif s.get("chart"):
            self._chart_body(slide, s, pv, top)
        elif s.get("bullets") or s.get("highlight"):
            self._bullets_body(slide, s, pv, top + 0.1)
        else:
            self.warn(self.cur.index, "appendix 页没有 image / columns+rows / chart / bullets 内容")

    # ================================================================================ 总控
    def build(self) -> None:
        slides = self.cfg.get("slides")
        if not isinstance(slides, list) or not slides:
            raise ValueError("slides.json 需要非空的 slides 数组")
        self.total = len(slides)
        for raw in slides:
            if not isinstance(raw, dict):
                raw = {"layout": "bullets", "title": str(raw)}
            self._render_slide(raw)

    def save(self, out_path: Path) -> Path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        self.prs.save(str(out_path))
        preview_path = Path(str(out_path) + ".preview.json")
        write_json(preview_path, self.stats.preview)
        return preview_path


def render_deck(json_path: Path) -> Dict[str, Any]:
    cfg = read_json(json_path)
    if not isinstance(cfg, dict):
        raise ValueError("slides.json 顶层必须是对象")
    r = DeckRenderer(cfg, json_path)
    r.build()
    out_rel = str(cfg.get("output") or "产出/路演.pptx")
    out_path = Path(out_rel) if Path(out_rel).is_absolute() else (r.root / out_rel).resolve()
    preview_path = r.save(out_path)
    st = r.stats
    durations = [s.get("duration") for s in cfg["slides"] if isinstance(s, dict) and isinstance(s.get("duration"), (int, float))]
    log(f"✓ {len(st.pages)} 页 → {out_path}" + (f"，{len(st.warnings)} 条警告" if st.warnings else ""))
    return {
        "ok": True, "output": str(out_path), "preview": str(preview_path), "pages": len(st.pages),
        "warnings": list(st.warnings), "anonymizeHits": st.anonymize_hits, "missingImages": st.missing_images,
        "durationSeconds": sum(durations) if durations else None,
        "perPage": [{"index": p.index, "layout": p.layout, "chars": p.chars, "minFont": p.min_font} for p in st.pages],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    force_utf8_stdout()
    ap = argparse.ArgumentParser(description="slides.json → 16:9 路演 pptx")
    ap.add_argument("slides_json", help="slides.json 路径")
    args = ap.parse_args(argv)
    path = Path(args.slides_json).resolve()
    try:
        result = render_deck(path)
    except FileNotFoundError as e:
        print(json.dumps({"ok": False, "error": f"文件不存在：{e}"}, ensure_ascii=False))
        return 1
    except Exception as e:
        print(json.dumps({"ok": False, "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc(limit=5)}, ensure_ascii=False))
        return 1
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
