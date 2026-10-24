# -*- coding: utf-8 -*-
"""python-pptx 1.0.2 路演 PPT 生成工具（16:9、主题色/字体、页面类型、原生图表、组合图、模板复用）。"""
from __future__ import annotations

import copy
import re
from typing import Sequence

from lxml import etree
from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LABEL_POSITION, XL_LEGEND_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.opc.constants import RELATIONSHIP_TYPE as RT
from pptx.oxml import parse_xml
from pptx.chart.axis import ValueAxis
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

# ------------------------------------------------------------ 主题 ----

THEME = {
    "primary": "1F3864", "secondary": "2E75B6", "accent": "C9A227", "dark": "262626", "gray": "7F7F7F",
    "light": "F2F2F2", "white": "FFFFFF", "positive": "2E8B57", "negative": "C0392B",
    "series": ["1F3864", "2E75B6", "C9A227", "8FAADC", "7F7F7F", "BDD7EE"],
}
LATIN_FONT = "Arial"
EA_FONT = "微软雅黑"

SLIDE_W, SLIDE_H = Inches(13.333), Inches(7.5)


def new_presentation(latin=LATIN_FONT, ea=EA_FONT, theme=THEME, template: str | None = None):
    """16:9 空白演示文稿，并改写主题色与主题字体（主题字体 = 母版/版式默认字体）。"""
    prs = Presentation(template) if template else Presentation()
    prs.slide_width, prs.slide_height = SLIDE_W, SLIDE_H
    set_theme_fonts_and_colors(prs, latin, ea, theme)
    return prs


def set_theme_fonts_and_colors(prs, latin, ea, theme):
    """python-pptx 没有主题 API：直接改 /ppt/theme/theme1.xml 的 clrScheme 与 fontScheme。"""
    theme_part = prs.slide_master.part.part_related_by(RT.THEME)
    root = etree.fromstring(theme_part.blob)
    ns = {"a": "http://schemas.openxmlformats.org/drawingml/2006/main"}
    color_map = {"dk2": theme["primary"], "lt2": theme["light"], "accent1": theme["series"][0],
                 "accent2": theme["series"][1], "accent3": theme["series"][2], "accent4": theme["series"][3],
                 "accent5": theme["series"][4], "accent6": theme["series"][5], "hlink": theme["secondary"]}
    for name, hexv in color_map.items():
        el = root.find(f".//a:clrScheme/a:{name}", ns)
        for child in list(el):
            el.remove(child)
        srgb = etree.SubElement(el, qn("a:srgbClr"))
        srgb.set("val", hexv)
    for tag in ("majorFont", "minorFont"):
        font_el = root.find(f".//a:fontScheme/a:{tag}", ns)
        font_el.find("a:latin", ns).set("typeface", latin)
        font_el.find("a:ea", ns).set("typeface", ea)
        # 针对简体中文脚本(Hans)的字体映射，PowerPoint 优先用它
        for f in font_el.findall("a:font", ns):
            if f.get("script") == "Hans":
                font_el.remove(f)
        hans = etree.SubElement(font_el, qn("a:font"))
        hans.set("script", "Hans")
        hans.set("typeface", ea)
    theme_part._blob = etree.tostring(root, xml_declaration=True, encoding="UTF-8", standalone=True)


# ------------------------------------------------------------ 文本 ----

def _set_ea_font(rpr, typeface):
    """a:rPr 子元素有固定顺序：... a:latin, a:ea, a:cs, a:sym, a:hlinkClick ..."""
    for old in rpr.findall(qn("a:ea")):
        rpr.remove(old)
    ea = etree.Element(qn("a:ea"))
    ea.set("typeface", typeface)
    latin = rpr.find(qn("a:latin"))
    if latin is not None:
        latin.addnext(ea)
        return ea
    for tag in ("a:cs", "a:sym", "a:hlinkClick", "a:hlinkMouseOver", "a:rtl", "a:extLst"):
        nxt = rpr.find(qn(tag))
        if nxt is not None:
            nxt.addprevious(ea)
            return ea
    rpr.append(ea)
    return ea


def style_run(run, size=18, bold=False, color="262626", latin=LATIN_FONT, ea=EA_FONT, italic=None):
    f = run.font
    f.size = Pt(size)
    f.bold = bold
    if italic is not None:
        f.italic = italic
    f.color.rgb = RGBColor.from_string(color)
    f.name = latin                       # 写入 a:latin
    _set_ea_font(run._r.get_or_add_rPr(), ea)  # 中文必须另写 a:ea，否则 PowerPoint 用主题 ea 字体
    return run


def add_text(slide, left, top, width, height, text, size=18, bold=False, color="262626", align=PP_ALIGN.LEFT,
             anchor=MSO_ANCHOR.TOP, latin=LATIN_FONT, ea=EA_FONT, line_spacing=1.1, word_wrap=True):
    """text 可为 str 或 [(文本, {样式覆盖}), ...] 的段落列表；返回文本框。"""
    tb = slide.shapes.add_textbox(left, top, width, height)
    tf = tb.text_frame
    tf.word_wrap = word_wrap
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.05)
    paras = text if isinstance(text, list) else [(text, {})]
    for i, (t, override) in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = override.get("align", align)
        p.line_spacing = line_spacing
        if override.get("space_before") is not None:
            p.space_before = Pt(override["space_before"])
        r = p.add_run()
        r.text = t
        style_run(r, override.get("size", size), override.get("bold", bold), override.get("color", color), latin, ea)
    return tb


def add_bullets(slide, left, top, width, height, items: Sequence[str | tuple[str, int]], size=18, color="262626",
                bullet_color="C9A227", line_spacing=1.2, space_after=8):
    """items: 文本 或 (文本, 层级0/1)。用 a:buChar 写真正的项目符号。"""
    tb = slide.shapes.add_textbox(left, top, width, height)
    tf = tb.text_frame
    tf.word_wrap = True
    for i, it in enumerate(items):
        text, lvl = (it, 0) if isinstance(it, str) else it
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.level = lvl
        p.line_spacing = line_spacing
        p.space_after = Pt(space_after)
        pPr = p._p.get_or_add_pPr()
        indent = Inches(0.3 + 0.3 * lvl)
        pPr.set("marL", str(int(indent)))
        pPr.set("indent", str(-int(Inches(0.25))))
        for old in pPr.findall(qn("a:buNone")) + pPr.findall(qn("a:buChar")) + pPr.findall(qn("a:buClr")):
            pPr.remove(old)
        bu_clr = etree.SubElement(pPr, qn("a:buClr"))
        etree.SubElement(bu_clr, qn("a:srgbClr")).set("val", bullet_color)
        bu_font = etree.SubElement(pPr, qn("a:buFont"))
        bu_font.set("typeface", "Arial")
        bu = etree.SubElement(pPr, qn("a:buChar"))
        bu.set("char", "■" if lvl == 0 else "–")
        r = p.add_run()
        r.text = text
        style_run(r, size - 2 * lvl, False, color)
    return tb


# ---------------------------------------------------------- 形状 ----

def add_rect(slide, left, top, width, height, fill="1F3864", line=None, shape=MSO_SHAPE.RECTANGLE, shadow=False):
    shp = slide.shapes.add_shape(shape, left, top, width, height)
    shp.fill.solid()
    shp.fill.fore_color.rgb = RGBColor.from_string(fill)
    if line:
        shp.line.color.rgb = RGBColor.from_string(line)
        shp.line.width = Pt(1)
    else:
        shp.line.fill.background()
    if not shadow:
        shp.shadow.inherit = False
    return shp


def add_card(slide, left, top, width, height, title, value, sub="", fill="F2F2F2", accent="1F3864"):
    """KPI 卡片：左侧色条 + 大数字 + 说明。"""
    add_rect(slide, left, top, width, height, fill=fill, shape=MSO_SHAPE.ROUNDED_RECTANGLE)
    add_rect(slide, left, top, Inches(0.08), height, fill=accent)
    add_text(slide, left + Inches(0.2), top + Inches(0.1), width - Inches(0.3), Inches(0.4), title, size=12, color="595959")
    add_text(slide, left + Inches(0.2), top + Inches(0.45), width - Inches(0.3), Inches(0.7), value, size=28, bold=True,
             color=accent)
    if sub:
        add_text(slide, left + Inches(0.2), top + height - Inches(0.45), width - Inches(0.3), Inches(0.35), sub, size=11,
                 color="7F7F7F")


def add_number_badge(slide, left, top, number, size=Inches(0.5), fill="C9A227"):
    """圆形数字徽标（无图标库时用形状 + 数字/符号代替图标）。"""
    c = add_rect(slide, left, top, size, size, fill=fill, shape=MSO_SHAPE.OVAL)
    tf = c.text_frame
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.CENTER
    r = p.add_run()
    r.text = str(number)
    style_run(r, 14, True, "FFFFFF")
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    return c


def add_process_flow(slide, left, top, width, steps: Sequence[str], height=Inches(0.9), gap=Inches(0.08)):
    """V 形箭头流程图。"""
    n = len(steps)
    w = (width - gap * (n - 1)) / n
    for i, s in enumerate(steps):
        shp = add_rect(slide, left + (w + gap) * i, top, w, height,
                       fill=THEME["series"][i % len(THEME["series"])], shape=MSO_SHAPE.CHEVRON if i else MSO_SHAPE.PENTAGON)
        tf = shp.text_frame
        tf.margin_left = Inches(0.25)
        p = tf.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        r = p.add_run()
        r.text = s
        style_run(r, 13, True, "FFFFFF")
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE


# ---------------------------------------------------------- 页框架 ----

def _frame(slide, title, subtitle=None, page_no=None, footer_text="智守康联 | 商业计划路演"):
    add_rect(slide, 0, 0, SLIDE_W, Inches(0.12), fill=THEME["primary"])
    add_text(slide, Inches(0.6), Inches(0.3), Inches(11), Inches(0.7), title, size=28, bold=True, color=THEME["primary"])
    add_rect(slide, Inches(0.6), Inches(1.02), Inches(0.9), Inches(0.05), fill=THEME["accent"])
    if subtitle:
        add_text(slide, Inches(0.6), Inches(1.1), Inches(11.5), Inches(0.5), subtitle, size=14, color=THEME["gray"])
    add_text(slide, Inches(0.6), Inches(7.05), Inches(6), Inches(0.3), footer_text, size=9, color=THEME["gray"])
    if page_no is not None:
        add_text(slide, Inches(12.0), Inches(7.05), Inches(0.8), Inches(0.3), str(page_no), size=9, color=THEME["gray"],
                 align=PP_ALIGN.RIGHT)
    return slide


def blank_slide(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])  # 6 = Blank


# ---------------------------------------------------------- 页面类型 ----

def cover_slide(prs, title, subtitle, team, date, logo_path=None):
    s = blank_slide(prs)
    add_rect(s, 0, 0, SLIDE_W, SLIDE_H, fill=THEME["primary"])
    add_rect(s, 0, Inches(5.2), SLIDE_W, Inches(0.06), fill=THEME["accent"])
    add_text(s, Inches(0.9), Inches(2.2), Inches(11), Inches(1.3), title, size=44, bold=True, color="FFFFFF")
    add_text(s, Inches(0.9), Inches(3.5), Inches(11), Inches(0.8), subtitle, size=22, color="D6DCE5")
    add_text(s, Inches(0.9), Inches(5.5), Inches(11), Inches(0.5), f"{team}    |    {date}", size=14, color="D6DCE5")
    if logo_path:
        s.shapes.add_picture(logo_path, Inches(11.3), Inches(0.5), height=Inches(0.9))
    return s


def agenda_slide(prs, items: Sequence[str], page_no=None):
    s = _frame(blank_slide(prs), "目录", page_no=page_no)
    for i, it in enumerate(items):
        col, row = divmod(i, 4)
        left = Inches(0.9 + col * 6.2)
        top = Inches(1.7 + row * 1.25)
        add_number_badge(s, left, top, i + 1)
        add_text(s, left + Inches(0.7), top + Inches(0.02), Inches(5), Inches(0.5), it, size=20, color=THEME["dark"])
    return s


def section_slide(prs, number, title, tagline=""):
    s = blank_slide(prs)
    add_rect(s, 0, 0, Inches(4.2), SLIDE_H, fill=THEME["primary"])
    add_text(s, Inches(0.6), Inches(2.6), Inches(3), Inches(1.5), f"{number:02d}", size=66, bold=True, color=THEME["accent"])
    add_text(s, Inches(4.8), Inches(2.8), Inches(8), Inches(1), title, size=36, bold=True, color=THEME["primary"])
    if tagline:
        add_text(s, Inches(4.8), Inches(3.8), Inches(8), Inches(0.8), tagline, size=16, color=THEME["gray"])
    return s


def bullets_slide(prs, title, items, subtitle=None, page_no=None, notes=None):
    s = _frame(blank_slide(prs), title, subtitle, page_no)
    add_bullets(s, Inches(0.8), Inches(1.7), Inches(11.5), Inches(5), items, size=20)
    if notes:
        set_notes(s, notes)
    return s


def image_text_slide(prs, title, image_path, bullets, subtitle=None, page_no=None, image_left=True, crop=None):
    """crop=(left, top, right, bottom) 以比例裁剪图片（0.1 = 裁掉 10%）。"""
    s = _frame(blank_slide(prs), title, subtitle, page_no)
    img_left = Inches(0.7) if image_left else Inches(7.0)
    txt_left = Inches(7.0) if image_left else Inches(0.7)
    pic = s.shapes.add_picture(image_path, img_left, Inches(1.7), width=Inches(5.8))
    if crop:
        pic.crop_left, pic.crop_top, pic.crop_right, pic.crop_bottom = crop
    if pic.height > Inches(5):  # 限高
        ratio = Inches(5) / pic.height
        pic.height = Inches(5)
        pic.width = int(pic.width * ratio)
    add_bullets(s, txt_left, Inches(1.8), Inches(5.6), Inches(5), bullets, size=18)
    return s


def chart_slide(prs, title, chart_builder, subtitle=None, page_no=None, takeaway=None, **kwargs):
    """chart_builder(slide, left, top, width, height, **kwargs) -> chart"""
    s = _frame(blank_slide(prs), title, subtitle, page_no)
    chart = chart_builder(s, Inches(0.7), Inches(1.6), Inches(8.3), Inches(5.2), **kwargs)
    if takeaway:
        add_rect(s, Inches(9.3), Inches(1.7), Inches(3.4), Inches(3.2), fill=THEME["light"], shape=MSO_SHAPE.ROUNDED_RECTANGLE)
        add_text(s, Inches(9.5), Inches(1.85), Inches(3.0), Inches(0.4), "核心结论", size=13, bold=True, color=THEME["accent"])
        add_text(s, Inches(9.5), Inches(2.3), Inches(3.0), Inches(2.5), takeaway, size=14, color=THEME["dark"])
    return s, chart


def table_slide(prs, title, header, rows, col_widths_in=None, subtitle=None, page_no=None, font_size=12,
                highlight_last_row=False):
    s = _frame(blank_slide(prs), title, subtitle, page_no)
    n_rows, n_cols = len(rows) + 1, len(header)
    width = Inches(sum(col_widths_in)) if col_widths_in else Inches(11.8)
    gt = s.shapes.add_table(n_rows, n_cols, Inches(0.7), Inches(1.7), width, Inches(0.4) * n_rows)
    table = gt.table
    # 关掉内置表格样式的条纹/首列，保证自定义底色生效
    tbl_pr = gt._element.graphic.graphicData.tbl.tblPr
    tbl_pr.set("firstRow", "1")
    tbl_pr.set("bandRow", "0")
    if col_widths_in:
        for j, w in enumerate(col_widths_in):
            table.columns[j].width = Inches(w)
    for j, h in enumerate(header):
        c = table.cell(0, j)
        c.text = ""
        c.fill.solid()
        c.fill.fore_color.rgb = RGBColor.from_string(THEME["primary"])
        r = c.text_frame.paragraphs[0].add_run()
        r.text = str(h)
        style_run(r, font_size, True, "FFFFFF")
        c.text_frame.paragraphs[0].alignment = PP_ALIGN.CENTER
        c.vertical_anchor = MSO_ANCHOR.MIDDLE
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            c = table.cell(i, j)
            c.text = ""
            c.fill.solid()
            is_last = highlight_last_row and i == len(rows)
            c.fill.fore_color.rgb = RGBColor.from_string("FFF4D6" if is_last else ("F2F2F2" if i % 2 == 0 else "FFFFFF"))
            r = c.text_frame.paragraphs[0].add_run()
            r.text = str(v)
            style_run(r, font_size, is_last, THEME["dark"])
            c.text_frame.paragraphs[0].alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.RIGHT
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
    return s, table


def team_slide(prs, title, members: Sequence[dict], page_no=None):
    """members: [{name, role, bio, photo(optional)}]"""
    s = _frame(blank_slide(prs), title, page_no=page_no)
    n = len(members)
    card_w = Inches(2.7)
    gap = (SLIDE_W - Inches(1.4) - card_w * n) / max(n - 1, 1)
    for i, m in enumerate(members):
        left = Inches(0.7) + (card_w + gap) * i
        add_rect(s, left, Inches(1.7), card_w, Inches(4.9), fill=THEME["light"], shape=MSO_SHAPE.ROUNDED_RECTANGLE)
        if m.get("photo"):
            pic = s.shapes.add_picture(m["photo"], left + Inches(0.75), Inches(1.9), width=Inches(1.2), height=Inches(1.2))
            pic.crop_left = pic.crop_right = 0.0
        else:
            add_number_badge(s, left + Inches(0.75), Inches(1.9), m["name"][0], size=Inches(1.2), fill=THEME["secondary"])
        add_text(s, left, Inches(3.25), card_w, Inches(0.4), m["name"], size=18, bold=True, color=THEME["primary"],
                 align=PP_ALIGN.CENTER)
        add_text(s, left, Inches(3.65), card_w, Inches(0.4), m["role"], size=12, color=THEME["accent"], align=PP_ALIGN.CENTER)
        add_text(s, left + Inches(0.15), Inches(4.1), card_w - Inches(0.3), Inches(2.3), m["bio"], size=11,
                 color="595959", align=PP_ALIGN.LEFT)
    return s


def financial_slide(prs, title, kpis: Sequence[tuple[str, str, str]], chart_builder, page_no=None, **kwargs):
    """上方 4 个 KPI 卡片 + 下方原生图表。"""
    s = _frame(blank_slide(prs), title, page_no=page_no)
    w = Inches(2.85)
    for i, (k, v, sub) in enumerate(kpis[:4]):
        add_card(s, Inches(0.7) + (w + Inches(0.17)) * i, Inches(1.5), w, Inches(1.3), k, v, sub,
                 accent=THEME["series"][i % 3])
    chart = chart_builder(s, Inches(0.7), Inches(3.0), Inches(11.9), Inches(3.9), **kwargs)
    return s, chart


def closing_slide(prs, title="感谢聆听", subtitle="期待与您共创智慧养老新生态", contact="contact@example.com"):
    s = blank_slide(prs)
    add_rect(s, 0, 0, SLIDE_W, SLIDE_H, fill=THEME["primary"])
    add_text(s, 0, Inches(2.6), SLIDE_W, Inches(1.2), title, size=48, bold=True, color="FFFFFF", align=PP_ALIGN.CENTER)
    add_text(s, 0, Inches(3.8), SLIDE_W, Inches(0.7), subtitle, size=20, color="D6DCE5", align=PP_ALIGN.CENTER)
    add_text(s, 0, Inches(5.0), SLIDE_W, Inches(0.5), contact, size=14, color=THEME["accent"], align=PP_ALIGN.CENTER)
    return s


def set_notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text


# ---------------------------------------------------------- 原生图表 ----

def set_chart_font(chart, size=12, latin=LATIN_FONT, ea=EA_FONT, color="404040"):
    """chart.font 只写 a:latin；中文轴标签/图例还要补 a:ea。"""
    chart.font.size = Pt(size)
    chart.font.name = latin
    chart.font.color.rgb = RGBColor.from_string(color)
    tx_pr = chart._chartSpace.find(qn("c:txPr"))
    for def_rpr in tx_pr.iter(qn("a:defRPr")):
        _set_ea_font(def_rpr, ea)


def _apply_series_colors(plot, colors, line=False):
    for i, s in enumerate(plot.series):
        c = RGBColor.from_string(colors[i % len(colors)])
        if line:
            s.format.line.color.rgb = c
            s.format.line.width = Pt(2.25)
            s.smooth = False
            s.marker.format.fill.solid()
            s.marker.format.fill.fore_color.rgb = c
        else:
            s.format.fill.solid()
            s.format.fill.fore_color.rgb = c


def _labels(plot, number_format='#,##0', position=None, size=10, color="404040"):
    plot.has_data_labels = True
    dl = plot.data_labels
    dl.number_format = number_format
    dl.number_format_is_linked = False
    dl.show_value = True
    if position is not None:
        dl.position = position
    dl.font.size = Pt(size)
    dl.font.color.rgb = RGBColor.from_string(color)


def add_bar_chart(slide, left, top, width, height, categories, series: dict[str, Sequence[float]], title=None,
                  stacked=False, number_format='#,##0', colors=None, legend=True, horizontal=False):
    cd = CategoryChartData()
    cd.categories = list(categories)
    for name, vals in series.items():
        cd.add_series(name, list(vals))
    if horizontal:
        ctype = XL_CHART_TYPE.BAR_STACKED if stacked else XL_CHART_TYPE.BAR_CLUSTERED
    else:
        ctype = XL_CHART_TYPE.COLUMN_STACKED if stacked else XL_CHART_TYPE.COLUMN_CLUSTERED
    chart = slide.shapes.add_chart(ctype, left, top, width, height, cd).chart
    set_chart_font(chart)
    plot = chart.plots[0]
    plot.gap_width = 80
    if stacked:
        plot.overlap = 100
    _apply_series_colors(plot, colors or THEME["series"])
    _labels(plot, number_format, None if stacked else XL_LABEL_POSITION.OUTSIDE_END,
            color="FFFFFF" if stacked else "404040")  # 堆叠柱标签落在深色柱内，用白字
    chart.has_legend = legend and len(series) > 1
    if chart.has_legend:
        chart.legend.position = XL_LEGEND_POSITION.BOTTOM
        chart.legend.include_in_layout = False
    chart.value_axis.has_major_gridlines = True
    chart.value_axis.major_gridlines.format.line.color.rgb = RGBColor.from_string("D9D9D9")
    chart.value_axis.format.line.fill.background()
    chart.value_axis.tick_labels.number_format = number_format
    chart.value_axis.tick_labels.number_format_is_linked = False
    if title:
        chart.has_title = True
        chart.chart_title.text_frame.text = title
        style_run(chart.chart_title.text_frame.paragraphs[0].runs[0], 14, True, THEME["dark"])
    return chart


def add_line_chart(slide, left, top, width, height, categories, series, title=None, number_format='#,##0', colors=None):
    cd = CategoryChartData()
    cd.categories = list(categories)
    for name, vals in series.items():
        cd.add_series(name, list(vals))
    chart = slide.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS, left, top, width, height, cd).chart
    set_chart_font(chart)
    plot = chart.plots[0]
    _apply_series_colors(plot, colors or THEME["series"], line=True)
    _labels(plot, number_format, XL_LABEL_POSITION.ABOVE)
    chart.has_legend = len(series) > 1
    if chart.has_legend:
        chart.legend.position = XL_LEGEND_POSITION.BOTTOM
        chart.legend.include_in_layout = False
    chart.value_axis.major_gridlines.format.line.color.rgb = RGBColor.from_string("D9D9D9")
    if title:
        chart.has_title = True
        chart.chart_title.text_frame.text = title
        style_run(chart.chart_title.text_frame.paragraphs[0].runs[0], 14, True, THEME["dark"])
    return chart


def add_pie_chart(slide, left, top, width, height, categories, values, title=None, donut=True, colors=None,
                  number_format='0.0%'):
    cd = CategoryChartData()
    cd.categories = list(categories)
    cd.add_series("占比", list(values))
    ctype = XL_CHART_TYPE.DOUGHNUT if donut else XL_CHART_TYPE.PIE
    chart = slide.shapes.add_chart(ctype, left, top, width, height, cd).chart
    set_chart_font(chart)
    plot = chart.plots[0]
    colors = colors or THEME["series"]
    for i, pt in enumerate(plot.series[0].points):  # 饼图逐点上色
        pt.format.fill.solid()
        pt.format.fill.fore_color.rgb = RGBColor.from_string(colors[i % len(colors)])
    plot.has_data_labels = True
    dl = plot.data_labels
    dl.show_percentage = True
    dl.show_value = False
    dl.show_category_name = True
    dl.number_format = number_format
    dl.number_format_is_linked = False
    dl.font.size = Pt(11)
    if not donut:
        dl.position = XL_LABEL_POSITION.OUTSIDE_END
    chart.has_legend = True
    chart.legend.position = XL_LEGEND_POSITION.RIGHT
    chart.legend.include_in_layout = False
    if title:
        chart.has_title = True
        chart.chart_title.text_frame.text = title
        style_run(chart.chart_title.text_frame.paragraphs[0].runs[0], 14, True, THEME["dark"])
    return chart


def add_combo_chart(slide, left, top, width, height, categories, bar_series: dict, line_series: dict,
                    bar_format='#,##0', line_format='0%', title=None, colors=None):
    """柱 + 折线(次坐标轴) 组合图。python-pptx 无原生 API：先建柱图，再把折线系列搬进新建的 c:lineChart。"""
    cd = CategoryChartData()
    cd.categories = list(categories)
    for name, vals in bar_series.items():
        cd.add_series(name, list(vals))
    for name, vals in line_series.items():
        cd.add_series(name, list(vals))
    gframe = slide.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, left, top, width, height, cd)
    chart = gframe.chart
    set_chart_font(chart)
    colors = colors or THEME["series"]
    plot_area = chart._chartSpace.chart.plotArea
    bar_chart = plot_area.find(qn("c:barChart"))
    sers = bar_chart.findall(qn("c:ser"))
    n_bar = len(bar_series)
    # 新 lineChart + 次坐标轴
    line_chart = etree.SubElement(plot_area, qn("c:lineChart"))
    etree.SubElement(line_chart, qn("c:grouping")).set("val", "standard")
    etree.SubElement(line_chart, qn("c:varyColors")).set("val", "0")
    for ser in sers[n_bar:]:
        bar_chart.remove(ser)
        # 去掉柱系列专有元素 invertIfNegative
        for inv in ser.findall(qn("c:invertIfNegative")):
            ser.remove(inv)
        line_chart.append(ser)
    etree.SubElement(line_chart, qn("c:marker")).set("val", "1")
    ax_ids = [el.get("val") for el in bar_chart.findall(qn("c:axId"))]
    sec_cat_id, sec_val_id = "50010", "50020"
    etree.SubElement(line_chart, qn("c:axId")).set("val", sec_cat_id)
    etree.SubElement(line_chart, qn("c:axId")).set("val", sec_val_id)
    # plotArea 里 chart 元素必须在 axis 元素之前：把 lineChart 移到第一个 c:catAx 之前
    first_ax = plot_area.find(qn("c:catAx"))
    first_ax.addprevious(line_chart)
    # 次分类轴(隐藏) + 次数值轴(右侧)
    sec_cat = parse_xml(  # 用 pptx.oxml.parse_xml 才能得到 python-pptx 的自定义元素类
        f'<c:catAx xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart"><c:axId val="{sec_cat_id}"/>'
        f'<c:scaling><c:orientation val="minMax"/></c:scaling><c:delete val="1"/><c:axPos val="b"/>'
        f'<c:majorTickMark val="out"/><c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/>'
        f'<c:crossAx val="{sec_val_id}"/><c:crosses val="autoZero"/><c:auto val="1"/><c:lblAlgn val="ctr"/>'
        f'<c:lblOffset val="100"/><c:noMultiLvlLbl val="0"/></c:catAx>')
    sec_val = parse_xml(
        f'<c:valAx xmlns:c="http://schemas.openxmlformats.org/drawingml/2006/chart"><c:axId val="{sec_val_id}"/>'
        f'<c:scaling><c:orientation val="minMax"/></c:scaling><c:delete val="0"/><c:axPos val="r"/>'
        f'<c:numFmt formatCode="{line_format}" sourceLinked="0"/><c:majorTickMark val="out"/>'
        f'<c:minorTickMark val="none"/><c:tickLblPos val="nextTo"/><c:crossAx val="{sec_cat_id}"/>'
        f'<c:crosses val="max"/><c:crossBetween val="between"/></c:valAx>')
    last_ax = plot_area.findall(qn("c:valAx"))[-1]
    last_ax.addnext(sec_cat)
    sec_cat.addnext(sec_val)
    # 样式
    bar_plot = chart.plots[0]
    bar_plot.gap_width = 80
    _apply_series_colors(bar_plot, colors)
    _labels(bar_plot, bar_format, XL_LABEL_POSITION.OUTSIDE_END)
    line_plot = chart.plots[1]
    _apply_series_colors(line_plot, [THEME["accent"]] + colors[2:], line=True)
    _labels(line_plot, line_format, XL_LABEL_POSITION.ABOVE, color=THEME["accent"])
    chart.has_legend = True
    chart.legend.position = XL_LEGEND_POSITION.BOTTOM
    chart.legend.include_in_layout = False
    # 注意：有两个 valAx 时 chart.value_axis 返回的是第 2 个(次轴)，主轴要显式取
    primary_val = ValueAxis(plot_area.findall(qn("c:valAx"))[0])
    primary_val.tick_labels.number_format = bar_format
    primary_val.tick_labels.number_format_is_linked = False
    primary_val.has_major_gridlines = True
    primary_val.major_gridlines.format.line.color.rgb = RGBColor.from_string("D9D9D9")
    if title:
        chart.has_title = True
        chart.chart_title.text_frame.text = title
        style_run(chart.chart_title.text_frame.paragraphs[0].runs[0], 14, True, THEME["dark"])
    return chart


# ---------------------------------------------------------- 模板复用 ----

def fill_placeholders(slide, mapping: dict[int, str]):
    """按占位符 idx 填文本，保留占位符原格式（只改第一个 run 的文本，删除其余 run）。"""
    for ph in slide.placeholders:
        idx = ph.placeholder_format.idx
        if idx not in mapping or not ph.has_text_frame:
            continue
        tf = ph.text_frame
        lines = mapping[idx].split("\n")
        first_p = tf.paragraphs[0]
        for extra in tf.paragraphs[1:]:
            extra._p.getparent().remove(extra._p)
        runs = first_p.runs
        if runs:
            runs[0].text = lines[0]
            for r in runs[1:]:
                r._r.getparent().remove(r._r)
        else:
            first_p.add_run().text = lines[0]
        for line in lines[1:]:
            p = tf.add_paragraph()
            p.text = line


def replace_text_everywhere(prs, mapping: dict[str, str]):
    """模板里 {{占位}} 风格的文本替换（跨 run 的占位符需先合并 run，这里按段落级合并后重写）。"""
    pattern = re.compile("|".join(re.escape(k) for k in mapping))
    for slide in prs.slides:
        for shape in slide.shapes:
            frames = []
            if shape.has_text_frame:
                frames.append(shape.text_frame)
            if shape.has_table:
                for row in shape.table.rows:
                    for cell in row.cells:
                        frames.append(cell.text_frame)
            for tf in frames:
                for p in tf.paragraphs:
                    full = "".join(r.text for r in p.runs)
                    if pattern.search(full):
                        new = pattern.sub(lambda m: mapping[m.group(0)], full)
                        if p.runs:
                            p.runs[0].text = new
                            for r in p.runs[1:]:
                                r._r.getparent().remove(r._r)


def delete_slide(prs, index):
    sld_id_lst = prs.slides._sldIdLst
    sld = list(sld_id_lst)[index]
    prs.part.drop_rel(sld.rId)
    sld_id_lst.remove(sld)


def duplicate_slide(prs, index):
    """复制一页（文本/形状/图片可用；图表部件会被两页共享，不建议复制图表页）。"""
    src = prs.slides[index]
    dst = prs.slides.add_slide(src.slide_layout)
    for shp in list(dst.shapes):
        shp._element.getparent().remove(shp._element)
    rid_map = {}
    for rel in src.part.rels.values():
        if rel.reltype in (RT.SLIDE_LAYOUT, RT.NOTES_SLIDE):
            continue
        if rel.is_external:
            new_rid = dst.part.rels.get_or_add_ext_rel(rel.reltype, rel.target_ref)
        else:
            new_rid = dst.part.rels.get_or_add(rel.reltype, rel.target_part)
        rid_map[rel.rId] = new_rid
    for shp in src.shapes:
        el = copy.deepcopy(shp._element)
        for node in el.iter():
            for attr in ("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed",
                         "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id",
                         "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}link"):
                if attr in node.attrib and node.attrib[attr] in rid_map:
                    node.attrib[attr] = rid_map[node.attrib[attr]]
        dst.shapes._spTree.append(el)
    return dst
