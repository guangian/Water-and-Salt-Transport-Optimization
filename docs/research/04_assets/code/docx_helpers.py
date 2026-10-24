# -*- coding: utf-8 -*-
"""python-docx 中文商业文档工具函数（python-docx 1.2.0 验证通过）。

覆盖：中文字体(w:eastAsia)、标题多级编号、TOC 域、页眉页脚页码域、封面、图片+SEQ 题注、
三线表/底纹表、分页分节(横向页)、段落格式、脚注、超链接/书签、列表编号重启、批注。
"""
from __future__ import annotations

import copy
from typing import Iterable, Sequence

from docx import Document
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.style import WD_STYLE_TYPE
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.opc.constants import RELATIONSHIP_TYPE as RT
from docx.opc.packuri import PackURI
from docx.opc.part import XmlPart
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn
from docx.shared import Cm, Mm, Pt, RGBColor

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
THEME_FONT_ATTRS = ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme")


# ------------------------------------------------- 按 schema 顺序插入子元素 ----
# Word 对 pPr/tcPr/tblPr/trPr 子元素顺序敏感（顺序错误可能报“发现无法读取的内容”），
# python-docx 只对它自己封装过的元素保证顺序，手工追加的元素要用下面的函数插入。
_ORDER = {
    "pPr": ["pStyle", "keepNext", "keepLines", "pageBreakBefore", "framePr", "widowControl", "numPr",
            "suppressLineNumbers", "pBdr", "shd", "tabs", "suppressAutoHyphens", "kinsoku", "wordWrap",
            "overflowPunct", "topLinePunct", "autoSpaceDE", "autoSpaceDN", "bidi", "adjustRightInd", "snapToGrid",
            "spacing", "ind", "contextualSpacing", "mirrorIndents", "suppressOverlap", "jc", "textDirection",
            "textAlignment", "textboxTightWrap", "outlineLvl", "divId", "cnfStyle", "rPr", "sectPr", "pPrChange"],
    "tcPr": ["cnfStyle", "tcW", "gridSpan", "hMerge", "vMerge", "tcBorders", "shd", "noWrap", "tcMar",
             "textDirection", "tcFitText", "vAlign", "hideMark", "headers", "cellIns", "cellDel", "cellMerge",
             "tcPrChange"],
    "tblPr": ["tblStyle", "tblpPr", "tblOverlap", "bidiVisual", "tblStyleRowBandSize", "tblStyleColBandSize",
              "tblW", "jc", "tblCellSpacing", "tblInd", "tblBorders", "shd", "tblLayout", "tblCellMar", "tblLook",
              "tblCaption", "tblDescription", "tblPrChange"],
    "trPr": ["cnfStyle", "divId", "gridBefore", "gridAfter", "wBefore", "wAfter", "cantSplit", "trHeight",
             "tblHeader", "tblCellSpacing", "jc", "hidden", "ins", "del", "trPrChange"],
}


def insert_ordered(parent, el):
    """把 el 插入 parent（pPr/tcPr/tblPr/trPr）中符合 schema 顺序的位置；同名旧元素先删除。"""
    order = _ORDER[parent.tag.split("}")[1]]
    name = el.tag.split("}")[1]
    for old in parent.findall(el.tag):
        parent.remove(old)
    successors = {qn("w:" + n) for n in order[order.index(name) + 1:]}
    for child in parent:
        if child.tag in successors:
            child.addprevious(el)
            return el
    parent.append(el)
    return el


# ---------------------------------------------------------------- 字体 ----

def _rfonts_of(rpr):
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts")
        rpr.insert(0, rfonts)  # rFonts 必须是 rPr 的第一个子元素
    # 主题字体属性优先级高于显式字体名，必须删掉，否则 Word 里仍显示主题字体
    for a in THEME_FONT_ATTRS:
        rfonts.attrib.pop(qn(a), None)
    return rfonts


def set_fonts(rpr, east_asia: str = "宋体", latin: str = "Times New Roman"):
    """在任意 rPr(run / style / docDefaults) 上设置 中文字体 + 西文字体。"""
    rfonts = _rfonts_of(rpr)
    rfonts.set(qn("w:ascii"), latin)
    rfonts.set(qn("w:hAnsi"), latin)
    rfonts.set(qn("w:eastAsia"), east_asia)
    rfonts.set(qn("w:cs"), latin)
    return rfonts


def set_run_font(run, east_asia="宋体", latin="Times New Roman", size: float | None = None,
                 bold: bool | None = None, italic: bool | None = None,
                 color: str | None = None):
    set_fonts(run._element.get_or_add_rPr(), east_asia, latin)
    if size is not None:
        run.font.size = Pt(size)
    if bold is not None:
        run.font.bold = bold
    if italic is not None:
        run.font.italic = italic
    if color:
        run.font.color.rgb = RGBColor.from_string(color)
    return run


def set_style_font(style, east_asia="宋体", latin="Times New Roman", size: float | None = None,
                   bold: bool | None = None, color: str | None = None):
    set_fonts(style.element.get_or_add_rPr(), east_asia, latin)
    if size is not None:
        style.font.size = Pt(size)
    if bold is not None:
        style.font.bold = bold
    if color:
        style.font.color.rgb = RGBColor.from_string(color)
    return style


def setup_chinese_document(doc, body_east_asia="宋体", heading_east_asia="黑体",
                           latin="Times New Roman", body_size=12, a4=True):
    """把 python-docx 默认模板改造成中文商业文档底座。"""
    # 1) docDefaults：全文默认字体 + 语言 zh-CN
    styles_el = doc.styles.element
    rpr_default = styles_el.find(qn("w:docDefaults")).find(qn("w:rPrDefault")).find(qn("w:rPr"))
    set_fonts(rpr_default, body_east_asia, latin)
    lang = rpr_default.find(qn("w:lang"))
    if lang is None:
        lang = OxmlElement("w:lang")
        rpr_default.append(lang)
    lang.set(qn("w:val"), "en-US")
    lang.set(qn("w:eastAsia"), "zh-CN")

    # 2) 正文样式：宋体小四、1.5 倍行距、首行缩进 2 字符、段后 0
    normal = doc.styles["Normal"]
    set_style_font(normal, body_east_asia, latin, size=body_size)
    pf = normal.paragraph_format
    pf.line_spacing = 1.5
    pf.space_after = Pt(0)
    pf.space_before = Pt(0)
    ind = normal.element.get_or_add_pPr().get_or_add_ind()
    ind.set(qn("w:firstLineChars"), "200")  # 以“字符”为单位的首行缩进，随字号自适应

    # 3) 标题样式：黑体、黑色、去掉主题色；标题不缩进
    sizes = {1: 22, 2: 16, 3: 14, 4: 12}
    for lvl in range(1, 5):
        st = doc.styles[f"Heading {lvl}"]
        set_style_font(st, heading_east_asia, latin, size=sizes[lvl], bold=True, color="000000")
        st.paragraph_format.first_line_indent = Pt(0)
        st.paragraph_format.space_before = Pt(12 if lvl > 1 else 18)
        st.paragraph_format.space_after = Pt(6)
        st.paragraph_format.line_spacing = 1.25
        st.paragraph_format.keep_with_next = True
        st.element.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "0")
    doc.styles["Heading 1"].paragraph_format.page_break_before = False

    # 4) 题注、引用、列表等样式字体
    set_style_font(doc.styles["Caption"], heading_east_asia, latin, size=10.5, bold=False, color="000000")
    doc.styles["Caption"].paragraph_format.alignment = WD_ALIGN_PARAGRAPH.CENTER
    doc.styles["Caption"].paragraph_format.keep_with_next = True  # 表题注与表同页（图题注在图后，见 add_figure）
    doc.styles["Caption"].element.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "0")
    for name in ("List Bullet", "List Bullet 2", "List Number", "List Number 2", "Quote", "Intense Quote",
                 "Header", "Footer", "Title", "Subtitle", "TOC Heading"):
        st = doc.styles[name]
        st.element.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "0")
    set_style_font(doc.styles["Title"], heading_east_asia, latin, size=26, bold=True, color="1F3864")
    set_style_font(doc.styles["Header"], body_east_asia, latin, size=9)
    set_style_font(doc.styles["Footer"], body_east_asia, latin, size=9)

    # 5) A4 + 常用页边距（默认模板是 Letter！）
    if a4:
        for sec in doc.sections:
            sec.page_width, sec.page_height = Mm(210), Mm(297)
            sec.left_margin = sec.right_margin = Cm(3.0)
            sec.top_margin = sec.bottom_margin = Cm(2.54)
            sec.header_distance = sec.footer_distance = Cm(1.5)
    return doc


# ------------------------------------------------------- 多级标题编号 ----

def add_heading_numbering(doc, levels=3, fmt: Sequence[str] | None = None,
                          east_asia="黑体", latin="Times New Roman"):
    """给 Heading 1..levels 挂多级编号（1 / 1.1 / 1.1.1）。返回 numId。"""
    numbering = doc.part.numbering_part.element
    fmt = fmt or ["%1", "%1.%2", "%1.%2.%3", "%1.%2.%3.%4"]
    abs_ids = [int(a.get(qn("w:abstractNumId"))) for a in numbering.findall(qn("w:abstractNum"))]
    abs_id = max(abs_ids, default=-1) + 1
    lvls = []
    for i in range(levels):
        lvls.append(
            f'<w:lvl w:ilvl="{i}"><w:start w:val="1"/><w:numFmt w:val="decimal"/>'
            f'<w:pStyle w:val="Heading{i + 1}"/><w:suff w:val="space"/>'
            f'<w:lvlText w:val="{fmt[i]}"/><w:lvlJc w:val="left"/>'
            f'<w:pPr><w:ind w:left="0" w:firstLine="0"/></w:pPr>'
            f'<w:rPr><w:rFonts w:ascii="{latin}" w:hAnsi="{latin}" w:eastAsia="{east_asia}"/></w:rPr></w:lvl>'
        )
    abstract = parse_xml(
        f'<w:abstractNum {nsdecls("w")} w:abstractNumId="{abs_id}">'
        f'<w:multiLevelType w:val="multilevel"/>{"".join(lvls)}</w:abstractNum>'
    )
    existing = numbering.findall(qn("w:abstractNum"))
    if existing:
        existing[-1].addnext(abstract)  # abstractNum 必须排在所有 num 之前
    else:
        numbering.insert(0, abstract)
    num_ids = [int(n.get(qn("w:numId"))) for n in numbering.findall(qn("w:num"))]
    num_id = max(num_ids, default=0) + 1
    numbering.append(parse_xml(
        f'<w:num {nsdecls("w")} w:numId="{num_id}"><w:abstractNumId w:val="{abs_id}"/></w:num>'))
    for i in range(levels):
        ppr = doc.styles[f"Heading {i + 1}"].element.get_or_add_pPr()
        numpr = ppr.get_or_add_numPr()
        numpr.get_or_add_ilvl().val = i
        numpr.get_or_add_numId().val = num_id
    # "TOC Heading" 基于 Heading 1，会继承编号；numId=0 表示“不编号”
    toc_ppr = doc.styles["TOC Heading"].element.get_or_add_pPr()
    toc_numpr = toc_ppr.get_or_add_numPr()
    toc_numpr.get_or_add_ilvl().val = 0
    toc_numpr.get_or_add_numId().val = 0
    return num_id


def restart_list_numbering(doc, paragraph, style_name="List Number"):
    """让某段成为新编号列表的第一项(从 1 重新开始)。后续同列表段落调用 continue_list_numbering。"""
    numbering = doc.part.numbering_part.element
    style_numpr = doc.styles[style_name].element.pPr.numPr
    base_num_id = style_numpr.numId.val
    abs_id = None
    for n in numbering.findall(qn("w:num")):
        if int(n.get(qn("w:numId"))) == base_num_id:
            abs_id = n.find(qn("w:abstractNumId")).get(qn("w:val"))
    num_id = max(int(n.get(qn("w:numId"))) for n in numbering.findall(qn("w:num"))) + 1
    numbering.append(parse_xml(
        f'<w:num {nsdecls("w")} w:numId="{num_id}"><w:abstractNumId w:val="{abs_id}"/>'
        f'<w:lvlOverride w:ilvl="0"><w:startOverride w:val="1"/></w:lvlOverride></w:num>'))
    continue_list_numbering(paragraph, num_id)
    return num_id


def continue_list_numbering(paragraph, num_id, ilvl=0):
    numpr = paragraph._p.get_or_add_pPr().get_or_add_numPr()
    numpr.get_or_add_ilvl().val = ilvl
    numpr.get_or_add_numId().val = num_id


# ---------------------------------------------------------------- 域 ----

def add_field(paragraph, instr: str, placeholder: str = "", dirty: bool = False,
              run_fmt: dict | None = None):
    """插入复杂域 { instr }，placeholder 为未更新前显示的缓存结果。返回结果 run。"""
    def _run():
        r = paragraph.add_run()
        if run_fmt:
            set_run_font(r, **run_fmt)
        return r

    r = _run()
    fld = OxmlElement("w:fldChar")
    fld.set(qn("w:fldCharType"), "begin")
    if dirty:
        fld.set(qn("w:dirty"), "true")
    r._r.append(fld)
    r = _run()
    it = OxmlElement("w:instrText")
    it.set(qn("xml:space"), "preserve")
    it.text = f" {instr} "
    r._r.append(it)
    r = _run()
    fld = OxmlElement("w:fldChar")
    fld.set(qn("w:fldCharType"), "separate")
    r._r.append(fld)
    result = _run()
    result.text = placeholder
    r = _run()
    fld = OxmlElement("w:fldChar")
    fld.set(qn("w:fldCharType"), "end")
    r._r.append(fld)
    return result


def _settings_insert(doc, el):
    """按 schema 顺序把 el 插入 settings.xml（放在 compat/rsids 等之前）。"""
    settings = doc.settings.element
    for tag in ("w:hdrShapeDefaults", "w:footnotePr", "w:endnotePr", "w:compat", "w:docVars", "w:rsids",
                "m:mathPr", "w:themeFontLang", "w:clrSchemeMapping", "w:shapeDefaults",
                "w:decimalSymbol", "w:listSeparator"):
        node = settings.find(qn(tag))
        if node is not None:
            node.addprevious(el)
            return el
    settings.append(el)
    return el


def enable_update_fields_on_open(doc):
    """Word 打开时弹窗“是否更新域”(TOC/页码/题注一起更新)。LibreOffice 忽略此设置。"""
    settings = doc.settings.element
    uf = settings.find(qn("w:updateFields"))
    if uf is None:
        uf = OxmlElement("w:updateFields")
        _settings_insert(doc, uf)
    uf.set(qn("w:val"), "true")


def add_toc(doc, levels="1-3", title="目  录", placeholder="（目录将在 Word 打开时自动更新；或右键→更新域）"):
    if title:
        p = doc.add_paragraph(style="TOC Heading")
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run(title)
        set_run_font(r, "黑体", size=16, bold=True, color="000000")
    p = doc.add_paragraph()
    p.paragraph_format.first_line_indent = Pt(0)
    add_field(p, f'TOC \\o "{levels}" \\h \\z \\u', placeholder, dirty=True)
    enable_update_fields_on_open(doc)
    return p


def add_page_number_footer(section, template="第 {PAGE} 页 / 共 {NUMPAGES} 页", size=9):
    footer = section.footer
    footer.is_linked_to_previous = False
    p = footer.paragraphs[0]
    p.text = ""
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    import re
    for token in re.split(r"(\{PAGE\}|\{NUMPAGES\}|\{SECTIONPAGES\})", template):
        if token in ("{PAGE}", "{NUMPAGES}", "{SECTIONPAGES}"):
            add_field(p, token[1:-1], "1", run_fmt=dict(size=size))
        elif token:
            set_run_font(p.add_run(token), size=size)
    return p


def set_header_text(section, text, size=9, align=WD_ALIGN_PARAGRAPH.CENTER, bottom_border=True):
    header = section.header
    header.is_linked_to_previous = False
    p = header.paragraphs[0]
    p.text = ""
    p.alignment = align
    set_run_font(p.add_run(text), size=size, color="595959")
    if bottom_border:
        pbdr = parse_xml(f'<w:pBdr {nsdecls("w")}><w:bottom w:val="single" w:sz="6" w:space="1" w:color="1F3864"/></w:pBdr>')
        insert_ordered(p._p.get_or_add_pPr(), pbdr)
    return p


def restart_page_numbering(section, start=1, fmt: str | None = None):
    sect_pr = section._sectPr
    pg = sect_pr.find(qn("w:pgNumType"))
    if pg is None:
        pg = OxmlElement("w:pgNumType")
        # pgNumType 需在 cols 之前
        cols = sect_pr.find(qn("w:cols"))
        (cols.addprevious(pg) if cols is not None else sect_pr.append(pg))
    pg.set(qn("w:start"), str(start))
    if fmt:
        pg.set(qn("w:fmt"), fmt)  # lowerRoman / upperRoman / decimal


# ------------------------------------------------------------- 封面 ----

def add_cover_page(doc, title, subtitle="", meta: Iterable[tuple[str, str]] = (), logo_path=None,
                   new_section_after=True):
    for _ in range(4):
        doc.add_paragraph()
    if logo_path:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.add_run().add_picture(logo_path, width=Cm(3))
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.first_line_indent = Pt(0)
    set_run_font(p.add_run(title), "黑体", size=30, bold=True, color="1F3864")
    if subtitle:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.first_line_indent = Pt(0)
        set_run_font(p.add_run(subtitle), "黑体", size=18, color="2F5597")
    for _ in range(8):
        doc.add_paragraph()
    for k, v in meta:
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.first_line_indent = Pt(0)
        set_run_font(p.add_run(f"{k}：{v}"), "宋体", size=14)
    if new_section_after:
        sec = doc.add_section(WD_SECTION.NEW_PAGE)
        return sec


# --------------------------------------------------- 图片 + 题注(SEQ) ----

class CaptionCounter:
    def __init__(self):
        self.counts = {}

    def next(self, label):
        self.counts[label] = self.counts.get(label, 0) + 1
        return self.counts[label]


def add_caption(doc, label, text, counter: CaptionCounter, chapter_style_ref: bool = False):
    """题注：图 1  说明 / 表 1  说明。chapter_style_ref=True 则为 图 2-1 形式(STYLEREF 1 + SEQ \\s 1)。"""
    n = counter.next(label)
    p = doc.add_paragraph(style="Caption")
    p.add_run(f"{label} ")
    if chapter_style_ref:
        add_field(p, "STYLEREF 1 \\s", "1")
        p.add_run("-")
        add_field(p, f"SEQ {label} \\* ARABIC \\s 1", str(n))
    else:
        add_field(p, f"SEQ {label} \\* ARABIC", str(n))
    p.add_run(f"  {text}")
    for r in p.runs:
        set_run_font(r, "黑体", size=10.5)
    if label == "图":  # 图注在图下方：不要与后文“同页”，否则可能把图注挤到下一页
        p.paragraph_format.keep_with_next = False
    return p


def add_figure(doc, image_path, caption, counter: CaptionCounter, width_cm=14.0):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.first_line_indent = Pt(0)
    p.paragraph_format.keep_with_next = True
    p.add_run().add_picture(image_path, width=Cm(width_cm))
    return add_caption(doc, "图", caption, counter)


# ------------------------------------------------------------- 表格 ----

def set_cell_shading(cell, fill_hex: str):
    tc_pr = cell._tc.get_or_add_tcPr()
    insert_ordered(tc_pr, parse_xml(f'<w:shd {nsdecls("w")} w:val="clear" w:color="auto" w:fill="{fill_hex}"/>'))


def set_cell_borders(cell, **edges):
    """edges: top/bottom/left/right = dict(val='single', sz=4, color='000000') 或 None(nil)。"""
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = tc_pr.find(qn("w:tcBorders"))
    if borders is None:
        borders = insert_ordered(tc_pr, OxmlElement("w:tcBorders"))
    for edge in ("top", "left", "bottom", "right"):
        if edge not in edges:
            continue
        spec = edges[edge]
        el = borders.find(qn(f"w:{edge}"))
        if el is None:
            el = OxmlElement(f"w:{edge}")
            borders.append(el)
        if spec is None:
            el.set(qn("w:val"), "nil")
        else:
            el.set(qn("w:val"), spec.get("val", "single"))
            el.set(qn("w:sz"), str(spec.get("sz", 4)))
            el.set(qn("w:space"), "0")
            el.set(qn("w:color"), spec.get("color", "000000"))


def set_table_borders(table, top=12, bottom=12, inside_h=None, inside_v=None, left=None, right=None,
                      color="000000"):
    """表级边框。参数为 sz(1/8 pt)；None 表示无边框。"""
    tbl_pr = table._tbl.tblPr
    parts = []
    for name, sz in (("top", top), ("left", left), ("bottom", bottom), ("right", right),
                     ("insideH", inside_h), ("insideV", inside_v)):
        if sz is None:
            parts.append(f'<w:{name} w:val="nil"/>')
        else:
            parts.append(f'<w:{name} w:val="single" w:sz="{sz}" w:space="0" w:color="{color}"/>')
    insert_ordered(tbl_pr, parse_xml(f'<w:tblBorders {nsdecls("w")}>{"".join(parts)}</w:tblBorders>'))


def set_col_widths(table, widths_cm: Sequence[float]):
    table.autofit = False
    tbl_pr = table._tbl.tblPr
    layout = tbl_pr.find(qn("w:tblLayout"))
    if layout is None:
        layout = insert_ordered(tbl_pr, OxmlElement("w:tblLayout"))
    layout.set(qn("w:type"), "fixed")
    for row in table.rows:
        for idx, w in enumerate(widths_cm):
            if idx < len(row.cells):
                row.cells[idx].width = Cm(w)
    for idx, w in enumerate(widths_cm):
        table.columns[idx].width = Cm(w)


def _format_cell_text(cell, text, bold=False, size=10.5, align=WD_ALIGN_PARAGRAPH.CENTER,
                      east_asia="宋体", color=None):
    cell.text = ""
    p = cell.paragraphs[0]
    p.alignment = align
    p.paragraph_format.first_line_indent = Pt(0)
    p.paragraph_format.line_spacing = 1.0
    p.paragraph_format.keep_with_next = True  # 短表不跨页（最后一行由 _release_last_row 解除）
    r = p.add_run(str(text))
    set_run_font(r, east_asia, size=size, bold=bold, color=color)
    cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER


def _release_last_row(table):
    for cell in table.rows[-1].cells:
        for p in cell.paragraphs:
            p.paragraph_format.keep_with_next = False
    for row in table.rows:  # 行内不拆分
        tr_pr = row._tr.get_or_add_trPr()
        if tr_pr.find(qn("w:cantSplit")) is None:
            insert_ordered(tr_pr, OxmlElement("w:cantSplit"))


def add_three_line_table(doc, header: Sequence[str], rows: Sequence[Sequence], widths_cm=None,
                         header_fill: str | None = None, font_size=10.5, first_col_left=True,
                         num_align_right=True):
    """学术/商业计划书常用三线表：顶线、底线粗，表头下细线，无竖线。"""
    table = doc.add_table(rows=1 + len(rows), cols=len(header))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    set_table_borders(table, top=12, bottom=12, inside_h=None, inside_v=None, left=None, right=None)
    for j, h in enumerate(header):
        c = table.cell(0, j)
        _format_cell_text(c, h, bold=True, size=font_size, east_asia="黑体")
        set_cell_borders(c, bottom=dict(sz=6))
        if header_fill:
            set_cell_shading(c, header_fill)
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            is_num = isinstance(v, (int, float)) or (isinstance(v, str) and v.replace(",", "").replace(".", "").replace("%", "").replace("-", "").isdigit())
            align = WD_ALIGN_PARAGRAPH.CENTER
            if j == 0 and first_col_left:
                align = WD_ALIGN_PARAGRAPH.LEFT
            elif is_num and num_align_right:
                align = WD_ALIGN_PARAGRAPH.RIGHT
            _format_cell_text(table.cell(i, j), v, size=font_size, align=align)
    _release_last_row(table)
    if widths_cm:
        set_col_widths(table, widths_cm)
    return table


def add_grid_table(doc, header, rows, widths_cm=None, header_fill="1F3864", header_font_color="FFFFFF",
                   zebra_fill="F2F2F2", font_size=10.5, border_color="BFBFBF"):
    """带底纹的全框线表格：深色表头 + 斑马纹。"""
    table = doc.add_table(rows=1 + len(rows), cols=len(header))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    set_table_borders(table, top=4, bottom=4, inside_h=4, inside_v=4, left=4, right=4, color=border_color)
    for j, h in enumerate(header):
        c = table.cell(0, j)
        _format_cell_text(c, h, bold=True, size=font_size, east_asia="黑体", color=header_font_color)
        set_cell_shading(c, header_fill)
    for i, row in enumerate(rows, start=1):
        for j, v in enumerate(row):
            c = table.cell(i, j)
            _format_cell_text(c, v, size=font_size,
                              align=WD_ALIGN_PARAGRAPH.LEFT if j == 0 else WD_ALIGN_PARAGRAPH.CENTER)
            if zebra_fill and i % 2 == 0:
                set_cell_shading(c, zebra_fill)
    _release_last_row(table)
    if widths_cm:
        set_col_widths(table, widths_cm)
    return table


def set_repeat_table_header(row):
    tr_pr = row._tr.get_or_add_trPr()
    el = OxmlElement("w:tblHeader")
    el.set(qn("w:val"), "true")
    insert_ordered(tr_pr, el)


# ----------------------------------------------------------- 分节 ----

def clear_page_number_restart(section):
    """add_section() 会深拷贝上一节的 sectPr，连 pgNumType(start=1) 一起复制 → 新节页码又从 1 开始。"""
    for pg in section._sectPr.findall(qn("w:pgNumType")):
        section._sectPr.remove(pg)


def add_landscape_section(doc):
    sec = doc.add_section(WD_SECTION.NEW_PAGE)
    clear_page_number_restart(sec)
    sec.orientation = WD_ORIENT.LANDSCAPE
    if sec.page_width < sec.page_height:
        sec.page_width, sec.page_height = sec.page_height, sec.page_width
    return sec


def add_portrait_section(doc):
    sec = doc.add_section(WD_SECTION.NEW_PAGE)
    clear_page_number_restart(sec)
    sec.orientation = WD_ORIENT.PORTRAIT
    if sec.page_width > sec.page_height:
        sec.page_width, sec.page_height = sec.page_height, sec.page_width
    return sec


# --------------------------------------------------------- 超链接 ----

def add_hyperlink(paragraph, text, url: str | None = None, anchor: str | None = None,
                  color="0563C1", underline=True, east_asia="宋体", size=None):
    hl = OxmlElement("w:hyperlink")
    if url:
        r_id = paragraph.part.relate_to(url, RT.HYPERLINK, is_external=True)
        hl.set(qn("r:id"), r_id)
    else:
        hl.set(qn("w:anchor"), anchor)
    hl.set(qn("w:history"), "1")
    r = OxmlElement("w:r")
    rpr = OxmlElement("w:rPr")
    set_fonts(rpr, east_asia)
    c = OxmlElement("w:color")
    c.set(qn("w:val"), color)
    rpr.append(c)
    if size:
        sz = OxmlElement("w:sz")
        sz.set(qn("w:val"), str(int(size * 2)))
        rpr.append(sz)
    if underline:
        u = OxmlElement("w:u")
        u.set(qn("w:val"), "single")
        rpr.append(u)
    r.append(rpr)
    t = OxmlElement("w:t")
    t.set(qn("xml:space"), "preserve")
    t.text = text
    r.append(t)
    hl.append(r)
    paragraph._p.append(hl)
    return hl


_bookmark_id = [100]


def add_bookmark(paragraph, name):
    _bookmark_id[0] += 1
    bid = str(_bookmark_id[0])
    start = OxmlElement("w:bookmarkStart")
    start.set(qn("w:id"), bid)
    start.set(qn("w:name"), name)
    end = OxmlElement("w:bookmarkEnd")
    end.set(qn("w:id"), bid)
    p = paragraph._p
    p.insert(1 if p.pPr is not None else 0, start)
    p.append(end)
    return name


# ----------------------------------------------------------- 脚注 ----

_FOOTNOTES_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    f'<w:footnotes {nsdecls("w", "r")}>'
    '<w:footnote w:type="separator" w:id="-1"><w:p><w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr>'
    '<w:r><w:separator/></w:r></w:p></w:footnote>'
    '<w:footnote w:type="continuationSeparator" w:id="0"><w:p><w:pPr><w:spacing w:after="0" w:line="240" w:lineRule="auto"/></w:pPr>'
    '<w:r><w:continuationSeparator/></w:r></w:p></w:footnote>'
    '</w:footnotes>'
)
FOOTNOTES_CT = "application/vnd.openxmlformats-officedocument.wordprocessingml.footnotes+xml"


class Footnotes:
    """python-docx 1.2.0 仍无脚注 API；本类手工创建 /word/footnotes.xml 部件。"""

    def __init__(self, doc):
        self.doc = doc
        self.part = XmlPart.load(PackURI("/word/footnotes.xml"), FOOTNOTES_CT,
                                 _FOOTNOTES_XML.encode("utf-8"), doc.part.package)
        doc.part.relate_to(self.part, RT.FOOTNOTES)
        self.next_id = 1
        styles = doc.styles
        if "Footnote Text" not in [s.name for s in styles]:
            st = styles.add_style("Footnote Text", WD_STYLE_TYPE.PARAGRAPH)
            st.base_style = styles["Normal"]
            set_style_font(st, "宋体", size=9)
            st.paragraph_format.line_spacing = 1.0
            st.paragraph_format.first_line_indent = Pt(0)
            st.element.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "0")
        # settings 里登记分隔符脚注（Word 生成文件的惯例，可选）
        settings = doc.settings.element
        if settings.find(qn("w:footnotePr")) is None:
            _settings_insert(doc, parse_xml(
                f'<w:footnotePr {nsdecls("w")}><w:footnote w:id="-1"/><w:footnote w:id="0"/></w:footnotePr>'))

    def add(self, paragraph, text):
        fid = self.next_id
        self.next_id += 1
        fn = parse_xml(
            f'<w:footnote {nsdecls("w")} w:id="{fid}"><w:p><w:pPr><w:pStyle w:val="FootnoteText"/></w:pPr>'
            f'<w:r><w:rPr><w:vertAlign w:val="superscript"/></w:rPr><w:footnoteRef/></w:r>'
            f'<w:r><w:t xml:space="preserve"> {text}</w:t></w:r></w:p></w:footnote>')
        self.part.element.append(fn)
        r = paragraph.add_run()
        rpr = r._r.get_or_add_rPr()
        va = OxmlElement("w:vertAlign")
        va.set(qn("w:val"), "superscript")
        rpr.append(va)
        ref = OxmlElement("w:footnoteReference")
        ref.set(qn("w:id"), str(fid))
        r._r.append(ref)
        return fid


# ----------------------------------------------------------- 段落 ----

def add_body_paragraph(doc, text, bold_spans: Sequence[str] = (), first_line_indent=True,
                       line_spacing=None, space_before=None, space_after=None, align=None):
    p = doc.add_paragraph()
    _add_inline(p, text, bold_spans)
    pf = p.paragraph_format
    if not first_line_indent:
        pf.first_line_indent = Pt(0)
        p._p.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "0")
    if line_spacing is not None:
        if isinstance(line_spacing, Pt):
            pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY
        pf.line_spacing = line_spacing
    if space_before is not None:
        pf.space_before = Pt(space_before)
    if space_after is not None:
        pf.space_after = Pt(space_after)
    if align is not None:
        p.alignment = align
    return p


def _add_inline(p, text, bold_spans=()):
    """把 text 中出现的 bold_spans 子串加粗，其余正常。"""
    if not bold_spans:
        set_run_font(p.add_run(text))
        return
    import re
    pattern = "(" + "|".join(re.escape(s) for s in bold_spans) + ")"
    for seg in re.split(pattern, text):
        if not seg:
            continue
        set_run_font(p.add_run(seg), bold=seg in bold_spans)


def add_quote_block(doc, text, fill="F2F2F2", bar_color="1F3864"):
    p = doc.add_paragraph()
    pf = p.paragraph_format
    pf.left_indent = Cm(0.8)
    pf.right_indent = Cm(0.8)
    pf.first_line_indent = Pt(0)
    p._p.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "0")
    pf.space_before = Pt(6)
    pf.space_after = Pt(6)
    ppr = p._p.get_or_add_pPr()
    insert_ordered(ppr, parse_xml(
        f'<w:pBdr {nsdecls("w")}><w:left w:val="single" w:sz="24" w:space="8" w:color="{bar_color}"/></w:pBdr>'))
    insert_ordered(ppr, parse_xml(f'<w:shd {nsdecls("w")} w:val="clear" w:color="auto" w:fill="{fill}"/>'))
    set_run_font(p.add_run(text), "楷体", size=11, color="404040")
    return p


def add_bullets(doc, items: Sequence[str], style="List Bullet"):
    out = []
    for it in items:
        p = doc.add_paragraph(style=style)
        set_run_font(p.add_run(it))
        out.append(p)
    return out


def add_numbered(doc, items: Sequence[str], restart=True, style="List Number"):
    out = []
    num_id = None
    for i, it in enumerate(items):
        p = doc.add_paragraph(style=style)
        set_run_font(p.add_run(it))
        if restart:
            if i == 0:
                num_id = restart_list_numbering(doc, p, style)
            else:
                continue_list_numbering(p, num_id)
        out.append(p)
    return out
