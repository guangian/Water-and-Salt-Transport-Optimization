# 04 Office 文档自动生成技术底座（python-docx / python-pptx / openpyxl / matplotlib / Markdown→docx）

| 项 | 内容 |
| --- | --- |
| 目标 | 为「自动生成商业竞赛材料（商业计划书 Word、路演 PPT、财务模型 Excel、商业图表）的 AI 系统」确定生成技术底座，并用可运行代码验证每一项结论 |
| 验证环境 | Ubuntu 24.04 (noble) x86_64，Python 3.12.3（venv：`/tmp/officetest/venv`），LibreOffice 24.2.7.2（apt，headless），pandoc 3.1.3，poppler-utils 24.02（pdftoppm/pdffonts/pdftotext），fonts-noto-cjk 1:20230817；无 Windows/MS Office，所有「Word/PowerPoint/Excel 打开后的行为」均标注为未在本机验证 |
| 库版本（`pip list` 实测） | python-docx 1.2.0 · python-pptx 1.0.2 · openpyxl 3.1.5 · matplotlib 3.11.2 · pandas 3.1.0 · pillow 13.0.0 · numpy 2.5.4 · lxml 6.1.3 · markdown-it-py 4.2.0 · formulas 1.3.4 · docxtpl 0.20.2 · PyMuPDF 2.0（MuPDF 1.28.5）· fonttools 4.67.0 · numpy-financial 1.0.0 · python3-uno（系统 python3，随 LibreOffice） |
| 验证方法 | 全部示例在 `/tmp/officetest` 实际运行；docx/pptx/xlsx 用 LibreOffice 转 PDF、`pdftoppm` 转 PNG 后逐页目视检查中文、编号、目录、图表；xlsx 公式结果用 numpy-financial 复算对照；一键复现脚本 `run_all.sh` 全流程 12.8 s 跑通（日志见附录 A） |
| 标记约定 | **已验证** = 本机实际运行并检查输出；**未验证** = 无对应软件（Windows/Word/PowerPoint/Excel），依据官方文档/规范/社区经验给出，需在 Windows 上复核 |
| 交付物位置 | 报告：本文件；全部源码 + 样图 + 运行日志：`/home/ubuntu/bizgold/docs/research/04_assets/`（`code/` 下 13 个文件可直接复用）；生成的 docx/pptx/xlsx/pdf 样例在 `/tmp/officetest/out/` |

目录

0. 结论速览（给架构设计者的 1 页纸）
1. 环境搭建与版本（已验证）
2. python-docx：商业计划书 Word 生成
3. python-pptx：路演 PPT 生成
4. openpyxl：财务模型 Excel 生成
5. matplotlib：商业图表
6. Markdown→docx 中间层方案对比与推荐流水线
7. 各库版本、坑点与规避（汇总）
8. 附录：一键复现日志、文件清单

---

## 0. 结论速览

1. **技术选型可行，四个库 + LibreOffice headless 即可覆盖全部需求**（已验证）：python-docx 1.2.0 生成计划书（含多级编号、TOC 域、页码域、SEQ 题注、三线表、横向页、脚注、超链接、批注）；python-pptx 1.0.2 生成 16:9 路演 PPT（10 类页面、原生柱/折/饼/组合图、备注、模板占位符替换）；openpyxl 3.1.5 生成 9 表财务模型（跨表公式、IF/NPV/IRR/SUMPRODUCT/COUNTIF、命名区域、数据验证、原生图表）；matplotlib 3.11.2 生成 11 类 300dpi 中文商业图。
2. **Linux 服务器上的中文字体是第一道坎**，解法（已验证）：`apt install fonts-noto-cjk`；matplotlib 必须**重建字体缓存**（装字体前已 import 过 matplotlib 则缓存是旧的，`fontManager` 看不到 Noto CJK）；LibreOffice 转 PDF 时用 `FONTCONFIG_FILE` 指向一份别名配置把 `宋体/黑体/微软雅黑/楷体` 映射到 `Noto Serif/Sans CJK SC`，否则 LibreOffice 把宋体和黑体都替换成同一种无衬线字体（pdffonts 实测）。
3. **docx 的目录（TOC）只能靠「域」**：python-docx 写 TOC 域 + `w:updateFields=true`（Word 打开时提示更新，未验证 Word 行为），服务端出 PDF 必须用 LibreOffice **UNO** 脚本 `lo_convert.py` 调 `DocumentIndexes.update()`（已验证：PDF 目录带真实页码；`soffice --convert-to pdf` 命令行**不会**更新目录，只显示占位文字）。
4. **openpyxl 写的公式没有缓存值**（`data_only=True` 读出 None，已验证），三种对策全部跑通：LibreOffice 重算（0.9 s，NPV/IRR 与 numpy-financial 一致）、`formulas` 库纯 Python 计算（0.8 s，支持 NPV/IRR/IF/SUMPRODUCT/COUNTIF/中文命名区域）、复制公式版并用计算值覆盖得到「数值版」（样式/图表/命名区域全部保留，pandas 可直接读）。交付给评委的文件另加 `fullCalcOnLoad=True` 让 Excel 打开即重算。
5. **Markdown→docx 推荐自写渲染器**（markdown-it-py 解析 + 本报告的 `docx_helpers`），而不是 pandoc/docxtpl：pandoc 出的 docx 没有中文字体设置、没有 SEQ 题注域、表格无样式、TOC 不带页码（已验证 XML）；docxtpl 适合固定版式的报名表/封面，不适合长文；自写渲染器 272 行即可支持标题层级+自动编号、段落首行缩进、粗斜体/代码/链接、嵌套列表、三线表+题注、图片+题注、引用块、分页标记、代码块（已验证）。
6. **需要在 Windows 上复核的点**（未验证）：Word 打开时 `updateFields` 弹窗更新目录；`docx2pdf`/`win32com` 转 PDF；PowerPoint 对手工拼的组合图 XML 的接受度（LibreOffice 渲染正常，XML 顺序按 OOXML 规范）；matplotlib 在 Windows 上要用英文字体名 `Microsoft YaHei`/`SimHei`。

---

## 1. 环境搭建与版本（已验证）

### 1.1 安装命令（Ubuntu 24.04，实际执行）

```bash
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
    python3.12-venv fonts-noto-cjk pandoc poppler-utils \
    libreoffice-writer libreoffice-calc libreoffice-impress python3-uno
mkdir -p /tmp/officetest && cd /tmp/officetest
python3 -m venv venv && ./venv/bin/pip install -U pip
./venv/bin/pip install python-docx python-pptx openpyxl matplotlib pandas pillow \
    markdown-it-py formulas pymupdf docxtpl
```

要点：
- `python3-uno` 必须用 **apt** 装且只能被**系统 python3** 导入（`/usr/bin/python3 -c "import uno"` 已验证可用；venv 里 `import uno` 不可用）。所以 `lo_convert.py` 用系统 python 跑，其余用 venv。
- `--no-install-recommends` 装 LibreOffice 三件套约 1 分钟，不含 Java，`--headless` 转换 docx/pptx/xlsx 均正常（已验证）。
- 字体：装 `fonts-noto-cjk` 后 `fc-list :lang=zh family` 出现 `Noto Sans/Serif/Sans Mono CJK SC/TC/JP/KR/HK`（之前只有 WenQuanYi Micro Hei 与 Droid Sans Fallback）。**但** `fc-match 宋体` 仍返回 `Noto Sans`（拉丁字体！）——fontconfig 不认识中文字体名，需要第 2.11 节的别名配置。

### 1.2 实测版本

| 组件 | 版本 | 备注 |
| --- | --- | --- |
| Python | 3.12.3 | venv |
| python-docx | 1.2.0 | 新增 `add_comment()`；仍无脚注 API |
| python-pptx | 1.0.2 | 无组合图/主题 API，需改 XML |
| openpyxl | 3.1.5 | `defined_names` 为 dict 风格 |
| matplotlib | 3.11.2 | 能枚举 .ttc 全部字体面 |
| pandas / numpy | 3.1.0 / 2.5.4 | |
| pillow | 13.0.0 | 拼接预览图 |
| lxml | 6.1.3 | 改 XML |
| markdown-it-py | 4.2.0 | Markdown 解析 |
| formulas | 1.3.4 | 纯 Python 公式引擎（依赖 schedula 1.6.15、XlsxWriter 3.3.0） |
| docxtpl | 0.20.2 | Jinja2 模板 |
| PyMuPDF | 2.0 (MuPDF 1.28.5) | PDF→PNG；`import fitz` 已弃用，用 `import pymupdf` |
| LibreOffice | 24.2.7.2 | headless 转换 + UNO |
| pandoc | 3.1.3 | 对比用 |
| poppler-utils | 24.02 | pdftoppm/pdffonts/pdftotext |

---

## 2. python-docx：商业计划书 Word 生成

全部功能封装在 `docx_helpers.py`（771 行，附 2.13），演示脚本 `01_docx_demo.py` 生成 7 页 A4 样例（封面、目录、正文 3 章、横向页、参考文献），经 `lo_convert.py` 更新目录后转 PDF，逐页检查通过。样图：`04_assets/docx_page_正文.png`、`docx_page_示例.png`、`docx_目录与横向页.png`。

### 2.1 中文字体（w:eastAsia）——已验证

结论：
1. `run.font.name = "Times New Roman"` 只写 `w:ascii`/`w:hAnsi`；中文必须另写 `w:eastAsia`：`rPr.rFonts.set(qn("w:eastAsia"), "宋体")`。
2. **默认模板的 Heading 1–9 样式带主题字体属性** `w:asciiTheme="majorHAnsi" w:eastAsiaTheme="majorEastAsia"…`（实测 `Heading 1` 的 rFonts 只有 4 个 *Theme 属性）。主题属性优先级高于显式字体名，不删掉它们，标题字体设置在 Word 里不生效。`set_fonts()` 统一 `attrib.pop()` 四个主题属性。
3. docDefaults（`w:docDefaults/w:rPrDefault/w:rPr/w:rFonts`）同样是主题字体，应一并改成 `宋体 + Times New Roman`，并把 `w:lang w:eastAsia` 设为 `zh-CN`（默认是 `en-US`）。
4. 默认模板页面是 **Letter（8.5×11 in）**，中文材料要改 A4（`Mm(210) × Mm(297)`），已在 `setup_chinese_document()` 处理。
5. 正文字体方案（样式级一次设置，不需要每个 run 设）：正文 宋体 小四 1.5 倍行距 首行缩进 2 字符；标题 黑体 加粗；题注 黑体 五号；引用块 楷体；西文一律 Times New Roman。
6. LibreOffice 渲染对照（pdffonts 实测，见 2.11）：加别名后 `宋体→NotoSerifCJKsc`、`黑体→NotoSansCJKsc`、`Times New Roman→LiberationSerif`，嵌入 PDF，无乱码。

核心代码（节选自 `docx_helpers.py`）：

```python
THEME_FONT_ATTRS = ("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme")

def set_fonts(rpr, east_asia="宋体", latin="Times New Roman"):
    rfonts = rpr.find(qn("w:rFonts"))
    if rfonts is None:
        rfonts = OxmlElement("w:rFonts"); rpr.insert(0, rfonts)   # rFonts 必须是 rPr 第一个子元素
    for a in THEME_FONT_ATTRS:
        rfonts.attrib.pop(qn(a), None)                              # 主题字体属性优先于显式字体名
    rfonts.set(qn("w:ascii"), latin); rfonts.set(qn("w:hAnsi"), latin)
    rfonts.set(qn("w:eastAsia"), east_asia); rfonts.set(qn("w:cs"), latin)

# 样式级：set_style_font(doc.styles["Normal"], "宋体", "Times New Roman", size=12)
# run 级：  set_run_font(run, "黑体", size=14, bold=True, color="1F3864")
# 首行缩进 2 字符（随字号自适应，优于 Pt(24)）：
normal.element.get_or_add_pPr().get_or_add_ind().set(qn("w:firstLineChars"), "200")
```

### 2.2 标题样式与多级编号——已验证

- python-docx 没有编号 API；做法：在 `numbering.xml` 新增一个 `w:abstractNum`（`multiLevelType=multilevel`，每级 `w:lvlText` = `%1` / `%1.%2` / `%1.%2.%3`，`w:pStyle` 指向 `Heading1..3`，`w:suff=space` 用空格而非制表符分隔编号与标题），再新增 `w:num` 指向它，最后把 `Heading n` 样式的 `pPr/numPr` 指向该 numId（`ppr.get_or_add_numPr()` 由 python-docx 保证元素顺序）。
- 坑 1：`w:abstractNum` 必须排在所有 `w:num` 之前，所以用 `existing[-1].addnext(abstract)` 而不是 `append`。
- 坑 2：**`TOC Heading` 样式基于 `Heading 1`，会继承编号**，首次生成时目录标题变成「1 目 录」、正文从 2 开始编号（实测截图）。修正：给 `TOC Heading` 写 `numPr numId=0`（0 = 不编号）。
- 坑 3：`add_heading(text, level)` 在 python-docx 中 level 0 = `Title` 样式；正文章节用 1–3。
- 效果：PDF 中为「1 执行摘要 / 2 项目概述 / 2.1 行业痛点 …」，目录条目同样带编号。

```python
def add_heading_numbering(doc, levels=3, fmt=None, east_asia="黑体", latin="Times New Roman"):
    numbering = doc.part.numbering_part.element
    fmt = fmt or ["%1", "%1.%2", "%1.%2.%3", "%1.%2.%3.%4"]
    abs_id = max([int(a.get(qn("w:abstractNumId"))) for a in numbering.findall(qn("w:abstractNum"))], default=-1) + 1
    lvls = "".join(
        f'<w:lvl w:ilvl="{i}"><w:start w:val="1"/><w:numFmt w:val="decimal"/><w:pStyle w:val="Heading{i+1}"/>'
        f'<w:suff w:val="space"/><w:lvlText w:val="{fmt[i]}"/><w:lvlJc w:val="left"/>'
        f'<w:pPr><w:ind w:left="0" w:firstLine="0"/></w:pPr>'
        f'<w:rPr><w:rFonts w:ascii="{latin}" w:hAnsi="{latin}" w:eastAsia="{east_asia}"/></w:rPr></w:lvl>'
        for i in range(levels))
    abstract = parse_xml(f'<w:abstractNum {nsdecls("w")} w:abstractNumId="{abs_id}">'
                         f'<w:multiLevelType w:val="multilevel"/>{lvls}</w:abstractNum>')
    numbering.findall(qn("w:abstractNum"))[-1].addnext(abstract)          # abstractNum 必须在 num 之前
    num_id = max([int(n.get(qn("w:numId"))) for n in numbering.findall(qn("w:num"))], default=0) + 1
    numbering.append(parse_xml(f'<w:num {nsdecls("w")} w:numId="{num_id}"><w:abstractNumId w:val="{abs_id}"/></w:num>'))
    for i in range(levels):
        numpr = doc.styles[f"Heading {i+1}"].element.get_or_add_pPr().get_or_add_numPr()
        numpr.get_or_add_ilvl().val = i; numpr.get_or_add_numId().val = num_id
    toc_numpr = doc.styles["TOC Heading"].element.get_or_add_pPr().get_or_add_numPr()   # 目录标题不编号
    toc_numpr.get_or_add_ilvl().val = 0; toc_numpr.get_or_add_numId().val = 0
    return num_id
```

### 2.3 自动目录（TOC）——三种方案的实测结论

| 方案 | 做法 | 实测结果 | 适用 |
| --- | --- | --- | --- |
| A. TOC 域 + `updateFields` | 写复杂域 `{ TOC \o "1-3" \h \z \u }`（`w:fldChar begin/separate/end`，`w:dirty=true`），`settings.xml` 加 `<w:updateFields w:val="true"/>` | XML 已验证生成、顺序正确（`updateFields` 位于 `savePreviewPicture` 之后、`footnotePr/compat` 之前）。**Word 打开时弹「是否更新域」—未在本机验证**；LibreOffice 命令行转换**忽略**它，PDF 里只显示占位文字 | 交付 .docx 给评委/用户在 Word 中打开 |
| B. 服务端 UNO 更新 | `lo_convert.py`：UNO 连接 headless LibreOffice → `doc.getDocumentIndexes().getByIndex(i).update()` → `getTextFields().refresh()` → 再 update 一次（目录自身会改变页码）→ 另存 docx + 导出 PDF | **已验证**：PDF 第 2 页为完整目录「1 执行摘要……1 / 2 项目概述……1 / … / 6 参考文献……5」（pdftotext 输出）；另存的 `*.updated.docx` 里目录在 `<w:sdt>` 内、条目样式 `TOC1/TOC2`、每条带内部超链接 | 服务端出 PDF / 出"已更新目录"的 docx |
| C. 手工目录 | 用 python 按标题列表生成静态段落（可加书签 + `w:hyperlink w:anchor` 内部链接，页码用 `PAGEREF 书签 \h` 域） | 未单独实现；页码仍是域，仍需更新，故不如 B | 不推荐，仅在禁用域时用 |

补充（已验证）：
- `soffice --headless --convert-to pdf` 单次约 1.3–1.6 s（含进程启动）；UNO 方式 1.4–1.6 s，多个文件可复用同一进程更快。
- LibreOffice 另存的 docx 把目录放进 `w:sdt`，**python-docx 的 `doc.paragraphs` 不遍历 sdt 内容**（实测 `paragraphs` 中只剩「目  录」一行），要读目录条目需直接查 XML。
- `NUMPAGES` 统计全文页数（含封面/目录）。本例页脚「第 1 页 / 共 7 页」从正文节开始计数但总数为 7。要「共 N 页」只算正文，Word 可用嵌套域 `{ = {NUMPAGES} - 2 }`（未验证 LibreOffice 支持），或页脚放到独立末节用 `SECTIONPAGES`。

```python
def add_field(paragraph, instr, placeholder="", dirty=False, run_fmt=None):
    """复杂域 { instr }，placeholder 为更新前显示的缓存结果"""
    def _run():
        r = paragraph.add_run()
        if run_fmt: set_run_font(r, **run_fmt)
        return r
    r = _run(); fld = OxmlElement("w:fldChar"); fld.set(qn("w:fldCharType"), "begin")
    if dirty: fld.set(qn("w:dirty"), "true")
    r._r.append(fld)
    r = _run(); it = OxmlElement("w:instrText"); it.set(qn("xml:space"), "preserve"); it.text = f" {instr} "; r._r.append(it)
    r = _run(); fld = OxmlElement("w:fldChar"); fld.set(qn("w:fldCharType"), "separate"); r._r.append(fld)
    result = _run(); result.text = placeholder
    r = _run(); fld = OxmlElement("w:fldChar"); fld.set(qn("w:fldCharType"), "end"); r._r.append(fld)
    return result

def add_toc(doc, levels="1-3", title="目  录", placeholder="（目录将在 Word 打开时自动更新；或右键→更新域）"):
    p = doc.add_paragraph(style="TOC Heading"); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    set_run_font(p.add_run(title), "黑体", size=16, bold=True, color="000000")
    p = doc.add_paragraph(); p.paragraph_format.first_line_indent = Pt(0)
    add_field(p, f'TOC \\o "{levels}" \\h \\z \\u', placeholder, dirty=True)
    enable_update_fields_on_open(doc)          # settings.xml: <w:updateFields w:val="true"/>
```

UNO 更新目录（节选自 `lo_convert.py`，完整见 2.13）：

```python
doc = desktop.loadComponentFromURL(uno.systemPathToFileUrl(path), "_blank", 0, (_prop("Hidden", True),))
idx = doc.getDocumentIndexes()
for i in range(idx.getCount()): idx.getByIndex(i).update()   # 目录
doc.getTextFields().refresh()                                  # PAGE/NUMPAGES/SEQ
for i in range(idx.getCount()): idx.getByIndex(i).update()   # 页码变化后再来一次
doc.storeToURL(uno.systemPathToFileUrl(out_docx), (_prop("FilterName", "MS Word 2007 XML"),))
doc.storeToURL(uno.systemPathToFileUrl(out_pdf), (_prop("FilterName", "writer_pdf_Export"),))
```

### 2.4 页眉页脚、页码域、分节、封面——已验证

- 封面 + 目录放在第 1 节（无页眉页脚），`doc.add_section()` 开第 2 节放正文：`footer.is_linked_to_previous = False` 后写「第 {PAGE} 页 / 共 {NUMPAGES} 页」（PAGE/NUMPAGES 为复杂域），`restart_page_numbering(sec, start=1)` 在 `sectPr` 插入 `<w:pgNumType w:start="1"/>`（必须放在 `w:cols` 之前）。
- 页眉文字 + 下边框线：段落 `w:pBdr/w:bottom`。
- **坑（实测踩到）：`doc.add_section()` 会深拷贝上一节的 `sectPr`，连 `pgNumType start=1` 一起复制**——后续横向节/纵向节页码又从 1 开始；LibreOffice 为维持奇偶页还自动插入空白页，导致 `NUMPAGES`=9 而 PDF 只有 7 页。修正：新节创建后 `clear_page_number_restart(sec)` 删除 `pgNumType`。
- 封面：`add_cover_page()` 用空段落垂直留白 + 居中大字（黑体 30 磅）+ 元信息行（团队/指导教师/日期），可选 logo；封面后 `add_page_break()` 再放目录。

### 2.5 图片插入与题注（SEQ 域）——已验证

- `paragraph.add_run().add_picture(path, width=Cm(14))` 居中；题注段落用 `Caption` 样式：`图 ` + 域 `SEQ 图 \* ARABIC` + 说明。Python 侧维护 `CaptionCounter` 把正确序号写进域的缓存结果，这样即使不更新域也显示正确；LibreOffice 转 PDF 时 SEQ 会被重新计算（PDF 中「图 1」「图 2」「表 1~3」正确）。
- **章节式编号（图 2-1）在 LibreOffice 下不可用（已验证）**：`{ STYLEREF 1 \s }-{ SEQ 图 \* ARABIC \s 1 }` 经 LibreOffice 渲染成「图 市场分析-1」「图 财务预测-2」（STYLEREF 取了标题文字而非编号，`\s 1` 不重置）。Word 中应正常（未验证）。跨平台稳妥做法：顺序编号「图 1」，或由 Python 计算章号写死文本「图 2-1」。
- 图与图注同页：图片段落 `keep_with_next=True`；图注段落 `keep_with_next=False`（否则图注会被"粘"到下一段而与图分离）。

```python
def add_caption(doc, label, text, counter, chapter_style_ref=False):
    n = counter.next(label)
    p = doc.add_paragraph(style="Caption"); p.add_run(f"{label} ")
    if chapter_style_ref:   # Word 专用：图 2-1
        add_field(p, "STYLEREF 1 \\s", "1"); p.add_run("-"); add_field(p, f"SEQ {label} \\* ARABIC \\s 1", str(n))
    else:                   # 跨平台：图 1
        add_field(p, f"SEQ {label} \\* ARABIC", str(n))
    p.add_run(f"  {text}")
    for r in p.runs: set_run_font(r, "黑体", size=10.5)
    if label == "图": p.paragraph_format.keep_with_next = False
    return p
```

### 2.6 表格：三线表 / 底纹表 / 列宽 / 合并 / 跨页——已验证

- 三线表：表级 `w:tblBorders` 顶线、底线 `sz=12`（1.5 磅），`insideH/insideV/left/right = nil`；表头行每个单元格 `w:tcBorders/w:bottom sz=6`；表头底纹 `w:shd fill="D9E2F3"`；数字列右对齐，首列左对齐。
- 底纹全框线表：表头深色（`1F3864`）白字 + 斑马纹（`F2F2F2`）+ 浅灰框线（`BFBFBF` sz=4）。
- 列宽：`table.autofit=False` + `w:tblLayout type=fixed` + **每个单元格** `cell.width`（只设 `column.width` 在 Word 中不生效）。
- 合并：`t.cell(4,1).merge(t.cell(4,3))`，合并后重新写文本与格式。
- 跨页：单元格段落 `keep_with_next=True`（最后一行解除）+ `w:cantSplit` 行内不拆分；`w:tblHeader` 让表头跨页重复。实测短表不再被拆到两页。
- **XML 子元素顺序**：手工 `append` 到 `tcPr/tblPr/trPr/pPr` 的元素（shd、tcBorders、tblBorders、tblLayout、cantSplit、pBdr）会排在 python-docx 已写入的 `vAlign/tblLook/jc/spacing` 之后，违反 OOXML 顺序（Word 可能报"无法读取的内容"）。`insert_ordered(parent, el)` 按 schema 顺序表插入；自检脚本扫描 `document.xml/header/footer/styles/footnotes` 的四类元素，结果 **misordered = 0**。

### 2.7 分页与分节（横向页）——已验证

`doc.add_page_break()` 分页；`add_landscape_section()`：`add_section(NEW_PAGE)` → 清除继承的 `pgNumType` → `orientation=LANDSCAPE` 并**交换宽高**（python-docx 只改 orientation 标记不会换纸张尺寸）；之后 `add_portrait_section()` 换回。页眉页脚默认 `is_linked_to_previous=True` 自动沿用。PDF 第 6 页为横向宽表、页码连续（第 4 页 / 共 7 页）。

### 2.8 段落格式——已验证

| 需求 | 代码 |
| --- | --- |
| 1.5 倍行距 | `pf.line_spacing = 1.5` |
| 固定值 20 磅 | `pf.line_spacing_rule = WD_LINE_SPACING.EXACTLY; pf.line_spacing = Pt(20)` |
| 首行缩进 2 字符 | `pPr.ind.set(qn("w:firstLineChars"), "200")`（Word/LibreOffice 都识别；`first_line_indent=Pt(24)` 只在小四时恰好等于 2 字） |
| 段前/段后 | `pf.space_before = Pt(6); pf.space_after = Pt(6)` |
| 两端对齐 | `p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY` |
| 标题与下文同页 | `pf.keep_with_next = True` |
| 不缩进（标题/题注/列表/封面） | 样式与段落上 `firstLineChars=0`（否则继承 Normal 的缩进） |
| 引用块 | 左侧 3 磅色条 `w:pBdr/w:left sz=24` + 底纹 `w:shd fill=F2F2F2` + 左右缩进 0.8 cm + 楷体 |

### 2.9 脚注 / 引用 / 超链接 / 书签 / 批注——已验证

- **脚注**：python-docx 1.2.0 仍无 API。`Footnotes` 类手工创建部件 `/word/footnotes.xml`（`XmlPart.load(PackURI("/word/footnotes.xml"), 内容类型, 初始 XML, package)` + `doc.part.relate_to(part, RT.FOOTNOTES)`），初始 XML 含 `separator`(id=-1) 与 `continuationSeparator`(id=0)；每次 `add()` 追加 `<w:footnote w:id=n>` 并在正文 run 插入 `<w:footnoteReference w:id=n/>`（上标）；`settings.xml` 登记 `w:footnotePr`。保存后 `[Content_Types].xml` 自动含 `footnotes+xml` Override，PDF 第 3 页底部出现「1 艾瑞咨询《2025 年中国智慧养老行业研究报告》…」。
- **参考文献/引用**：正文用脚注或「[1]」上标 + 文末「参考文献」编号列表（`add_numbered`，从 1 重新编号）。
- **超链接**：外部 `paragraph.part.relate_to(url, RT.HYPERLINK, is_external=True)` 得 rId → `<w:hyperlink r:id>`；内部 `<w:hyperlink w:anchor="书签名">` + `add_bookmark(heading_paragraph, "exec_summary")`（`bookmarkStart` 插在 `pPr` 之后）。默认模板**没有 `Hyperlink` 字符样式**，颜色/下划线直接写在 run 上。
- **批注**：`doc.add_comment(runs=p.runs[0], text="请补充官网上线时间。", author="评审老师", initials="PS")`（1.2.0 新 API，生成 `word/comments.xml`），可用于「AI 评审意见」回写到计划书。

### 2.10 项目符号 / 编号列表——已验证

- 默认模板自带 `List Bullet/List Bullet 2/3`、`List Number/List Number 2/3` 样式，直接 `add_paragraph(style=...)`。
- **坑：`List Number` 全文共用一个 numId，第二个列表会从上一个列表末尾继续编号**。`restart_list_numbering()` 新建 `w:num`（引用同一 abstractNum + `w:lvlOverride/w:startOverride=1`）并写到首段 `numPr`，后续段 `continue_list_numbering()` 使用同一 numId。
- 列表段落要把 `firstLineChars` 设 0（否则继承正文的 2 字符缩进，项目符号位置错乱）。

### 2.11 docx→pdf 转换方案

| 平台 | 方案 | 状态 | 要点 |
| --- | --- | --- | --- |
| Linux 服务器 | `soffice --headless --convert-to pdf --outdir OUT file.docx` | **已验证**（1.3 s/文件） | 不更新目录域；需 `FONTCONFIG_FILE` 字体别名；并发时用 `-env:UserInstallation=file:///tmp/lo_profile_N` 隔离配置目录 |
| Linux 服务器 | UNO（`lo_convert.py`，系统 python3 + python3-uno） | **已验证** | 更新目录/域后导出；同时产出 `*.updated.docx`；xlsx 可 `calculateAll()` 重算；pptx 直接导出 |
| Windows | `pip install docx2pdf` → `docx2pdf.convert("a.docx", "a.pdf")`（内部调用 Word COM） | 未验证 | 需本机装 Word；不更新目录 |
| Windows | pywin32：`word = win32com.client.Dispatch("Word.Application"); d = word.Documents.Open(path); d.TablesOfContents(1).Update(); d.Fields.Update(); d.ExportAsFixedFormat(pdf, 17); d.Close(False)` | 未验证 | 能更新目录再导出；路径用绝对路径，中文路径 OK（Python 3 str） |

字体别名配置 `fonts.conf`（已验证：`FONTCONFIG_FILE=/tmp/officetest/fonts.conf fc-match 宋体` → `Noto Serif CJK SC`，`pdffonts` 显示 PDF 内嵌 `NotoSerifCJKsc-Regular/Bold`、`NotoSansCJKsc-Regular/Bold`、`LiberationSerif`）。不配别名时 LibreOffice 用自带替换表把宋体与黑体都替成 `NotoSansCJKsc`，正文失去衬线体：

```xml
<?xml version="1.0"?>
<!DOCTYPE fontconfig SYSTEM "fonts.dtd">
<!-- Linux 服务器上把 Windows 中文字体名映射到 Noto CJK，供 LibreOffice 转 PDF 时使用。
     用法：FONTCONFIG_FILE=/tmp/officetest/fonts.conf soffice (headless) ...
     注意 XML 注释内不能出现两个连续连字符 -->
<fontconfig>
  <include ignore_missing="yes">/etc/fonts/fonts.conf</include>
  <alias binding="same"><family>宋体</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>SimSun</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>新宋体</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>仿宋</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>FangSong</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>楷体</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>KaiTi</family><prefer><family>Noto Serif CJK SC</family></prefer></alias>
  <alias binding="same"><family>黑体</family><prefer><family>Noto Sans CJK SC</family></prefer></alias>
  <alias binding="same"><family>SimHei</family><prefer><family>Noto Sans CJK SC</family></prefer></alias>
  <alias binding="same"><family>微软雅黑</family><prefer><family>Noto Sans CJK SC</family></prefer></alias>
  <alias binding="same"><family>Microsoft YaHei</family><prefer><family>Noto Sans CJK SC</family></prefer></alias>
  <alias binding="same"><family>等线</family><prefer><family>Noto Sans CJK SC</family></prefer></alias>
  <alias binding="same"><family>DengXian</family><prefer><family>Noto Sans CJK SC</family></prefer></alias>
  <alias binding="same"><family>Times New Roman</family><prefer><family>Liberation Serif</family><family>Tinos</family></prefer></alias>
  <alias binding="same"><family>Arial</family><prefer><family>Liberation Sans</family><family>Arimo</family></prefer></alias>
  <alias binding="same"><family>Calibri</family><prefer><family>Carlito</family></prefer></alias>
</fontconfig>
```

注意：该 XML 的注释里不能出现 `--`（第一次写成 `soffice --headless` 导致 fontconfig 报 `not well-formed`，已踩坑）。

### 2.12 生成结果校验（已验证）

- `[Content_Types].xml` 含 footnotes/comments Override；`document.xml.rels` 含 footnotes 与外部超链接关系；`numbering.xml` 中 abstractNum 全部位于 num 之前；`settings.xml` 顺序 `… savePreviewPicture, updateFields, footnotePr, compat …`。
- 四类属性容器子元素顺序自检：0 处错位。
- LibreOffice 打开无告警，PDF 7 页 A4，字体全部嵌入。

### 2.13 代码

`docx_helpers.py`（完整，已验证）：

```python
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
```

`01_docx_demo.py`（生成 7 页样例）：

```python
# -*- coding: utf-8 -*-
"""python-docx 1.2.0 全功能演示：生成 out/01_商业计划书示例.docx"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from docx import Document  # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH  # noqa: E402
from docx.shared import Pt  # noqa: E402

import docx_helpers as H  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
CHART1 = os.path.join(OUT, "charts", "01_市场规模堆叠柱.png")
CHART2 = os.path.join(OUT, "charts", "08_现金流组合.png")

doc = Document()
H.setup_chinese_document(doc, body_east_asia="宋体", heading_east_asia="黑体", latin="Times New Roman", body_size=12)
H.add_heading_numbering(doc, levels=3)
fn = H.Footnotes(doc)
cap = H.CaptionCounter()

# ---------- 第 1 节：封面 + 目录（无页码） ----------
sec_cover = doc.sections[0]
H.add_cover_page(doc, "智守康联：居家养老 AI 健康监测平台", "商业计划书",
                 meta=[("参赛团队", "智守康联团队"), ("指导教师", "张三 教授"), ("日期", "2026 年 10 月")],
                 new_section_after=False)
doc.add_page_break()
H.add_toc(doc, levels="1-3")

# ---------- 第 2 节：正文（页码从 1 开始，页眉页脚） ----------
sec_body = doc.add_section()
H.set_header_text(sec_body, "智守康联商业计划书  |  2026 全国高校商业精英挑战赛")
H.add_page_number_footer(sec_body, "第 {PAGE} 页 / 共 {NUMPAGES} 页")
H.restart_page_numbering(sec_body, start=1)

h = doc.add_heading("执行摘要", level=1)
H.add_bookmark(h, "exec_summary")
H.add_body_paragraph(doc, "智守康联聚焦独居老人的居家健康风险，提供“硬件 + 订阅服务”的一体化解决方案。"
                     "2026 年目标覆盖 3 万户家庭，实现营业收入 1,800 万元，净利润 260 万元。",
                     bold_spans=["3 万户家庭", "1,800 万元"])
p = H.add_body_paragraph(doc, "本方案中的市场数据主要来自艾瑞咨询与国家统计局公开报告")
fn.add(p, "艾瑞咨询《2025 年中国智慧养老行业研究报告》，2025 年 6 月。")
p.add_run("，行业规模测算采用自下而上法。")
H.set_run_font(p.runs[-1])

doc.add_heading("项目概述", level=1)
doc.add_heading("行业痛点", level=2)
H.add_bullets(doc, ["独居老人意外跌倒后平均 4.5 小时才被发现", "现有穿戴设备佩戴依从性低于 40%",
                    "社区养老服务人手不足，人均服务 200 户以上"])
doc.add_heading("解决方案", level=2)
H.add_numbered(doc, ["毫米波雷达无感监测，无需佩戴", "边缘 AI 模型 3 分钟内识别跌倒/长时间静止",
                     "家属 App + 社区大屏双端联动"])
H.add_quote_block(doc, "“我们不是在卖设备，而是在卖一份让子女安心的确定性。” —— 团队创始人")

doc.add_heading("市场分析", level=1)
doc.add_heading("市场规模", level=2)
H.add_body_paragraph(doc, "如图 1 所示，智慧养老市场规模自 2021 年的 160 亿元增长至 2026 年的 790 亿元，"
                     "年复合增长率达 37.6%。")
H.add_figure(doc, CHART1, "智慧养老市场规模（2021—2026E，亿元）", cap, width_cm=14)
doc.add_heading("竞争格局", level=2)
H.add_body_paragraph(doc, "主要竞争对手对比如表 1 所示（三线表）。")
H.add_caption(doc, "表", "主要竞争对手对比", cap)
H.add_three_line_table(doc, ["公司", "核心产品", "2025 年营收（万元）", "市场份额"],
                       [["A 公司", "穿戴手环", "12,500", "28.5%"], ["B 公司", "摄像头监护", "9,200", "21.0%"],
                        ["本项目", "毫米波雷达", "1,350", "6.5%"]], widths_cm=[3, 4, 4.5, 3], header_fill="D9E2F3")
doc.add_paragraph()
H.add_body_paragraph(doc, "带底纹与合并单元格的全框线表格示例见表 2。")
H.add_caption(doc, "表", "三年财务摘要（万元）", cap)
t = H.add_grid_table(doc, ["指标", "2024A", "2025E", "2026E"],
                     [["营业收入", "480", "980", "1,800"], ["毛利率", "33%", "38%", "42%"], ["净利润", "-120", "60", "260"],
                      ["合并行示例", "", "", ""]], widths_cm=[4, 3.5, 3.5, 3.5])
merged = t.cell(4, 1).merge(t.cell(4, 3))
H._format_cell_text(merged, "此单元格由 3 列合并而成", size=10.5)
H.set_repeat_table_header(t.rows[0])

doc.add_heading("财务预测", level=1)
doc.add_heading("现金流", level=2)
H.add_body_paragraph(doc, "现金流预测见图 2，项目预计在 2027 年实现累计现金流转正。")
H.add_figure(doc, CHART2, "现金流预测（万元）", cap, width_cm=14)
doc.add_heading("段落格式示例", level=2)
H.add_body_paragraph(doc, "本段使用固定行距 20 磅、段前 6 磅、段后 6 磅、首行缩进 2 字符。" * 3,
                     line_spacing=Pt(20), space_before=6, space_after=6)
H.add_body_paragraph(doc, "本段不缩进且两端对齐，1.5 倍行距。" * 4, first_line_indent=False,
                     align=WD_ALIGN_PARAGRAPH.JUSTIFY)
p = H.add_body_paragraph(doc, "更多信息见团队官网：")
H.add_hyperlink(p, "https://example.com/zhishou", url="https://example.com/zhishou")
p.add_run("；返回 ")
H.set_run_font(p.runs[-1])
H.add_hyperlink(p, "执行摘要", anchor="exec_summary")
p.add_run("（文档内部链接）。")
H.set_run_font(p.runs[-1])
# python-docx 1.2.0 新增：批注（评审意见）
doc.add_comment(runs=p.runs[0], text="请补充官网上线时间。", author="评审老师", initials="PS")

# ---------- 第 3 节：横向页 ----------
sec_land = H.add_landscape_section(doc)
doc.add_heading("附录：横向页宽表", level=1)
H.add_caption(doc, "表", "五年利润表预测（横向页）", cap)
H.add_three_line_table(doc, ["项目", "2024A", "2025E", "2026E", "2027E", "2028E", "CAGR"],
                       [["营业收入", "480", "980", "1,800", "3,000", "4,500", "75%"],
                        ["营业成本", "320", "610", "1,040", "1,650", "2,400", "66%"],
                        ["净利润", "-120", "60", "260", "520", "900", "—"]],
                       widths_cm=[4, 3, 3, 3, 3, 3, 3])

# ---------- 第 4 节：恢复纵向 ----------
H.add_portrait_section(doc)
doc.add_heading("参考文献", level=1)
H.add_numbered(doc, ["艾瑞咨询. 2025 年中国智慧养老行业研究报告[R]. 2025.",
                     "国家统计局. 2024 年国民经济和社会发展统计公报[R]. 2025."])

path = os.path.join(OUT, "01_商业计划书示例.docx")
doc.save(path)
print("saved", path)
```

`lo_convert.py`（UNO：目录更新 / xlsx 重算 / 导出 PDF，系统 python3 运行）：

```python
#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""用 LibreOffice UNO 在服务端完成：更新目录/域 → 另存 docx → 导出 PDF；xlsx 强制重算后另存。

必须用系统 python3（apt install python3-uno），不是 venv 的 python。
用法：
  python3 lo_convert.py in.docx  out_dir   # 生成 out_dir/in.updated.docx + out_dir/in.pdf
  python3 lo_convert.py in.xlsx  out_dir   # 重算公式后生成 out_dir/in.recalc.xlsx + out_dir/in.pdf
  python3 lo_convert.py in.pptx  out_dir   # 生成 out_dir/in.pdf
"""
import os
import subprocess
import sys
import time

import uno  # noqa: E402  (python3-uno)
from com.sun.star.beans import PropertyValue  # noqa: E402

PROFILE = "file:///tmp/officetest/lo_profile"  # 独立用户配置，避免与桌面版冲突、可并发多实例
PIPE = "officetest_pipe"


def _prop(name, value):
    p = PropertyValue()
    p.Name, p.Value = name, value
    return p


def start_office():
    cmd = ["soffice", "--headless", "--invisible", "--nologo", "--norestore",
           f"-env:UserInstallation={PROFILE}", f"--accept=pipe,name={PIPE};urp;StarOffice.ComponentContext"]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    local_ctx = uno.getComponentContext()
    resolver = local_ctx.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local_ctx)
    for _ in range(60):
        try:
            ctx = resolver.resolve(f"uno:pipe,name={PIPE};urp;StarOffice.ComponentContext")
            return proc, ctx
        except Exception:
            time.sleep(0.5)
    proc.kill()
    raise RuntimeError("LibreOffice 启动超时")


def convert(path, out_dir):
    path = os.path.abspath(path)
    out_dir = os.path.abspath(out_dir)  # systemPathToFileUrl 需要绝对路径，否则 IOException
    os.makedirs(out_dir, exist_ok=True)
    stem, ext = os.path.splitext(os.path.basename(path))
    ext = ext.lower()
    proc, ctx = start_office()
    try:
        desktop = ctx.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)
        doc = desktop.loadComponentFromURL(uno.systemPathToFileUrl(path), "_blank", 0,
                                           (_prop("Hidden", True), _prop("ReadOnly", False)))
        pdf_filter = {"docx": "writer_pdf_Export", "xlsx": "calc_pdf_Export", "pptx": "impress_pdf_Export"}[ext[1:]]
        if ext == ".docx":
            idx = doc.getDocumentIndexes()
            for i in range(idx.getCount()):
                idx.getByIndex(i).update()          # 目录/图表目录
            doc.getTextFields().refresh()           # PAGE/NUMPAGES/SEQ 等域
            # 目录更新后页码可能变化，再更新一次
            for i in range(idx.getCount()):
                idx.getByIndex(i).update()
            out_docx = os.path.join(out_dir, f"{stem}.updated.docx")
            doc.storeToURL(uno.systemPathToFileUrl(out_docx), (_prop("FilterName", "MS Word 2007 XML"),))
            print("docx (TOC updated):", out_docx, "| indexes:", idx.getCount())
        elif ext == ".xlsx":
            doc.calculateAll()
            out_xlsx = os.path.join(out_dir, f"{stem}.recalc.xlsx")
            doc.storeToURL(uno.systemPathToFileUrl(out_xlsx), (_prop("FilterName", "Calc MS Excel 2007 XML"),))
            print("xlsx (recalculated):", out_xlsx)
        out_pdf = os.path.join(out_dir, f"{stem}.pdf")
        doc.storeToURL(uno.systemPathToFileUrl(out_pdf), (_prop("FilterName", pdf_filter),))
        print("pdf:", out_pdf)
        doc.close(True)
    finally:
        try:
            desktop.terminate()
        except Exception:
            pass
        proc.wait(timeout=20)


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(sys.argv[1])))
```

---

## 3. python-pptx：路演 PPT 生成

封装在 `pptx_helpers.py`（670 行，附 3.9），`02_pptx_demo.py` 生成 15 页 16:9 样例并演示模板占位符替换；LibreOffice 转 PDF→PNG 全部页面检查通过（`04_assets/pptx_全部页面缩略.png`、`pptx_组合图.png`）。

### 3.1 16:9 版式、空白母版、自定义主题色/字体——已验证

- `prs.slide_width, prs.slide_height = Inches(13.333), Inches(7.5)`（默认模板是 4:3 的 10×7.5 in）。**改尺寸不会缩放已有版式里的占位符**（模板复用测试中默认版式的占位符仍按 4:3 坐标摆放），自建页面一律用 `slide_layouts[6]`（Blank）+ 绝对坐标。
- python-pptx 没有主题 API。`set_theme_fonts_and_colors()` 取 `prs.slide_master.part.part_related_by(RT.THEME)`（普通 `Part`，只有 `.blob`），lxml 改 `a:clrScheme`（dk2/lt2/accent1–6/hlink）与 `a:fontScheme`（`majorFont/minorFont` 的 `a:latin`、`a:ea`，并追加 `<a:font script="Hans" typeface="微软雅黑"/>`），再写回 `theme_part._blob`。保存后 theme1.xml 实测：`majorFont: Arial / 微软雅黑`，`accent1=1F3864`。好处：所有未显式设字体的文本（含图表默认文本）都落到主题字体。
- 主题色：商务蓝金 `1F3864 / 2E75B6 / C9A227 / 8FAADC / 7F7F7F / BDD7EE`。

### 3.2 中文字体设置——已验证

`run.font.name` 只写 `a:latin`；中文字形由 `a:ea` 决定，必须手工补 `<a:ea typeface="微软雅黑"/>`，且 `a:rPr` 子元素有固定顺序（`… a:latin, a:ea, a:cs, a:sym, a:hlinkClick …`），`_set_ea_font()` 插在 `a:latin` 之后。图表同理：`chart.font.name` 只产生 `a:latin`，`set_chart_font()` 遍历 `c:txPr` 下所有 `a:defRPr` 补 `a:ea`（实测 chart1.xml 含 `<a:ea typeface="微软雅黑"/>`）。

```python
def _set_ea_font(rpr, typeface):
    for old in rpr.findall(qn("a:ea")): rpr.remove(old)
    ea = etree.Element(qn("a:ea")); ea.set("typeface", typeface)
    latin = rpr.find(qn("a:latin"))
    if latin is not None: latin.addnext(ea); return ea
    for tag in ("a:cs", "a:sym", "a:hlinkClick", "a:hlinkMouseOver", "a:rtl", "a:extLst"):
        nxt = rpr.find(qn(tag))
        if nxt is not None: nxt.addprevious(ea); return ea
    rpr.append(ea); return ea

def style_run(run, size=18, bold=False, color="262626", latin="Arial", ea="微软雅黑", italic=None):
    f = run.font; f.size = Pt(size); f.bold = bold; f.color.rgb = RGBColor.from_string(color); f.name = latin
    _set_ea_font(run._r.get_or_add_rPr(), ea)
```

### 3.3 常用页面类型生成函数——已验证（15 页样例）

| 函数 | 内容 | 要点 |
| --- | --- | --- |
| `cover_slide(prs, title, subtitle, team, date, logo)` | 深色整页底 + 金色分隔线 + 44 磅标题 | 全部矩形/文本框绝对定位 |
| `agenda_slide(prs, items)` | 圆形数字徽标 + 条目，两列自动排布 | `add_number_badge` |
| `section_slide(prs, n, title, tagline)` | 左侧 1/3 深色色块 + 大号章节号 | |
| `bullets_slide(prs, title, items)` | 两级项目符号（`a:buChar` ■ / –，`marL/indent` 控制悬挂） | python-pptx 文本框默认**无**项目符号，必须写 `a:buChar` |
| `image_text_slide(prs, title, img, bullets, crop=...)` | 左图右文/右图左文，图片限高，支持按比例裁剪 | `pic.crop_left/right/top/bottom = 0.03` |
| `chart_slide(prs, title, builder, takeaway=...)` | 原生图表 + 右侧「核心结论」卡片 | builder = `add_bar_chart/add_line_chart/add_pie_chart/add_combo_chart` |
| `table_slide(prs, title, header, rows, col_widths_in, highlight_last_row)` | 深色表头、斑马纹、末行高亮、首列左对齐数字右对齐 | 关闭内置表格样式条纹 `tblPr.bandRow=0` 自定义底色才生效 |
| `team_slide(prs, title, members)` | 卡片 + 头像（无照片时用首字徽标） | |
| `financial_slide(prs, title, kpis, builder)` | 4 个 KPI 卡片 + 下方图表 | |
| `closing_slide(prs)` | 结尾页 | |
| `set_notes(slide, text)` | 演讲者备注 | `slide.notes_slide.notes_text_frame.text`（实测可读回） |
| `add_process_flow(slide, ...)` | V 形箭头流程（PENTAGON + CHEVRON） | 形状代替图标 |

页框架 `_frame()` 统一：顶部色条、标题、金色短线、副标题、页脚文字、页码。

### 3.4 原生图表——已验证

- `CategoryChartData` + `slide.shapes.add_chart(XL_CHART_TYPE.*, ...)`。
- 系列颜色：`series.format.fill.solid(); fore_color.rgb = RGBColor(...)`（柱/饼）；折线用 `series.format.line.color.rgb`、`line.width`、`marker.format.fill`、`smooth=False`。
- 饼/环：逐点上色 `series.points[i].format.fill`；标签 `show_percentage/show_category_name`、`number_format='0.0%'`、`number_format_is_linked=False`。
- 数据标签：`plot.has_data_labels=True; dl.number_format; dl.position = XL_LABEL_POSITION.OUTSIDE_END/ABOVE; dl.font.size/color`。堆叠柱的标签落在深色柱内，用白字。
- 图例 `chart.legend.position = XL_LEGEND_POSITION.BOTTOM; include_in_layout=False`；网格线颜色 `value_axis.major_gridlines.format.line.color.rgb`；坐标轴数字格式 `value_axis.tick_labels.number_format`。
- **组合图（柱 + 折线次坐标轴）python-pptx 无 API**，`add_combo_chart()` 的做法（已验证 LibreOffice 渲染正确，PowerPoint 未验证）：先建含全部系列的柱图 → 在 `c:plotArea` 新建 `c:lineChart`，把折线系列 `c:ser` 移入（删除柱系列专有的 `c:invertIfNegative`）→ 加 `c:marker` 与两个新 `c:axId` → **`lineChart` 必须排在第一个 `c:catAx` 之前**（plotArea 里图表元素先于坐标轴元素）→ 追加隐藏的次分类轴 `c:catAx(delete=1)` 与右侧次数值轴 `c:valAx(axPos=r, crosses=max, numFmt 0%)`。
- **坑：有两个 `c:valAx` 时 `chart.value_axis` 返回的是第 2 个**（python-pptx 为 XY 图设计），主轴要 `ValueAxis(plot_area.findall(qn("c:valAx"))[0])` 显式取；并且新建轴元素要用 `pptx.oxml.parse_xml`（得到 python-pptx 自定义元素类）而不是 `lxml.etree.fromstring`（否则 `.get_or_add_numFmt` 等方法缺失，实测报 `AttributeError`）。

### 3.5 形状与图标、图片裁剪、备注——已验证

- 形状：`add_shape(MSO_SHAPE.ROUNDED_RECTANGLE/OVAL/CHEVRON/PENTAGON)`，`fill.solid()`、`line.fill.background()` 去边框、`shadow.inherit=False` 去阴影（LibreOffice 对圆角矩形仍画淡阴影，PowerPoint 未验证）。
- 图标：python-pptx **不能直接插 SVG**（需 PNG 回退 + `asvg:svgBlip` 扩展，未实现）。可行路线：① 形状 + 数字/符号徽标（本样例）；② 用 matplotlib/PIL 预渲染 PNG 图标；③ 图标字体（Windows 的 Segoe MDL2 Assets）不跨平台，不推荐。
- 图片裁剪：`pic.crop_left = 0.03`（比例），实测写入 `a:srcRect` 并可读回。
- 备注：`slide.notes_slide.notes_text_frame.text = "..."`。

### 3.6 模板 .pptx 复用与占位符替换——已验证

- 读取：`Presentation("模板.pptx")` → 遍历 `slide.placeholders` 得 `(idx, type, name)`（实测 Title Slide：`0 CENTER_TITLE / 1 SUBTITLE`；Title and Content：`0 TITLE / 1 OBJECT`）。
- `fill_placeholders(slide, {idx: text})`：只改第一个 run 的文本、删多余 run → 保留模板字体/颜色；多行用 `\n` 拆成多个段落。
- `replace_text_everywhere(prs, {"{{项目名称}}": "智守康联"})`：按段落合并 run 后替换，覆盖文本框与表格；**跨 run 的占位符**（PowerPoint 编辑时常把 `{{` 与 `}}` 拆开）靠段落级合并解决。
- 删除页 `delete_slide(prs, i)`（`_sldIdLst.remove` + `drop_rel`）；复制页 `duplicate_slide(prs, i)`（深拷贝 shape XML + 复制图片等关系并重映射 `r:embed`；图表部件会被两页共享，复制图表页不建议）。实测复制含图片的第 5 页再删除，文件正常。
- 推荐模式：公司模板只提供**母版/版式/主题**，页面由代码从 Blank 版式生成；需要「原样复用的页」（如固定的致谢页）才用占位符替换。

### 3.7 字体嵌入限制

- python-pptx **无法嵌入字体**（需要 `presentation.xml` 的 `p:embeddedFontLst` + 混淆格式的 `.fntdata` 部件，由 PowerPoint 保存时生成）。LibreOffice 导出 PPTX 可选「嵌入字体」但对 CJK 大字体文件体积巨大，且 PowerPoint 读取 LibreOffice 嵌入字体兼容性不稳定（未验证）。
- 实务建议：PPT 文本用目标机器必有的字体（Windows：微软雅黑/等线；macOS：PingFang SC），主题字体设为微软雅黑；对外发送用 PDF（本流程 LibreOffice 导出 PDF 已嵌入 Noto CJK）；若必须保证 PPT 视觉一致，把关键视觉元素（标题图、数据图）渲染成 PNG 插入。

### 3.8 把 PPT 渲染成预览图——已验证

| 方法 | 命令/代码 | 结果 |
| --- | --- | --- |
| LibreOffice → PDF → pdftoppm | `soffice --headless --convert-to pdf x.pptx` → `pdftoppm -png -r 100 x.pdf out/slide` | 15 页 PDF（959.98×540 pt = 16:9），`slide-01.png…slide-15.png`，2.5 s |
| LibreOffice → PDF → PyMuPDF | `pymupdf.open(pdf)[i].get_pixmap(dpi=100).save(png)` | 同上，无需 poppler，纯 pip |
| `soffice --convert-to png x.pptx` | | **只输出第 1 页**（实测只得到 1 个 PNG），不可用于全页预览 |
| UNO 导出 | `lo_convert.py x.pptx out/` | 产出 PDF，再转 PNG |

### 3.9 代码

`pptx_helpers.py`（完整）：

```python
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
```

`02_pptx_demo.py`：

```python
# -*- coding: utf-8 -*-
"""python-pptx 1.0.2 演示：生成 out/02_路演示例.pptx 与模板复用 out/02_模板填充.pptx"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from pptx import Presentation  # noqa: E402
from pptx.util import Inches  # noqa: E402

import pptx_helpers as P  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
CHART_IMG = os.path.join(OUT, "charts", "09_竞争矩阵.png")
CANVAS_IMG = os.path.join(OUT, "charts", "11_商业模式画布.png")

prs = P.new_presentation()
years = ["2024A", "2025E", "2026E", "2027E", "2028E"]

P.cover_slide(prs, "智守康联：居家养老 AI 健康监测平台", "商业计划路演 | 2026 全国高校商业精英挑战赛", "智守康联团队", "2026 年 10 月")
P.agenda_slide(prs, ["项目概述", "市场分析", "产品与技术", "商业模式", "财务预测", "团队介绍", "融资计划"], page_no=2)
P.section_slide(prs, 1, "项目概述", "为 3,000 万独居老人提供无感健康守护")
P.bullets_slide(prs, "行业痛点与解决方案", [
    "独居老人意外跌倒后平均 4.5 小时才被发现", ("现有穿戴设备依从性低于 40%，夜间摘除率高", 1),
    "毫米波雷达无感监测 + 边缘 AI 3 分钟预警", ("无需佩戴、无摄像头隐私顾虑", 1), ("家属 App + 社区大屏双端联动", 1),
], subtitle="要点页：两级项目符号", page_no=4, notes="演讲备注：强调 4.5 小时这个数字来自民政部 2024 年抽样调查。")
s = P.image_text_slide(prs, "竞争格局", CHART_IMG, ["本项目定位“低价高能”象限", "A 公司价格高、依赖穿戴设备",
                                                 "C 公司仅提供摄像头方案，隐私争议大"], subtitle="图文页：matplotlib 图片 + 要点，图片四边各裁 3%",
                       page_no=5, crop=(0.03, 0.03, 0.03, 0.03))
s, chart = P.chart_slide(prs, "市场规模：五年 CAGR 37.6%", P.add_bar_chart, subtitle="数据图表页：原生堆叠柱状图（可在 PPT 内编辑数据）",
                         page_no=6, takeaway="智慧养老硬件 + 订阅市场 2026 年将达 790 亿元，软件订阅占比从 19% 提升至 38%。",
                         categories=["2021", "2022", "2023", "2024", "2025", "2026E"],
                         series={"硬件": [120, 150, 190, 240, 300, 370], "软件订阅": [30, 55, 90, 140, 210, 300],
                                 "服务": [10, 18, 30, 50, 80, 120]}, stacked=True)
P.chart_slide(prs, "用户增长曲线", P.add_line_chart, subtitle="原生折线图（标记 + 数据标签）", page_no=7,
              categories=years, series={"付费家庭数（户）": [2500, 8000, 30000, 65000, 110000],
                                        "活跃设备数（台）": [3000, 9500, 36000, 78000, 130000]})
P.chart_slide(prs, "收入结构（2026E）", P.add_pie_chart, subtitle="原生环形图：逐点上色 + 百分比标签", page_no=8,
              categories=["设备销售", "年度订阅", "数据增值", "政府采购"], values=[0.42, 0.35, 0.08, 0.15],
              takeaway="订阅收入占比 35%，经常性收入提升估值倍数。")
P.chart_slide(prs, "营收与毛利率", P.add_combo_chart, subtitle="组合图：柱（营收，万元）+ 折线（毛利率，次坐标轴）", page_no=9,
              categories=years, bar_series={"营业收入": [480, 980, 1800, 3000, 4500]},
              line_series={"毛利率": [0.33, 0.38, 0.42, 0.45, 0.47]})
P.table_slide(prs, "五年财务预测摘要", ["指标（万元）", "2024A", "2025E", "2026E", "2027E", "2028E"],
              [["营业收入", "480", "980", "1,800", "3,000", "4,500"], ["营业成本", "320", "610", "1,040", "1,650", "2,400"],
               ["毛利", "160", "370", "760", "1,350", "2,100"], ["费用合计", "280", "310", "500", "830", "1,200"],
               ["净利润", "-120", "60", "260", "520", "900"]], col_widths_in=[2.8, 1.8, 1.8, 1.8, 1.8, 1.8],
              subtitle="表格页：深色表头 + 斑马纹 + 末行高亮", page_no=10, highlight_last_row=True)
P.team_slide(prs, "核心团队", [
    {"name": "李华", "role": "CEO / 创始人", "bio": "清华大学生物医学工程硕士，前某头部医疗器械公司产品总监，8 年智能硬件经验。"},
    {"name": "王敏", "role": "CTO", "bio": "毫米波雷达算法专家，发表 SCI 论文 6 篇，主导 3 代产品研发。"},
    {"name": "张强", "role": "COO", "bio": "原连锁养老机构区域总经理，熟悉民政采购与社区渠道。"},
    {"name": "赵云", "role": "CFO", "bio": "注册会计师，曾任职四大会计师事务所，负责两轮融资。"},
], page_no=11)
P.financial_slide(prs, "财务亮点", [("2026E 营业收入", "1,800 万", "同比 +84%"), ("2026E 净利润", "260 万", "净利率 14%"),
                                 ("累计现金流转正", "2027 年", "无需追加融资"), ("本轮融资", "1,500 万", "出让 15% 股权")],
                  P.add_bar_chart, page_no=12, categories=years,
                  series={"营业收入": [480, 980, 1800, 3000, 4500], "净利润": [-120, 60, 260, 520, 900]})
s = P.image_text_slide(prs, "商业模式画布", CANVAS_IMG, ["九宫格由 matplotlib 生成后插入", "也可用形状 + 文本框原生绘制"],
                       subtitle="图文页（右图）", page_no=13, image_left=False)
s = P.bullets_slide(prs, "实施路径", ["形状与“图标”示例：无图标库时用 形状 + 数字徽标 / V 形箭头流程"], page_no=14)
P.add_process_flow(s, Inches(0.8), Inches(3.0), Inches(11.7), ["需求调研", "产品研发", "试点运营", "市场推广", "融资扩张"])
for i, t in enumerate(["3 个社区试点", "6 项专利", "民政采购入围", "3 万户覆盖"]):
    P.add_number_badge(s, Inches(1.0) + Inches(3.0) * i, Inches(4.6), i + 1)
    P.add_text(s, Inches(1.65) + Inches(3.0) * i, Inches(4.62), Inches(2.3), Inches(0.5), t, size=16)
P.closing_slide(prs)

# 复制一页（文本/形状/图片）再删掉它，验证 duplicate/delete
dup = P.duplicate_slide(prs, 4)
P.delete_slide(prs, len(prs.slides._sldIdLst) - 1)

path = os.path.join(OUT, "02_路演示例.pptx")
prs.save(path)
print("saved", path, "slides:", len(prs.slides))

# ---------------- 模板 .pptx 复用与占位符替换 ----------------
tpl = Presentation()  # 用默认模板里的版式当“公司模板”
tpl.slide_width, tpl.slide_height = P.SLIDE_W, P.SLIDE_H
s1 = tpl.slides.add_slide(tpl.slide_layouts[0])  # Title Slide: idx0 标题, idx1 副标题
s1.shapes.title.text = "{{项目名称}}"
s1.placeholders[1].text = "{{副标题}}"
s2 = tpl.slides.add_slide(tpl.slide_layouts[1])  # Title and Content: idx0 标题, idx1 正文
s2.shapes.title.text = "要点"
s2.placeholders[1].text = "待填充"
tpl_path = os.path.join(OUT, "02_模板.pptx")
tpl.save(tpl_path)

prs2 = Presentation(tpl_path)
for s in prs2.slides:
    print("layout:", s.slide_layout.name, [(ph.placeholder_format.idx, ph.placeholder_format.type, ph.name) for ph in s.placeholders])
P.replace_text_everywhere(prs2, {"{{项目名称}}": "智守康联", "{{副标题}}": "居家养老 AI 健康监测平台"})
P.fill_placeholders(prs2.slides[1], {0: "三大核心优势", 1: "无感监测\n3 分钟预警\n双端联动"})
out2 = os.path.join(OUT, "02_模板填充.pptx")
prs2.save(out2)
print("saved", out2)
```

---

## 4. openpyxl：财务模型 Excel 生成

`03_xlsx_demo.py` 生成 9 张工作表的 DCF 财务模型（公式版），`03b_xlsx_values.py` 演示公式值问题的三种对策并产出数值版。LibreOffice 重算后转 PDF 检查（`04_assets/xlsx_全部工作表.png`、`xlsx_利润表与图表.png`）。

### 4.1 工作簿结构——已验证

| 工作表 | 内容 | 关键公式 |
| --- | --- | --- |
| 假设 | 11 行逐年假设（蓝字、关键项黄底）+ 全局参数（税率/折现率/永续增长率/初始投资/营运资金比，均注册为**命名区域**）+ 情景下拉（数据验证）+ 情景系数 | `=IF(C23="乐观",1.2,IF(C23="悲观",0.8,1))` |
| 收入 | 付费家庭数×情景系数、设备收入、订阅收入、合计、增速、订阅占比 | `='假设'!C5*情景系数`、`=IF(B8=0,"",C8/B8-1)` |
| 成本 | 硬件/订阅成本、三项费用、折旧（SUMPRODUCT 分摊历年资本开支） | `=SUMPRODUCT('假设'!$C$14:E$14/'假设'!$C$15:E$15)` |
| 利润表 | 收入→毛利→EBITDA→EBIT→所得税（亏损不缴）→净利润→净利率 | `=-IF(C12>0,C12*税率,0)` |
| 现金流量表 | 净利润+折旧−资本开支−营运资金增加=FCF，期初投资在 B 列，累计现金流 | `=-初始投资`、`=SUM($B$9:C9)` |
| 资产负债表 | 简化：货币资金/营运资金/固定资产净值/初始投资长期资产 vs 实收资本/留存收益，**检验行** | `=ROUND(C9-C12,2)` → 实测 5 年均为 0 |
| 估值 | NPV、Gordon 终值、EV、项目 NPV、IRR、回收期、情景加权 NPV | `=NPV(C4,'现金流量表'!C9:G9)`、`=IRR('现金流量表'!B9:G9)`、`=COUNTIF(...,"<0")+1`、`=SUMPRODUCT(B16:B18,C16:C18)` |
| 敏感性 | 折现率一维敏感性（公式实时算） | `=NPV(A5,'现金流量表'!$C$9:$G$9)-初始投资` |
| 图表 | 引用利润表/现金流量表的数据区 + 原生折线图 + 柱状图 | |

验证数据（LibreOffice 重算 ⇄ numpy-financial）：FCF = [-800, -404.8, -213.1, 272.2, 1529.0, 3554.4]；NPV@12% = **1851.02**（npf.npv 1851.02）；IRR = **38.84%**（npf.irr 0.3884）；回收期 4 年；资产负债表检验行 [0,0,0,0,0]。

### 4.2 公式写入规则（LibreOffice 24.2 实测 + Excel 规范）

| 写法 | LibreOffice 重算结果 | Excel（规范/经验，未验证） | 结论 |
| --- | --- | --- | --- |
| `=SUM('数据 表'!A1:A2)`（表名加单引号） | 30 ✔ | ✔ | **所有表名一律加单引号**（含空格/中文/标点时必须，不含时无害） |
| `=SUM(数据 表!A1:A2)`（有空格未加引号） | #VALUE! | 文件修复/错误 | 禁止 |
| `=sum(...)` 小写函数名 | 30 ✔（回写为大写） | ✔（Excel 大小写不敏感） | 可用，但建议统一大写 |
| `=IF(x>5;1;0)` 分号分隔 | #VALUE!（被改写成 `=IF(...;1;0,TRUE())`） | 文件损坏提示 | **只能用英文逗号**，openpyxl 原样写入不做转换 |
| `=IF(x>5,1,0)` | 1 ✔ | ✔ | |
| `=XLOOKUP(...)` 无前缀 | #NAME? | #NAME?（Excel 要求 `_xlfn.` 前缀） | 2010 以后新增函数必须写 `_xlfn.XLOOKUP`、`_xlfn.IFS`、`_xlfn.TEXTJOIN`… |
| `=_xlfn.XLOOKUP(...)` | #NAME?（**LibreOffice 24.2 没有 XLOOKUP**，24.8 起支持） | ✔（Excel 365） | 服务端若用 LibreOffice 重算，避免 XLOOKUP，用 INDEX/MATCH |
| `=_xlfn.IFS(...)` | 大 ✔ | ✔ | |
| `=TEXT(A1,"0.00")&"万元"` | `10.00万元` ✔ | ✔ | |
| `=数据表!A1` 不存在的表 | #NAME? | #REF!/#NAME? | 生成时先校验表名 |
| `openpyxl.utils.FORMULAE` | `NPV/IRR/SUMPRODUCT/IF/COUNTIF` ∈ FORMULAE；`XLOOKUP/IFS/TEXTJOIN` ∉ | | 不在集合里的函数名就是要加 `_xlfn.` 的函数（实测输出见附录 A） |

其它：相对/绝对引用与 Excel 一致（`$C$9`）；中文命名区域 `税率`、`情景系数` 在 LibreOffice 与 `formulas` 库中均可直接参与公式（已验证）；`wb.defined_names[name] = DefinedName(name, attr_text="'假设'!$C$18")` 是 3.1 的写法（旧版 `wb.defined_names.append()` 已不可用）。

### 4.3 数字格式与颜色约定——已验证

```python
NUM  = '#,##0;[Red](#,##0);"-"'      # 千分位，负数红色加括号，零显示 -
NUM1 = '#,##0.0;[Red](#,##0.0);"-"'
PCT  = '0.0%'
MULT = '0.0"x"'
WAN  = '#,##0"万元"'
F_INPUT   = Font(name="Arial", size=10, color="0000FF")   # 蓝：手工输入
F_FORMULA = Font(name="Arial", size=10, color="000000")   # 黑：本表公式
F_LINK    = Font(name="Arial", size=10, color="008000")   # 绿：跨表引用
FILL_INPUT = PatternFill("solid", fgColor="FFF2CC")       # 黄底：核心假设
FILL_TOTAL = PatternFill("solid", fgColor="D9E2F3")       # 合计行
TOP_BOTTOM = Border(top=Side(style="thin"), bottom=Side(style="double"))   # 净利润行上单下双线
```

`put(ws, row, col, value, kind="input|formula|link", fmt=...)` 一个函数统一字体/格式/对齐/边框；中文标签用微软雅黑，数字用 Arial（LibreOffice 下回退到 Noto Sans CJK / Liberation Sans，PDF 实测正常）。

### 4.4 冻结窗格、命名区域、数据验证、原生图表、打印设置——已验证

- 冻结：`ws.freeze_panes = "C5"`（冻结前 4 行与前 2 列）。
- 数据验证：下拉 `DataValidation(type="list", formula1='"基准,乐观,悲观"')`；范围 `DataValidation(type="decimal", operator="between", formula1="0", formula2="1")` + `dv.add("C8:G10")`；LibreOffice 往返后 2 条验证仍在。
- 批注：`cell.comment = Comment("2024 年实际数…", "财务模型")`。
- 原生图表：`LineChart/BarChart` + `Reference` + `add_data(titles_from_data=True)` + `set_categories`；系列颜色 `series.graphicalProperties.line.solidFill = "1F3864"`、线宽 EMU、`marker.symbol="circle"`；**`series.smooth = False` 必须显式写**——不写 `<smooth>` 时 LibreOffice 按平滑曲线渲染（实测）；`y_axis.delete = False / x_axis.delete = False` 防止 openpyxl 3.1 默认隐藏坐标轴；图表数据源引用公式单元格时显示依赖缓存值（见 4.5）。
- 打印/PDF：`page_setup.orientation="landscape"`、`sheet_properties.pageSetUpPr.fitToPage=True`、`fitToWidth=1, fitToHeight=0`——LibreOffice 转 PDF 从 16 页碎片变成每表 1 页共 9 页。
- 隐藏网格线 `ws.sheet_view.showGridLines = False`；列宽 `column_dimensions["A"].width = 24`。

### 4.5 公式计算值为空：问题与对策——全部已验证

问题：openpyxl 只写公式字符串，不写 `<v>` 缓存值。`load_workbook(..., data_only=True)` 读到 `None`（实测）；pandas `read_excel` 同样为 NaN；原生图表若引用这些单元格，在未重算的查看器（如某些预览组件）里是空图。Excel/WPS 打开会自动计算，但「另存」之前文件里仍无值。

| 对策 | 做法 | 实测 | 保留样式/图表/命名区域/数据验证 | 适用 |
| --- | --- | --- | --- | --- |
| A. LibreOffice 重算 | `soffice --headless --convert-to xlsx --outdir out 模型.xlsx`（LibreOffice 对没有缓存值的单元格一律计算）或 UNO `doc.calculateAll()` 后 `storeToURL(..., "Calc MS Excel 2007 XML")` | 0.9 s；NPV/IRR 与 numpy-financial 一致；检验行全 0 | 是（实测命名区域 6 个、数据验证 2 条、图表 2 个、冻结窗格均在），但文件被 LibreOffice 重写，样式细节可能微调 | 服务端读值、出 PDF |
| B. `formulas` 库 | `xl = formulas.ExcelModel().loads(path).finish(); sol = xl.calculate(); sol["'[模型.xlsx]估值'!C10"].value` | 0.8 s；支持 NPV/IRR/IF/SUMPRODUCT/COUNTIF/跨表/中文命名区域，结果与 A 完全一致 | 只读计算，不改文件（也可 `xl.write(dirpath=...)` 写出带值的副本，未测） | 无 LibreOffice 的环境、单元测试 |
| C. Python 自算 + 数值版 | 以公式版为模板 `load_workbook(公式版)`，把每个 `str.startswith("=")` 的单元格替换为 A/B 得到的值后另存「数值版」 | 240 个单元格，<0.1 s；pandas 直接读到利润表数据框 | 全保留（同一文件改值） | 交付评委的「只读数值版」、供 PPT/Word 取数 |
| 补充 | `wb.calculation = CalcProperties(fullCalcOnLoad=True)` | 写入 `workbook.xml` `calcPr fullCalcOnLoad=1` | | 让 Excel 打开时强制全量重算（Excel 行为未验证） |

推荐流水线：**公式版（交付、可编辑）→ LibreOffice 重算得到值 → 数值版（评委只读 + 系统取数）**，PPT/Word 中的财务数字一律从数值版读取，保证三份材料口径一致。

### 4.6 代码

`03_xlsx_demo.py`（公式版，完整）：

```python
# -*- coding: utf-8 -*-
"""openpyxl 3.1.5 财务模型演示：生成 out/03_财务模型.xlsx（公式版）与 out/03_财务模型_数值版.xlsx。

颜色约定（投行/咨询建模惯例）：蓝色字=手工输入假设，黑色字=本表公式，绿色字=跨表引用，黄色底=关键假设。
"""
import os

import numpy_financial as npf
from openpyxl import Workbook
from openpyxl.chart import BarChart, LineChart, Reference
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import FORMULAE, get_column_letter
from openpyxl.workbook.defined_name import DefinedName
from openpyxl.workbook.properties import CalcProperties
from openpyxl.worksheet.datavalidation import DataValidation

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

YEARS = [2024, 2025, 2026, 2027, 2028]
N = len(YEARS)
COLS = [get_column_letter(3 + i) for i in range(N)]  # C..G 为年份列

# ---------------------------------------------------------------- 样式 ----
F_INPUT = Font(name="Arial", size=10, color="0000FF")          # 蓝：输入
F_FORMULA = Font(name="Arial", size=10, color="000000")        # 黑：公式
F_LINK = Font(name="Arial", size=10, color="008000")           # 绿：跨表引用
F_HEADER = Font(name="微软雅黑", size=10, bold=True, color="FFFFFF")
F_TITLE = Font(name="微软雅黑", size=14, bold=True, color="1F3864")
F_LABEL = Font(name="微软雅黑", size=10)
F_BOLD = Font(name="微软雅黑", size=10, bold=True)
FILL_HEADER = PatternFill("solid", fgColor="1F3864")
FILL_INPUT = PatternFill("solid", fgColor="FFF2CC")
FILL_TOTAL = PatternFill("solid", fgColor="D9E2F3")
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(top=THIN, bottom=THIN, left=THIN, right=THIN)
TOP_BOTTOM = Border(top=Side(style="thin", color="000000"), bottom=Side(style="double", color="000000"))
NUM = '#,##0;[Red](#,##0);"-"'            # 负数红色加括号，零显示 -
NUM1 = '#,##0.0;[Red](#,##0.0);"-"'
PCT = "0.0%"
MULT = '0.0"x"'
WAN = '#,##0"万元"'
RIGHT = Alignment(horizontal="right", vertical="center")
CENTER = Alignment(horizontal="center", vertical="center")


def title(ws, text, sub=None):
    ws["A1"] = text
    ws["A1"].font = F_TITLE
    if sub:
        ws["A2"] = sub
        ws["A2"].font = Font(name="微软雅黑", size=9, color="7F7F7F")
    ws.sheet_view.showGridLines = False


def header_row(ws, row, labels, start_col=1):
    for j, lab in enumerate(labels):
        c = ws.cell(row=row, column=start_col + j, value=lab)
        c.font = F_HEADER
        c.fill = FILL_HEADER
        c.alignment = CENTER
        c.border = BORDER


def label(ws, row, text, bold=False, indent=0):
    c = ws.cell(row=row, column=1, value=text)
    c.font = F_BOLD if bold else F_LABEL
    c.alignment = Alignment(indent=indent)
    return c


def put(ws, row, col, value, kind="formula", fmt=NUM, bold=False, fill=None):
    """kind: input / formula / link"""
    c = ws.cell(row=row, column=col, value=value)
    c.font = Font(name="Arial", size=10, bold=bold,
                  color={"input": "0000FF", "formula": "000000", "link": "008000"}[kind])
    c.number_format = fmt
    c.alignment = RIGHT
    c.border = BORDER
    if fill:
        c.fill = fill
    return c


def year_row(ws, row, text, cells, kind="formula", fmt=NUM, bold=False, fill=None, indent=0):
    label(ws, row, text, bold, indent)
    for i, v in enumerate(cells):
        put(ws, row, 3 + i, v, kind, fmt, bold, fill)


def set_widths(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


wb = Workbook()

# ================================================================ 假设页 ====
ws = wb.active
ws.title = "假设"
title(ws, "关键假设（蓝字 = 可修改输入；黄底 = 核心假设）", "单位：万元 / 户 / %")
header_row(ws, 4, ["项目", "单位"] + [f"{y}E" if y > 2024 else f"{y}A" for y in YEARS])
assump = [
    # (标签, 单位, 值列表, 格式, 名称)
    ("付费家庭数", "户", [2500, 8000, 30000, 65000, 110000], "#,##0", None),
    ("客单价（设备）", "元/户", [1200, 1150, 1100, 1050, 1000], "#,##0", None),
    ("年订阅费", "元/户/年", [360, 360, 400, 400, 420], "#,##0", None),
    ("订阅转化率", "%", [0.6, 0.65, 0.7, 0.72, 0.75], PCT, None),
    ("硬件成本率", "%", [0.62, 0.58, 0.55, 0.52, 0.50], PCT, None),
    ("订阅服务成本率", "%", [0.35, 0.32, 0.30, 0.28, 0.26], PCT, None),
    ("销售费用率", "%", [0.25, 0.20, 0.16, 0.14, 0.12], PCT, None),
    ("研发费用", "万元", [150, 180, 220, 260, 300], "#,##0", None),
    ("管理费用", "万元", [80, 100, 130, 160, 200], "#,##0", None),
    ("资本开支", "万元", [200, 150, 250, 300, 300], "#,##0", None),
    ("折旧年限", "年", [5, 5, 5, 5, 5], "0", None),
]
ROW_A = {}
r = 5
for lab, unit, vals, fmt, _ in assump:
    label(ws, r, lab)
    ws.cell(row=r, column=2, value=unit).font = F_LABEL
    for i, v in enumerate(vals):
        put(ws, r, 3 + i, v, "input", fmt, fill=FILL_INPUT if lab in ("付费家庭数", "客单价（设备）") else None)
    ROW_A[lab] = r
    r += 1
r += 1
label(ws, r, "全局参数", bold=True)
r += 1
globals_ = [("所得税率", 0.25, PCT, "税率"), ("折现率（WACC）", 0.12, PCT, "折现率"), ("永续增长率", 0.03, PCT, "永续增长率"),
            ("初始投资（2023 末）", 800, "#,##0", "初始投资"), ("营运资金占收入比", 0.10, PCT, "营运资金比")]
for lab, v, fmt, name in globals_:
    label(ws, r, lab)
    put(ws, r, 3, v, "input", fmt, fill=FILL_INPUT)
    # 命名区域：公式里可写 税率 而不是 假设!$C$20
    wb.defined_names[name] = DefinedName(name, attr_text=f"'假设'!$C${r}")
    ROW_A[lab] = r
    r += 1
label(ws, r, "情景选择")
put(ws, r, 3, "基准", "input", "@", fill=FILL_INPUT)
dv = DataValidation(type="list", formula1='"基准,乐观,悲观"', allow_blank=False, showDropDown=False)
dv.error, dv.errorTitle = "请选择 基准/乐观/悲观", "无效输入"
ws.add_data_validation(dv)
dv.add(ws.cell(row=r, column=3))
ROW_A["情景选择"] = r
r += 1
label(ws, r, "情景收入系数")
put(ws, r, 3, f'=IF(C{r - 1}="乐观",1.2,IF(C{r - 1}="悲观",0.8,1))', "formula", "0.00")
wb.defined_names["情景系数"] = DefinedName("情景系数", attr_text=f"'假设'!$C${r}")
ROW_A["情景系数"] = r
dv2 = DataValidation(type="decimal", operator="between", formula1="0", formula2="1", allow_blank=False)
dv2.error, dv2.errorTitle = "比率应在 0~1 之间", "范围错误"
ws.add_data_validation(dv2)
dv2.add(f"C{ROW_A['订阅转化率']}:G{ROW_A['订阅服务成本率']}")
ws.cell(row=ROW_A["付费家庭数"], column=3).comment = Comment("2024 年实际数，来源：公司后台。", "财务模型")
ws.freeze_panes = "C5"
set_widths(ws, [24, 12, 12, 12, 12, 12, 12])

# ================================================================ 收入页 ====
wr = wb.create_sheet("收入")
title(wr, "收入预测", "单位：万元；绿色 = 引用假设页")
header_row(wr, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wr, 5, "付费家庭数（户）", [f"='假设'!{c}{ROW_A['付费家庭数']}*情景系数" for c in COLS], "link", "#,##0")
year_row(wr, 6, "设备收入", [f"={c}5*'假设'!{c}{ROW_A['客单价（设备）']}/10000" for c in COLS], "link")
year_row(wr, 7, "订阅收入", [f"={c}5*'假设'!{c}{ROW_A['订阅转化率']}*'假设'!{c}{ROW_A['年订阅费']}/10000" for c in COLS], "link")
year_row(wr, 8, "营业收入合计", [f"=SUM({c}6:{c}7)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wr, 9, "收入增速", ["" if i == 0 else f"=IF({COLS[i - 1]}8=0,\"\",{c}8/{COLS[i - 1]}8-1)" for i, c in enumerate(COLS)],
         "formula", PCT)
year_row(wr, 10, "订阅收入占比", [f"=IF({c}8=0,0,{c}7/{c}8)" for c in COLS], "formula", PCT)
wr.freeze_panes = "C5"
set_widths(wr, [24, 4, 12, 12, 12, 12, 12])

# ================================================================ 成本页 ====
wc = wb.create_sheet("成本")
title(wc, "成本与费用", "单位：万元")
header_row(wc, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wc, 5, "硬件成本", [f"='收入'!{c}6*'假设'!{c}{ROW_A['硬件成本率']}" for c in COLS], "link")
year_row(wc, 6, "订阅服务成本", [f"='收入'!{c}7*'假设'!{c}{ROW_A['订阅服务成本率']}" for c in COLS], "link")
year_row(wc, 7, "营业成本合计", [f"=SUM({c}5:{c}6)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wc, 8, "销售费用", [f"='收入'!{c}8*'假设'!{c}{ROW_A['销售费用率']}" for c in COLS], "link")
year_row(wc, 9, "研发费用", [f"='假设'!{c}{ROW_A['研发费用']}" for c in COLS], "link")
year_row(wc, 10, "管理费用", [f"='假设'!{c}{ROW_A['管理费用']}" for c in COLS], "link")
year_row(wc, 11, "期间费用合计", [f"=SUM({c}8:{c}10)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
# 折旧：SUMPRODUCT 对历年资本开支按折旧年限分摊（简化：直线法，当年起计）
dep = []
for i, c in enumerate(COLS):
    rng = f"'假设'!$C${ROW_A['资本开支']}:{c}${ROW_A['资本开支']}"
    life = f"'假设'!$C${ROW_A['折旧年限']}:{c}${ROW_A['折旧年限']}"
    dep.append(f"=SUMPRODUCT({rng}/{life})")
year_row(wc, 12, "折旧摊销（SUMPRODUCT）", dep, "link")
wc.freeze_panes = "C5"
set_widths(wc, [24, 4, 12, 12, 12, 12, 12])

# ================================================================ 利润表 ====
wi = wb.create_sheet("利润表")
title(wi, "利润表", "单位：万元")
header_row(wi, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wi, 5, "营业收入", [f"='收入'!{c}8" for c in COLS], "link", bold=True)
year_row(wi, 6, "减：营业成本", [f"=-'成本'!{c}7" for c in COLS], "link")
year_row(wi, 7, "毛利", [f"={c}5+{c}6" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wi, 8, "毛利率", [f"=IF({c}5=0,0,{c}7/{c}5)" for c in COLS], "formula", PCT)
year_row(wi, 9, "减：期间费用", [f"=-'成本'!{c}11" for c in COLS], "link")
year_row(wi, 10, "EBITDA", [f"={c}7+{c}9" for c in COLS], "formula", bold=True)
year_row(wi, 11, "减：折旧摊销", [f"=-'成本'!{c}12" for c in COLS], "link")
year_row(wi, 12, "EBIT", [f"={c}10+{c}11" for c in COLS], "formula", bold=True)
year_row(wi, 13, "减：所得税（IF 亏损不缴）", [f"=-IF({c}12>0,{c}12*税率,0)" for c in COLS], "formula")
year_row(wi, 14, "净利润", [f"={c}12+{c}13" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wi, 15, "净利率", [f"=IF({c}5=0,0,{c}14/{c}5)" for c in COLS], "formula", PCT)
for c in COLS:
    wi[f"{c}14"].border = TOP_BOTTOM
wi.freeze_panes = "C5"
set_widths(wi, [26, 4, 12, 12, 12, 12, 12])

# ================================================================ 现金流量表 ====
wf = wb.create_sheet("现金流量表")
title(wf, "现金流量表（自由现金流口径）", "单位：万元")
header_row(wf, 4, ["项目", "2023"] + [str(y) for y in YEARS])
year_row(wf, 5, "净利润", [f"='利润表'!{c}14" for c in COLS], "link")
year_row(wf, 6, "加：折旧摊销", [f"='成本'!{c}12" for c in COLS], "link")
year_row(wf, 7, "减：资本开支", [f"=-'假设'!{c}{ROW_A['资本开支']}" for c in COLS], "link")
year_row(wf, 8, "减：营运资金增加", [f"=-('收入'!{c}8-{'0' if i == 0 else f'收入!{COLS[i - 1]}8'})*营运资金比"
                                 if i else f"=-'收入'!{c}8*营运资金比" for i, c in enumerate(COLS)], "link")
year_row(wf, 9, "自由现金流（FCF）", [f"=SUM({c}5:{c}8)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
put(wf, 9, 2, "=-初始投资", "link", NUM, bold=True, fill=FILL_TOTAL)  # B9 = 期初投资
year_row(wf, 10, "累计现金流", [f"=SUM($B$9:{c}9)" for c in COLS], "formula", bold=True)
put(wf, 10, 2, "=B9", "formula", NUM, bold=True)
wf.freeze_panes = "C5"
set_widths(wf, [26, 12, 12, 12, 12, 12, 12])

# ================================================================ 资产负债表 ====
wbsh = wb.create_sheet("资产负债表")
title(wbsh, "资产负债表（简化）", "单位：万元；检验行应为 0")
header_row(wbsh, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wbsh, 5, "货币资金", [f"='现金流量表'!{c}10+初始投资" for c in COLS], "link")
year_row(wbsh, 6, "营运资金", [f"='收入'!{c}8*营运资金比" for c in COLS], "link")
year_row(wbsh, 7, "固定资产净值", [f"=SUM('假设'!$C${ROW_A['资本开支']}:{c}${ROW_A['资本开支']})-SUM('成本'!$C$12:{c}$12)" for c in COLS], "link")
year_row(wbsh, 8, "初始投资形成的长期资产", ["=初始投资"] * N, "link")
year_row(wbsh, 9, "资产合计", [f"=SUM({c}5:{c}8)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wbsh, 10, "实收资本", ["=初始投资"] * N, "link")
year_row(wbsh, 11, "留存收益", [f"=SUM('利润表'!$C$14:{c}$14)" for c in COLS], "link")
year_row(wbsh, 12, "所有者权益合计", [f"=SUM({c}10:{c}11)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wbsh, 13, "检验：资产-权益（应为 0）", [f"=ROUND({c}9-{c}12,2)" for c in COLS], "formula", NUM)
wbsh.freeze_panes = "C5"
set_widths(wbsh, [26, 4, 12, 12, 12, 12, 12])

# ================================================================ 估值页 ====
wv = wb.create_sheet("估值")
title(wv, "DCF 估值", "单位：万元")
label(wv, 4, "折现率（WACC）")
put(wv, 4, 3, "=折现率", "link", PCT)
label(wv, 5, "永续增长率")
put(wv, 5, 3, "=永续增长率", "link", PCT)
label(wv, 6, "预测期 FCF 现值（NPV）")
put(wv, 6, 3, f"=NPV(C4,'现金流量表'!C9:G9)", "link")
label(wv, 7, "终值（Gordon）")
put(wv, 7, 3, "='现金流量表'!G9*(1+C5)/(C4-C5)", "link")
label(wv, 8, "终值现值")
put(wv, 8, 3, f"=C7/(1+C4)^{N}", "formula")
label(wv, 9, "企业价值（EV）", bold=True)
put(wv, 9, 3, "=C6+C8", "formula", bold=True, fill=FILL_TOTAL)
label(wv, 10, "项目 NPV（含初始投资）", bold=True)
put(wv, 10, 3, "=C6-初始投资", "formula", bold=True, fill=FILL_TOTAL)
label(wv, 11, "项目 IRR", bold=True)
put(wv, 11, 3, "=IRR('现金流量表'!B9:G9)", "link", PCT, bold=True, fill=FILL_TOTAL)
label(wv, 12, "投资回收期（年，静态）")
put(wv, 12, 3, "=COUNTIF('现金流量表'!C10:G10,\"<0\")+1", "link", "0")
label(wv, 14, "情景加权 NPV（SUMPRODUCT）", bold=True)
header_row(wv, 15, ["情景", "概率", "NPV"])
for i, (name, p, mult) in enumerate([("悲观", 0.25, 0.7), ("基准", 0.5, 1.0), ("乐观", 0.25, 1.3)]):
    wv.cell(row=16 + i, column=1, value=name).font = F_LABEL
    put(wv, 16 + i, 2, p, "input", PCT)
    put(wv, 16 + i, 3, f"=$C$10*{mult}", "formula")
label(wv, 19, "加权 NPV", bold=True)
put(wv, 19, 3, "=SUMPRODUCT(B16:B18,C16:C18)", "formula", bold=True, fill=FILL_TOTAL)
set_widths(wv, [28, 12, 14])

# ================================================================ 敏感性 ====
wsn = wb.create_sheet("敏感性")
title(wsn, "敏感性分析", "上表：公式实时计算（折现率一维）；下表：Python 预计算数值（二维）")
header_row(wsn, 4, ["折现率", "NPV（含初始投资）", "相对基准"])
for i, rate in enumerate([0.08, 0.10, 0.12, 0.14, 0.16]):
    put(wsn, 5 + i, 1, rate, "input", PCT)
    put(wsn, 5 + i, 2, f"=NPV(A{5 + i},'现金流量表'!$C$9:$G$9)-初始投资", "link")
    put(wsn, 5 + i, 3, f"=IF('估值'!$C$10=0,0,B{5 + i}/'估值'!$C$10-1)", "link", PCT)
set_widths(wsn, [16, 20, 14, 14, 14, 14])

# ================================================================ 图表页 ====
wch = wb.create_sheet("图表")
title(wch, "图表（原生 Excel 图表，数据引用利润表/现金流量表）")
header_row(wch, 3, ["年份", "营业收入", "净利润", "自由现金流"])
for i, (y, c) in enumerate(zip(YEARS, COLS)):
    wch.cell(row=4 + i, column=1, value=str(y)).font = F_LABEL
    put(wch, 4 + i, 2, f"='利润表'!{c}5", "link")
    put(wch, 4 + i, 3, f"='利润表'!{c}14", "link")
    put(wch, 4 + i, 4, f"='现金流量表'!{c}9", "link")
line = LineChart()
line.title = "营业收入与净利润（万元）"
line.y_axis.title = "万元"
line.x_axis.title = "年份"
line.height, line.width = 8, 16
line.add_data(Reference(wch, min_col=2, min_row=3, max_col=3, max_row=3 + N), titles_from_data=True)
line.set_categories(Reference(wch, min_col=1, min_row=4, max_row=3 + N))
for s, color in zip(line.series, ["1F3864", "C9A227"]):
    s.graphicalProperties.line.solidFill = color
    s.graphicalProperties.line.width = 28000  # EMU
    s.marker.symbol = "circle"
    s.marker.graphicalProperties.solidFill = color
    s.smooth = False  # 不写 <smooth> 时 LibreOffice 会按平滑曲线渲染，显式关闭
line.y_axis.majorGridlines.spPr = None
line.y_axis.delete = False  # openpyxl 3.1 默认轴可能被隐藏，显式打开
line.x_axis.delete = False
wch.add_chart(line, "F3")
bar = BarChart()
bar.type = "col"
bar.title = "自由现金流（万元）"
bar.height, bar.width = 8, 16
bar.add_data(Reference(wch, min_col=4, min_row=3, max_row=3 + N), titles_from_data=True)
bar.set_categories(Reference(wch, min_col=1, min_row=4, max_row=3 + N))
bar.series[0].graphicalProperties.solidFill = "2E75B6"
bar.series[0].invertIfNegative = False
bar.dataLabels = None
bar.y_axis.delete = False
bar.x_axis.delete = False
bar.legend = None
wch.add_chart(bar, "F21")
set_widths(wch, [10, 14, 14, 14])

# 打印/PDF 导出设置：横向、缩放到一页宽（LibreOffice 转 PDF 时也生效）
for sh in wb.worksheets:
    sh.page_setup.orientation = "landscape"
    sh.page_setup.paperSize = sh.PAPERSIZE_A4
    sh.sheet_properties.pageSetUpPr.fitToPage = True
    sh.page_setup.fitToWidth = 1
    sh.page_setup.fitToHeight = 0
    sh.print_options.horizontalCentered = True

# 让 Excel 打开时强制全量重算（python 写的公式没有缓存值）
wb.calculation = CalcProperties(fullCalcOnLoad=True)
wb.properties.creator = "bizgold-office-engine"

path = os.path.join(OUT, "03_财务模型.xlsx")
wb.save(path)
print("saved", path)
print("FORMULAE check:", {f: (f in FORMULAE) for f in ["NPV", "IRR", "SUMPRODUCT", "IF", "XLOOKUP", "IFS", "TEXTJOIN", "COUNTIF"]})
```

`03b_xlsx_values.py`（三种对策 + 数值版）：

```python
# -*- coding: utf-8 -*-
"""公式缓存值为空的三种对策，全部实际运行：
A) LibreOffice 重算（命令行 / UNO lo_convert.py）
B) formulas 库纯 Python 计算
C) Python 自算后写『数值版』：复制公式版工作簿，把每个公式单元格替换为计算值（样式不变）
"""
import os
import shutil
import subprocess
import sys
import time

from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SRC = os.path.join(OUT, "03_财务模型.xlsx")
ENV = dict(os.environ, FONTCONFIG_FILE=os.path.join(HERE, "fonts.conf"))

# ---------- A) LibreOffice 重算 ----------
t0 = time.time()
recalc_dir = os.path.join(OUT, "xlsx_recalc")
subprocess.run(["soffice", "--headless", "--convert-to", "xlsx", "--outdir", recalc_dir, SRC], check=True,
               capture_output=True, env=ENV)
recalc_path = os.path.join(recalc_dir, os.path.basename(SRC))
wb_vals = load_workbook(recalc_path, data_only=True)
print(f"A) LibreOffice 重算 {time.time() - t0:.1f}s  NPV={wb_vals['估值']['C10'].value:.2f}  IRR={wb_vals['估值']['C11'].value:.4f}"
      f"  检验行={[wb_vals['资产负债表'][f'{c}13'].value for c in 'CDEFG']}")

# ---------- B) formulas 库 ----------
try:
    import formulas
    t0 = time.time()
    xl = formulas.ExcelModel().loads(SRC).finish()
    sol = xl.calculate()
    def get(sheet, cell):
        key = f"'[{os.path.basename(SRC)}]{sheet.upper()}'!{cell}"
        return sol[key].value[0][0] if key in sol else None
    print(f"B) formulas {formulas.__version__} 计算 {time.time() - t0:.1f}s  NPV={get('估值','C10')}  IRR={get('估值','C11')}"
          f"  净利润2026={get('利润表','E14')}  情景选择={get('假设','C23')} 情景系数={get('假设','C24')}  加权NPV={get('估值','C19')}")
except Exception as e:  # noqa: BLE001
    print("B) formulas 库失败:", type(e).__name__, str(e)[:300])

# ---------- C) 数值版 ----------
# 做法：以公式版为模板，用 LibreOffice(或 formulas/自算) 得到的值覆盖公式单元格 → 样式/图表/命名区域全部保留，
# 适合给评委看的“只读数值版”，也适合 pandas 读取。
t0 = time.time()
wb_f = load_workbook(SRC)              # 公式版（含样式）
n = 0
for ws in wb_f.worksheets:
    ws_v = wb_vals[ws.title]
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str) and c.value.startswith("="):
                c.value = ws_v[c.coordinate].value
                n += 1
wb_f.calculation.fullCalcOnLoad = False
dst = os.path.join(OUT, "03_财务模型_数值版.xlsx")
wb_f.save(dst)
chk = load_workbook(dst, data_only=True)
print(f"C) 数值版写入 {n} 个单元格 {time.time() - t0:.1f}s  NPV={chk['估值']['C10'].value:.2f}  (pandas 可直接读)")
import pandas as pd
df = pd.read_excel(dst, sheet_name="利润表", header=3, index_col=0).dropna(how="all")
print(df.iloc[:, 1:].head(6).round(1).to_string())
```

---

## 5. matplotlib：商业图表

`bizcharts.py`（449 行）+ `04_mpl_demo.py` 生成 11 张 300dpi PNG，逐张目视检查中文无乱码、布局无遮挡（样图：`04_assets/01_市场规模堆叠柱.png`、`04_瀑布图.png`、`06_龙卷风.png`、`08_现金流组合.png`、`11_商业模式画布.png`）。

### 5.1 中文字体自动选择——已验证

```python
FONT_CANDIDATES = {
    "Windows": ["Microsoft YaHei", "SimHei", "DengXian", "SimSun", "KaiTi"],
    "Darwin":  ["PingFang SC", "Hiragino Sans GB", "Heiti SC", "STHeiti", "Songti SC"],
    "Linux":   ["Noto Sans CJK SC", "Source Han Sans SC", "Source Han Sans CN", "WenQuanYi Micro Hei",
                "WenQuanYi Zen Hei", "AR PL UMing CN", "Droid Sans Fallback"],
}
def setup_chinese_font(prefer=None, verbose=True) -> str:
    names = {f.name for f in fm.fontManager.ttflist}
    chosen = next((c for c in candidates if c in names), None)
    if chosen is None:
        fm._load_fontmanager(try_read_cache=False)     # 重建缓存（刚 apt 装完字体时必需）
        ...
    plt.rcParams["font.sans-serif"] = [chosen, "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False          # 否则负号显示为方框
```

实测要点：
1. **字体缓存**：本机装 fonts-noto-cjk 之前 `~/.cache/matplotlib/fontlist-v3.11.0.json` 已存在，首次 `fontManager` 只看到 `Droid Sans Fallback / WenQuanYi Micro Hei`；`fm._load_fontmanager(try_read_cache=False)` 重建后出现 `Noto Sans CJK SC/TC/JP/KR/HK`、`Noto Serif CJK *`、`Noto Sans Mono CJK *`（并自动写回缓存，之后新进程直接可见）。替代做法：删除 `fontlist-*.json`。
2. **.ttc 多字体面**：Ubuntu 的 Noto CJK 是 `NotoSansCJK-Regular.ttc`（JP/KR/SC/TC/HK 10 个面）。matplotlib 3.11.2 能枚举全部面（`findfont("Noto Sans CJK SC")` → 该 ttc；bold 自动命中 `NotoSansCJK-Bold.ttc`）。若旧版 matplotlib 只看到 `Noto Sans CJK JP`，可用 fontTools 抽出 SC 面：`TTCollection(path).fonts[idx].save("NotoSansCJKsc-Regular.otf")` 再 `fm.fontManager.addfont(...)`（已验证可行）。
3. 用 `fc-list :lang=zh family` 检查系统字体；无字体时函数抛出带 apt 命令与 fc-list 清单的 `RuntimeError`。
4. Windows 上 matplotlib 识别的是**英文族名**（`Microsoft YaHei`、`SimHei`、`SimSun`、`KaiTi`、`DengXian`），不是「微软雅黑」（未在 Windows 验证，依据 font_manager 读取 name 表英文名的实现）；Word/PPT 的 XML 里中英文名都可用。
5. 输出 PNG 用 `savefig(dpi=300, bbox_inches="tight", facecolor="white")`；若要输出 PDF/SVG 嵌入 Word 不适用（python-docx 只接受位图），PNG 即可；CFF 的 OTF（Noto CJK）在 matplotlib PDF 后端需 `pdf.fonttype=3`（未验证，仅 PNG 路线已验证）。

### 5.2 统一配色——已验证

```python
PALETTES = {
    "blue_gold":  {"primary": "#1F3864", "secondary": "#2E75B6", "accent": "#C9A227", "accent2": "#E2B93B",
                   "positive": "#2E8B57", "negative": "#C0392B", "neutral": "#7F7F7F", "light": "#D6DCE5",
                   "series": ["#1F3864", "#2E75B6", "#C9A227", "#8FAADC", "#7F7F7F", "#BDD7EE", "#E2B93B", "#404040"]},
    "consulting": {"primary": "#051C2C", "secondary": "#2251FF", "accent": "#00A9F4", "accent2": "#AAE6F0",
                   "positive": "#00A651", "negative": "#E5133B", "neutral": "#8C8C8C", "light": "#E6E6E6",
                   "series": ["#051C2C", "#2251FF", "#00A9F4", "#AAE6F0", "#8C8C8C", "#C4C4C4", "#034B6F", "#1E6B52"]},
}
apply_style()   # 去上/右边框、浅灰横向网格、标题左对齐加粗、图例无框、savefig.dpi=300、prop_cycle=series
```

与 Word（题注黑体、表头底纹 D9E2F3、深色表头 1F3864）和 PPT 主题（accent1=1F3864、强调金 C9A227）同一套色值，三件套视觉统一。

### 5.3 函数签名（全部已实现并出图）

| 函数 | 签名 | 用途 |
| --- | --- | --- |
| 市场规模堆叠柱 | `market_size_stacked_bar(years, segments: dict[str, list], title, unit="亿元", cagr_label=True, source=None)` | 自动算 CAGR 并标注 |
| 份额条形 | `share_barh(labels, values, title, highlight="本项目", unit="%", source=None)` | 本项目高亮金色 |
| 增长折线 | `growth_line(x, series: dict, title, unit, annotate_last=True, yoy_bars=None, source=None)` | 可叠加同比柱（次轴） |
| 瀑布图 | `waterfall(labels, deltas, title, total_label="净利润", unit="万元")` | 首项起始值，自动合计柱 |
| 漏斗 | `funnel(stages, values, title)` | 最窄保留 22% 宽度，右侧标转化率 |
| 敏感性龙卷风 | `tornado(factors, low, high, base, title, unit, low_label="-10%", high_label="+10%")` | 按摆幅排序，标签贴柱外端 |
| 盈亏平衡 | `breakeven(fixed_cost, price, var_cost, max_units, title, unit, qty_unit)` | 标注平衡点、盈亏区填色 |
| 现金流柱线组合 | `cashflow_combo(periods, inflow, outflow, title, unit)` | 流入/流出柱 + 累计净现金流线 |
| 竞争矩阵散点 | `competition_matrix(names, x, y, size=None, highlight=None, xlabel, ylabel, quadrant_labels)` | 气泡大小=规模，四象限标签 |
| 甘特图 | `gantt(tasks=[(名称, 开始期, 持续, 阶段)], total_periods=12, period_label="月")` | 按阶段配色 |
| 商业模式画布 | `business_model_canvas(blocks: dict[九宫格名, list[str]], title)` | 10×6 网格九宫格 |
| 通用 | `setup_chinese_font()`, `use_palette(name)`, `apply_style(base_font_size)`, `save_fig(fig, path, dpi=300)` | |

出图过程中修过的问题（均已重新出图确认）：龙卷风标签与柱重叠（改为贴各自柱外端）、现金流累计值标签压线（改 `annotate` 偏移 + 白底框）、增长折线末端标注与右轴刻度重叠（横轴改序号并留白）、漏斗窄柱内白字不可见（最小宽度 + 文字外置）。

### 5.4 代码

`bizcharts.py`（完整）：

```python
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
```

`04_mpl_demo.py`：

```python
# -*- coding: utf-8 -*-
"""生成 11 张商业样图到 out/charts/，并打印所用字体。"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import bizcharts as bc  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "out", "charts")
os.makedirs(OUT, exist_ok=True)

font = bc.setup_chinese_font()
bc.apply_style()

years = [2021, 2022, 2023, 2024, 2025, 2026]
bc.save_fig(bc.market_size_stacked_bar(
    years, {"硬件": [120, 150, 190, 240, 300, 370], "软件订阅": [30, 55, 90, 140, 210, 300], "服务": [10, 18, 30, 50, 80, 120]},
    title="智慧养老市场规模（亿元）", source="艾瑞咨询、团队测算"), f"{OUT}/01_市场规模堆叠柱.png")

bc.save_fig(bc.share_barh(["A公司", "B公司", "本项目", "C公司", "其他"], [28.5, 21.0, 6.5, 15.0, 29.0],
                          title="细分市场份额（2025）", highlight="本项目", source="公司年报"), f"{OUT}/02_份额条形.png")

bc.save_fig(bc.growth_line(years, {"营业收入": [320, 480, 690, 980, 1350, 1800], "毛利": [96, 158, 248, 372, 540, 756]},
                           title="收入与毛利增长（万元）", unit="万元", yoy_bars=[0, 0.5, 0.44, 0.42, 0.38, 0.33],
                           source="财务模型"), f"{OUT}/03_增长折线.png")

bc.save_fig(bc.waterfall(["营业收入", "原材料", "人工", "营销费用", "研发费用", "管理费用", "税金"],
                         [1800, -620, -310, -220, -180, -120, -90], title="2026E 利润瀑布（万元）", source="财务模型"),
            f"{OUT}/04_瀑布图.png")

bc.save_fig(bc.funnel(["曝光", "访问", "注册", "试用", "付费"], [1_200_000, 360_000, 90_000, 27_000, 8_100],
                      title="用户转化漏斗", source="后台数据"), f"{OUT}/05_漏斗.png")

bc.save_fig(bc.tornado(["客单价", "获客成本", "续费率", "人工成本", "折现率"],
                       low=[1800, 3100, 2200, 2900, 3000], high=[3600, 2300, 3200, 2500, 2400], base=2700,
                       source="敏感性分析"), f"{OUT}/06_龙卷风.png")

bc.save_fig(bc.breakeven(fixed_cost=600, price=0.45, var_cost=0.18, max_units=4000, title="盈亏平衡分析",
                         unit="万元", qty_unit="台"), f"{OUT}/07_盈亏平衡.png")

bc.save_fig(bc.cashflow_combo(["2024", "2025", "2026", "2027", "2028"], [300, 700, 1400, 2300, 3300],
                              [900, 950, 1200, 1700, 2200], title="现金流预测（万元）"), f"{OUT}/08_现金流组合.png")

bc.save_fig(bc.competition_matrix(["A公司", "B公司", "C公司", "本项目", "D公司"], [8, 6, 3, 4, 7], [7, 5, 3, 8, 4],
                                  size=[90, 60, 30, 20, 50], highlight="本项目"), f"{OUT}/09_竞争矩阵.png")

bc.save_fig(bc.gantt([("需求调研", 1, 2, "筹备期"), ("产品研发", 2, 4, "研发期"), ("试点运营", 5, 3, "试运营"),
                      ("市场推广", 7, 4, "推广期"), ("融资与扩张", 9, 4, "推广期")], total_periods=12),
            f"{OUT}/10_甘特图.png")

bc.save_fig(bc.business_model_canvas({
    "重要伙伴": ["社区养老服务中心", "三甲医院老年科", "硬件代工厂"], "关键业务": ["健康监测设备研发", "AI 预警算法"],
    "核心资源": ["专利 6 项", "医疗数据库", "研发团队"], "价值主张": ["24h 无感健康监测", "异常 3 分钟预警", "家属远程可视"],
    "客户关系": ["订阅制服务", "社区驻点"], "渠道通路": ["民政采购", "线上商城", "养老机构"],
    "客户细分": ["独居老人家庭", "养老机构", "社区政府"], "成本结构": ["硬件 BOM 成本", "研发人员薪酬", "获客成本"],
    "收入来源": ["设备销售", "年度订阅费", "数据增值服务"]}), f"{OUT}/11_商业模式画布.png")

print("charts written:", sorted(os.listdir(OUT)))
print("font:", font)
```

---

## 6. Markdown→docx 中间层方案对比与推荐流水线

### 6.1 方案对比（pandoc 与 docxtpl 均实际运行）

| 方案 | 实测情况 | 优点 | 缺点 | 结论 |
| --- | --- | --- | --- | --- |
| 直接 python-docx 编程 | 第 2 节 | 控制力最强 | 代理需要调用 API 而非写文本；章节内容与排版耦合 | 作为渲染器的底层 |
| pandoc 3.1.3（`pandoc ch.md -o x.docx --toc --resource-path out`） | 转换成功；XML 检查：TOC 域在 `w:sdt` 内但**无页码**、无 `updateFields`；图注为 `ImageCaption` 样式**无 SEQ 域**；样式表无 `w:eastAsia` 字体（只有 `lang=en-US`）；表格样式 `Table`（无三线/底纹）；标题无编号 | 零代码、Markdown 支持最全（脚注、公式、引用）、可用 `--reference-doc` 套用样式 | 中文字体/三线表/题注域/多级编号都要靠 reference.docx + 后处理；外部二进制；细节不可控 | 备选/应急 |
| docxtpl 0.20.2（Jinja2 模板） | `{{ project }}`、`{%tr for m in members %}` 行循环实测成功（注意 `{%tr for%}` 与 `{%tr endfor%}` 必须各占独立表行，写在同一行内报 `TemplateSyntaxError: unknown tag 'endfor'`） | 版式由人在 Word 里定，适合**报名表、封面、固定结构附件** | 不适合长篇自由内容；跨 run 的标签易碎 | 固定表单用 |
| 自写渲染器 `md2docx.py`（markdown-it-py 4.2.0 + docx_helpers） | 272 行；样章渲染：3 级标题自动编号、首行缩进段落、粗/斜体/代码/链接、嵌套列表、三线表+题注、图+题注、引用块、分页、代码块全部正确（`04_assets/md2docx_渲染结果.png`） | 代理只写 Markdown；排版规则集中在渲染器；可随时加约定（题注、分页、宽度） | 需维护；Markdown 子集 | **推荐** |

### 6.2 推荐流水线（适合 AI 代理逐章生成）

```
代理（LLM）逐章输出 Markdown ──┐
   · # / ## / ### 标题（渲染器自动多级编号，不要手写 1.1）
   · 段落 / **粗体** / *斜体* / `代码` / [链接](url)
   · - / 1. 列表（支持嵌套）
   · 表：标题行 + GFM 表格（渲染为三线表 + "表 N" 题注）
   · ![题注](charts/xxx.png "width=13cm")（图 + "图 N" 题注）
   · > 评审提示/引言（引用块）
   · <!-- pagebreak -->（强制分页）
图表（bizcharts.py）→ PNG 300dpi ──┤
财务数值版（openpyxl 数值版）→ 数字 ──┤
                                   ▼
            md2docx.build_document(chapters, out.docx, base_dir, 封面信息, 页眉)
                 ├─ setup_chinese_document + add_heading_numbering
                 ├─ 封面 → 目录(TOC 域, updateFields) → 新节(页眉/页码从 1)
                 └─ MarkdownDocxRenderer.render(每章 md)
                                   ▼
            lo_convert.py out.docx  →  out.pdf（目录已更新）+ out.updated.docx
```

要点：
- **一个 `MarkdownDocxRenderer` 实例贯穿全书**：图/表题注计数器、列表编号、脚注都在实例上，分章 `render()` 多次即可连续编号。
- 代理不写章节号、不写「图 1」，由渲染器生成；正文引用图表用「如下图所示/见上表」或由系统在生成后回填编号。
- 图片路径相对 `base_dir`；markdown-it 会把中文路径百分号编码，渲染器已 `unquote`（首次实测没解码导致「缺失图片」，已修复）。
- 表格题注约定：表格前一段以 `表：` 或 `Table:` 开头。
- 每章渲染完成后做 XML 顺序自检（2.6）并用 LibreOffice 转一次 PDF 做"能打开"检查。

### 6.3 渲染器核心代码

`md2docx.py`（完整，已验证）：

```python
# -*- coding: utf-8 -*-
"""Markdown → docx 渲染器（markdown-it-py 4.2.0 + python-docx 1.2.0 + docx_helpers）。

面向 AI 代理逐章生成：代理只输出 Markdown（含图片引用、表格、题注约定），渲染器负责排版。
支持：# 标题层级（自动多级编号）、段落（首行缩进）、**粗体**/*斜体*/`代码`/[链接]、无序/有序/嵌套列表、
GFM 表格（前一段以“表：”开头即为题注，否则取表格上方 `Table:` 行）、图片 ![题注](路径 "width=12cm")、
> 引用块、分页标记 <!-- pagebreak -->、``` 代码块、--- 忽略。
"""
from __future__ import annotations

import os
import re
from typing import Iterator

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Pt
from markdown_it import MarkdownIt
from markdown_it.token import Token

import docx_helpers as H

PAGEBREAK_RE = re.compile(r"^\s*(<!--\s*pagebreak\s*-->|\\pagebreak|---\s*pagebreak\s*---)\s*$", re.I)
WIDTH_RE = re.compile(r"width\s*=\s*([\d.]+)\s*(cm|in)?", re.I)
TABLE_CAPTION_RE = re.compile(r"^(表|Table)\s*[:：]\s*(.+)$")


class MarkdownDocxRenderer:
    def __init__(self, doc: Document | None = None, base_dir: str = ".", heading_offset: int = 0,
                 table_style: str = "three_line", image_width_cm: float = 14.0):
        self.doc = doc or H.setup_chinese_document(Document())
        if doc is None:
            H.add_heading_numbering(self.doc, levels=3)
        self.base_dir = base_dir
        self.heading_offset = heading_offset
        self.table_style = table_style
        self.image_width_cm = image_width_cm
        self.captions = H.CaptionCounter()
        self.footnotes: H.Footnotes | None = None
        self.md = MarkdownIt("commonmark", {"html": True}).enable(["table", "strikethrough"])
        self._pending_table_caption: str | None = None
        self._list_stack: list[tuple[str, int | None]] = []  # (type, num_id)

    # ------------------------------------------------------------ 入口 ----
    def render(self, md_text: str):
        tokens = self.md.parse(md_text)
        it = iter(range(len(tokens)))
        i = 0
        while i < len(tokens):
            i = self._dispatch(tokens, i)
        return self.doc

    def save(self, path: str):
        self.doc.save(path)
        return path

    # -------------------------------------------------------- 块级分发 ----
    def _dispatch(self, tokens: list[Token], i: int) -> int:
        t = tokens[i]
        typ = t.type
        if typ == "heading_open":
            level = int(t.tag[1]) + self.heading_offset
            inline = tokens[i + 1]
            p = self.doc.add_heading("", level=min(level, 9))
            self._inline(p, inline.children or [], east_asia="黑体")
            return i + 3
        if typ == "paragraph_open":
            inline = tokens[i + 1]
            text = inline.content.strip()
            # 分页标记（单独一段）
            if PAGEBREAK_RE.match(text):
                self.doc.add_page_break()
                return i + 3
            m = TABLE_CAPTION_RE.match(text)
            if m and i + 3 < len(tokens) and tokens[i + 3].type == "table_open":
                self._pending_table_caption = m.group(2)
                return i + 3
            # 纯图片段落 → 图 + 题注
            children = inline.children or []
            if children and all(c.type in ("image", "softbreak") for c in children):
                for c in children:
                    if c.type == "image":
                        self._image(c)
                return i + 3
            in_list = bool(self._list_stack)
            if in_list:
                p = self._list_paragraph()
            else:
                p = self.doc.add_paragraph()
            self._inline(p, children)
            return i + 3
        if typ == "html_block":
            if PAGEBREAK_RE.match(t.content.strip()):
                self.doc.add_page_break()
            return i + 1
        if typ == "bullet_list_open":
            self._list_stack.append(("bullet", None))
            return i + 1
        if typ == "ordered_list_open":
            self._list_stack.append(("ordered", None))
            return i + 1
        if typ in ("bullet_list_close", "ordered_list_close"):
            self._list_stack.pop()
            return i + 1
        if typ in ("list_item_open", "list_item_close"):
            return i + 1
        if typ == "blockquote_open":
            j = i + 1
            texts = []
            while tokens[j].type != "blockquote_close":
                if tokens[j].type == "inline":
                    texts.append(tokens[j].content)
                j += 1
            H.add_quote_block(self.doc, "\n".join(texts))
            return j + 1
        if typ == "table_open":
            return self._table(tokens, i)
        if typ == "fence" or typ == "code_block":
            self._code_block(t.content)
            return i + 1
        if typ == "hr":
            return i + 1
        return i + 1

    # ------------------------------------------------------------ 列表 ----
    def _list_paragraph(self):
        depth = len(self._list_stack)
        kind, num_id = self._list_stack[-1]
        suffix = "" if depth == 1 else f" {min(depth, 3)}"
        style = ("List Bullet" if kind == "bullet" else "List Number") + suffix
        p = self.doc.add_paragraph(style=style)
        if kind == "ordered":
            if num_id is None:  # 每个有序列表从 1 重新编号
                num_id = H.restart_list_numbering(self.doc, p, style)
                self._list_stack[-1] = (kind, num_id)
            else:
                H.continue_list_numbering(p, num_id)
        return p

    # ------------------------------------------------------------ 行内 ----
    def _inline(self, p, children: list[Token], east_asia="宋体"):
        bold = italic = code = strike = False
        link: str | None = None
        for c in children:
            if c.type == "text":
                if link:
                    H.add_hyperlink(p, c.content, url=link, east_asia=east_asia)
                else:
                    r = p.add_run(c.content)
                    H.set_run_font(r, "等线" if code else east_asia, "Consolas" if code else "Times New Roman",
                                   bold=bold or None, italic=italic or None)
                    if strike:
                        r.font.strike = True
            elif c.type == "strong_open":
                bold = True
            elif c.type == "strong_close":
                bold = False
            elif c.type == "em_open":
                italic = True
            elif c.type == "em_close":
                italic = False
            elif c.type == "s_open":
                strike = True
            elif c.type == "s_close":
                strike = False
            elif c.type == "code_inline":
                r = p.add_run(c.content)
                H.set_run_font(r, "等线", "Consolas", size=10.5)
            elif c.type == "link_open":
                link = c.attrs.get("href")
            elif c.type == "link_close":
                link = None
            elif c.type in ("softbreak",):
                p.add_run(" ")
            elif c.type == "hardbreak":
                p.add_run().add_break()
            elif c.type == "image":
                self._image(c)
            elif c.type == "html_inline":
                if "<br" in c.content:
                    p.add_run().add_break()

    # ------------------------------------------------------------ 图片 ----
    def _image(self, tok: Token):
        from urllib.parse import unquote
        src = unquote(tok.attrs.get("src", ""))  # markdown-it 会把中文路径百分号编码，必须解码
        path = src if os.path.isabs(src) else os.path.join(self.base_dir, src)
        alt = "".join(c.content for c in (tok.children or [])) or tok.content
        width = self.image_width_cm
        title = tok.attrs.get("title") or ""
        m = WIDTH_RE.search(title)
        if m:
            width = float(m.group(1)) * (2.54 if (m.group(2) or "cm").lower() == "in" else 1)
        if not os.path.exists(path):
            p = self.doc.add_paragraph()
            H.set_run_font(p.add_run(f"[缺失图片：{src}]"), color="C00000")
            return
        H.add_figure(self.doc, path, alt, self.captions, width_cm=width)

    # ------------------------------------------------------------ 表格 ----
    def _table(self, tokens: list[Token], i: int) -> int:
        rows: list[list[str]] = []
        aligns: list[str | None] = []
        j = i + 1
        cur: list[str] | None = None
        in_head = False
        while tokens[j].type != "table_close":
            tt = tokens[j]
            if tt.type == "thead_open":
                in_head = True
            elif tt.type == "thead_close":
                in_head = False
            elif tt.type == "tr_open":
                cur = []
            elif tt.type == "tr_close":
                rows.append(cur or [])
            elif tt.type in ("th_open", "td_open"):
                if in_head:
                    style = tt.attrs.get("style", "")
                    aligns.append(style.split(":")[-1].strip() if style else None)
            elif tt.type == "inline" and cur is not None:
                cur.append(tt.content)
            j += 1
        if not rows:
            return j + 1
        header, body = rows[0], rows[1:]
        if self._pending_table_caption:
            H.add_caption(self.doc, "表", self._pending_table_caption, self.captions)
            self._pending_table_caption = None
        if self.table_style == "grid":
            H.add_grid_table(self.doc, header, body)
        else:
            H.add_three_line_table(self.doc, header, body, header_fill="D9E2F3")
        self.doc.add_paragraph()  # 表后空行，避免与下文粘连
        return j + 1

    # ---------------------------------------------------------- 代码块 ----
    def _code_block(self, content: str):
        from docx.oxml import parse_xml
        from docx.oxml.ns import nsdecls
        p = self.doc.add_paragraph()
        p.paragraph_format.first_line_indent = Pt(0)
        p._p.get_or_add_pPr().get_or_add_ind().set("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}firstLineChars", "0")
        H.insert_ordered(p._p.get_or_add_pPr(), parse_xml(f'<w:shd {nsdecls("w")} w:val="clear" w:color="auto" w:fill="F5F5F5"/>'))
        lines = content.rstrip("\n").split("\n")
        for k, line in enumerate(lines):
            r = p.add_run(line)
            H.set_run_font(r, "等线", "Consolas", size=9)
            if k < len(lines) - 1:
                r.add_break()


# ---------------------------------------------------------- 组装流水线 ----

def build_document(chapters_md: list[str], out_path: str, base_dir: str, title: str, subtitle: str,
                   meta: list[tuple[str, str]], header_text: str, with_toc: bool = True):
    """代理逐章产出的 Markdown 列表 → 完整 docx（封面 + 目录 + 页眉页脚 + 正文）。"""
    doc = H.setup_chinese_document(Document())
    H.add_heading_numbering(doc, levels=3)
    H.add_cover_page(doc, title, subtitle, meta=meta, new_section_after=False)
    if with_toc:
        doc.add_page_break()
        H.add_toc(doc, levels="1-3")
    sec = doc.add_section()
    H.set_header_text(sec, header_text)
    H.add_page_number_footer(sec, "第 {PAGE} 页 / 共 {NUMPAGES} 页")
    H.restart_page_numbering(sec, start=1)
    renderer = MarkdownDocxRenderer(doc, base_dir=base_dir)
    for md in chapters_md:
        renderer.render(md)
    renderer.save(out_path)
    return out_path
```

`05_md2docx_demo.py`（样章 + pandoc 对比 + docxtpl 对比）：

```python
# -*- coding: utf-8 -*-
"""Markdown→docx 流水线演示 + pandoc / docxtpl 对比。"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from md2docx import build_document  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

CH1 = """# 市场分析

## 市场规模与增长

智慧养老市场正处于**高速成长期**。根据艾瑞咨询数据，2021—2026 年市场规模的年复合增长率达到 *37.6%*，其中软件订阅占比持续提升（参见 `图 1`）。

![智慧养老市场规模（亿元）](charts/01_市场规模堆叠柱.png "width=13cm")

核心驱动因素包括：

- 人口老龄化加速：60 岁以上人口占比 **21.1%**
- 政策支持：民政部“十四五”养老服务规划明确智慧养老方向
  - 地方配套补贴覆盖 28 个省份
  - 社区适老化改造专项资金
- 技术成熟：毫米波雷达成本三年下降 60%

## 竞争格局

表：主要竞争对手对比（2025 年）

| 公司 | 核心产品 | 营收（万元） | 市场份额 |
|:-----|:--------|------------:|--------:|
| A 公司 | 穿戴手环 | 12,500 | 28.5% |
| B 公司 | 摄像头监护 | 9,200 | 21.0% |
| 本项目 | 毫米波雷达 | 1,350 | 6.5% |

> 评审提示：对比表应突出本项目差异化，而非罗列参数。

实施步骤如下：

1. 完成 3 个社区试点（2026 Q1）
2. 取得医疗器械二类备案
3. 启动 A 轮融资，详见 [融资计划](https://example.com/financing)

<!-- pagebreak -->

# 财务预测

## 收入预测

| 年份 | 2024A | 2025E | 2026E |
|------|------:|------:|------:|
| 营业收入 | 480 | 980 | 1,800 |
| 净利润 | -120 | 60 | 260 |

![现金流预测（万元）](charts/08_现金流组合.png)

```
假设：客单价 1,200 元 → 1,000 元（年降 4%）
      订阅转化率 60% → 75%
```
"""

path = build_document([CH1], os.path.join(OUT, "05_md渲染示例.docx"), base_dir=OUT, title="智守康联商业计划书",
                      subtitle="（Markdown 渲染流水线示例）", meta=[("团队", "智守康联"), ("日期", "2026-10")],
                      header_text="智守康联商业计划书")
print("saved", path)

# ---- pandoc 对比：同一份 Markdown ----
md_path = os.path.join(OUT, "05_chapter.md")
open(md_path, "w", encoding="utf-8").write(CH1)
r = subprocess.run(["pandoc", md_path, "-o", os.path.join(OUT, "05_pandoc.docx"), "--toc", "--resource-path", OUT],
                   capture_output=True, text=True)
print("pandoc rc:", r.returncode, r.stderr.strip()[:200])

# ---- docxtpl 对比：Jinja2 模板填充 ----
from docx import Document  # noqa: E402
from docxtpl import DocxTemplate  # noqa: E402

tpl_doc = Document()
tpl_doc.add_heading("{{ project }} 报名信息表", 1)
tpl_doc.add_paragraph("团队：{{ team }}；指导教师：{{ advisor }}")
t = tpl_doc.add_table(rows=4, cols=3)  # {%tr for %} 与 {%tr endfor %} 必须各占独立的一行
t.style = "Table Grid"
t.cell(0, 0).text, t.cell(0, 1).text, t.cell(0, 2).text = "姓名", "角色", "学院"
t.cell(1, 0).text = "{%tr for m in members %}"
t.cell(2, 0).text, t.cell(2, 1).text, t.cell(2, 2).text = "{{ m.name }}", "{{ m.role }}", "{{ m.college }}"
t.cell(3, 0).text = "{%tr endfor %}"
tpl_path = os.path.join(OUT, "05_docxtpl_模板.docx")
tpl_doc.save(tpl_path)
dt = DocxTemplate(tpl_path)
dt.render({"project": "智守康联", "team": "智守康联团队", "advisor": "张三",
           "members": [{"name": "李华", "role": "CEO", "college": "经管学院"}, {"name": "王敏", "role": "CTO", "college": "电子学院"}]})
dt.save(os.path.join(OUT, "05_docxtpl_输出.docx"))
d = Document(os.path.join(OUT, "05_docxtpl_输出.docx"))
print("docxtpl:", d.paragraphs[0].text, "| 表格行数:", len(d.tables[0].rows), [c.text for c in d.tables[0].rows[1].cells])
```

---

## 7. 各库版本、坑点与规避（汇总）

| # | 库/工具（版本） | 坑点 | 规避 | 状态 |
| --- | --- | --- | --- | --- |
| 1 | python-docx 1.2.0 | `font.name` 不设中文字体 | `rFonts` 写 `w:eastAsia`，并删 4 个 `*Theme` 属性（Heading 样式自带） | 已验证 |
| 2 | python-docx | 默认模板是 Letter、`lang=en-US`、无 Hyperlink/脚注样式 | `setup_chinese_document()` 改 A4、zh-CN；脚注样式自建 | 已验证 |
| 3 | python-docx | 无脚注 API | 手工 `footnotes.xml` 部件（`Footnotes` 类） | 已验证 |
| 4 | python-docx | 无编号 API；`TOC Heading` 继承 Heading 1 的编号 | 手写 abstractNum/num；TOC Heading `numId=0` | 已验证 |
| 5 | python-docx | `List Number` 跨列表连续编号 | 新建 `w:num` + `startOverride` | 已验证 |
| 6 | python-docx | `add_section()` 复制上一节 `pgNumType`，页码重启、LibreOffice 多出空白页 | 新节删除 `pgNumType` | 已验证 |
| 7 | python-docx | 手工 append 的 pPr/tcPr/tblPr/trPr 子元素顺序错 | `insert_ordered()` 按 schema 顺序插入 + 自检脚本 | 已验证 |
| 8 | python-docx | 横向页只改 orientation 不换纸张 | 交换 `page_width/page_height` | 已验证 |
| 9 | python-docx | 目录域不更新 | Word：`updateFields`（未验证）；服务端：UNO `DocumentIndexes.update()`（已验证）；`--convert-to pdf` 不更新 | 已验证 |
| 10 | python-docx | `STYLEREF`/`SEQ \s` 章节式题注在 LibreOffice 下错误 | 用顺序编号或 Python 写死「图 2-1」 | 已验证 |
| 11 | python-docx | `doc.paragraphs` 不含 `w:sdt` 内段落（LibreOffice 保存的目录） | 需要时直接查 XML | 已验证 |
| 12 | python-docx | 只设 `column.width` 列宽不生效 | 每个 cell 设 `width` + `tblLayout fixed` | 经验（样例用此法，PDF 列宽正确） |
| 13 | python-pptx 1.0.2 | 文本/图表中文字体：`font.name` 只写 `a:latin` | 补 `a:ea`（run、`c:txPr/a:defRPr`）；主题 fontScheme 一并改 | 已验证 |
| 14 | python-pptx | 无主题 API | 改 theme1.xml blob | 已验证 |
| 15 | python-pptx | 无组合图 API | 手工 `c:lineChart` + 次坐标轴；lineChart 置于 catAx 前 | LibreOffice 已验证，PowerPoint 未验证 |
| 16 | python-pptx | 两个 valAx 时 `chart.value_axis` 返回次轴；`etree.fromstring` 造的元素缺 python-pptx 方法 | 显式 `ValueAxis(valAx_lst[0])`；用 `pptx.oxml.parse_xml` | 已验证 |
| 17 | python-pptx | 文本框默认无项目符号 | 写 `a:buChar/a:buClr/a:buFont` + `marL/indent` | 已验证 |
| 18 | python-pptx | 改 `slide_width/height` 不缩放模板版式 | 自建页用 Blank 版式 + 绝对坐标；模板要本来就是 16:9 | 已验证 |
| 19 | python-pptx | 不能嵌入字体、不能直接插 SVG、无复制页 API | 用通用字体/PDF 分发；PNG 图标；`duplicate_slide()`（图表页不复制） | 已验证（复制含图片页） |
| 20 | python-pptx | 表格自定义底色被内置样式条纹覆盖 | `tblPr.bandRow="0"` | 已验证 |
| 21 | openpyxl 3.1.5 | 公式无缓存值（`data_only` 为 None、pandas NaN、图表空） | LibreOffice 重算 / `formulas` / 数值版 / `fullCalcOnLoad` | 已验证 |
| 22 | openpyxl | 分号分隔符、未加引号的表名 → #VALUE!/损坏 | 逗号 + 表名一律 `'表名'!` | 已验证（LibreOffice） |
| 23 | openpyxl | 新函数需 `_xlfn.` 前缀；LibreOffice 24.2 无 XLOOKUP | 查 `FORMULAE`；服务端重算避免 XLOOKUP | 已验证 |
| 24 | openpyxl | `defined_names.append` 已移除 | `wb.defined_names[name] = DefinedName(...)` | 已验证 |
| 25 | openpyxl | 折线不写 `smooth` → LibreOffice 画平滑曲线；坐标轴默认可能隐藏 | `series.smooth=False`；`axis.delete=False` | 已验证 |
| 26 | openpyxl | 转 PDF 时表格被切成多页 | `fitToPage/fitToWidth=1/landscape` | 已验证 |
| 27 | matplotlib 3.11.2 | 新装字体不可见（旧缓存） | `fm._load_fontmanager(try_read_cache=False)` 或删 `fontlist-*.json` | 已验证 |
| 28 | matplotlib | 负号成方框 | `axes.unicode_minus=False` | 已验证 |
| 29 | matplotlib | Windows 需英文字体名 `Microsoft YaHei/SimHei` | 候选表按平台排列 | 未验证（无 Windows） |
| 30 | LibreOffice 24.2.7 | fontconfig 不认 `宋体/黑体/微软雅黑`，宋体黑体都被替成 Noto Sans | `FONTCONFIG_FILE=fonts.conf` 别名（注释里不能有 `--`） | 已验证 |
| 31 | LibreOffice | `--convert-to png` 只出第 1 页；`systemPathToFileUrl` 需绝对路径（相对路径 IOException） | 走 PDF→pdftoppm/PyMuPDF；`os.path.abspath` | 已验证 |
| 32 | LibreOffice | 并发/桌面冲突 | `-env:UserInstallation=file:///tmp/lo_profile_x` 独立配置目录 | 已用于 lo_convert.py |
| 33 | pandoc 3.1.3 | 无中文字体/题注域/编号，TOC 无页码 | 仅做应急，需 `--reference-doc` | 已验证 XML |
| 34 | docxtpl 0.20.2 | `{%tr for%}`/`{%tr endfor%}` 须各占一行 | 模板按文档写 | 已验证 |
| 35 | markdown-it-py 4.2.0 | 图片 src 被百分号编码（中文路径） | `urllib.parse.unquote` | 已验证 |
| 36 | PyMuPDF 2.0 | `import fitz` 弃用警告 | `import pymupdf` | 已验证 |
| 37 | 通用 | 中文文件名/路径 | Linux 全链路（python-docx/pptx/openpyxl/matplotlib/soffice/UNO）已验证可用；Windows 下 Python 3 str 路径亦可，避免依赖 `cv2.imread` 等不支持非 ASCII 路径的 C 库 | Linux 已验证 |
| 38 | 通用 | Python 文件 CRLF 行尾在 Linux 上 shebang 失效（`bash\r`） | 保存为 LF | 已踩坑 |

---

## 8. 附录

### 附录 A：一键复现日志（`run_all.sh`，12.8 s）

```bash
#!/usr/bin/env bash
# 一键复现：生成全部示例并转 PDF/PNG。前置：apt 装好 libreoffice-{writer,calc,impress} python3-uno fonts-noto-cjk poppler-utils
set -e
cd "$(dirname "$0")"
export FONTCONFIG_FILE="$PWD/fonts.conf"
PY=./venv/bin/python
echo "== versions =="; $PY --version; $PY -m pip list 2>/dev/null | grep -iE "^(python-docx|python-pptx|openpyxl|matplotlib|pandas|pillow|markdown-it-py|formulas|PyMuPDF|docxtpl|lxml|numpy|fonttools) " ; soffice --version | head -1; pandoc --version | head -1
echo "== 4 matplotlib =="; $PY 04_mpl_demo.py
echo "== 1 docx =="; $PY 01_docx_demo.py; /usr/bin/python3 lo_convert.py "out/01_商业计划书示例.docx" out/pdf
echo "== 2 pptx =="; $PY 02_pptx_demo.py; soffice --headless --convert-to pdf --outdir out/pptx_render "out/02_路演示例.pptx" >/dev/null; pdftoppm -png -r 50 "out/pptx_render/02_路演示例.pdf" out/pptx_render/slide; ls out/pptx_render/slide-*.png | wc -l
echo "== 3 xlsx =="; $PY 03_xlsx_demo.py; $PY 03b_xlsx_values.py 2>&1 | grep -E "^(A|B|C)\)"
echo "== 5 md2docx =="; $PY 05_md2docx_demo.py; /usr/bin/python3 lo_convert.py "out/05_md渲染示例.docx" out/md_render
echo "ALL DONE"
```

运行输出（去掉警告行）：

```text
== versions ==
Python 3.12.3
docxtpl           0.20.2
fonttools         4.67.0
formulas          1.3.4
lxml              6.1.3
markdown-it-py    4.2.0
matplotlib        3.11.2
numpy             2.5.4
openpyxl          3.1.5
pandas            3.1.0
pillow            13.0.0
pymupdf           2.0
python-docx       1.2.0
python-pptx       1.0.2
LibreOffice 24.2.7.2 420(Build:2)
pandoc 3.1.3
== 4 matplotlib ==
[bizcharts] 中文字体: Noto Sans CJK SC (Linux)
charts written: ['01_市场规模堆叠柱.png', '02_份额条形.png', '03_增长折线.png', '04_瀑布图.png', '05_漏斗.png', '06_龙卷风.png', '07_盈亏平衡.png', '08_现金流组合.png', '09_竞争矩阵.png', '10_甘特图.png', '11_商业模式画布.png', '_contact.png', '_contact2.png']
font: Noto Sans CJK SC
== 1 docx ==
saved /tmp/officetest/out/01_商业计划书示例.docx
docx (TOC updated): /tmp/officetest/out/pdf/01_商业计划书示例.updated.docx | indexes: 1
pdf: /tmp/officetest/out/pdf/01_商业计划书示例.pdf
== 2 pptx ==
saved /tmp/officetest/out/02_路演示例.pptx slides: 15
layout: Title Slide [(0, <PP_PLACEHOLDER_TYPE.CENTER_TITLE: 3>, 'Title 1'), (1, <PP_PLACEHOLDER_TYPE.SUBTITLE: 4>, 'Subtitle 2')]
layout: Title and Content [(0, <PP_PLACEHOLDER_TYPE.TITLE: 1>, 'Title 1'), (1, <PP_PLACEHOLDER_TYPE.OBJECT: 7>, 'Content Placeholder 2')]
saved /tmp/officetest/out/02_模板填充.pptx
Warning: failed to launch javaldx - java may not function correctly
15
== 3 xlsx ==
saved /tmp/officetest/out/03_财务模型.xlsx
FORMULAE check: {'NPV': True, 'IRR': True, 'SUMPRODUCT': True, 'IF': True, 'XLOOKUP': False, 'IFS': False, 'TEXTJOIN': False, 'COUNTIF': True}
A) LibreOffice 重算 0.9s  NPV=1851.02  IRR=0.3884  检验行=[0, 0, 0, 0, 0]
B) formulas 1.3.4 计算 0.8s  NPV=1851.0235082204508  IRR=0.3884136785442349  净利润2026=705.4499999999999  情景选择=基准 情景系数=1  加权NPV=1851.0235082204508
C) 数值版写入 240 个单元格 0.0s  NPV=1851.02  (pandas 可直接读)
== 5 md2docx ==
saved /tmp/officetest/out/05_md渲染示例.docx
pandoc rc: 0 
docxtpl: 智守康联 报名信息表 | 表格行数: 3 ['李华', 'CEO', '经管学院']
docx (TOC updated): /tmp/officetest/out/md_render/05_md渲染示例.updated.docx | indexes: 1
pdf: /tmp/officetest/out/md_render/05_md渲染示例.pdf
ALL DONE
```

### 附录 B：交付文件清单

| 路径 | 内容 |
| --- | --- |
| `04_assets/code/docx_helpers.py` | python-docx 工具库（字体、编号、域、目录、页眉页脚、封面、题注、表格、分节、超链接、书签、脚注、段落、列表、XML 顺序） |
| `04_assets/code/pptx_helpers.py` | python-pptx 工具库（主题、字体、10 类页面、4 类图表含组合图、形状、模板占位符、复制/删除页） |
| `04_assets/code/bizcharts.py` | matplotlib 商业图表库（字体自动选择、配色、11 类图） |
| `04_assets/code/md2docx.py` | Markdown→docx 渲染器 + `build_document()` 流水线 |
| `04_assets/code/lo_convert.py` | LibreOffice UNO：目录更新、xlsx 重算、PDF 导出（系统 python3） |
| `04_assets/code/fonts.conf` | Linux 字体别名（宋体/黑体/微软雅黑→Noto CJK） |
| `04_assets/code/01_docx_demo.py … 05_md2docx_demo.py`、`run_all.sh` | 演示与复现脚本 |
| `04_assets/run_all.log` | 全流程运行日志 |
| `04_assets/*.png` | 样图：5 张 matplotlib 图、docx 3 张页面渲染、pptx 2 张、xlsx 2 张、md2docx 1 张 |
| `/tmp/officetest/out/` | 生成的 `01_商业计划书示例.docx/.pdf`、`02_路演示例.pptx/.pdf`、`03_财务模型.xlsx`、`03_财务模型_数值版.xlsx`、`05_md渲染示例.docx/.pdf` 等（临时目录） |

### 附录 C：样图索引

- `04_assets/docx_page_正文.png`：python-docx 正文页（编号标题、粗体、脚注、列表、引用块、页眉线、页码）
- `04_assets/docx_page_示例.png`：图 + 题注、三线表、底纹合并表
- `04_assets/docx_目录与横向页.png`：UNO 更新后的目录页、横向宽表页
- `04_assets/pptx_全部页面缩略.png`：15 页路演 PPT；`pptx_组合图.png`：柱 + 折线次轴组合图
- `04_assets/xlsx_全部工作表.png`：9 张工作表；`xlsx_利润表与图表.png`：颜色约定与原生图表
- `04_assets/md2docx_渲染结果.png`：Markdown 样章渲染
- `04_assets/01_市场规模堆叠柱.png`、`04_瀑布图.png`、`06_龙卷风.png`、`08_现金流组合.png`、`11_商业模式画布.png`：matplotlib 样图
