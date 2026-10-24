# bizkit 工具包

金策 BizGold 的 Python 产出工具。所有命令都在 `agent/toolkit/` 下以模块方式运行：

```bash
cd agent/toolkit
python -m bizkit.<模块> ...        # Windows 用 py -3 / 项目 venv 的 python
```

约定：每个 CLI 把一段 JSON 结果打印到 stdout（`ok` 字段表示是否成功），人类可读的进度写到 stderr/`[bizkit]` 前缀行；
路径一律相对**工作区根**（含 `任务书.md` / `CLAUDE.md` / `project.json` 的目录），绝对路径原样使用。
运行环境：Python 3.10+，依赖见 `requirements.txt`；测试 `python -m pytest -q tests/`。

---

## A1 模块：charts / docx_render / preread / pdf / package

### `bizkit.charts` — charts.json → 300 dpi PNG

```bash
python -m bizkit.charts 04_图表/charts.json              # 全部
python -m bizkit.charts 04_图表/charts.json --only fig_a,fig_b
```

退出码：0 全部成功；2 部分失败（失败项在 `failed`，其余照常输出）；1 charts.json 本身读不了。
输出：`outputDir/<id>.png` + `outputDir/index.json`（id → file/title/source/type/widthMm/heightMm/px）。
图内不画标题（题注由 Word 排版），资料来源写在图左下角；所有图型都保证黑白打印可辨（填充纹理 / 线型 / 线末标签）。

支持 18 种 `type`：`bar` `grouped_bar` `stacked_bar` `hbar` `line` `area` `pie_as_bar` `donut`（>4 类自动退化为 `pie_as_bar` 并警告）
`waterfall` `funnel` `tornado` `breakeven` `cashflow_combo` `scatter_matrix` `positioning_map` `gantt` `radar` `bmc`。

```jsonc
{
  "theme": "bizgold", "dpi": 300, "outputDir": "04_图表",
  "charts": [
    { "id": "fig_market_size", "type": "stacked_bar", "title": "汝瓷工艺品市场规模", "width": 150, "height": 85,
      "yLabel": "市场规模", "unit": "亿元", "source": "资料来源：中国陶瓷工业协会（2025）", "note": "2026E 为预测值",
      "data": { "categories": ["2022", "2023", "2024", "2025E", "2026E"],
                "series": [ { "name": "线下", "values": [18.2, 19.5, 20.1, 21.0, 21.8] },
                            { "name": "线上", "values": [5.1, 7.2, 9.6, 12.5, 15.9] } ] } },
    { "id": "fig_funnel", "type": "funnel", "title": "传播转化漏斗", "unit": "万",
      "data": { "stages": [ { "name": "曝光", "value": 5200 }, { "name": "点击", "value": 420 }, { "name": "下单", "value": 4.7 } ] } },
    { "id": "fig_gantt", "type": "gantt", "title": "执行排期",
      "data": { "tasks": [ { "name": "产品焕新", "start": "2026-04", "end": "2026-07", "group": "产品" } ],
                "milestones": [ { "name": "首发上市", "date": "2026-08" } ] } },
    { "id": "fig_bmc", "type": "bmc", "title": "商业模式画布",
      "data": { "keyPartners": ["…"], "keyActivities": ["…"], "keyResources": ["…"], "valuePropositions": ["…"],
                "customerRelationships": ["…"], "channels": ["…"], "customerSegments": ["…"],
                "costStructure": ["…"], "revenueStreams": ["…"] } }
  ]
}
```

各图型的 `data` 字段在 `bizkit/charts.py` 顶部 `DATA_SPECS` 一处定义（校验失败时会原样附在错误里），摘要：
`waterfall {items:[{name, value, total?}]}`、`tornado {base, items:[{name, low, high}]}`、
`breakeven {fixedCost, unitPrice, unitVarCost, maxUnits?, unitName?}`、`cashflow_combo {categories, values, cumulative?}`、
`scatter_matrix {xLabel, yLabel, points:[{name, x, y, size?, self?}], quadrants?, xMid?, yMid?}`、
`positioning_map {xAxis:{low, high}, yAxis:{low, high}, points:[…]}`、`radar {axes:[…≥3], series:[{name, values}], max?}`、
`hbar {categories, values, sort?, highlight?}`、`donut {categories, values, centerText?}`。
图表级可选字段（与 `data` 同级）：`unit`（写进坐标轴标题与数值标签）、`xLabel/yLabel`、`valueLabels`、`legend`、
`decimals`（0–4）、`width/height`（40–400 mm，未给时按图型取默认尺寸）、`source`、`note`。
字段缺失或类型错误时返回中文提示，例如 `[fig_x] stacked_bar 需要 data.series（数组）`。

结果 JSON：

```json
{ "ok": true, "outputDir": ".../04_图表", "index": ".../04_图表/index.json", "count": 19,
  "rendered": [ { "id": "fig_market_size", "file": "fig_market_size.png", "type": "stacked_bar", "px": [1772, 1004], "warnings": [] } ],
  "failed": [ { "id": "fig_bad", "type": "line", "error": "[fig_bad] line 需要 data.series（数组）" } ],
  "warnings": [ "[fig_motivation] donut 有 6 类（>4），已自动改为 pie_as_bar" ], "generatedAt": "2026-03-01" }
```

### `bizkit.docx_render` — book.json + Markdown → Word（可选 PDF）

```bash
python -m bizkit.docx_render 05_正文/book.json            # book.json.pdf 为 true 时顺带转 PDF
python -m bizkit.docx_render 05_正文/book.json --no-pdf
```

```jsonc
{
  "cover": { "title": "宋韵落天青，一釉生万象", "subtitle": "弘宝汝瓷品牌策划书", "teamName": "青瓷工作室",
             "company": "指导企业：汝州市弘宝汝瓷坊有限公司", "group": "本科组", "date": "2026 年 3 月" },
  "abstract": "05_正文/00_概要.md",
  "chapters": [ { "file": "05_正文/01_引言.md" }, { "file": "05_正文/02_品牌现状分析.md" } ],
  "appendices": [ { "file": "05_正文/A1_调研问卷.md" } ],
  "references": "01_调研/参考文献.json",
  "style": { "bodyFont": "宋体", "bodySize": 12, "headingFont": "黑体", "latinFont": "Times New Roman",
             "lineSpacing": 1.5, "firstLineIndent": 2, "toc": true, "tocDepth": 2,
             "header": "弘宝汝瓷品牌策划书", "tableStyle": "banded", "h1Align": "center" },
  "anonymize": { "forbidden": ["某某大学"], "whitelist": ["河南财经"] },
  "output": "产出/弘宝汝瓷品牌策划书.docx",
  "pdf": true
}
```

Markdown 子集：`#`–`###` 标题（作者写的 “1.2 / 第一章” 编号会被剥掉并统一重编）、段落、`**粗体**`/`*斜体*`、
`- ` / `1. ` 列表、`> ` 引用、`| |` 表格（上一行 `表：题注`，表后一行 `{资料来源：…}`）、`![题注](路径)`、
`[1]` / `[2-4]` 文中引用（上标）、`[^k]: 内容` 脚注（归为章末注释）、`---` 或 `<!-- pagebreak -->` 分页、`{资料来源：…}` 来源行。
图宽上限 15.5 cm；图/表题注按章自动编号（图 2-1、表 A-1）；一级标题另起一页；封面不含学校信息。
目录是真正的 TOC 域（Word 打开时 `updateFields` 自动刷新），同时写入静态目录文字——若本机能转 PDF，页码从 PDF 回填（`tocPageNumbers: "pdf"`），否则为估算值。

参考文献 JSON 为数组，按 `type` 输出 GB/T 7714 格式：`report`→[R/OL]、`web`→[EB/OL]、`gov`→[A/OL]、`paper`→[J]、
`book`→[M]、`annual-report`→[R]、`news`→[N]、`standard`→[S]、`thesis`→[D]、`database`→[DB/OL]。
字段（见 contracts §7.2）：`id type author title publisher year url accessed note`，另可选 `place date journal volume issue pages
newspaper edition number`；有 `url` 即按联机文献加 `/OL` 与访问日期。

结果 JSON：

```json
{ "ok": true, "output": ".../产出/弘宝汝瓷品牌策划书.docx", "chapters": 3, "figures": 19, "tables": 6, "tableRows": 38,
  "words": 3209, "pagesEstimate": 14.4, "pages": 25, "pdf": ".../产出/弘宝汝瓷品牌策划书.pdf", "converter": "libreoffice",
  "tocPageNumbers": "pdf", "anonymizeHits": [ { "file": "05_正文/02_品牌现状分析.md", "line": 80, "term": "某某大学", "count": 1 } ],
  "missingImages": [], "warnings": [] }
```

`pagesEstimate = 字数/650 + 图数×0.45 + 表行数/40`；缺图处插入红色「【缺图：路径】」占位并在 `warnings` 提示，不中断渲染。

### `bizkit.preread` — 资料预读摘要

```bash
python -m bizkit.preread inputs/企业资料.pdf                 # → inputs/_preread/企业资料.pdf.txt
python -m bizkit.preread inputs/销售数据.xlsx --out 01_调研/预读
```

支持 pdf（文本 + 表格，pdfplumber）、docx（段落 + 表格）、xlsx/csv（每个 sheet 的形状、列类型、前 N 行、数值列统计）、
txt/md/json；图片只记录尺寸。结果：

```json
{ "ok": true, "file": ".../示例企业资料.pdf", "kind": "pdf", "out": ".../_preread/示例企业资料.pdf.txt",
  "words": 486, "pages": 2, "tables": 1, "warnings": [] }
```

### `bizkit.pdf` — docx → pdf

```bash
python -m bizkit.pdf 产出/策划书.docx [--out 产出/策划书.pdf] [--prefer libreoffice]
python -m bizkit.pdf --pages 产出/策划书.pdf          # 只数页数
python -m bizkit.pdf --detect                          # 列出本机可用转换器
```

转换器顺序：Word COM（会先更新域）→ WPS COM → docx2pdf → LibreOffice（独立 profile，无头）。
结果 `{ "ok": true, "pdf": "...", "converter": "libreoffice", "pages": 25, "tried": [], "seconds": 4.2 }`（`tried` 列出之前失败的转换器及原因）；
全部失败时 `{ "ok": false, "reason": "所有转换器均失败", "tried": [ { "converter": "libreoffice", "error": "..." } ] }`。
`docx_render` 在 `pdf: true` 时内部调用它，转换失败只记一条警告，docx 照常产出。

### `bizkit.package` — 提交包

```bash
python -m bizkit.package --track brand --workspace .      # → 08_提交包/ + <命名>.zip
```

命名规则取 `tracks/<track>/naming.json`（先找工作区下的，再找 `agent/tracks/` 下的），都没有时用内置默认
`{组别}-{学校}-{团队名}-{作品名}`；`{学校}` 等需手填的占位符
保留为「学校名称」字样并在 README.md 中提示替换。复制并重命名 `产出/*.docx|pdf|pptx`、`03_财务/财务模型.xlsx`、
`06_路演/{讲稿,答辩题库,视频分镜}.md`、`01_调研/调研工具包/`、最新一轮 `07_评审/`；生成 `声明与授权/` 三份模板
（注明以组委会当年附件为准）、`数据真实性清单.md`（`numbers.json` 中 `synthetic` 条目 + 正文「【示例数据」标记）、
`README.md`（文件清单 + 提交前检查）。

```json
{ "ok": true, "dir": ".../08_提交包", "zip": ".../08_提交包/本科组-学校名称-青瓷工作室-宋韵天青……策划书.zip",
  "baseName": "本科组-学校名称-青瓷工作室-宋韵天青——弘宝汝瓷品牌策划书", "naming": "default",
  "files": [ { "src": "产出/弘宝汝瓷品牌策划书.docx", "dst": "本科组-学校名称-….docx", "size": 2491268, "note": "品牌策划书（Word 源文件）" } ],
  "missing": [ "03_财务/财务模型.xlsx" ], "placeholders": [ "学校名称" ], "pdfPages": 25, "syntheticCount": 2, "warnings": [] }
```

### 示例与测试

`examples/ws_demo/` 是一个覆盖全部 18 种图型、3 章 + 1 附录 + 7 条参考文献的迷你工作区（弘宝汝瓷品牌策划书）；
`examples/make_demo_inputs.py` 生成 `inputs/` 下的示例 PDF / xlsx。端到端测试：

```bash
cd agent/toolkit && python -m pytest -q tests/test_a1.py
```

---

## A2 模块：pptx_render

### `bizkit.pptx_render` — slides.json → 16:9 路演 PPT（契约 §7.3）

```bash
python -m bizkit.pptx_render 06_路演/slides.json
```

从空白 16:9 母版用绝对坐标生成，不依赖模板。主题取 `theme.py`（海军蓝主色、金色强调、白底；cover / section / closing
深蓝底白字 + 金线），字体名按 `fonts.doc_fonts()` 写入主题与每个 run（中文微软雅黑 `a:ea`、西文 Arial `a:latin`），
由打开它的机器解析，Linux 生成机不装这些字体也没关系（本机预览用 LibreOffice 会回退到 Noto CJK）。
路径约定：`image` 相对 **slides.json 所在目录**（`../04_图表/fig_x.png`），`output` 相对工作区根。

```jsonc
{ "title": "汝瓷文创创业计划路演", "subtitle": "青瓷工作室 · 创业组", "aspect": "16:9", "theme": "bizgold",
  "minutes": 10, "footer": "汝瓷文创创业计划 · 青瓷工作室",
  "anonymize": { "forbidden": ["某某大学"], "whitelist": [] },
  "slides": [ { "layout": "cover", "title": "…", "subtitle": "…", "note": "讲稿", "duration": 20 }, … ],
  "output": "产出/路演.pptx" }
```

页型与字段（每页通用 `note` = 演讲者备注，`duration` 秒会附在备注末尾）：

| layout | 必填 | 可选 | 版式 |
|---|---|---|---|
| `cover` | title | subtitle | 深蓝整页 + 金色竖线，标题 44pt（放不下时在标点处平衡断两行） |
| `agenda` | title, items[] (≤8) | | 金色圆形序号 + 条目（>4 项两列），右侧深蓝面板放演示标题 |
| `section` | title | kicker | 深蓝整页，自动章节序号 01/02… |
| `bullets` | title, bullets[] (≤6，每条 ≤40 字) | highlight, highlightLabel | 金色方块项目符号 24→18pt 自适应；有 highlight 时右侧大数字卡片 |
| `image` | title, image | caption, bullets[] (≤3) | 图左文右（无 bullets 时整宽），按比例适配，caption 14pt 灰字 |
| `two_column` | title, left{heading,bullets}, right{heading,bullets} | | 两张浅灰卡片，左蓝右金色条 |
| `chart` | title, chart{type, categories[], series[{name,values}], unit} | insight, decimals | python-pptx 原生图表；insight 放图下金边结论框 |
| `table` | title, columns[] (≤7), rows[][] (≤8) | | 表头海军蓝白字、斑马纹、自动列宽、数字列右对齐、18→14pt |
| `kpi` | title, kpis[{value,label}] (2–4) | | 大数字卡片 44→28pt（必须一行放下） |
| `timeline` | title, steps[{when,what}] (3–6) | | 水平时间轴，蓝金交替序号圆点 |
| `team` | title, members[{role,duty}] (≤6) | | 卡片 + 角色首字徽标，>3 人两行 |
| `quote` | title, quote | by | 大号金色引号 + 36→24pt 引言 |
| `closing` | title | subtitle, contact | 深蓝整页 |
| `appendix` | title | image / bullets / columns+rows / chart 字段同上 | 右上角金色「备用页」角标，内容按字段自动选版式 |

`chart.type`：`column_clustered`（竖柱）、`bar`（横条，第一类在最上）、`line`（带标记折线）、`pie`（≤4 类，超出自动改竖柱）、
`stacked`（堆叠柱，深色段白字标签）；别名 `column / grouped_bar / hbar / stacked_bar / donut` 也接受。
系列色取 `theme.series`，数据标签与轴标签 14pt，图例在底部（单系列不显示图例），`unit` 作数值轴标题。

每页固定元素：顶部标题（结论句，32pt，放不下缩到 28pt，再长换两行）、金色短线、页脚（`footer` + 页码 `n / N`）、演讲者备注。

渲染前校验（超出先警告，再截断或缩字号，保证不溢出；**警告应清零——精简文字而不是接受缩字号**）：

| 检查 | 处理 |
|---|---|
| 标题 > 28 字 | 警告 |
| 正文 > 120 字（表格页按行列数约束，不计） | 警告 |
| 要点 > 6 条 / 单条 > 40 字 | 警告；多余条目截去；放不下时字号 24→18pt，仍放不下按行数均分截断加「…」 |
| 表格 > 7 列 / > 8 行 | 警告并截去；字号 18→14pt |
| 饼图 > 4 类、未知 chart.type、系列值个数与类目不符 | 警告；改簇状柱 / 对齐补空 |
| kpi 不在 2–4、timeline 不在 3–6、team > 6、agenda > 8 | 警告；多余截去 |
| 缺图 / 未知 layout | 画「【缺图：路径】」占位框并计入 `missingImages` / 按 bullets 渲染 |
| `anonymize.forbidden` 命中（含备注） | 替换为「XX」，记入 `anonymizeHits`；`whitelist` 中的词受保护 |

输出：`<output>`（pptx）与 `<output>.preview.json`——数组，每页
`{ index, layout, title, texts[], image, chart, table, kpis, steps, members, quote, note }`，无关字段为 `null`，
`image` 为相对工作区根的路径（`04_图表/fig_market_size.png`），前端用它拼 HTML 预览。

结果 JSON：

```json
{ "ok": true, "output": ".../产出/路演.pptx", "preview": ".../产出/路演.pptx.preview.json", "pages": 20,
  "warnings": [ "第 5 页：要点 7 条 > 6，已保留前 6 条" ],
  "anonymizeHits": [ { "page": 17, "field": "members[0].duty", "term": "某某大学", "count": 1 } ],
  "missingImages": [], "durationSeconds": 540,
  "perPage": [ { "index": 1, "layout": "cover", "chars": 10, "minFont": 16 } ] }
```

`perPage.chars` 为该页正文字数（去空白字符数），`minFont` 为该页正文 / 标签 / 表格 / 图表的最小字号（不含页脚页码），
下限 14pt；正文要点 ≥ 18pt。本机预览：`soffice --headless --convert-to pdf 产出/路演.pptx` 再 `pdftoppm -png -r 80`。

示例与测试：`examples/ws_demo2/`（汝瓷文创创业计划路演，20 页覆盖全部 14 种页型与 5 种图表类型，
图表 PNG 由 `04_图表/charts.json` 经 `bizkit.charts` 生成）：

```bash
cd agent/toolkit && python -m pytest -q tests/test_pptx.py
```
