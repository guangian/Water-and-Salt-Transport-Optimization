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

---

## A2B 模块：finmodel

### `bizkit.finmodel` — assumptions.json → 财务模型.xlsx / numbers.json / 财务摘要.md / 财务图表.charts.json（契约 §7.4 / §7.5）

```bash
python -m bizkit.finmodel 03_财务/assumptions.json                 # 输出到 assumptions.json 所在目录
python -m bizkit.finmodel 03_财务/assumptions.json --out 03_财务   # --out 相对工作区根
```

纯 Python 计算引擎（不依赖 Excel 重算），openpyxl 写表。`mode` 二选一：

| mode | 适用 | 工作表 |
|---|---|---|
| `startup` | 创业组：三表 + 盈亏平衡 + 单位经济 + 敏感性 + 估值 + 融资 | 说明 / 假设 / 驱动 / 收入 / 成本 / 人力 / 利润表 / 现金流量表 / 资产负债表 / 盈亏平衡 / 单位经济 / 敏感性 / 估值 / 融资 / 检查 |
| `budget` | 品牌 / 营销组：传播预算 + 排期 + KPI（可选 ROI） | 说明 / 预算总表 / 排期 / KPI（/ ROI） |

**startup 假设**（金额一律以 `currency` 元填写，`unit` 决定表内换算，默认 10000 = 万元；`years` 默认 5，`monthlyFirstYear` 首年按 12 个月、其后按年）：

```jsonc
{
  "mode": "startup", "company": "弘宝汝瓷文创", "unit": 10000, "startYear": 2026, "years": 5, "monthlyFirstYear": true,
  "scenarios": { "base": 1.0, "optimistic": 1.2, "pessimistic": 0.8 },            // 需求量倍数
  "revenueStreams": [
    { "name": "文创礼品线上零售", "driver": "volume_price", "unitName": "件", "cogsRatio": 0.42, "source": "单价=天猫同类 TOP10 均价 [12]",
      "volume": { "y1Monthly": [260, 290, …], "yoyGrowth": [0, 0.85, 0.55, 0.38, 0.26] },   // 或 volume.y1 = 首年合计（自动 S 形爬坡）
      "price":  { "base": 298, "yoyChange": [0, 0.03, 0.03, 0.03, 0.03] } },
    { "name": "企业定制", "driver": "clients_arpu", "clients": { "y1Monthly": [0, 0, 1, …], "yoyGrowth": […] }, "arpu": { "base": 60000 }, "cogsRatio": 0.55 }
  ],
  "costs": {
    "fixed": [ { "name": "场地租金", "monthly": 12000, "yoyChange": 0.05, "source": "…" } ],
    "variablePctOfRevenue": [ { "name": "平台佣金与支付", "pct": 0.06 } ],
    "marketing": { "pctOfRevenue": [0.16, 0.14, 0.11, 0.09, 0.08], "minAnnual": 300000 },   // 取 max(占比, 年度下限)
    "headcount": [ { "role": "品牌/运营", "y": [3, 4, 5, 6, 7], "monthlySalary": 8000, "socialRate": 0.28 } ],
    "capex": [ { "name": "设备与装修", "year": 1, "month": 1, "amount": 400000, "depreciationYears": 5 } ],
    "otherOpexPctOfRevenue": 0.03 },
  "workingCapital": { "receivableDays": 30, "inventoryDays": 60, "payableDays": 45 },
  "tax": { "incomeTaxRate": 0.25 },                                                 // 亏损结转：按累计利润计税
  "financing": [ { "round": "种子轮", "year": 1, "month": 1, "amount": 2000000, "preMoney": 8000000,
                   "use": [ { "item": "产品研发", "pct": 0.35 }, … ] } ],
  "openingCash": 500000,
  "valuation": { "methods": ["dcf", "ps_comparables", "vc_method"], "wacc": 0.18, "terminalGrowth": 0.03, "psMultiple": 3.0, "exitPE": 20, "exitYear": 5, "targetIrr": 0.40 },
  "unitEconomics": { "cac": 120, "repeatRate": 0.35, "grossMarginRef": 0.58, "lifetimeOrders": 3.2, "lifetimeMonths": 24 },
  "sensitivity": { "drivers": ["price", "volume", "cogsRatio", "cac", "fixedCost"], "range": [-0.3, -0.15, 0.15, 0.3] }
}
```

收入驱动 `driver`：`volume_price`（`volume` × `price`）、`clients_arpu`（`clients` × `arpu`，年合同额）、
`subscribers_arpu`（`subscribers{y1Monthly | y1, yoyGrowth, churnMonthly, opening}` 为新增订阅，按月滚动活跃数，收入 = 活跃订阅·月 × `arpu`（元/订阅/月））、
`funnel`（`traffic` × `conversion` × `aov`）。数量类字段给 `y1Monthly[12]` 或 `y1` 合计（自动 S 形爬坡）+ `yoyGrowth[years]`；
价格类（`price / arpu / aov / conversion`）给 `base` + `yoyChange[years]`。
敏感性 `drivers` 可选：`price volume cogsRatio cac fixedCost conversion marketing salary churn receivableDays`。

计算口径：毛利 → 营销 / 人力（人数 × 月薪 × (1+社保)）/ 固定 / 变动 / 其他费用 → EBITDA → 直线折旧 → EBIT → 所得税（累计计税基数取 running max，
即亏损结转、不退税）→ 净利润；营运资本按 DSO / DIO / DPO 天数折算期末余额，间接法现金流；融资按轮次注入实收资本。
盈亏平衡月 = 累计 EBIT 首次转正的月份（年度期间内线性插值），盈亏平衡销量按主收入流边际贡献倒推；
单位经济 LTV = 客单价 × 毛利率 × 生命周期订单数，回收期 = CAC ÷ 月毛利贡献；三情景用需求倍数重算全模型；
敏感性为单因素 ±range 对第 3 年净利润与最低现金的影响；估值 DCF（FCF = CFO + CFI，Gordon 终值）/ PS 可比 / VC 法三者取区间，
融资一致性 `equity = amount / (preMoney + amount)`。

**「公式 + 数值」方案**：假设 / 驱动 / 收入 / 成本 / 人力 五张表写 Python 数值（蓝字黄底 = 可改输入）；
利润表 / 现金流量表 / 资产负债表 / 检查写 **Excel 公式**（黑字 = 本表公式，绿字 = 跨表引用），列布局五张表完全一致
（A 科目 | B 期初 | C–N 首年 12 个月 | O 首年合计 `=SUM(C:N)` | P–S 第 2–5 年 | 右侧灰字「数值校验」列 = Python 算的年度值），
评委改任一蓝字，三表与「检查」表同步重算；「检查」表的 Master Check `=SUMPRODUCT(ABS(...))` 对公式结果与 Python 值逐项求差，
`<0.0001` 显示「通过」。所有表冻结首行首列、千分位 / 百分比格式、横向一页宽打印；利润表与现金流量表带 openpyxl 原生图表。

**budget 假设**：

```jsonc
{ "mode": "budget", "unit": 10000, "totalCap": 3000000, "startYear": 2026, "headcount": 6, "contingencyPct": 0.08,
  "phases": [ { "name": "第一阶段：声量引爆（2026.4–6）",                  // 或显式 "start": "2026-04", "end": "2026-06"
                "items": [ { "item": "抖音信息流投放", "category": "媒介投放", "unitPrice": 35, "quantity": 16000, "unit": "CPM", "source": "巨量引擎 2025 刊例 [18]" } ] } ],
  "kpis": [ { "metric": "曝光量", "baseline": "—", "target": "5000 万", "tool": "巨量云图", "cycle": "周", "phase": "第一阶段" } ],
  "roi": { "incrementalRevenue": 6800000, "grossMargin": 0.55 } }
```

「预算总表」每行 `金额 = 单价 × 数量` 为公式，阶段小计 / 合计 / 机动金 / 总计 / 上限余量 / 人均预算均为公式；「排期」为阶段 × 月份矩阵；
超过 `totalCap` 时 `warnings` 给出「超过上限」并在摘要中标红提示。

**输出**（都写到 `--out` 或 assumptions.json 所在目录）：

- `财务模型.xlsx`；
- `numbers.json`：`{ generatedAt, unit, mode, company, entries: [ { key, label, value, unit, source, usedIn } ] }`，
  `source` 形如 `03_财务/财务模型.xlsx!利润表!O2`，可直接用 openpyxl 定位核对。startup 固定 key：
  `revenue_y1..y5 net_profit_y1..y5 gross_margin_y1..y5 breakeven_month min_cash financing_amount pre_money post_money equity_offered cac ltv ltv_cac payback_months`
  （金额万元、比率 %、月数「月」）；budget：`budget_total budget_contingency budget_per_person kpi_1..n roi`；
- `财务摘要.md`：五年主表、盈亏平衡、单位经济、三情景、敏感性排名、估值区间、融资三元组、勾稽检查与「财务逻辑」段落；
- `财务图表.charts.json`：`bizkit.charts` 可直接渲染（`outputDir` = `04_图表`），startup 6 张
  `fig_fin_revenue_mix`(stacked_bar) `fig_fin_revenue_profit`(line) `fig_fin_cashflow`(cashflow_combo) `fig_fin_breakeven`(breakeven)
  `fig_fin_tornado`(tornado) `fig_fin_use_of_funds`(pie_as_bar)；budget 2 张 `fig_fin_budget_mix`(pie_as_bar) `fig_fin_budget_phases`(stacked_bar)。

结果 JSON（最后一行 stdout；`checks` 为三表勾稽的最大绝对差，单位元）：

```json
{ "ok": true, "mode": "startup",
  "files": { "xlsx": ".../03_财务/财务模型.xlsx", "numbers": ".../03_财务/numbers.json", "summary": ".../03_财务/财务摘要.md", "charts": ".../03_财务/财务图表.charts.json" },
  "checks": { "balanceSheet": 0.0, "cashflow": 0.0, "retained": 0.0 },
  "warnings": [ "客户生命周期默认 24 个月（可用 unitEconomics.lifetimeMonths 覆盖）" ],
  "breakevenMonth": 17, "minCash": -58.1, "charts": 6 }
```

验收规则（技能 fin-model §4，写进摘要「检查」段，不满足只警告不报错）：三表勾稽差 < 1 元；盈亏平衡 12–36 个月；
融资前现金最低点 + 融资额 ≥ 0；出让比例 = 融资额 ÷ 投后（±1 pp）；估值方法 ≥ 2 种；资金用途合计 = 100%；
首年亏损或微利（净利率 ≤ 10%）；LTV:CAC ≥ 3。另有软提示进 `warnings`：某年净利率 > 30%、单轮出让 > 20%、LTV:CAC > 5、
收入接近等比增长、条款 preMoney 偏离三法区间、融资后现金仍为负、预算项缺来源等。

示例与测试：`examples/ws_demo2/03_财务/assumptions.json`（弘宝汝瓷文创，两条收入流，盈亏平衡 17 个月、最低现金 −58.1 万、出让 20%）、
`assumptions_budget.json`（品牌传播预算 300 万上限，3 阶段 12 项，ROI 25%）。测试覆盖手算对照、三表勾稽、
xlsx 重新打开核对公式 / 配色 / 冻结窗格、numbers.json 的 source 可定位，并在装有 LibreOffice 时 headless 重算校验公式值与 Python 一致：

```bash
cd agent/toolkit && python -m pytest -q tests/test_finmodel.py
```

---

## A2C 模块：checks

### `bizkit.checks` — 机器检查器：工作区 → 07_评审/第N轮/机器检查.json + .md（契约 §7.7；标准 §6.6 / §6.7）

```bash
python -m bizkit.checks --track brand --workspace . --round 2          # 全量检查 → 07_评审/第2轮/机器检查.json|md
python -m bizkit.checks --workspace .                                  # 不给 --round → 07_评审/即时/；赛道从 project.json 推断
python -m bizkit.checks --track case --workspace . --variant retail     # 案例赛道变体：合并 requiredEvidenceByVariant / coverFieldsByVariant
python -m bizkit.checks --track brand --workspace . --chapter 05_正文/02_品牌现状分析.md   # 单章快检，只打印 JSON 不写文件
python -m bizkit.checks --workspace . --no-write                        # 只打印不落盘
```

纯规则 / 正则 / 表格启发式，不调用模型。阈值来自 `tracks/<id>/checks.json`（经 `track_defaults.load_track_config()` 补齐两套字段名），
硬门槛名称来自 `tracks/<id>/rubric.json.hardGates`，封面必填字段来自 `naming.json.coverFields`（变体可覆盖）。
**容错优先**：输入文件缺失 → 对应项 `status=skipped`；单项内部异常 → 该项 skipped 且 `summary.degraded=true`；
整个进程任何意外都以 `{"ok": false, "error": …}` 结束，不抛堆栈。stdout 最后一行恒为结果 JSON。

**读什么**：`project.json`、`05_正文/*.md` + `book.json`（封面 / 附录清单）、`01_调研/参考文献.json`、`03_财务/numbers.json`（账本）与 `财务模型.xlsx`
（「检查」「预算总表」「敏感性」页）、`04_图表/index.json`（图表来源）、`06_路演/讲稿.md`、`答辩题库.md`、`slides.json`（或 `产出/路演.pptx.preview.json`）、
`产出/路演.render.json`（`pptx_render` 的字号统计）、`产出/*.pdf`（真实页数，pypdf）。

**检查项**（每项 `{id, value, threshold, status: pass|warn|fail|skipped, detail, locations?[{file, line, text, …}]}`，`locations` 带文件 / 行号或页码 / 原句）：

| 组 | id | 规则要点 |
|---|---|---|
| 结构 / 篇幅 | `struct_coverage` `page_count` `word_count` | 一级标题 / 文件名 ↔ `requiredElements`；页数优先读 PDF，否则 字数/650 + 图×0.45 + 表行/40；超 `pageLimit` fail、偏离 `pageTarget±tol` warn |
| 匿名 | `anonymity_hits` | `anonymityRegex` − `anonymityWhitelist`（等长占位遮蔽，命中位置映射回原句），扫正文 / 路演 md / slides.json / 封面 |
| 数据 | `numeric_density` | 每千字数值 token（含 % / 万 / 亿 / 元 / 倍），去表格分隔线与图片路径；全文 + `numericDensityChapters` 逐章阈值 |
| 来源 | `source_count` `source_domains` `source_official_ratio` `source_freshness` `source_pairing` `chart_count` | 参考文献条数 / 独立域名 / 官方占比（`OFFICIAL_DOMAIN_HINTS` + `OFFICIAL_TYPES`）/ 近 N 年；图表后 3 行内「来源」或 index.json 有 source |
| 用语 | `forbidden_terms` `forbidden_terms_adlaw` `weak_terms` | 禁语（needSource 词同句有「来源 / 据《 / [n]」放过）；广告法词减白名单，「第一」排除序数与列举，「100%」只认宣传搭配；弱词每百字 |
| 示例数据 | `synthetic_marker_coverage` | numbers.json `synthetic:true` 的值，在正文 / 讲稿 / 题库每次出现 ±3 行（表格整表、PPT 同页）内须有「示例数据」 |
| 路演 | `slide_metrics` `deck_font_min` `script_timing` `sentence_length` `qa_bank` | 正文页数区间（不含 appendix）、每页 ≤ 120 字、结论句标题 ≥ 80%、含数字页 ≥ 80%、图表页 ≥ 8；render.json 最小字号；讲稿字数 ÷ 190；句长中位数 ≤ 25 / 最长 ≤ 40；`### Qn【维度】` 题数、每维度 ≥ 4、主答 + 结论 / 证据 / 行动 |
| 赛道专项（`requiredEvidence` 开关） | `has_primary_research` `has_tam_sam_som` `has_competitor_matrix` `has_unit_economics` `has_sensitivity` `three_statement_check` `financing_consistency` `revenue_shape` `budget_recalc` `budget_in_range` `kpi_table` `positioning_statement` `big_idea` `kol_tiers` `offline_event` `risk_matrix` `brief_response_table` `appendix_questionnaire` `honor_ratio` `title_pattern` `ip_list` `param_compare_table` `milestone_gantt` `equity_table` `team_matrix` `short_video_script` `platform_matrix` `ratio_table_with_formula` `comprehensive_income_forecast` `annual_report_page_refs` `contingency` `itinerary_timeline` `itinerary_grid` `venue_capacity_check` `vendor_quote_table` `regional_feature` … | 竞品矩阵 = ≥ 3 个非本品实体列 × ≥ 5 行或题注含「竞品」；KPI 表 = 指标 / 基线 / 目标 / 工具 / 周期 ≥ 3 列 × ≥ 5 行；预算复算 = Σ行 vs 合计 vs 账本 ≤ 1%、机动金 5–10%；三表勾稽读 numbers.json 或 xlsx「检查」页；`brief_response_table` 仅 `project.subject.isAssignedCase` 且赛道 `if_assigned|true` 时检查，否则 skipped |
| 跨材料 | `cross_doc_consistency` | 账本（`numberKeys`）↔ 正文 / PPT / 讲稿 / 题库：同句出现 label 关键词时比对同类数字（±1%，万 / 亿换算，只看同句 / 同表行，泛化词只看紧邻数字，情景 / 阶段 / 机动金句不报） |

硬门槛按 `rubric.hardGates` 逐个输出 `{id, passed, skipped?, detail, …}`：`anonymity_hits` `page_limit` `word_min` `cover_fields` `cross_doc_consistency`
`financing_consistency` `forbidden_terms_adlaw` `budget_in_range`；对应输入缺失时 `passed: true, skipped: true`，不误伤。

结果 JSON（同时写入 `机器检查.json`；`.md` 为检查项表 + 硬门槛表 + 跨材料不一致清单 + 问题定位 + 一句总评）：

```json
{ "track": "brand", "variant": null, "round": 1, "generatedAt": "…", "workspace": "…", "config": { "source": "…/tracks/brand/checks.json", "warning": null },
  "checks": [ { "id": "anonymity_hits", "value": 1, "threshold": { "pass": 0 }, "status": "fail", "detail": "1 处疑似院校 / 人员信息：05_正文/02_品牌现状分析.md:15「大学」",
                "locations": [ { "file": "05_正文/02_品牌现状分析.md", "line": 15, "term": "大学", "text": "本次问卷由华中某大学市场营销专业的同学协助发放。" } ] }, "…" ],
  "hardGates": [ { "id": "anonymity_hits", "passed": false, "detail": "…", "hits": [ "…" ] }, { "id": "cover_fields", "passed": true, "required": [ "作品名称", "团队名称" ], "missing": [] }, "…" ],
  "crossDoc": { "consistent": false, "checked": 4, "keys": [ "budget_total", "sample_size", "kpi_exposure", "roi" ],
                "mismatches": [ { "key": "budget_total", "label": "推广预算总额（含机动金）", "unit": "万元", "expected": 300, "book": 300.0, "deck": 300.0, "script": 280.0, "qa": 300.0,
                                  "locations": [ { "file": "06_路演/讲稿.md", "line": 29, "text": "传播预算 280 万元，分六块。", "found": "280万元", "value": 280.0, "doc": "script" } ] } ] },
  "summary": { "pass": 21, "warn": 5, "fail": 8, "skipped": 1, "degraded": false, "hardGatePassed": false },
  "verdict": "硬门槛未过（anonymity_hits、cross_doc_consistency、forbidden_terms_adlaw），不可提交；另有 8 项 fail、5 项 warn。",
  "errors": [], "files": { "json": "…/07_评审/第1轮/机器检查.json", "md": "…/07_评审/第1轮/机器检查.md" }, "ok": true }
```

单章模式（`--chapter`）只跑 `numeric_density source_pairing forbidden_terms forbidden_terms_adlaw anonymity_hits weak_terms` 六项，
结果多一个 `"mode": "chapter"`，不写文件；写作智能体可在每章落盘后立即调用。

示例与测试：`examples/ws_demo3/`（云岫茶坊，品牌赛道命题案例的迷你工作区，**故意**含一处匿名违规「华中某大学」、一处禁语「蓝海」、
一处广告法用语「最佳」、一处缺「示例数据」标记的样本量、讲稿里一处与账本不一致的预算 280 vs 300；同时含竞品矩阵 / KPI 表 / 预算表 /
风险矩阵 / 命题回应表 / 定位陈述 / 26 题问卷），`tests/test_checks.py` 断言这些缺陷逐一命中、缺文件全部 skipped 不崩、
对 `examples/ws_demo`（brand，读 PDF 真页数）与 `ws_demo2`（startup，读 xlsx「检查」页、发现 PPT 两处数字不一致）整体运行无 errors、单章模式与 CLI 入口、案例赛道变体合并：

```bash
cd agent/toolkit && python -m pytest -q tests/test_checks.py
```
