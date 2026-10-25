# 金策 BizGold · 接口与数据契约（v1.0）

> 所有模块（服务端、前端、工具包、智能体规范）以本文件为准。字段一律 camelCase（JSON）；文件路径一律相对项目工作区根 `workspace/<projectId>/`；编码 UTF-8；时间 ISO-8601。

## 1. 设置 `~/.bizgold/settings.json`

```jsonc
{
  "api": {
    "baseUrl": "https://llm.goaichat.top",        // Anthropic 格式根地址（不带 /v1）
    "apiKey": "",                                  // 为空则读环境变量 BIZGOLD_API_KEY / ANTHROPIC_AUTH_TOKEN
    "fallbackBaseUrl": "https://model.kimitk.top", // 连不上时切换（服务端健康检查决定）
    "models": {                                    // 角色 → 模型 id（网关 /v1/models 可查）
      "main": "deepseek-v4-pro",
      "researcher": "kimi-k3",
      "finance": "deepseek-v4-pro",
      "writer": "deepseek-v4-pro",
      "pitch": "deepseek-v4-pro",
      "judgeA": "kimi-k3",
      "judgeB": "glm-5.3",
      "judgeC": "qwen3.8-max",
      "inspector": "deepseek-v4-pro",
      "fast": "deepseek-v4-flash"                  // 后台辅助调用（ANTHROPIC_DEFAULT_HAIKU_MODEL）
    },
    "contextTokens": 200000,
    "disableThinking": false
  },
  "search": {
    "providers": ["bocha", "tavily", "serper", "bing", "duckduckgo"],  // 尝试顺序；无 key 的 API 自动跳过
    "bochaKey": "", "tavilyKey": "", "serperKey": "",
    "maxResults": 8
  },
  "python": { "path": "", "pipIndex": "https://pypi.tuna.tsinghua.edu.cn/simple" },
  "run": { "maxTurns": 1500, "totalTimeoutMin": 360, "idleTimeoutMin": 45, "bashTimeoutMs": 600000,
           "reviewRounds": 3, "targetScore": 85 },
  "paths": { "workspaceDir": "" },                 // 为空 = <repo>/workspace
  "ui": { "theme": "light", "language": "zh-CN" }
}
```

`GET /api/settings` 返回时 `apiKey` 以 `"sk-****last4"` 形式脱敏，并附 `apiKeySource: "settings" | "env" | "none"`。`PUT /api/settings` 中 `apiKey` 若为脱敏串则不改动。

## 2. 项目 `workspace/<id>/project.json`

```jsonc
{
  "id": "p_20261024_153012_ab12",
  "name": "宋韵天青——弘宝汝瓷品牌策划书",      // 作品名
  "teamName": "青瓷工作室",                      // 团队名（匿名制下代替院校）
  "track": "brand",                              // brand | startup | digital | case | mice
  "group": "本科组",                              // 组别：研究生组/本科组/高职组/命题组/创业组/创客组/在华留学生组
  "subject": {                                   // 策划对象 / 项目
    "type": "brand",                             // brand | startup | case-company | event
    "name": "弘宝汝瓷",
    "company": "汝州市弘宝汝瓷坊有限公司",
    "isAssignedCase": false,                     // 是否命题企业/命题组
    "caseVariant": "accounting",                 // 仅 case 赛道：accounting（会计案例）| retail（流通业）
    "caseTitle": "",                             // 命题题目（如有）
    "brief": "……用户填写的项目/品牌/命题简介……",
    "region": "河南汝州",
    "industry": "非遗工艺品 / 文创"
  },
  "team": [ { "role": "队长/总策划", "major": "市场营销", "strength": "品牌传播" } ],  // 匿名信息，无姓名院校
  "params": {
    "researchMode": "synthetic",                 // synthetic（示例数据，标注）| provided（仅用用户资料）
    "pageTarget": 60,                            // 正文目标页数（≤80）
    "pitchMinutes": 10,                          // 讲稿时长
    "reviewRounds": 3,
    "targetScore": 85,
    "language": "zh-CN",
    "models": { "main": "deepseek-v4-pro", "...": "…" }   // 快照自设置
  },
  "inputs": [ { "name": "企业资料.pdf", "path": "inputs/企业资料.pdf", "size": 123456, "kind": "pdf", "preread": "inputs/_preread/企业资料.pdf.txt" } ],
  "status": "idle",                              // idle | running | paused | completed | failed | stopped
  "currentStage": "write",
  "stages": {                                    // 9 阶段
    "parse":    { "status": "done",    "startedAt": "…", "endedAt": "…", "note": "" },
    "research": { "status": "done" },
    "strategy": { "status": "done" },
    "finance":  { "status": "done" },
    "charts":   { "status": "done" },
    "write":    { "status": "running" },
    "pitch":    { "status": "pending" },
    "review":   { "status": "pending", "round": 0 },
    "package":  { "status": "pending" }
  },
  "session": { "sessionId": "…", "lastPromptAt": "…", "mode": "run" },   // mode: run | modify | resume
  "usage": { "inputTokens": 0, "outputTokens": 0, "cacheReadTokens": 0, "turns": 0, "costEstimateUsd": 0, "durationMs": 0 },
  "review": { "rounds": [ { "round": 1, "score": 78.5, "byJudge": { "investor": 76, "professor": 80, "executive": 79.5 },
                             "hardGatePassed": false, "path": "07_评审/第1轮" } ], "latestScore": 78.5 },
  "artifacts": { "book": "产出/弘宝汝瓷品牌策划书.docx", "pdf": "产出/弘宝汝瓷品牌策划书.pdf", "deck": "产出/路演.pptx",
                 "model": "03_财务/财务模型.xlsx", "script": "06_路演/讲稿.md", "qa": "06_路演/答辩题库.md", "package": "08_提交包/提交包.zip" },
  "createdAt": "…", "updatedAt": "…"
}
```

阶段 `status`：`pending | running | done | failed | skipped`。

## 3. REST API（前缀 `/api`，JSON；错误统一 `{ "error": { "code": "…", "message": "…" } }`）

| 方法 | 路径 | 说明 |
|---|---|---|
| GET | `/health` | `{ ok, version, python: {found, path, version, missing: [pkg]}, pdfConverter: "word"|"wps"|"libreoffice"|null, fonts: {cjk: string[]}, gitBash: bool, workspaceDir }` |
| GET / PUT | `/settings` | 见 §1 |
| POST | `/settings/test` | body `{ baseUrl?, apiKey? }` → `{ ok, latencyMs, models: string[], error? }`（请求 `/v1/models`，失败时尝试 fallback） |
| GET | `/tracks` | `[ { id, name, docName, pageLimit, wordMin, scoring: {written, pitch}, elements: string[], groups: string[], assignedCases: [{company, title}] } ]` |
| GET | `/knowledge` / `/knowledge/:slug` | 知识库目录与 Markdown 内容 |
| GET | `/projects` | 项目列表（摘要） |
| POST | `/projects` | body = project.json 的 `name/teamName/track/group/subject/team/params` → 创建工作区并返回完整 project |
| GET / DELETE | `/projects/:id` | 详情 / 删除（删除需 `?confirm=1`） |
| PATCH | `/projects/:id` | 修改元数据/参数（运行中仅允许改 `name`） |
| POST | `/projects/:id/files` | multipart 上传到 `inputs/`，返回 inputs 数组（自动预读） |
| DELETE | `/projects/:id/files?name=` | 删除输入文件 |
| POST | `/projects/:id/run` | body `{ mode: "run"|"resume", fromStage?: string }` → `{ ok }`；运行中返回 409 |
| POST | `/projects/:id/stop` | 中止（abort） |
| POST | `/projects/:id/chat` | body `{ message }` → 修改模式（resume 会话），事件走 SSE |
| GET | `/projects/:id/events?since=<seq>` | **SSE**：先回放 `events.jsonl` 中 seq > since 的事件，再实时推送；`event: agent`，`data: AgentEvent`；心跳 `: ping` 每 15s |
| GET | `/projects/:id/artifacts` | 文件树 `{ name, path, type: "dir"|"file", size, mtime, children? }`（排除 `.claude-home`、`toolkit`、`node_modules`、`.cache`） |
| GET | `/projects/:id/file?path=` | 原样下载/内联（按扩展名设置 Content-Type；路径必须在工作区内） |
| GET | `/projects/:id/preview?path=` | `{ kind: "markdown"|"image"|"slides"|"table"|"pdf"|"text"|"binary", content?, url?, sheets? }`；md 中相对图片改写为 `/api/projects/:id/file?path=` |
| GET | `/projects/:id/review` | `{ latest: ReviewRound, rounds: ReviewRound[], machine: MachineCheckReport }` |
| GET | `/projects/:id/numbers` | `numbers.json` 内容 + 一致性检查结果 |
| POST | `/projects/:id/package` | 重新打包 → `{ path }` |
| POST | `/projects/:id/open-folder` | 在系统文件管理器打开工作区（本地应用便利功能） |

## 4. 事件 `AgentEvent`（SSE 与 `events.jsonl` 每行一个）

```jsonc
{ "seq": 123, "ts": "2026-10-24T15:31:02.123Z", "projectId": "p_…",
  "type": "status" | "stage" | "text" | "text_delta" | "thinking" | "tool_use" | "tool_result" | "file" | "subagent" | "usage" | "review" | "error" | "done",
  "agent": { "id": "main" | "<parent_tool_use_id>", "name": "主代理" | "行业研究员" | "…" },   // 子代理用 Agent 工具的 tool_use_id 标识，名称来自其 input.subagent_type/description
  "data": { … } }
```

| type | data |
|---|---|
| `status` | `{ status: project.status, message }` |
| `stage` | `{ stage, status, note, round? }` |
| `text` | `{ text }`（完整助手文本块） |
| `text_delta` | `{ text }`（增量，前端拼接；服务端已按 16ms 合并） |
| `thinking` | `{ text, delta: bool }` |
| `tool_use` | `{ toolUseId, name, input }`（input 截断至 4KB；Write/Edit 只保留路径与前 300 字） |
| `tool_result` | `{ toolUseId, name, isError, output }`（截断 4KB） |
| `file` | `{ path, op: "write"|"edit"|"delete", size }`（来自 hooks / 文件监听） |
| `subagent` | `{ toolUseId, name, status: "start"|"end", summary? }` |
| `usage` | `{ inputTokens, outputTokens, cacheReadTokens, turns, costEstimateUsd }` |
| `review` | `{ round, score, byJudge, hardGatePassed }` |
| `error` | `{ message, recoverable: bool }` |
| `done` | `{ subtype: "success"|"error"|"stopped", result, usage }` |

## 5. 任务书 `任务书.md`（服务端渲染后放在工作区根，主代理第一步阅读）

```markdown
# 任务书
- 赛道：品牌策划赛道（brand） | 组别：本科组 | 语言：中文
- 作品名：…… | 团队名：……（匿名制：全文不得出现院校/姓名/联系方式）
- 策划对象：品牌 弘宝汝瓷（汝州市弘宝汝瓷坊有限公司，河南汝州，非遗工艺品/文创）
- 是否命题企业：否 | 命题题目：—
- 项目简介：……
- 团队（匿名）：队长/总策划（市场营销，品牌传播）；……
- 参数：调研模式=示例数据 | 正文目标 60 页（上限 80）| 讲稿 10 分钟 | 评审轮数 3 | 目标分 85
- 用户资料（已预读为文本，路径见下）：
  - inputs/企业资料.pdf → inputs/_preread/企业资料.pdf.txt（2.3 万字）
- 赛道规范：tracks/brand/spec.md | 量规：tracks/brand/rubric.json | 命名：tracks/brand/naming.json
- 开始指令：按 CLAUDE.md 执行流程完成全部阶段，直至 08_提交包 生成。
```

## 6. 进度协议

- 主代理在每个阶段开始/结束调用 MCP 工具 `mcp__bizgold__update_progress({ stage, status: "running"|"done"|"failed", note?, round? })`。
- 兜底：主代理每完成一步在对话中输出完整清单：

```
1. ✅ 读题解析
2. ✅ 调研
3. 🔄 策略
4. ⏳ 财务
5. ⏳ 图表
6. ⏳ 撰写
7. ⏳ 路演
8. ⏳ 评审（第 0/3 轮）
9. ⏳ 打包
```
服务端用正则 `^\s*(\d)\.\s*(✅|🔄|⏳|❌)\s*(读题解析|调研|策略|财务|图表|撰写|路演|评审|打包)` 解析。

## 7. 工具包数据文件

### 7.1 `04_图表/charts.json`

```jsonc
{ "outputDir": "04_图表", "dpi": 300, "theme": "bizgold",
  "charts": [
    { "id": "fig_market_size", "type": "stacked_bar", "title": "中国汝瓷工艺品市场规模（2019—2025）",
      "data": { "categories": ["2019","2020","2021","2022","2023","2024","2025E"],
                "series": [ { "name": "线下", "values": [12.1, 10.8, 13.2, 12.9, 15.4, 17.0, 18.6] },
                            { "name": "线上", "values": [3.2, 4.9, 6.8, 8.1, 10.2, 12.5, 14.9] } ] },
      "xLabel": "", "yLabel": "市场规模（亿元）", "unit": "亿元", "valueLabels": true,
      "source": "资料来源：中国陶瓷工业协会《2025 年陶瓷行业运行报告》；团队整理", "note": "2025E 为估计值",
      "width": 160, "height": 90 }  // 毫米；省略用默认 160×90
  ] }
```
`type` ∈ `bar | grouped_bar | stacked_bar | hbar | line | area | pie_as_bar | waterfall | funnel | tornado | breakeven | cashflow_combo | scatter_matrix | gantt | radar | bmc | positioning_map | donut`（`donut` 仅允许 ≤4 类）。每种 type 的 `data` 结构在 `biz-charts` 技能中给出示例。输出 `04_图表/<id>.png`，并生成 `04_图表/index.json`（id → 文件、标题、来源）。

### 7.2 `05_正文/book.json`

```jsonc
{ "track": "brand", "docType": "品牌策划书",
  "cover": { "title": "宋韵落天青，一釉生万象", "subtitle": "弘宝汝瓷品牌策划书", "teamName": "青瓷工作室",
             "company": "指导企业：汝州市弘宝汝瓷坊有限公司", "group": "本科组", "date": "2026 年 3 月", "extraLines": [] },
  "abstract": "05_正文/00_概要.md",            // 概要（1 页）
  "chapters": [ { "file": "05_正文/01_引言.md", "numbered": true }, { "file": "05_正文/02_品牌现状分析.md" } ],
  "appendices": [ { "file": "05_正文/A1_调研问卷.md" } ],
  "references": "01_调研/参考文献.json",       // 条目数组 → GB/T 7714 顺序编码
  "style": { "bodyFont": "宋体", "bodySize": 12, "headingFont": "黑体", "latinFont": "Times New Roman",
             "lineSpacing": 1.5, "firstLineIndent": 2, "pageNumbers": true, "toc": true, "tocDepth": 2,
             "header": "弘宝汝瓷品牌策划书", "watermark": "" },
  "anonymize": { "forbidden": ["XX大学", "XX学院"], "whitelist": ["河南财经", "弘宝"] },
  "output": "产出/弘宝汝瓷品牌策划书.docx", "pdf": true }
```

章节 Markdown 约定（渲染器支持的子集）：`#`/`##`/`###` 标题（渲染时按章节自动编号 1、1.1、1.1.1）；段落；`**粗体**`；有序/无序列表；表格（GFM，首行表头；单元格支持 `**`）；图片 `![图题](../04_图表/fig_x.png)` → 居中 + 自动“图 N-M 图题”；表格前一行 `表：表题` 或 `Table: 表题` → 自动“表 N-M 表题”；`> 引用`；`---` 分页；`[^n]` 脚注；行内 `{来源：…}` 渲染为小号灰色来源行；`<!-- pagebreak -->`。参考文献引用写 `[12]`。

参考文献 `01_调研/参考文献.json`：`[{ "id": 1, "type": "report|web|paper|book|annual-report|gov|news", "author": "…", "title": "…", "publisher": "…", "year": 2025, "url": "…", "accessed": "2026-10-24", "note": "" }]`。

### 7.3 `06_路演/slides.json`

```jsonc
{ "title": "弘宝汝瓷品牌策划路演", "subtitle": "青瓷工作室 · 本科组", "aspect": "16:9", "theme": "bizgold",
  "minutes": 10, "footer": "弘宝汝瓷品牌策划 · 青瓷工作室",
  "slides": [
    { "layout": "cover", "title": "宋韵落天青，一釉生万象", "subtitle": "弘宝汝瓷品牌策划", "note": "开场 20 秒：……" },
    { "layout": "agenda", "title": "目录", "items": ["现状诊断", "定位与核心价值", "品牌元素", "传播执行", "预算与效果", "风险与保护"] },
    { "layout": "section", "title": "一、现状诊断", "kicker": "问题在哪" },
    { "layout": "bullets", "title": "结论句标题：汝瓷市场 6 年增长 2.3 倍，但弘宝线上份额不足 3%",
      "bullets": ["要点一（≤30 字）", "要点二", "要点三"], "highlight": "18.6 亿元", "highlightLabel": "2025E 市场规模", "note": "讲稿备注" },
    { "layout": "image", "title": "……", "image": "../04_图表/fig_market_size.png", "caption": "资料来源：……", "bullets": ["读图一句话"] },
    { "layout": "two_column", "title": "……", "left": { "heading": "现状", "bullets": ["…"] }, "right": { "heading": "目标", "bullets": ["…"] } },
    { "layout": "chart", "title": "……", "chart": { "type": "bar"|"line"|"pie"|"column_clustered", "categories": ["2026","2027","2028"], "series": [{ "name": "收入", "values": [1.2, 2.8, 4.9] }], "unit": "百万元" }, "insight": "一句结论" },
    { "layout": "table", "title": "……", "columns": ["维度", "竞品 A", "竞品 B", "弘宝"], "rows": [["价格带", "…", "…", "…"]], "note": "" },
    { "layout": "kpi", "title": "……", "kpis": [ { "value": "384", "label": "有效问卷" }, { "value": "72%", "label": "愿意尝试新品" } ] },
    { "layout": "timeline", "title": "……", "steps": [ { "when": "2026 Q2", "what": "…" } ] },
    { "layout": "team", "title": "团队与分工", "members": [ { "role": "总策划", "duty": "…" } ] },
    { "layout": "quote", "title": "品牌主张", "quote": "一釉生万象", "by": "Big Idea" },
    { "layout": "closing", "title": "谢谢", "subtitle": "恳请评委指正", "contact": "" },
    { "layout": "appendix", "title": "备用页：敏感性分析", "image": "../04_图表/fig_tornado.png" }
  ],
  "output": "产出/路演.pptx" }
```
每页可带 `note`（演讲者备注）与 `duration`（秒）。渲染器限制：每页正文 ≤120 字、要点 ≤6 条、正文字号 ≥ 18pt（要点页 20–24）、标题 ≥ 28pt。

### 7.4 `03_财务/assumptions.json`（创业计划 / 需要财务模型的赛道）

```jsonc
{ "mode": "startup",                       // startup（三表模型）| budget（策划案预算表）
  "company": "弘宝汝瓷文创", "currency": "元", "unit": 10000,   // 展示单位：万元
  "startYear": 2026, "years": 5, "monthlyFirstYear": true,
  "scenarios": { "base": 1.0, "optimistic": 1.2, "pessimistic": 0.8 },   // 需求乘数
  "revenueStreams": [
    { "name": "文创礼品线上零售", "driver": "volume_price",
      "volume": { "y1Monthly": [200,220,260,300,340,380,420,460,520,600,700,800], "yoyGrowth": [0, 0.8, 0.5, 0.35, 0.25] },
      "price": { "base": 268, "yoyChange": [0, 0.03, 0.03, 0.03, 0.03] },
      "cogsRatio": 0.42, "source": "单价来自天猫同类 TOP10 均价；首年销量来自渠道访谈" },
    { "name": "企业定制", "driver": "clients_arpu", "clients": { "y1": 12, "yoyGrowth": [0, 1.0, 0.6, 0.4, 0.3] }, "arpu": { "base": 50000 }, "cogsRatio": 0.55, "source": "…" }
  ],
  "costs": {
    "fixed": [ { "name": "场地租金", "monthly": 15000, "yoyChange": 0.05 }, { "name": "SaaS 与软件", "monthly": 3000 } ],
    "variablePctOfRevenue": [ { "name": "平台佣金与支付", "pct": 0.06 }, { "name": "物流", "pct": 0.04 } ],
    "marketing": { "pctOfRevenue": [0.18, 0.15, 0.12, 0.10, 0.09], "minAnnual": 300000 },
    "headcount": [ { "role": "品牌/运营", "y": [3,5,7,9,11], "monthlySalary": 9000, "socialRate": 0.30 },
                   { "role": "设计/产品", "y": [2,3,4,5,6], "monthlySalary": 11000, "socialRate": 0.30 } ],
    "capex": [ { "name": "设备与装修", "year": 1, "amount": 400000, "depreciationYears": 5 } ],
    "otherOpexPctOfRevenue": 0.03 },
  "workingCapital": { "receivableDays": 30, "inventoryDays": 60, "payableDays": 45 },
  "tax": { "incomeTaxRate": 0.25, "vatIncluded": false },
  "financing": [ { "round": "种子轮", "year": 1, "month": 1, "amount": 3000000, "preMoney": 12000000, "use": [ { "item": "产品研发", "pct": 0.35 }, { "item": "品牌营销", "pct": 0.40 }, { "item": "团队与运营", "pct": 0.25 } ] } ],
  "openingCash": 500000,
  "valuation": { "methods": ["dcf", "ps_comparables", "vc_method"], "wacc": 0.18, "terminalGrowth": 0.03, "psMultiple": 3.0, "exitPE": 20, "exitYear": 5, "targetIrr": 0.40 },
  "unitEconomics": { "cac": 85, "repeatRate": 0.35, "grossMarginRef": 0.58, "lifetimeOrders": 3.2 },
  "sensitivity": { "drivers": ["price", "volume", "cogsRatio", "cac", "fixedCost"], "range": [-0.3, -0.15, 0.15, 0.3] }
}
```
`budget` 模式：
```jsonc
{ "mode": "budget", "currency": "元", "totalCap": 3000000,
  "phases": [ { "name": "第一阶段：声量引爆（2026.4–6）", "items": [ { "item": "抖音信息流投放", "unitPrice": 35, "quantity": 20000, "unit": "CPM", "source": "巨量引擎 2025 刊例" } ] } ],
  "contingencyPct": 0.08,
  "kpis": [ { "metric": "曝光量", "baseline": "—", "target": "5000 万", "tool": "巨量云图", "cycle": "周" } ] }
```

### 7.5 `03_财务/numbers.json`（数字账本；工具包生成，MCP `record_numbers` 可追加）

```jsonc
{ "generatedAt": "…", "unit": "万元",
  "entries": [
    { "key": "revenue_y1", "label": "第 1 年营业收入", "value": 286.4, "unit": "万元", "source": "03_财务/财务模型.xlsx!利润表!C5", "usedIn": [] },
    { "key": "revenue_y3", "label": "第 3 年营业收入", "value": 1180.2, "unit": "万元", "source": "…" },
    { "key": "gross_margin_y3", "label": "第 3 年毛利率", "value": 0.57, "unit": "%", "source": "…" },
    { "key": "breakeven_month", "label": "盈亏平衡月份", "value": 19, "unit": "月", "source": "…" },
    { "key": "min_cash", "label": "现金最低点", "value": -212.5, "unit": "万元", "source": "…" },
    { "key": "financing_amount", "label": "融资额", "value": 300, "unit": "万元", "source": "…" },
    { "key": "post_money", "label": "投后估值", "value": 1500, "unit": "万元" },
    { "key": "equity_offered", "label": "出让比例", "value": 0.20, "unit": "%" },
    { "key": "som_y1", "label": "首年可获得市场 SOM", "value": 420, "unit": "万元", "source": "01_调研/数据/市场规模.csv" },
    { "key": "sample_size", "label": "问卷有效样本量", "value": 412, "unit": "份", "source": "01_调研/调研工具包/示例结果.md", "synthetic": true },
    { "key": "budget_total", "label": "推广预算总额", "value": 300, "unit": "万元" }
  ] }
```
固定 key 集合（质检比对用）：`revenue_y1..y5, net_profit_y1..y5, gross_margin_y1..y5, breakeven_month, min_cash, financing_amount, pre_money, post_money, equity_offered, cac, ltv, ltv_cac, payback_months, tam, sam, som_y1, sample_size, budget_total, kpi_*`。

### 7.6 评审 `07_评审/第N轮/评审_<judge>.json`

```jsonc
{ "round": 1, "judge": "investor", "model": "kimi-k3", "track": "brand", "rubric": "tracks/brand/rubric.json",
  "scores": [ { "dimension": "现状分析", "weight": 16, "score": 4, "cap": 5, "evidence": "…引用正文位置…", "issues": ["…"], "fixes": ["…具体改法…"] } ],
  "total": 81.2,                            // Σ score/5*weight
  "officialEstimate": { "written": 49, "pitch": 33, "total": 82 },   // 按官方评分表口径的模拟分
  "hardGates": [ { "gate": "anonymity", "passed": true, "detail": "" } ],
  "topFixes": [ { "priority": 1, "where": "02_品牌现状分析.md §2.3", "what": "补 384 份问卷的样本结构表", "gain": "+3" } ],
  "verdict": "接近一等奖水平；修复前三项后可达 85+" }
```
汇总 `07_评审/第N轮/评审汇总.json`：`{ round, score（三评委中位数，缺席 1 人时取均值）, mean, byJudge, byDimension, minDimension, minJudgeDimension?（任一评委任一维度最低分，可选）, flags?（三评委 flags 并集，可选，如 synthetic_unmarked / anonymity / brief_not_answered）, machine: {…}, hardGatePassed, officialEstimate: {written, pitch, total}, fixes: [{ id, owner: "writer|pitch|finance|main", file, where, what, accept, gain }], converged: bool, remaining: string[] }`。`rubric` 字段固定为 `"tracks/<赛道>/rubric.json"`。

### 7.7 机器检查 `07_评审/第N轮/机器检查.json`

```jsonc
{ "track": "brand", "generatedAt": "…",
  "checks": [ { "id": "struct_coverage", "value": 1.0, "threshold": { "pass": 1.0, "warn": 0.85 }, "status": "pass"|"warn"|"fail", "detail": "10/10 要素命中" } ],
  "hardGates": [ { "id": "anonymity_hits", "passed": true, "hits": [] }, { "id": "page_limit", "passed": true, "pages": 63, "limit": 80 } ],
  "crossDoc": { "consistent": true, "mismatches": [ { "key": "financing_amount", "book": 300, "deck": 300, "script": 280 } ] },
  "summary": { "pass": 18, "warn": 3, "fail": 1 } }
```
检查项 id 与阈值见 `02_国金级材料标准.md §6.7`，赛道阈值在 `tracks/<id>/checks.json`。

2026.10 新增三项文风检查（只追加到 `checks[]`，不进 `hardGates`；单章快检 `--chapter` 也输出）：

| id | 定义 | value / threshold | 状态规则 | locations 条目 |
|---|---|---|---|---|
| `cliche_terms` | 套话 / 万能模板句密度：正文（不含附录）+ 讲稿，词表 `checks.json.clicheTerms`（字面词、含 `…` 的模板句如「随着…发展」、`{term, regex:true}`），`clicheWhitelist` 先遮蔽 | 每万字命中次数 / `{pass: 3, fail: 9, per: "万字"}`（`clichePer10k` 可改） | ≤ 3 pass；≥ 9（3 倍）fail；其间 warn | `{file, line, term, text}`；另有 `hits`、`words` |
| `synthetic_research_markers` | 示例数据标注（章节级）：`params.researchMode == "synthetic"` 或账本含 `synthetic:true` 时，正文 / 讲稿 / 题库中引用一手调研结论的句子（调研显示 / 受访者… / 问卷·访谈·样本 + 计数或百分比 / 账本示例数值带同类单位）所在**文件**、PPT **页**内必须出现「示例数据」标记 | 缺标记的文件 / 页数 / `{pass: 0}` | 0 pass；≥ 1 fail；非 synthetic 模式 skipped | 触发句 `{file, line, term, text}`；另有 `files`、`missingFiles[]`。与既有 `synthetic_marker_coverage`（账本示例数值 ±3 行内须有标记）互补 |
| `unsourced_numbers` | 无源关键数字：正文散文句中带 亿元 / 万元 / 万 / 亿 / 元 / % / 个百分点 / 倍 / 万人 等单位的数字，同句或同段须有 `[n]`、`[^n]`、「示例数据」、来源 / 据《 / 据…局·协会·报告、测算 / 模型 / 假设 / 估算、命题、一手调研口径（问卷 / 访谈 / 受访者）、「见第 N 章 / 见表」指向，或数值能在 `numbers.json` 找到（允许千分位 / 小数位 / 万↔亿↔元 / 20%↔0.2）；读图 / 读表句与紧跟图表之后的段落由图表来源覆盖；年份、章节号、图表号、页码、文献号本身、表格单元格、区间后段（25—35）不计 | 无源 / 关键数字总数 / `{pass: 0.10, warn: 0.25}` | ≤ 10% pass；≤ 25% warn；否则 fail | 样例 ≤ 20 条 `{file, line, term, text(句子前 60 字)}`；另有 `unsourced`、`total` |

### 7.8 赛道文件

- `rubric.json`：`{ "track", "dimensions": [ { "id", "name", "weight", "official": "品牌现状分析准确性 10", "anchors": { "1": "…", "3": "…", "5": "…" }, "bandAnchors"?: { "60": "…", "75": "…", "85": "…", "92": "…" }, "evidenceRequired"?: ["…"], "machineCapRules": [ { "check": "has_primary_research", "when": "fail", "cap": 2 } ] } ], "targetScore": 85, "minDimension": 3, "bands"?: {…}, "scoringRules"?: [...], "judgeFocus"?: {…} }`（带 ? 的为 2026.2 起的可选键，只供评委阅读；`checks.py` 只读 `hardGates`）
- `checks.json`：`{ "pageLimit": 80, "wordMin": 0, "chartMin": 20, "sourceMin": 15, "numericDensity": 15, "requiredElements": [ { "id": "summary", "name": "概要", "aliases": ["执行摘要", "摘要"] } ], "forbiddenTerms": ["…"], "anonymityRegex": "…" }`
- `naming.json`：`{ "fileNameRule": "{抽签编号或学校}-{团队名}-{作品名}", "coverFields": ["作品名称", "团队名称", "指导企业名称"], "coverForbidden": ["院校", "姓名", "联系方式"], "declarations": ["作品原创性声明", "作品版权声明", "企业授权书"], "deckMinPages": 0, "videoMinutes": 10, "qaMinutes": 5 }`
- `outline.json`：默认章节树 `[ { "no": "1", "title": "引言", "file": "01_引言.md", "pages": [1,2], "mustAnswer": ["…"], "charts": ["…"], "scoring": ["内容完整程度"] } ]`

## 8. 子代理任务书格式（主代理 → 子代理，写入 `任务/<role>_<n>.md` 后在 Agent 调用中引用）

```markdown
# 子任务：行业研究员 #2 —— 竞品与渠道
- 背景：（两三句，含品牌/项目名与赛道）
- 目标：回答以下问题……（编号）
- 输入：任务书.md、00_任务解析/任务解析.md、inputs/_preread/*
- 输出：01_调研/02_竞品与渠道.md（结构：结论摘要 → 证据表（数值+口径+年份+来源 URL）→ 对策划的含义）；数据 CSV 写到 01_调研/数据/；新增参考文献追加到 01_调研/参考文献.json
- 约束：每条数据带来源；无法核实的标注【待核实】；不得编造 URL；检索次数 ≤ 12；完成后返回 ≤300 字摘要
```

## 9. 环境变量（服务端读取）

`BIZGOLD_PORT`（默认 5178）、`BIZGOLD_API_BASE_URL`、`BIZGOLD_API_KEY`（或 `ANTHROPIC_AUTH_TOKEN`）、`BIZGOLD_MODEL`、`BIZGOLD_WORKSPACE_DIR`、`BIZGOLD_PYTHON`、`BIZGOLD_SETTINGS_DIR`（默认 `~/.bizgold`）、`BOCHA_API_KEY` / `TAVILY_API_KEY` / `SERPER_API_KEY`。环境变量优先级高于设置文件。
