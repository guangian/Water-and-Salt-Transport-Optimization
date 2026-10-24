# 金策 BizGold —— 全国高校商业精英挑战赛 国金级智能体 · 架构设计

> 版本 v1.0（2026-10-24）。本文件是整个项目的总设计与各模块开发者的共同依据；接口字段以 `docs/contracts.md` 为准。

## 0. 一句话定位

**把一支顶级商赛团队（行业研究员 + 财务建模师 + 主笔 + 路演设计师 + 三位评委 + 质检员）装进一个本地 Web 应用**：用户给出赛道、命题/品牌/项目信息与资料，系统端到端产出 **国赛金奖级** 的全套参赛材料——策划书/商业计划书（Word + PDF）、财务模型（Excel，带公式）、路演 PPT、10 分钟讲稿、答辩题库、调研工具包、合规报告与提交包——并用官方评分表反复自评、迭代到目标分。

设计原则（来自对参考项目 Mrite 的拆解与两份调研报告）：

1. **成熟内核 + 领域规范**：引擎用 Claude Agent SDK（Claude Code 同源能力：Read/Write/Edit/Bash/Grep/Glob/Agent/Skill），领域知识全部写进 `CLAUDE.md`、赛道规范、子代理与技能文件——模型负责内容，规范负责水准。
2. **内容与排版分离**：模型只写结构化内容（Markdown 章节、JSON 幻灯片、JSON 财务假设），确定性的 Python 工具包负责渲染 Word / PPT / Excel / 图表，保证任何模型、任何一次运行都得到专业版式。
3. **一本数字账本**：财务模型先产出 `numbers.json`，正文、PPT、讲稿、题库全部引用同一账本；质检员机器比对跨材料一致性（评委最常抓的穿帮点）。
4. **评分表驱动**：章节大纲与官方评分项一一对应；生成后由三位不同模型的评委按量规打分，分数低于目标线自动进入修改轮。
5. **合规零扣分**：匿名、页数、字体、命名、禁语、AI 使用规定全部做成机器检查项与硬门槛。
6. **诚实的数据**：真实的一手调研只能由学生完成。系统生成问卷/访谈提纲/执行计划，并（在"示例数据"模式下）生成**明确标注为示例的模拟结果**供占位，所有示例数据在产出物中标记“【示例数据，须替换为真实调研结果】”，并在《数据真实性清单》中逐条列出。

## 1. 支持的赛道（v1）

| 赛道 id | 名称 | 书面材料 | 官方评分映射 | v1 完成度 |
|---|---|---|---|---|
| `brand` | 品牌策划赛道（含在华留学生赛道） | 品牌策划书 10 要素，≤80 页 | 策划书 60（6×10）+ 展示答辩 40（8×5） | 完整 |
| `startup` | 创业计划赛道（命题组/创业组/创客组） | 商业计划书 ≥3000 字 ≤80 页 宋体小四 | 创新 20 / 商业 30 / 团队 10 + 陈述 20 / 答辩 20 | 完整 |
| `digital` | 数字营销赛道 | 数字营销方案 9 要素 + 短视频脚本 | 方案 60（6×10）+ 展示答辩 40 | 完整 |
| `case` | 流通业经营模拟 / 会计与商业管理案例 | 案例分析报告（8 / 12 要素） | 书面 60 或 70 + 展示 | 规范版 |
| `mice` | 商务会奖旅游策划 | 会奖策划书 10 要素，预算 200–400 万 | 书面 60 + 演讲答辩 40 | 规范版 |

每个赛道一个目录 `agent/tracks/<id>/`：`spec.md`（章节大纲 + 每章必答 + 篇幅 + 评分映射 + 文风）、`rubric.json`（评委量规）、`checks.json`（机器检查阈值）、`naming.json`（命名/封面/匿名规则）、`outline.json`（默认章节树）。新增赛道只需新增目录。

## 2. 总体架构

```
┌───────────────────────────── 浏览器（React + Vite + Tailwind）──────────────────────────────┐
│ 工作台 │ 新建项目向导 │ 项目详情：运行控制台 · 产出物 · 评审 · 修改对话 · 数字账本 │ 设置 │ 知识库 │
└───────────────────────────────▲──────────────────────────────────────────────────────────┘
                                │ REST + SSE（/api/*）
┌───────────────────────────────┴──────────────────────────────────────────────────────────┐
│ apps/server（Node 22 + TypeScript + Fastify）                                             │
│  ProjectStore（workspace/<id>/project.json, events.jsonl）  SettingsStore（~/.bizgold）    │
│  AgentRunner ── Claude Agent SDK query() ── claude 二进制 ──► 第三方 Anthropic 兼容网关     │
│     ├─ 进程内 MCP 工具：web_search / web_fetch / update_progress / record_numbers           │
│     ├─ hooks：PostToolUse 记录文件写入；Stop 兜底                                            │
│     └─ 事件归一化 → SSE 推送 + 落盘                                                        │
│  PythonRunner（toolkit CLI）  PreviewService（md/图片/pptx→html/docx→pdf）  Packager（zip）  │
└───────────────────────────────▲──────────────────────────────────────────────────────────┘
                                │ 复制到每个项目工作区
┌───────────────────────────────┴──────────────────────────────────────────────────────────┐
│ agent/（智能体大脑，随工作区分发）                                                           │
│  CLAUDE.md（总流程规范） tracks/<id>/（赛道规范+量规） .claude/agents/（子代理）              │
│  .claude/skills/（技能） toolkit/bizkit（Python 渲染与检查工具包） knowledge/（知识库）        │
└──────────────────────────────────────────────────────────────────────────────────────────┘
```

### 2.1 目录结构

```
bizgold/
├── package.json                # pnpm workspace 根；scripts: dev / build / start / check
├── pnpm-workspace.yaml
├── apps/
│   ├── server/                 # Fastify 服务 + AgentRunner
│   │   ├── src/index.ts        # 启动：静态托管 web/dist，注册路由
│   │   ├── src/routes/*.ts     # settings / projects / events / artifacts / review / tracks / knowledge
│   │   ├── src/agent/runner.ts # 封装 query()，事件归一化，进度协议，会话恢复
│   │   ├── src/agent/mcp.ts    # 进程内 MCP 工具（web_search / web_fetch / update_progress / record_numbers）
│   │   ├── src/agent/prompt.ts # 任务书与启动提示词构建
│   │   ├── src/services/*.ts   # store / python / preview / packager / search / fonts
│   │   └── package.json
│   └── web/                    # React 前端
│       ├── src/pages/*.tsx     # Dashboard / NewProject / Project / Settings / Knowledge
│       ├── src/components/*    # 控制台、进度、文件树、预览器、评分雷达、数字账本……
│       ├── src/api.ts          # fetch + SSE 客户端
│       └── package.json
├── agent/                      # 智能体大脑（见第 4 章）
│   ├── CLAUDE.md
│   ├── tracks/{brand,startup,digital,case,mice}/
│   ├── .claude/agents/*.md
│   ├── .claude/skills/*/SKILL.md
│   ├── .claude/settings.json   # hooks（可选）、permissions
│   ├── toolkit/bizkit/*.py     # Python 工具包 + requirements.txt
│   └── knowledge/*.md          # 精炼知识（评分表、要素清单、数据源、写作范例、评委 50 问……）
├── docs/                       # 调研报告、架构、契约、用户手册、测试记录
├── scripts/                    # setup.ps1 / setup.sh / check-env.* / run-e2e.*
└── workspace/                  # 运行期项目（.gitignore）
```

### 2.2 一次项目运行的流程（9 个阶段）

| # | 阶段 id | 名称 | 执行者 | 产出（工作区相对路径） |
|---|---|---|---|---|
| 1 | `parse` | 读题解析 | 主代理 | `00_任务解析/任务解析.md`（命题解读、评分映射、材料清单、约束与假设） |
| 2 | `research` | 调研 | 行业研究员 ×N（并行）+ 消费者洞察 | `01_调研/*.md`（每条数据带来源）、`01_调研/数据/*.csv`、`01_调研/参考文献.json`、`01_调研/调研工具包/`（问卷、访谈提纲、执行计划、示例结果） |
| 3 | `strategy` | 策略 | 主代理 | `02_策略/策略框架.md`（定位逻辑链 / 商业模式画布 / Big Idea / 章节结论句） |
| 4 | `finance` | 财务 | 财务建模师 | `03_财务/assumptions.json` → `财务模型.xlsx` + `numbers.json` + `财务摘要.md` + `建模说明.md` |
| 5 | `charts` | 图表 | 主代理（工具包） | `04_图表/*.png` + `04_图表/charts.json` |
| 6 | `write` | 撰写 | 主笔 | `05_正文/NN_*.md` → `05_正文/book.json` → `产出/<作品名>.docx`（+pdf） |
| 7 | `pitch` | 路演 | 路演设计师 | `06_路演/slides.json` → `产出/路演.pptx`；`讲稿.md`；`答辩题库.md`；`视频分镜.md` |
| 8 | `review` | 评审 | 评委团 ×3 + 质检员 | `07_评审/第N轮/评审_*.json|md`、`机器检查.json`、`质检整改清单.md`、`评审汇总.json`、`修改清单.md`；未达标按 owner 派修改任务（≤ 最大轮数） |
| 9 | `package` | 打包 | 主代理（工具包） | `08_提交包/`（按赛道命名规则重命名、声明/授权书模板、数据真实性清单、README）+ zip |

进度由主代理通过 MCP 工具 `update_progress(stage, status, note)` 上报（服务端另有基于“✅/🔄/⏳ 列表”文本的兜底解析）。

### 2.3 模型分工（默认，可在设置中改）

| 角色 | 默认模型 | 理由 |
|---|---|---|
| 主代理 / 主笔 / 路演 | `deepseek-v4-pro` | 长中文文档、指令遵循稳定 |
| 行业研究员 / 消费者洞察 | `kimi-k3` | 检索归纳与长上下文 |
| 财务建模师 | `deepseek-v4-pro` | 数值与 JSON 严谨 |
| 评委 A（投资人）| `kimi-k3` | 多模型评委团，降低“自评自夸” |
| 评委 B（学者）| `glm-5.3` | 同上 |
| 评委 C（企业高管）| `qwen3.8-max` | 同上 |
| 质检员 | `deepseek-v4-pro` | 机器检查解释 |

子代理模型通过 `.claude/agents/*.md` 的 frontmatter `model:` 指定；服务端在启动前把设置里的模型名写入工作区的 agents 文件（模板占位符替换）。

## 3. 服务端设计（apps/server）

### 3.1 关键组件

- **SettingsStore**：`~/.bizgold/settings.json`（Windows：`%USERPROFILE%\.bizgold\settings.json`），字段见契约。API 密钥只存在该文件或环境变量 `BIZGOLD_API_KEY` / `ANTHROPIC_AUTH_TOKEN`，**绝不写入仓库目录**。
- **ProjectStore**：`workspace/<projectId>/project.json`（元数据、输入、阶段状态、评分、会话 id）、`events.jsonl`（归一化事件，供 SSE 回放）、`inputs/`（用户上传）、以及 2.2 的产出目录。
- **AgentRunner**：
  - 准备工作区：复制 `agent/`（CLAUDE.md、.claude、toolkit、knowledge、tracks/<track>）；渲染 `任务书.md`（赛道、组别、命题/品牌/项目信息、团队信息、参数、上传资料预读摘要）；渲染 agents 模板中的模型占位符。
  - 预读输入：PDF→文本（pdfplumber，经 PythonRunner）、docx→文本、xlsx/csv→列摘要、图片→保留路径并提示用 Read 查看。
  - 启动 `query()`：`systemPrompt: {type:'preset', preset:'claude_code', append}`、`settingSources:['project']`、`permissionMode:'bypassPermissions'`、`allowDangerouslySkipPermissions:true`、`maxTurns` 可配置（默认 1500）、`includePartialMessages:true`、`persistSession:true`、`mcpServers:{bizgold}`、`disallowedTools:['WebSearch','WebFetch','AskUserQuestion']`（内置 WebFetch 绕过来源记录，统一走 MCP `web_fetch`）、`env` 见技术手册 §1.11/§4（`...process.env`、`ANTHROPIC_BASE_URL`、`ANTHROPIC_AUTH_TOKEN`、`ANTHROPIC_MODEL`、`ANTHROPIC_DEFAULT_HAIKU_MODEL=<同主模型或快速模型>`、`CLAUDE_CODE_SUBAGENT_MODEL`、`CLAUDE_CONFIG_DIR=<workspace>/.claude-home`、关遥测/自动更新/自动记忆、`CLAUDE_CODE_MAX_CONTEXT_TOKENS`、Windows 下 `CLAUDE_CODE_GIT_BASH_PATH` 自动探测）、`stderr` 收集。
  - 事件归一化：SDK 消息 → `AgentEvent`（契约 §4），含子代理归属（`parent_tool_use_id`）；增量文本 16ms 合并；落盘 + 广播。
  - 超时与保护：总时长（默认 6h）、无输出超时（默认 45min）、工具超时；`abortController` 停止；进程清理。
  - 修改模式：`resume: sessionId` + 新提示词；禁止修改 `inputs/`。
  - 失败恢复：`project.json` 记录 `sessionId` 与阶段状态；“续跑”从当前阶段继续（提示词带已完成清单）。
- **MCP 工具**（`createSdkMcpServer` + `tool()`）：
  - `web_search(query, count?, site?, freshness?)`：按设置的提供商顺序尝试（博查 Bocha API → Tavily → Serper → Bing 网页解析 → DuckDuckGo lite），返回 `[{title,url,snippet,date?}]`；请求节流与缓存（内存 + `workspace/<id>/.cache/search.json`）。
  - `web_fetch(url, max_chars?)`：抓取并抽正文为 Markdown（去导航/脚本，保留表格），返回 `{title,url,content,fetched_at}`；拒绝内网地址；PDF 链接自动抽文本。
  - `update_progress(stage, status, note?)`：写 `project.json.stages` 并广播。
  - `record_numbers(entries)`：合并写入 `03_财务/numbers.json`（质检比对用）。
- **PythonRunner**：定位 Python（设置 > `python` > `python3`），首次运行 `pip install -r toolkit/requirements.txt`（国内镜像可配置），执行工具包 CLI，统一 UTF-8。
- **PreviewService**：Markdown（含图片重写为 /api/file 链接）、图片、`slides.json`→前端 HTML 渲染、`xlsx`→JSON（SheetJS）、docx→pdf（Windows 优先 Word COM / WPS，其次 LibreOffice；均无则仅提供 docx）、PDF 页数（pypdf）。
- **Packager**：按 `naming.json` 重命名、复制声明模板、生成 README 与数据真实性清单、打 zip。

### 3.2 安全与隔离

- 文件路由只允许访问 `workspace/<id>/` 内路径（规范化后前缀校验）。
- 工作区 `CLAUDE_CONFIG_DIR` 独立，不读用户 `~/.claude`。
- 不记录密钥到日志；`events.jsonl` 过滤 `Authorization`/`sk-` 模式。
- 搜索/抓取工具禁止访问 localhost/内网段。

## 4. 智能体大脑（agent/）

### 4.1 `CLAUDE.md` 总规范的结构（对标 Mrite 华为杯模板的颗粒度）

1. ⛔ 最高优先级规则：全自动不提问；Python 脚本 UTF-8 头；进度协议；数字溯源铁律（每个数字来自 `numbers.json`/调研数据文件/带来源的检索结果）；匿名铁律；示例数据标注铁律；中文命名；不读二进制；同一操作最多重试 3 次；不得结束回复直到提交包生成。
2. 工作区结构与文件命名。
3. 执行流程 Step 0–8（每步：输入 → 动作 → 工具命令 → 验收标准 → 产出）。
4. 子代理调度规范（何时并行、给子代理的任务书格式、结果回收与核验）。
5. 写作规范（国金级文风：结论先行、数值四要素、禁语表、句长、加粗克制、图注四段、章末小结）。
6. 图表与版式规范（工具包调色板、图型选择、每章图表配额）。
7. 评审与迭代协议（量规、目标线、修改清单格式、最大轮数、收敛判定）。
8. 合规与提交（页数/字体/命名/封面/声明/AI 规定）。
9. 常见错误速修表。

赛道差异全部放在 `tracks/<id>/spec.md`，CLAUDE.md 通过“读取 `tracks/当前赛道/spec.md`”引用。

### 4.2 子代理（`.claude/agents/*.md`）

| 文件 | 角色 | 工具 | 输出 |
|---|---|---|---|
| `industry-researcher.md` | 行业研究员 | web_search / web_fetch / Read / Write / Bash | 带来源的调研笔记 + CSV 数据 + 参考文献条目 |
| `consumer-insight.md` | 消费者洞察与调研设计 | 同上 | 用户画像、问卷、访谈提纲、调研计划、示例结果（标注） |
| `financial-modeler.md` | 财务建模师 | Read / Write / Bash | assumptions.json → 模型 xlsx + numbers.json + 财务摘要 |
| `chief-writer.md` | 主笔 | Read / Write / Edit / Bash | 各章 Markdown（按 spec 大纲，引用图表与账本） |
| `pitch-designer.md` | 路演设计师 | Read / Write / Bash | slides.json、讲稿、题库、分镜 |
| `judge-investor.md` / `judge-professor.md` / `judge-executive.md` | 评委团 | Read / Grep / Glob（只读） | rubric 打分 JSON + 修改意见（主代理落盘） |
| `chart-maker.md` | 制图员 | Read / Write / Edit / Bash | charts.json + 渲染 + 抽查 |
| `qa-inspector.md` | 质检员 | Read / Bash | 运行机器检查、跨材料一致性、合规报告 |

### 4.3 技能（`.claude/skills/*/SKILL.md`）

`docx-report`（book.json 规范与渲染命令）、`pptx-deck`（slides.json 规范、页型、设计规则）、`fin-model`（assumptions.json 规范、三表逻辑、估值）、`biz-charts`（charts.json 规范、选图）、`market-sizing`（TAM/SAM/SOM 双法）、`survey-design`（问卷/样本量/质控）、`brand-creative`（Big Idea、品牌元素、ICC 准则）、`pitch-qa`（讲稿节奏、答辩 50 问）、`compliance-check`（检查命令与修复）、`web-research`（检索策略、来源分级、引用格式 GB/T 7714）。

### 4.4 工具包 `toolkit/bizkit`（Python 3.10+）

| 模块 | 功能 | CLI |
|---|---|---|
| `fonts.py` | 中文字体发现（Windows 微软雅黑/宋体/黑体；Linux Noto/文泉驿）与 matplotlib/docx/pptx 字体名映射 | — |
| `theme.py` | 配色（主色海军蓝 `#14305C`、金 `#C9A227`、辅助灰蓝；图表 6 色序列；语义色）与字号体系 | — |
| `charts.py` | 柱/堆叠柱/分组柱/折线/面积/排序条形/瀑布/漏斗/龙卷风/盈亏平衡/现金流柱线/竞争矩阵散点/甘特/雷达/商业模式画布/定位图；300 dpi；统一 `来源：` 脚注 | `python -m bizkit.charts charts.json` |
| `docx_render.py` | `book.json`（封面、目录、章节 Markdown 文件列表、附录、参考文献）→ docx：样式（宋体小四/黑体标题/Times New Roman 西文）、多级编号、自动目录域、页眉页脚页码、图注表注 SEQ 域、三线表/底纹表、分页、引用；匿名清洗 | `python -m bizkit.docx_render book.json` |
| `pptx_render.py` | `slides.json` → 16:9 pptx（主题、14 种页型、原生图表、图片、备注、页码、匿名） | `python -m bizkit.pptx_render slides.json` |
| `finmodel.py` | `assumptions.json` → `财务模型.xlsx`（假设/驱动/收入/成本/人力/利润表/现金流量表/资产负债表/盈亏平衡/单位经济/敏感性/估值/检查 页；公式 + 同步写入计算值）+ `numbers.json`；`budget` 模式生成策划案预算表 | `python -m bizkit.finmodel assumptions.json` |
| `checks.py` | 契约 §6 的机器检查项（结构覆盖、匿名、数值密度、来源数、图表数、禁语、跨材料一致性、页数/字数估算、PPT 指标、讲稿时长、题库覆盖）→ `机器检查.json` + Markdown | `python -m bizkit.checks --track brand --workspace .` |
| `pdf.py` | docx→pdf（Word COM / WPS / LibreOffice）、页数统计 | `python -m bizkit.pdf in.docx` |
| `package.py` | 命名、声明模板、README、zip | `python -m bizkit.package --track brand` |
| `preread.py` | PDF/docx/xlsx/csv 预读为文本摘要 | `python -m bizkit.preread <file>` |
| `search.py` / `fetch.py` | 与服务端同逻辑的命令行检索/抓取（供脚本与离线调试） | `python -m bizkit.search "关键词"` |

## 5. 前端设计（apps/web）

- 视觉：深海军蓝 + 金色点缀的专业风格，浅色为默认、支持深色；Inter/系统中文字体；lucide 图标；recharts 雷达/条形。
- 页面：
  1. **工作台**：项目卡片（赛道徽标、状态、最近评分、更新时间）、快捷新建、运行中项目实时摘要、环境健康（API/Python/字体/PDF 转换器）。
  2. **新建项目向导**（三步）：① 赛道与组别（卡片选择，显示评分权重与材料要求摘要）② 项目信息（作品名/团队名、命题企业下拉或自填品牌/项目、项目简介、上传资料拖拽区、团队成员匿名信息）③ 生成参数（模型、调研模式、页数目标、讲稿时长、评审轮数、目标分、搜索提供商）→ 开始运行。
  3. **项目详情**：
     - 运行控制台：9 阶段进度（✅/🔄/⏳ + 用时）、实时流（主代理与各子代理分泳道，工具调用折叠卡片，思考折叠）、Token/费用估算、停止/续跑；
     - 产出物：左树右预览（md/图片/slides 预览/表格/PDF），下载与“打开所在文件夹”；
     - 评审：最新轮雷达图 + 维度条、三位评委意见、机器检查项通过/失败表、硬门槛状态、历史轮次折线；
     - 修改对话：对话框发指令（“把第三章的竞品表扩到 5 家”），复用会话；
     - 数字账本：关键数字表（名称/值/单位/来源文件/被引用位置）与一致性标记。
  4. **设置**：API 地址/密钥/模型（拉取 `/v1/models`、连通性测试）、角色模型映射、搜索 API、Python 路径与依赖检查、字体检查、输出目录、主题。
  5. **知识库**：赛道规范、官方评分表、获奖作品特征、评委 50 问、数据源清单、使用手册（渲染 `agent/knowledge/*.md` 与 `tracks/*/spec.md`）。

## 6. 质量闭环

1. **生成前**：任务解析必须列出评分映射表（评分项 → 章节 → 关键证据），缺项即补入大纲。
2. **生成中**：主笔每章写完后运行 `bizkit.checks --chapter`（数值密度、来源配对、禁语）；财务模型 `Checks` 页必须为 0。
3. **生成后**：三评委独立打分（不同模型）→ 取均值；质检员机器检查；合并《修改清单》（按“得分提升/工作量”排序）。
4. **迭代**：目标线默认 85（国金线），每维 ≥3；最大 3 轮；每轮只改清单中的问题，重新渲染与检查；收敛（提升 <2 分）即停止。
5. **交付**：评审报告与机器检查随提交包一起给用户；《数据真实性清单》列出需替换的示例数据。

## 7. 运行环境

- Node ≥ 20（用户机 22.21，VM 22.23）、pnpm 10；Python ≥ 3.10（用户机 Anaconda 3.12）；Windows 需 Git for Windows（Bash 工具）；可选：Word/WPS 或 LibreOffice（PDF）。
- 安装：`pnpm install`（含 `@anthropic-ai/claude-agent-sdk` 平台二进制）；`python -m pip install -r agent/toolkit/requirements.txt`。
- 启动：`pnpm dev`（server 5178 + web 5173）或 `pnpm build && pnpm start`（单端口 5178 托管前端）。
- 密钥：环境变量或设置页；仓库内无任何密钥。

## 8. 测试与验收

- 单元：工具包渲染（docx/pptx/xlsx 可打开、字体/样式正确）、检查器阈值、搜索解析器、事件归一化。
- 集成：VM 上用真实网关跑完整流程（品牌策划 + 创业计划各 1 例），记录时长、Token、评分曲线。
- 用户机：Windows 下完整复跑一例；核对 Word/PPT/Excel 打开效果、PDF 页数、匿名检查。
- 验收标准：9 阶段全部 ✅；提交包齐全；评审总分 ≥85；机器硬门槛全部通过；跨材料一致性 100%。
