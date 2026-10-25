# 金策 BizGold —— 全国高校商业精英挑战赛 国金级材料全自动生成规范 v1.1

> 适用：全国高校商业精英挑战赛（CUBEC）各策划书 / 计划书类赛道（品牌策划、创业计划、数字营销、案例分析、会奖旅游……）。目标：**全国总决赛一等奖（“国金”）水准**——章节与官方评分项一一对应、每个数字可追溯、财务闭环、创意可一句话复述、合规零扣分，并经三位评委按量规打分达标。
> 本文件只规定**流程、质量门槛与写作风格**；赛道差异一律见 `tracks/<赛道>/spec.md`（章节大纲、每章必答、评分映射）；Word / PPT / Excel / 图表由 `toolkit/bizkit` 渲染，你只写内容（Markdown / JSON）。

---

## ⛔ 最高优先级规则（违反任何一条即视为任务失败）

1. **全自动**：禁止向用户提问、确认或等待；遇到决策点自行选择最合理方案，并把依据写进 `00_任务解析/任务解析.md` 的“关键假设”表。
2. **三份文件先读后做**：`任务书.md` → 本文件 → `tracks/<赛道>/spec.md`。内容冲突时：任务书 > 赛道规范 > 本文件；流程冲突时：本文件优先。
3. **数字溯源铁律**：正文、PPT、讲稿、题库里出现的每一个数字只能来自三处：① `03_财务/numbers.json`（财务、预算、KPI、调研关键数字）② `01_调研/数据/*.csv` 或调研笔记中带来源的条目 ③ `任务书.md` 与 `inputs/_preread/` 用户资料。禁止凭印象、估算或“看起来合理”地填数。外部数据写足四要素：**数值 + 口径 + 年份 + 来源**，并在参考文献中有对应条目。
4. **示例数据铁律**（`researchMode = synthetic` 时）：问卷、访谈、实地观察等一手调研结果由你模拟生成时，所在段落开头必须带标记 **【示例数据，须替换为真实调研结果】**；`numbers.json` 对应条目 `synthetic: true`；PPT 页 `note`、讲稿段落末尾、题库该题“证据”后同样加 `（示例数据，须替换）`；全部写入《数据真实性清单》与《团队确认清单》。禁止把模拟结果包装成真实调研。`researchMode = provided` 时只使用用户资料与公开数据，不得虚构一手调研。
5. **匿名铁律**：全部产出（封面、页眉页脚、正文、图片、PPT 母版与备注、文件属性、讲稿）不得出现院校名称、校徽、院系、指导教师姓名、学生姓名、学号、电话、邮箱。团队只以“团队名”出现，署名用“XX 品牌策划工作室 / XX 团队”。机器检查命中即修复。例外：会计与商业管理案例赛道（`subject.caseVariant = accounting`）细则要求封面写选手与指导教师姓名——封面对应字段一律用占位 `【提交前由团队填写：选手姓名】`、`【提交前由团队填写：指导教师姓名】`，正文、页眉页脚、PPT、讲稿仍不得出现任何姓名；打包 README 提醒填写。
6. **来源真实**：参考文献只写真实存在、可访问的来源（`web_fetch` 成功抓到，或国家统计局 / 部委 / 上市公司年报 / 行业协会 / 头部咨询机构的已知出版物）。不得编造 URL、报告名、数据；无法核实的数据标 **【待核实】** 并降级为定性表述。
7. **进度协议**：每个阶段开始与结束调用 `mcp__bizgold__update_progress({stage, status, note, round?})`，同时在对话中输出完整进度清单（格式见第五章；评审阶段 round 从 0 起）。阶段 id：`parse / research / strategy / finance / charts / write / pitch / review / package`。
8. **不读二进制**：PDF / 图片 / Office 文件不要直接 `Read`；用 `inputs/_preread/<文件名>.txt`，没有预读文本就运行 `python -m bizkit.preread <file> --out inputs/_preread`。确需看图时才用 `Read` 读 PNG/JPG。
9. **Windows 兼容**：你写的每个 `.py` 脚本文件前两行固定为
   ```python
   import sys, io
   sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
   ```
   `python -c` 一行命令不需要这两行，但命令前必须有 `PYTHONUTF8=1`。生成 Word / PPT / Excel / 图表**一律用工具包 CLI**（`python -m bizkit.xxx`），不要手写 python-docx / python-pptx 代码（模块缺失时按规则 14 的降级表处理）。调用 Python 统一用 `"${BIZGOLD_PYTHON:-python}"`，并先 `export PYTHONPATH=toolkit PYTHONUTF8=1`（Windows Git Bash 同样写法）。
10. **不得结束回复**：直到 `08_提交包/*.zip` 生成、进度清单全部 ✅、并输出交付摘要。中途任何失败都要修复或降级后继续。
11. **子代理**：允许并行派发（研究员 ≤ 4、主笔 ≤ 2、评委 = 3）。每个子代理的任务书先写入 `任务/<角色>_<编号>.md`，`Agent` 的 prompt 只给文件路径与一句话目标。派发、异步回收、核对与重派规则见第四章。
12. **重试上限**：同一命令失败最多重试 3 次；第 3 次仍失败则降级（换图型、减少检索、缩小范围）并在《任务解析》“异常记录”登记。**反常结果先当 bug 查**：盈亏平衡 < 6 个月、毛利率 > 80%、SOM > SAM 的 10%、样本有效率 100%、ROI > 10 倍、KOL 单价与刊例差 3 倍以上——先查假设与脚本，不得直接写进正文。
13. **文件纪律**：产出文件名用中文；不产生无关文件；临时文件放 `tmp/`，打包前删除；不修改 `inputs/`、`toolkit/`、`tracks/`、`knowledge/`。
14. **环境自检与降级**：Step 0 第一个动作是运行工具包自检（Step 0 动作 0），结果写进 `00_任务解析/环境自检.md`。缺失的模块按「七、7.1 降级表」处理：**降级是预案，不是失败**，不得因模块缺失停止流程或向用户提问；每次走降级路径都要在《任务解析》异常记录登记一行（模块 / 替代做法 / 对交付物的影响）。规则 9 的“一律用工具包 CLI”在该模块缺失时让位于 7.1 的替代做法。

---

## 项目结构（工作区）

```
<工作区>/
├── 任务书.md                     ← 服务端生成：赛道、组别、对象、参数、资料清单（先读）
├── CLAUDE.md                     ← 本文件
├── project.json                  ← 服务端维护，只读
├── inputs/                       ← 用户资料；_preread/ 下是预读文本（只读）
├── tracks/<赛道>/                 ← spec.md 章节大纲与评分映射 / rubric.json / checks.json / naming.json / outline.json
├── knowledge/                    ← 评分表与扣分细则、国金写作标准、财务模型标准、市场测算与数据源、路演与答辩、命题企业与获奖特征、来源索引
├── toolkit/bizkit/               ← 渲染与检查工具包（python -m bizkit.*）；README.md 为 CLI 手册
├── .claude/agents|skills/        ← 子代理与技能
├── 任务/                         ← 子代理任务书（编号规则见第四章）
├── 00_任务解析/  任务解析.md  环境自检.md
├── 01_调研/  00_调研总览.md  01_宏观与行业.md  02_消费者与渠道.md  03_竞品与标杆.md  04_对象企业与政策.md  05_消费者洞察.md
│            数据/*.csv  数据/关键数字_*.json  数据/关键数字_汇总.json  参考文献_A…D.json  参考文献_消费者.json  参考文献.json  文献编号映射.json
│            调研工具包/（问卷.md 访谈提纲.md 执行计划.md 示例结果.md）
├── 02_策略/策略框架.md
├── 03_财务/  assumptions.json  财务模型.xlsx  numbers.json  财务摘要.md  建模说明.md  财务图表.charts.json
├── 04_图表/  图表需求单.md  charts.json  *.png  index.json
├── 05_正文/  00_概要.md  01_*.md … A1_*.md  新增文献_<批号>.json（合并后删除）  book.json
├── 06_路演/  slides.json  讲稿.md  答辩题库.md  视频分镜.md
├── 07_评审/第N轮/  机器检查.json|md  质检整改清单.md  评审_investor|professor|executive.json|md  评审汇总.json  修改清单.md
├── 产出/     <作品名>.docx  <作品名>.pdf  路演.pptx  路演.pptx.preview.json
├── 08_提交包/ 按命名规则重命名的全部交付物 + 声明模板 + 数据真实性清单.md + 团队确认清单.md + README.md + 提交包.zip
└── tmp/
```

---

## 一、执行流程（Step 0–8；每步“输入 → 动作 → 验收 → 产出”）

> 提速与上下文铁律：① 任何文件只 `Read` 一次，需要回看先看自己的笔记（《任务解析》《调研总览》）；② > 300 行的文件先 `Grep` 定位再按 `offset/limit` 分段读所需部分；③ `knowledge/` 只按需定点读（见 Step 0），`knowledge/评审量规总表.md` 与 `tracks/<赛道>/rubric.json` 重复，不读；④ 已有结果直接复用；⑤ 每完成一步重新输出进度清单，并附一行“上下文自评：已读文件 N 个，下一步派子代理 / 写文件而非对话输出”；⑥ 一次 `Bash` 不超过 10 分钟、单个脚本 ≤ 5 分钟，长任务拆开；⑦ 大段内容（策略框架、charts.json、book.json）写文件，不在对话里复述。

### Step 0 读题解析（parse）

输入（按此顺序与方式读，总量控制在 25 K 字以内）：
1. `任务书.md`（全文）。
2. `tracks/<赛道>/spec.md`（全文）与 `rubric.json`（全文）；`checks.json`、`naming.json`、`outline.json` 只在用到字段时读。
3. `knowledge/评分表与扣分细则.md`：先 `Grep -n "^### "` 列出小节，只读本赛道评分表小节 + “跨赛道通用规则速查”与“扣分细则”小节；其他赛道不读。
4. `knowledge/命题企业与获奖特征.md`：仅命题 / 指定企业赛道读，且先 `Grep` 企业名或命题关键词，只读命中的条目与“高等级获奖作品的共性特征”小节。
5. `inputs/_preread/*.txt`：每个文件先读前 120 行，再 `Grep` 任务书中的关键词（品牌名 / 产品 / 区域 / 目标）命中行 ±20 行；> 600 行的文件写 `tmp/预读索引_<文件>.md`（每 100 行一句摘要）后只读索引。
6. `knowledge/国金写作标准.md`、`财务模型标准.md`、`市场测算与数据源.md`、`路演与答辩.md` 由对应子代理或对应 Step 按需读，Step 0 不读。

动作：
0. **环境自检**（必须最先做）：
   ```bash
   export PYTHONPATH=toolkit PYTHONUTF8=1
   for m in preread charts docx_render pdf package fonts pptx_render finmodel checks search; do
     "${BIZGOLD_PYTHON:-python}" -c "import importlib; importlib.import_module('bizkit.$m')" 2>/dev/null \
       && echo "OK      bizkit.$m" || echo "MISSING bizkit.$m"
   done
   "${BIZGOLD_PYTHON:-python}" -m bizkit.fonts 2>/dev/null | head -5
   ```
   把输出写入 `00_任务解析/环境自检.md`（表：模块 / 状态 / 降级做法 / 影响的 Step）。`charts`、`docx_render`、`package` 任一缺失视为环境错误：重试导入 1 次后仍缺失，记录异常，继续按 7.1 降级。
1. 写 `00_任务解析/任务解析.md`，七部分：
   1. **对象解读**：品牌 / 项目 / 命题的事实清单（主体、产品、区域、阶段、已有资产、命题目标 KGI / KPI），逐条标来源（任务书 / 资料第几页）。命题企业须逐条列出命题要求。
   2. **评分映射表**（固定表头：官方评分项 / 分值 / 对应章节 / 该章必须出现的证据（表 / 图 / 数字 / 一手资料）/ 目标得分）。覆盖率必须 100%。
   3. **材料清单与格式约束**：页数上限与目标页数（目标 = 任务书 pageTarget，且 ≤ 上限 × 90%）、字体字号、封面字段、命名规则、PPT 页数区间（`checks.json.pitch.deckPages`）、讲稿时长、视频要求、匿名与 AI 规定。
   4. **关键假设**：每条格式 `假设 → 依据 → 若不成立的影响`，≥ 5 条。
   5. **研究问题清单**：12–20 个具体问题，按研究员 A–D 分组（A 宏观与行业 / B 消费者与渠道 / C 竞品与标杆 / D 对象企业与政策），每个问题注明期望产出（数据表 / 事实 / 引用）。
   6. **异常记录**：留空表（时间 / 环节 / 现象 / 处理 / 对交付的影响）。
   7. **页数预算表**：读 `tracks/<赛道>/outline.json` 的 `chapters[].pages` 与 `appendices`，按 `pageTarget` 等比缩放：章节合计 = pageTarget × 0.85（取整到 0.5 页），附录合计 = pageTarget × 0.15 且 ≤ spec 附录上限；每章目标页 = 该章区间中值 ÷ 全部中值之和 × 章节合计，不得低于区间下限的 80%。表字段：文件名 / 目标页 / 目标字数（650 字/页，减去图表占页：图 0.45 页、表 40 行 1 页）/ 必配图表数。每个一级标题另起一页，预算已含该空白。该表是 Step 5 任务书与超页裁剪的唯一依据。

验收：环境自检文件存在；评分映射 100% 覆盖；研究问题 ≥ 12；页数预算表合计 = pageTarget ± 1 页。

### Step 1 调研（research）

动作：
1. 并行派发 4 名 **行业研究员**（`industry-researcher`），任务书 `任务/研究员_A.md … 研究员_D.md`（每份写明预算：检索 ≤ 12 次、抓取 ≤ 8 次、≤ 15 分钟；`web_search` 返回空时换更具体的关键词或直接 `web_fetch` 已知权威站点，不要反复重试同一查询），输出文件固定：A `01_调研/01_宏观与行业.md`（规模、增速、结构、政策、趋势，TAM 口径数据）；B `02_消费者与渠道.md`（人群、动机、决策链、渠道格局、价格带）；C `03_竞品与标杆.md`（3–5 家直接竞品 + 1–2 个跨界标杆，维度 ≥ 6，含传播案例）；D `04_对象企业与政策.md`（企业资产、产品线、现有传播，只做事实不做赞美；创业项目的技术路线与知识产权；命题区域产业政策）。任务书六段格式见第四章：背景、目标问题（抄自任务解析第 5 部分）、输入（只列必需文件）、输出文件与结构、约束（每条数据带来源；检索 ≤ 12 次；不编造 URL；文献写 `参考文献_<字母>.json`）、回报格式。
2. 同一回复里派发 1 名 **消费者洞察**（`consumer-insight`，任务书 `任务/消费者洞察_1.md`，输出 `01_调研/05_消费者洞察.md` 与 `调研工具包/`）：用户画像（2–3 个，含一句洞察）、问卷（≥ 25 题）、访谈提纲（≥ 12 问）、调研执行计划（样本量 n = Z²p(1−p)/e²，95% / 5% → 384；配额与质控；质量门槛 α ≥ 0.7、KMO ≥ 0.6、有效率 ≥ 80%）、`synthetic` 模式下生成**标注的示例结果**（样本结构表 + 5–8 个关键发现 + 交叉表 1–2 张 + α / KMO / 有效率三项指标，CSV 落盘）。
3. 回收与合并（全部研究员与消费者洞察回收并通过第四章核对之后）：
   a. **文献合并**：读 `01_调研/参考文献_A.json … _D.json`、`参考文献_消费者.json`（文件可能是数组，也可能是 `{"references": [...]}` 对象，两种都要接受；条目字段：`id type author title publisher year url accessed fetched`），按 URL 去重（同 URL 保留 `fetched: true` 者），按 A → B → C → D → 消费者的顺序顺延编号写 `01_调研/参考文献.json`，并写 `01_调研/文献编号映射.json`（`{"A1": 1, "A2": 2, "B1": 9, …}`）；用一段 Python 把 `01_调研/*.md` 与 `01_调研/数据/*.csv` 中的临时号替换为正式号（`[A1]` → `[1]`）。之后任何环节新增文献只能走「`新增文献_*.json` → 主代理合并」，不得直接编辑 `参考文献.json`。
   b. 写 `01_调研/00_调研总览.md`——关键事实 20 条（每条四要素 + 正式文献号）、数据表索引（CSV 文件名 → 字段 → 用途）、信息空白与【待核实】清单、对策略的 5 条启示。
   c. 核对 `参考文献.json`：条数 ≥ `checks.json.sourceMin`、独立域名 ≥ 8、官方 / 年报 / 协会来源占比 ≥ 40%、近 3 年 ≥ 70%、`fetched: true` 占比 ≥ 60%。不达标**最多补派 1 名研究员一次**（任务书只列缺口，检索 ≤ 10 次、抓取 ≤ 8 次、≤ 12 分钟）；补派后仍不达标就在调研总览“信息空白”里登记差距并**继续下一步**——这些指标是评审扣分项，不是阻塞门槛，调研阶段总时长控制在 40 分钟以内。
   d. **关键数字汇总**：把 `01_调研/数据/关键数字_*.json` 合并为 `01_调研/数据/关键数字_汇总.json`（同 key 冲突取官方 / 更近年份，并在总览注明）；**此时不记账**，记账在 Step 3 验收后执行。

验收：5 份调研笔记 + 调研总览存在；CSV ≥ 4 个；`参考文献.json` 达标且无临时号残留（`grep -rn "\[[A-D]\d" 01_调研` 为空）；示例数据全部带标记。

### Step 2 策略（strategy）

动作：写 `02_策略/策略框架.md`（≤ 3000 字正文 + 附表）：
1. 问题树（≤ 3 个核心问题）与洞察句（人群 + 场景 + 冲突）。
2. **定位逻辑链**（品牌 / 营销类，细则列明的 5 步）：宏微观环境预测 → 市场细分 → 目标市场 → 竞争对手异同（竞品矩阵）→ 独特性与风格；定位陈述句：“对于（目标人群），X 是（品类）里（差异点）的品牌，因为（支持点）。”核心价值须给出“企业内部认可 + 目标顾客认可”的双验证证据。
3. 商业模式（创业 / 案例类）：商业模式画布九格、收入公式、单位经济目标、护城河。
4. **Big Idea**：按 `.claude/skills/brand-creative/SKILL.md` 的“10 → 3 → 1”流程发散收敛（5 个张力源各 2 个候选 → 7 项清单打分 → 前 3 对照 → 选 1 并说明为何不选另两个），候选表作为附表；最终给一句创意主张 + 视觉锤 + 语言钉，并说明如何延展到 ≥ 3 个渠道与 ≥ 1 个线下场景；给出与竞品传播主题的差异矩阵。
5. **各章结论句**：按 `outline.json` 每章一句（带数字），后续写作围绕结论句展开。
6. **数字需求单**：列出财务阶段必须产出的 key（契约 §7.5 固定 key 集合）与策划案预算上限。
7. 标题：按 `outline.json.titlePattern`（品牌：“对仗主标题（8–14 字，含文化 / 情绪意象）——品牌名 + 品牌策划书”；创业：“项目名——一句话定位 商业计划书”；其他赛道见各自 outline）。
8. 命题 / 指定企业：写“命题回应表”（命题要求 / 本方案对应章节 / 关键举措 / 可量化目标），供引言章原样使用。

验收：定位陈述句与 Big Idea 各一句；候选表 ≥ 10 条；章结论句数量 = 章节数；数字需求单完整。

### Step 3 财务（finance）

动作：写 `任务/建模师_1.md`，派发 **财务建模师**（`financial-modeler`）：
- 赛道专项技能：案例赛道先读 `case-analysis`，会奖赛道先读 `mice-program`，数字营销先读 `digital-video`（主笔与路演同样适用）。
- 创业计划 / 案例类：`mode = startup`，`assumptions.json` 每条假设带 `source`（调研数据 / 行业对标 / 任务书），运行 `python -m bizkit.finmodel 03_财务/assumptions.json`，要求“检查”页全部为 0、出让比例 = 融资额 ÷ 投后、盈亏平衡在 12–36 个月之间（否则重新审视假设并记录理由）、收入曲线非等差等比、悲观情景现金最低点 + 融资额 ≥ 0。
- 品牌 / 数字营销 / 会奖：`mode = budget`，预算 = Σ 单价 × 数量（单价有来源：平台刊例、同类案例、供应商报价区间）、机动金 5–10%、总额 ≤ 命题上限（会奖 200–400 万）、KPI 表 ≥ 5 行（指标 / 基线 / 目标 / 工具 / 周期）、ROI 或品牌资产指标。
- 产出 `numbers.json`（固定 key）、`财务模型.xlsx`（budget 模式同名文件，内含预算总表 / 排期 / KPI 页）、`财务摘要.md`、`建模说明.md`（假设来源表、对标理由、最敏感因子、评委追问要点）、`财务图表.charts.json`。
- **记账合并（主代理，建模师回收并核对后立即做，顺序不可颠倒）**：调用 `mcp__bizgold__record_numbers(entries = 01_调研/数据/关键数字_汇总.json 的全部条目)`。条目 key 若不在契约 §7.5 固定 key 集合内，统一加前缀 `research_`（如 `research_price_median`）；示例数据保留 `synthetic: true`；`source` 必须是 CSV 路径或带文献号的调研笔记位置。然后重新 `Read 03_财务/numbers.json`，确认同时含财务 key 与 `sample_size`、`tam`（budget 模式为 `budget_total`）。此后任何人不得直接编辑 `numbers.json`：财务数字改动只能改 `assumptions.json` 重跑 finmodel，调研数字只能通过 `record_numbers` 追加。

验收：`numbers.json` 覆盖数字需求单全部 key 且含调研 key；检查页为 0；预算 ≤ 上限；`财务摘要.md` 抽 5 个数字与 `numbers.json` 一致。

### Step 4 图表（charts）

动作：
1. 主代理写 `04_图表/图表需求单.md`：按 `outline.json` 每章的必备图表列出（id / 图型 / 结论句 / 数据来源文件（CSV 或 numbers.json key）/ 用于哪一章）；把 `03_财务/财务图表.charts.json` 的条目直接并入需求单；总数 ≥ `checks.json.chartMin`。
2. 写 `任务/制图_1.md`，派发 **制图员**（`chart-maker`）：按需求单读 CSV / numbers.json，写 `04_图表/charts.json` 并运行 `python -m bizkit.charts`，每图带 `source`，`Read` 抽查 3 张图确认中文与数值无误；失败项用 `--only` 修复重跑。
3. 回收核对：`04_图表/index.json` 条数 ≥ chartMin、全部文件存在、无中文乱码（主代理再抽看 1 张）。品牌视觉以“色板 / 字体 / 图形元素说明表 + 应用场景清单”呈现（表格，不画 LOGO）。

验收：index.json 达标；需求单中的每一项在 index 中有对应 id。

### Step 5 撰写（write）

前置：`00_任务解析/任务解析.md` 第 7 部分「页数预算表」存在；`04_图表/index.json` 与 `03_财务/numbers.json` 就绪。

动作：
1. **分批**：章节清单 = `outline.json.chapters`（含 `00_概要`）+ `appendices`。去掉 `00_概要` 后按顺序每 3 章一批，附录单独一批，`00_概要` 单独放在最后一批（它要引用各章结论与数字）；批数 = ⌈(章数 − 1) ÷ 3⌉ + 2。**同时最多 2 批并行**：先派第 1、2 批，全部回收并核对后再派第 3、4 批，依此类推。
2. 每批写任务书 `任务/主笔_<批号>.md`（六段格式，≤ 1500 字），必须包含：本批章节的文件名与一级标题（照抄 outline）、每章目标页与目标字数（来自页数预算表）、每章结论句（来自 `02_策略/策略框架.md`）、每章必答（来自 `outline.mustAnswer`）、本章可用图表 id 与标题（来自 `04_图表/index.json`，按章分配）、可用数字 key（来自 numbers.json）、命题回应表（引言批）、本批禁止触碰的文件（其他章、`book.json`、`01_调研/参考文献.json`）。prompt：“读 `任务/主笔_<批号>.md` 并完成；只写本批章节文件；完成后返回 ≤ 300 字摘要。”
3. 主笔每章写完运行 `bizkit.checks --chapter`（模块缺失时跳过，由第 5 步统一检查）。主笔新增引用一律写到 `05_正文/新增文献_<批号>.json`（字段同参考文献.json，`id` 用 `W<批号>-<序号>` 临时号），正文暂用 `[W1-3]` 这类临时号引用。
4. **回收核对**（每批，按第四章）：文件齐全且标题与 outline 一致；每章首段为加粗结论句；每章 ≥ 1 图 1 表；章末有「本章小结」；示例数据段带标记；每章字数在目标字数 ± 15%；无 `[待补]` / `TODO`；必答问题逐条有回应（摘要中列出未答项则补调研或降级为【待核实】）。不合格带具体问题重派该批一次。
5. **合并与渲染（主代理）**：
   - 合并全部 `05_正文/新增文献_*.json` 进 `01_调研/参考文献.json`（按 URL 去重、顺延编号），用 Python 把 `05_正文/*.md` 中的临时号替换为正式号，删除 `新增文献_*.json`；
   - 写 `05_正文/book.json`（契约 §7.2；封面字段来自任务书与 `naming.json.coverFields`；`anonymize.forbidden` 填任务书中出现的院校词；`pdf: true`）；
   - 运行 `"${BIZGOLD_PYTHON:-python}" -m bizkit.docx_render 05_正文/book.json`，读输出 JSON：匿名命中 → 回到 md 源头修复后重渲染；记录页数（有 PDF 用 `bizkit.pdf --pages`，否则用估计值）。
6. **页数控制**：目标 = pageTarget ± 10%，上限 `checks.json.pageLimit` 不可超。超页：对照页数预算表找超支最多的章，先删重复图、再合并表、再压缩附录（问卷缩小字号、CSV 明细移出正文），每次调整后重渲染；不足：补读图句、与基准对照、敏感性分析，不放大图片凑页。调整 ≤ 2 轮，仍不达标登记异常并继续。

验收：全部章节与附录 md 存在；`grep -rn "\[W[0-9]" 05_正文` 为空（无临时号残留）；docx 生成；页数在目标区间且 ≤ 上限；章检查无 fail（或已登记降级）。

### Step 6 路演（pitch）

动作：写 `任务/路演_1.md`（含：各章结论句、`numbers.json` 关键 key 与值、`04_图表/index.json` 可用图、`checks.json.pitch` 的页数区间 / 字数区间 / 题库维度、示例数据 key 列表），派发 **路演设计师**（`pitch-designer`；数字营销赛道另读 `digital-video` 技能）：
1. `06_路演/slides.json`：正文页数 = `checks.json.pitch.deckPages` 区间（备用页 ≤ `deckPagesAppendixMax`，layout `appendix` 放最后）；每页标题为结论句（含数字或判断），正文 ≤ 120 字、要点 ≤ 6 条；标题 ≥ 32 pt、正文 ≥ 24 pt、图注 ≥ 16 pt、全稿 ≤ 3 色（溢出只能删字不能缩字）；数字全部来自 `numbers.json` 与正文；图来自 `04_图表/`；每页 `note` 写讲稿段落（引用示例数据的页在 `note` 末尾加 `（示例数据，须替换）`）；运行 `python -m bizkit.pptx_render 06_路演/slides.json`。
2. `讲稿.md`：字数区间 = `checks.json.pitch.scriptWords`（10 分钟 ≈ 1 750 字），前 150 字内出现场景 / 数字钩子与一句定位；每段首句结论；段落时间分配见 `pitch-qa` 技能；每段标注 `[0:00–0:45]` 与讲者角色；引用示例数据的段落末尾加 `（示例数据，须替换）`。
3. `答辩题库.md`：≥ `checks.json.pitch.qaMin` 题，维度 = `checks.json.pitch.qaDimensions` 各 ≥ 4；格式 `### Q1【维度】问题` + `主答：角色` + `结论（10 秒）/ 证据（30 秒）/ 行动（10 秒）`，单答 ≤ 150 字；尖锐题 ≥ 8（含“示例数据如何替换为真实调研”）。
4. `视频分镜.md`：10 分钟陈述视频分镜（时间轴、画面、出镜人角色、PPT 页码、台词要点、字幕），开头 3–5 名队员真人出镜，结尾致谢。

验收：pptx 生成且警告 = 0；讲稿字数在区间；题库数量与维度达标；分镜含出镜提示；示例数据标注已传播到三份文件。

### Step 7 评审（review）

> 一轮 = 质检 → 三评委 → 汇总 → 判定 →（未达标）修改 → 重渲染。轮数上限 `reviewRounds`（默认 3）。主代理在本步**只做调度、落盘、汇总与判定**：不亲自改稿、不替评委打分、不改评委给出的分数。

动作（第 N 轮，N 从 1 起）：
1. **质检**：`mkdir -p 07_评审/第N轮`；写 `任务/质检_第N轮.md`（轮次、赛道、作品 docx / pdf 路径、上一轮修改清单路径），派发 **质检员**（`qa-inspector`）。回收后确认 `07_评审/第N轮/机器检查.json` 存在且含 `summary`，并读 `质检整改清单.md`；若 `bizkit.checks` 不可用，按规则 14 用 `tmp/checks_lite.py` 生成同名文件并置 `summary.degraded = true`。
2. **评委任务书**：按文末「附 A 评委任务书模板」写三份 `任务/评委_investor_第N轮.md`、`任务/评委_professor_第N轮.md`、`任务/评委_executive_第N轮.md`。三份只有 `judge` 不同。
3. **并行派发三位评委**（`judge-investor` / `judge-professor` / `judge-executive`）：在同一条回复里连发三个 `Agent` 调用，prompt 固定为“读 `任务/评委_<judge>_第N轮.md` 并完成评审；回复只包含一个 ```json 块和一个 ```markdown 块，不要有其他文字。”按第四章「异步回收」等待三份结果。
4. **落盘**（评委没有 Write 工具）：对每位评委的回复——
   - 取第一个 ```json 块写入 `07_评审/第N轮/评审_<judge>.json`，取第一个 ```markdown 块写入 `07_评审/第N轮/评审_<judge>.md`（用 Write 原样写，不改分数、不改措辞）；
   - 写入前补齐 / 覆盖四个字段：`round = N`、`judge = "<judge>"`、`track = "<赛道>"`、`model`（取任务书“角色模型”表中该评委的模型名；无则填 `"n/a"`）；
   - 校验：`"${BIZGOLD_PYTHON:-python}" -c "import json,sys;d=json.load(open(sys.argv[1],encoding='utf-8'));assert len(d['scores'])==int(sys.argv[2]),'维度数不符';assert 0<=d['total']<=100;print('ok')" 07_评审/第N轮/评审_<judge>.json <rubric 维度数>`；
   - 解析失败：先修复明显格式问题（尾逗号、中文引号、`//` 注释、未闭合括号）再存；仍失败或维度数不符 → **重派该评委一次**（同一任务书，prompt 追加“上次输出的 JSON 无法解析：<错误原文>，请严格重新输出”）；第二次仍失败则该评委本轮缺席：`评审汇总.json.byJudge.<judge> = null`，在《任务解析》异常记录登记。
5. **汇总** `07_评审/第N轮/评审汇总.json`（契约 §7.6）：
   - `byJudge`：三评委 `total`；`score` = 三评委 **中位数**（官方“去掉最高分和最低分后取平均”在 3 位评委时即中位数；缺席 1 人时取两人均值），另记 `mean`；
   - `byDimension`：每个维度三评委均值；`minDimension` = 其中最小值；
   - `machine`：抄机器检查 `summary`；`hardGatePassed`：rubric.hardGates 全部通过（机器检查为准，评委 `hardGates` 为辅，任一来源判失败即失败）；
   - `officialEstimate`：三评委 `officialEstimate.{written,pitch,total}` 各取中位数；
   - `fixes`：合并三评委 `topFixes` + 机器检查 fail 项 + `质检整改清单.md`，按“同文件同位置”去重，按 `gain ÷ 工作量` 排序，取前 12 条；每条字段 `{ id, owner: "writer|pitch|finance|main", file, where, what, accept, gain }`；
   - `converged`：N ≥ 2 且 `score − 上一轮 score < 2`；`remaining`：本轮仍未解决的问题（字符串数组，判定 b 时必填）。
   同时写人读版 `修改清单.md`（表：序号 / owner / 文件与位置 / 现状 / 改法 / 验收 / 预期提升），并调用 `update_progress(stage="review", status="running", round=N, note="第N轮 中位数 xx，最低维度 x.x，硬门槛 通过|未过")`。
6. **判定**（按顺序）：
   a. `score ≥ targetScore` 且 `minDimension ≥ 3` 且 `hardGatePassed` → 评审结束，`update_progress(stage="review", status="done", round=N)`，进入 Step 8；
   b. 否则若 `N ≥ reviewRounds` 或 `converged` → 停止迭代，`remaining` 写清剩余问题，进入 Step 8（交付摘要必须原样列出 `remaining`）；
   c. 否则进入修改轮。
7. **修改轮**（只改 `修改清单.md` 中的条目）：
   - 按 `owner` 拆任务书：`任务/修改_第N轮_主笔.md`（writer 项，> 6 条时按章拆成两份 `_主笔1 / _主笔2` 并行）、`任务/修改_第N轮_路演.md`（pitch 项）、`任务/修改_第N轮_财务.md`（finance 项）。每份按「附 B 修改任务书模板」写：条目表 + 可改文件白名单 + “只改清单项，不得重写整章 / 整份”，并写明“先读 `.claude/skills/revision-round/SKILL.md`”。分别派发 `chief-writer` / `pitch-designer` / `financial-modeler`。`owner = main` 的条目（book.json 封面字段、charts.json 单图数据、参考文献条目、匿名词替换）由主代理直接改。
   - 执行顺序与重渲染（有财务项时从第一步起，否则从受影响步起）：财务 → `bizkit.finmodel` → 受影响图表 `bizkit.charts 04_图表/charts.json --only <ids>` → 正文 → `bizkit.docx_render 05_正文/book.json` → `bizkit.pptx_render 06_路演/slides.json`。任何一步失败按规则 12 处理。
   - 回收核对：逐条对照 `accept` 验收；未完成项带原文重派一次；仍未完成的保留到下一轮清单并标注“第 N 轮未完成”。
   - N = N + 1，回到动作 1。

验收：每轮目录含 `机器检查.json|md`、`质检整改清单.md`、三份 `评审_<judge>.json|md`（缺席者有异常登记）、`评审汇总.json`、`修改清单.md`；`project.json.review` 由 `update_progress(round=N)` 同步；最后一轮 `评审汇总.json.remaining` 与交付摘要一致。

### Step 8 打包（package）

动作：
1. **提交前自检表**（逐项打勾，写入交付摘要）：匿名 0 命中 / 页数 ≤ 上限 / 封面字段 = naming.json / 命名规则占位已说明 / PPT 页数在区间 / 讲稿字数在区间 / 题库 ≥ qaMin / 示例数据条数 = 数据真实性清单条数 / 参考文献无临时号 / numbers.json 与正文抽查 5 个一致 / zip 内文件数 = README 清单数。
2. 按 `.claude/skills/ai-disclosure/SKILL.md` 写 `08_提交包/团队确认清单.md` 与 `08_提交包/内部_AI辅助使用说明.md`（逐项勾选：替换示例数据 / 填学校名 / 补设计稿 / 真人出镜录视频 / 签名指印 / 企业授权书盖章 / 全文通读修改（建议 ≥ 30% 段落重写为团队口吻）/ 4 年同品牌查重自查 / 删除所有“示例数据”标记）。
3. `python -m bizkit.package --track <赛道> --workspace .`；核对 `08_提交包/README.md` 清单；删除 `tmp/`。
4. 输出最终进度清单（全 ✅）与**交付摘要**（固定 6 段）：① 文件清单 ② 各轮分数表（轮次 / 中位数 / 最低维度 / 硬门槛）③ 最终官方口径估分（written / pitch / total）④ 剩余问题（原样抄 `评审汇总.json.remaining`，无则写“无”）⑤ 示例数据清单条数与所在文件 ⑥ 团队待办（抄团队确认清单）。

---

## 二、写作规范（国金级文风，评委打分的直接依据）

### 2.1 总原则
- **结论先行**：章、节、段首句都是结论（含数字或判断），后面给证据与原因。
- **数值四要素**：数值 + 口径 + 年份 + 来源，如“2025 年中国汝瓷工艺品零售规模 33.5 亿元（线上 + 线下，零售口径），资料来源：中国陶瓷工业协会 [3]”。
- **一图表两读图**：每张图表后至少两句：一句形态或数值、一句原因或含义；做不到就删图。
- **一手与二手分开**：一手调研写样本量、时间、地点、方法、有效率；示例数据带标记。
- **对比才有意义**：份额、增速、价格、转化率必须与基准（竞品 / 行业均值 / 上一年）并列出现。
- **可执行**：每个动作有目标、渠道、内容、时间、负责人（角色）、预算、KPI。

### 2.2 禁语与红线
禁语（正文 0 次）：颠覆、独一无二、蓝海、众所周知、前景广阔、填补空白（无证据）、据统计（无来源）、预计将（无依据）、我们学校、老师说、史无前例；广告法绝对化用语：最佳、最好、最优、最先进、第一、第一品牌、全网第一、行业第一、最高级、顶级、国家级、世界级、唯一、独家、极致、首选、100%（“最低点 / 最近三年 / 最终 / 最小 / 最大值”属正常用法，不算）。
弱评价词（每百字 ≤ 0.1）：较好、较高、较大、显著、良好、合理、有效、明显、一定程度。
红线：不复述题目原文；企业荣誉 / 现有做法合计 ≤ 5% 篇幅；不写“AI 生成”之类字样；不出现任何身份信息。

### 2.3 句式与版面
- 句长 ≤ 35 字（讲稿 ≤ 25 字）；段落 3–6 句；主语用“本方案 / 本团队 / 本项目”。
- 术语首次出现给解释；百分比与百分点不混用；金额统一单位（万元或亿元，全文一致）。
- 加粗只用于章结论句与关键数字，每页 ≤ 2 处。
- 图注 / 表注：图题在图下、表题在表上，编号由渲染器生成；图下必写“资料来源：…”。
- 引用：正文用 `[n]`，参考文献按 GB/T 7714-2015 顺序编码由渲染器生成。

### 2.4 篇幅与密度
- 正文页数 = `pageTarget`，上限见 `checks.json.pageLimit`（除封面目录外含附录）。
- 数值密度：策划案全文 ≥ 15 个 / 千字，BP ≥ 20 个 / 千字，市场与财务章 ≥ 25 个 / 千字。
- 图表配额：见 `checks.json.chartMin`（策划案 ≥ 20、BP ≥ 15、案例 ≥ 12），每章 ≥ 1 图 1 表。
- 参考文献：见 `checks.json.sourceMin`；独立域名 ≥ 8；官方 / 年报 / 协会 ≥ 40%；近 3 年 ≥ 70%。

---

## 三、工具包速查（详见 `toolkit/README.md` 与各技能）

| 用途 | 命令 | 输入 | 输出 |
|---|---|---|---|
| 预读资料 | `python -m bizkit.preread <file> --out inputs/_preread` | pdf/docx/xlsx/csv | 文本摘要 |
| 图表 | `python -m bizkit.charts 04_图表/charts.json [--only id,id]` | charts.json | png + index.json |
| Word | `python -m bizkit.docx_render 05_正文/book.json` | book.json + 章节 md | docx（+pdf） |
| PPT | `python -m bizkit.pptx_render 06_路演/slides.json` | slides.json | pptx + preview.json |
| 财务 | `python -m bizkit.finmodel 03_财务/assumptions.json --out 03_财务` | assumptions.json | 财务模型.xlsx + numbers.json + 财务摘要 + 图表 json |
| 检查 | `python -m bizkit.checks --track <赛道> --workspace . [--round N] [--chapter md]` | 工作区 | 机器检查.json/md |
| PDF | `python -m bizkit.pdf 产出/x.docx` / `--pages 产出/x.pdf` | docx / pdf | pdf / 页数 |
| 打包 | `python -m bizkit.package --track <赛道> --workspace .` | 工作区 | 08_提交包 + zip |
| 检索 | MCP `mcp__bizgold__web_search` / `mcp__bizgold__web_fetch`（首选）；命令行 `python -m bizkit.search "关键词"` 为备用 | | |
| 记账 | MCP `mcp__bizgold__record_numbers(entries)` 追加非财务关键数字（样本量、SOM、KPI） | | numbers.json |

命令前统一：`export PYTHONPATH=toolkit PYTHONUTF8=1`，用 `"${BIZGOLD_PYTHON:-python}"` 调用。

---

## 四、子代理调度规范

- **任务书**：路径 `任务/<角色>_<编号>.md`，编号规则固定——研究员 `研究员_A…D`、消费者洞察 `消费者洞察_1`、建模师 `建模师_1`、制图 `制图_1`、主笔 `主笔_<批号>`、路演 `路演_1`、质检 `质检_第N轮`、评委 `评委_<judge>_第N轮`、修改 `修改_第N轮_<角色>`。六段：背景 / 目标问题（编号）/ 输入文件（**只列必需文件，并写“不要读其他文件”**）/ 输出文件与结构 / 约束 / 回报格式。≤ 1500 字；子代理需要的事实（对象、命题目标、结论句、数字 key 与值、角色模型表）直接抄进任务书，不让子代理再翻全库。
- **调用**：`Agent` 的 prompt 固定一句：“读 `任务/xxx.md` 并完成；只写规定的输出文件；完成后返回 ≤ 300 字摘要（含输出文件路径、数据条数、待核实项）”。需要并行时在**同一条回复**里连发多个 `Agent` 调用（研究员 ≤ 4、主笔 ≤ 2、评委 = 3）。
- **异步回收**：`Agent` 可能立即返回 `status: "async_launched"`（后台运行），真正结果随后以 `task_notification` 送达。规则：① 收到 `async_launched` 后不要重复派发同一任务；② 等到该任务的通知到达（或 `task_progress` 显示完成）再做回收核对；③ 等待期间可以做不依赖其结果的工作（写下一份任务书、整理进度清单），但**不得进入下一个 Step**；④ 超过 30 分钟无通知 → `ls -la` 查输出文件，已生成且非空则按文件回收，否则重派一次。
- **回收核对清单（每个子代理必做）**：输出文件存在且非空；文件内无 `TODO` / `[待补]` / `XXX`；每条数据带来源或临时文献号；示例数据带标记；摘要所列文件与实际一致；本角色的专项项（研究员：CSV 最后一列为 `source`；主笔：首段加粗结论句 + 章末小结；建模师：检查页为 0；制图：index.json 齐全；路演：渲染警告 = 0；评委：JSON 可解析且维度数正确）。任一不满足 → 带具体问题重派同一角色一次（prompt 追加“上次缺：…”）；第二次仍不满足 → 降级（缩小范围）并登记异常。主代理**不得**自己“顺手补数”或替子代理编造结果。
- 研究员之间不重复主题；两源数据不一致时取官方 / 更近年份，并在总览注明。
- 评委只读：不得修改任何文件；三位并行；主代理不得替评委打分、不得修改评委给出的分数与评语。
- 子代理默认会加载本 CLAUDE.md；评委与质检员的任务书必须自足（它们可能看不到本文件）。

---

## 五、进度清单格式（每完成一步重新输出）

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
标记放在编号后、名称前；失败用 ❌ 并在下一行写一句原因与处理。

---

## 六、合规与提交（扣分细则：院校信息每处扣 3 分、人数不符扣 3 分、超时扣 1–3 分）

- 书面：A4 竖版 PDF；除封面目录外 ≤ 页数上限；创业计划 ≥ 3000 字、宋体小四 / Times New Roman 小四；封面字段按 `naming.json`（品牌策划：作品名 / 团队名 / 指导企业名；创业计划：组别 / 团队名 / 项目名），其余信息一律不出现。
- 命名：按 `naming.json.fileNameRule`，学校名由用户在提交时填写（打包时用“学校名称”占位并在 README 提醒）。
- 声明：作品原创性声明、作品版权声明、企业授权书（模板由打包生成，提醒手写签名 + 指印、企业盖章）。
- AI 规定：细则禁止“完全使用 AI 工具创作”。本系统定位为**辅助工作台**：所有示例数据必须由团队替换为真实调研，正文须经团队修改与确认；《团队确认清单》《数据真实性清单》、交付摘要三处口径一致。
- 视频：10 分钟陈述视频开头 3–5 名队员真人出镜；分镜已给提示。
- 查重：标题与正文不得套用网传模板段落；同校同品牌 4 年只能用一次（README 提醒）。

---

## 七、常见错误速修

| 现象 | 原因 | 处理 |
|---|---|---|
| 图表中文变方块 | 字体未找到 | 工具包已自动选字体；若仍异常运行 `python -m bizkit.fonts` 查看并在 charts.json 顶层设 `"font": "<可用字体名>"` |
| docx 里出现“【缺图：…】” | 图片路径错 | 图片路径相对章节 md 所在目录（通常 `../04_图表/xxx.png`）；先 `ls 04_图表` |
| finmodel “检查”页不为 0 | 假设不自洽（负库存、负现金未融资、折旧超资本开支） | 读 `财务摘要.md` 的检查说明，调整假设后重跑；不要手改 xlsx |
| 页数超上限 | 图表过多 / 附录过长 | 按页数预算表找超支章：删重复图、两图并排说明改为一图、问卷原件缩小字号、CSV 明细移出正文 |
| 匿名检查命中 | 文中有院校 / 姓名 | 替换为团队名 / 角色；检查页眉、封面、PPT 备注、图片文字 |
| 评委分数低于 3 的维度 | 缺证据表 / 缺一手资料 / 数字无来源 | 对照 rubric 的 5 分锚定补证据，不是改措辞 |
| web_search 返回空 | 提供商被限流 | 换关键词（更具体 + 年份 + 机构名）、改用 `web_fetch` 直接抓权威站点（stats.gov.cn、协会官网、上市公司公告） |
| PDF 转换失败 | 本机无 Word / WPS / LibreOffice | 保留 docx 并在交付摘要提示用户用 Word “另存为 PDF”；页数用估算值 |
| 讲稿超时 | 字数 > 区间 | 删过渡句与形容词，保留数字与结论；每段首句结论 |
| `Agent` 返回 async_launched 后没结果 | 子代理在后台运行 | 等 task_notification；30 分钟无果查文件再决定重派（第四章） |

### 7.1 工具包模块缺失时的降级表（规则 14）

| 缺失模块 | 影响 | 替代做法（唯一允许的例外） | 交付影响登记 |
|---|---|---|---|
| `bizkit.finmodel` | Step 3 无法生成 xlsx / numbers.json / 财务摘要 | 建模师仍写 `assumptions.json`（契约 §7.4）；然后写 `tmp/finmodel_lite.py`（允许使用 pandas / openpyxl）：按假设逐年计算收入、成本、利润、现金流（startup：三表 + 盈亏平衡月 + 现金最低点 + 出让比例；budget：预算 = Σ 单价 × 数量 + 机动金 + KPI 表），把结果写成 `03_财务/财务模型.xlsx`（每张表第一行写公式说明）、`03_财务/numbers.json`（契约 §7.5 固定 key）、`03_财务/财务摘要.md`（含“检查”表：负现金未融资 / 负库存 / 出让比例误差 三项结果）与 `03_财务/财务图表.charts.json`。脚本运行 ≤ 5 分钟 | “财务模型为脚本生成，Excel 内无活公式；提交前请人工复核” |
| `bizkit.pptx_render` | Step 6 无 pptx | 仍完整产出 `slides.json`、`讲稿.md`、`答辩题库.md`、`视频分镜.md`；再写 `tmp/pptx_lite.py`（允许 python-pptx）按 slides.json 生成“标题 + 要点 + 图片 + 备注”的 16:9 素版 `产出/路演.pptx`（正文 ≥ 24 pt，≤ 3 色，母版不含任何院校信息，备注写 `note`）；不做复杂版式 | “PPT 为素版，视觉需团队用模板二次排版；内容与数字以 slides.json 为准” |
| `bizkit.checks` | Step 5 逐章检查、Step 7 机器检查、硬门槛判定全部缺失 | 主代理写一次 `tmp/checks_lite.py`（正则 + 计数，≤ 150 行）并在 Step 5/7 复用，输出与契约 §7.7 同结构的 `机器检查.json|md`（`summary.degraded = true`）。最少覆盖：`anonymity_hits`（用 checks.json 的 anonymityRegex + whitelist 扫 05_正文、06_路演、03_财务、book.json、slides.json）、`forbidden_terms` / `forbidden_terms_adlaw`（forbiddenTerms + advertisingLawTerms − whitelist）、`numeric_density`（每千字数字个数）、`chart_count`（`04_图表/index.json` 条数与正文 `![` 引用数）、`source_count`（参考文献.json 条数、独立域名数）、`struct_coverage`（requiredElements 名称或别名是否出现在一级标题）、`page_limit`（有 pdf 用 `bizkit.pdf --pages`，否则按 650 字/页估算）、`cross_doc_consistency`（numbers.json 每个 value 在正文 / slides.json / 讲稿中出现的写法是否一致：同一 key 的数字出现两种以上数值即 fail）。覆盖不到的检查项写 `"status": "skipped"`，评委对 skipped 项不施加 cap | “机器检查为简化版，`source_pairing / url_alive / font_check` 未执行，提交前请人工核对” |
| `bizkit.search` | 无命令行检索备用 | 只用 MCP `mcp__bizgold__web_search / web_fetch`；MCP 也不可用时研究员只用已知权威出版物并全部标【待核实】，调研总览首行注明“本次调研未联网” | “参考文献未经在线核实” |
| `bizkit.pdf` 转换失败 | 无 PDF、无精确页数 | 保留 docx；页数用 docx_render 输出的估计值；交付摘要提示“用 Word 另存为 PDF 后核对页数 ≤ 上限” | 同左 |

---

## 八、新项目使用

1. 服务端已生成 `任务书.md`、复制资料并预读。
2. 你收到“开始”指令后：读三份文件 → Step 0 … Step 8 → 输出交付摘要。
3. 修改模式（用户后续指令）：只改指定内容，改后重新渲染相关产出并运行 `checks`，不得改 `inputs/`。

---

## 附 A 评委任务书模板（`任务/评委_<judge>_第N轮.md`；`{…}` 由主代理填写）

```markdown
# 评审任务：{judge} 第 {N} 轮

## 固定值（原样写入 JSON）
round = {N}；judge = "{judge}"；track = "{赛道}"；rubric = "tracks/{赛道}/rubric.json"；model = "n/a"（主代理落盘时覆盖）

## 必读（按顺序；只读这些，不要读 knowledge/、tracks/ 其他文件）
1. `tracks/{赛道}/rubric.json`（维度、锚定、machineCapRules、hardGates、officialWeight）
2. `07_评审/第{N}轮/机器检查.json`（只看 `summary` 与 `checks` 中 fail 项；缺失则跳过 cap）
3. `05_正文/00_概要.md`
4. 正文各章（逐章读，大于 400 行的分两段读）：{按 outline 顺序列出 05_正文/NN_*.md 全路径}
5. 附录：{A1…} 只 `Glob` 核对存在并读前 30 行
{N ≥ 2 时追加} 6. `07_评审/第{N-1}轮/修改清单.md` 与 `评审_{judge}.json`：先核对上一轮你提出的 topFixes 是否落实，未落实的在本轮 issues 中重复并标 “【第 {N-1} 轮遗留】”

## 选读（只在打对应维度时读）
`03_财务/财务摘要.md`（财务 / 预算维度）、`06_路演/slides.json`、`06_路演/讲稿.md`、`06_路演/答辩题库.md`（pitchDimensions）

## 评分口径
- 每维 1–5 分（0.5 步长），锚定以 rubric 为准；机器检查 fail 且有 capRule 时 ≤ cap。
- total = Σ(score ÷ 5 × weight)；officialEstimate.written = Σ(score ÷ 5 × officialWeight)，pitch = Σ pitchDimensions(score ÷ 5 × weight) × {官方展示总分} ÷ 100，total = written + pitch。
- evidence 必须写“文件名 + 小节标题 + 你看到的原句或数字”；issues ≤ 3 条 / 维；fixes 写到“在哪个文件哪一节加什么表 / 图 / 数”。
- hardGates 逐项写 `{ "gate": "...", "passed": true|false, "where": "..." }`。
- topFixes 5 条，字段 `{ where, what, gain }`，按 gain ÷ 工作量排序。

## 回复格式
只输出两个代码块：先 ```json（契约 §7.6 完整结构），再 ```markdown（维度表 + 2–3 句评语 + 5 条修改建议）。不要寒暄、不要解释。
```

## 附 B 修改任务书模板（`任务/修改_第N轮_<角色>.md`）

```markdown
# 修改任务：第 {N} 轮 — {角色}

## 条目（只改这些；按序号逐条完成并在回报中逐条确认）
| 序号 | 文件与位置 | 现状（评委 / 机器原话） | 改法 | 验收（主代理将按此核对） | 预期提升 |
|---|---|---|---|---|---|
| {id} | {file}#{where} | {现状} | {what} | {accept} | {gain} |

## 可改文件白名单
{文件列表}。白名单之外的文件一律不动；需要新增数字先查 `03_财务/numbers.json`，账本没有的数字不能写（可改为定性表述并标【待核实】）。

## 规则
- 只改清单项，不重写整章 / 整份；保留原有结论句、图表引用与文献号。
- 新增引用写 `05_正文/新增文献_修改{N}.json`（临时号 `R{N}-1…`），正文用临时号引用（主代理合并）。
- 改完每个文件运行对应检查（正文：`bizkit.checks --chapter`；路演：`bizkit.pptx_render`；财务：`bizkit.finmodel` 检查页为 0）。
- 回报 ≤ 300 字：逐条“序号 → 已完成 / 未完成（原因）”，列出改动文件。
```
