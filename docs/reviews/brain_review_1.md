# 金策 BizGold「大脑」评审报告 #1

| 项 | 内容 |
|---|---|
| 评审对象 | `agent/CLAUDE.md` v1.0、`agent/tracks/{brand,startup,digital,case,mice}/*`、`agent/.claude/agents/*.md`（9）、`agent/.claude/skills/*/SKILL.md`（10）、`agent/.claude/settings.json`、`agent/knowledge/*.md`、`docs/architecture.md`、`docs/contracts.md` |
| 评审依据 | `docs/research/01_赛事情报.md`（§2、§3、§5.4、§6.2）、`02_国金级材料标准.md`（全文）、`03_agent_sdk技术手册.md`（§1.9、§2.6、§3.3、§3.4）、`agent/toolkit/README.md`（工具包实际实现状态）、数模 CLAUDE.md 参考样本 |
| 评审角色 | 商赛国赛评委 × Claude Code / Agent SDK 工程师 |
| 评审方式 | 只读走查 + 文件间交叉比对 + 按 CLAUDE.md 模拟主代理执行；未运行模型、未改动任何现有文件 |
| 严重度定义 | **高** = 不改则流程跑不通或必然触发硬门槛失败 / 扣分；**中** = 可跑通但会造成返工、分数损失或 deepseek/kimi 级模型明显歧义；**低** = 文案、命名与可维护性 |

---

## ① 总评

骨架扎实：九阶段、数字账本、内容/排版分离、三评委量规＋硬门槛，赛道文件对官方评分项覆盖率高。但当前跑不通完整流程：CLAUDE.md 依赖 4 个尚未实现的工具模块且无降级；撰写分批、修改轮执行者、子代理异步回收、参考文献并发写入对弱模型都是歧义点；Step 0 读 60 KB 知识吃掉主代理近半上下文。质量上限卡在创意发散、一手调研标准与官方计分口径未落地。先修 8 个高危项，再谈国金。

---

## ② 问题清单

### 2.1 清单总表

> 高危项的**完整替换文本**在 2.2 节，按编号对应；中/低项的改法直接写在表内（可粘贴）。

| 序号 | 严重度 | 文件与位置 | 问题 | 建议的具体改法 |
|---|---|---|---|---|
| H1 | 高 | `agent/CLAUDE.md` Step 7；`agents/judge-*.md` 「输出」段；`docs/contracts.md` §7.6 汇总行 | ① 修改轮**没有执行者**：第 4 点只写“按清单修改”，主代理自己改稿会撑爆上下文，派谁、任务书长什么样、改后重渲染顺序均未定义；② 评委任务书 `任务/评委_<judge>_第N轮.md` 无模板，评委 JSON 要求 `model` 字段但评委不知道自己的模型名；③ 汇总用“三评委均值”，官方规则是“去掉最高最低取平均”（3 人即中位数）；④ “每个维度 ≥ 3”未说明是谁的分（单评委/均值）；⑤ “提升 < 2 分视为收敛”第 1 轮无参照；⑥ ```json 解析失败的重派流程缺失；⑦ 评委 N ≥ 2 轮不知道上一轮改了什么 | 整节替换 + 两份任务书模板 + 评委 agent 输出段替换，见 **2.2-H1** |
| H2 | 高 | `agent/CLAUDE.md` 规则 9、Step 3/5/6/7、§三命令表；`agents/financial-modeler.md` 第 17 行；`agents/pitch-designer.md` 第 15 行；`agents/qa-inspector.md`；`skills/fin-model`、`pptx-deck` | `toolkit/bizkit/` 实际只有 `charts / docx_render / preread / pdf / package / fonts / theme`（README 标题即“A1 模块”），而大脑无条件调用 `bizkit.finmodel`、`bizkit.pptx_render`、`bizkit.checks`、`bizkit.search`。按规则 9“一律用工具包 CLI、不要手写 python-docx/pptx”，模型会在 Step 3 卡死（重试 3 次后不知如何降级），Step 5 的逐章 checks、Step 6 的 PPT、Step 7 的机器检查全部失效 → 硬门槛无法判定 → 评审循环无法收敛 | 新增规则 14「环境自检与降级」+ Step 0 自检动作 + §七降级表，见 **2.2-H2**；同时把 4 个模块列为工具包 A2 的阻塞项 |
| H3 | 高 | `agent/CLAUDE.md` Step 5；`agents/chief-writer.md` 「收尾」段 | ① “可按章分 2 批并行，每批不超过 3 章”最多覆盖 6 章，品牌 8 章 + 4 附录、创业 13 章 + 5 附录无法分完；② “最后一批负责 book.json”在并行下无“最后”；③ 两批主笔都可能改 `参考文献.json`（竞态）；④ 概要需引用各章数字却可能排在第一批；⑤ 无页数预算，`outline` 区间之和（品牌 48–60、创业 53–71）与 `pageTarget` 无关，附录区间为 `[0,0]`，渲染器 H1 新起页会再吃 10+ 页 | 整节替换 + 收尾段替换 + Step 0 新增「页数预算表」，见 **2.2-H3** |
| H4 | 高 | `agents/industry-researcher.md` 第 20 行；`agents/consumer-insight.md`；`CLAUDE.md` Step 1 第 3 点 | 4 名研究员并行“读入 → 追加 → 写回” `01_调研/参考文献.json`，必然出现覆盖丢失与 `id` 重号（两名研究员同时读到 max id = 7，都写 8），后续正文 `[n]` 全部错位，`source_pairing`/GB/T 7714 直接失分 | 每人只写 `参考文献_<字母>.json`（临时号 `A1…`），主代理回收时合并编号并全局替换，见 **2.2-H4** |
| H5 | 高 | `CLAUDE.md` Step 1 第 3 点 vs Step 3；`docs/contracts.md` §7.5 | Step 1 结束就 `record_numbers` 写入 `numbers.json`，Step 3 的 finmodel 随后“生成” `numbers.json`（契约原文“工具包生成，MCP 可追加”），顺序上调研数字会被覆盖；`关键数字_*.json` 的 key 也无命名约束，与固定 key 集合冲突时无规则 → `cross_doc_consistency` 硬门槛与样本量引用失效 | 记账移到 Step 3 验收后、增加前缀与二次核对，见 **2.2-H5** |
| H6 | 高 | `CLAUDE.md` 规则 11、§四；SDK 手册 §2.6/§3.3 | ① SDK 子代理**默认后台化**（`Agent` 立即返回 `async_launched`，结果经 `task_notification` 送达），CLAUDE.md 完全没写如何等待与回收，弱模型会误判“子代理失败”而重复派发或直接跳到下一步；② 回收核对项分散在规则 11 与 §四，缺统一清单与重派上限；③ 任务书编号规则在 Step 1（A–D）、§四（`<n>`）、研究员 agent（`研究员_*`）、契约 §8（`#2`）四处不一致 | §四整节替换，见 **2.2-H6** |
| H7 | 高 | `tracks/case/spec.md` 第 12/40 行；`tracks/case/checks.json`；`tracks/case/rubric.json`；`docs/contracts.md` §2 `subject`；`CLAUDE.md` 规则 5 | ① 会计案例细则要求封面写选手与指导教师姓名，而 CLAUDE.md 规则 5“全部产出不得出现指导教师姓名、学生姓名”为最高优先级 → 模型要么违规要么被封面卡住；② case `anonymityWhitelist` 含 `中国科学技术大学`（真实院校名进白名单，命中院校信息每处扣 3 分的风险）；③ `subject.caseVariant` 在 spec 中使用但契约未定义；④ 会计 70 分 9 项（含“行业数据聚类分析 5”“问题识别 5”）与 rubric 7 维无映射；流通业“体现地区性特点（10 分）”在 spec / checks 中无任何要求 | 规则 5 加例外、checks/rubric/契约补丁，见 **2.2-H7** |
| H8 | 高 | `CLAUDE.md` Step 0 输入行、“提速铁律”；`agents/judge-*.md` 读法；`knowledge/评审量规总表.md` | Step 0 要求读 `评分表与扣分细则.md`（36 KB，覆盖全部赛道）+ `命题企业与获奖特征.md`（28 KB）+ 全部 `inputs/_preread/*` + spec + rubric，加上自动注入的 CLAUDE.md（26.5 KB），主代理在写第一个文件前已消耗约 60–100 K token（deepseek/kimi 128 K 窗口）；`knowledge/` 还有一份 24 KB 的 `评审量规总表.md` 与 rubric.json 重复；每个子代理默认再加载一遍 CLAUDE.md；评委 maxTurns 40 逐章读 80 页正文无分段策略 | 改为“按赛道/按命题企业定点读取 + 禁止重复读取 + 大文件分段”，见 **2.2-H8** |
| M1 | 中 | `docs/architecture.md` 第 126 行、§4.2、§4.4 | `disallowedTools:['WebSearch','AskQuestion']` 工具名错误（SDK 为 `AskUserQuestion`，错名不生效，模型仍可能提问→违反规则 1）；§4.2 写评委工具“Read / Bash”，agent 文件实为 `Read, Grep, Glob`；§4.4 写“10 种页型”，契约 §7.3 与 `pptx-deck` 技能为 14 种 | 改为 `disallowedTools:['WebSearch','WebFetch','AskUserQuestion']`（WebFetch 走 MCP，避免内置 WebFetch 绕过 `web_fetch` 的来源记录）；§4.2/§4.4 文字对齐实际文件 |
| M2 | 中 | `agents/*.md` frontmatter（9 个） | ① `model: {{MODEL_X}}` 未加引号，YAML 中 `{` 开启 flow mapping，占位符未被替换（本地调试、单测、模板渲染失败）时整份 agent 解析失败并被静默丢弃；② 无 `skills:` 预加载，靠正文“先 Read SKILL.md”——可行但多耗 1 轮且弱模型常跳过；③ 评委与质检员的任务书已自足，仍加载 26.5 KB CLAUDE.md | 每个 agent 改为：`model: "{{MODEL_X}}"`；按角色加 `skills: [fin-model]` / `[pptx-deck, pitch-qa]` / `[web-research, market-sizing]` / `[docx-report]` / `[survey-design]` / `[biz-charts]`；三位评委与 `qa-inspector` 加 `omitClaudeMd: true`（SDK ≥ 0.3.271；低于该版本忽略即可） |
| M3 | 中 | `tracks/*/rubric.json` `officialMapping`；`agents/judge-*.md` 打分规则第 3 条 | `officialEstimate` 的换算只有一句散文，评委各算各的：brand 的 `structure` 维度映射到“展示：内容完整程度 5”（把书面维度算进展示分）；startup 映射引用了不存在的维度 `product`；三评委 `officialEstimate` 不可比 | 每个维度加 `"officialWeight": <分值>`（书面部分之和 = 60 或 70），评委按 `written = Σ score÷5×officialWeight`、`pitch = Σ pitchDimensions(score÷5×weight)×官方展示总分÷100` 计算；startup 删除 `product`，`structure/style` 的 `officialWeight` 设 0（仅内部维度）；写进 judge agent 第 3 条 |
| M4 | 中 | `tracks/brand/checks.json` 第 30–33 行；其余 4 赛道同款字段 | ① `advertisingLawTerms` 含单字 `最`，“最低点 / 最近三年 / 最终 / 最高价”全部误报，硬门槛 `forbidden_terms_adlaw` 永远过不了；② `forbiddenTerms` 缺 CLAUDE.md 2.2 列出的 `据统计`、`填补空白`、`预计将`；③ `weakTermsPer100: 0.1` 没有词表，机器无法执行；④ rubric capRule 用 `risk_matrix`，checks 字段是 `risk_matrix_min`；⑤ `anonymityRegex` 含 `导师`，创业计划“创业导师 / 行业导师”为合理表述 | 替换为：`"forbiddenTerms": ["颠覆","独一无二","蓝海","众所周知","前景广阔","填补空白","我们学校","老师说","据统计","预计将","史无前例"]`、`"advertisingLawTerms": ["最佳","最好","最优","最先进","第一","第一品牌","全网第一","行业第一","最高级","顶级","国家级","世界级","唯一","独家","极致","首选","100%"]`、`"advertisingLawWhitelist": ["最低点","最近","最终","最小","最大值","最高值","首选渠道"]`、`"weakTerms": ["较好","较高","较大","显著","良好","合理","有效","明显","一定程度"]`；rubric 改 `"check": "risk_matrix_min"`；`anonymityWhitelist` 加 `"创业导师","行业导师","导师制"` |
| M5 | 中 | `CLAUDE.md` Step 6 第 1、3 点 | 页数“18–20 页”与题库六维度是品牌口径；创业 `checks.json.pitch.deckPages=[20,24]`、维度为“技术与壁垒 / 商业模式与财务 / 运营与执行”，案例为“财务分析 / 战略逻辑…”；主代理按 CLAUDE.md 验收会与 checks 冲突 | 两处改为引用：“页数区间 = `checks.json.pitch.deckPages`（+ 备用页 ≤ `deckPagesAppendixMax`）”、“题库 ≥ `pitch.qaMin`，维度 = `pitch.qaDimensions` 各 ≥ 4 题” |
| M6 | 中 | `skills/fin-model/SKILL.md` budget 模式；`agents/financial-modeler.md` 第 19 行；`CLAUDE.md` 项目结构树 | ① budget 模式产出写 `预算表.xlsx`，契约与结构树是 `财务模型.xlsx`，打包命名会漏；② “在 `财务摘要.md` 末尾追加建模说明”——该文件由 finmodel 生成，重跑即丢失 | 统一文件名为 `财务模型.xlsx`（budget 模式内含“预算表 / KPI / ROI”三页）；建模说明单独写 `03_财务/建模说明.md`，并在契约 §7 与结构树登记 |
| M7 | 中 | `agents/qa-inspector.md` 第 16–17 行；`CLAUDE.md` 结构树 / Step 7；`docs/contracts.md` §7.7 | 质检员产出 `07_评审/第N轮/质检整改清单.md` 未在任何契约与 CLAUDE.md 中登记，主代理不知道要读它；第 16 行用裸 `python` 而非 `"${BIZGOLD_PYTHON:-python}"` | 结构树 `07_评审/第N轮/` 加 `质检整改清单.md`；H1 替换文本的汇总步骤已把它并入 `fixes`；命令改用环境变量写法 |
| M8 | 中 | `agents/consumer-insight.md` 第 15 行 | 输出文件 `01_调研/0X_消费者洞察.md` 的 `X` 未定义，与研究员 `01_…04_` 可能重号 | 固定为 `01_调研/05_消费者洞察.md`；研究员 A–D 固定为 `01_宏观与行业 / 02_消费者与渠道 / 03_竞品与标杆 / 04_对象企业与政策` 并写入 Step 1 |
| M9 | 中 | `skills/pptx-deck/SKILL.md`；`agents/pitch-designer.md` 第 15 行 | 正文字号下限 18 pt；02 报告 §3 的国赛投影标准是正文 ≥ 24 pt、≤ 3 色、每页 ≤ 6 条；评委隔 8 米看 18 pt 等于没写 | 改为“标题 ≥ 32 pt、正文 ≥ 24 pt、图注 ≥ 16 pt、全稿 ≤ 3 色（主色/强调色/中性灰）；溢出只能删字不能缩字” |
| M10 | 中 | `skills/survey-design/SKILL.md` | 缺国金一手调研的硬指标：Cronbach α ≥ 0.7、KMO ≥ 0.6、有效率 ≥ 80%、访谈 8–15 人、实地观察 ≥ 1 次、交叉表 + 卡方、Kano / 联合分析任选其一；rubric `status` 维度 5 分锚定要求“样本 ≥ 384 + 访谈 ≥ 8 + 实地”，技能却不要求 | 在技能“执行计划”段加“质量门槛表”（指标 / 阈值 / 计算方法 / 未达标处理），示例结果必须报告 α、KMO、有效率三个数（示例数据标记） |
| M11 | 中 | `tracks/brand/outline.json`、`tracks/startup/outline.json` `pages`；`tracks/*/spec.md` 页数表 | 章节页数区间与 `pageTarget` 脱钩，附录为 `[0,0]`；创业 13 章 + 5 附录各起新页，上限 80 页在 `pageTarget=70` 时很容易被突破 | 由 H3 的「页数预算表」在 Step 0 按 `pageTarget` 等比缩放；outline 的 `pages` 改为 `pageShare`（占比），附录给 `[2,12]` 区间 |
| M12 | 中 | `tracks/brand/spec.md` 第 47、61 行 | 附录写 A1–A6，`outline.json` 只有 A1–A4（A5 内容样稿并入 A4，A6 由渲染器生成——spec 未说明）；写作要点第 9 条“数字营销留学生组”属于 digital 赛道 | spec 附录行改为“A1–A4（内容样稿并入 A4；参考文献由渲染器自动生成，不占附录编号）”；第 9 条移到 `tracks/digital/spec.md` |
| M13 | 中 | `tracks/startup/rubric.json` `officialMapping`；`tracks/case/rubric.json` `hardGates` | startup 映射引用不存在的 `product` 维；case 的 `hardGates` 缺 `cover_fields`（会计封面字段最多，最该查） | startup：`创新性 20（innovation 15 + 产品技术方案并入 innovation）`；case `hardGates` 加 `"cover_fields"`，`regional_feature`（流通业变体） |
| M14 | 中 | `CLAUDE.md` Step 6、§六“AI 规定”；`agents/pitch-designer.md`；`skills/pptx-deck`、`pitch-qa` | 示例数据标记只要求出现在正文段落与 numbers.json，PPT 页、讲稿、题库引用同一数字时**不带标记**，用户照读即把模拟数据当真实数据讲给评委；§六只在交付摘要“提示”，缺可核对的《团队确认清单》 | Step 6 加一条：“凡引用 `synthetic: true` 的数字，slides 的 `note` 与讲稿段落末尾加 `（示例数据，须替换）`，题库该题 `证据` 后加同样标注”；Step 8 交付物加 `08_提交包/团队确认清单.md`（替换示例数据 / 填学校名 / 真人出镜 / 签名指印 / 企业盖章 / 全文通读修改 ≥ 30% 段落 / 4 年查重自查，逐项勾选） |
| M15 | 中 | `skills/brand-creative/SKILL.md`；`CLAUDE.md` Step 2 第 4 点 | Big Idea 只规定“一句主张 + 视觉锤 + 语言钉”，没有**发散方法**；弱模型会直接产出第一个想到的常规创意（“国潮 × 年轻化”），rubric `creative` 维很难过 3 分 | 技能加“10 → 3 → 1”流程：按 5 个张力源（品类冲突 / 文化符号 / 消费场景 / 对手空白 / 数据反直觉）各出 2 个候选 → 用 02 报告 §2.4 创意检查清单 7 项打分（满分 35）→ 保留前 3 写一页对照 → 选 1 并说明为何另两个不选；把候选表写进 `02_策略/策略框架.md` 附表 |
| M16 | 中 | `tracks/*/outline.json` 引言 `mustAnswer`；`tracks/*/checks.json requiredEvidence`；rubric `fit` 维 | 命题企业“逐条回应命题要求”只在 Step 0 出现，正文无强制位置，机器无检查，`fit` 维没有 cap | 引言 `mustAnswer` 加“命题回应表（命题要求 / 本方案对应章节 / 关键举措 / 可量化目标）”；checks `requiredEvidence` 加 `"brief_response_table": true`；rubric `fit` 加 `machineCapRules: [{check: brief_response_table, when: fail, cap: 2}]` |
| M17 | 中 | `CLAUDE.md` 规则 12 与 §七；数模 CLAUDE.md 样本 | 缺“反常结果先当 bug 查”“单脚本 ≤ 5 分钟”“禁止重复读取同一文件”“提交前自检表”这类对弱模型最有效的护栏 | 见 ③ 补充内容 B |
| M18 | 中 | `CLAUDE.md` Step 2 / Step 4 | 策略框架与 `charts.json`（20 张图 × 数据数组）全部由主代理亲自写，是主代理上下文第二大消耗点；charts.json 还要求主代理手里有全部 CSV 数据 | Step 4 改为派发 `chart-maker` 子代理（新 agent，见 ③-A），主代理只写“图表需求单”（id / 图型 / 结论句 / 数据来源文件）；Step 2 策略框架限 ≤ 3000 字 |
| L1 | 低 | `CLAUDE.md` 规则 9 | “每个 Python 脚本前两行固定为 …”对 `python -c` 一行命令无法执行，弱模型会纠结 | 改为“每个 `.py` 脚本文件前两行…；`python -c` 一行命令不需要，但命令前必须有 `PYTHONUTF8=1`” |
| L2 | 低 | `docs/contracts.md` §7.6 | `rubric` 字段含义不明（文件名？版本？） | 定义为 `"rubric": "tracks/<赛道>/rubric.json"`，并在 rubric.json 加 `"version": "2026.1"` |
| L3 | 低 | `knowledge/赛事总览.md` | 含主办方联系人姓名、电话、邮箱；主代理可能把它抄进《任务解析》的赛事信息栏，污染匿名检查 | 删除联系方式段，或在文件头加“以下联系方式不得写入任何产出” |
| L4 | 低 | `knowledge/*.md`（8 个） | 全部是自动摘录，正文带 `[S07][R12]` 文献号但无对应 URL 表，子代理无法回查 | 每个文件末尾附“引用索引”（编号 → 标题 → URL），或在文件头注明“编号对应 docs/research/01、02 的来源表” |
| L5 | 低 | `CLAUDE.md` 2.4 vs `tracks/brand/checks.json` | 参考文献下限 CLAUDE.md 写“策划案 ≥ 15 条”，checks 也是 15；02 报告 §1.7 国金样本中位数约 25 条 | `sourceMin` 提到 20（品牌）/ 25（创业），`sourceOfficialRatioMin` 保持 0.4 |
| L6 | 低 | `agents/judge-*.md` 第 12 行 | 三位评委都“逐章读 `05_正文/*.md`”，附录（问卷全文、CSV 明细）也会被读进上下文 | 改为“正文章节逐章读；附录只 `Glob` 核对存在并读前 30 行” |
| L7 | 低 | `CLAUDE.md` Step 8 | 交付摘要“最终得分与官方口径模拟分”没有要求写“剩余问题”与“轮次历史” | 交付摘要固定 6 段：文件清单 / 各轮分数表 / 最终官方口径估分 / 剩余问题（来自 `评审汇总.json.remaining`）/ 示例数据清单条数 / 团队待办 |
| L8 | 低 | `agent/.claude/settings.json` | `permissions.allow` 含 `Skill` 与 `Agent`，但无 `Bash` 细粒度限制，子代理（不受 bypassPermissions 影响）拿到的是全量 Bash；与架构“只读评委”目标一致性靠 agent `tools` 字段保证，可接受 | 可加 `"deny": ["Bash(rm -rf:*)", "Bash(git:*)"]` 作保险 |

### 2.2 高危项的完整替换文本

#### 2.2-H1 评审与修改轮（替换 `agent/CLAUDE.md` 「### Step 7 评审（review）」整节）

~~~markdown
### Step 7 评审（review）

> 一轮 = 质检 → 三评委 → 汇总 → 判定 →（未达标）修改 → 重渲染。轮数上限 `reviewRounds`（默认 3）。主代理在本步**只做调度、落盘、汇总与判定**：不亲自改稿、不替评委打分、不改评委给出的分数。

动作（第 N 轮，N 从 1 起）：
1. **质检**：`mkdir -p 07_评审/第N轮`；写 `任务/质检_第N轮.md`（轮次、赛道、作品 docx / pdf 路径、上一轮修改清单路径），派发 **质检员**（`qa-inspector`）。回收后确认 `07_评审/第N轮/机器检查.json` 存在且含 `summary`；若 `bizkit.checks` 不可用，按规则 14 用 `tmp/checks_lite.py` 生成同名文件并置 `summary.degraded = true`。
2. **评委任务书**：按本文末「附 A 评委任务书模板」写三份 `任务/评委_investor_第N轮.md`、`任务/评委_professor_第N轮.md`、`任务/评委_executive_第N轮.md`。三份只有 `judge` 不同。
3. **并行派发三位评委**（`judge-investor` / `judge-professor` / `judge-executive`）：在同一条回复里连发三个 `Agent` 调用，prompt 固定为“读 `任务/评委_<judge>_第N轮.md` 并完成评审；回复只包含一个 ```json 块和一个 ```markdown 块，不要有其他文字。”按 §四「异步回收」等待三份结果。
4. **落盘**（评委没有 Write 工具）：对每位评委的回复——
   - 取第一个 ```json 块写入 `07_评审/第N轮/评审_<judge>.json`，取第一个 ```markdown 块写入 `07_评审/第N轮/评审_<judge>.md`（用 Write 原样写，不改分数、不改措辞）；
   - 写入前补齐 / 覆盖四个字段：`round = N`、`judge = "<judge>"`、`track = "<赛道>"`、`model`（取 `project.json.models.<judge>`；无此字段填 `"n/a"`）；
   - 校验：`"${BIZGOLD_PYTHON:-python}" -c "import json,sys;d=json.load(open(sys.argv[1],encoding='utf-8'));assert len(d['scores'])==int(sys.argv[2]),'维度数不符';assert 0<=d['total']<=100;print('ok')" 07_评审/第N轮/评审_<judge>.json <rubric 维度数>`；
   - 解析失败：先修复明显格式问题（尾逗号、中文引号、`//` 注释、未闭合括号）再存；仍失败或维度数不符 → **重派该评委一次**（同一任务书，prompt 追加“上次输出的 JSON 无法解析：<错误原文>，请严格重新输出”）；第二次仍失败则该评委本轮缺席：`评审汇总.json.byJudge.<judge> = null`，在《任务解析》异常记录登记。
5. **汇总** `07_评审/第N轮/评审汇总.json`（契约 §7.6，字段如下）：
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
   - 按 `owner` 拆任务书：`任务/修改_第N轮_主笔.md`（writer 项，> 6 条时按章拆成两份 `_主笔1 / _主笔2` 并行）、`任务/修改_第N轮_路演.md`（pitch 项）、`任务/修改_第N轮_财务.md`（finance 项）。每份按「附 B 修改任务书模板」写：条目表 + 可改文件白名单 + “只改清单项，不得重写整章 / 整份”。分别派发 `chief-writer` / `pitch-designer` / `financial-modeler`。`owner = main` 的条目（book.json 封面字段、charts.json 单图数据、参考文献条目、匿名词替换）由主代理直接改。
   - 执行顺序与重渲染（有财务项时从第一步起，否则从受影响步起）：财务 → `bizkit.finmodel` → 受影响图表 `bizkit.charts 04_图表/charts.json --only <ids>` → 正文 → `bizkit.docx_render 05_正文/book.json` → `bizkit.pptx_render 06_路演/slides.json`。任何一步失败按规则 12 处理。
   - 回收核对：逐条对照 `accept` 验收；未完成项带原文重派一次；仍未完成的保留到下一轮清单并标注“第 N 轮未完成”。
   - N = N + 1，回到动作 1。

验收：每轮目录含 `机器检查.json|md`、`质检整改清单.md`、三份 `评审_<judge>.json|md`（缺席者有异常登记）、`评审汇总.json`、`修改清单.md`；`project.json.review` 由 `update_progress(round=N)` 同步；最后一轮 `评审汇总.json.remaining` 与交付摘要一致。
~~~

**附 A 评委任务书模板**（追加到 CLAUDE.md 文末，或放 `tracks/_templates/评委任务书.md`；`{…}` 由主代理填写）

~~~markdown
# 评审任务：{judge} 第 {N} 轮

## 固定值（原样写入 JSON）
round = {N}；judge = "{judge}"；track = "{赛道}"；rubric = "tracks/{赛道}/rubric.json"；model = "n/a"（主代理落盘时覆盖）

## 必读（按顺序；只读这些，不要读 knowledge/、tracks/ 其他文件）
1. `tracks/{赛道}/rubric.json`（维度、锚定、machineCapRules、hardGates、officialWeight）
2. `07_评审/第{N}轮/机器检查.json`（只看 `summary` 与 `results` 中 fail 项；缺失则跳过 cap）
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
~~~

**附 B 修改任务书模板**（`任务/修改_第N轮_<角色>.md`）

~~~markdown
# 修改任务：第 {N} 轮 — {角色}

## 条目（只改这些；按序号逐条完成并在回报中逐条确认）
| 序号 | 文件与位置 | 现状（评委 / 机器原话） | 改法 | 验收（主代理将按此核对） | 预期提升 |
|---|---|---|---|---|---|
| {id} | {file}#{where} | {现状} | {what} | {accept} | {gain} |

## 可改文件白名单
{文件列表}。白名单之外的文件一律不动；需要新增数字先查 `03_财务/numbers.json`，账本没有的数字不能写（可改为定性表述并标【待核实】）。

## 规则
- 只改清单项，不重写整章 / 整份；保留原有结论句、图表引用与文献号。
- 新增引用写 `05_正文/新增文献_修改{N}.json`（临时号 `R{N}-1…`），正文用临时号引用。
- 改完每个文件运行对应检查（正文：`bizkit.checks --chapter`；路演：`bizkit.pptx_render`；财务：`bizkit.finmodel` 检查页为 0）。
- 回报 ≤ 300 字：逐条“序号 → 已完成 / 未完成（原因）”，列出改动文件。
~~~

**三份 `agents/judge-*.md` 「输出」段最后一句替换为**

~~~markdown
你没有 Write 工具。最终回复**只包含**两个代码块：先一个 ```json 块（契约 §7.6 完整结构，`round / judge / track / rubric` 按任务书固定值填写，`model` 填 "n/a"），再一个 ```markdown 块（人读版）。代码块外不得有任何文字；JSON 中不得有注释、尾逗号或中文引号。主代理负责落盘并覆盖 `model` 字段。
~~~

**`docs/contracts.md` §7.6 汇总行替换为**

~~~text
汇总 `07_评审/第N轮/评审汇总.json`：`{ round, score（三评委中位数，缺席 1 人时取均值）, mean, byJudge, byDimension, minDimension, machine: {…}, hardGatePassed, officialEstimate: {written, pitch, total}, fixes: [{ id, owner, file, where, what, accept, gain }], converged: bool, remaining: string[] }`。
~~~

#### 2.2-H2 环境自检与降级（新增 `agent/CLAUDE.md` 规则 14 + Step 0 动作 + §七新增小节）

**在 ⛔ 规则 13 之后追加规则 14：**

~~~markdown
14. **环境自检与降级**：Step 0 第一个动作是运行工具包自检（见 Step 0 动作 0），把结果写进 `00_任务解析/环境自检.md`。缺失的模块按「七、7.1 降级表」处理：**降级是预案，不是失败**，不得因模块缺失停止流程或向用户提问；每次走降级路径都要在《任务解析》异常记录登记一行（模块 / 替代做法 / 对交付物的影响）。规则 9 的“一律用工具包 CLI”在该模块缺失时让位于 7.1 的替代做法。
~~~

**Step 0「动作」列表之前插入：**

~~~markdown
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
~~~

**§七 末尾新增小节：**

~~~markdown
### 7.1 工具包模块缺失时的降级表（规则 14）

| 缺失模块 | 影响 | 替代做法（唯一允许的例外） | 交付影响登记 |
|---|---|---|---|
| `bizkit.finmodel` | Step 3 无法生成 xlsx / numbers.json / 财务摘要 | 建模师仍写 `assumptions.json`（契约 §7.4）；然后写 `tmp/finmodel_lite.py`（允许使用 pandas / openpyxl）：按假设逐年计算收入、成本、利润、现金流（startup：三表 + 盈亏平衡月 + 现金最低点 + 出让比例；budget：预算 = Σ 单价 × 数量 + 机动金 + KPI 表），把结果写成 `03_财务/财务模型.xlsx`（每张表第一行写公式说明）、`03_财务/numbers.json`（契约 §7.5 固定 key）、`03_财务/财务摘要.md`（含“检查”表：负现金未融资 / 负库存 / 出让比例误差 三项结果）与 `03_财务/财务图表.charts.json`。脚本运行 ≤ 5 分钟 | “财务模型为脚本生成，Excel 内无活公式；提交前请用模板重建或人工复核” |
| `bizkit.pptx_render` | Step 6 无 pptx | 仍完整产出 `slides.json`、`讲稿.md`、`答辩题库.md`、`视频分镜.md`；再写 `tmp/pptx_lite.py`（允许 python-pptx）按 slides.json 生成“标题 + 要点 + 图片 + 备注”的 16:9 素版 `产出/路演.pptx`（正文 ≥ 24 pt，≤ 3 色，母版不含任何院校信息，备注写 `note`）；不做复杂版式 | “PPT 为素版，视觉需团队用模板二次排版；内容与数字以 slides.json 为准” |
| `bizkit.checks` | Step 5 逐章检查、Step 7 机器检查、硬门槛判定全部缺失 | 主代理写一次 `tmp/checks_lite.py`（正则 + 计数，≤ 150 行）并在 Step 5/7 复用，输出与契约 §7.7 同结构的 `机器检查.json|md`（`summary.degraded = true`）。最少覆盖：`anonymity_hits`（用 checks.json 的 anonymityRegex + whitelist 扫 05_正文、06_路演、03_财务、book.json、slides.json）、`forbidden_terms` / `forbidden_terms_adlaw`（forbiddenTerms + advertisingLawTerms − whitelist）、`numeric_density`（每千字数字个数）、`chart_count`（`04_图表/index.json` 条数与正文 `![` 引用数）、`source_count`（参考文献.json 条数、独立域名数）、`struct_coverage`（requiredElements 名称或别名是否出现在一级标题）、`page_limit`（有 pdf 用 `bizkit.pdf --pages`，否则按 650 字/页估算）、`cross_doc_consistency`（numbers.json 每个 value 在正文 / slides.json / 讲稿中出现的写法是否一致：同一 key 的数字出现两种以上数值即 fail）。覆盖不到的检查项写 `"status": "skipped"`，评委对 skipped 项不施加 cap | “机器检查为简化版，`source_pairing / url_alive / font_check` 未执行，提交前请人工核对” |
| `bizkit.search` | 无命令行检索备用 | 只用 MCP `mcp__bizgold__web_search / web_fetch`；MCP 也不可用时研究员只用已知权威出版物并全部标【待核实】，调研总览首行注明“本次调研未联网” | “参考文献未经在线核实” |
| `bizkit.pdf` 转换失败 | 无 PDF、无精确页数 | 保留 docx；页数用 docx_render 输出的估计值；交付摘要提示“用 Word 另存为 PDF 后核对页数 ≤ 上限” | 同左 |
~~~

#### 2.2-H3 撰写分批与页数预算（替换 `agent/CLAUDE.md` 「### Step 5 撰写（write）」整节；Step 0 新增第 7 部分；替换 `chief-writer.md` 「收尾」段）

**Step 0「动作」第 6 部分之后追加：**

~~~markdown
7. **页数预算表**：读 `tracks/<赛道>/outline.json` 的 `chapters[].pages` 与 `appendices`，按 `pageTarget` 等比缩放：章节合计 = pageTarget × 0.85（取整到 0.5 页），附录合计 = pageTarget × 0.15 且 ≤ spec 附录上限；每章目标页 = 该章区间中值 ÷ 全部中值之和 × 章节合计，不得低于区间下限的 80%。表字段：文件名 / 目标页 / 目标字数（650 字/页，减去图表占页：图 0.45 页、表 40 行 1 页）/ 必配图表数。每个一级标题另起一页，预算已含该空白。该表是 Step 5 任务书与超页裁剪的唯一依据。
~~~

**Step 5 整节：**

~~~markdown
### Step 5 撰写（write）

前置：`00_任务解析/任务解析.md` 第 7 部分「页数预算表」存在；`04_图表/index.json` 与 `03_财务/numbers.json` 就绪。

动作：
1. **分批**：章节清单 = `outline.json.chapters`（含 `00_概要`）+ `appendices`。去掉 `00_概要` 后按顺序每 3 章一批，附录单独一批，`00_概要` 单独放在最后一批（它要引用各章结论与数字）；批数 = ⌈(章数 − 1) ÷ 3⌉ + 2。**同时最多 2 批并行**：先派第 1、2 批，全部回收并核对后再派第 3、4 批，依此类推。
2. 每批写任务书 `任务/主笔_<批号>.md`（六段格式，≤ 1500 字），必须包含：本批章节的文件名与一级标题（照抄 outline）、每章目标页与目标字数（来自页数预算表）、每章结论句（来自 `02_策略/策略框架.md`）、每章必答（来自 `outline.mustAnswer`）、本章可用图表 id 与标题（来自 `04_图表/index.json`，按章分配）、可用数字 key（来自 numbers.json）、本批禁止触碰的文件（其他章、`book.json`、`01_调研/参考文献.json`）。prompt：“读 `任务/主笔_<批号>.md` 并完成；只写本批章节文件；完成后返回 ≤ 300 字摘要。”
3. 主笔每章写完运行 `bizkit.checks --chapter`（模块缺失时跳过，由第 5 步统一检查）。主笔新增引用一律写到 `05_正文/新增文献_<批号>.json`（字段同参考文献.json，`id` 用 `W<批号>-<序号>` 临时号），正文暂用 `[W1-3]` 这类临时号引用。
4. **回收核对**（每批，按 §四）：文件齐全且标题与 outline 一致；每章首段为加粗结论句；每章 ≥ 1 图 1 表；章末有「本章小结」；示例数据段带标记；每章字数在目标字数 ± 15%；无 `[待补]` / `TODO`；必答问题逐条有回应（摘要中列出未答项则补调研或降级为【待核实】）。不合格带具体问题重派该批一次。
5. **合并与渲染（主代理）**：
   - 合并全部 `05_正文/新增文献_*.json` 进 `01_调研/参考文献.json`（按 URL 去重、顺延编号），用 Python 把 `05_正文/*.md` 中的临时号替换为正式号，删除 `新增文献_*.json`；
   - 写 `05_正文/book.json`（契约 §7.2；封面字段来自任务书；`anonymize.forbidden` 填任务书中出现的院校词；`pdf: true`）；
   - 运行 `"${BIZGOLD_PYTHON:-python}" -m bizkit.docx_render 05_正文/book.json`，读输出 JSON：匿名命中 → 修复后重渲染；记录页数（有 PDF 用 `bizkit.pdf --pages`，否则用估计值）。
6. **页数控制**：目标 = pageTarget ± 10%，上限 `checks.json.pageLimit` 不可超。超页：对照页数预算表找超支最多的章，先删重复图、再合并表、再压缩附录（问卷缩小字号、CSV 明细移出正文），每次调整后重渲染；不足：补读图句、与基准对照、敏感性分析，不放大图片凑页。调整 ≤ 2 轮，仍不达标登记异常并继续。

验收：全部章节与附录 md 存在；`grep -rn "\[W[0-9]" 05_正文` 为空（无临时号残留）；docx 生成；页数在目标区间且 ≤ 上限；章检查无 fail（或已登记降级）。
~~~

**`agents/chief-writer.md` 「## 收尾（最后一批负责）」整段替换为：**

~~~markdown
## 收尾
不写 `book.json`、不运行 `docx_render`、不改 `01_调研/参考文献.json`（主代理在全部批次回收后统一完成）。本批新增引用写入 `05_正文/新增文献_<批号>.json`（字段同参考文献.json，`id` 用 `W<批号>-<序号>`），正文中用 `[W<批号>-<序号>]` 临时引用。
完成后返回 ≤ 300 字摘要：章节文件清单、每章字数与估页、图表引用数、新增文献数、自检结果、未能回答的必答问题。
~~~

#### 2.2-H4 参考文献并发写入（替换 `agents/industry-researcher.md` 第 20 行 + `CLAUDE.md` Step 1 第 3 点）

**`industry-researcher.md` 第 20 行（“追加 `01_调研/参考文献.json`…”）替换为：**

~~~markdown
- 写 `01_调研/参考文献_<主题字母>.json`（**只写自己这一个文件，不读不写 `参考文献.json`**；JSON 数组；`id` 用 `<主题字母><序号>`，如 `A1`、`A2`）：`{ "id": "A1", "type": "report|web|paper|book|annual-report|gov|news", "author": "机构或作者", "title": "标题", "publisher": "发布方", "year": 2025, "url": "…", "accessed": "YYYY-MM-DD", "fetched": true }`（`fetched` = 是否用 web_fetch 成功打开过原文）。调研笔记与 CSV 的 `source` 列用同一临时号引用（如 `[A1]`）；正式编号由主代理统一分配。
~~~

（`consumer-insight.md` 若追加文献，同样改为 `参考文献_消费者.json`，临时号 `C1…`。）

**`CLAUDE.md` Step 1 第 3 点替换为：**

~~~markdown
3. 回收与合并（全部研究员与消费者洞察回收并通过 §四核对之后）：
   a. **文献合并**：读 `01_调研/参考文献_A.json … _D.json`、`参考文献_消费者.json`，按 URL 去重（同 URL 保留 `fetched: true` 者），按 A → B → C → D → 消费者的顺序顺延编号写 `01_调研/参考文献.json`，并写 `01_调研/文献编号映射.json`（`{"A1": 1, "A2": 2, "B1": 9, …}`）；用一段 Python 把 `01_调研/*.md` 与 `01_调研/数据/*.csv` 中的临时号替换为正式号（`[A1]` → `[1]`）。之后任何环节新增文献只能走「`新增文献_*.json` → 主代理合并」，不得直接编辑 `参考文献.json`。
   b. 写 `01_调研/00_调研总览.md`——关键事实 20 条（每条四要素 + 正式文献号）、数据表索引（CSV 文件名 → 字段 → 用途）、信息空白与【待核实】清单、对策略的 5 条启示。
   c. 核对 `参考文献.json`：条数 ≥ `checks.json.sourceMin`、独立域名 ≥ 8、官方 / 年报 / 协会来源占比 ≥ 40%、近 3 年 ≥ 70%、`fetched: true` 占比 ≥ 60%；不达标补派 1 名研究员（任务书只列缺口）。
   d. **关键数字汇总**：把 `01_调研/数据/关键数字_*.json` 合并为 `01_调研/数据/关键数字_汇总.json`（同 key 冲突取官方 / 更近年份，并在总览注明）；**此时不记账**，记账在 Step 3 验收后执行（见 Step 3）。
~~~

#### 2.2-H5 数字账本的写入顺序（替换 `CLAUDE.md` Step 3 「产出」与「验收」两行）

~~~markdown
- 产出 `numbers.json`（固定 key）、`财务模型.xlsx`、`财务摘要.md`、`建模说明.md`、`财务图表.charts.json`。
- **记账合并（主代理，建模师回收并核对后立即做，顺序不可颠倒）**：调用 `mcp__bizgold__record_numbers(entries = 01_调研/数据/关键数字_汇总.json 的全部条目)`。条目 key 若不在契约 §7.5 固定 key 集合内，统一加前缀 `research_`（如 `research_price_median`）；示例数据保留 `synthetic: true`；`source` 必须是 CSV 路径或带文献号的调研笔记位置。然后重新 `Read 03_财务/numbers.json`，确认同时含财务 key 与 `sample_size`、`tam`（budget 模式为 `budget_total`）。此后任何人不得直接编辑 `numbers.json`：财务数字改动只能改 `assumptions.json` 重跑 finmodel，调研数字只能通过 `record_numbers` 追加。

验收：`numbers.json` 覆盖数字需求单全部 key 且含调研 key；检查页为 0；预算 ≤ 上限；`财务摘要.md` 抽 5 个数字与 `numbers.json` 一致。
~~~

#### 2.2-H6 子代理调度与异步回收（替换 `CLAUDE.md` 「## 四、子代理调度规范」整节；规则 11 末尾加“回收规则见第四章”）

~~~markdown
## 四、子代理调度规范

- **任务书**：路径 `任务/<角色>_<编号>.md`，编号规则固定——研究员 `研究员_A…D`、消费者洞察 `消费者洞察_1`、建模师 `建模师_1`、主笔 `主笔_<批号>`、路演 `路演_1`、质检 `质检_第N轮`、评委 `评委_<judge>_第N轮`、修改 `修改_第N轮_<角色>`。六段：背景 / 目标问题（编号）/ 输入文件（**只列必需文件，并写“不要读其他文件”**）/ 输出文件与结构 / 约束 / 回报格式。≤ 1500 字；子代理需要的事实（对象、命题目标、结论句、数字 key）直接抄进任务书，不让子代理再翻全库。
- **调用**：`Agent` 的 prompt 固定一句：“读 `任务/xxx.md` 并完成；只写规定的输出文件；完成后返回 ≤ 300 字摘要（含输出文件路径、数据条数、待核实项）”。需要并行时在**同一条回复**里连发多个 `Agent` 调用（研究员 ≤ 4、主笔 ≤ 2、评委 = 3）。
- **异步回收**：`Agent` 可能立即返回 `status: "async_launched"`（后台运行），真正结果随后以 `task_notification` 送达。规则：① 收到 `async_launched` 后不要重复派发同一任务；② 等到该任务的通知到达（或 `task_progress` 显示完成）再做回收核对；③ 等待期间可以做不依赖其结果的工作（写下一份任务书、整理进度清单），但**不得进入下一个 Step**；④ 超过 30 分钟无通知 → `ls -la` 查输出文件，已生成且非空则按文件回收，否则重派一次。
- **回收核对清单（每个子代理必做）**：输出文件存在且非空；文件内无 `TODO` / `[待补]` / `XXX`；每条数据带来源或临时文献号；示例数据带标记；摘要所列文件与实际一致；本角色的专项项（研究员：CSV 最后一列为 `source`；主笔：首段加粗结论句 + 章末小结；建模师：检查页为 0；路演：渲染警告 = 0；评委：JSON 可解析且维度数正确）。任一不满足 → 带具体问题重派同一角色一次（prompt 追加“上次缺：…”）；第二次仍不满足 → 降级（缩小范围）并登记异常。主代理**不得**自己“顺手补数”或替子代理编造结果。
- 研究员之间不重复主题；两源数据不一致时取官方 / 更近年份，并在总览注明。
- 评委只读：不得修改任何文件；三位并行；主代理不得替评委打分、不得修改评委给出的分数与评语。
- 子代理默认会加载本 CLAUDE.md；评委与质检员的任务书必须自足（`omitClaudeMd: true` 时它们看不到本文件）。
~~~

#### 2.2-H7 案例赛道的匿名例外、契约与评分映射

**`CLAUDE.md` 规则 5 末尾追加一句：**

~~~markdown
例外：会计与商业管理案例赛道（`subject.caseVariant = accounting`）细则要求封面写选手与指导教师姓名——封面对应字段一律用占位 `【提交前由团队填写：选手姓名】`、`【提交前由团队填写：指导教师姓名】`，正文、页眉页脚、PPT、讲稿仍不得出现任何姓名；该两个占位串写入 `checks.json.anonymityWhitelist`，打包 README 提醒填写。
~~~

**`tracks/case/checks.json` 补丁（替换 / 新增字段）：**

~~~json
{
  "anonymityWhitelist": ["大学生", "University of", "科大讯飞", "【提交前由团队填写：选手姓名】", "【提交前由团队填写：指导教师姓名】"],
  "coverFieldsByVariant": {
    "accounting": ["团队名称", "组别", "【提交前由团队填写：选手姓名】", "【提交前由团队填写：指导教师姓名】"],
    "retail": ["团队名称", "作品名称", "案例企业名称"]
  },
  "requiredEvidenceByVariant": {
    "accounting": { "problem_identification": true, "industry_cluster_analysis": true, "ratio_table_with_formula": true, "annual_report_page_refs": true, "three_statement_check": true },
    "retail": { "regional_feature": true, "analysis_framework_diagram": true, "novel_method": true, "word_min": 3000, "deck_min_pages": 20 }
  }
}
~~~

说明：删除白名单中的 `中国科学技术大学`（真实院校名不得进白名单，命中后由主代理按规则 5 替换为“案例企业母校”之类的中性表述或直接删除）。

**`tracks/case/rubric.json` 新增字段：**

~~~json
"officialMapByVariant": {
  "accounting": { "total": 70, "items": [
    { "name": "内容完整程度", "weight": 5, "from": ["structure"] },
    { "name": "问题识别", "weight": 5, "from": ["environment"] },
    { "name": "案例企业财务指标分析", "weight": 10, "from": ["finanalysis"] },
    { "name": "行业数据聚类分析", "weight": 5, "from": ["finanalysis", "environment"] },
    { "name": "论述和决策确认的清晰程度", "weight": 10, "from": ["strategy"] },
    { "name": "结论的正确有效性", "weight": 10, "from": ["strategy", "forecast"] },
    { "name": "财务与会计知识的运用", "weight": 10, "from": ["finanalysis", "forecast"] },
    { "name": "对问题的切题回应度", "weight": 10, "from": ["implementation", "governance"] },
    { "name": "报告和阐释形式", "weight": 5, "from": ["structure"] } ],
    "pitchTotal": 30 },
  "retail": { "total": 60, "items": [
    { "name": "分析框架", "weight": 20, "from": ["strategy", "structure"] },
    { "name": "分析对象（含地区性特点）", "weight": 10, "from": ["environment"] },
    { "name": "分析视角", "weight": 10, "from": ["strategy"] },
    { "name": "分析方法", "weight": 10, "from": ["finanalysis"] },
    { "name": "结论", "weight": 10, "from": ["forecast", "implementation"] } ],
    "pitchTotal": 40 }
}
~~~

**`docs/contracts.md` §2 `subject` 补字段：** `"caseVariant": "accounting" | "retail"`（仅 case 赛道；服务端由命题库写入，任务书同步输出）。

**`tracks/case/spec.md` 流通业要素表追加一行：** “地区性特点｜必须有一节分析案例企业在参赛队所在省份的门店 / 供应链 / 消费者布局与当地政策，配 1 张区域数据表（来源：省统计局 / 企业年报分地区数据）｜评分项“分析对象 10 分”原文要求”。

#### 2.2-H8 上下文预算（替换 `CLAUDE.md` 「提速铁律」一行与 Step 0 「输入」一行；评委读法见 L6）

**「提速铁律」替换为：**

~~~markdown
> 提速与上下文铁律：① 任何文件只 `Read` 一次，需要回看先看自己的笔记（《任务解析》《调研总览》）；② > 300 行的文件先 `Grep` 定位再按 `offset/limit` 分段读所需部分；③ `knowledge/` 只按需定点读（见 Step 0），`knowledge/评审量规总表.md` 与 `tracks/<赛道>/rubric.json` 重复，不读；④ 已有结果直接复用；⑤ 每完成一步重新输出进度清单；⑥ 一次 `Bash` 不超过 10 分钟、单个脚本 ≤ 5 分钟，长任务拆开；⑦ 大段内容（策略框架、charts.json、book.json）写文件，不在对话里复述。
~~~

**Step 0「输入」行替换为：**

~~~markdown
输入（按此顺序与方式读，总量控制在 25 K 字以内）：
1. `任务书.md`（全文）。
2. `tracks/<赛道>/spec.md`（全文）与 `rubric.json`（全文）；`checks.json`、`naming.json`、`outline.json` 只在用到字段时读。
3. `knowledge/评分表与扣分细则.md`：先 `Grep -n "^### "` 列出小节，只读本赛道评分表小节 + “通用扣分细则”小节；其他赛道不读。
4. `knowledge/命题企业与获奖特征.md`：仅命题 / 指定企业赛道读，且先 `Grep` 企业名或命题关键词，只读命中的条目与“获奖标题范式”小节。
5. `inputs/_preread/*.txt`：每个文件先读前 120 行与 `Grep` 任务书中的关键词（品牌名 / 产品 / 区域 / 目标）命中行 ±20 行；> 600 行的文件写 `tmp/预读索引_<文件>.md`（每 100 行一句摘要）后只读索引。
6. `knowledge/国金写作标准.md`、`财务模型标准.md`、`市场测算与数据源.md`、`路演与答辩.md` 由对应子代理或对应 Step 按需读，Step 0 不读。
~~~

---

## ③ 你会补充的内容

### A. 新增 / 改造的 agents 与 skills

| 名称 | 类型 | 作用 | 要点 |
|---|---|---|---|
| `chart-maker` | agent（新） | 承接 Step 4：按主代理写的“图表需求单”读 CSV / numbers.json，写 `04_图表/charts.json` 并渲染、看图抽查 | `tools: Read, Write, Edit, Bash, Grep, Glob`，`skills: [biz-charts]`，maxTurns 50；输出 index.json 核对 ≥ chartMin；解决 M18 |
| `creative-ideation` | skill（新） | Big Idea 的 10 → 3 → 1 发散收敛，含 7 项检查清单评分表、与竞品传播主题的差异矩阵、视觉锤 / 语言钉生成规则、商标初筛（中国商标网类别检索提示） | 解决 M15；被 Step 2 与 `brand-creative` 引用 |
| `revision-round` | skill（新） | 修改轮操作规程：读清单 → 定位 → 最小改动 → 重检查 → 逐条回报；给出“只改不重写”的 Edit 用法示例与禁止事项 | 配合 H1 的修改任务书，供 `chief-writer / pitch-designer / financial-modeler` 共用 |
| `case-analysis` | skill（新） | 会计案例：年报页码引用、比率表（公式列）、杜邦分解、现金流质量、行业聚类（同业 5–8 家 2–3 个指标分组）、问题识别树；流通业：地区性分析框架与前沿方法清单（RFM、DEA、GIS 热力、情景规划） | 解决 H7 的质量面 |
| `mice-program` | skill（新） | 会奖：主题 → 日程网格（半天粒度）→ 供应商与单价来源 → 200–400 万预算表 → 应急预案 → 宣传视频四平台脚本 | 当前 mice 只有 spec 无技能支撑 |
| `digital-video` | skill（新） | 数字营销：短视频脚本（预选赛 40 分）、投放矩阵、数据回传与归因、A/B 方案 | 当前 digital 只有 spec |
| `ai-disclosure` | skill（新） | 生成 `08_提交包/团队确认清单.md` 与内部《AI 辅助使用说明》（不进提交物），说明哪些内容必须由团队重写 / 替换 / 核实 | 解决 M14、⑤合规 |
| `survey-design` | skill（改） | 加质量门槛表（α、KMO、有效率、访谈数、实地）、进阶分析二选一、示例结果必须报告三项指标 | 解决 M10 |
| `pptx-deck` | skill（改） | 字号 / 色数 / 每页信息量改为国赛投影标准；加“示例数据标注传播”规则 | 解决 M9、M14 |

### B. 新增的 CLAUDE.md 护栏（借鉴数模 CLAUDE.md 的颗粒度）

1. **反常结果先当 bug 查**：盈亏平衡 < 6 个月、毛利率 > 80%、SOM > SAM 10%、样本有效率 100%、ROI > 10 倍——先查假设与脚本，不得直接写进正文。
2. **禁止重复读取**：同一文件第二次 `Read` 前必须说明原因；评委 / 主笔任务书里直接给结论句与数字，减少回读。
3. **单脚本 ≤ 5 分钟、单 Bash ≤ 10 分钟**；图表渲染分批 `--only`。
4. **提交前自检表**（Step 8 前，主代理逐项打勾写入交付摘要）：匿名 0 命中 / 页数 ≤ 上限 / 封面字段 = naming.json / 命名规则 / PPT 页数区间 / 讲稿字数区间 / 题库 ≥ 40 / 示例数据条数 = 数据真实性清单条数 / 参考文献无临时号 / numbers.json 与正文抽查 5 个一致 / zip 内文件数 = README 清单数。
5. **进度清单重输出**时附一行“上下文自评：已读文件 N 个，建议下一步派子代理 / 写文件而非对话输出”。
6. **AI 使用说明三处一致**：交付摘要、`团队确认清单.md`、README 对“示例数据 / 须团队重写”口径一致。

### C. 新增的机器检查项（给 `bizkit.checks` 的需求，同时写入各 `checks.json`）

| check id | 口径 | 阈值 | 用途 |
|---|---|---|---|
| `source_pairing` | “图 / 表 N”编号后 3 行内出现“资料来源 / 数据来源” | ≥ 95% pass / < 80% fail | 02 §6.7 |
| `url_alive` | 参考文献 URL HEAD/GET 2xx–3xx | ≥ 90% / < 75% fail | 02 §6.7；研究员 `fetched` 字段可先替代 |
| `synthetic_marker_coverage` | `numbers.json` 中 `synthetic: true` 的 value 在正文 / slides / 讲稿 / 题库每次出现处 ±3 行内有“示例数据”字样 | 100% | ⑤合规 |
| `brief_response_table` | 引言含“命题回应表”且行数 ≥ 命题要求条数 | pass/fail | M16 |
| `honor_ratio` | 含“荣誉 / 获奖 / 认证 / 称号”的句子字数 ÷ 全文 | ≤ 5% | 规则 2.2 红线 |
| `title_pattern` | 标题符合 `outline.titlePattern`（破折号、字数区间） | pass/fail | spec 第 7 点 |
| `weak_terms` | `weakTerms` 词表每百字频次 | ≤ 0.1 | M4 |
| `deck_font_min` / `deck_color_count` | PPT 正文最小字号 / 颜色数 | ≥ 24 pt / ≤ 3 | M9 |
| `script_sentence_len` | 讲稿句长中位数 / 最长句 | ≤ 25 / ≤ 40 字 | 2.3 |
| `regional_feature` | 流通业变体正文含“所在省 / 本地 / 区域”章节且配区域数据表 | pass/fail | H7 |
| `similarity_template` | 与 `knowledge` 中“网传模板句”词表的 8-gram 重合率 | < 2% | 查重（简化版） |
| `appendix_questionnaire` | 附录含问卷全文且题数 ≥ 25 | pass/fail | rubric `structure` 5 分锚定 |

### D. 知识库补件

- `knowledge/单价库.md`：KOL 分层刊例区间、线下活动 / 物料 / 投放单价区间，每条带来源与年份（budget 模式“单价有来源”的现成依据）。
- `knowledge/获奖标题库.md`：近 4 届一等奖标题 30 条 + 结构拆解（已在 `命题企业与获奖特征.md` 中的抽出来单列，便于 Grep）。
- `knowledge/竞品传播主题库.md`：常见品类头部品牌近两年传播主题（避免 Big Idea 撞车）。

---

## ④ 模拟走查（按现行 CLAUDE.md，以 deepseek/kimi 级主代理视角逐步执行）

> 场景：品牌策划赛道，`researchMode = synthetic`，`pageTarget = 60`，`reviewRounds = 3`，`targetScore = 85`。记号：✅ 可无歧义执行；⚠ 有歧义但能猜；⛔ 卡死或必然出错。

| Step | 主代理动作（按文件顺序） | 歧义 / 卡点 | 预估主代理上下文增量 |
|---|---|---|---|
| 启动 | 自动注入 CLAUDE.md（26.5 KB）；读 `任务书.md`、`tracks/brand/spec.md`；`update_progress(parse, running)`；输出进度清单 | ✅ 规则 2 顺序清楚。⚠ 规则 7 要求“开始与结束都调用”，但进度清单格式的“评审（第 0/3 轮）”第一次输出时 round 填什么未说明 | ~25 K token |
| 0 读题 | 读 `rubric.json`、`评分表与扣分细则.md`（36 KB 全部赛道）、`命题企业与获奖特征.md`（28 KB）、全部 `inputs/_preread/*.txt`；写 `任务解析.md` 六部分 | ⛔ H8：知识文件全读后上下文已 60–100 K，部分模型在此触发压缩，后续忘记 CLAUDE.md 细则。⚠ “评分映射 100% 覆盖”没有给表头模板（官方项 / 章节 / 证据 / 目标分）→ 弱模型写成散文。⚠ 页数目标“≤ 上限 90%”与 pageTarget 由谁定冲突（任务书已给 60）。⛔ H2：没有任何环境自检，模块缺失到 Step 3 才发现 | +40–70 K |
| 1 调研 | 写 4 份 `任务/研究员_?.md`（编号 A–D 还是 1–4？§四说 `<n>`，契约 §8 说 `#2`）；同一回复派 4 个 `Agent` + 1 个 `consumer-insight`；等待；读 5 份摘要；写 `00_调研总览.md`；核对参考文献指标；`record_numbers` | ⛔ H6：`Agent` 返回 `async_launched`，CLAUDE.md 无说明 → 典型失败是主代理立即“回收”发现文件不存在、按规则 11 重派，派出 8–10 个研究员。⛔ H4：4 人并发追加 `参考文献.json`，id 重号。⚠ consumer-insight 输出文件 `0X_` 未定（M8）。⚠ “官方/年报/协会 ≥ 40%”需要主代理逐条判断 type，无字段映射（`type ∈ {gov, annual-report, report?}` 哪些算官方）。⛔ H5：此时写入 `numbers.json` 的调研数字会在 Step 3 被覆盖 | +15 K（读总览与摘要） |
| 2 策略 | 读 `00_调研总览.md`（可能 3–5 K 字）、`outline.json`；写 `02_策略/策略框架.md` 七部分 | ⚠ M15：Big Idea 无发散方法，产出常规创意。⚠ 第 5 点“各章结论句”需要先知道章节数 → 必须读 outline（未在输入列出）。⚠ 第 7 点标题范式对数字营销 / 会奖没有说明。✅ 其余可执行 | +10 K |
| 3 财务 | 写 `任务/建模师_1.md`；派 `financial-modeler`；回收读 `财务摘要.md`、`numbers.json` | ⛔ H2：`python -m bizkit.finmodel` → `No module named bizkit.finmodel`；建模师按规则 12 重试 3 次后不知如何降级（规则 9 禁止手写），很可能回复“无法完成”；主代理随后也无预案。⚠ budget 模式文件名 `预算表.xlsx` vs `财务模型.xlsx`（M6）。⚠ “建模说明追加到财务摘要末尾”重跑即丢（M6） | +8 K |
| 4 图表 | 读 `biz-charts/SKILL.md`、`财务图表.charts.json`、各 CSV；亲自写 `04_图表/charts.json`（≥ 20 图含数据数组）；运行 charts；`Read` 2 张 PNG | ⚠ M18：主代理必须把 CSV 数据搬进 JSON，20 张图 ≈ 400–600 行 JSON，上下文大增且易写错 `data` 结构（charts.py 校验失败 → exit 2 → 逐图修）。⚠ 品牌视觉“色板 / 字体 / 图形元素说明表”没有对应图型，模型可能试图用 `bmc` 乱画。✅ 其余清楚 | +20–30 K |
| 5 撰写 | 按“2 批 ≤ 3 章”写 `任务/主笔_1.md`、`_2.md`；派 2 个 `chief-writer`；回收；（谁写 book.json？）；渲染；读页数 | ⛔ H3：8 章 + 4 附录分不完；两批都认为自己不是“最后一批”或都写 book.json；两批同时改 `参考文献.json`。⚠ 概要若在第 1 批，引用不到后续章节数字。⚠ 页数 ± 10% 的调整“删重复图 / 合并表”没有定位依据（无页数预算表）。⛔ H2：`bizkit.checks --chapter` 不存在 → 主笔每章自检失败 3 次 → 浪费 maxTurns 80 的大半 | +15 K |
| 6 路演 | 写 `任务/路演_1.md`；派 `pitch-designer`；回收核对 pptx 警告、讲稿字数、题库维度 | ⛔ H2：`bizkit.pptx_render` 不存在 → 无 pptx。⚠ M5：CLAUDE.md 的页数 / 维度与 `checks.json.pitch` 不一致，主代理按 CLAUDE.md 验收会与机器检查冲突。⚠ M14：示例数据不带标记进入讲稿 | +6 K |
| 7 评审 | 派 `qa-inspector`（⛔ `bizkit.checks` 不存在 → 机器检查缺失，硬门槛无法判定）；写 3 份评委任务书（无模板，内容自拟）；并行派 3 评委；从回复中抠 ```json / ```markdown 落盘；写 `评审汇总.json` 与 `修改清单.md`；判定 | ⚠ 评委 JSON 的 `model` 字段评委不知道；`rubric` 字段含义不明（L2）。⚠ “每个维度 ≥ 3”是均值还是任一评委（H1）。⚠ 第 1 轮“提升 < 2 分”无参照。⛔ H1：未达标后“按清单修改”——主代理没有被告知派谁、任务书怎么写、改后重渲染顺序；最可能的行为是主代理自己 Edit 正文（上下文暴涨、改坏引用）或重新派 `chief-writer` 重写整章（分数震荡不收敛）。⚠ 3 评委用均值而官方是去高去低 | 每轮 +15–25 K（三份评委 JSON+MD 约 6–8 K 字） |
| 8 打包 | `bizkit.package`；核对 README；删 `tmp/`；输出进度清单与交付摘要 | ✅ 模块存在。⚠ 交付摘要无固定结构（L7），剩余问题可能被省略。⚠ “填写学校名”与命名规则占位靠 README 提醒，`团队确认清单` 不存在（M14） | +5 K |

**走查结论**：按现行文件，一次完整运行会在 Step 1（异步回收 / 文献竞态）出现第一次失真，在 Step 3 第一次卡死（finmodel 缺失），即使人工补齐工具包，也会在 Step 5（分批）与 Step 7（修改轮）再次停摆或不收敛；主代理在 Step 0 结束时的上下文占用已接近 deepseek/kimi 窗口的一半。修复 H1–H8 后，流程在纸面上可无歧义走通；质量上限的提升主要靠 M10 / M15 / M16 与 ③-A 的新技能。

---

## 附：五赛道对官方细则的覆盖核对（完整性维度速查）

| 赛道 | 官方评分 / 硬性要求 | 大脑覆盖 | 缺口 |
|---|---|---|---|
| 品牌策划 | 6×10 书面 + 8×5 展示；≤ 80 页；封面 3 字段；4 年查重；视频 3–5 人出镜 | spec / rubric / checks 齐全，映射清晰 | `fit` 维无机器检查（M16）；附录 A5/A6 与 outline 不一致（M12）；`最` 误报（M4） |
| 创业计划 | 创新 20 / 商业 30 / 团队 10 + 陈述 20 / 答辩 20；≥ 3000 字；PPT ≥ 20 页（命题组省赛扣分）；惠州命题方向 | rubric 10 维、`word_min` / `financing_consistency` 硬门槛、deck [20,24] | `officialMapping` 引用不存在的 `product`（M13）；团队章无“股权表 / 分工互补表”作为 requiredEvidence（建议加 `equity_table`、`team_matrix`——`equity_table` 已有，补 `team_matrix`） |
| 数字营销 | 预选赛短视频 40 + 方案 60；6×10 | spec 含短视频规则，rubric 6 维映射明确 | 无技能支撑（③-A `digital-video`）；视频评分无 rubric `pitchDimensions` 对应项 |
| 会计 / 流通业案例 | 会计 70（9 项）/ 30；≤ 50 页含附录；封面含姓名；流通业 60（含地区性 10）/ 40，≥ 3000 字、PPT ≥ 20 页 | spec 两变体、checks 有 `three_statement_check` 等 | H7 全部；`caseVariant` 未入契约；无 `case-analysis` 技能 |
| 会奖旅游 | 6×10 + 40；预算 200–400 万硬门槛；宣传视频四平台 | `budget_in_range` 硬门槛、spec 10 要素 | 无技能支撑（③-A `mice-program`）；日程网格 / 供应商表未列为 requiredEvidence（建议加 `itinerary_grid`、`vendor_quote_table`） |

（评审完）
