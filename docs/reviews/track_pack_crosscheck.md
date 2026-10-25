# 专题包总审报告（crosscheck）

- 范围：`agent/tracks/{brand,startup,digital,case,mice}/专题/00–06`（30 文件）、五个 `spec.md` §5 专题包段、五个赛道技能 `SKILL.md` frontmatter 与分工。
- 对照基准：`agent/CLAUDE.md`、`docs/contracts.md` §7、`bizkit/checks.py` check id、`bizkit/charts.py` 18 图型、`bizkit/pptx_render.py` 14 版式、各赛道 `checks.json / rubric.json / outline.json / naming.json / track.json`。
- 方法：`/tmp/bizgold_review3/cross/scan.py` 自动扫描（未知标识符、不存在的路径、无访问日期 URL、匿名正则、阈值与配置不一致、文件大小 / 行数）+ 逐文件人工通读。扫描输出 `cross/scan1.txt`（修前）、`cross/scan3.txt`（修后），完整 diff 见 `cross/final.diff`。
- 结果：改动 25 个文件（24 个专题 md + `tracks/mice/spec.md` 1 行），无新增 / 删除文件，未触碰 CLAUDE.md、agents、toolkit、rubric / checks / outline、其他技能。
- 验证：`cd agent/toolkit && .venv-bizgold/bin/python -m pytest -q tests` → **94 passed**；`tracks/*/*.json` 全部 `json.load` 通过；`git status` 仅 25 个 M、无 untracked；`git` 状态未做任何改变（无 add / commit / stash / checkout）。

## 一、跨包共性问题 → 修改

| # | 问题 | 修改（文件） |
|---|---|---|
| C1 | `chart` 版式只支持原生柱 / 条 / 堆叠 / 折线 / 饼环（`pptx_render.py`），各包 06 骨架却把 funnel / tornado / waterfall / scatter_matrix / gantt / radar / cashflow_combo / pie_as_bar 标成 `chart` → 渲染时会失败或降级 | 五个 `06_路演与答辩题库.md`：头注加入版式限制说明，对应页改为 `image（xx 图）`，原生图标明子类型（`chart（hbar）`、`chart（line / column_clustered）`、`chart（pie，≤ 4 类）`） |
| C2 | 06 页时间轴合计与「讲稿 9:15–9:30」规则不一致（startup 页表合计 10:40、digital 9:50、case / mice 10:00） | startup / digital / case / mice `06`：逐页重排到合计 9:30（brand 9:25），标题写明「合计 9:30 / 留余量」 |
| C3 | 单答时限五包不一（60 秒 / 90 秒 / 1 分钟 / 45 秒 / 仅 150 字） | 五包 `02` §演讲 + `06` §题库统一为「单答 ≤ 60 秒 / 150 字」 |
| C4 | 索引段落标题五包不一致（「文件用途（一行一个）」「文件一句话」「哪个阶段 / 哪个角色读哪一份」「谁在什么阶段读哪一份」），部分带无意义「大小」列 | 五个 `00_索引.md` 统一为 `## 文件一句话用途` / `## 谁在哪个阶段读哪一份` / `## 使用约定`，去掉大小列（startup、case 原缺 `## 使用约定`，已补） |
| C5 | 引用的细则 PDF 文件名被「…」截断，路径不存在（case 01/02/03/06、mice 01/02） | 替换为 `docs/research/raw/` 下真实文件名（`sdssfw_1449_【山东省赛细则】2026年全国高校商业精英挑战赛会计与商业管理案例竞赛山东省赛细则.pdf`、`sdssfw_1443_【山东省竞赛细则】2026年全国高校商业精英挑战赛商务会奖旅游策划竞赛实施细则.pdf`） |
| C6 | 单文件超出简报 3–8 KB（brand 06 8.8 KB、startup 02 / 06、digital 06、case 01 8.2–8.4 KB） | 压缩赘述（不删规则与数字），全部 ≤ 8 192 B；索引均 ≤ 60 行 |
| C7 | 个别 URL 无访问日期 | brand `01`（cnipa）、digital `05`（巨量算数）补「访问 2026-10-25」；其余扫描命中均已在同段 / 段首统一标注日期（误报） |

## 二、分包问题 → 修改

### brand（`agent/tracks/brand/专题/`）
| # | 问题 | 修改 |
|---|---|---|
| B1 | `00` / `03` 引用 gitignore 的运行目录 `workspace/p_20261025_010541_a2ey`，仓库内不存在 | `00_索引.md` 使用约定、`03` §4 标题 / 头注 ⑤ 改为「内部运行样本（不随仓库分发）」，只保留对照结论 |
| B2 | `02` §4 写 `rubric.json.machineCapRules`，但该顶层键为 null，cap 规则在 `dimensions[].machineCapRules`；且漏 `forbidden_terms` cap | `02` §4 标题改为 `rubric.json.dimensions[].machineCapRules`（check id 见 `checks.json.requiredEvidence`），补 `forbidden_terms → 规范 ≤ 3` |
| B3 | `06` 题库标题「20 题」实际只有 18 题，团队与资源维度为空（`qaDimensions` 要求六维各 ≥ 4） | 补 Q10、Q18（团队与资源），补充题方向注明该维度再补 ≥ 2 |
| B4 | `00` 对 06 的描述（20 页 / 18 题）与正文（19 页 + 备用 / 20 题）不符 | `00` 行改为「19 页 PPT 骨架 + 备用页…评委最常问 20 题」 |
| B5 | `05` 写 `tracks/brand/track.json.financeMode`（非法路径写法） | 改为 `tracks/brand/track.json` 的 `financeMode` |
| B6 | `06` Q4 示例绑定具体命题「电解质软糖」 | 改为通用表述「认知 vs 行为的差值」 |

### startup（`agent/tracks/startup/专题/`）
| # | 问题 | 修改 |
|---|---|---|
| S1 | `02` §1 / `03` 附录表引用「附录 A6 迭代创新说明」，但 `outline.json` 只有 A1–A5 | `02`：注明「`outline.json` 默认 A1–A5，需要时在 `book.json` 追加 `A6_迭代创新说明.md`」；`03` 附录行改为「A1–A5（+A6）」。未改 outline.json（非错误，A6 仅非首次参赛需要） |
| S2 | `02` §4 扣分规则来源写成 §3.17 且遗漏 PDF < 3 000 字项 | 改为「PPT < 20 页或省赛 PDF < 3 000 字扣 3 分（`naming.json.deckMinPages / anonymityPenalty`）…（`01_赛事情报.md` §3.17）」 |
| S3 | `02` §5「单答 > 90 秒」与其余文件 60 秒冲突 | 统一 60 秒（见 C3） |
| S4 | `06` 页表合计 10:40，p9 / p13 漏斗与 p15 / p16 图型用 `chart` | 重排到 9:30；p9 / p13 → `image（funnel 图）`，p15 → `chart（line / column_clustered）`，p16 → `chart（pie，≤ 4 类）` |
| S5 | `00` 缺「使用约定」，索引标题不统一 | 补 `## 使用约定`（本包只补赛道差异；外部 URL 仅用于核规则，不直接进参考文献） |
| S6（保留，报告说明） | `05` Step 2 写「不做 Big Idea 10→3→1」，与 CLAUDE.md Step 2 第 4 条字面冲突 | 未改：startup `checks.json.requiredEvidence` 无 `big_idea`，创业计划书用价值主张替代创意主张是赛道合理差异；建议下一轮在 CLAUDE.md Step 2 写明「品牌 / 数字赛道必做，其他赛道按专题 05」 |

### digital（`agent/tracks/digital/专题/`）
| # | 问题 | 修改 |
|---|---|---|
| D1 | `06` 页表合计 9:50，p14 甘特 / p4 / p15 / p17 版式标注不清 | 重排到 9:30；p14 → `chart（hbar）或 image`，甘特缩略 `image`；头注加版式限制 |
| D2 | `04` 一处把纯表格写成图 | 改为「表（不画图）」 |
| D3 | `00` 带「大小」列、标题不统一 | 统一（见 C4） |
| D4 | `02` §5 单答 90 秒 | 统一 60 秒 |

### case（`agent/tracks/case/专题/`）
| # | 问题 | 修改 |
|---|---|---|
| K1 | 细则 PDF 路径截断（01/02/03/06） | 见 C5 |
| K2 | `03` 引用 `sem.ysu.edu.cn` 页面（院校域名 + 同句出现院校语境） | 删除该 URL，保留结论 |
| K3 | `06` 骨架版式：p5 scatter_matrix、p7 waterfall、p8 cashflow_combo、p17 tornado 标 `chart` | 改 `image（xx 图）`；p4 `chart（line）`、p12 `table`、p16 `chart（grouped_bar）或 table`；时间轴重排到 9:30 |
| K4 | `02` 缺单答时限 | 补「单答 ≤ 60 秒」 |
| K5 | `00` 缺「使用约定」、标题后无空行 | 补 3 条使用约定（含：科大讯飞示例数字须重新核实并经 `record_numbers` 登记后才能进正文） |

### mice（`agent/tracks/mice/专题/` + `tracks/mice/spec.md`）
| # | 问题 | 修改 |
|---|---|---|
| M1（关键） | `02` §三、`05` #27 与 Step 8 行、`spec.md` §4 要求封面写占位 `【提交前由团队填写：参赛者姓名 / 指导教师姓名】`；但 mice `checks.json.anonymityRegex` 含「指导教师」且 `anonymityWhitelist` 无占位句（case 赛道才有），CLAUDE.md 匿名例外也只给 `caseVariant = accounting` → 照做必命中 anonymity、评审循环无法通过 | `02` §三封面行、`05` Step 8 行与 #27、`spec.md` §4 第 2 条统一改为：系统封面只含 `naming.json.coverFields` 四项，不生成姓名或占位句；参赛者 / 指导教师姓名由团队提交前手工补填（细则七（一）），README 提醒。**这是对 spec.md 的唯一改动**（明显与 checks.json / CLAUDE.md 冲突） |
| M2 | `02` §四把方案 PDF 命名写成宣传视频规则「院校名称-团队名称-队长姓名-组别」 | 改为 `naming.json.fileNameRule`（`{学校名称}-{团队名称}-{作品名称}`）与 `promoVideoNameRule`（`…-队长姓名-组别-宣传视频`）分开，占位项列明，补 PDF ≤ 30 MB |
| M3 | `06` 骨架：p6 radar、p13 gantt、p15 pie_as_bar、p16 tornado、p18 funnel 标 `chart`；页表合计 10:00 | 改 `image（xx 图）` / `chart（hbar）或 image`；重排到 9:30，头注加版式限制 |
| M4 | `02` §二「≤ 45 秒」、`06` 只写 150 字 | 统一「单答 ≤ 60 秒 / 150 字」 |
| M5 | 细则 PDF 路径截断（01/02） | 见 C5 |
| M6 | `00` 标题 / 大小列 | 统一（见 C4） |

## 三、接线核对（无需修改的结论）
- 五个 `00_索引.md` 的 Step → 角色 → 文件 均与 CLAUDE.md Step 0 规则一致：研究员读 `05` 调研段、策略读 `01`、建模师读 `05` 财务段（mice 另读 `01` §5、`04` 第 7 章行）、主笔读 `03` / `04`、路演读 `06`、评委 / 质检读 `02`。
- 五个 `spec.md` §5 专题包表与索引一致，未改动（mice 仅改 §4 一行，见 M1）。
- 五个赛道技能 frontmatter：`name` 与目录名一致、`description` 均为可触发的场景描述（133–209 字）；技能引用的专题路径均存在；分工明确为「技能 = 操作步骤 / 写法，专题 = 标准、评分、图表、自检」，无内容重复。case-analysis 技能第 11 条的封面占位写法与 CLAUDE.md 会计例外、case `checks.json` 白名单一致。
- 阈值核对（deckPages、scriptWords 1 650–1 950、qaMin 40、chartMin 20/15/18/12/15、sourceMin 20/25/18/12/12、pageLimit 80 / case 50、sourceDomainsMin、officialRatio、contingency 5–10%、预案 ≥ 8）全部与各赛道 `checks.json` / `rubric.json` 一致。
- 扫描残留项均为误报：`fig1_roadmap` 等为图 id，`budget_contingency` 为 finmodel budget 键，`risk_social` 为 startup rubric 维度，`content_5707286` 等为 URL 片段，「教授」为评委角色，「我们学校 / 老师说」为禁语引用，「40 页版 / 第 10–14 页」非页数上限。

## 四、保留项与风险（供下一轮决策）
1. **院校域名 URL**：brand `03`（`wxy.wfu.edu.cn`、`news.jiangnan.edu.cn`）、startup `02` / `03`（`cupk.edu.cn`、`hebust.edu.cn`、`cafuc.edu.cn`、`gdufe.edu.cn`）引用的是细则转载页与获奖经验帖，域名不等于文中出现院校名，机器匿名正则不命中，故保留；索引已注明「仅用于核规则，不直接进参考文献」。若希望零院校痕迹，下一轮可改为引用组委会原始公告或写作「公开经验帖」。
2. **mice 封面姓名**：细则要求封面写姓名，现按「系统不生成、团队手工补填」处理。若希望与 case 赛道一致（生成占位句），需同时改 mice `checks.json.anonymityWhitelist`（加占位句）与 CLAUDE.md 匿名例外条款，本轮未越权。
3. **startup A6 附录**：`outline.json` 未加 A6，非首次参赛团队需在 `book.json` 手工追加。
4. **示例数据耦合**：brand `03` §4 的「样本作品体检」与 case 包的科大讯飞数字来源于单次运行 / 单一年报，已在索引标注「须重新核实」；它们是包内最容易过期的内容。

## 五、成熟度评分（1–5）与下一轮建议
| 专题包 | 评分 | 理由 | 下一轮建议 |
|---|---|---|---|
| brand | 4 | 方法论与评分映射最完整，cap 速查与 check id 对齐；题库补齐 20 题 | 把 `03` §4 的样本体检换成脱敏的通用「国银→国金」差距表；为 `04` 的图 id 清单与 `outline.json` 章节号建立一一映射 |
| startup | 4 | 评分表逐项「国金 / 国银 / 失分」颗粒最细，06 骨架可直接落 slides.json | 决定 A6 是否进 outline.json；在 CLAUDE.md 明确 Big Idea 步骤的赛道适用范围；把 `02` §6 的融资一致容差（±1 pp）同步写进 `05` 自检表对应条 |
| digital | 4 | 推理链 → 图表库 → 自检闭环清晰；行业基准值全部标注「须换当期来源」 | `01` 来源段中的第三方报告（飞瓜、FACT+S）链接稳定性差，建议改为「报告名 + 发布方 + 年份」并保留镜像说明；补一条直播合规（第十九条）自检 |
| case | 3.5 | 可复算链路（比率 → 杜邦 → 综合收益表）扎实，页数 / 来源阈值与 50 页 pageLimit 对齐 | 科大讯飞示例数字需要每届刷新机制（建议放到 `docs/research/` 而非专题）；流通业变体仅一节，缺独立评分映射与图表行；`06` 题库按 accounting / retail 两变体拆分 |
| mice | 3.5 | 预算三层法、日程四连表、合规底线具体可执行；本轮修正了与 checks.json 冲突的封面姓名规则 | 决定封面姓名策略（见四.2）并同步 checks / CLAUDE.md；`01` 肇庆场地示例数据改为「字段模板 + 核实步骤」；补「人群管理方案」互借的具体页 / 表引用 |

通用建议：① 给 `scan.py` 这类一致性扫描固化为 `toolkit/tests` 的一个只读测试（检查专题包引用的路径、check id、图型、版式存在），防止下次漂移；② 五包 `06` 的页时间轴可改为「秒数」单列由脚本求和，避免手工累加错误；③ 专题包的外部 URL 建议集中到每文件末尾「来源」段并统一日期，便于下一届整体刷新。
