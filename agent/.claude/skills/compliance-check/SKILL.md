---
name: compliance-check
description: 用 bizkit.checks 做机器检查与合规自检：结构覆盖、匿名、页数字数字体、数值密度、来源、图表、禁语与广告法、套话密度、示例数据标注、无源关键数字、财务勾稽、融资一致性、PPT 指标、讲稿时长、题库、跨材料数字一致性；阅读检查结果并修复。每章写完、评审轮次、打包前使用。
---

# compliance-check：评委一眼能抓到的硬伤，机器先抓

## 命令
```bash
export PYTHONPATH=toolkit PYTHONUTF8=1
"${BIZGOLD_PYTHON:-python}" -m bizkit.checks --track <赛道> --workspace . --chapter 05_正文/02_xxx.md   # 单章快检
"${BIZGOLD_PYTHON:-python}" -m bizkit.checks --track <赛道> --workspace . --round N                  # 全量 → 07_评审/第N轮/机器检查.json|md
```
阈值来自 `tracks/<赛道>/checks.json`（缺省用 `bizkit/track_defaults.py`）。

## 结果怎么读
- `hardGates`：任一 `passed:false` = 不可提交。匿名（`anonymity_hits`，正则 + 白名单）、页数上限、封面字段、跨材料一致性、广告法用词、（创业计划）字数 ≥ 3000、（会奖）预算区间。
- `checks[]`：`fail` 必修；`warn` 尽量修；`skipped` 表示文件缺失（如尚未生成 PPT）。
- `crossDoc.mismatches`：列出账本 key 在正文 / PPT / 讲稿中的不同数值——以 `numbers.json` 为准修改文档；若账本错了，先改 assumptions 重跑 finmodel 再改文档。

## 修复指南
| 项 | 修法 |
|---|---|
| struct_coverage < 100% | 对照 `checks.json.requiredElements` 补章节或改标题为标准名（可用别名） |
| anonymity_hits | 替换为团队名 / 角色；检查封面、页眉、PPT 备注、图片文字、文件属性 |
| numeric_density 低 | 补数值四要素、基线对照、读图句；删空话 |
| source_count / domains / official_ratio 低 | 补派研究员找 A/B 级来源；不要凑无关文献 |
| source_pairing 低 | 每张图表后补 `{资料来源：…}` 或图内 source |
| chart_count 低 | 按 outline 的必备图表补齐，用 charts.json 生成 |
| forbidden_terms | 逐条改写（“颠覆”→“替代”；“据统计”→“据国家统计局《…》”） |
| cliche_terms warn / fail | 按 `locations` 逐处改写为数字或删句（“赋能渠道”→“覆盖 1,200 家便利店”；“随着…的发展”→直接给增速与来源）；专业搭配加进 `checks.json.clicheWhitelist` |
| synthetic_research_markers fail | 在 `missingFiles` 列出的章节 / 讲稿 / 题库 / PPT 页补标记：正文段首 `【示例数据，须替换为真实调研结果】`，PPT `note` 末尾、讲稿段末、题库“证据”后 `（示例数据，须替换）`；不是把调研句删掉 |
| unsourced_numbers warn / fail | 对 `locations` 里的数字句补 `[n]` / 「资料来源：」/ 「测算 / 假设」字样，或把数字记入 `numbers.json`（record_numbers）后引用；读图句以“读图：”开头并紧跟图表 |
| has_primary_research fail | 补样本量 / 方法 / 时间 / 有效率句子（示例数据带标记） |
| has_competitor_matrix fail | 补 3–5 家 × ≥ 5 维表 |
| budget_recalc fail | 预算表每行 单价 × 数量 = 金额，合计与账本一致，机动金 5–10% |
| kpi_table fail | KPI 表含 指标 / 基线 / 目标 / 工具 / 周期 五列 ≥ 5 行 |
| three_statement_check / financing_consistency fail | 改 assumptions 重跑 finmodel |
| slide_metrics | 页数区间、每页 ≤ 120 字、结论句标题 ≥ 80%、含数字页 ≥ 80% |
| script_timing / sentence_length | 字数区间、句长 ≤ 25 |
| qa_bank | ≥ 40 题、六维度各 ≥ 4、三段式、主答人 |

## 文风三项（2026.10 起，正文 + 讲稿；`--chapter` 单章快检同样输出，不进硬门槛）
| 项 | 统计口径 | 阈值 | 配置键 |
|---|---|---|---|
| `cliche_terms` | 正文（不含附录）+ 讲稿中套话 / 模板句命中次数 ÷ 字数 × 10 000；`value` = 每万字次数，`locations` 给文件 / 行 / 原句 | ≤ 3 pass；3–9 warn；≥ 9 fail | `clicheTerms`（字面词；含 `…` 的模板句如「随着…发展」「助力…高质量发展」；`{"term": 正则, "regex": true}`）、`clicheWhitelist`、`clichePer10k`（默认 3） |
| `synthetic_research_markers` | `researchMode == synthetic`（或账本有 `synthetic:true`）时，引用一手调研结论的句子（调研显示 / 受访者认为 / 问卷·访谈·样本 + 计数或百分比 / 账本示例数值）所在的**文件**（正文章节、讲稿、题库）或 **PPT 页** 必须含「示例数据」；`value` = 缺标记的文件 / 页数，`missingFiles` 列出 | 0 pass；≥ 1 fail；非 synthetic 模式 skipped | 无（标记写法见 CLAUDE.md 示例数据铁律）；与 `synthetic_marker_coverage`（示例数值 ±3 行内须有标记）互补 |
| `unsourced_numbers` | 正文散文句里带 亿元 / 万元 / 万 / 亿 / 元 / % / 个百分点 / 倍 / 万人 的数字，同句或同段无 `[n]` / 示例数据 / 来源 / 测算·模型·假设 / 一手调研口径 / 「见第 N 章」指向，且不在 `numbers.json`（允许千分位、小数位、万↔亿↔元、20%↔0.2）的条数 ÷ 关键数字总数；表格单元格、年份、章节号、图表号、页码、文献号、读图句不计 | ≤ 10% pass；≤ 25% warn；> 25% fail | 无 |

## 打包前最终清单
PDF ≤ 页数上限、文件命名占位已说明、封面三字段、声明模板、数据真实性清单、README 待办（替换示例数据 / 填学校名 / 签名指印 / 录视频 / 授权书盖章）。
