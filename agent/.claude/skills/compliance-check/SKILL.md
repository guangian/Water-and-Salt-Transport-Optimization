---
name: compliance-check
description: 用 bizkit.checks 做机器检查与合规自检：结构覆盖、匿名、页数字数字体、数值密度、来源、图表、禁语与广告法、财务勾稽、融资一致性、PPT 指标、讲稿时长、题库、跨材料数字一致性；阅读检查结果并修复。每章写完、评审轮次、打包前使用。
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
| has_primary_research fail | 补样本量 / 方法 / 时间 / 有效率句子（示例数据带标记） |
| has_competitor_matrix fail | 补 3–5 家 × ≥ 5 维表 |
| budget_recalc fail | 预算表每行 单价 × 数量 = 金额，合计与账本一致，机动金 5–10% |
| kpi_table fail | KPI 表含 指标 / 基线 / 目标 / 工具 / 周期 五列 ≥ 5 行 |
| three_statement_check / financing_consistency fail | 改 assumptions 重跑 finmodel |
| slide_metrics | 页数区间、每页 ≤ 120 字、结论句标题 ≥ 80%、含数字页 ≥ 80% |
| script_timing / sentence_length | 字数区间、句长 ≤ 25 |
| qa_bank | ≥ 40 题、六维度各 ≥ 4、三段式、主答人 |

## 打包前最终清单
PDF ≤ 页数上限、文件命名占位已说明、封面三字段、声明模板、数据真实性清单、README 待办（替换示例数据 / 填学校名 / 签名指印 / 录视频 / 授权书盖章）。
