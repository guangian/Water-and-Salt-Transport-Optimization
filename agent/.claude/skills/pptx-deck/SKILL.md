---
name: pptx-deck
description: 用 slides.json + bizkit.pptx_render 生成 16:9 国金级路演 PPT（14 种页型、原生图表、演讲者备注、匿名）。设计路演页结构、写 slides.json、渲染或修改 PPT 时使用。
---

# pptx-deck：slides.json → 路演 PPT

## 1. 逐页结构（10 分钟 = 18–20 页正文 + 2–4 页备用；8 分钟 12–14 页；15 分钟 22–25 页）

**品牌策划 / 数字营销 / 会奖（策划案版）**
1 封面（标题范式 + 团队名）→ 2 钩子页（一个场景或一个数字，kpi 或 quote）→ 3 目录（agenda）→ 4 命题 / 品牌一句话与目标（bullets + highlight）→ 5–7 现状诊断（市场规模图 / 消费者洞察 kpi / 竞品矩阵 table）→ 8 问题树 → 9 定位逻辑链（two_column：细分 vs 目标市场）→ 10 定位图（image）→ 11 定位陈述与核心价值（quote）→ 12 Big Idea（quote + 视觉说明）→ 13 品牌元素系统（table：元素 / 设计 / 寓意）→ 14 传播策略总图（image 或 table：阶段 × 渠道）→ 15 KOL 与内容样稿（two_column）→ 16 线下事件（bullets/image）→ 17 执行排期（timeline）→ 18 预算与 KPI（chart + kpi）→ 19 资产保护与风险（table）→ 20 预期效果与回扣定位（kpi + 一句话）→ 21 致谢（closing）→ 备用页：敏感性 tornado、问卷样本结构、详细预算表、竞品案例。

**创业计划（BP 版）**
1 封面 → 2 痛点钩子 → 3 目录 → 4 一句话定位 + 最硬证据（kpi）→ 5 痛点与为什么是现在 → 6 产品 / 方案（image）→ 7 技术与创新（table 参数对比）→ 8 壁垒 → 9 市场规模 TAM/SAM/SOM（image funnel）→ 10 用户与一手调研（kpi）→ 11 竞争矩阵（table）→ 12 商业模式画布（image bmc）→ 13 单位经济（kpi）→ 14 营销与渠道（funnel/table）→ 15 牵引力 / 已完成事项（timeline）→ 16 团队（team）→ 17 财务预测（chart：收入与净利）→ 18 融资计划与资金用途（kpi + chart）→ 19 风险与对策（table）→ 20 社会价值与愿景（quote）→ 21 致谢 → 备用：敏感性、三表摘要、估值方法、知识产权清单。

## 2. 页型字段（契约 §7.3；路径相对 slides.json 所在目录）

| layout | 必填 | 可选 |
|---|---|---|
| cover | title | subtitle, note |
| agenda | title, items[] | |
| section | title | kicker |
| bullets | title, bullets[] (≤6, 每条 ≤40 字) | highlight, highlightLabel |
| image | title, image | caption, bullets[] (≤3) |
| two_column | title, left{heading,bullets}, right{heading,bullets} | |
| chart | title, chart{type: bar|column_clustered|line|pie|stacked, categories[], series[{name,values}], unit} | insight |
| table | title, columns[], rows[][] (≤7 列 ≤8 行) | note |
| kpi | title, kpis[{value,label}] (2–4) | |
| timeline | title, steps[{when,what}] (3–6) | |
| team | title, members[{role,duty}] | |
| quote | title, quote | by |
| closing | title | subtitle, contact |
| appendix | title | image / bullets / table 字段同上 |
每页通用：`note`（演讲者备注 = 该页讲稿段）、`duration`（秒）。页数区间、备用页上限、题库维度以 `tracks/<赛道>/checks.json.pitch` 为准（主代理会抄进任务书）。顶层：`title, subtitle, aspect:"16:9", theme:"bizgold", minutes, footer, anonymize{forbidden[]}, output`。

## 3. 设计规则（渲染器会校验并警告）
- 标题 = 结论句（含数字或判断动词），≤ 28 字；正文 ≤ 120 字；要点 ≤ 6 条；字号标题 ≥ 32 pt、正文 ≥ 24 pt、图注 ≥ 16 pt（国赛投影标准，评委隔几米看）；全稿 ≤ 3 色（主色 / 强调色 / 中性灰）；溢出只能删字不能缩字。
- 示例数据传播：凡引用 `numbers.json` 中 `synthetic: true` 的数字，该页 `note` 末尾加 `（示例数据，须替换）`。
- 每页至少 1 个数字（kpi / highlight / 图表 / 表格）；图表 ≥ 8 页；数字与账本一致。
- 一页一个信息层级：不要图 + 表 + 六条要点堆在一页。
- 颜色由主题控制（海军蓝 / 金），不要在文字里写颜色要求。
- 匿名：`anonymize.forbidden` 填院校词；备注里也不要出现姓名。
- 备用页统一放最后并用 `appendix` 页型（渲染器加“备用页”角标，讲时跳过）。

## 4. 命令与结果
```bash
export PYTHONPATH=toolkit PYTHONUTF8=1
"${BIZGOLD_PYTHON:-python}" -m bizkit.pptx_render 06_路演/slides.json
```
输出 `产出/路演.pptx` 与 `产出/路演.pptx.preview.json`（前端预览用）。结果 JSON：`pages`、`warnings[]`（文字溢出 / 字号过小 / 缺图）、`perPage[{chars,minFont}]`。警告必须清零：精简文字而不是缩字号。

## 5. 与讲稿的配合
每页 `note` 就是讲稿该段；`duration` 之和 ≈ 分钟 × 60；封面与致谢各 ≤ 20 秒；核心页（定位 / Big Idea / 预算或财务）各 ≥ 45 秒。
