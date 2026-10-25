---
name: biz-charts
description: 用 charts.json + bizkit.charts 批量生成 300 dpi 商务图表（18 种图型：柱 / 堆叠 / 条形 / 折线 / 面积 / 排序条形 / 瀑布 / 漏斗 / 龙卷风 / 盈亏平衡 / 现金流柱线 / 竞争矩阵 / 定位图 / 甘特 / 雷达 / 商业模式画布 / 环形）。做图表阶段、给正文或 PPT 配图、修复图表问题时使用。
---

# biz-charts：结论句决定图型

## 1. 选图表
| 你要说的结论 | 图型 `type` | data 要点 |
|---|---|---|
| 规模随时间增长、构成变化 | `stacked_bar` / `bar` / `area` | categories（年份）+ series |
| A 比 B 大、排名 | `hbar`（sort:true, highlight:"本品牌"） | categories + values |
| 占比（≤ 4 类） | `donut`；> 4 类用 `pie_as_bar` | categories + values |
| 趋势 / 多条线对比 | `line` | categories + series（≤ 4 条） |
| 从哪到哪（预算拆解、利润桥） | `waterfall` | items[{name,value,total?}] |
| 转化 / 漏斗 | `funnel` | stages[{name,value}] |
| 哪个参数最敏感 | `tornado` | base + items[{name,low,high}] |
| 卖多少才保本 | `breakeven` | fixedCost, unitPrice, unitVarCost |
| 现金流与累计 | `cashflow_combo` | categories + values (+cumulative) |
| 竞争格局 / 四象限 | `scatter_matrix` | xLabel,yLabel, points[{name,x,y,size?,self?}], quadrants[] |
| 品牌定位（语义两极） | `positioning_map` | xAxis{low,high}, yAxis{low,high}, points[] |
| 排期 / 里程碑 | `gantt` | tasks[{name,start,end,group}], milestones[] |
| 多维评分（品牌资产、壁垒） | `radar` | axes[≥3], series[{name,values}] |
| 商业模式 | `bmc` | 九格各数组 |
| 并列对比多组 | `grouped_bar` | categories + series |

硬规则：同屏定性色 ≤ 4；纵轴从 0 开始（折线除外）；> 8 类改横向条形；三个数以内写进句子不画图；图内不写标题（题注由 Word / PPT 承担）；每张图 `source` 必填（“资料来源：机构《报告》年份 [文献号]”或“资料来源：本方案财务模型 / 示例调研数据”）。

## 2. charts.json 模板（契约 §7.1）
```jsonc
{ "theme": "bizgold", "dpi": 300, "outputDir": "04_图表",
  "charts": [
    { "id": "fig_market_size", "type": "stacked_bar", "title": "汝瓷工艺品市场规模（2021—2026E）", "unit": "亿元", "yLabel": "市场规模",
      "source": "资料来源：中国陶瓷工业协会《2025 年陶瓷行业运行报告》[3]", "note": "2026E 为预测值",
      "data": { "categories": ["2021","2022","2023","2024","2025","2026E"], "series": [ { "name": "线下", "values": [17.1,18.2,19.5,20.1,21.0,21.8] }, { "name": "线上", "values": [4.0,5.1,7.2,9.6,12.5,15.9] } ] } },
    { "id": "fig_positioning", "type": "positioning_map", "title": "汝瓷品牌定位图",
      "source": "资料来源：团队竞品监测（2026 年 3 月）", "data": { "xAxis": { "low": "传统礼器", "high": "日常生活" }, "yAxis": { "low": "大众价位", "high": "高端收藏" },
      "points": [ { "name": "弘宝汝瓷", "x": 0.3, "y": 0.6, "self": true }, { "name": "竞品 A", "x": -0.5, "y": 0.8 }, { "name": "竞品 B", "x": 0.6, "y": -0.3 }, { "name": "竞品 C", "x": -0.2, "y": -0.6 } ] } },
    { "id": "fig_tornado", "type": "tornado", "title": "预算 ±20% 对核心 KPI 的影响", "unit": "万元", "source": "资料来源：本方案预算模型",
      "data": { "base": 300, "items": [ { "name": "KOL 投放单价", "low": 262, "high": 338 }, { "name": "信息流 CPM", "low": 280, "high": 320 } ] } } ] }
```
id 命名：`fig_<英文短名>`；文件名即 `<id>.png`。宽度：正文整幅 150–160 mm，半幅 80 mm；PPT 用图可用 160×90。

## 3. 命令
```bash
export PYTHONPATH=toolkit PYTHONUTF8=1
"${BIZGOLD_PYTHON:-python}" -m bizkit.charts 04_图表/charts.json            # 全部
"${BIZGOLD_PYTHON:-python}" -m bizkit.charts 04_图表/charts.json --only fig_a,fig_b
```
结果 JSON：`rendered[]`、`failed[]`（错误信息含缺失字段）、`warnings[]`；输出 `04_图表/index.json`。退出码 2 = 部分失败，修复后用 `--only` 重跑失败项。

## 4. 正文引用
`![汝瓷工艺品市场规模（2021—2026E）](../04_图表/fig_market_size.png)` → 渲染为“图 2-1 …”。图后两句读图：一句数值形态（“线上占比从 19% 升至 42%”），一句含义（“说明年轻消费者的入口在电商与内容平台”）。

## 5. 排错
- 中文方块：运行 `python -m bizkit.fonts` 查看可用字体，charts.json 顶层加 `"font": "Noto Sans CJK SC"`（或 Windows 的 `Microsoft YaHei`）。
- `donut` 超过 4 类自动改 `pie_as_bar`。
- 数值标签重叠：减少类别或改 `hbar`。
