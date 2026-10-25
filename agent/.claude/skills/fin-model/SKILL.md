---
name: fin-model
description: 用 assumptions.json + bizkit.finmodel 生成带公式的财务模型 Excel、numbers.json 数字账本与财务摘要。创业计划 / 案例类用 startup 模式（收入驱动、三表、盈亏平衡、单位经济、三情景、敏感性、DCF/PS/VC 估值、融资一致性）；策划类用 budget 模式（预算明细、机动金、KPI、ROI）。做财务阶段、调整假设或被评委质疑数字时使用。
---

# fin-model：假设 → 模型 → 账本

## 1. 口径与原则
- 单位：`assumptions.unit` 为展示换算（10000 = 万元）；输入用元；账本统一万元。
- 每条假设有 `source`（调研 CSV / 文献号 / 刊例 / 当地薪酬 / 行业对标），没有来源的假设不能进模型。
- 收入必须由驱动因子推出：`volume_price`（销量 × 单价）、`clients_arpu`（客户数 × ARPU）、`subscribers_arpu`（订阅数 × ARPU，含流失率）、`funnel`（流量 × 转化率 × 客单价）。
- 首年按月（S 型爬坡，12 个值），后续按年增长率；相邻年增速差异 > 5 个百分点；营销费率随规模递减；毛利率对标行业。
- 三大穿帮零容忍：等差 / 等比增长；永远盈利（首年应亏损或微利）；融资额 ÷ 投后 ≠ 出让比例。

## 2. startup 模式 assumptions.json（契约 §7.4 精简版）
```jsonc
{ "mode": "startup", "company": "弘宝汝瓷文创", "currency": "元", "unit": 10000, "startYear": 2026, "years": 5, "monthlyFirstYear": true,
  "scenarios": { "base": 1.0, "optimistic": 1.2, "pessimistic": 0.8 },
  "revenueStreams": [
    { "name": "文创礼品线上零售", "driver": "volume_price", "volume": { "y1Monthly": [200,220,260,300,340,380,420,460,520,600,700,800], "yoyGrowth": [0,0.8,0.5,0.35,0.25] },
      "price": { "base": 268, "yoyChange": [0,0.03,0.03,0.03,0.03] }, "cogsRatio": 0.42, "source": "单价=天猫同类 TOP10 均价 [12]；首年销量=渠道访谈 + SOM 的 8%" },
    { "name": "企业定制", "driver": "clients_arpu", "clients": { "y1": 12, "yoyGrowth": [0,1.0,0.6,0.4,0.3] }, "arpu": { "base": 50000 }, "cogsRatio": 0.55, "source": "…" } ],
  "costs": { "fixed": [ { "name": "场地租金", "monthly": 15000, "yoyChange": 0.05 } ],
             "variablePctOfRevenue": [ { "name": "平台佣金与支付", "pct": 0.06 } ],
             "marketing": { "pctOfRevenue": [0.18,0.15,0.12,0.10,0.09], "minAnnual": 300000 },
             "headcount": [ { "role": "品牌/运营", "y": [3,5,7,9,11], "monthlySalary": 9000, "socialRate": 0.30 } ],
             "capex": [ { "name": "设备与装修", "year": 1, "amount": 400000, "depreciationYears": 5 } ],
             "otherOpexPctOfRevenue": 0.03 },
  "workingCapital": { "receivableDays": 30, "inventoryDays": 60, "payableDays": 45 },
  "tax": { "incomeTaxRate": 0.25 },
  "financing": [ { "round": "种子轮", "year": 1, "month": 1, "amount": 3000000, "preMoney": 12000000, "use": [ { "item": "产品研发", "pct": 0.35 }, { "item": "品牌营销", "pct": 0.40 }, { "item": "团队与运营", "pct": 0.25 } ] } ],
  "openingCash": 500000,
  "valuation": { "methods": ["dcf", "ps_comparables", "vc_method"], "wacc": 0.18, "terminalGrowth": 0.03, "psMultiple": 3.0, "exitPE": 20, "exitYear": 5, "targetIrr": 0.40 },
  "unitEconomics": { "cac": 85, "repeatRate": 0.35, "grossMarginRef": 0.58, "lifetimeOrders": 3.2 },
  "sensitivity": { "drivers": ["price", "volume", "cogsRatio", "cac", "fixedCost"], "range": [-0.3, -0.15, 0.15, 0.3] } }
```
产出：`财务模型.xlsx`（说明 / 假设 / 驱动 / 收入 / 成本 / 人力 / 利润表 / 现金流量表 / 资产负债表 / 盈亏平衡 / 单位经济 / 敏感性 / 估值 / 融资 / 检查）、`numbers.json`、`财务摘要.md`、`财务图表.charts.json`。

## 3. budget 模式（品牌 / 数字营销 / 会奖）
```jsonc
{ "mode": "budget", "currency": "元", "totalCap": 3000000,
  "phases": [ { "name": "第一阶段：声量引爆（2026.4–6）", "items": [ { "item": "抖音信息流投放", "unitPrice": 35, "quantity": 20000, "unit": "CPM", "source": "巨量引擎 2025 刊例 [18]" } ] } ],
  "contingencyPct": 0.08,
  "kpis": [ { "metric": "曝光量", "baseline": "—", "target": "5000 万", "tool": "巨量云图", "cycle": "周" } ],
  "roi": { "incrementalRevenue": 4200000, "grossMargin": 0.55, "brandAssetNote": "品牌搜索指数目标 +120%" } }
```
产出：`财务模型.xlsx`（budget 模式内含 预算总表 / 排期 / KPI 页，文件名与 startup 模式一致以便打包）、`numbers.json`（budget_total、budget_contingency、kpi_*、roi）、`财务摘要.md`、`财务图表.charts.json`。另由建模师手写 `建模说明.md`（独立文件）。会奖赛道总额必须在 200–400 万之间。

## 4. 命令与验收
```bash
export PYTHONPATH=toolkit PYTHONUTF8=1
"${BIZGOLD_PYTHON:-python}" -m bizkit.finmodel --example budget    # 或 startup：打印完整可运行模板，字段即全部可用字段
"${BIZGOLD_PYTHON:-python}" -m bizkit.finmodel 03_财务/assumptions.json --out 03_财务
```
字段拿不准就看 `--example` 的模板照抄改数（未列出的键会被忽略），**不要去读 finmodel.py 源码**——它有 2700 行，读完只会浪费上下文。
读 `03_财务/财务摘要.md`：
- 检查页全部 0（资产 = 负债 + 权益；现金流期末 = 资产负债表货币资金；净利 → 留存收益）。
- 盈亏平衡 12–36 个月；现金最低点 + 融资额 ≥ 0；出让比例 = 融资额 ÷ 投后（±1 pp）。
- 估值 ≥ 2 种方法，区间合理；资金用途合计 = 100%。
- 账本 key 齐全：`revenue_y1..y5 net_profit_y1..y5 gross_margin_y1..y5 breakeven_month min_cash financing_amount pre_money post_money equity_offered cac ltv ltv_cac payback_months`（budget：`budget_total budget_contingency kpi_* roi`）。
不达标 → 改假设重跑（最多 3 轮），不要手改 xlsx。

## 5. 评委常问与模型应答
“销量为什么能翻倍？”→ 指向 SOM 与渠道容量假设来源；“毛利率怎么来的？”→ 行业对标文献号；“最坏情况？”→ 悲观情景现金最低点与应对；“估值依据？”→ 三法区间与取值理由；“钱花在哪？”→ 资金用途表与里程碑对应。
