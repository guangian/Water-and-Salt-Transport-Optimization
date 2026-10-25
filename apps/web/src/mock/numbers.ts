import type { NumbersData } from '../types';

export const MOCK_NUMBERS_BRAND: NumbersData = {
  generatedAt: '2026-10-18T10:40:00.000Z',
  unit: '万元',
  entries: [
    { key: 'budget_total', label: '推广预算总额', value: 292, unit: '万元', source: '03_财务/财务模型.xlsx!预算表!F10', usedIn: ['05_正文/07_财务预算.md', '06_路演/slides.json#10', '06_路演/讲稿.md'], consistency: 'ok' },
    { key: 'budget_phase1', label: '第一阶段预算（声量引爆）', value: 124, unit: '万元', source: '03_财务/财务模型.xlsx!预算表!F3', usedIn: ['05_正文/07_财务预算.md'], consistency: 'ok' },
    { key: 'budget_phase2', label: '第二阶段预算（种草转化）', value: 60, unit: '万元', source: '03_财务/财务模型.xlsx!预算表!F5', usedIn: ['05_正文/07_财务预算.md'], consistency: 'ok' },
    { key: 'budget_phase3', label: '第三阶段预算（体验沉淀）', value: 42, unit: '万元', source: '03_财务/财务模型.xlsx!预算表!F7', usedIn: ['05_正文/07_财务预算.md'], consistency: 'ok' },
    { key: 'kpi_exposure', label: '曝光量目标', value: 5000, unit: '万次', source: '03_财务/财务模型.xlsx!KPI!C2', usedIn: ['05_正文/07_财务预算.md', '06_路演/讲稿.md P9'], consistency: 'mismatch' },
    { key: 'kpi_engagement', label: '互动量目标', value: 120, unit: '万次', source: '03_财务/财务模型.xlsx!KPI!C3', usedIn: ['05_正文/07_财务预算.md'], consistency: 'ok' },
    { key: 'kpi_online_gmv', label: '线上 GMV 目标', value: 480, unit: '万元', source: '03_财务/财务模型.xlsx!KPI!C4', usedIn: ['05_正文/00_概要.md', '05_正文/07_财务预算.md', '06_路演/slides.json#7'], consistency: 'ok' },
    { key: 'roi_y1', label: '首年传播 ROI', value: 1.6, unit: '倍', source: '03_财务/财务模型.xlsx!效果!D8', usedIn: ['05_正文/07_财务预算.md', '06_路演/答辩题库.md'], consistency: 'ok' },
    { key: 'tam', label: '汝瓷工艺品市场规模（2025E）', value: 33.5, unit: '亿元', source: '01_调研/数据/市场规模.csv', usedIn: ['05_正文/02_品牌现状分析.md', '06_路演/slides.json#4'], consistency: 'ok' },
    { key: 'online_share_2025', label: '线上渠道占比（2025E）', value: 0.445, unit: '%', source: '01_调研/数据/市场规模.csv', usedIn: ['05_正文/02_品牌现状分析.md'], consistency: 'ok' },
    { key: 'som_y1', label: '首年可获得市场 SOM', value: 480, unit: '万元', source: '02_策略/策略框架.md §3', usedIn: ['05_正文/03_品牌定位与核心价值.md'], consistency: 'ok' },
    { key: 'sample_size', label: '问卷有效样本量', value: 412, unit: '份', source: '01_调研/调研工具包/示例结果.md', usedIn: ['05_正文/02_品牌现状分析.md', '06_路演/slides.json#9'], synthetic: true, consistency: 'ok' },
    { key: 'premium_willing', label: '愿意溢价 20% 的比例', value: 0.72, unit: '%', source: '01_调研/调研工具包/示例结果.md', usedIn: ['05_正文/02_品牌现状分析.md', '06_路演/slides.json#9'], synthetic: true, consistency: 'ok' },
    { key: 'conversion_rate', label: '到店转化率', value: 0.032, unit: '%', source: '03_财务/assumptions.json', usedIn: ['05_正文/07_财务预算.md'], synthetic: true, consistency: 'ok' },
    { key: 'competitor_price_main', label: '东道汝窑主力价格带', value: 1280, unit: '元', source: '01_调研/数据/竞品价格.csv', usedIn: ['05_正文/02_品牌现状分析.md', '06_路演/slides.json#8'], consistency: 'ok' },
    { key: 'search_growth', label: '“茶席美学”搜索量两年增幅', value: 1.8, unit: '%', source: '01_调研/03_消费者洞察.md', usedIn: ['05_正文/02_品牌现状分析.md'], synthetic: true, consistency: 'unchecked' },
  ],
  consistency: {
    consistent: false,
    mismatches: [{ key: 'kpi_exposure', label: '曝光量目标', book: '5000 万', deck: '5000 万', script: '4800 万' }],
  },
};

export const MOCK_NUMBERS_STARTUP: NumbersData = {
  generatedAt: '2026-10-24T16:40:00.000Z',
  unit: '万元',
  entries: [
    { key: 'revenue_y1', label: '第 1 年营业收入', value: 406, unit: '万元', source: '03_财务/财务模型.xlsx!利润表!C5', usedIn: ['05_正文/00_执行摘要.md'], consistency: 'ok' },
    { key: 'revenue_y3', label: '第 3 年营业收入', value: 1900, unit: '万元', source: '03_财务/财务模型.xlsx!利润表!E5', usedIn: ['05_正文/00_执行摘要.md'], consistency: 'ok' },
    { key: 'revenue_y5', label: '第 5 年营业收入', value: 3940, unit: '万元', source: '03_财务/财务模型.xlsx!利润表!G5', usedIn: [], consistency: 'unchecked' },
    { key: 'gross_margin_y3', label: '第 3 年毛利率', value: 0.62, unit: '%', source: '03_财务/财务模型.xlsx!利润表!E8', usedIn: [], consistency: 'unchecked' },
    { key: 'net_profit_y3', label: '第 3 年净利润', value: 352.5, unit: '万元', source: '03_财务/财务模型.xlsx!利润表!E14', usedIn: [], consistency: 'unchecked' },
    { key: 'breakeven_month', label: '盈亏平衡月份', value: 19, unit: '月', source: '03_财务/财务模型.xlsx!盈亏平衡!B12', usedIn: [], consistency: 'unchecked' },
    { key: 'min_cash', label: '现金最低点', value: -212.5, unit: '万元', source: '03_财务/财务模型.xlsx!现金流量表!C20', usedIn: [], consistency: 'unchecked' },
    { key: 'financing_amount', label: '融资额', value: 300, unit: '万元', source: '03_财务/assumptions.json', usedIn: ['05_正文/00_执行摘要.md'], consistency: 'ok' },
    { key: 'pre_money', label: '投前估值', value: 1200, unit: '万元', source: '03_财务/assumptions.json', usedIn: [], consistency: 'unchecked' },
    { key: 'post_money', label: '投后估值', value: 1500, unit: '万元', source: '03_财务/财务模型.xlsx!估值!C9', usedIn: [], consistency: 'unchecked' },
    { key: 'equity_offered', label: '出让比例', value: 0.2, unit: '%', source: '03_财务/财务模型.xlsx!估值!C10', usedIn: [], consistency: 'unchecked' },
    { key: 'cac', label: '客户获取成本', value: 4.2, unit: '万元', source: '03_财务/assumptions.json', usedIn: [], synthetic: true, consistency: 'unchecked' },
    { key: 'ltv', label: '客户生命周期价值', value: 26.8, unit: '万元', source: '03_财务/财务模型.xlsx!单位经济!C6', usedIn: [], synthetic: true, consistency: 'unchecked' },
    { key: 'ltv_cac', label: 'LTV / CAC', value: 6.4, unit: '倍', source: '03_财务/财务模型.xlsx!单位经济!C7', usedIn: [], synthetic: true, consistency: 'unchecked' },
    { key: 'tam', label: '惠州电子信息制造业产值（2025）', value: 5420, unit: '亿元', source: '01_调研/数据/电子信息产业规模.csv', usedIn: ['05_正文/02_市场分析.md'], consistency: 'ok' },
    { key: 'sam', label: '可服务市场 SAM', value: 38.6, unit: '亿元', source: '01_调研/01_行业规模与政策.md §2', usedIn: ['05_正文/02_市场分析.md'], consistency: 'ok' },
    { key: 'som_y1', label: '首年可获得市场 SOM', value: 420, unit: '万元', source: '02_策略/策略框架.md §3', usedIn: [], consistency: 'unchecked' },
    { key: 'sample_size', label: '企业访谈数', value: 12, unit: '家', source: '01_调研/03_客户访谈与需求.md', usedIn: ['05_正文/02_市场分析.md'], consistency: 'ok' },
  ],
  consistency: { consistent: true, mismatches: [] },
};
