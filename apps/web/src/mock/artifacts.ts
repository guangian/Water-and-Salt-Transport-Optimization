import type { ArtifactNode, Preview, Project, SheetData, SlidesDeck } from '../types';

// ---------- 工具：树构造 ----------

interface FileSpec {
  path: string;
  size: number;
  mtime?: string;
}

export function buildTree(files: FileSpec[], rootName = '工作区'): ArtifactNode {
  const root: ArtifactNode = { name: rootName, path: '', type: 'dir', children: [] };
  for (const f of files) {
    const parts = f.path.split('/');
    let node = root;
    let acc = '';
    parts.forEach((part, i) => {
      acc = acc ? `${acc}/${part}` : part;
      const isLast = i === parts.length - 1;
      node.children ??= [];
      let child = node.children.find((c) => c.name === part);
      if (!child) {
        child = isLast
          ? { name: part, path: acc, type: 'file', size: f.size, mtime: f.mtime ?? '2026-10-18T13:20:00.000Z' }
          : { name: part, path: acc, type: 'dir', children: [], mtime: f.mtime };
        node.children.push(child);
      }
      node = child;
    });
  }
  const sortRec = (n: ArtifactNode) => {
    if (!n.children) return;
    n.children.sort((a, b) => (a.type === b.type ? a.name.localeCompare(b.name, 'zh-CN') : a.type === 'dir' ? -1 : 1));
    n.children.forEach(sortRec);
  };
  sortRec(root);
  return root;
}

// ---------- 工具：SVG 图表（模拟 bizkit.charts 输出） ----------

const PALETTE = ['#14305C', '#C9A227', '#5B6B82', '#2E8B6E', '#B5533C', '#7F97BF'];

export function svgBarChart(opts: {
  title: string;
  categories: string[];
  series: { name: string; values: number[] }[];
  unit?: string;
  source?: string;
  stacked?: boolean;
}): string {
  const W = 960;
  const H = 540;
  const pad = { l: 80, r: 40, t: 70, b: 90 };
  const plotW = W - pad.l - pad.r;
  const plotH = H - pad.t - pad.b;
  const n = opts.categories.length;
  const maxVal = Math.max(
    ...opts.categories.map((_, i) =>
      opts.stacked ? opts.series.reduce((s, se) => s + (se.values[i] ?? 0), 0) : Math.max(...opts.series.map((se) => se.values[i] ?? 0)),
    ),
  );
  const niceMax = Math.ceil(maxVal / 5) * 5 || 1;
  const groupW = plotW / n;
  const barW = opts.stacked ? groupW * 0.5 : (groupW * 0.7) / opts.series.length;
  let bars = '';
  opts.categories.forEach((cat, i) => {
    let yAcc = pad.t + plotH;
    opts.series.forEach((se, j) => {
      const v = se.values[i] ?? 0;
      const h = (v / niceMax) * plotH;
      const x = opts.stacked ? pad.l + i * groupW + (groupW - barW) / 2 : pad.l + i * groupW + groupW * 0.15 + j * barW;
      const y = opts.stacked ? yAcc - h : pad.t + plotH - h;
      bars += `<rect x="${x.toFixed(1)}" y="${y.toFixed(1)}" width="${barW.toFixed(1)}" height="${h.toFixed(1)}" fill="${PALETTE[j % PALETTE.length]}" rx="2"/>`;
      if (!opts.stacked || j === opts.series.length - 1) {
        const labelV = opts.stacked ? opts.series.reduce((s, s2) => s + (s2.values[i] ?? 0), 0) : v;
        const labelText = Number.isInteger(labelV) ? String(labelV) : labelV.toFixed(1);
        bars += `<text x="${(x + barW / 2).toFixed(1)}" y="${(y - 6).toFixed(1)}" font-size="14" text-anchor="middle" fill="#334155">${labelText}</text>`;
      }
      if (opts.stacked) yAcc = y;
    });
    bars += `<text x="${(pad.l + i * groupW + groupW / 2).toFixed(1)}" y="${pad.t + plotH + 24}" font-size="15" text-anchor="middle" fill="#334155">${cat}</text>`;
  });
  let grid = '';
  for (let k = 0; k <= 5; k++) {
    const y = pad.t + plotH - (plotH * k) / 5;
    grid += `<line x1="${pad.l}" y1="${y}" x2="${W - pad.r}" y2="${y}" stroke="#E2E8F0" stroke-width="1"/>`;
    grid += `<text x="${pad.l - 10}" y="${y + 5}" font-size="13" text-anchor="end" fill="#64748B">${((niceMax * k) / 5).toFixed(0)}</text>`;
  }
  const legend = opts.series
    .map(
      (se, j) =>
        `<g transform="translate(${pad.l + j * 140}, ${H - 48})"><rect width="14" height="14" rx="2" fill="${PALETTE[j % PALETTE.length]}"/><text x="20" y="12" font-size="14" fill="#334155">${se.name}</text></g>`,
    )
    .join('');
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Microsoft YaHei, PingFang SC, Noto Sans CJK SC, sans-serif">
<rect width="${W}" height="${H}" fill="#FFFFFF"/>
<text x="${pad.l}" y="40" font-size="22" font-weight="700" fill="#14305C">${opts.title}</text>
${opts.unit ? `<text x="${pad.l}" y="${pad.t - 12}" font-size="13" fill="#64748B">单位：${opts.unit}</text>` : ''}
${grid}${bars}${legend}
${opts.source ? `<text x="${pad.l}" y="${H - 16}" font-size="12" fill="#94A3B8">${opts.source}</text>` : ''}
</svg>`;
  return `data:image/svg+xml;charset=utf-8,${encodeURIComponent(svg)}`;
}

export const MOCK_IMAGES: Record<string, string> = {
  '04_图表/fig_market_size.png': svgBarChart({
    title: '中国汝瓷工艺品市场规模（2019—2025E）',
    categories: ['2019', '2020', '2021', '2022', '2023', '2024', '2025E'],
    series: [
      { name: '线下', values: [12.1, 10.8, 13.2, 12.9, 15.4, 17.0, 18.6] },
      { name: '线上', values: [3.2, 4.9, 6.8, 8.1, 10.2, 12.5, 14.9] },
    ],
    unit: '亿元',
    source: '资料来源：中国陶瓷工业协会《2025 年陶瓷行业运行报告》；团队整理。2025E 为估计值。',
    stacked: true,
  }),
  '04_图表/fig_competitor_price.png': svgBarChart({
    title: '主要竞品茶器礼盒价格带对比',
    categories: ['弘宝汝瓷', '东道汝窑', '朱氏汝瓷', '玉松汝瓷', '景德镇 A'],
    series: [
      { name: '入门款', values: [268, 399, 298, 329, 459] },
      { name: '主力款', values: [680, 1280, 760, 890, 1580] },
    ],
    unit: '元',
    source: '资料来源：天猫、京东旗舰店 2026-10 公开售价；团队整理。',
  }),
  '04_图表/fig_budget_split.png': svgBarChart({
    title: '年度传播预算分配（合计 300 万元）',
    categories: ['内容制作', '抖音投放', '小红书种草', '线下快闪', '公关与 KOL', '监测与应急'],
    series: [{ name: '预算（万元）', values: [54, 78, 60, 42, 42, 24] }],
    unit: '万元',
    source: '资料来源：03_财务/财务模型.xlsx!预算表。',
  }),
  '04_图表/fig_revenue_forecast.png': svgBarChart({
    title: '智检云眸 2027—2031 年收入预测（基准情景）',
    categories: ['2027', '2028', '2029', '2030', '2031'],
    series: [
      { name: 'SaaS 订阅', values: [286, 742, 1380, 2210, 3120] },
      { name: '边缘一体机', values: [120, 310, 520, 690, 820] },
    ],
    unit: '万元',
    source: '资料来源：03_财务/财务模型.xlsx!利润表；示例数据，须替换为真实调研结果。',
    stacked: true,
  }),
};

// 一个极简但合法的 PDF（Helvetica 文本），用于 Mock 的 PDF 预览
const MINI_PDF = `%PDF-1.4
1 0 obj << /Type /Catalog /Pages 2 0 R >> endobj
2 0 obj << /Type /Pages /Kids [3 0 R] /Count 1 >> endobj
3 0 obj << /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >> endobj
4 0 obj << /Length 190 >> stream
BT /F1 26 Tf 72 740 Td (BizGold - Sample PDF Preview) Tj ET
BT /F1 14 Tf 72 700 Td (Mock mode: the real server converts docx to PDF.) Tj ET
BT /F1 14 Tf 72 680 Td (Pages: 63 / Limit: 80) Tj ET
endstream
endobj
5 0 obj << /Type /Font /Subtype /Type1 /BaseFont /Helvetica >> endobj
trailer << /Root 1 0 R >>
%%EOF`;

export const MOCK_PDF_URL = `data:application/pdf;base64,${typeof btoa === 'function' ? btoa(MINI_PDF) : ''}`;

// ---------- slides.json（12 页） ----------

export const MOCK_SLIDES: SlidesDeck = {
  title: '弘宝汝瓷品牌策划路演',
  subtitle: '青瓷工作室 · 本科组',
  aspect: '16:9',
  theme: 'bizgold',
  minutes: 10,
  footer: '弘宝汝瓷品牌策划 · 青瓷工作室',
  slides: [
    { layout: 'cover', title: '宋韵落天青，一釉生万象', subtitle: '弘宝汝瓷品牌策划 · 青瓷工作室', note: '开场 20 秒：从“天青色等烟雨”切入，说明汝瓷的文化稀缺性与弘宝的机会。', duration: 25 },
    { layout: 'agenda', title: '目录', items: ['现状诊断', '定位与核心价值', '品牌元素', '传播执行', '预算与效果', '风险与保护'], duration: 20 },
    { layout: 'section', title: '一、现状诊断', kicker: '问题在哪', duration: 10 },
    {
      layout: 'bullets',
      title: '汝瓷市场 6 年增长 2.3 倍，但弘宝线上份额不足 3%',
      bullets: ['线上渠道 2019—2025E 复合增速 29%，线下仅 7%', '弘宝 92% 收入来自汝州本地门店与批发', '品牌认知停留在“工艺品”，缺少生活方式叙事', '25–40 岁新中产对“茶席美学”搜索量两年增长 180%'],
      highlight: '18.6 亿元',
      highlightLabel: '2025E 线下市场规模',
      note: '强调“增长在线上、我们在线下”的错位。',
      duration: 60,
    },
    { layout: 'image', title: '市场规模：线上份额从 21% 升至 44%', image: '../04_图表/fig_market_size.png', caption: '资料来源：中国陶瓷工业协会《2025 年陶瓷行业运行报告》；团队整理', bullets: ['读图一句话：增量来自线上，弘宝尚未入场'], duration: 40 },
    {
      layout: 'two_column',
      title: '定位：从“汝州工艺品”到“宋韵生活方式品牌”',
      left: { heading: '现状', bullets: ['区域性非遗工艺品', '以器物工艺为卖点', '礼品渠道为主'] },
      right: { heading: '目标', bullets: ['面向新中产的宋韵生活方式', '以“茶席美学”为场景', '线上内容驱动 + 线下体验'] },
      duration: 50,
    },
    {
      layout: 'chart',
      title: '三年收入目标：线上占比提升至 35%',
      chart: { type: 'column_clustered', categories: ['2026', '2027', '2028'], series: [{ name: '线下', values: [1420, 1560, 1690] }, { name: '线上', values: [120, 480, 910] }], unit: '万元' },
      insight: '线上收入三年增长 6.6 倍，带动总收入复合增速 19%',
      duration: 45,
    },
    {
      layout: 'table',
      title: '竞品对比：弘宝的价格带优势与内容短板',
      columns: ['维度', '东道汝窑', '朱氏汝瓷', '弘宝汝瓷'],
      rows: [
        ['主力价格带', '1280 元', '760 元', '680 元'],
        ['线上粉丝量', '42 万', '9.8 万', '1.2 万'],
        ['内容更新频率', '日更', '周更', '月更'],
        ['非遗传承人背书', '有', '有', '有（省级）'],
      ],
      note: '数据为 2026-10 公开信息，团队整理。',
      duration: 45,
    },
    {
      layout: 'kpi',
      title: '调研关键发现',
      kpis: [
        { value: '412', label: '有效问卷（示例数据）' },
        { value: '72%', label: '愿意为“宋韵美学”溢价 20%' },
        { value: '3.4 次', label: '目标人群年均茶器购买频次' },
        { value: '58%', label: '通过短视频首次了解汝瓷' },
      ],
      duration: 40,
    },
    {
      layout: 'timeline',
      title: '传播执行：三阶段 12 个月',
      steps: [
        { when: '2027 Q1', what: '品牌焕新发布 + 天青色视觉系统上线' },
        { when: '2027 Q2', what: '抖音 / 小红书“一釉生万象”内容季' },
        { when: '2027 Q3', what: '郑州 / 杭州茶博会快闪体验店' },
        { when: '2027 Q4', what: '年度礼盒上市 + 会员体系沉淀' },
      ],
      duration: 45,
    },
    { layout: 'quote', title: '品牌主张', quote: '一釉生万象', by: 'Big Idea · 天青色不是颜色，是一种生活的留白', duration: 20 },
    { layout: 'closing', title: '谢谢', subtitle: '恳请评委指正', contact: '', duration: 10 },
  ],
  output: '产出/路演.pptx',
};

// ---------- Markdown 章节 ----------

export const MOCK_CHAPTER_MD = `# 品牌现状分析

## 外部环境：政策与文化红利叠加，增量市场在线上

**结论先行**：汝瓷工艺品市场 2019—2025E 从 15.3 亿元增至 33.5 亿元，增长 2.2 倍；其中线上渠道占比由 21% 升至 44%，是全部增量的 71% 来源。{来源：中国陶瓷工业协会《2025 年陶瓷行业运行报告》，2026-01}

![中国汝瓷工艺品市场规模（2019—2025E）](../04_图表/fig_market_size.png)

| 年份 | 线下（亿元） | 线上（亿元） | 线上占比 | 同比 |
|---|---|---|---|---|
| 2021 | 13.2 | 6.8 | 34.0% | +27.4% |
| 2022 | 12.9 | 8.1 | 38.6% | +5.0% |
| 2023 | 15.4 | 10.2 | 39.8% | +21.9% |
| 2024 | 17.0 | 12.5 | 42.4% | +15.2% |
| 2025E | 18.6 | 14.9 | 44.5% | +13.6% |

表：汝瓷工艺品市场规模与渠道结构

### 政策面
- 《“十五五”非物质文化遗产保护规划》提出“推动传统工艺高质量传承发展”，汝瓷烧制技艺列入国家级非遗代表性项目名录 [3]。
- 河南省 2025 年出台《关于促进文旅文创融合战略实施意见》，对非遗品牌化给予最高 200 万元奖补 [4]。

### 消费面
- 25–40 岁新中产对“茶席美学”“宋式生活”相关内容的搜索量两年增长 **180%**（示例数据，须替换为真实调研结果）[5]。
- 示例问卷（n=412）显示，72% 的受访者愿意为“有文化叙事的汝瓷器物”支付 20% 以上溢价。

## 内部环境：工艺与传承是资产，品牌叙事与渠道是短板

> 弘宝汝瓷拥有省级非遗传承人 2 名、天青釉配方专利 1 项，但 2025 年线上收入仅占 2.8%，品牌内容月更不足 2 条。

| 维度 | 优势（S） | 劣势（W） |
|---|---|---|
| 产品 | 天青釉发色稳定，开片均匀 | SKU 偏传统，缺少日用茶席场景产品 |
| 渠道 | 汝州门店 + 礼品批发稳定 | 线上近乎空白，无私域运营 |
| 品牌 | 省级传承人背书 | 品牌主张缺位，视觉系统陈旧 |
| 组织 | 家族作坊决策快 | 缺少品牌与内容人才 |

表：弘宝汝瓷内部环境 SW 分析

### 本章小结
弘宝汝瓷处在“赛道增长、自身掉队”的典型错位之中：增量在线上、在年轻新中产、在生活方式叙事，而弘宝的资产集中在工艺与本地渠道。**品牌策划的核心命题是：用一个可持续的生活方式叙事，把工艺资产转译为线上可传播、可购买的品牌资产。**
`;

export const MOCK_TASK_MD = `# 任务书
- 赛道：品牌策划赛道（brand） | 组别：本科组 | 语言：中文
- 作品名：宋韵天青——弘宝汝瓷品牌策划书 | 团队名：青瓷工作室（匿名制：全文不得出现院校/姓名/联系方式）
- 策划对象：品牌 弘宝汝瓷（汝州市弘宝汝瓷坊有限公司，河南汝州，非遗工艺品/文创）
- 是否命题企业：否 | 命题题目：—
- 项目简介：弘宝汝瓷是汝州本地非遗传承品牌，以天青釉茶器与文创礼品为主，线下渠道为主、线上份额不足 3%……
- 团队（匿名）：队长/总策划（市场营销，品牌传播）；调研与数据（统计学，问卷与数据分析）；创意与视觉（视觉传达设计，品牌元素设计）；财务预算（会计学，预算编制）
- 参数：调研模式=示例数据 | 正文目标 60 页（上限 80）| 讲稿 10 分钟 | 评审轮数 3 | 目标分 85
- 用户资料（已预读为文本，路径见下）：
  - inputs/弘宝汝瓷企业资料.pdf → inputs/_preread/弘宝汝瓷企业资料.pdf.txt（2.3 万字）
  - inputs/2025 销售数据.xlsx → inputs/_preread/2025 销售数据.xlsx.txt（3 个 sheet）
- 赛道规范：tracks/brand/spec.md | 量规：tracks/brand/rubric.json | 命名：tracks/brand/naming.json
- 开始指令：按 CLAUDE.md 执行流程完成全部阶段，直至 08_提交包 生成。
`;

export const MOCK_SCRIPT_MD = `# 路演讲稿（10 分钟 · 约 2600 字）

> 节奏：开场 0:25 → 现状 2:00 → 定位 2:00 → 元素与传播 3:00 → 预算与效果 1:30 → 风险与收尾 1:05

## P1 开场（0:00–0:25）
各位评委好。“天青色等烟雨”，等的是一窑好瓷。今天我们要讲的，是如何让一个汝州的非遗作坊，被更多年轻人等到。

## P4 现状诊断（0:45–1:45）
过去六年，汝瓷工艺品市场增长了 2.2 倍，其中 **71% 的增量来自线上**。而弘宝汝瓷 92% 的收入仍来自汝州本地门店与礼品批发，线上份额不足 3%。问题不在产品，在叙事与渠道。

## P6 定位（2:25–3:15）
我们把弘宝从“汝州工艺品”重新定位为“宋韵生活方式品牌”——面向 25–40 岁新中产，以茶席美学为场景……

## P12 收尾（9:40–10:00）
一釉生万象。谢谢各位评委，恳请指正。
`;

export const MOCK_QA_MD = `# 答辩题库（评委视角 · 30 问）

## 一、商业与财务（10 问）
1. **300 万元预算的 ROI 如何测算？** —— 预算分 6 类，按渠道 CPM/CPE 刊例测算曝光 5000 万、互动 120 万、到店 + 成交转化 3.2%，对应新增线上收入 480 万元，首年 ROI 1.6（见 03_财务/财务模型.xlsx!预算表）。
2. **为什么线上占比目标是 35% 而不是更高？** —— 参考东道汝窑用 5 年把线上做到 41%，35% 是三年内克制且可达的目标……

## 二、品牌与创意（10 问）
11. **“一釉生万象”是否过于抽象？** —— 它是品牌主张，不是广告语；落地口号为“把宋朝的天色，搬到你的茶席”……

## 三、调研与数据（10 问）
21. **问卷样本量 412 份是否足够？** —— 置信水平 95%、误差 ±5% 对应最小样本 385 份；**本方案为示例数据，须替换为真实调研结果**……
`;

export const MOCK_FIXES_MD = `# 修改清单（第 3 轮）

| 优先级 | 位置 | 修改内容 | 预计提升 | 状态 |
|---|---|---|---|---|
| 1 | 02_品牌现状分析.md §2.3 | 补 412 份问卷的样本结构表（年龄/城市/收入） | +3 | 已完成 |
| 2 | 06_品牌推广与传播.md §6.4 | 补抖音/小红书投放的 CPM 刊例来源 | +2 | 已完成 |
| 3 | 讲稿.md P9 | 统一“曝光 5000 万”口径（讲稿写了 4800 万） | +1 | 待处理 |
| 4 | 08_品牌资产保护.md | 增加商标近似查询结果截图说明 | +1 | 已完成 |
`;

export const MOCK_README_MD = `# 提交包说明

- 文件命名规则：全国预选赛抽签编号-团队名称-作品名称（例：03-20-青瓷工作室-宋韵天青——弘宝汝瓷品牌策划书.pdf）
- 已包含：策划书 PDF（63 页）、作品原创性声明、作品版权声明、企业授权书模板、数据真实性清单
- **数据真实性清单**：共 6 条示例数据，提交前须替换为真实调研结果（见《数据真实性清单.md》）
- 匿名检查：0 处命中；禁语检查：0 处命中
`;

// ---------- xlsx sheets ----------

export const MOCK_SHEETS: SheetData[] = [
  {
    name: '假设',
    rows: [
      ['项目', '取值', '单位', '来源'],
      ['年度传播预算上限', 3000000, '元', '任务书'],
      ['抖音信息流 CPM', 35, '元/CPM', '巨量引擎 2025 刊例'],
      ['小红书薯条 CPE', 1.8, '元/互动', '小红书聚光 2025 刊例'],
      ['到店转化率', '3.2%', '%', '示例数据，须替换'],
      ['客单价（线上）', 268, '元', '天猫同类 TOP10 均价'],
      ['应急预备金比例', '8%', '%', '方法论默认'],
    ],
  },
  {
    name: '预算表',
    rows: [
      ['阶段', '项目', '单价', '数量', '单位', '小计（元）', '来源'],
      ['第一阶段：声量引爆（2027.1–3）', '品牌视觉系统与 TVC 制作', 540000, 1, '项', 540000, '供应商报价均值'],
      ['第一阶段：声量引爆（2027.1–3）', '抖音信息流投放', 35, 20000, 'CPM', 700000, '巨量引擎 2025 刊例'],
      ['第二阶段：种草转化（2027.4–9）', '小红书种草（KOC 200 篇）', 1500, 200, '篇', 300000, '蒲公英平台均价'],
      ['第二阶段：种草转化（2027.4–9）', '小红书薯条投放', 1.8, 166000, '互动', 300000, '聚光平台刊例'],
      ['第三阶段：体验沉淀（2027.10–12）', '郑州 / 杭州茶博会快闪店', 210000, 2, '场', 420000, '展位 + 搭建报价'],
      ['全年', '公关与 KOL 合作', 420000, 1, '项', 420000, '示例数据，须替换'],
      ['全年', '监测与应急预备金', 240000, 1, '项', 240000, '预算 8%'],
      ['合计', '', '', '', '', 2920000, ''],
    ],
  },
  {
    name: 'KPI',
    rows: [
      ['指标', '基线', '目标', '监测工具', '周期'],
      ['曝光量', '—', '5000 万', '巨量云图', '周'],
      ['互动量', '—', '120 万', '聚光 / 云图', '周'],
      ['线上 GMV', '42 万元', '480 万元', '生意参谋', '月'],
      ['品牌搜索指数', '120', '1500', '巨量算数', '月'],
      ['会员数', '0', '3 万', '微信 CRM', '月'],
    ],
  },
];

export const MOCK_STARTUP_SHEETS: SheetData[] = [
  {
    name: '利润表',
    rows: [
      ['科目（万元）', '2027', '2028', '2029', '2030', '2031'],
      ['营业收入', 406, 1052, 1900, 2900, 3940],
      ['营业成本', 170.5, 420.8, 722, 1044, 1379],
      ['毛利', 235.5, 631.2, 1178, 1856, 2561],
      ['毛利率', 0.58, 0.6, 0.62, 0.64, 0.65],
      ['销售费用', 73.1, 157.8, 228, 290, 354.6],
      ['研发费用', 180, 260, 340, 420, 480],
      ['管理费用', 96, 128, 160, 198, 236],
      ['净利润', -121.6, 60.3, 352.5, 735.5, 1163.4],
    ],
  },
  {
    name: '现金流量表',
    rows: [
      ['科目（万元）', '2027', '2028', '2029', '2030', '2031'],
      ['经营活动现金流', -98.2, 88.6, 401.2, 790.4, 1210.8],
      ['投资活动现金流', -140, -60, -80, -90, -100],
      ['筹资活动现金流', 300, 0, 0, 0, 0],
      ['期末现金', 111.8, 140.4, 461.6, 1162, 2272.8],
    ],
  },
  {
    name: '敏感性',
    rows: [
      ['驱动因素', '-30%', '-15%', '基准', '+15%', '+30%'],
      ['单价', -412, -83, 352.5, 788, 1223],
      ['销量', -395, -74, 352.5, 779, 1206],
      ['成本率', 622, 487, 352.5, 218, 83],
      ['CAC', 410, 381, 352.5, 324, 295],
    ],
  },
];

// ---------- 各项目文件清单 ----------

const BRAND_FILES: FileSpec[] = [
  { path: '任务书.md', size: 2_140, mtime: '2026-10-18T09:16:02.000Z' },
  { path: 'inputs/弘宝汝瓷企业资料.pdf', size: 2_341_002, mtime: '2026-10-18T09:15:40.000Z' },
  { path: 'inputs/2025 销售数据.xlsx', size: 184_320, mtime: '2026-10-18T09:15:41.000Z' },
  { path: 'inputs/产品图册.jpg', size: 1_204_112, mtime: '2026-10-18T09:15:41.000Z' },
  { path: 'inputs/_preread/弘宝汝瓷企业资料.pdf.txt', size: 70_230 },
  { path: 'inputs/_preread/2025 销售数据.xlsx.txt', size: 9_840 },
  { path: '00_任务解析/任务解析.md', size: 8_420, mtime: '2026-10-18T09:22:10.000Z' },
  { path: '01_调研/01_市场规模与趋势.md', size: 14_210 },
  { path: '01_调研/02_竞品与渠道.md', size: 12_870 },
  { path: '01_调研/03_消费者洞察.md', size: 11_300 },
  { path: '01_调研/数据/市场规模.csv', size: 1_240 },
  { path: '01_调研/数据/竞品价格.csv', size: 980 },
  { path: '01_调研/参考文献.json', size: 6_410 },
  { path: '01_调研/调研工具包/问卷.md', size: 7_800 },
  { path: '01_调研/调研工具包/访谈提纲.md', size: 3_200 },
  { path: '01_调研/调研工具包/执行计划.md', size: 2_900 },
  { path: '01_调研/调研工具包/示例结果.md', size: 5_600 },
  { path: '02_策略/策略框架.md', size: 9_950 },
  { path: '03_财务/assumptions.json', size: 3_110 },
  { path: '03_财务/财务模型.xlsx', size: 61_440 },
  { path: '03_财务/numbers.json', size: 2_420 },
  { path: '03_财务/财务摘要.md', size: 4_100 },
  { path: '04_图表/charts.json', size: 9_200 },
  { path: '04_图表/fig_market_size.png', size: 182_400 },
  { path: '04_图表/fig_competitor_price.png', size: 164_200 },
  { path: '04_图表/fig_budget_split.png', size: 150_900 },
  { path: '04_图表/index.json', size: 1_100 },
  { path: '05_正文/00_概要.md', size: 2_600 },
  { path: '05_正文/01_引言.md', size: 4_300 },
  { path: '05_正文/02_品牌现状分析.md', size: 16_800 },
  { path: '05_正文/03_品牌定位与核心价值.md', size: 15_200 },
  { path: '05_正文/04_品牌元素设计.md', size: 12_400 },
  { path: '05_正文/05_品牌推广与传播.md', size: 18_900 },
  { path: '05_正文/06_品牌资产保护与风险管理.md', size: 9_700 },
  { path: '05_正文/07_财务预算.md', size: 8_100 },
  { path: '05_正文/A1_调研问卷.md', size: 7_800 },
  { path: '05_正文/book.json', size: 2_300 },
  { path: '06_路演/slides.json', size: 11_900 },
  { path: '06_路演/讲稿.md', size: 9_600 },
  { path: '06_路演/答辩题库.md', size: 14_300 },
  { path: '06_路演/视频分镜.md', size: 5_200 },
  { path: '07_评审/第1轮/评审_investor.json', size: 6_100 },
  { path: '07_评审/第1轮/评审_professor.json', size: 6_400 },
  { path: '07_评审/第1轮/评审_executive.json', size: 6_000 },
  { path: '07_评审/第1轮/评审汇总.json', size: 3_300 },
  { path: '07_评审/第1轮/机器检查.json', size: 5_900 },
  { path: '07_评审/第1轮/修改清单.md', size: 2_800 },
  { path: '07_评审/第2轮/评审汇总.json', size: 3_300 },
  { path: '07_评审/第2轮/机器检查.json', size: 5_900 },
  { path: '07_评审/第2轮/修改清单.md', size: 2_400 },
  { path: '07_评审/第3轮/评审_investor.json', size: 6_200 },
  { path: '07_评审/第3轮/评审_professor.json', size: 6_500 },
  { path: '07_评审/第3轮/评审_executive.json', size: 6_100 },
  { path: '07_评审/第3轮/评审汇总.json', size: 3_400 },
  { path: '07_评审/第3轮/机器检查.json', size: 6_000 },
  { path: '07_评审/第3轮/修改清单.md', size: 1_900 },
  { path: '08_提交包/03-20-青瓷工作室-宋韵天青——弘宝汝瓷品牌策划书.pdf', size: 8_120_400, mtime: '2026-10-18T13:29:50.000Z' },
  { path: '08_提交包/作品原创性声明.docx', size: 24_100 },
  { path: '08_提交包/作品版权声明.docx', size: 23_800 },
  { path: '08_提交包/企业授权书模板.docx', size: 25_300 },
  { path: '08_提交包/数据真实性清单.md', size: 3_100 },
  { path: '08_提交包/README.md', size: 1_400 },
  { path: '08_提交包/提交包.zip', size: 9_410_200, mtime: '2026-10-18T13:30:12.000Z' },
  { path: '产出/弘宝汝瓷品牌策划书.docx', size: 2_940_100, mtime: '2026-10-18T12:40:10.000Z' },
  { path: '产出/弘宝汝瓷品牌策划书.pdf', size: 8_120_400, mtime: '2026-10-18T12:41:30.000Z' },
  { path: '产出/路演.pptx', size: 4_210_300, mtime: '2026-10-18T13:02:10.000Z' },
];

const STARTUP_FILES: FileSpec[] = [
  { path: '任务书.md', size: 2_300, mtime: '2026-10-24T15:31:02.000Z' },
  { path: 'inputs/惠州市十五五规划建议（节选）.pdf', size: 1_120_440 },
  { path: 'inputs/仲恺电子信息企业访谈纪要.docx', size: 88_210 },
  { path: '00_任务解析/任务解析.md', size: 9_100 },
  { path: '01_调研/01_行业规模与政策.md', size: 15_300 },
  { path: '01_调研/02_竞品与定价.md', size: 13_100 },
  { path: '01_调研/03_客户访谈与需求.md', size: 10_900 },
  { path: '01_调研/数据/电子信息产业规模.csv', size: 1_600 },
  { path: '01_调研/参考文献.json', size: 7_200 },
  { path: '01_调研/调研工具包/访谈提纲.md', size: 3_400 },
  { path: '02_策略/策略框架.md', size: 10_800 },
  { path: '03_财务/assumptions.json', size: 4_900 },
  { path: '03_财务/财务模型.xlsx', size: 88_400 },
  { path: '03_财务/numbers.json', size: 3_100 },
  { path: '03_财务/财务摘要.md', size: 4_600 },
  { path: '04_图表/charts.json', size: 8_800 },
  { path: '04_图表/fig_revenue_forecast.png', size: 176_300 },
  { path: '04_图表/fig_market_size.png', size: 180_000 },
  { path: '05_正文/00_执行摘要.md', size: 2_900 },
  { path: '05_正文/01_项目与产品.md', size: 14_600 },
  { path: '05_正文/02_市场分析.md', size: 16_100 },
];

export const PROJECT_FILES: Record<string, FileSpec[]> = {
  p_20261018_091533_a1f3: BRAND_FILES,
  [`p_20261024_153012_ab12`]: STARTUP_FILES,
  p_20261021_201145_c7d9: [
    { path: '任务书.md', size: 2_000 },
    { path: 'inputs/门店经营数据 2025Q1-Q3.csv', size: 41_020 },
    { path: '00_任务解析/任务解析.md', size: 7_600 },
    { path: '01_调研/01_平台与人群.md', size: 12_100 },
    { path: '01_调研/02_竞品账号拆解.md', size: 10_300 },
    { path: '02_策略/策略框架.md', size: 8_900 },
    { path: '03_财务/assumptions.json', size: 2_800 },
    { path: '03_财务/财务模型.xlsx', size: 48_000 },
    { path: '04_图表/charts.json', size: 6_200 },
  ],
  p_20261022_143301_e4b8: [
    { path: '任务书.md', size: 1_900 },
    { path: '00_任务解析/任务解析.md', size: 6_900 },
  ],
  p_20261024_101209_f0e2: [
    { path: '任务书.md', size: 2_100 },
    { path: 'inputs/科大讯飞 2024 年年度报告.pdf', size: 9_820_112 },
  ],
};

export function artifactsFor(project: Project): ArtifactNode {
  const files = PROJECT_FILES[project.id] ?? [{ path: '任务书.md', size: 2_000 }];
  return buildTree(files, project.name);
}

// ---------- 预览 ----------

const GENERIC_MD = (title: string) => `# ${title}

> 本文件由金策 BizGold Mock 模式生成，用于演示预览器。真实运行时此处为智能体产出的实际内容。

## 结论摘要
- 结论一：数值四要素（数值、口径、年份、来源）齐备。
- 结论二：每章末尾有小结句，呼应评分项。

## 证据表
| 指标 | 数值 | 口径 | 年份 | 来源 |
|---|---|---|---|---|
| 示例指标 | 18.6 亿元 | 线下零售额 | 2025E | 中国陶瓷工业协会 |
`;

export function previewFor(project: Project, path: string): Preview {
  const lower = path.toLowerCase();
  if (MOCK_IMAGES[path]) return { kind: 'image', url: MOCK_IMAGES[path] };
  if (lower.endsWith('.png') || lower.endsWith('.jpg') || lower.endsWith('.jpeg')) {
    return { kind: 'image', url: MOCK_IMAGES['04_图表/fig_market_size.png'] };
  }
  if (path === '06_路演/slides.json' || lower.endsWith('.pptx') || lower.endsWith('.preview.json')) {
    return { kind: 'slides', slides: project.track === 'brand' ? MOCK_SLIDES : { ...MOCK_SLIDES, title: `${project.subject.name}路演`, footer: `${project.subject.name} · ${project.teamName}` } };
  }
  if (lower.endsWith('.xlsx')) return { kind: 'table', sheets: project.track === 'startup' ? MOCK_STARTUP_SHEETS : MOCK_SHEETS };
  if (lower.endsWith('.csv')) {
    return {
      kind: 'table',
      sheets: [
        {
          name: 'csv',
          rows: [
            ['year', 'offline', 'online', 'source'],
            [2021, 13.2, 6.8, '中国陶瓷工业协会'],
            [2022, 12.9, 8.1, '中国陶瓷工业协会'],
            [2023, 15.4, 10.2, '中国陶瓷工业协会'],
            [2024, 17.0, 12.5, '中国陶瓷工业协会'],
          ],
        },
      ],
    };
  }
  if (lower.endsWith('.pdf')) return { kind: 'pdf', url: MOCK_PDF_URL };
  if (lower.endsWith('.docx') || lower.endsWith('.zip')) {
    return { kind: 'binary', size: 2_940_100, message: lower.endsWith('.docx') ? 'Word 文档不支持在线预览，请下载查看；同名 PDF 可直接预览。' : '压缩包请下载后查看。' };
  }
  if (lower.endsWith('.json')) {
    if (path.endsWith('numbers.json')) {
      return { kind: 'text', content: JSON.stringify({ generatedAt: '2026-10-18T10:40:00.000Z', unit: '万元', entries: [{ key: 'budget_total', label: '推广预算总额', value: 292, unit: '万元', source: '03_财务/财务模型.xlsx!预算表!F10' }] }, null, 2) };
    }
    if (path.endsWith('参考文献.json')) {
      return {
        kind: 'text',
        content: JSON.stringify(
          [
            { id: 1, type: 'report', author: '中国陶瓷工业协会', title: '2025 年陶瓷行业运行报告', publisher: '中国陶瓷工业协会', year: 2026, url: 'https://www.ccia.org.cn/report/2025', accessed: '2026-10-18' },
            { id: 2, type: 'gov', author: '文化和旅游部', title: '“十五五”非物质文化遗产保护规划', year: 2026, url: 'https://www.mct.gov.cn/', accessed: '2026-10-18' },
          ],
          null,
          2,
        ),
      };
    }
    return { kind: 'text', content: JSON.stringify({ note: 'Mock JSON 文件', path }, null, 2) };
  }
  if (lower.endsWith('.txt')) return { kind: 'text', content: `（预读文本节选）\n\n汝州市弘宝汝瓷坊有限公司成立于 2008 年，现有省级非遗传承人 2 名……\n\n2025 年营业收入 1540 万元，其中线上 42 万元（2.8%）。` };
  if (lower.endsWith('.md')) {
    if (path === '任务书.md') return { kind: 'markdown', content: MOCK_TASK_MD };
    if (path.includes('品牌现状分析') || path.includes('02_市场分析')) return { kind: 'markdown', content: MOCK_CHAPTER_MD };
    if (path.endsWith('讲稿.md')) return { kind: 'markdown', content: MOCK_SCRIPT_MD };
    if (path.endsWith('答辩题库.md')) return { kind: 'markdown', content: MOCK_QA_MD };
    if (path.endsWith('修改清单.md')) return { kind: 'markdown', content: MOCK_FIXES_MD };
    if (path.endsWith('README.md')) return { kind: 'markdown', content: MOCK_README_MD };
    if (path.endsWith('数据真实性清单.md')) {
      return {
        kind: 'markdown',
        content: `# 数据真实性清单\n\n| # | 数据 | 出现位置 | 当前值 | 性质 | 替换要求 |\n|---|---|---|---|---|---|\n| 1 | 问卷有效样本量 | 02_品牌现状分析 §2.3、路演 P9 | 412 份 | 示例数据 | 替换为真实问卷结果 |\n| 2 | 愿意溢价比例 | 02_品牌现状分析 §2.3 | 72% | 示例数据 | 替换为真实问卷结果 |\n| 3 | 到店转化率 | 07_财务预算 §7.2 | 3.2% | 示例数据 | 以试点投放数据替换 |\n`,
      };
    }
    return { kind: 'markdown', content: GENERIC_MD(path.split('/').pop()?.replace(/\.md$/, '') ?? path) };
  }
  return { kind: 'binary', message: '不支持预览的文件类型。' };
}

export function mockFileUrl(project: Project, path: string): string {
  if (MOCK_IMAGES[path]) return MOCK_IMAGES[path];
  const p = previewFor(project, path);
  if (p.kind === 'image' && p.url) return p.url;
  if (p.kind === 'pdf' && p.url) return p.url;
  if (p.kind === 'markdown' || p.kind === 'text') return `data:text/plain;charset=utf-8,${encodeURIComponent(p.content ?? '')}`;
  if (p.kind === 'slides' && p.slides) return `data:application/json;charset=utf-8,${encodeURIComponent(JSON.stringify(p.slides, null, 2))}`;
  if (p.kind === 'table' && p.sheets) return `data:text/csv;charset=utf-8,${encodeURIComponent(p.sheets[0].rows.map((r) => r.join(',')).join('\n'))}`;
  return `data:text/plain;charset=utf-8,${encodeURIComponent(`Mock 模式下无法下载二进制文件：${path}`)}`;
}
