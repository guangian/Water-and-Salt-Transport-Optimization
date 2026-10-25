import type { JudgeId, JudgeReview, JudgeScore, MachineCheckReport, ReviewData, ReviewFix, ReviewRound, TrackId } from '../types';

export const JUDGE_META: Record<JudgeId, { name: string; persona: string; defaultModel: string }> = {
  investor: { name: '评委 A · 投资人', persona: '关注商业逻辑、财务可信度与增长路径', defaultModel: 'kimi-k3' },
  professor: { name: '评委 B · 学者', persona: '关注理论框架、调研方法与论证严谨', defaultModel: 'glm-5.3' },
  executive: { name: '评委 C · 企业高管', persona: '关注落地可行性、预算与执行细节', defaultModel: 'qwen3.8-max' },
};

const BRAND_DIMENSIONS: { dimension: string; weight: number }[] = [
  { dimension: '现状分析准确性', weight: 12 },
  { dimension: '定位与核心价值', weight: 12 },
  { dimension: '品牌元素创意性', weight: 10 },
  { dimension: '推广传播可行性', weight: 12 },
  { dimension: '资产保护与风险', weight: 8 },
  { dimension: '财务预算合理性', weight: 8 },
  { dimension: '内容完整与逻辑', weight: 10 },
  { dimension: '数据与证据质量', weight: 10 },
  { dimension: '路演材料专业度', weight: 10 },
  { dimension: '创新性', weight: 8 },
];

const STARTUP_DIMENSIONS: { dimension: string; weight: number }[] = [
  { dimension: '创新性与技术含量', weight: 20 },
  { dimension: '商业模式完整性', weight: 10 },
  { dimension: '目标市场与客户定位', weight: 8 },
  { dimension: '盈利模式与财务预测', weight: 12 },
  { dimension: '行业调研与一手资料', weight: 10 },
  { dimension: '竞争分析', weight: 10 },
  { dimension: '团队架构', weight: 10 },
  { dimension: '路演材料专业度', weight: 10 },
  { dimension: '数据与证据质量', weight: 10 },
];

export function dimensionsFor(track: TrackId): { dimension: string; weight: number }[] {
  return track === 'startup' || track === 'case' ? STARTUP_DIMENSIONS : BRAND_DIMENSIONS;
}

function scoresFromPattern(track: TrackId, pattern: number[], evidencePrefix: string): JudgeScore[] {
  return dimensionsFor(track).map((d, i) => {
    const score = pattern[i % pattern.length];
    const issues = score >= 4.5 ? [] : score >= 4 ? [`${d.dimension}：论证充分但案例可再具体`] : [`${d.dimension}：证据链不完整，缺少数值支撑`, `${d.dimension}：与评分项表述未一一对应`];
    const fixes = score >= 4.5 ? [] : score >= 4 ? [`补充 1–2 个带来源的数据点`] : [`按“结论句 + 证据表 + 含义”重写该节`, `补齐来源与年份口径`];
    return { dimension: d.dimension, weight: d.weight, score, cap: 5, evidence: `${evidencePrefix} · ${d.dimension}相关章节`, issues, fixes };
  });
}

function totalOf(scores: JudgeScore[]): number {
  return Math.round(scores.reduce((s, x) => s + (x.score / x.cap) * x.weight, 0) * 10) / 10;
}

function judge(track: TrackId, round: number, id: JudgeId, pattern: number[], verdict: string, topFixes: ReviewFix[], model?: string): JudgeReview {
  const scores = scoresFromPattern(track, pattern, `05_正文`);
  const total = totalOf(scores);
  return {
    round,
    judge: id,
    model: model ?? JUDGE_META[id].defaultModel,
    track,
    rubric: `tracks/${track}/rubric.json`,
    scores,
    total,
    officialEstimate: { written: Math.round(total * 0.6 * 10) / 10, pitch: Math.round(total * 0.4 * 10) / 10, total: Math.round(total) },
    hardGates: [
      { gate: 'anonymity', passed: true, detail: '未发现院校/姓名/联系方式' },
      { gate: 'page_limit', passed: round > 1, detail: round > 1 ? '63 / 80 页' : '估算 84 页，超出 80 页上限' },
      { gate: 'required_elements', passed: true, detail: '10/10 要素齐备' },
    ],
    topFixes,
    verdict,
  };
}

function machine(track: TrackId, round: number): MachineCheckReport {
  const r = round;
  return {
    track,
    generatedAt: `2026-10-18T1${r + 0}:20:00.000Z`,
    checks: [
      { id: 'struct_coverage', name: '结构覆盖（必备要素）', value: 1.0, threshold: { pass: 1.0, warn: 0.85 }, status: 'pass', detail: '10/10 要素命中' },
      { id: 'word_count', name: '正文字数', value: 28640 + r * 1200, threshold: { pass: 15000, warn: 10000 }, status: 'pass', detail: `${(28640 + r * 1200).toLocaleString()} 字（估算 ${60 + r} 页）` },
      { id: 'page_estimate', name: '页数估算', value: r === 1 ? 84 : 63, threshold: { pass: 80, warn: 76 }, status: r === 1 ? 'fail' : 'warn', detail: r === 1 ? '估算 84 页，超出上限' : '估算 63 页，距上限留有余量；附录 11 页' },
      { id: 'numeric_density', name: '数值密度（每千字数值数）', value: 14.2 + r, threshold: { pass: 15, warn: 10 }, status: 14.2 + r >= 15 ? 'pass' : 'warn', detail: `${(14.2 + r).toFixed(1)} / 千字` },
      { id: 'source_count', name: '参考文献条数', value: 12 + r * 4, threshold: { pass: 15, warn: 10 }, status: 12 + r * 4 >= 15 ? 'pass' : 'warn', detail: `${12 + r * 4} 条，GB/T 7714 格式` },
      { id: 'source_pairing', name: '数据与来源配对率', value: 0.86 + r * 0.04, threshold: { pass: 0.95, warn: 0.85 }, status: 0.86 + r * 0.04 >= 0.95 ? 'pass' : 'warn', detail: `${Math.round((0.86 + r * 0.04) * 100)}% 的数值句带来源` },
      { id: 'chart_count', name: '图表数量', value: 16 + r * 2, threshold: { pass: 20, warn: 14 }, status: 16 + r * 2 >= 20 ? 'pass' : 'warn', detail: `${16 + r * 2} 幅图 + 9 张表` },
      { id: 'figure_caption', name: '图注四段式', value: 1.0, threshold: { pass: 1.0, warn: 0.9 }, status: 'pass', detail: '所有图表含“图 N-M 标题 + 来源”' },
      { id: 'forbidden_terms', name: '禁语检查', value: r === 1 ? 3 : 0, threshold: { pass: 0, warn: 0 }, status: r === 1 ? 'fail' : 'pass', detail: r === 1 ? '命中 3 处：“众所周知”“业内领先”“绝对”' : '0 处命中' },
      { id: 'anonymity_hits', name: '匿名检查', value: 0, threshold: { pass: 0, warn: 0 }, status: 'pass', detail: '未发现院校/姓名/联系方式' },
      { id: 'synthetic_marking', name: '示例数据标注', value: 1.0, threshold: { pass: 1.0, warn: 0.9 }, status: 'pass', detail: '6/6 处示例数据已标注' },
      { id: 'cross_doc', name: '跨材料一致性', value: r === 3 ? 0.95 : 0.88, threshold: { pass: 1.0, warn: 0.9 }, status: r === 3 ? 'warn' : 'fail', detail: r === 3 ? '20 个关键数字中 1 个不一致' : '20 个关键数字中 3 个不一致' },
      { id: 'summary_present', name: '概要页', value: true, status: 'pass', detail: '概要 1 页，含目标/策略/预算/效果' },
      { id: 'toc_present', name: '目录', value: true, status: 'pass', detail: '自动目录域，二级' },
      { id: 'chapter_conclusions', name: '章末小结', value: 7 / 7, threshold: { pass: 1.0, warn: 0.8 }, status: 'pass', detail: '7/7 章有小结句' },
      { id: 'deck_pages', name: '路演页数', value: 12, threshold: { pass: 10, warn: 8 }, status: 'pass', detail: '12 页，含封面与结尾' },
      { id: 'deck_text_density', name: '路演单页文字量', value: 96, threshold: { pass: 120, warn: 150 }, status: 'pass', detail: '最大单页 96 字（上限 120）' },
      { id: 'script_duration', name: '讲稿时长估算', value: 9.8, threshold: { pass: 10, warn: 11 }, status: 'pass', detail: '约 2600 字 → 9.8 分钟' },
      { id: 'qa_coverage', name: '题库覆盖评分维度', value: 1.0, threshold: { pass: 1.0, warn: 0.8 }, status: 'pass', detail: '30 问覆盖全部 6 个评分项' },
      { id: 'budget_sum', name: '预算合计校验', value: true, status: 'pass', detail: '预算表合计 292 万元 = 各项之和' },
      { id: 'font_check', name: '字体检查（docx）', value: null, status: 'skipped', detail: '需 Word/WPS 环境；当前未安装 PDF 转换器' },
      { id: 'pdf_pages', name: 'PDF 实际页数', value: null, status: 'skipped', detail: '未生成 PDF（无转换器）' },
    ],
    hardGates: [
      { id: 'anonymity_hits', passed: true, hits: [] },
      { id: 'page_limit', passed: r > 1, pages: r === 1 ? 84 : 63, limit: 80 },
      { id: 'required_elements', passed: true, detail: '10/10' },
      { id: 'forbidden_terms', passed: r > 1, hits: r === 1 ? ['众所周知', '业内领先', '绝对'] : [] },
      { id: 'cover_fields', passed: true, detail: '作品名称 / 团队名称 / 指导企业名称' },
      { id: 'file_naming', passed: true, detail: '03-20-青瓷工作室-宋韵天青——弘宝汝瓷品牌策划书.pdf' },
    ],
    crossDoc: {
      consistent: false,
      mismatches:
        r === 3
          ? [{ key: 'kpi_exposure', label: '曝光量目标', book: '5000 万', deck: '5000 万', script: '4800 万' }]
          : [
              { key: 'kpi_exposure', label: '曝光量目标', book: '5000 万', deck: '4500 万', script: '4800 万' },
              { key: 'budget_total', label: '推广预算总额', book: 300, deck: 292, script: 300 },
              { key: 'sample_size', label: '问卷有效样本量', book: 412, deck: 384, script: 412 },
            ],
    },
    summary: r === 1 ? { pass: 14, warn: 3, fail: 3, skipped: 2 } : r === 2 ? { pass: 16, warn: 3, fail: 1, skipped: 2 } : { pass: 18, warn: 2, fail: 0, skipped: 2 },
  };
}

function round(track: TrackId, n: number, judges: JudgeReview[], fixes: ReviewFix[], converged: boolean): ReviewRound {
  const byJudge: Partial<Record<JudgeId, number>> = {};
  for (const j of judges) byJudge[j.judge] = j.total;
  const score = Math.round((judges.reduce((s, j) => s + j.total, 0) / judges.length) * 10) / 10;
  const m = machine(track, n);
  return {
    round: n,
    score,
    byJudge,
    hardGatePassed: m.hardGates.every((g) => g.passed),
    converged,
    fixes,
    judges,
    machine: m,
    path: `07_评审/第${n}轮`,
    at: `2026-10-18T1${n + 0}:25:00.000Z`,
  };
}

export const MOCK_REVIEW_BRAND: ReviewData = (() => {
  const track: TrackId = 'brand';
  const r1 = round(
    track,
    1,
    [
      judge(track, 1, 'investor', [4, 4, 3, 4, 3, 4, 4, 3, 4, 4], '商业逻辑成立，但预算 ROI 与转化链条缺少量化，第 7 章预算须与第 6 章传播动作一一对应。', [
        { priority: 1, where: '07_财务预算.md §7.2', what: '补预算—动作—KPI 对照表', gain: '+3' },
        { priority: 2, where: '02_品牌现状分析.md §2.3', what: '补问卷样本结构', gain: '+2' },
        { priority: 3, where: '全文', what: '清除 3 处禁语', gain: '+1' },
      ]),
      judge(track, 1, 'professor', [4, 4, 4, 4, 4, 4, 4, 4, 4, 4], '框架完整，STP 逻辑清楚；调研方法说明不足，应补充抽样与问卷信效度。', [
        { priority: 1, where: '02_品牌现状分析.md §2.3', what: '补抽样方法与样本结构表', gain: '+3' },
        { priority: 2, where: '01_调研/参考文献.json', what: '来源增至 15 条以上', gain: '+2' },
      ]),
      judge(track, 1, 'executive', [4, 4, 4, 4, 4, 4, 4, 3, 4, 4], '传播方案可执行，但页数超限、部分数据无来源，提交前必须压缩至 80 页以内。', [
        { priority: 1, where: '05_正文/A1_调研问卷.md', what: '附录问卷压缩为 4 页', gain: '+2（合规）' },
        { priority: 2, where: '06_品牌推广与传播.md §6.4', what: '补刊例来源', gain: '+2' },
      ]),
    ],
    [
      { priority: 1, where: '05_正文/A1_调研问卷.md', what: '压缩附录至 80 页以内（硬门槛）', gain: '合规' },
      { priority: 2, where: '02_品牌现状分析.md §2.3', what: '补 412 份问卷的样本结构表', gain: '+3' },
      { priority: 3, where: '07_财务预算.md §7.2', what: '补预算—动作—KPI 对照表', gain: '+3' },
      { priority: 4, where: '全文', what: '清除“众所周知/业内领先/绝对”等禁语', gain: '+1' },
    ],
    false,
  );
  const r2 = round(
    track,
    2,
    [
      judge(track, 2, 'investor', [4, 4.5, 4, 4, 4, 4.5, 4, 4, 4, 4], '预算与 KPI 对照清晰，ROI 1.6 可信；建议补充悲观情景。', [
        { priority: 1, where: '07_财务预算.md §7.4', what: '补三情景预算效果', gain: '+2' },
      ]),
      judge(track, 2, 'professor', [4.5, 4, 4, 4, 4, 4, 4.5, 4, 4, 4.5], '样本结构表补齐后调研可信度明显提升；参考文献 20 条、格式规范。', [
        { priority: 1, where: '03_品牌定位与核心价值.md §3.2', what: '定位图加入竞品坐标来源', gain: '+1' },
      ]),
      judge(track, 2, 'executive', [4, 4, 4, 4.5, 4, 4, 4, 4, 4.5, 4], '页数达标；传播排期与预算一致；建议把快闪店的落地供应商与报价写入附录。', [
        { priority: 1, where: '06_品牌推广与传播.md §6.5', what: '补快闪店报价与供应商', gain: '+1' },
      ]),
    ],
    [
      { priority: 1, where: '07_财务预算.md §7.4', what: '补三情景预算效果', gain: '+2' },
      { priority: 2, where: '讲稿.md / slides.json', what: '统一预算总额口径（292 万 vs 300 万）', gain: '+1' },
    ],
    false,
  );
  const r3 = round(
    track,
    3,
    [
      judge(track, 3, 'investor', [4.5, 4.5, 4, 4.5, 4, 4.5, 4.5, 4, 4.5, 4.5], '接近一等奖水平：财务与传播闭环完整，三情景分析到位；讲稿曝光口径与正文不一致是唯一硬伤。', [
        { priority: 1, where: '06_路演/讲稿.md P9', what: '曝光量统一为 5000 万', gain: '+1' },
        { priority: 2, where: '07_财务预算.md', what: '补预算执行的月度现金流', gain: '+1' },
        { priority: 3, where: '08_品牌资产保护.md', what: '补商标近似查询结果', gain: '+0.5' },
      ]),
      judge(track, 3, 'professor', [4.5, 4.5, 4.5, 4, 4.5, 4, 4.5, 4.5, 4.5, 4.5], '论证严谨、数据来源完整，品牌元素设计与定位高度一致，可作为范例。', [
        { priority: 1, where: '04_品牌元素设计.md §4.3', what: '增加 ICC 准则条款引用', gain: '+0.5' },
      ]),
      judge(track, 3, 'executive', [4.5, 4, 4.5, 4.5, 4, 4.5, 4.5, 4, 4.5, 4.5], '可直接落地执行；建议补应急预案的触发阈值与负责人分工。', [
        { priority: 1, where: '06_品牌资产保护与风险管理.md §6.3', what: '补应急预案触发阈值', gain: '+1' },
      ]),
    ],
    [
      { priority: 1, where: '06_路演/讲稿.md P9', what: '曝光量统一为 5000 万（跨材料一致性）', gain: '+1' },
      { priority: 2, where: '06_品牌资产保护与风险管理.md §6.3', what: '补应急预案触发阈值与分工', gain: '+1' },
      { priority: 3, where: '07_财务预算.md', what: '补预算执行的月度现金流', gain: '+1' },
    ],
    true,
  );
  return { latest: r3, rounds: [r1, r2, r3], machine: r3.machine ?? null };
})();

/** 为运行中的项目在评审阶段动态生成一轮评审结果 */
export function makeReviewRound(track: TrackId, n: number, base: number): ReviewRound {
  const p = (delta: number) => dimensionsFor(track).map((_, i) => Math.min(5, Math.max(3, Math.round(((base + delta) / 100) * 5 * 2 + (i % 3 === 0 ? -0.5 : 0)) / 2)));
  const judges = [
    judge(track, n, 'investor', p(-2), '财务模型三表勾稽正确，盈亏平衡在第 19 个月；客户获取成本假设偏乐观，建议补充试点客户的实际数据。', [
      { priority: 1, where: '05_正文/05_财务预测.md §5.3', what: '补 CAC 的试点数据来源', gain: '+2' },
      { priority: 2, where: '05_正文/04_商业模式.md', what: '补订阅续费率敏感性', gain: '+1' },
    ]),
    judge(track, n, 'professor', p(1), '选题契合惠州“人工智能+”命题方向，调研访谈 12 家企业形成一手资料；建议补技术路线的对比表。', [
      { priority: 1, where: '05_正文/01_项目与产品.md §1.4', what: '补技术路线对比表', gain: '+2' },
    ]),
    judge(track, n, 'executive', p(0), '落地路径清晰，仲恺首批 30 家工厂的推进计划可执行；建议补数据安全与产线停机风险的应对。', [
      { priority: 1, where: '05_正文/08_风险分析.md', what: '补数据安全与停机风险应对', gain: '+1' },
    ]),
  ];
  return round(track, n, judges, [
    { priority: 1, where: '05_正文/05_财务预测.md §5.3', what: '补 CAC 的试点数据来源', gain: '+2' },
    { priority: 2, where: '05_正文/01_项目与产品.md §1.4', what: '补技术路线对比表', gain: '+2' },
    { priority: 3, where: '05_正文/08_风险分析.md', what: '补数据安全与停机风险应对', gain: '+1' },
  ], false);
}
