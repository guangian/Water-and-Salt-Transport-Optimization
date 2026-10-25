import type { ModelMap, Project, StageId, StageState } from '../types';
import { emptyStages, STAGE_IDS } from '../lib/stages';

const DEFAULT_MODELS: ModelMap = {
  main: 'deepseek-v4-pro',
  researcher: 'kimi-k3',
  finance: 'deepseek-v4-pro',
  writer: 'deepseek-v4-pro',
  pitch: 'deepseek-v4-pro',
  judgeA: 'kimi-k3',
  judgeB: 'glm-5.3',
  judgeC: 'qwen3.8-max',
  inspector: 'deepseek-v4-pro',
  fast: 'deepseek-v4-flash',
};

/** 构造按顺序完成的阶段表：前 doneCount 个 done，第 doneCount+1 个为 current 状态。 */
function stagesUpTo(
  start: string,
  minutesEach: number[],
  doneCount: number,
  current?: { status: StageState['status']; note?: string },
): Record<StageId, StageState> {
  const stages = emptyStages();
  let t = Date.parse(start);
  STAGE_IDS.forEach((id, i) => {
    if (i < doneCount) {
      const dur = (minutesEach[i] ?? 10) * 60_000;
      stages[id] = { status: 'done', startedAt: new Date(t).toISOString(), endedAt: new Date(t + dur).toISOString() };
      if (id === 'review') stages[id].round = 3;
      t += dur;
    } else if (i === doneCount && current) {
      stages[id] = { status: current.status, startedAt: new Date(t).toISOString(), note: current.note };
      if (current.status === 'failed') stages[id].endedAt = new Date(t + 4 * 60_000).toISOString();
    }
  });
  return stages;
}

export const RUNNING_PROJECT_ID = 'p_20261024_153012_ab12';
const RUNNING_START = new Date(Date.now() - 91 * 60_000).toISOString();

export const MOCK_PROJECTS: Project[] = [
  {
    id: 'p_20261018_091533_a1f3',
    name: '宋韵天青——弘宝汝瓷品牌策划书',
    teamName: '青瓷工作室',
    track: 'brand',
    group: '本科组',
    subject: {
      type: 'brand',
      name: '弘宝汝瓷',
      company: '汝州市弘宝汝瓷坊有限公司',
      isAssignedCase: false,
      caseTitle: '',
      brief:
        '弘宝汝瓷是汝州本地非遗传承品牌，以天青釉茶器与文创礼品为主，线下渠道为主、线上份额不足 3%。希望面向 25–40 岁新中产与茶文化爱好者重塑品牌定位，建立“宋韵生活方式”品牌资产，并设计 300 万元以内的年度传播预算。',
      region: '河南汝州',
      industry: '非遗工艺品 / 文创',
    },
    team: [
      { role: '队长 / 总策划', major: '市场营销', strength: '品牌传播与策划' },
      { role: '调研与数据', major: '统计学', strength: '问卷设计与数据分析' },
      { role: '创意与视觉', major: '视觉传达设计', strength: '品牌元素设计' },
      { role: '财务预算', major: '会计学', strength: '预算编制与成本控制' },
    ],
    params: { researchMode: 'synthetic', pageTarget: 60, pitchMinutes: 10, reviewRounds: 3, targetScore: 85, language: 'zh-CN', models: DEFAULT_MODELS },
    inputs: [
      { name: '弘宝汝瓷企业资料.pdf', path: 'inputs/弘宝汝瓷企业资料.pdf', size: 2_341_002, kind: 'pdf', preread: 'inputs/_preread/弘宝汝瓷企业资料.pdf.txt', prereadStatus: 'done', prereadChars: 23_410 },
      { name: '2025 销售数据.xlsx', path: 'inputs/2025 销售数据.xlsx', size: 184_320, kind: 'xlsx', preread: 'inputs/_preread/2025 销售数据.xlsx.txt', prereadStatus: 'done', prereadChars: 3_280 },
      { name: '产品图册.jpg', path: 'inputs/产品图册.jpg', size: 1_204_112, kind: 'image', prereadStatus: 'done' },
    ],
    status: 'completed',
    currentStage: null,
    stages: stagesUpTo('2026-10-18T09:16:00.000Z', [6, 38, 12, 22, 9, 61, 27, 54, 5], 9),
    session: { sessionId: 'sess_7f2a9c', lastPromptAt: '2026-10-18T13:20:00.000Z', mode: 'run' },
    usage: { inputTokens: 4_812_300, outputTokens: 612_480, cacheReadTokens: 2_903_100, turns: 486, costEstimateUsd: 6.84, durationMs: 4 * 3600_000 + 14 * 60_000 },
    review: {
      rounds: [
        { round: 1, score: 78.5, byJudge: { investor: 76.0, professor: 80.0, executive: 79.5 }, hardGatePassed: false, path: '07_评审/第1轮' },
        { round: 2, score: 84.2, byJudge: { investor: 83.0, professor: 85.5, executive: 84.1 }, hardGatePassed: true, path: '07_评审/第2轮' },
        { round: 3, score: 87.6, byJudge: { investor: 86.4, professor: 88.9, executive: 87.5 }, hardGatePassed: true, path: '07_评审/第3轮' },
      ],
      latestScore: 87.6,
    },
    artifacts: {
      book: '产出/弘宝汝瓷品牌策划书.docx',
      pdf: '产出/弘宝汝瓷品牌策划书.pdf',
      deck: '产出/路演.pptx',
      model: '03_财务/财务模型.xlsx',
      script: '06_路演/讲稿.md',
      qa: '06_路演/答辩题库.md',
      package: '08_提交包/提交包.zip',
    },
    createdAt: '2026-10-18T09:15:33.000Z',
    updatedAt: '2026-10-18T13:30:12.000Z',
  },
  {
    id: RUNNING_PROJECT_ID,
    name: '智检云眸——电子信息产业 AI 质检与预测性维护服务平台',
    teamName: '云眸智造团队',
    track: 'startup',
    group: '命题组',
    subject: {
      type: 'startup',
      name: '智检云眸',
      company: '惠州智检云眸科技有限公司（拟设立）',
      isAssignedCase: true,
      caseTitle: '2.2 电子信息产业 AI 质检与预测性维护',
      brief:
        '面向惠州仲恺高新区电子信息制造集群的中小企业，提供“视觉质检模型 + 设备预测性维护”的订阅制 SaaS 与边缘一体机。首批聚焦 PCB/SMT 产线外观缺陷检测与贴片机关键部件寿命预测，目标三年服务 300 家工厂。',
      region: '广东惠州',
      industry: '工业人工智能 / 电子信息制造服务',
    },
    team: [
      { role: '队长 / CEO', major: '工商管理', strength: '战略与融资' },
      { role: 'CTO', major: '计算机科学', strength: '工业视觉算法' },
      { role: '产品与运营', major: '工业工程', strength: '产线流程与客户成功' },
      { role: '财务', major: '金融学', strength: '财务建模与估值' },
      { role: '市场', major: '市场营销', strength: 'B2B 渠道与品牌' },
    ],
    params: { researchMode: 'synthetic', pageTarget: 55, pitchMinutes: 10, reviewRounds: 3, targetScore: 85, language: 'zh-CN', models: DEFAULT_MODELS },
    inputs: [
      { name: '惠州市十五五规划建议（节选）.pdf', path: 'inputs/惠州市十五五规划建议（节选）.pdf', size: 1_120_440, kind: 'pdf', preread: 'inputs/_preread/惠州市十五五规划建议（节选）.pdf.txt', prereadStatus: 'done', prereadChars: 18_200 },
      { name: '仲恺电子信息企业访谈纪要.docx', path: 'inputs/仲恺电子信息企业访谈纪要.docx', size: 88_210, kind: 'docx', preread: 'inputs/_preread/仲恺电子信息企业访谈纪要.docx.txt', prereadStatus: 'done', prereadChars: 6_540 },
    ],
    status: 'running',
    currentStage: 'write',
    // 演示用：前 5 个阶段共 91 分钟，恰好在页面打开时进入「撰写」阶段
    stages: stagesUpTo(RUNNING_START, [5, 41, 11, 26, 8], 5, { status: 'running' }),
    session: { sessionId: 'sess_c01d4e', lastPromptAt: RUNNING_START, mode: 'run' },
    usage: { inputTokens: 2_104_880, outputTokens: 238_120, cacheReadTokens: 1_422_000, turns: 212, costEstimateUsd: 2.97, durationMs: 91 * 60_000 },
    review: { rounds: [], latestScore: null },
    artifacts: { model: '03_财务/财务模型.xlsx' },
    createdAt: new Date(Date.parse(RUNNING_START) - 48_000).toISOString(),
    updatedAt: new Date().toISOString(),
  },
  {
    id: 'p_20261021_201145_c7d9',
    name: '一口鲜甜，渔你相约——渔惑脆肉鲩火锅数字营销方案',
    teamName: '鲜甜计划小组',
    track: 'digital',
    group: '高职组',
    subject: {
      type: 'brand',
      name: '渔惑脆肉鲩火锅',
      company: '珠海市渔惑餐饮管理有限公司',
      isAssignedCase: true,
      caseTitle: '渔惑脆肉鲩火锅数字营销',
      brief: '围绕珠海本地特色脆肉鲩火锅，制定抖音 + 小红书 + 哔哩哔哩三平台数字营销方案与 3 分钟创意短视频脚本，目标 90 天内门店到店转化提升 30%。',
      region: '广东珠海',
      industry: '餐饮连锁',
    },
    team: [
      { role: '队长 / 策略', major: '电子商务', strength: '内容运营' },
      { role: '短视频', major: '数字媒体艺术', strength: '剪辑与分镜' },
      { role: '投放与数据', major: '市场营销', strength: '投放优化' },
    ],
    params: { researchMode: 'provided', pageTarget: 45, pitchMinutes: 10, reviewRounds: 2, targetScore: 85, language: 'zh-CN', models: DEFAULT_MODELS },
    inputs: [
      { name: '门店经营数据 2025Q1-Q3.csv', path: 'inputs/门店经营数据 2025Q1-Q3.csv', size: 41_020, kind: 'csv', preread: 'inputs/_preread/门店经营数据 2025Q1-Q3.csv.txt', prereadStatus: 'done', prereadChars: 1_900 },
    ],
    status: 'stopped',
    currentStage: 'charts',
    stages: stagesUpTo('2026-10-21T20:12:00.000Z', [4, 36, 9, 14], 4, { status: 'failed', note: '用户手动停止（图表渲染进行中）' }),
    session: { sessionId: 'sess_9a1b22', lastPromptAt: '2026-10-21T20:12:00.000Z', mode: 'run' },
    usage: { inputTokens: 1_280_400, outputTokens: 142_300, cacheReadTokens: 760_000, turns: 131, costEstimateUsd: 1.62, durationMs: 68 * 60_000 },
    review: { rounds: [], latestScore: null },
    artifacts: {},
    createdAt: '2026-10-21T20:11:45.000Z',
    updatedAt: '2026-10-21T21:20:03.000Z',
  },
  {
    id: 'p_20261022_143301_e4b8',
    name: '洛阳·牡丹盛会——卓越新能源 2027 年度盛典会奖旅游策划',
    teamName: '洛水会奖工作室',
    track: 'mice',
    group: '本科组',
    subject: {
      type: 'event',
      name: '卓越新能源 2027 年度盛典',
      company: '卓越新能源科技（虚拟客户）',
      isAssignedCase: true,
      caseTitle: '洛阳四天三夜公司年会（200–300 人，预算 200–400 万元）',
      brief: '为一家新能源企业 260 人团队设计洛阳四天三夜年会：欢迎晚宴、总结大会、龙门石窟与隋唐洛阳城研学、牡丹主题颁奖晚宴，预算 320 万元。',
      region: '河南洛阳',
      industry: '会奖旅游 / 新能源',
    },
    team: [
      { role: '队长 / 行程设计', major: '会展经济与管理', strength: '行程与供应商统筹' },
      { role: '预算', major: '财务管理', strength: '报价与成本' },
      { role: '活动策划', major: '旅游管理', strength: '主题与视觉' },
    ],
    params: { researchMode: 'synthetic', pageTarget: 40, pitchMinutes: 10, reviewRounds: 3, targetScore: 85, language: 'zh-CN', models: DEFAULT_MODELS },
    inputs: [],
    status: 'failed',
    currentStage: 'research',
    stages: stagesUpTo('2026-10-22T14:34:00.000Z', [5], 1, { status: 'failed', note: '网关连续 3 次返回 502（kimi-k3），请检查 API 设置后续跑' }),
    session: { sessionId: 'sess_3c7e10', lastPromptAt: '2026-10-22T14:34:00.000Z', mode: 'run' },
    usage: { inputTokens: 96_200, outputTokens: 8_140, cacheReadTokens: 40_000, turns: 14, costEstimateUsd: 0.11, durationMs: 9 * 60_000 },
    review: { rounds: [], latestScore: null },
    artifacts: {},
    createdAt: '2026-10-22T14:33:01.000Z',
    updatedAt: '2026-10-22T14:43:50.000Z',
  },
  {
    id: 'p_20261024_101209_f0e2',
    name: '智启未来——科大讯飞未来三年经营提升商业计划书',
    teamName: '星火策略顾问组',
    track: 'case',
    group: '研究生组',
    subject: {
      type: 'case-company',
      name: '科大讯飞',
      company: '科大讯飞股份有限公司（002230）',
      isAssignedCase: true,
      caseTitle: '基于 2022–2024 年业绩报告的未来三年经营提升计划（2026 届命题）',
      brief: '以商业策略顾问身份，基于 2022–2024 年年报进行财务比率分析、行业聚类分析，提出未来三年策略计划、实施计划与综合收益表预测。',
      region: '安徽合肥',
      industry: '人工智能 / 软件与信息技术服务',
    },
    team: [
      { role: '队长 / 财务分析', major: '会计学', strength: '财务比率与估值' },
      { role: '战略', major: '企业管理', strength: '战略分析框架' },
      { role: '行业研究', major: '金融学', strength: '行业数据与聚类' },
    ],
    params: { researchMode: 'provided', pageTarget: 45, pitchMinutes: 10, reviewRounds: 3, targetScore: 85, language: 'zh-CN', models: DEFAULT_MODELS },
    inputs: [
      { name: '科大讯飞 2024 年年度报告.pdf', path: 'inputs/科大讯飞 2024 年年度报告.pdf', size: 9_820_112, kind: 'pdf', prereadStatus: 'pending' },
    ],
    status: 'idle',
    currentStage: null,
    stages: emptyStages(),
    usage: { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, turns: 0, costEstimateUsd: 0, durationMs: 0 },
    review: { rounds: [], latestScore: null },
    artifacts: {},
    createdAt: '2026-10-24T10:12:09.000Z',
    updatedAt: '2026-10-24T10:12:09.000Z',
  },
];

export { DEFAULT_MODELS };
