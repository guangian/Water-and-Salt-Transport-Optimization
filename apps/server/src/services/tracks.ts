import path from 'node:path';
import fsp from 'node:fs/promises';
import { AGENT_DIR } from '../config.js';
import type { TrackId, TrackInfo } from '../types.js';
import { readJsonOrNull } from '../util/fs.js';

/**
 * 内置赛道摘要（docs/research/01_赛事情报.md §2 材料要求、§3.0 评分、§4.4 2026 命题企业）。
 * agent/tracks/<id>/track.json 存在时优先用文件内容（由规范撰写者维护），这里只是兜底。
 */
export const BUILTIN_TRACKS: TrackInfo[] = [
  {
    id: 'brand',
    name: '品牌策划赛道',
    docName: '品牌策划书',
    pageLimit: 80,
    wordMin: 0,
    scoring: { written: 60, pitch: 40 },
    elements: ['目录', '概要', '引言', '品牌现状分析', '品牌定位与核心价值', '品牌元素设计', '品牌推广与传播', '品牌资产保护和风险管理', '财务预算', '附录和参考资料'],
    groups: ['研究生组', '本科组', '高职组', '在华留学生组'],
    assignedCases: [
      { company: '蚌埠市伟伟冰米酒餐饮管理有限公司', title: '伟伟冰酒酿品牌视觉符号策划方案' },
      { company: '桂林长发小寨博物馆有限公司', title: '长发小寨品牌方案（瑶族长发文化）' },
      { company: '河北康贝尔食品有限公司', title: '马大姐电解质软糖项目策划' },
      { company: '湖南白云湖生态茶业股份有限公司', title: '城步峒茶（“世界峒茶”）品牌策划' },
      { company: '南京龙韵狮风科技有限公司', title: '中华龙狮舞品牌策划方案' },
      { company: '山东树上村庄旅游开发有限公司', title: '树上村庄品牌策划' },
      { company: '山西生旺油脂股份有限公司', title: '完香牌产品策划' },
      { company: '天士力医药集团股份有限公司', title: '品牌与 Z 世代叙事构建' },
      { company: '渭南市永信一餐饮有限公司', title: '永信一清真菜（非遗老字号）品牌策划' },
      { company: '新化县桃花源农业开发有限公司', title: '“蒙洱茶”品牌策划' },
    ],
    notes: '策划书 60（6 项×10：现状分析准确性/定位正确性/元素创意性/推广可行性/资产保护与风险管理全面性/财务预算合理性）+ 展示答辩 40（8 项×5）；A4 竖版 PDF，除封面目录外 ≤80 页；封面含作品名/团队名/指导企业名，全文不得出现院校与人员信息；命名“抽签编号-团队名称-作品名称”。',
  },
  {
    id: 'startup',
    name: '创业计划赛道',
    docName: '商业计划书',
    pageLimit: 80,
    wordMin: 3000,
    scoring: { written: 60, pitch: 40 },
    elements: ['执行摘要', '项目与产品/服务', '市场分析与行业调研', '商业模式与盈利模式', '营销策略与渠道', '竞争分析', '团队架构与组织', '财务预测与融资计划', '风险分析与应对', '附录（调查问卷、证明材料）'],
    groups: ['命题组', '创业组', '创客组', '在华留学生组'],
    assignedCases: [
      { company: '惠州市（命题组）', title: '制造业与服务业协同发展（惠州产业数字化升级）' },
      { company: '惠州市（命题组）', title: '人工智能赋能“2+1”现代产业体系' },
      { company: '惠州市（命题组）', title: '对外贸易高质量发展助力高水平对外开放' },
    ],
    notes: '商业计划书 60（创新性 20 / 商业性 30 / 团队架构 10）+ 展示答辩 40（陈述 20 / 答辩 20）；不少于 3000 字，除封面目录外 ≤80 页，中文宋体小四、英文 Times New Roman 小四；封面仅含组别、团队名称、项目名称；命名“参赛组别-学校名称-团队名称-商业计划书名称”。',
  },
  {
    id: 'digital',
    name: '数字营销赛道',
    docName: '数字营销方案',
    pageLimit: 80,
    wordMin: 0,
    scoring: { written: 60, pitch: 40 },
    elements: ['目录', '概要', '引言', '目标企业分析', '市场分析', '数字营销策略制定', '执行计划与预算', '数字营销效果', '附录和参考资料'],
    groups: ['研究生组', '本科组', '高职组'],
    assignedCases: [],
    notes: '方案 60（分析全面性与深度/用户洞察与策略定位/策略可行性与渠道匹配/数据逻辑与评估机制/执行计划与财务预算合理性/创新性与价值，各 10）+ 展示答辩 40；另需 ≤3 分钟创意短视频在 ≥3 个平台发布；命名“学校名称-团队名称-作品名称-作品形式”。',
  },
  {
    id: 'case',
    name: '案例分析赛道（流通业经营模拟 / 会计与商业管理案例）',
    docName: '案例分析报告',
    pageLimit: 50,
    wordMin: 3000,
    scoring: { written: 60, pitch: 40 },
    elements: ['目录', '概要', '引言', '案例简介', '企业内外部环境分析', '发展瓶颈分析', '结论与启示', '附录与参考资料'],
    groups: ['研究生组', '本科组', '高职组'],
    assignedCases: [{ company: '科大讯飞股份有限公司（002230）', title: '基于 2022–2024 年业绩报告的未来三年商业计划（会计与商业管理案例 2026 命题）' }],
    notes: '流通业经营模拟：书面报告 60（分析框架 20 / 对象 10 / 视角 10 / 方法 10 / 结论 10）+ 陈述答辩 40，报告 ≥3000 字、配套 PPT ≥20 页；会计与商业管理案例：商业计划书 70 + 视频展示 30，12 要素，除封面目录外 ≤50 页。',
  },
  {
    id: 'mice',
    name: '商务会奖旅游策划赛道',
    docName: '商务会奖旅游策划方案',
    pageLimit: 80,
    wordMin: 0,
    scoring: { written: 60, pitch: 40 },
    elements: ['目录', '客户公司介绍', '方案目的', '方案可行性分析', '方案主题', '具体行程设计与规划', '突发事件预案', '费用预算', '方案特色', '方案预计效果'],
    groups: ['本科组', '高职组'],
    assignedCases: [
      { company: '自选客户公司（广东省肇庆市）', title: '200–300 人四天三夜公司年会，预算 200–400 万元' },
      { company: '自选客户公司（四川省南充市）', title: '200–300 人四天三夜公司年会，预算 200–400 万元' },
      { company: '自选客户公司（河南省洛阳市）', title: '200–300 人四天三夜公司年会，预算 200–400 万元' },
    ],
    notes: '书面策划 60（内容原创性与结构完整性/合理性与逻辑性/活动安排与目标一致性/新颖与创意性/排版设计与可读性/财务预算合理性，各 10）+ 演讲答辩 40（8 项×5）；需 10 分钟展示视频与 1–3 分钟宣传视频。',
  },
];

export class TrackService {
  private cache: { ts: number; tracks: TrackInfo[] } | null = null;

  constructor(private readonly tracksDir = path.join(AGENT_DIR, 'tracks')) {}

  async list(): Promise<TrackInfo[]> {
    if (this.cache && Date.now() - this.cache.ts < 10000) return this.cache.tracks;
    const out: TrackInfo[] = [];
    for (const builtin of BUILTIN_TRACKS) {
      const fromFile = await readJsonOrNull<Partial<TrackInfo>>(path.join(this.tracksDir, builtin.id, 'track.json'));
      out.push(fromFile ? { ...builtin, ...fromFile, id: builtin.id } : builtin);
    }
    // 规范目录里新增的赛道（无内置摘要）也列出来
    const names = await fsp.readdir(this.tracksDir).catch(() => [] as string[]);
    for (const n of names) {
      if (out.some((t) => t.id === n)) continue;
      const fromFile = await readJsonOrNull<TrackInfo>(path.join(this.tracksDir, n, 'track.json'));
      if (fromFile) out.push({ ...fromFile, id: n as TrackId });
    }
    this.cache = { ts: Date.now(), tracks: out };
    return out;
  }

  async get(id: string): Promise<TrackInfo | undefined> {
    return (await this.list()).find((t) => t.id === id);
  }
}
