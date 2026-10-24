import type { AgentEvent, AgentRef, Project, StageId, StageStatus } from '../types';
import { STAGE_IDS, STAGES } from '../lib/stages';

type DistributiveOmit<T, K extends PropertyKey> = T extends unknown ? Omit<T, K> : never;
export type ScriptEvent = DistributiveOmit<AgentEvent, 'seq' | 'ts' | 'projectId'>;

export interface Step {
  /** 与上一条事件的间隔（毫秒），默认 300–800 随机 */
  delay?: number;
  ev: ScriptEvent;
}

const MAIN: AgentRef = { id: 'main', name: '主代理' };

let toolCounter = 0;
function tid(prefix = 'toolu'): string {
  toolCounter += 1;
  return `${prefix}_${toolCounter.toString(36).padStart(4, '0')}${Math.random().toString(36).slice(2, 6)}`;
}

function rand(min: number, max: number): number {
  return Math.round(min + Math.random() * (max - min));
}

// ---------- 原子事件 ----------

function ev<T extends ScriptEvent['type']>(type: T, agent: AgentRef, data: Extract<ScriptEvent, { type: T }>['data'], delay?: number): Step {
  return { delay, ev: { type, agent, data } as ScriptEvent };
}

/** 把一段文本拆成 3–6 个增量再补一条完整 text */
function text(agent: AgentRef, content: string, chunks = 4): Step[] {
  const steps: Step[] = [];
  const n = Math.max(1, Math.min(chunks, Math.ceil(content.length / 18)));
  const size = Math.ceil(content.length / n);
  for (let i = 0; i < n; i++) {
    steps.push(ev('text_delta', agent, { text: content.slice(i * size, (i + 1) * size) }, rand(300, 600)));
  }
  steps.push(ev('text', agent, { text: content }, rand(300, 500)));
  return steps;
}

function thinking(agent: AgentRef, content: string): Step[] {
  return [ev('thinking', agent, { text: content, delta: false }, rand(400, 800))];
}

function tool(agent: AgentRef, name: string, input: unknown, output: string, isError = false): Step[] {
  const toolUseId = tid();
  return [
    ev('tool_use', agent, { toolUseId, name, input }, rand(300, 700)),
    ev('tool_result', agent, { toolUseId, name, isError, output }, rand(400, 800)),
  ];
}

function file(agent: AgentRef, path: string, op: 'write' | 'edit' | 'delete', size: number): Step {
  return ev('file', agent, { path, op, size }, rand(300, 500));
}

function stage(id: StageId, status: StageStatus, note?: string, round?: number): Step {
  return ev('stage', MAIN, { stage: id, status, note, round }, rand(300, 600));
}

function progress(id: StageId, status: 'running' | 'done' | 'failed', note?: string, round?: number): Step[] {
  return [
    ...tool(MAIN, 'mcp__bizgold__update_progress', { stage: id, status, note, round }, `{"ok":true,"stage":"${id}","status":"${status}"}`),
    stage(id, status, note, round),
  ];
}

function usage(scale: number): Step {
  return ev(
    'usage',
    MAIN,
    {
      inputTokens: Math.round(180_000 * scale),
      outputTokens: Math.round(22_000 * scale),
      cacheReadTokens: Math.round(110_000 * scale),
      turns: Math.round(18 * scale),
      costEstimateUsd: Math.round(0.26 * scale * 100) / 100,
    },
    300,
  );
}

function checklist(done: number, running: number | null, round?: string): string {
  return STAGES.map((s, i) => {
    const mark = i < done ? '✅' : i === running ? '🔄' : '⏳';
    const name = s.id === 'review' && round ? `评审（第 ${round} 轮）` : s.name;
    return `${i + 1}. ${mark} ${name}`;
  }).join('\n');
}

interface SubagentSpec {
  name: string;
  subagentType: string;
  description: string;
}

/** 子代理生命周期：Agent 工具调用 → subagent start → 子代理事件 → subagent end → Agent 工具结果 */
function subagent(spec: SubagentSpec, body: (ref: AgentRef) => Step[], summary: string): { start: Step[]; body: Step[]; end: Step[] } {
  const toolUseId = tid('agent');
  const ref: AgentRef = { id: toolUseId, name: spec.name };
  return {
    start: [
      ev('tool_use', MAIN, { toolUseId, name: 'Agent', input: { subagent_type: spec.subagentType, description: spec.description, prompt: `见 任务/${spec.subagentType}.md` } }, rand(300, 600)),
      ev('subagent', MAIN, { toolUseId, name: spec.name, status: 'start' }, 300),
    ],
    body: body(ref),
    end: [
      ev('subagent', MAIN, { toolUseId, name: spec.name, status: 'end', summary }, rand(300, 500)),
      ev('tool_result', MAIN, { toolUseId, name: 'Agent', isError: false, output: summary }, 300),
    ],
  };
}

/** 交错两个子代理的事件流，模拟并行 */
function interleave(a: Step[], b: Step[]): Step[] {
  const out: Step[] = [];
  let i = 0;
  let j = 0;
  while (i < a.length || j < b.length) {
    if (i < a.length) out.push(a[i++]);
    if (i < a.length && Math.random() < 0.4) out.push(a[i++]);
    if (j < b.length) out.push(b[j++]);
  }
  return out;
}

// ---------- 阶段脚本 ----------

interface Ctx {
  p: Project;
  S: string; // 策划对象
  trackName: string;
  docName: string;
  outName: string;
  trackId: string;
}

function ctxOf(p: Project): Ctx {
  const names: Record<string, [string, string]> = {
    brand: ['品牌策划赛道', '品牌策划书'],
    startup: ['创业计划赛道', '商业计划书'],
    digital: ['数字营销赛道', '数字营销方案'],
    case: ['会计与商业管理案例竞赛', '案例分析商业计划书'],
    mice: ['商务会奖旅游策划竞赛', '会奖旅游策划书'],
  };
  const [trackName, docName] = names[p.track] ?? ['赛道', '策划书'];
  return { p, S: p.subject.name || p.name, trackName, docName, outName: `${p.subject.name || '作品'}${docName}`, trackId: p.track };
}

function parseStage(c: Ctx): Step[] {
  return [
    ...progress('parse', 'running'),
    ...text(MAIN, `开始读题解析。先阅读任务书与 ${c.trackName} 规范，建立「评分项 → 章节 → 关键证据」映射表，再列出材料清单与硬约束。`),
    ...tool(MAIN, 'Read', { file_path: '任务书.md' }, `# 任务书\n- 赛道：${c.trackName}（${c.trackId}） | 组别：${c.p.group}\n- 作品名：${c.p.name}\n- 策划对象：${c.S}（${c.p.subject.company}）\n…（共 42 行）`),
    ...tool(MAIN, 'Read', { file_path: `tracks/${c.trackId}/spec.md` }, `# ${c.trackName}规范\n## 章节大纲与篇幅配额\n| 章 | 标题 | 页数 | 必答问题 | 对应评分项 |\n…（共 188 行）`),
    ...tool(MAIN, 'Read', { file_path: `tracks/${c.trackId}/rubric.json` }, `{"track":"${c.trackId}","dimensions":[{"id":"status","name":"现状分析准确性","weight":12,…}],"targetScore":85,"minDimension":3}`),
    ...tool(MAIN, 'Glob', { pattern: 'inputs/_preread/*.txt' }, c.p.inputs.map((i) => i.preread ?? i.path).join('\n') || '（无用户资料）'),
    ...thinking(MAIN, `评分项与章节的对应关系：现状分析 ↔ 第 2 章，需要外部/内部双视角与一手数据；定位 ↔ 第 3 章，必须给出竞争对手识别与异同点；预算 ↔ 第 7 章，须与传播动作一一对应。硬约束：匿名、≤80 页、引用标注出处。用户资料已预读 ${c.p.inputs.length} 份，其中企业经营数据可直接作为内部环境分析的证据。`),
    ...tool(MAIN, 'Write', { file_path: '00_任务解析/任务解析.md', content: `# 任务解析\n\n## 1. 命题解读\n${c.S}……` }, '文件已写入（8.4 KB）'),
    file(MAIN, '00_任务解析/任务解析.md', 'write', 8_420),
    ...text(MAIN, `任务解析完成：评分映射 6 项 → 7 章；材料清单 8 项（${c.docName}、财务模型、路演 PPT、讲稿、题库、调研工具包、声明、提交包）；约束 5 条。\n\n${checklist(0, 0)}`),
    ...progress('parse', 'done', '评分映射 6→7，材料清单 8 项'),
    usage(0.4),
  ];
}

function researchStage(c: Ctx): Step[] {
  const r1 = subagent(
    { name: '行业研究员 #1', subagentType: 'industry-researcher', description: '行业研究员 #1 —— 市场规模与趋势' },
    (ref) => [
      ...text(ref, `任务：测算 ${c.S} 所在行业 2019—2025 年市场规模与渠道结构，每条数据带来源。先检索行业协会与统计年鉴。`),
      ...tool(ref, 'mcp__bizgold__web_search', { query: `${c.p.subject.industry} 市场规模 2025 报告`, count: 8 }, `[{"title":"2025 年行业运行报告","url":"https://www.example-association.org.cn/report/2025","snippet":"2025 年市场规模预计达 33.5 亿元，同比增长 13.6%……","date":"2026-01-15"},{"title":"行业统计年鉴 2025","url":"https://data.stats.gov.cn/","snippet":"……"}]（共 8 条）`),
      ...tool(ref, 'mcp__bizgold__web_fetch', { url: 'https://www.example-association.org.cn/report/2025', max_chars: 6000 }, `{"title":"2025 年行业运行报告","content":"## 一、总体规模\\n2025 年市场规模预计 33.5 亿元，其中线上渠道 14.9 亿元，占比 44.5%……（已抽取 5,820 字）","fetched_at":"2026-10-24T15:40:12Z"}`),
      ...text(ref, `核心数据：2025E 市场规模 33.5 亿元（+13.6%），线上占比由 2019 年 21% 升至 44.5%，增量的 71% 来自线上。竞品价格带、渠道结构已整理为证据表。`),
      ...tool(ref, 'Write', { file_path: '01_调研/01_市场规模与趋势.md', content: '# 市场规模与趋势\n\n## 结论摘要\n…' }, '文件已写入（14.2 KB）'),
      file(ref, '01_调研/01_市场规模与趋势.md', 'write', 14_210),
      ...tool(ref, 'Write', { file_path: '01_调研/数据/市场规模.csv', content: 'year,offline,online,source\n2021,13.2,6.8,…' }, '文件已写入（1.2 KB）'),
      file(ref, '01_调研/数据/市场规模.csv', 'write', 1_240),
      ...tool(ref, 'Edit', { file_path: '01_调研/参考文献.json', old_string: ']', new_string: ',{"id":7,…}]' }, '已追加 7 条参考文献'),
      file(ref, '01_调研/参考文献.json', 'edit', 6_410),
    ],
    `完成市场规模测算：2025E 33.5 亿元，线上占比 44.5%；输出 01_市场规模与趋势.md、市场规模.csv，新增 7 条参考文献。`,
  );
  const r2 = subagent(
    { name: '消费者洞察', subagentType: 'consumer-insight', description: '消费者洞察与调研设计' },
    (ref) => [
      ...text(ref, `任务：刻画 ${c.S} 的目标人群画像，设计问卷与访谈提纲，并在「示例数据」模式下生成标注为示例的模拟结果。`),
      ...tool(ref, 'mcp__bizgold__web_search', { query: `${c.p.subject.industry} 消费者 画像 2025 年轻人 偏好`, count: 6 }, `[{"title":"2025 新中产消费趋势报告","url":"https://report.example.com/2025","snippet":"25–40 岁人群对文化叙事类产品的溢价意愿达 68%……"}]（共 6 条）`),
      ...text(ref, `目标人群：25–40 岁新中产，一二线城市占 63%，首次触达渠道以短视频为主（58%）。问卷 28 题、访谈提纲 12 问已设计；示例结果（n=412）全部标注“【示例数据，须替换为真实调研结果】”。`),
      ...tool(ref, 'Write', { file_path: '01_调研/03_消费者洞察.md', content: '# 消费者洞察\n…' }, '文件已写入（11.3 KB）'),
      file(ref, '01_调研/03_消费者洞察.md', 'write', 11_300),
      ...tool(ref, 'Write', { file_path: '01_调研/调研工具包/问卷.md', content: '# 问卷\n…' }, '文件已写入（7.8 KB）'),
      file(ref, '01_调研/调研工具包/问卷.md', 'write', 7_800),
      ...tool(ref, 'Write', { file_path: '01_调研/调研工具包/示例结果.md', content: '# 示例结果（示例数据，须替换）\n…' }, '文件已写入（5.6 KB）'),
      file(ref, '01_调研/调研工具包/示例结果.md', 'write', 5_600),
    ],
    `完成消费者洞察：人群画像 3 类、问卷 28 题、访谈提纲 12 问、示例结果 n=412（已标注）。`,
  );
  return [
    ...progress('research', 'running'),
    ...text(MAIN, `调研阶段：并行派出行业研究员与消费者洞察研究员，任务书已写入 任务/ 目录。要求每条数据带来源，检索次数 ≤12。`),
    ...tool(MAIN, 'Write', { file_path: '任务/industry-researcher_1.md', content: '# 子任务：行业研究员 #1 —— 市场规模与趋势\n…' }, '文件已写入'),
    ...r1.start,
    ...r2.start,
    ...interleave(r1.body, r2.body),
    ...r1.end,
    ...r2.end,
    ...text(MAIN, `调研完成：3 份调研笔记、2 个数据 CSV、18 条参考文献；示例数据 6 处已标注并登记。\n\n${checklist(1, 1)}`),
    ...progress('research', 'done', '3 份笔记 / 18 条参考文献'),
    usage(1.6),
  ];
}

function strategyStage(c: Ctx): Step[] {
  return [
    ...progress('strategy', 'running'),
    ...text(MAIN, `策略阶段：基于调研结论推导定位逻辑链（趋势 → 洞察 → 竞争空隙 → 定位 → 核心价值），并写出每章结论句。`),
    ...tool(MAIN, 'Read', { file_path: '01_调研/01_市场规模与趋势.md' }, '# 市场规模与趋势\n## 结论摘要\n…（共 212 行）'),
    ...tool(MAIN, 'Read', { file_path: '01_调研/03_消费者洞察.md' }, '# 消费者洞察\n…（共 168 行）'),
    ...thinking(MAIN, `竞争空隙：主要竞品要么走高价收藏路线，要么走低价日用路线，“有文化叙事的中价位生活方式”是空白。${c.S} 的工艺与传承资产可以转译为这一定位。Big Idea 需要同时支撑视觉、内容与活动。`),
    ...text(MAIN, `定位：面向 25–40 岁新中产的生活方式品牌；核心价值：传承可感知、审美可日常、价格可触达。Big Idea 与各章结论句已写入策略框架。`),
    ...tool(MAIN, 'Write', { file_path: '02_策略/策略框架.md', content: '# 策略框架\n\n## 1. 定位逻辑链\n…' }, '文件已写入（9.9 KB）'),
    file(MAIN, '02_策略/策略框架.md', 'write', 9_950),
    ...progress('strategy', 'done', '定位逻辑链 + Big Idea'),
    usage(2.0),
  ];
}

function financeStage(c: Ctx): Step[] {
  const f = subagent(
    { name: '财务建模师', subagentType: 'financial-modeler', description: '财务建模师 —— assumptions → 模型 + numbers.json' },
    (ref) => [
      ...text(ref, `任务：为 ${c.S} 建立${c.trackId === 'startup' || c.trackId === 'case' ? '五年三表财务模型（含盈亏平衡、敏感性与估值）' : '分阶段预算表与效果测算'}，所有假设标注来源。`),
      ...tool(ref, 'Write', { file_path: '03_财务/assumptions.json', content: `{"mode":"${c.trackId === 'startup' || c.trackId === 'case' ? 'startup' : 'budget'}",…}` }, '文件已写入（3.1 KB）'),
      file(ref, '03_财务/assumptions.json', 'write', 3_110),
      ...tool(ref, 'Bash', { command: 'python -m bizkit.finmodel 03_财务/assumptions.json', description: '渲染财务模型' }, `✅ 03_财务/财务模型.xlsx 已生成（13 个 sheet，公式 412 个）\n✅ Checks 页：0 个错误\n✅ 03_财务/numbers.json：18 条关键数字\n✅ 03_财务/财务摘要.md`),
      file(ref, '03_财务/财务模型.xlsx', 'write', 61_440),
      file(ref, '03_财务/numbers.json', 'write', 2_420),
      file(ref, '03_财务/财务摘要.md', 'write', 4_100),
      ...tool(ref, 'mcp__bizgold__record_numbers', { entries: [{ key: 'budget_total', label: '推广预算总额', value: 292, unit: '万元', source: '03_财务/财务模型.xlsx!预算表!F10' }] }, '{"ok":true,"merged":18}'),
      ...text(ref, `模型完成：Checks=0；关键数字 18 条已登记到数字账本，正文与路演须统一引用。`),
    ],
    `财务模型完成：财务模型.xlsx（13 sheet，Checks=0）、numbers.json 18 条、财务摘要.md。`,
  );
  return [...progress('finance', 'running'), ...f.start, ...f.body, ...f.end, ...text(MAIN, `财务阶段完成，数字账本 18 条。\n\n${checklist(3, 3)}`), ...progress('finance', 'done', 'Checks=0，账本 18 条'), usage(2.6)];
}

function chartsStage(c: Ctx): Step[] {
  return [
    ...progress('charts', 'running'),
    ...text(MAIN, `图表阶段：根据各章结论句与数据文件编写 charts.json（18 幅），用工具包统一渲染 300 dpi。`),
    ...tool(MAIN, 'Write', { file_path: '04_图表/charts.json', content: '{"outputDir":"04_图表","dpi":300,"theme":"bizgold","charts":[…18 项]}' }, '文件已写入（9.2 KB）'),
    file(MAIN, '04_图表/charts.json', 'write', 9_200),
    ...tool(MAIN, 'Bash', { command: 'python -m bizkit.charts 04_图表/charts.json', description: '渲染图表' }, `渲染 18 幅图表 → 04_图表/\n  fig_market_size.png (stacked_bar) ✓\n  fig_competitor_price.png (grouped_bar) ✓\n  fig_budget_split.png (hbar) ✓\n  …\n字体：Microsoft YaHei；来源脚注已统一\n索引：04_图表/index.json`),
    file(MAIN, '04_图表/fig_market_size.png', 'write', 182_400),
    file(MAIN, '04_图表/fig_competitor_price.png', 'write', 164_200),
    file(MAIN, '04_图表/fig_budget_split.png', 'write', 150_900),
    ...text(MAIN, `图表渲染完成 18/18，${c.S} 的所有图表统一配色与来源脚注。`),
    ...progress('charts', 'done', '18 幅图表'),
    usage(2.9),
  ];
}

function writeStage(c: Ctx): Step[] {
  const chapters: [string, number][] = [
    ['05_正文/00_概要.md', 2_600],
    ['05_正文/01_引言.md', 4_300],
    ['05_正文/02_现状分析.md', 16_800],
    ['05_正文/03_定位与核心价值.md', 15_200],
    ['05_正文/04_方案设计.md', 12_400],
    ['05_正文/05_推广与执行.md', 18_900],
    ['05_正文/06_风险管理.md', 9_700],
    ['05_正文/07_财务预算.md', 8_100],
  ];
  const w = subagent(
    { name: '主笔', subagentType: 'chief-writer', description: '主笔 —— 按大纲撰写各章' },
    (ref) => {
      const steps: Step[] = [...text(ref, `开始撰写《${c.outName}》。按 spec 大纲逐章写作：结论先行、数值四要素、引用图表与账本数字、章末小结。`)];
      chapters.forEach(([path, size], i) => {
        steps.push(...tool(ref, 'Write', { file_path: path, content: `# ${path.split('/').pop()?.replace(/^\d+_|\.md$/g, '')}\n\n**结论先行**：…` }, `文件已写入（${(size / 1024).toFixed(1)} KB）`));
        steps.push(file(ref, path, 'write', size));
        if (i === 2 || i === 5) {
          steps.push(...tool(ref, 'Bash', { command: `python -m bizkit.checks --chapter ${path} --track ${c.trackId}`, description: '章节自检' }, `数值密度 16.8/千字 ✓  来源配对 94% ✓  禁语 0 ✓  字数 ${Math.round(size / 2.1)}`));
          steps.push(...text(ref, `第 ${i} 章自检通过：数值密度 16.8/千字，来源配对 94%，无禁语。继续下一章。`, 2));
        }
      });
      steps.push(...tool(ref, 'Write', { file_path: '05_正文/book.json', content: `{"track":"${c.trackId}","docType":"${c.docName}","cover":{…},"chapters":[…8 章],"output":"产出/${c.outName}.docx","pdf":true}` }, '文件已写入（2.3 KB）'));
      steps.push(file(ref, '05_正文/book.json', 'write', 2_300));
      steps.push(...tool(ref, 'Bash', { command: 'python -m bizkit.docx_render 05_正文/book.json', description: '渲染 Word' }, `✅ 产出/${c.outName}.docx（估算 61 页，正文 28,640 字，图 18 / 表 9）\n✅ 匿名清洗：0 处命中\n⚠ PDF 转换器未找到（Word/WPS/LibreOffice），已跳过 PDF`));
      steps.push(file(ref, `产出/${c.outName}.docx`, 'write', 2_940_100));
      steps.push(...text(ref, `《${c.outName}》初稿完成：8 章 + 附录，估算 61 页；docx 已生成。`));
      return steps;
    },
    `初稿完成：8 章 + 附录，估算 61 页；产出/${c.outName}.docx 已生成，PDF 因无转换器跳过。`,
  );
  return [...progress('write', 'running'), ...w.start, ...w.body, ...w.end, ...text(MAIN, `撰写阶段完成。\n\n${checklist(5, 5)}`), ...progress('write', 'done', '61 页 docx'), usage(4.1)];
}

function pitchStage(c: Ctx): Step[] {
  const d = subagent(
    { name: '路演设计师', subagentType: 'pitch-designer', description: '路演设计师 —— slides.json / 讲稿 / 题库' },
    (ref) => [
      ...text(ref, `任务：${c.p.params.pitchMinutes} 分钟路演，12 页；结论句标题，每页 ≤120 字；讲稿按页配时；题库覆盖全部评分维度。`),
      ...tool(ref, 'Write', { file_path: '06_路演/slides.json', content: '{"title":"…","aspect":"16:9","theme":"bizgold","slides":[…12 页]}' }, '文件已写入（11.9 KB）'),
      file(ref, '06_路演/slides.json', 'write', 11_900),
      ...tool(ref, 'Bash', { command: 'python -m bizkit.pptx_render 06_路演/slides.json', description: '渲染 PPT' }, `✅ 产出/路演.pptx（12 页，16:9，主题 bizgold）\n   页型：cover/agenda/section/bullets/image/two_column/chart/table/kpi/timeline/quote/closing\n✅ 备注已写入；匿名检查 0 处`),
      file(ref, '产出/路演.pptx', 'write', 4_210_300),
      ...tool(ref, 'Write', { file_path: '06_路演/讲稿.md', content: `# 路演讲稿（${c.p.params.pitchMinutes} 分钟）\n…` }, '文件已写入（9.6 KB）'),
      file(ref, '06_路演/讲稿.md', 'write', 9_600),
      ...tool(ref, 'Write', { file_path: '06_路演/答辩题库.md', content: '# 答辩题库（30 问）\n…' }, '文件已写入（14.3 KB）'),
      file(ref, '06_路演/答辩题库.md', 'write', 14_300),
      ...text(ref, `路演材料完成：12 页 PPT、约 2600 字讲稿（9.8 分钟）、30 问题库。`),
    ],
    `路演完成：路演.pptx 12 页、讲稿 9.8 分钟、题库 30 问。`,
  );
  return [...progress('pitch', 'running'), ...d.start, ...d.body, ...d.end, ...progress('pitch', 'done', '12 页 / 讲稿 / 题库'), usage(4.9)];
}

function judgeBody(ref: AgentRef, judge: string, score: number, model: string, verdict: string, round: number): Step[] {
  return [
    ...text(ref, `（${model}）按量规逐维评分，只读 05_正文、产出 与 06_路演。`, 2),
    ...tool(ref, 'Read', { file_path: '05_正文/02_现状分析.md' }, '# 现状分析\n…（共 310 行）'),
    ...tool(ref, 'Read', { file_path: '06_路演/讲稿.md' }, '# 路演讲稿\n…（共 96 行）'),
    ...text(ref, `总分 ${score.toFixed(1)}。${verdict}`),
    ...tool(ref, 'Write', { file_path: `07_评审/第${round}轮/评审_${judge}.json`, content: `{"round":${round},"judge":"${judge}","total":${score},…}` }, '文件已写入（6.1 KB）'),
    file(ref, `07_评审/第${round}轮/评审_${judge}.json`, 'write', 6_100),
  ];
}

function reviewRound(c: Ctx, round: number, scores: [number, number, number], passed: boolean): Step[] {
  const [a, b, d] = scores;
  const avg = Math.round(((a + b + d) / 3) * 10) / 10;
  const j1 = subagent({ name: '评委 A · 投资人', subagentType: 'judge-investor', description: '评委 A（投资人）' }, (ref) => judgeBody(ref, 'investor', a, c.p.params.models.judgeA ?? 'kimi-k3', passed ? '财务闭环完整，三情景分析到位。' : '预算与 KPI 之间缺少量化对照；页数估算超限。', round), `评委 A 总分 ${a.toFixed(1)}`);
  const j2 = subagent({ name: '评委 B · 学者', subagentType: 'judge-professor', description: '评委 B（学者）' }, (ref) => judgeBody(ref, 'professor', b, c.p.params.models.judgeB ?? 'glm-5.3', passed ? '论证严谨、来源完整。' : '调研方法说明不足，应补抽样与样本结构。', round), `评委 B 总分 ${b.toFixed(1)}`);
  const j3 = subagent({ name: '评委 C · 企业高管', subagentType: 'judge-executive', description: '评委 C（企业高管）' }, (ref) => judgeBody(ref, 'executive', d, c.p.params.models.judgeC ?? 'qwen3.8-max', passed ? '可直接落地执行。' : '传播方案可执行，但部分数据无来源。', round), `评委 C 总分 ${d.toFixed(1)}`);
  const q = subagent(
    { name: '质检员', subagentType: 'qa-inspector', description: '质检员 —— 机器检查与跨材料一致性' },
    (ref) => [
      ...tool(ref, 'Bash', { command: `python -m bizkit.checks --track ${c.trackId} --workspace . --round ${round}`, description: '机器检查' }, passed ? `检查 22 项：pass 18 / warn 2 / fail 0 / skipped 2\n硬门槛 6/6 通过\n跨材料一致性：20 个关键数字，1 个不一致（kpi_exposure：讲稿 4800 万 vs 正文 5000 万）` : `检查 22 项：pass 14 / warn 3 / fail 3 / skipped 2\n硬门槛 4/6 通过（page_limit ✗ 84>80；forbidden_terms ✗ 3 处）\n跨材料一致性：3 个不一致`),
      ...tool(ref, 'Write', { file_path: `07_评审/第${round}轮/机器检查.json`, content: '{…}' }, '文件已写入（5.9 KB）'),
      file(ref, `07_评审/第${round}轮/机器检查.json`, 'write', 5_900),
    ],
    passed ? '机器检查 18/22 通过，硬门槛全过，1 处跨材料不一致' : '机器检查 3 项失败，硬门槛 2 项未过',
  );
  return [
    ...progress('review', 'running', `第 ${round} 轮评审`, round),
    ...text(MAIN, `第 ${round} 轮评审：三位评委（${c.p.params.models.judgeA ?? 'kimi-k3'} / ${c.p.params.models.judgeB ?? 'glm-5.3'} / ${c.p.params.models.judgeC ?? 'qwen3.8-max'}）独立打分，质检员并行运行机器检查。`),
    ...j1.start,
    ...j2.start,
    ...j3.start,
    ...q.start,
    ...interleave(interleave(j1.body, j2.body), interleave(j3.body, q.body)),
    ...j1.end,
    ...j2.end,
    ...j3.end,
    ...q.end,
    ev('review', MAIN, { round, score: avg, byJudge: { investor: a, professor: b, executive: d }, hardGatePassed: passed }, 300),
    ...tool(MAIN, 'Write', { file_path: `07_评审/第${round}轮/评审汇总.json`, content: `{"round":${round},"score":${avg},…}` }, '文件已写入'),
    file(MAIN, `07_评审/第${round}轮/评审汇总.json`, 'write', 3_300),
    ...tool(MAIN, 'Write', { file_path: `07_评审/第${round}轮/修改清单.md`, content: '# 修改清单\n…' }, '文件已写入'),
    file(MAIN, `07_评审/第${round}轮/修改清单.md`, 'write', 2_800),
    ...text(
      MAIN,
      passed
        ? `第 ${round} 轮评审：均分 ${avg.toFixed(1)} ≥ 目标 ${c.p.params.targetScore}，硬门槛全部通过，提升 <2 分视为收敛，进入打包。`
        : `第 ${round} 轮评审：均分 ${avg.toFixed(1)}（目标 ${c.p.params.targetScore}），硬门槛 2 项未过（页数 84>80、禁语 3 处），按修改清单进入修改轮。`,
    ),
  ];
}

function fixRound(c: Ctx): Step[] {
  return [
    ...text(MAIN, `修改轮：只改清单中的问题——压缩附录至 4 页、补样本结构表、清除禁语、统一预算口径，然后重新渲染与检查。`),
    ...tool(MAIN, 'Edit', { file_path: '05_正文/A1_调研问卷.md', old_string: '（28 题完整版）', new_string: '（精简版，完整问卷见提交包附件）' }, '已编辑'),
    file(MAIN, '05_正文/A1_调研问卷.md', 'edit', 3_900),
    ...tool(MAIN, 'Edit', { file_path: '05_正文/02_现状分析.md', old_string: '## 2.3 消费者洞察', new_string: '## 2.3 消费者洞察\n\n表：问卷样本结构（n=412，示例数据）…' }, '已编辑'),
    file(MAIN, '05_正文/02_现状分析.md', 'edit', 17_400),
    ...tool(MAIN, 'Grep', { pattern: '众所周知|业内领先|绝对', path: '05_正文' }, '05_正文/01_引言.md:12\n05_正文/04_方案设计.md:88\n05_正文/05_推广与执行.md:140'),
    ...tool(MAIN, 'Edit', { file_path: '05_正文/01_引言.md', old_string: '众所周知，', new_string: '据《2025 年行业运行报告》，' }, '已编辑'),
    file(MAIN, '05_正文/01_引言.md', 'edit', 4_320),
    ...tool(MAIN, 'Bash', { command: 'python -m bizkit.docx_render 05_正文/book.json', description: '重新渲染' }, `✅ 产出/${c.outName}.docx（估算 63 页）\n✅ 匿名清洗 0 处`),
    file(MAIN, `产出/${c.outName}.docx`, 'write', 2_960_800),
  ];
}

function reviewStage(c: Ctx): Step[] {
  return [
    ...reviewRound(c, 1, [79.6, 82.8, 81.8], false),
    usage(6.0),
    ...fixRound(c),
    ...reviewRound(c, 2, [85.4, 87.2, 86.0], true),
    ...progress('review', 'done', '第 2 轮 86.2 分，收敛', 2),
    usage(7.4),
  ];
}

function packageStage(c: Ctx): Step[] {
  return [
    ...progress('package', 'running'),
    ...tool(MAIN, 'Bash', { command: `python -m bizkit.package --track ${c.trackId}`, description: '打包提交' }, `命名规则：{抽签编号或学校}-{团队名}-{作品名}\n✅ 08_提交包/XX-${c.p.teamName}-${c.p.name}.docx\n✅ 08_提交包/作品原创性声明.docx\n✅ 08_提交包/作品版权声明.docx\n✅ 08_提交包/企业授权书模板.docx\n✅ 08_提交包/数据真实性清单.md（6 条示例数据）\n✅ 08_提交包/README.md\n✅ 08_提交包/提交包.zip（9.4 MB）`),
    file(MAIN, '08_提交包/数据真实性清单.md', 'write', 3_100),
    file(MAIN, '08_提交包/README.md', 'write', 1_400),
    file(MAIN, '08_提交包/提交包.zip', 'write', 9_410_200),
    ...text(MAIN, `全部阶段完成。提交包已生成：${c.docName}（docx）、财务模型、路演 PPT、讲稿、题库、声明模板、数据真实性清单（6 条示例数据须替换为真实调研结果）。\n\n${checklist(9, null, '2/' + c.p.params.reviewRounds)}`),
    ...progress('package', 'done', '提交包.zip 9.4 MB'),
    usage(7.6),
    ev('done', MAIN, { subtype: 'success', result: '提交包已生成', usage: { turns: 137, costEstimateUsd: 1.98 } }, 400),
    ev('status', MAIN, { status: 'completed', message: '运行完成' }, 200),
  ];
}

const STAGE_BUILDERS: Record<StageId, (c: Ctx) => Step[]> = {
  parse: parseStage,
  research: researchStage,
  strategy: strategyStage,
  finance: financeStage,
  charts: chartsStage,
  write: writeStage,
  pitch: pitchStage,
  review: reviewStage,
  package: packageStage,
};

/** 单个阶段的脚本（用于为已完成阶段补历史事件） */
export function buildStageScript(project: Project, stage: StageId): Step[] {
  return STAGE_BUILDERS[stage](ctxOf(project));
}

export function buildRunScript(project: Project, fromStage: StageId, mode: 'run' | 'resume'): Step[] {
  const c = ctxOf(project);
  const startIdx = STAGE_IDS.indexOf(fromStage);
  const steps: Step[] = [ev('status', MAIN, { status: 'running', message: mode === 'resume' ? `从「${STAGES[startIdx].name}」阶段续跑` : '开始运行' }, 200)];
  if (mode === 'resume') {
    const doneNames = STAGES.slice(0, startIdx).map((s) => s.name).join('、');
    steps.push(...text(MAIN, `续跑：已完成阶段 ${doneNames || '无'}，从「${STAGES[startIdx].name}」继续；先核对已有产出物再开始。`));
    steps.push(...tool(MAIN, 'Glob', { pattern: '**/*.md' }, '00_任务解析/任务解析.md\n01_调研/01_市场规模与趋势.md\n…'));
  } else {
    steps.push(...text(MAIN, `收到任务书：${project.name}（${c.trackName} · ${project.group}）。按 CLAUDE.md 执行 9 个阶段，直至提交包生成。`));
  }
  for (let i = startIdx; i < STAGE_IDS.length; i++) steps.push(...STAGE_BUILDERS[STAGE_IDS[i]](c));
  return steps;
}

/** 修改模式脚本（修改对话） */
export function buildModifyScript(project: Project, message: string): Step[] {
  const c = ctxOf(project);
  const target = /竞品|第三章|3 章/.test(message) ? '05_正文/03_定位与核心价值.md' : /预算|财务/.test(message) ? '05_正文/07_财务预算.md' : /讲稿|路演|PPT|ppt/.test(message) ? '06_路演/讲稿.md' : '05_正文/02_现状分析.md';
  return [
    ev('status', MAIN, { status: 'running', message: '修改模式' }, 200),
    ...text(MAIN, `收到修改指令：“${message}”。计划：定位到 ${target}，按指令修改，保持数字与账本一致，然后重新渲染并自检。`),
    ...tool(MAIN, 'Read', { file_path: target }, `# ${target.split('/').pop()}\n…（共 240 行）`),
    ...thinking(MAIN, `修改需要新增内容但不能突破页数上限（当前 63/80）；新增表格约 0.5 页，可接受。所有新增数字必须带来源。`),
    ...tool(MAIN, 'Edit', { file_path: target, old_string: '（原段落）', new_string: '（修改后的段落，含新增表格与来源）' }, '已编辑'),
    file(MAIN, target, 'edit', 16_900),
    ...tool(MAIN, 'Bash', { command: `python -m bizkit.checks --chapter ${target} --track ${c.trackId}`, description: '章节自检' }, '数值密度 17.1/千字 ✓  来源配对 96% ✓  禁语 0 ✓'),
    ...tool(MAIN, 'Bash', { command: target.startsWith('06_') ? 'python -m bizkit.pptx_render 06_路演/slides.json' : 'python -m bizkit.docx_render 05_正文/book.json', description: '重新渲染' }, target.startsWith('06_') ? '✅ 产出/路演.pptx（12 页）' : `✅ 产出/${c.outName}.docx（估算 64 页）`),
    file(MAIN, target.startsWith('06_') ? '产出/路演.pptx' : `产出/${c.outName}.docx`, 'write', 2_980_000),
    ...text(MAIN, `已完成修改：${target} 已更新并重新渲染，章节自检通过（数值密度 17.1/千字，来源配对 96%）。相关数字与数字账本一致，跨材料一致性未受影响。`),
    ev('done', MAIN, { subtype: 'success', result: '修改完成' }, 300),
    ev('status', MAIN, { status: project.status === 'running' ? 'running' : 'completed', message: '修改完成' }, 200),
  ];
}
