import type { Project, StageId, TrackInfo } from '../types.js';
import { STAGE_IDS, STAGE_NAMES } from '../types.js';

const SUBJECT_TYPE_NAMES: Record<Project['subject']['type'], string> = {
  brand: '品牌',
  startup: '创业项目',
  'case-company': '案例企业',
  event: '会奖活动',
};

function cjkWordCount(text: string): number {
  const cjk = (text.match(/[\u3400-\u4dbf\u4e00-\u9fff]/g) ?? []).length;
  const latin = (text.match(/[A-Za-z0-9][A-Za-z0-9.%,]*/g) ?? []).length;
  return cjk + latin;
}

export function formatWordCount(n: number): string {
  if (n >= 10000) return `${(n / 10000).toFixed(1)} 万字`;
  return `${n} 字`;
}

/** 契约 §5 任务书 */
export function renderTaskSheet(p: Project, track: TrackInfo | undefined, prereadWords: Record<string, number> = {}): string {
  const s = p.subject;
  const subjectBits = [s.company, s.region, s.industry].filter(Boolean).join('，');
  const subjectLine = `${SUBJECT_TYPE_NAMES[s.type] ?? '对象'} ${s.name || '（未填写）'}${subjectBits ? `（${subjectBits}）` : ''}`;
  const team = p.team.length ? p.team.map((m) => `${m.role || '成员'}（${[m.major, m.strength].filter(Boolean).join('，') || '—'}）`).join('；') : '（未填写，按 3–5 人匿名团队处理）';
  const research = p.params.researchMode === 'provided' ? '仅用用户资料' : '示例数据（须逐条标注【示例数据，须替换为真实调研结果】）';
  const inputs = p.inputs.length
    ? p.inputs.map((i) => `  - ${i.path}${i.preread ? ` → ${i.preread}${prereadWords[i.name] !== undefined ? `（${formatWordCount(prereadWords[i.name]!)}）` : ''}` : '（二进制/图片，需要时用 Read 查看）'}`).join('\n')
    : '  - （无上传资料）';
  const trackName = track ? `${track.name}（${p.track}）` : p.track;
  const limit = track?.pageLimit ?? 80;
  const lines = [
    '# 任务书',
    `- 赛道：${trackName} | 组别：${p.group || '本科组'} | 语言：${p.params.language === 'en' ? '英文' : '中文'}`,
    `- 作品名：${p.name} | 团队名：${p.teamName || '（未填写）'}（匿名制：全文不得出现院校/姓名/联系方式）`,
    `- 策划对象：${subjectLine}`,
    `- 是否命题企业：${s.isAssignedCase ? '是' : '否'} | 命题题目：${s.caseTitle || '—'}`,
    `- 项目简介：${s.brief ? s.brief.replace(/\r?\n/g, ' ') : '（未填写，请基于公开信息与调研补全）'}`,
    `- 团队（匿名）：${team}`,
    `- 参数：调研模式=${research} | 正文目标 ${p.params.pageTarget} 页（上限 ${limit}）| 讲稿 ${p.params.pitchMinutes} 分钟 | 评审轮数 ${p.params.reviewRounds} | 目标分 ${p.params.targetScore}`,
    '- 用户资料（已预读为文本，路径见下）：',
    inputs,
    `- 赛道规范：tracks/${p.track}/spec.md | 量规：tracks/${p.track}/rubric.json | 命名：tracks/${p.track}/naming.json | 检查阈值：tracks/${p.track}/checks.json`,
    track ? `- 赛道摘要：${track.docName}；评分 书面 ${track.scoring.written} + 展示答辩 ${track.scoring.pitch}；必备要素：${track.elements.join('、')}` : '',
    track?.notes ? `- 赛道要点：${track.notes}` : '',
    '- 开始指令：按 CLAUDE.md 执行流程完成全部阶段，直至 08_提交包 生成。',
    '',
  ].filter((l) => l !== '');
  return lines.join('\n') + '\n';
}

export function prereadWordCount(text: string): number {
  return cjkWordCount(text);
}

export interface PromptContext {
  mode: 'run' | 'resume' | 'modify' | 'smoke';
  project: Project;
  fromStage?: StageId;
  message?: string;
}

export const RUN_PROMPT = '请阅读 任务书.md 与 CLAUDE.md，按执行流程完成全部阶段，直至 08_提交包 生成。全程自动执行，不要向用户提问。';

export const SMOKE_PROMPT =
  '这是一次冒烟测试，请严格按顺序执行：' +
  '1) 调用 mcp__bizgold__update_progress 把 stage=parse 设为 running；' +
  '2) 用 mcp__bizgold__web_search 搜索「全国高校商业精英挑战赛 2026 品牌策划」；' +
  '3) 对第一条结果调用 mcp__bizgold__web_fetch；' +
  '4) 把搜索与抓取的要点（标题、URL、不超过 300 字摘要）写入 00_任务解析/冒烟.md；' +
  '5) 调用 update_progress 把 parse 设为 done；' +
  '6) 最后输出契约 §6 格式的 9 行进度清单（1. ✅ 读题解析 … 9. ⏳ 打包）。不要做其他事情。';

/** 续跑/修改提示词：带已完成阶段清单与继续点 */
export function buildPrompt(ctx: PromptContext): string {
  const { project: p } = ctx;
  if (ctx.mode === 'smoke') return SMOKE_PROMPT;
  if (ctx.mode === 'run') return RUN_PROMPT;
  const done = STAGE_IDS.filter((s) => p.stages[s]?.status === 'done').map((s) => STAGE_NAMES[s]);
  const failed = STAGE_IDS.filter((s) => p.stages[s]?.status === 'failed').map((s) => STAGE_NAMES[s]);
  if (ctx.mode === 'resume') {
    const from = ctx.fromStage ?? firstUnfinishedStage(p) ?? 'package';
    return [
      '这是一次续跑。请先阅读 任务书.md 与 CLAUDE.md，并检查工作区现有产出，不要重做已完成阶段。',
      `已完成阶段：${done.length ? done.join('、') : '无'}。${failed.length ? `上次失败阶段：${failed.join('、')}。` : ''}`,
      `请从「${STAGE_NAMES[from]}」阶段继续，按 CLAUDE.md 执行流程完成剩余阶段，直至 08_提交包 生成。`,
      '每个阶段开始/结束仍需调用 mcp__bizgold__update_progress 上报进度。',
      '注意：子代理（Agent 工具）是前台同步执行的——调用返回时结果就在工具结果里，不存在“等待通知”；不要为了等待而结束回复。若上次是在等待子代理时结束的，先 ls 检查该子代理的输出文件是否已生成：已生成就按文件回收，没生成就用原任务书重新派发一次。',
    ].join('\n');
  }
  // modify
  return [
    '用户对已生成的材料提出修改要求：',
    '"""',
    (ctx.message ?? '').trim(),
    '"""',
    '约束：不得修改 inputs/ 目录；只改与要求相关的文件；改完后重新运行相应的渲染（docx_render / pptx_render / charts / finmodel）与机器检查（checks），',
    '若影响数字账本（numbers.json）须同步更新正文、PPT 与讲稿中的引用；完成后用 update_progress 上报受影响阶段，并输出 ≤200 字的修改说明。',
  ].join('\n');
}

export function firstUnfinishedStage(p: Project): StageId | undefined {
  return STAGE_IDS.find((s) => p.stages[s]?.status !== 'done' && p.stages[s]?.status !== 'skipped');
}

/** systemPrompt.append：简短运行环境说明（领域规范都在 CLAUDE.md，这里只讲环境事实） */
export function buildSystemAppend(opts: { workspaceDir: string; pythonCmd: string; isWindows: boolean; projectId: string }): string {
  const py = opts.pythonCmd;
  return [
    '## 运行环境（金策 BizGold 服务端）',
    `- 当前工作区：${opts.workspaceDir}（即 cwd）。所有产出写在工作区内的相对路径下；不要在工作区外创建文件。`,
    `- Python 解释器：${py}（环境变量 BIZGOLD_PYTHON 同值）。工具包位于 toolkit/，调用方式：在工作区根执行 \`${quote(py)} -m bizkit.<模块> ...\`，PYTHONPATH 已指向 toolkit/，无需 cd。`,
    `- 自定义工具：mcp__bizgold__web_search（联网检索）、mcp__bizgold__web_fetch（抓取网页/PDF 正文）、mcp__bizgold__update_progress（上报阶段进度，每个阶段开始与结束都必须调用）、mcp__bizgold__record_numbers（登记关键数字到 03_财务/numbers.json）。内置 WebSearch 不可用。`,
    `- 本次运行无人值守：不要向用户提问或等待确认；遇到可选项自行决定并在产出中说明假设。${opts.isWindows ? ' 当前为 Windows：Bash 工具运行在 Git Bash 中，路径用正斜杠，命令里避免 PowerShell 专有语法。' : ''}`,
    '- 所有文本文件 UTF-8 编码；Python 脚本首行加 `# -*- coding: utf-8 -*-`；中文回答。',
  ].join('\n');
}

function quote(p: string): string {
  return /\s/.test(p) ? `"${p}"` : p;
}
