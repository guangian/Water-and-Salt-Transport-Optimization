import type { Project, StageId, StageState, StageStatus } from '../types';

export interface StageDef {
  id: StageId;
  index: number;
  name: string;
  short: string;
  description: string;
  executor: string;
}

export const STAGES: StageDef[] = [
  { id: 'parse', index: 1, name: '读题解析', short: '读题', description: '命题解读、评分映射、材料清单', executor: '主代理' },
  { id: 'research', index: 2, name: '调研', short: '调研', description: '行业研究员并行检索，带来源数据', executor: '行业研究员 × N' },
  { id: 'strategy', index: 3, name: '策略', short: '策略', description: '定位逻辑链、商业模式、Big Idea', executor: '主代理' },
  { id: 'finance', index: 4, name: '财务', short: '财务', description: 'assumptions → 模型 xlsx + numbers.json', executor: '财务建模师' },
  { id: 'charts', index: 5, name: '图表', short: '图表', description: '300dpi 图表与 charts.json', executor: '工具包' },
  { id: 'write', index: 6, name: '撰写', short: '撰写', description: '各章 Markdown → docx/pdf', executor: '主笔' },
  { id: 'pitch', index: 7, name: '路演', short: '路演', description: 'slides.json → pptx、讲稿、题库', executor: '路演设计师' },
  { id: 'review', index: 8, name: '评审', short: '评审', description: '三评委打分 + 机器检查', executor: '评委团 + 质检员' },
  { id: 'package', index: 9, name: '打包', short: '打包', description: '命名、声明、数据真实性清单、zip', executor: '工具包' },
];

export const STAGE_IDS: StageId[] = STAGES.map((s) => s.id);

export const STAGE_BY_ID: Record<StageId, StageDef> = STAGES.reduce(
  (acc, s) => {
    acc[s.id] = s;
    return acc;
  },
  {} as Record<StageId, StageDef>,
);

export function stageName(id: StageId | null | undefined): string {
  return id ? STAGE_BY_ID[id]?.name ?? id : '—';
}

export function emptyStages(): Record<StageId, StageState> {
  return STAGE_IDS.reduce(
    (acc, id) => {
      acc[id] = { status: 'pending' };
      return acc;
    },
    {} as Record<StageId, StageState>,
  );
}

/** 阶段用时（毫秒）。运行中则按 now 计算。 */
export function stageDurationMs(stage: StageState | undefined, now: number = Date.now()): number | null {
  if (!stage?.startedAt) return null;
  const start = Date.parse(stage.startedAt);
  if (Number.isNaN(start)) return null;
  const end = stage.endedAt ? Date.parse(stage.endedAt) : stage.status === 'running' ? now : NaN;
  if (Number.isNaN(end)) return null;
  return Math.max(0, end - start);
}

/** 当前阶段（优先取 running；否则取最后一个 done 的下一个）。 */
export function currentStageId(stages: Record<StageId, StageState>): StageId | null {
  const running = STAGE_IDS.find((id) => stages[id]?.status === 'running');
  if (running) return running;
  const failed = STAGE_IDS.find((id) => stages[id]?.status === 'failed');
  if (failed) return failed;
  let last: StageId | null = null;
  for (const id of STAGE_IDS) {
    if (stages[id]?.status === 'done' || stages[id]?.status === 'skipped') last = id;
  }
  if (!last) return null;
  const idx = STAGE_IDS.indexOf(last);
  return idx + 1 < STAGE_IDS.length ? STAGE_IDS[idx + 1] : null;
}

/** 整体进度百分比：done=1、running=0.5。 */
export function stageProgressPercent(stages: Record<StageId, StageState>): number {
  let score = 0;
  for (const id of STAGE_IDS) {
    const st = stages[id]?.status;
    if (st === 'done' || st === 'skipped') score += 1;
    else if (st === 'running') score += 0.5;
  }
  return Math.round((score / STAGE_IDS.length) * 100);
}

export function countStages(stages: Record<StageId, StageState>): Record<StageStatus, number> {
  const out: Record<StageStatus, number> = { pending: 0, running: 0, done: 0, failed: 0, skipped: 0 };
  for (const id of STAGE_IDS) out[stages[id]?.status ?? 'pending'] += 1;
  return out;
}

const CHECKLIST_RE = /^\s*(\d)\.\s*(✅|🔄|⏳|❌)\s*(读题解析|调研|策略|财务|图表|撰写|路演|评审|打包)/u;
const EMOJI_STATUS: Record<string, StageStatus> = { '✅': 'done', '🔄': 'running', '⏳': 'pending', '❌': 'failed' };
const NAME_TO_ID: Record<string, StageId> = STAGES.reduce(
  (acc, s) => {
    acc[s.name] = s.id;
    return acc;
  },
  {} as Record<string, StageId>,
);

/** 兜底解析助手文本中的「✅/🔄/⏳ 清单」（契约 §6）。 */
export function parseChecklist(text: string): Partial<Record<StageId, StageStatus>> {
  const out: Partial<Record<StageId, StageStatus>> = {};
  for (const line of text.split(/\r?\n/)) {
    const m = CHECKLIST_RE.exec(line);
    if (!m) continue;
    const id = NAME_TO_ID[m[3]];
    const status = EMOJI_STATUS[m[2]];
    if (id && status) out[id] = status;
  }
  return out;
}

export function projectProgress(project: Pick<Project, 'stages'>): number {
  return stageProgressPercent(project.stages);
}
