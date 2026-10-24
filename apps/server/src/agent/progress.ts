import type { StageId, StageStatus } from '../types.js';

/**
 * 进度兜底解析（契约 §6）：主代理在对话中输出的
 *   1. ✅ 读题解析
 *   8. ⏳ 评审（第 0/3 轮）
 * 清单。update_progress 工具是主通道，这里只做兜底，所以解析结果由调用方按“只升不降”合并。
 */

const NAME_TO_STAGE: Record<string, StageId> = {
  读题解析: 'parse',
  调研: 'research',
  策略: 'strategy',
  财务: 'finance',
  图表: 'charts',
  撰写: 'write',
  路演: 'pitch',
  评审: 'review',
  打包: 'package',
};

const EMOJI_TO_STATUS: Record<string, StageStatus> = {
  '✅': 'done',
  '🔄': 'running',
  '⏳': 'pending',
  '❌': 'failed',
};

export const PROGRESS_LINE_RE = /^\s*(\d)\.\s*(✅|🔄|⏳|❌)\s*(读题解析|调研|策略|财务|图表|撰写|路演|评审|打包)(?:\s*[（(]\s*第\s*(\d+)\s*\/\s*(\d+)\s*轮\s*[)）])?/u;

export interface ProgressItem {
  stage: StageId;
  status: StageStatus;
  round?: number;
  totalRounds?: number;
}

export function parseProgressChecklist(text: string): ProgressItem[] {
  if (!text || !/[✅🔄⏳❌]/u.test(text)) return [];
  const out: ProgressItem[] = [];
  const seen = new Set<StageId>();
  for (const line of text.split(/\r?\n/)) {
    const m = PROGRESS_LINE_RE.exec(line);
    if (!m) continue;
    const stage = NAME_TO_STAGE[m[3] ?? ''];
    const status = EMOJI_TO_STATUS[m[2] ?? ''];
    if (!stage || !status || seen.has(stage)) continue;
    seen.add(stage);
    const item: ProgressItem = { stage, status };
    if (m[4] !== undefined) item.round = Number(m[4]);
    if (m[5] !== undefined) item.totalRounds = Number(m[5]);
    out.push(item);
  }
  return out;
}

const RANK: Record<StageStatus, number> = { pending: 0, skipped: 1, running: 2, failed: 3, done: 4 };

/** 兜底解析只允许“推进”阶段状态：done 不会被 pending/running 覆盖，running 不会退回 pending */
export function shouldApplyFallback(current: StageStatus | undefined, next: StageStatus): boolean {
  if (!current) return true;
  if (current === next) return false;
  return RANK[next] > RANK[current];
}
