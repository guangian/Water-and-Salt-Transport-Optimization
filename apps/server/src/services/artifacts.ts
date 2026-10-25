import path from 'node:path';
import fsp from 'node:fs/promises';
import type { Project, ReviewRoundSummary } from '../types.js';
import { readJsonOrNull, toPosixRelative } from '../util/fs.js';
import type { NumbersLedger } from '../agent/mcp.js';

/** 产出物路径约定（架构 §2.2）→ project.artifacts 的键 */
export async function scanArtifacts(ws: string): Promise<Record<string, string>> {
  const out: Record<string, string> = {};
  const pick = async (dir: string, pred: (name: string) => boolean, key: string) => {
    const abs = path.join(ws, dir);
    const names = await fsp.readdir(abs).catch(() => [] as string[]);
    const hits = names.filter(pred).filter((n) => !n.startsWith('~$') && !n.startsWith('.'));
    if (hits.length === 0) return;
    // 多个候选取最近修改的
    let best: { name: string; mtime: number } | null = null;
    for (const n of hits) {
      const st = await fsp.stat(path.join(abs, n)).catch(() => null);
      if (st && (!best || st.mtimeMs > best.mtime)) best = { name: n, mtime: st.mtimeMs };
    }
    if (best) out[key] = toPosixRelative(ws, path.join(abs, best.name));
  };
  await pick('产出', (n) => n.endsWith('.docx'), 'book');
  await pick('产出', (n) => n.endsWith('.pdf'), 'pdf');
  await pick('产出', (n) => n.endsWith('.pptx'), 'deck');
  await pick('03_财务', (n) => n.endsWith('.xlsx'), 'model');
  await pick('06_路演', (n) => n === '讲稿.md', 'script');
  await pick('06_路演', (n) => n === '答辩题库.md', 'qa');
  await pick('08_提交包', (n) => n.endsWith('.zip'), 'package');
  return out;
}

export interface ReviewSummaryFile {
  round?: number;
  score?: number;
  byJudge?: Record<string, number>;
  hardGatePassed?: boolean;
  machine?: unknown;
  fixes?: unknown[];
  converged?: boolean;
}

export interface MachineCheckReport {
  track?: string;
  generatedAt?: string;
  checks?: unknown[];
  hardGates?: unknown[];
  crossDoc?: unknown;
  summary?: { pass: number; warn: number; fail: number };
}

export type ReviewRoundData = ReviewSummaryFile & { path: string; judges?: Record<string, unknown>[]; machineReport?: MachineCheckReport };

export interface ReviewData {
  latest: ReviewRoundData | null;
  rounds: ReviewRoundData[];
  machine: MachineCheckReport | null;
  judges: Record<string, unknown>[];
}

/** 读取 07_评审/第N轮/{评审汇总.json, 机器检查.json, 评审_*.json}；缺则为空 */
export async function readReview(ws: string): Promise<ReviewData> {
  const dir = path.join(ws, '07_评审');
  const names = await fsp.readdir(dir).catch(() => [] as string[]);
  const roundDirs = names
    .map((n) => ({ n, m: /^第(\d+)轮$/.exec(n) }))
    .filter((x): x is { n: string; m: RegExpExecArray } => !!x.m)
    .map((x) => ({ name: x.n, round: Number(x.m[1]) }))
    .sort((a, b) => a.round - b.round);
  // 每轮目录里的 评审_<judge>.json（三评委原始打分），前端评审页按 round.judges 渲染
  const readJudges = async (roundDir: string): Promise<Record<string, unknown>[]> => {
    const out: Record<string, unknown>[] = [];
    const files = await fsp.readdir(path.join(dir, roundDir)).catch(() => [] as string[]);
    for (const f of files.sort()) {
      if (/^评审_.+\.json$/.test(f) && f !== '评审汇总.json') {
        const j = await readJsonOrNull<Record<string, unknown>>(path.join(dir, roundDir, f));
        if (j) out.push({ judge: f.slice(3, -5), ...j });
      }
    }
    return out;
  };
  const rounds: ReviewData['rounds'] = [];
  for (const r of roundDirs) {
    const sum = await readJsonOrNull<ReviewSummaryFile>(path.join(dir, r.name, '评审汇总.json'));
    if (!sum) continue;
    const machine = await readJsonOrNull<MachineCheckReport>(path.join(dir, r.name, '机器检查.json'));
    const judges = await readJudges(r.name);
    rounds.push({ ...sum, round: sum.round ?? r.round, path: `07_评审/${r.name}`, judges, ...(machine ? { machineReport: machine } : {}) });
  }
  const latestDir = roundDirs.at(-1);
  const machine = latestDir ? await readJsonOrNull<MachineCheckReport>(path.join(dir, latestDir.name, '机器检查.json')) : null;
  const judges = latestDir ? await readJudges(latestDir.name) : [];
  return { latest: rounds.at(-1) ?? null, rounds, machine, judges };
}

export function toReviewRounds(data: ReviewData): ReviewRoundSummary[] {
  return data.rounds.map((r) => ({
    round: r.round ?? 0,
    score: typeof r.score === 'number' ? r.score : 0,
    byJudge: r.byJudge ?? {},
    hardGatePassed: Boolean(r.hardGatePassed),
    path: r.path,
  }));
}

/** 同步评审结果到 project（返回是否有变化） */
export async function syncReviewIntoProject(ws: string, p: Project): Promise<{ changed: boolean; latest: ReviewRoundSummary | null }> {
  const data = await readReview(ws);
  const rounds = toReviewRounds(data);
  const latest = rounds.at(-1) ?? null;
  const before = JSON.stringify(p.review);
  p.review = { rounds, latestScore: latest ? latest.score : null };
  return { changed: before !== JSON.stringify(p.review), latest };
}

export interface NumbersReport {
  ledger: NumbersLedger | null;
  consistency: { checked: number; consistent: boolean; issues: { key: string; message: string }[] };
}

/**
 * numbers.json + 一致性检查：同 key 重复、缺来源、非数值、占位符示例数据未标注。
 * 跨材料比对（正文/PPT/讲稿）由工具包 checks 完成，这里只做账本自洽性。
 */
export async function readNumbers(ws: string): Promise<NumbersReport> {
  const ledger = await readJsonOrNull<NumbersLedger>(path.join(ws, '03_财务', 'numbers.json'));
  const issues: NumbersReport['consistency']['issues'] = [];
  const entries = Array.isArray(ledger?.entries) ? ledger!.entries : [];
  const seen = new Map<string, number>();
  for (const e of entries) {
    if (!e || typeof e.key !== 'string') {
      issues.push({ key: '?', message: '条目缺少 key' });
      continue;
    }
    seen.set(e.key, (seen.get(e.key) ?? 0) + 1);
    if (typeof e.value !== 'number' && (typeof e.value !== 'string' || !/^-?[\d.,]+%?$/.test(e.value))) issues.push({ key: e.key, message: `value 非数值：${String(e.value)}` });
    if (!e.source && !e.synthetic) issues.push({ key: e.key, message: '缺少来源（source）' });
    if (typeof e.unit === 'string' && e.unit === '%' && typeof e.value === 'number' && e.value > 1.5 && e.value <= 100) issues.push({ key: e.key, message: `单位为 % 但值为 ${e.value}，请统一用小数（0.57）或百分数并在 label 说明` });
  }
  for (const [k, n] of seen) if (n > 1) issues.push({ key: k, message: `key 重复 ${n} 次` });
  return { ledger, consistency: { checked: entries.length, consistent: issues.length === 0, issues } };
}
