import path from 'node:path';
import fsp from 'node:fs/promises';
import crypto from 'node:crypto';
import { EventEmitter } from 'node:events';
import type {
  AgentEvent,
  AgentEventInput,
  ModelMap,
  Project,
  ProjectCreateInput,
  ProjectSummary,
  StageId,
  StageState,
  TrackId,
} from '../types.js';
import { ApiError, STAGE_IDS } from '../types.js';
import { ensureDir, exists, readJsonOrNull, writeJsonAtomic } from '../util/fs.js';
import { createLogger, redactDeep } from '../util/log.js';

const log = createLogger('store');

export const PROJECT_ID_RE = /^p_\d{8}_\d{6}_[a-z0-9]{4}$/;
const TRACKS: TrackId[] = ['brand', 'startup', 'digital', 'case', 'mice'];

function pad(n: number, w = 2): string {
  return String(n).padStart(w, '0');
}

export function newProjectId(now = new Date()): string {
  const d = `${now.getFullYear()}${pad(now.getMonth() + 1)}${pad(now.getDate())}`;
  const t = `${pad(now.getHours())}${pad(now.getMinutes())}${pad(now.getSeconds())}`;
  const rand = crypto.randomBytes(3).toString('base64url').toLowerCase().replace(/[^a-z0-9]/g, 'x').slice(0, 4).padEnd(4, 'x');
  return `p_${d}_${t}_${rand}`;
}

export function emptyStages(): Record<StageId, StageState> {
  const stages = {} as Record<StageId, StageState>;
  for (const id of STAGE_IDS) stages[id] = id === 'review' ? { status: 'pending', round: 0 } : { status: 'pending' };
  return stages;
}

function str(v: unknown, fallback = ''): string {
  return typeof v === 'string' ? v : fallback;
}

function num(v: unknown, fallback: number): number {
  return typeof v === 'number' && Number.isFinite(v) ? v : fallback;
}

/** 把客户端输入整理成完整 project.json（契约 §2），缺省值补齐，未知字段忽略 */
export function buildProject(input: ProjectCreateInput, defaults: { reviewRounds: number; targetScore: number; models: Partial<ModelMap> }, now = new Date()): Project {
  const track = str(input.track) as TrackId;
  if (!TRACKS.includes(track)) throw new ApiError(400, 'bad_request', `未知赛道：${String(input.track)}`);
  const name = str(input.name).trim();
  if (!name) throw new ApiError(400, 'bad_request', '作品名不能为空');
  const subj = (input.subject ?? {}) as Partial<Project['subject']>;
  const params = (input.params ?? {}) as Partial<Project['params']>;
  const iso = now.toISOString();
  return {
    id: newProjectId(now),
    name,
    teamName: str(input.teamName).trim(),
    track,
    group: str(input.group, '本科组'),
    subject: {
      type: (['brand', 'startup', 'case-company', 'event'] as const).includes(subj.type as never) ? (subj.type as Project['subject']['type']) : track === 'startup' ? 'startup' : track === 'case' ? 'case-company' : track === 'mice' ? 'event' : 'brand',
      name: str(subj.name),
      company: str(subj.company),
      isAssignedCase: Boolean(subj.isAssignedCase),
      caseTitle: str(subj.caseTitle),
      brief: str(subj.brief),
      region: str(subj.region),
      industry: str(subj.industry),
    },
    team: Array.isArray(input.team)
      ? input.team.map((m) => ({ role: str(m?.role), major: str(m?.major), strength: str(m?.strength) }))
      : [],
    params: {
      researchMode: params.researchMode === 'provided' ? 'provided' : 'synthetic',
      pageTarget: Math.min(80, Math.max(5, num(params.pageTarget, 60))),
      pitchMinutes: Math.max(3, num(params.pitchMinutes, 10)),
      reviewRounds: Math.max(0, num(params.reviewRounds, defaults.reviewRounds)),
      targetScore: Math.min(100, Math.max(0, num(params.targetScore, defaults.targetScore))),
      language: str(params.language, 'zh-CN'),
      models: { ...defaults.models, ...(params.models ?? {}) },
    },
    inputs: [],
    status: 'idle',
    currentStage: '',
    stages: emptyStages(),
    session: { sessionId: '', lastPromptAt: '', mode: '' },
    usage: { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, turns: 0, costEstimateUsd: 0, durationMs: 0 },
    review: { rounds: [], latestScore: null },
    artifacts: {},
    createdAt: iso,
    updatedAt: iso,
  };
}

export function summarize(p: Project): ProjectSummary {
  const stages = {} as Record<StageId, StageState['status']>;
  for (const id of STAGE_IDS) stages[id] = p.stages[id]?.status ?? 'pending';
  return {
    id: p.id,
    name: p.name,
    teamName: p.teamName,
    track: p.track,
    group: p.group,
    status: p.status,
    currentStage: p.currentStage,
    latestScore: p.review?.latestScore ?? null,
    stages,
    createdAt: p.createdAt,
    updatedAt: p.updatedAt,
  };
}

/**
 * 项目存储：workspace/<id>/project.json + events.jsonl。
 * 写 project.json 原子；events.jsonl 追加写，seq 单调递增（进程内缓存，首次从文件尾部恢复）。
 * 同一项目的事件写入通过 per-project 队列串行化，避免并发追加交错。
 */
export class ProjectStore {
  private seqCache = new Map<string, number>();
  private writeQueues = new Map<string, Promise<void>>();
  private bus = new EventEmitter();

  constructor(readonly workspaceDir: string) {
    this.bus.setMaxListeners(0);
  }

  dirOf(id: string): string {
    if (!PROJECT_ID_RE.test(id)) throw new ApiError(400, 'bad_id', '项目 id 格式不正确');
    return path.join(this.workspaceDir, id);
  }

  private fileOf(id: string): string {
    return path.join(this.dirOf(id), 'project.json');
  }

  private eventsFileOf(id: string): string {
    return path.join(this.dirOf(id), 'events.jsonl');
  }

  async list(): Promise<ProjectSummary[]> {
    await ensureDir(this.workspaceDir);
    const names = await fsp.readdir(this.workspaceDir).catch(() => [] as string[]);
    const out: ProjectSummary[] = [];
    for (const n of names) {
      if (!PROJECT_ID_RE.test(n)) continue;
      const p = await readJsonOrNull<Project>(path.join(this.workspaceDir, n, 'project.json'));
      if (p) out.push(summarize(p));
    }
    out.sort((a, b) => b.updatedAt.localeCompare(a.updatedAt));
    return out;
  }

  async get(id: string): Promise<Project> {
    const p = await readJsonOrNull<Project>(this.fileOf(id));
    if (!p) throw new ApiError(404, 'not_found', `项目不存在：${id}`);
    if (!p.stages) p.stages = emptyStages();
    return p;
  }

  async exists(id: string): Promise<boolean> {
    return PROJECT_ID_RE.test(id) && (await exists(this.fileOf(id)));
  }

  async save(p: Project): Promise<Project> {
    p.updatedAt = new Date().toISOString();
    await writeJsonAtomic(this.fileOf(p.id), p);
    return p;
  }

  /** 读-改-写（串行化，避免并发覆盖） */
  async mutate(id: string, fn: (p: Project) => void | Promise<void>): Promise<Project> {
    return this.enqueue(`project:${id}`, async () => {
      const p = await this.get(id);
      await fn(p);
      return this.save(p);
    });
  }

  async create(p: Project): Promise<Project> {
    const dir = this.dirOf(p.id);
    await ensureDir(path.join(dir, 'inputs'));
    await writeJsonAtomic(this.fileOf(p.id), p);
    await fsp.writeFile(this.eventsFileOf(p.id), '', { flag: 'a' });
    log.info('project created', { id: p.id, track: p.track });
    return p;
  }

  async remove(id: string): Promise<void> {
    const dir = this.dirOf(id);
    await fsp.rm(dir, { recursive: true, force: true, maxRetries: 3, retryDelay: 200 });
    this.seqCache.delete(id);
    log.info('project removed', { id });
  }

  // ---------- 事件 ----------

  private async lastSeq(id: string): Promise<number> {
    const cached = this.seqCache.get(id);
    if (cached !== undefined) return cached;
    let seq = 0;
    try {
      const text = await fsp.readFile(this.eventsFileOf(id), 'utf8');
      const idx = text.lastIndexOf('\n', text.length - 2);
      const last = text.slice(idx + 1).trim();
      if (last) seq = (JSON.parse(last) as AgentEvent).seq ?? 0;
    } catch {
      seq = 0;
    }
    this.seqCache.set(id, seq);
    return seq;
  }

  /**
   * 追加事件：分配 seq/ts，脱敏后落盘并广播。
   * persist=false 时只分配 seq 并广播（用于 text_delta / thinking 增量）：
   * 增量随后总会以完整 text/thinking 块再次到达，落盘只会让长跑的回放文件膨胀几十倍。
   */
  async appendEvent(id: string, input: AgentEventInput, opts: { persist?: boolean } = {}): Promise<AgentEvent> {
    return this.enqueue(`events:${id}`, async () => {
      const seq = (await this.lastSeq(id)) + 1;
      const ev: AgentEvent = redactDeep({ seq, ts: new Date().toISOString(), projectId: id, ...input });
      if (opts.persist !== false) {
        await fsp.appendFile(this.eventsFileOf(id), JSON.stringify(ev) + '\n', 'utf8');
      }
      this.seqCache.set(id, seq);
      this.bus.emit(id, ev);
      return ev;
    });
  }

  async readEvents(id: string, since = 0): Promise<AgentEvent[]> {
    let text = '';
    try {
      text = await fsp.readFile(this.eventsFileOf(id), 'utf8');
    } catch {
      return [];
    }
    const out: AgentEvent[] = [];
    for (const line of text.split('\n')) {
      if (!line.trim()) continue;
      try {
        const ev = JSON.parse(line) as AgentEvent;
        if (ev.seq > since) out.push(ev);
      } catch {
        /* 跳过损坏行 */
      }
    }
    return out;
  }

  subscribe(id: string, listener: (ev: AgentEvent) => void): () => void {
    this.bus.on(id, listener);
    return () => this.bus.off(id, listener);
  }

  private enqueue<T>(key: string, task: () => Promise<T>): Promise<T> {
    const prev = this.writeQueues.get(key) ?? Promise.resolve();
    const run = prev.then(task, task);
    const settled = run.then(
      () => undefined,
      () => undefined,
    );
    this.writeQueues.set(key, settled);
    settled.then(() => {
      if (this.writeQueues.get(key) === settled) this.writeQueues.delete(key);
    });
    return run;
  }
}
