import type { BizgoldApi, EventHandlers } from '../api';
import { ApiError } from '../lib/errors';
import type {
  AgentEvent,
  ArtifactKey,
  InputKind,
  NumbersData,
  Project,
  ProjectInput,
  ReviewData,
  Settings,
  StageId,
} from '../types';
import { currentStageId, emptyStages, STAGE_IDS } from '../lib/stages';
import { loadJson, saveJson } from '../lib/storage';
import { MOCK_TRACKS } from './tracks';
import { MOCK_HEALTH, MOCK_MODELS, MOCK_SETTINGS } from './settings';
import { DEFAULT_MODELS, MOCK_PROJECTS, RUNNING_PROJECT_ID } from './projects';
import { artifactsFor, buildTree, mockFileUrl, previewFor, PROJECT_FILES } from './artifacts';
import { makeReviewRound, MOCK_REVIEW_BRAND } from './review';
import { MOCK_NUMBERS_BRAND, MOCK_NUMBERS_STARTUP } from './numbers';
import { knowledgeContent, MOCK_KNOWLEDGE_DOCS } from './knowledge';
import { buildModifyScript, buildRunScript, buildStageScript, type Step } from './script';

const LS_CREATED = 'bizgold.mock.createdProjects';
const LS_DELETED = 'bizgold.mock.deletedIds';
const LS_SETTINGS = 'bizgold.mock.settings';

function clone<T>(v: T): T {
  return JSON.parse(JSON.stringify(v)) as T;
}

function sleep(ms: number): Promise<void> {
  return new Promise((r) => setTimeout(r, ms));
}

function latency(): Promise<void> {
  return sleep(120 + Math.random() * 220);
}

function kindOf(name: string): InputKind {
  const ext = name.split('.').pop()?.toLowerCase() ?? '';
  if (ext === 'pdf') return 'pdf';
  if (ext === 'docx' || ext === 'doc') return 'docx';
  if (ext === 'xlsx' || ext === 'xls') return 'xlsx';
  if (ext === 'csv') return 'csv';
  if (['png', 'jpg', 'jpeg', 'webp', 'gif'].includes(ext)) return 'image';
  if (['txt', 'md'].includes(ext)) return 'text';
  return 'other';
}

const ARTIFACT_RULES: [RegExp, ArtifactKey][] = [
  [/^产出\/.*\.docx$/, 'book'],
  [/^产出\/.*\.pdf$/, 'pdf'],
  [/^产出\/路演\.pptx$/, 'deck'],
  [/^03_财务\/财务模型\.xlsx$/, 'model'],
  [/^06_路演\/讲稿\.md$/, 'script'],
  [/^06_路演\/答辩题库\.md$/, 'qa'],
  [/^08_提交包\/提交包\.zip$/, 'package'],
];

interface FileSpec {
  path: string;
  size: number;
  mtime?: string;
}

class MockStore {
  projects = new Map<string, Project>();
  reviews = new Map<string, ReviewData>();
  numbers = new Map<string, NumbersData>();
  extraFiles = new Map<string, FileSpec[]>();
  settings: Settings;
  runStartedAt = new Map<string, number>();

  constructor() {
    const deleted = new Set(loadJson<string[]>(LS_DELETED, []));
    for (const p of MOCK_PROJECTS) if (!deleted.has(p.id)) this.projects.set(p.id, clone(p));
    for (const p of loadJson<Project[]>(LS_CREATED, [])) {
      const copy = clone(p);
      // 页面刷新后 Mock 运行器不复存在：运行中的项目回退为“已暂停”，可续跑
      if (copy.status === 'running') {
        copy.status = 'paused';
        for (const id of STAGE_IDS) if (copy.stages[id].status === 'running') copy.stages[id] = { ...copy.stages[id], status: 'pending', startedAt: undefined };
      }
      this.projects.set(copy.id, copy);
    }
    this.reviews.set('p_20261018_091533_a1f3', clone(MOCK_REVIEW_BRAND));
    this.numbers.set('p_20261018_091533_a1f3', clone(MOCK_NUMBERS_BRAND));
    this.numbers.set(RUNNING_PROJECT_ID, clone(MOCK_NUMBERS_STARTUP));
    this.settings = loadJson<Settings>(LS_SETTINGS, clone(MOCK_SETTINGS));
  }

  persist(): void {
    const created = [...this.projects.values()].filter((p) => !MOCK_PROJECTS.some((s) => s.id === p.id));
    saveJson(LS_CREATED, created);
    const deleted = MOCK_PROJECTS.filter((s) => !this.projects.has(s.id)).map((s) => s.id);
    saveJson(LS_DELETED, deleted);
  }

  get(id: string): Project {
    const p = this.projects.get(id);
    if (!p) throw new ApiError(`项目不存在：${id}`, 'not_found', 404);
    return p;
  }

  files(id: string): FileSpec[] {
    const base = PROJECT_FILES[id] ?? [{ path: '任务书.md', size: 2_000 }];
    const extra = this.extraFiles.get(id) ?? [];
    const map = new Map<string, FileSpec>();
    for (const f of [...base, ...extra]) map.set(f.path, f);
    return [...map.values()];
  }

  /** 事件副作用：同步 project.json 的阶段/状态/用量/评分/产出物 */
  apply(id: string, ev: AgentEvent): void {
    const p = this.projects.get(id);
    if (!p) return;
    switch (ev.type) {
      case 'stage': {
        const prev = p.stages[ev.data.stage] ?? { status: 'pending' };
        const next = { ...prev, status: ev.data.status, note: ev.data.note ?? prev.note, round: ev.data.round ?? prev.round };
        if (ev.data.status === 'running') {
          next.startedAt = prev.status === 'running' && prev.startedAt ? prev.startedAt : ev.ts;
          next.endedAt = undefined;
        }
        if (ev.data.status === 'done' || ev.data.status === 'failed' || ev.data.status === 'skipped') next.endedAt = ev.ts;
        p.stages[ev.data.stage] = next;
        p.currentStage = ev.data.status === 'running' ? ev.data.stage : currentStageId(p.stages);
        break;
      }
      case 'status': {
        p.status = ev.data.status;
        if (ev.data.status === 'completed') p.currentStage = null;
        break;
      }
      case 'usage': {
        const started = this.runStartedAt.get(id);
        p.usage = { ...ev.data, durationMs: started ? Date.now() - started + (p.usage.durationMs ?? 0) * 0 : p.usage.durationMs };
        break;
      }
      case 'review': {
        const summary = { round: ev.data.round, score: ev.data.score, byJudge: ev.data.byJudge, hardGatePassed: ev.data.hardGatePassed, path: `07_评审/第${ev.data.round}轮` };
        const idx = p.review.rounds.findIndex((r) => r.round === ev.data.round);
        if (idx >= 0) p.review.rounds[idx] = summary;
        else p.review.rounds.push(summary);
        p.review.latestScore = ev.data.score;
        const round = makeReviewRound(p.track, ev.data.round, ev.data.score);
        round.score = ev.data.score;
        round.byJudge = ev.data.byJudge;
        round.hardGatePassed = ev.data.hardGatePassed;
        round.converged = ev.data.hardGatePassed;
        const existing = this.reviews.get(id) ?? { latest: null, rounds: [], machine: null };
        const rounds = existing.rounds.filter((r) => r.round !== ev.data.round).concat(round).sort((a, b) => a.round - b.round);
        this.reviews.set(id, { latest: round, rounds, machine: round.machine ?? null });
        break;
      }
      case 'file': {
        const list = this.extraFiles.get(id) ?? [];
        const i = list.findIndex((f) => f.path === ev.data.path);
        const spec = { path: ev.data.path, size: ev.data.size ?? 0, mtime: ev.ts };
        if (ev.data.op === 'delete') {
          if (i >= 0) list.splice(i, 1);
        } else if (i >= 0) list[i] = spec;
        else list.push(spec);
        this.extraFiles.set(id, list);
        for (const [re, key] of ARTIFACT_RULES) if (re.test(ev.data.path)) p.artifacts[key] = ev.data.path;
        if (ev.data.path === '03_财务/numbers.json' && !this.numbers.has(id)) {
          this.numbers.set(id, clone(p.track === 'startup' || p.track === 'case' ? MOCK_NUMBERS_STARTUP : MOCK_NUMBERS_BRAND));
        }
        break;
      }
      case 'done': {
        if (ev.data.subtype === 'stopped') p.status = 'stopped';
        if (ev.data.subtype === 'error') p.status = 'failed';
        break;
      }
      default:
        break;
    }
    p.updatedAt = ev.ts;
    if (!MOCK_PROJECTS.some((s) => s.id === id)) this.persist();
  }
}

class MockRunner {
  private logs = new Map<string, AgentEvent[]>();
  private seqs = new Map<string, number>();
  private listeners = new Map<string, Set<(ev: AgentEvent) => void>>();
  private timers = new Map<string, ReturnType<typeof setTimeout>>();
  private modes = new Map<string, 'run' | 'resume' | 'modify'>();

  constructor(private store: MockStore) {}

  log(id: string): AgentEvent[] {
    let l = this.logs.get(id);
    if (!l) {
      l = [];
      this.logs.set(id, l);
    }
    return l;
  }

  isRunning(id: string): boolean {
    return this.timers.has(id);
  }

  mode(id: string): 'run' | 'resume' | 'modify' | undefined {
    return this.modes.get(id);
  }

  private nextSeq(id: string): number {
    const n = (this.seqs.get(id) ?? 0) + 1;
    this.seqs.set(id, n);
    return n;
  }

  emit(id: string, step: Step['ev'], ts = new Date().toISOString(), notify = true): AgentEvent {
    const ev = { ...step, seq: this.nextSeq(id), ts, projectId: id } as AgentEvent;
    this.log(id).push(ev);
    this.store.apply(id, ev);
    if (notify) for (const fn of this.listeners.get(id) ?? []) fn(ev);
    return ev;
  }

  /** 为已完成阶段补历史事件（时间戳均匀分布在阶段起止之间） */
  seedHistory(project: Project, stages: StageId[]): void {
    for (const stage of stages) {
      const st = project.stages[stage];
      const start = st.startedAt ? Date.parse(st.startedAt) : Date.now() - 3600_000;
      const end = st.endedAt ? Date.parse(st.endedAt) : start + 600_000;
      const steps = buildStageScript(project, stage);
      steps.forEach((s, i) => {
        const ts = new Date(start + ((end - start) * (i + 1)) / (steps.length + 1)).toISOString();
        const ev = { ...s.ev, seq: this.nextSeq(project.id), ts, projectId: project.id } as AgentEvent;
        this.log(project.id).push(ev);
      });
    }
  }

  start(id: string, steps: Step[], mode: 'run' | 'resume' | 'modify', initialDelay = 400): void {
    this.stopTimer(id);
    this.modes.set(id, mode);
    let i = 0;
    const tick = () => {
      if (i >= steps.length) {
        this.timers.delete(id);
        this.modes.delete(id);
        return;
      }
      const step = steps[i++];
      this.emit(id, step.ev);
      const delay = step.delay ?? 300 + Math.random() * 500;
      this.timers.set(id, setTimeout(tick, delay));
    };
    this.timers.set(id, setTimeout(tick, initialDelay));
  }

  private stopTimer(id: string): void {
    const t = this.timers.get(id);
    if (t) clearTimeout(t);
    this.timers.delete(id);
  }

  stop(id: string): void {
    if (!this.isRunning(id)) return;
    this.stopTimer(id);
    this.modes.delete(id);
    const p = this.store.projects.get(id);
    const running = p ? STAGE_IDS.find((s) => p.stages[s].status === 'running') : undefined;
    if (running) this.emit(id, { type: 'stage', agent: { id: 'main', name: '主代理' }, data: { stage: running, status: 'failed', note: '用户手动停止' } });
    this.emit(id, { type: 'done', agent: { id: 'main', name: '主代理' }, data: { subtype: 'stopped', result: '已由用户停止' } });
    this.emit(id, { type: 'status', agent: { id: 'main', name: '主代理' }, data: { status: 'stopped', message: '已停止' } });
  }

  subscribe(id: string, since: number, handlers: EventHandlers): () => void {
    const fn = handlers.onEvent;
    handlers.onStatus?.('connecting');
    const t = setTimeout(() => {
      handlers.onStatus?.('open');
      for (const ev of this.log(id)) if (ev.seq > since) fn(ev);
      let set = this.listeners.get(id);
      if (!set) {
        set = new Set();
        this.listeners.set(id, set);
      }
      set.add(fn);
    }, 60);
    return () => {
      clearTimeout(t);
      this.listeners.get(id)?.delete(fn);
      handlers.onStatus?.('closed');
    };
  }
}

export function createMockApi(): BizgoldApi {
  const store = new MockStore();
  const runner = new MockRunner(store);

  // 运行中的演示项目：补齐前 5 个阶段的历史事件，再从「撰写」阶段实时续播
  const running = store.projects.get(RUNNING_PROJECT_ID);
  if (running) {
    runner.seedHistory(running, ['parse', 'research', 'strategy', 'finance', 'charts']);
    store.runStartedAt.set(running.id, Date.now() - running.usage.durationMs);
    // 跳过开头的 status / 「收到任务书」文本，直接从「撰写」阶段的 update_progress 开始续播
    const steps = buildRunScript(running, 'write', 'run');
    const firstTool = steps.findIndex((s) => s.ev.type === 'tool_use');
    runner.start(running.id, steps.slice(Math.max(0, firstTool)), 'run', 1500);
  }

  const maskKey = (k: string): string => (k && !k.startsWith('sk-****') ? `sk-****${k.slice(-4)}` : k);

  const api: BizgoldApi = {
    async health() {
      await latency();
      return clone(MOCK_HEALTH);
    },
    async getSettings() {
      await latency();
      return clone({ ...store.settings, api: { ...store.settings.api, apiKey: maskKey(store.settings.api.apiKey), apiKeySource: store.settings.api.apiKey ? 'settings' : 'none' } });
    },
    async putSettings(settings) {
      await latency();
      const next = clone(settings);
      if (next.api.apiKey.startsWith('sk-****') || next.api.apiKey === '') next.api.apiKey = store.settings.api.apiKey;
      for (const k of ['bochaKey', 'tavilyKey', 'serperKey'] as const) {
        if (next.search[k].startsWith('sk-****')) next.search[k] = store.settings.search[k];
      }
      store.settings = next;
      saveJson(LS_SETTINGS, next);
      return api.getSettings();
    },
    async testSettings(body) {
      await sleep(700 + Math.random() * 500);
      const url = body.baseUrl || store.settings.api.baseUrl;
      if (/bad|invalid|错误/.test(url)) return { ok: false, latencyMs: 0, models: [], error: `无法连接 ${url}：ENOTFOUND；备用地址同样失败` };
      return { ok: true, latencyMs: Math.round(180 + Math.random() * 240), models: [...MOCK_MODELS], baseUrlUsed: url };
    },
    async getTracks() {
      await latency();
      return clone(MOCK_TRACKS);
    },
    async getKnowledge() {
      await latency();
      return clone(MOCK_KNOWLEDGE_DOCS);
    },
    async getKnowledgeDoc(slug) {
      await latency();
      const doc = knowledgeContent(slug);
      if (!doc) throw new ApiError(`知识库文档不存在：${slug}`, 'not_found', 404);
      return doc;
    },
    async getProjects() {
      await latency();
      return clone([...store.projects.values()].sort((a, b) => Date.parse(b.updatedAt) - Date.parse(a.updatedAt)));
    },
    async createProject(input) {
      await latency();
      const now = new Date();
      const pad = (n: number) => String(n).padStart(2, '0');
      const id = `p_${now.getFullYear()}${pad(now.getMonth() + 1)}${pad(now.getDate())}_${pad(now.getHours())}${pad(now.getMinutes())}${pad(now.getSeconds())}_${Math.random().toString(16).slice(2, 6)}`;
      const project: Project = {
        id,
        name: input.name,
        teamName: input.teamName,
        track: input.track,
        group: input.group,
        subject: input.subject,
        team: input.team,
        params: { ...input.params, models: { ...DEFAULT_MODELS, ...input.params.models } },
        inputs: [],
        status: 'idle',
        currentStage: null,
        stages: emptyStages(),
        usage: { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, turns: 0, costEstimateUsd: 0, durationMs: 0 },
        review: { rounds: [], latestScore: null },
        artifacts: {},
        createdAt: now.toISOString(),
        updatedAt: now.toISOString(),
      };
      store.projects.set(id, project);
      store.persist();
      return clone(project);
    },
    async getProject(id) {
      await latency();
      return clone(store.get(id));
    },
    async deleteProject(id) {
      await latency();
      store.get(id);
      runner.stop(id);
      store.projects.delete(id);
      store.persist();
    },
    async patchProject(id, patch) {
      await latency();
      const p = store.get(id);
      if (p.status === 'running') {
        if (patch.name) p.name = patch.name;
      } else {
        Object.assign(p, clone(patch));
      }
      p.updatedAt = new Date().toISOString();
      store.persist();
      return clone(p);
    },
    async uploadFiles(id, files) {
      const p = store.get(id);
      await sleep(500 + files.length * 350);
      const added: ProjectInput[] = files.map((f) => {
        const kind = kindOf(f.name);
        const textual = kind === 'pdf' || kind === 'docx' || kind === 'xlsx' || kind === 'csv' || kind === 'text';
        return {
          name: f.name,
          path: `inputs/${f.name}`,
          size: f.size,
          kind,
          preread: textual ? `inputs/_preread/${f.name}.txt` : undefined,
          prereadStatus: 'done',
          prereadChars: textual ? Math.round(f.size / 40) : undefined,
        };
      });
      for (const a of added) {
        const i = p.inputs.findIndex((x) => x.name === a.name);
        if (i >= 0) p.inputs[i] = a;
        else p.inputs.push(a);
      }
      p.updatedAt = new Date().toISOString();
      store.persist();
      return clone(p.inputs);
    },
    async deleteFile(id, name) {
      await latency();
      const p = store.get(id);
      p.inputs = p.inputs.filter((x) => x.name !== name);
      store.persist();
      return clone(p.inputs);
    },
    async run(id, body) {
      await latency();
      const p = store.get(id);
      if (runner.isRunning(id) || p.status === 'running') throw new ApiError('项目正在运行中', 'conflict', 409);
      let from: StageId = 'parse';
      if (body.mode === 'resume') {
        from = body.fromStage ?? STAGE_IDS.find((s) => p.stages[s].status !== 'done' && p.stages[s].status !== 'skipped') ?? 'package';
        p.stages[from] = { status: 'pending' };
      } else {
        p.stages = emptyStages();
        p.review = { rounds: [], latestScore: null };
        p.artifacts = {};
        p.usage = { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, turns: 0, costEstimateUsd: 0, durationMs: 0 };
        store.reviews.delete(id);
        store.extraFiles.delete(id);
      }
      p.status = 'running';
      p.currentStage = from;
      p.session = { sessionId: p.session?.sessionId ?? `sess_${Math.random().toString(36).slice(2, 8)}`, lastPromptAt: new Date().toISOString(), mode: body.mode };
      p.updatedAt = new Date().toISOString();
      store.runStartedAt.set(id, Date.now());
      store.persist();
      runner.start(id, buildRunScript(p, from, body.mode), body.mode);
      return { ok: true };
    },
    async stop(id) {
      await latency();
      store.get(id);
      if (!runner.isRunning(id)) throw new ApiError('项目当前没有在运行', 'not_running', 409);
      runner.stop(id);
      return { ok: true };
    },
    async chat(id, message) {
      await latency();
      const p = store.get(id);
      if (runner.isRunning(id)) throw new ApiError('项目正在运行中，请等待完成或先停止后再发送修改指令', 'conflict', 409);
      if (!message.trim()) throw new ApiError('指令不能为空', 'bad_request', 400);
      p.session = { ...p.session, lastPromptAt: new Date().toISOString(), mode: 'modify' };
      runner.start(id, buildModifyScript(p, message.trim()), 'modify');
      return { ok: true };
    },
    async getArtifacts(id) {
      await latency();
      const p = store.get(id);
      return buildTree(store.files(id), p.name) ?? artifactsFor(p);
    },
    async getPreview(id, path) {
      await sleep(150 + Math.random() * 250);
      const p = store.get(id);
      return clone(previewFor(p, path));
    },
    fileUrl(id, path) {
      const p = store.projects.get(id);
      if (!p) return '#';
      return mockFileUrl(p, path);
    },
    async getReview(id) {
      await latency();
      store.get(id);
      return clone(store.reviews.get(id) ?? { latest: null, rounds: [], machine: null });
    },
    async getNumbers(id) {
      await latency();
      store.get(id);
      return clone(store.numbers.get(id) ?? { generatedAt: '', unit: '万元', entries: [], consistency: { consistent: true, mismatches: [] } });
    },
    async repackage(id) {
      const p = store.get(id);
      await sleep(1200);
      const list = store.extraFiles.get(id) ?? [];
      if (!list.some((f) => f.path === '08_提交包/提交包.zip') && !(PROJECT_FILES[id] ?? []).some((f) => f.path === '08_提交包/提交包.zip')) {
        list.push({ path: '08_提交包/提交包.zip', size: 9_410_200, mtime: new Date().toISOString() });
        store.extraFiles.set(id, list);
      }
      p.artifacts.package = '08_提交包/提交包.zip';
      p.updatedAt = new Date().toISOString();
      return { path: '08_提交包/提交包.zip' };
    },
    async openFolder(id) {
      await latency();
      store.get(id);
      return { ok: true };
    },
    subscribeEvents(id, since, handlers) {
      return runner.subscribe(id, since, handlers);
    },
  };
  return api;
}
