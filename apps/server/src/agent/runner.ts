import path from 'node:path';
import { query, type Options } from '@anthropic-ai/claude-agent-sdk';
import { IS_WINDOWS } from '../config.js';
import type { AgentEventInput, Project, ProjectStatus, Settings, StageId, StageStatus } from '../types.js';
import { ApiError, STAGE_BUDGET_MIN, STAGE_IDS, STAGE_NAMES } from '../types.js';
import type { StageBudget } from './mcp.js';
import { createLogger, errorMessage, redact } from '../util/log.js';
import type { ProjectStore } from '../services/store.js';
import type { SettingsStore } from '../services/settings.js';
import { detectGitBash, type PythonRunner } from '../services/python.js';
import { SearchService } from '../services/search.js';
import type { TrackService } from '../services/tracks.js';
import { scanArtifacts, syncReviewIntoProject } from '../services/artifacts.js';
import { pickBaseUrl } from '../services/gateway.js';
import { EventNormalizer, MAIN_AGENT } from './events.js';
import { createBizgoldMcpServer } from './mcp.js';
import { buildPrompt, buildSystemAppend } from './prompt.js';
import { shouldApplyFallback } from './progress.js';
import { prepareWorkspace, writePythonShims } from './workspace.js';

const log = createLogger('runner');
/** 主代理提前收尾时自动续跑的最大次数（每次续跑都是完整的新 query，带已完成阶段清单） */
const AUTO_CONTINUE_MAX = 6;

export interface StartOptions {
  mode: 'run' | 'resume' | 'modify';
  fromStage?: StageId;
  message?: string;
  /** 冒烟：短提示词 + 小 maxTurns，用于验证链路 */
  smoke?: boolean;
  maxTurns?: number;
  /** 服务端自动续跑的次数（主代理提前收尾时由 runner 自动 resume，用户手动启动时为 0） */
  autoContinue?: number;
}

interface ActiveRun {
  projectId: string;
  ws: string;
  abort: AbortController;
  startedAt: number;
  lastActivity: number;
  stopRequested: boolean;
  timeoutReason: string | null;
  chain: Promise<void>;
  normalizer: EventNormalizer;
  stderrTail: string[];
  timers: NodeJS.Timeout[];
  mode: StartOptions['mode'];
  smoke: boolean;
  autoContinue: number;
}

export interface RunnerDeps {
  store: ProjectStore;
  settings: SettingsStore;
  python: PythonRunner;
  tracks: TrackService;
}

/**
 * 封装 query()：准备工作区 → 构造 Options → 消费消息流 → 归一化落盘/广播 → 更新 project.json。
 * 每个项目同一时刻最多一个运行；stop() 通过 AbortController 终止子进程。
 */
export class AgentRunner {
  private active = new Map<string, ActiveRun>();
  private shuttingDown = false;

  constructor(private readonly deps: RunnerDeps) {
    const onExit = () => this.shutdown();
    process.once('SIGINT', onExit);
    process.once('SIGTERM', onExit);
    process.once('beforeExit', onExit);
  }

  isRunning(projectId: string): boolean {
    return this.active.has(projectId);
  }

  get activeCount(): number {
    return this.active.size;
  }

  /** 进程退出：中止所有运行（子进程随 abort 结束） */
  shutdown(): void {
    if (this.shuttingDown) return;
    this.shuttingDown = true;
    for (const run of this.active.values()) {
      run.stopRequested = true;
      run.abort.abort();
    }
  }

  /**
   * 启动时清理上次进程异常退出留下的 running 状态：子进程早已不在，标记为 paused（可续跑），
   * 正在进行的阶段记为 failed 并留注释；同时追加一条 status 事件让前端回放时能看到中断原因。
   */
  async recoverStale(): Promise<string[]> {
    const { store } = this.deps;
    const recovered: string[] = [];
    for (const summary of await store.list()) {
      if (summary.status !== 'running' || this.active.has(summary.id)) continue;
      await store.mutate(summary.id, (p) => {
        p.status = 'paused';
        for (const s of STAGE_IDS) {
          const st = p.stages[s];
          if (st?.status === 'running') {
            st.status = 'failed';
            st.endedAt = new Date().toISOString();
            st.note = '服务重启，运行被中断';
          }
        }
      });
      await store.appendEvent(summary.id, { type: 'status', agent: MAIN_AGENT, data: { status: 'paused', message: '服务重启，上次运行被中断；可从未完成阶段续跑' } });
      recovered.push(summary.id);
    }
    if (recovered.length) log.warn(`recovered ${recovered.length} interrupted project(s)`, recovered);
    return recovered;
  }

  async stop(projectId: string): Promise<boolean> {
    const run = this.active.get(projectId);
    if (!run) return false;
    run.stopRequested = true;
    run.normalizer.markAborted();
    run.abort.abort();
    return true;
  }

  async start(projectId: string, opts: StartOptions): Promise<void> {
    if (this.active.has(projectId)) throw new ApiError(409, 'already_running', '项目正在运行中');
    const { store, settings: settingsStore, python, tracks } = this.deps;
    const project = await store.get(projectId);
    const settings = await settingsStore.get();
    const { key: apiKey, source } = await settingsStore.apiKey();
    if (!apiKey) throw new ApiError(400, 'no_api_key', '未配置 API 密钥：请在设置页填写，或设置环境变量 BIZGOLD_API_KEY');
    if (opts.mode === 'modify' && !opts.message?.trim()) throw new ApiError(400, 'bad_request', '修改模式需要 message');
    if (opts.mode === 'resume' && !project.session.sessionId) {
      log.warn(`project ${projectId} has no session to resume, falling back to a fresh run with resume prompt`);
    }

    const ws = store.dirOf(projectId);
    const track = await tracks.get(project.track);
    const pythonPath = (await python.interpreter()) ?? (IS_WINDOWS ? 'python.exe' : 'python3');

    // 新一轮 run：重置阶段/用量/会话；resume/modify 保留
    if (opts.mode === 'run') {
      await store.mutate(projectId, (p) => {
        for (const s of STAGE_IDS) p.stages[s] = s === 'review' ? { status: 'pending', round: 0 } : { status: 'pending' };
        p.currentStage = '';
        p.usage = { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, turns: 0, costEstimateUsd: 0, durationMs: 0 };
        p.session = { sessionId: '', lastPromptAt: '', mode: 'run' };
      });
    }
    const fresh = await store.get(projectId);
    const prep = await prepareWorkspace(ws, fresh, settings.api.models, track);

    const abort = new AbortController();
    const run: ActiveRun = {
      projectId,
      ws,
      abort,
      startedAt: Date.now(),
      lastActivity: Date.now(),
      stopRequested: false,
      timeoutReason: null,
      chain: Promise.resolve(),
      normalizer: undefined as unknown as EventNormalizer,
      stderrTail: [],
      timers: [],
      mode: opts.mode,
      smoke: opts.smoke === true,
      autoContinue: opts.autoContinue ?? 0,
    };
    run.normalizer = new EventNormalizer((ev) => this.enqueue(run, ev), {
      workspaceDir: ws,
      pricing: settings.api.pricing,
      onSession: (sid) => {
        void store.mutate(projectId, (p) => {
          p.session.sessionId = sid;
          p.session.lastPromptAt = new Date().toISOString();
          p.session.mode = opts.mode;
        });
      },
    });
    this.active.set(projectId, run);

    await store.mutate(projectId, (p) => {
      p.status = 'running';
      p.session.mode = opts.mode;
      p.session.lastPromptAt = new Date().toISOString();
    });
    await this.enqueue(run, {
      type: 'status',
      agent: MAIN_AGENT,
      data: { status: 'running', message: `开始${opts.smoke ? '冒烟' : opts.mode === 'run' ? '运行' : opts.mode === 'resume' ? '续跑' : '修改'}（模型 ${settings.api.models.main}，密钥来源 ${source}）`, warnings: prep.warnings },
    });
    for (const w of prep.warnings) await this.enqueue(run, { type: 'error', agent: MAIN_AGENT, data: { message: w, recoverable: true } });

    // 网关地址选择（主/备用）
    const gw = await pickBaseUrl(settings.api.baseUrl, settings.api.fallbackBaseUrl, apiKey);
    if (gw.note) await this.enqueue(run, { type: 'status', agent: MAIN_AGENT, data: { status: 'running', message: gw.note } });

    const gitBash = await detectGitBash();
    const search = new SearchService(async () => (await settingsStore.get()).search, ws);
    const mcp = createBizgoldMcpServer({
      projectId,
      workspaceDir: ws,
      search,
      python,
      updateStage: (stage, status, note, round) => this.updateStage(run, stage, status, note, round, 'tool'),
      budget: () => this.budgetSnapshot(run),
      onFileWritten: (rel, size) => void this.enqueue(run, { type: 'file', agent: MAIN_AGENT, data: { path: rel, op: 'write', size } }),
    });

    const sessionId = fresh.session.sessionId;
    const useResume = (opts.mode === 'resume' || opts.mode === 'modify') && Boolean(sessionId);
    const promptMode = opts.smoke ? 'smoke' : opts.mode;
    let prompt = buildPrompt({ mode: promptMode, project: fresh, fromStage: opts.fromStage, message: opts.message });
    if (!useResume && opts.mode === 'modify') prompt = '请先阅读 任务书.md 与 CLAUDE.md 了解项目背景与工作区现状。\n' + prompt;

    const maxTurns = opts.maxTurns ?? (opts.smoke ? 15 : settings.run.maxTurns);
    await writePythonShims(ws, pythonPath).catch((e) => log.warn('python shims failed', errorMessage(e)));
    const env = this.buildEnv(settings, apiKey, gw.baseUrl, ws, pythonPath, gitBash);
    const options: Options = {
      cwd: ws,
      model: settings.api.models.main,
      systemPrompt: {
        type: 'preset',
        preset: 'claude_code',
        append: buildSystemAppend({ workspaceDir: ws, pythonCmd: pythonPath, isWindows: IS_WINDOWS, projectId }),
      },
      settingSources: ['project'],
      permissionMode: 'bypassPermissions',
      allowDangerouslySkipPermissions: true,
      maxTurns,
      includePartialMessages: true,
      persistSession: true,
      strictMcpConfig: true,
      forwardSubagentText: true,
      mcpServers: { bizgold: mcp },
      disallowedTools: ['WebSearch', 'WebFetch', 'AskUserQuestion'],
      abortController: abort,
      env,
      stderr: (data: string) => {
        const line = redact(data.trimEnd());
        if (!line) return;
        run.stderrTail.push(line);
        if (run.stderrTail.length > 60) run.stderrTail.shift();
        log.debug(`[cli ${projectId}] ${line}`);
      },
      ...(useResume ? { resume: sessionId } : {}),
      ...(settings.api.disableThinking ? { thinking: { type: 'disabled' as const } } : {}),
    };

    // 超时保护：总时长 + 空闲
    const totalMs = Math.max(1, settings.run.totalTimeoutMin) * 60_000;
    const idleMs = Math.max(1, settings.run.idleTimeoutMin) * 60_000;
    run.timers.push(
      setTimeout(() => this.timeout(run, `总时长超过 ${settings.run.totalTimeoutMin} 分钟`), totalMs),
      setInterval(() => {
        if (Date.now() - run.lastActivity > idleMs) this.timeout(run, `${settings.run.idleTimeoutMin} 分钟无输出`);
      }, 30_000),
    );

    log.info(`start ${opts.mode}${opts.smoke ? ' (smoke)' : ''} project=${projectId} model=${settings.api.models.main} resume=${useResume} maxTurns=${maxTurns}`);
    void this.loop(run, prompt, options);
  }

  private buildEnv(settings: Settings, apiKey: string, baseUrl: string, ws: string, pythonPath: string, gitBash: string | null): Options['env'] {
    const m = settings.api.models;
    const env: Record<string, string | undefined> = {
      // 子代理前台同步执行：SDK 0.3.x 默认把 Agent 调用后台化（async_launched），而单次 query 模式下主代理一旦“等待通知”结束回复，整个进程就结束、子代理被杀；
      // 关掉后台任务后 Agent 调用阻塞到子代理完成并直接返回结果，同一条回复里的多个 Agent 调用仍然并行。
      CLAUDE_CODE_DISABLE_BACKGROUND_TASKS: '1',
      ...process.env,
      ANTHROPIC_BASE_URL: baseUrl,
      ANTHROPIC_AUTH_TOKEN: apiKey,
      // 两种凭据并存会触发 CLI 警告且可能发错头，这里只用 Bearer
      ANTHROPIC_API_KEY: undefined,
      ANTHROPIC_MODEL: m.main,
      ANTHROPIC_DEFAULT_HAIKU_MODEL: m.fast || m.main,
      CLAUDE_CODE_SUBAGENT_MODEL: m.main,
      CLAUDE_CONFIG_DIR: path.join(ws, '.claude-home'),
      CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC: '1',
      DISABLE_TELEMETRY: '1',
      DISABLE_ERROR_REPORTING: '1',
      DISABLE_AUTOUPDATER: '1',
      DISABLE_BUG_COMMAND: '1',
      CLAUDE_CODE_DISABLE_AUTO_MEMORY: '1',
      CLAUDE_CODE_MAX_CONTEXT_TOKENS: String(settings.api.contextTokens || 200000),
      BASH_DEFAULT_TIMEOUT_MS: String(settings.run.bashTimeoutMs || 600000),
      BASH_MAX_TIMEOUT_MS: String(Math.max(settings.run.bashTimeoutMs || 600000, 1_800_000)),
      PYTHONUTF8: '1',
      PYTHONIOENCODING: 'utf-8',
      BIZGOLD_PYTHON: pythonPath,
      BIZGOLD_WORKSPACE: ws,
    };
    // 工作区 .bin 里的 python / python3 垫片指向配置的解释器：主代理和子代理偶尔会直接敲 `python3 -c`，
    // 在 Windows 上 python3 往往是应用商店别名或根本不存在，这里保证它们都落到装好依赖的那一个解释器。
    const shimDir = path.join(ws, '.bin');
    const sep = IS_WINDOWS ? ';' : ':';
    env.PATH = `${shimDir}${sep}${process.env.PATH ?? ''}`;
    if (IS_WINDOWS) env.Path = env.PATH;
    if (settings.api.disableThinking) env.CLAUDE_CODE_DISABLE_THINKING = '1';
    // 技术手册 §4.5：网关对 beta 字段（context_management / defer_loading 等）返回 400 时的首选开关。
    // 本机实测 llm.goaichat.top 不需要（smoke 不加该开关即通过），因此默认关闭、设置可开。
    if (settings.api.disableExperimentalBetas) env.CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS = '1';
    if (IS_WINDOWS && gitBash) env.CLAUDE_CODE_GIT_BASH_PATH = gitBash;
    return env;
  }

  private timeout(run: ActiveRun, reason: string): void {
    if (run.timeoutReason || run.stopRequested) return;
    run.timeoutReason = reason;
    log.warn(`project ${run.projectId} timeout: ${reason}`);
    void this.enqueue(run, { type: 'error', agent: MAIN_AGENT, data: { message: `运行超时：${reason}，已中止`, recoverable: true } });
    run.abort.abort();
  }

  private async loop(run: ActiveRun, prompt: string, options: Options): Promise<void> {
    const { projectId } = run;
    let threw: unknown = null;
    try {
      for await (const msg of query({ prompt, options })) {
        run.lastActivity = Date.now();
        run.normalizer.handle(msg);
      }
    } catch (e) {
      threw = e;
    }
    run.normalizer.flush();
    for (const t of run.timers) clearTimeout(t);
    await run.chain.catch(() => undefined);

    const result = run.normalizer.lastResult;
    let finalStatus: ProjectStatus;
    if (run.stopRequested && !run.timeoutReason) {
      finalStatus = 'stopped';
      if (!result) await this.enqueue(run, { type: 'done', agent: MAIN_AGENT, data: { subtype: 'stopped', result: '用户已停止', usage: null } });
    } else if (run.timeoutReason) {
      finalStatus = 'failed';
      if (!result) await this.enqueue(run, { type: 'done', agent: MAIN_AGENT, data: { subtype: 'error', result: `运行超时：${run.timeoutReason}`, usage: null } });
    } else if (threw && !result) {
      finalStatus = 'failed';
      const msg = errorMessage(threw);
      const tail = run.stderrTail.slice(-8).join('\n');
      await this.enqueue(run, { type: 'error', agent: MAIN_AGENT, data: { message: `运行异常：${msg}${tail ? `\n--- stderr ---\n${tail}` : ''}`, recoverable: false } });
      await this.enqueue(run, { type: 'done', agent: MAIN_AGENT, data: { subtype: 'error', result: msg, usage: null } });
    } else if (result?.subtype === 'success') {
      finalStatus = 'completed';
    } else if (result?.subtype === 'stopped') {
      finalStatus = 'stopped';
    } else {
      finalStatus = 'failed';
      if (threw) log.warn(`query ended with error after result: ${errorMessage(threw)}`);
    }

    let shouldAutoContinue = false;
    try {
      await this.deps.store.mutate(projectId, async (p) => {
        const allDone = STAGE_IDS.every((s) => p.stages[s]?.status === 'done' || p.stages[s]?.status === 'skipped');
        const prematurely = finalStatus === 'completed' && !allDone;
        // 成功结束但阶段未全部完成（模型提前收尾/冒烟/对半成品做修改）→ paused，便于续跑
        p.status = prematurely ? 'paused' : finalStatus;
        p.usage.durationMs += Date.now() - run.startedAt;
        // 提前收尾时进行中的阶段并没有做完，回退为 pending 让续跑从它开始；真正失败才标 failed
        for (const s of STAGE_IDS) if (p.stages[s]?.status === 'running') p.stages[s]!.status = finalStatus === 'completed' ? 'pending' : 'failed';
        p.artifacts = { ...p.artifacts, ...(await scanArtifacts(run.ws)) };
        await syncReviewIntoProject(run.ws, p);
        // 完整运行（非冒烟、非修改）提前收尾：主代理常在“等待子代理 / 下一步”时结束回复，服务端自动续跑（上限 AUTO_CONTINUE_MAX 次）
        shouldAutoContinue = prematurely && !run.smoke && run.mode !== 'modify' && run.autoContinue < AUTO_CONTINUE_MAX && !!p.session.sessionId;
      });
      const p = await this.deps.store.get(projectId);
      await this.enqueue(run, {
        type: 'status',
        agent: MAIN_AGENT,
        data: { status: p.status, message: shouldAutoContinue ? `主代理提前收尾（阶段未全部完成），服务端自动续跑（第 ${run.autoContinue + 1}/${AUTO_CONTINUE_MAX} 次）` : `运行结束：${p.status}` },
      });
    } catch (e) {
      log.error('finalize failed', errorMessage(e));
    } finally {
      await run.chain.catch(() => undefined);
      this.active.delete(projectId);
      log.info(`end project=${projectId} status=${finalStatus} elapsed=${Math.round((Date.now() - run.startedAt) / 1000)}s`);
    }
    if (shouldAutoContinue) {
      setTimeout(() => {
        this.start(projectId, { mode: 'resume', autoContinue: run.autoContinue + 1 }).catch((e) => log.warn(`auto-continue failed: ${errorMessage(e)}`));
      }, 3000);
    }
  }

  /** 事件串行化：落盘 + 广播 + 同步 project.json（保持与消息流一致的顺序） */
  private enqueue(run: ActiveRun, ev: AgentEventInput): Promise<void> {
    const task = async () => {
      try {
        await this.persistAndApply(run, ev);
      } catch (e) {
        log.error('event handling failed', errorMessage(e));
      }
    };
    run.chain = run.chain.then(task, task);
    return run.chain;
  }

  private async persistAndApply(run: ActiveRun, ev: AgentEventInput): Promise<void> {
    const { store } = this.deps;
    const { projectId } = run;
    if (ev.type === 'stage' && ev.data.source === 'fallback') {
      const stage = ev.data.stage as StageId;
      const status = ev.data.status as StageStatus;
      const p = await store.get(projectId);
      if (!shouldApplyFallback(p.stages[stage]?.status, status)) return;
      await this.updateStage(run, stage, status, undefined, typeof ev.data.round === 'number' ? ev.data.round : undefined, 'fallback');
      return;
    }
    if (ev.type === 'usage') {
      const d = ev.data as { inputTokens?: number; outputTokens?: number; cacheReadTokens?: number; turns?: number; costEstimateUsd?: number };
      await store.mutate(projectId, (p) => {
        // result.modelUsage 是会话累计值：直接赋值而不是累加
        p.usage.inputTokens = d.inputTokens ?? p.usage.inputTokens;
        p.usage.outputTokens = d.outputTokens ?? p.usage.outputTokens;
        p.usage.cacheReadTokens = d.cacheReadTokens ?? p.usage.cacheReadTokens;
        p.usage.turns = d.turns ?? p.usage.turns;
        p.usage.costEstimateUsd = d.costEstimateUsd ?? p.usage.costEstimateUsd;
      });
    }
    const isDelta = ev.type === 'text_delta' || (ev.type === 'thinking' && ev.data.delta === true);
    await store.appendEvent(projectId, ev, { persist: !isDelta });
  }

  /** 当前阶段用时 / 预算（分钟），供 MCP 工具在回执里提醒主代理；没有运行中的阶段时返回 null */
  private async budgetSnapshot(run: ActiveRun): Promise<StageBudget | null> {
    const p = await this.deps.store.get(run.projectId);
    const stage = (Object.keys(p.stages) as StageId[]).find((s) => p.stages[s]?.status === 'running') ?? p.currentStage;
    const st = stage ? p.stages[stage] : undefined;
    if (!stage || !st) return null;
    const startedMs = st.startedAt ? Date.parse(st.startedAt) : run.startedAt;
    const elapsedMin = Math.max(0, Math.round((Date.now() - startedMs) / 60_000));
    const totalMin = Math.round(((p.usage?.durationMs ?? 0) + (Date.now() - run.startedAt)) / 60_000);
    return { stage, elapsedMin, budgetMin: STAGE_BUDGET_MIN[stage], totalMin };
  }

  /** update_progress 工具 / 兜底解析共用：写 project.json.stages 并广播 stage 事件；返回给主代理看的用时/预算提示 */
  private async updateStage(run: ActiveRun, stage: StageId, status: StageStatus, note: string | undefined, round: number | undefined, source: 'tool' | 'fallback'): Promise<string> {
    const { store } = this.deps;
    const now = new Date().toISOString();
    let reviewLatest: { round: number; score: number; byJudge: Record<string, number>; hardGatePassed: boolean } | null = null;
    let budgetNote = '';
    await store.mutate(run.projectId, async (p) => {
      const cur = p.stages[stage] ?? { status: 'pending' };
      const totalMin = Math.round(((p.usage?.durationMs ?? 0) + (Date.now() - run.startedAt)) / 60_000);
      const budgetMin = STAGE_BUDGET_MIN[stage];
      if (status === 'running') {
        budgetNote = `全程已用 ${totalMin} 分钟；本阶段预算 ${budgetMin} 分钟，超预算时按 CLAUDE.md 的降级规则收尾，不要为单项指标反复补派。`;
      } else if (status === 'done' || status === 'failed') {
        const startedMs = cur.startedAt ? Date.parse(cur.startedAt) : Date.now();
        const usedMin = Math.max(0, Math.round((Date.now() - startedMs) / 60_000));
        const over = usedMin - budgetMin;
        budgetNote = `本阶段用时 ${usedMin} 分钟（预算 ${budgetMin}${over > 0 ? `，超出 ${over}` : ''}）；全程已用 ${totalMin} 分钟。`;
      }
      cur.status = status;
      if (status === 'running') {
        cur.startedAt = cur.startedAt ?? now;
        delete cur.endedAt;
        p.currentStage = stage;
      } else if (status === 'done' || status === 'failed' || status === 'skipped') {
        cur.startedAt = cur.startedAt ?? now;
        cur.endedAt = now;
        if (p.currentStage === stage && status === 'done') {
          const next = STAGE_IDS[STAGE_IDS.indexOf(stage) + 1];
          p.currentStage = next ?? stage;
        }
      }
      if (note !== undefined) cur.note = note;
      if (round !== undefined) cur.round = round;
      p.stages[stage] = cur;
      if (status === 'done') {
        p.artifacts = { ...p.artifacts, ...(await scanArtifacts(run.ws)) };
        if (stage === 'review' || stage === 'package') {
          const r = await syncReviewIntoProject(run.ws, p);
          if (r.changed && r.latest) reviewLatest = r.latest;
        }
      }
    });
    await store.appendEvent(run.projectId, {
      type: 'stage',
      agent: MAIN_AGENT,
      data: { stage, status, note: note ?? '', ...(round !== undefined ? { round } : {}), source, name: STAGE_NAMES[stage] },
    });
    if (reviewLatest) {
      const r = reviewLatest as { round: number; score: number; byJudge: Record<string, number>; hardGatePassed: boolean };
      await store.appendEvent(run.projectId, { type: 'review', agent: MAIN_AGENT, data: { round: r.round, score: r.score, byJudge: r.byJudge, hardGatePassed: r.hardGatePassed } });
    }
    return budgetNote;
  }
}
