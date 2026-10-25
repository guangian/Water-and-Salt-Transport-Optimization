/**
 * 契约类型（docs/contracts.md §1、§2、§4）。字段 camelCase，时间 ISO-8601，路径相对工作区根。
 */

export type TrackId = 'brand' | 'startup' | 'digital' | 'case' | 'mice';

export const STAGE_IDS = [
  'parse',
  'research',
  'strategy',
  'finance',
  'charts',
  'write',
  'pitch',
  'review',
  'package',
] as const;
export type StageId = (typeof STAGE_IDS)[number];

export const STAGE_NAMES: Record<StageId, string> = {
  parse: '读题解析',
  research: '调研',
  strategy: '策略',
  finance: '财务',
  charts: '图表',
  write: '撰写',
  pitch: '路演',
  review: '评审',
  package: '打包',
};

export type StageStatus = 'pending' | 'running' | 'done' | 'failed' | 'skipped';
export type ProjectStatus = 'idle' | 'running' | 'paused' | 'completed' | 'failed' | 'stopped';
export type RunMode = 'run' | 'resume' | 'modify';

// ---------- 设置（契约 §1） ----------

export interface ModelMap {
  main: string;
  researcher: string;
  finance: string;
  writer: string;
  pitch: string;
  judgeA: string;
  judgeB: string;
  judgeC: string;
  inspector: string;
  fast: string;
}

/** 价格表：USD / 百万 token；缺省 0（第三方网关的 total_cost_usd 不可信，自行计价） */
export interface ModelPrice {
  input: number;
  output: number;
  cacheRead?: number;
}

export interface Settings {
  api: {
    baseUrl: string;
    apiKey: string;
    fallbackBaseUrl: string;
    models: ModelMap;
    contextTokens: number;
    disableThinking: boolean;
    /** 网关拒绝 beta 字段（400）时置 true：子进程加 CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1 */
    disableExperimentalBetas: boolean;
    pricing: Record<string, ModelPrice>;
  };
  search: {
    providers: SearchProviderId[];
    bochaKey: string;
    tavilyKey: string;
    serperKey: string;
    maxResults: number;
  };
  python: { path: string; pipIndex: string };
  run: {
    maxTurns: number;
    totalTimeoutMin: number;
    idleTimeoutMin: number;
    bashTimeoutMs: number;
    reviewRounds: number;
    targetScore: number;
  };
  paths: { workspaceDir: string };
  ui: { theme: 'light' | 'dark'; language: string };
}

export type SearchProviderId = 'bocha' | 'tavily' | 'serper' | 'so360' | 'bing' | 'duckduckgo';

export type ApiKeySource = 'settings' | 'env' | 'none';

/** GET /api/settings 的返回：apiKey 脱敏 + 来源 */
export interface PublicSettings extends Omit<Settings, 'api'> {
  api: Settings['api'] & { apiKeySource: ApiKeySource };
}

// ---------- 项目（契约 §2） ----------

export interface ProjectSubject {
  type: 'brand' | 'startup' | 'case-company' | 'event';
  name: string;
  company: string;
  isAssignedCase: boolean;
  caseTitle: string;
  brief: string;
  region: string;
  industry: string;
}

export interface TeamMember {
  role: string;
  major: string;
  strength: string;
}

export interface ProjectParams {
  researchMode: 'synthetic' | 'provided';
  pageTarget: number;
  pitchMinutes: number;
  reviewRounds: number;
  targetScore: number;
  language: string;
  models: Partial<ModelMap>;
}

export interface ProjectInput {
  name: string;
  path: string;
  size: number;
  kind: string;
  preread?: string;
}

export interface StageState {
  status: StageStatus;
  startedAt?: string;
  endedAt?: string;
  note?: string;
  round?: number;
}

export interface ProjectUsage {
  inputTokens: number;
  outputTokens: number;
  cacheReadTokens: number;
  turns: number;
  costEstimateUsd: number;
  durationMs: number;
}

export interface ReviewRoundSummary {
  round: number;
  score: number;
  byJudge: Record<string, number>;
  hardGatePassed: boolean;
  path: string;
}

export interface Project {
  id: string;
  name: string;
  teamName: string;
  track: TrackId;
  group: string;
  subject: ProjectSubject;
  team: TeamMember[];
  params: ProjectParams;
  inputs: ProjectInput[];
  status: ProjectStatus;
  currentStage: StageId | '';
  stages: Record<StageId, StageState>;
  session: { sessionId: string; lastPromptAt: string; mode: RunMode | '' };
  usage: ProjectUsage;
  review: { rounds: ReviewRoundSummary[]; latestScore: number | null };
  artifacts: Record<string, string>;
  createdAt: string;
  updatedAt: string;
}

export type ProjectCreateInput = Pick<Project, 'name' | 'teamName' | 'track'> &
  Partial<Pick<Project, 'group' | 'subject' | 'team' | 'params'>>;

export interface ProjectSummary {
  id: string;
  name: string;
  teamName: string;
  track: TrackId;
  group: string;
  status: ProjectStatus;
  currentStage: StageId | '';
  latestScore: number | null;
  stages: Record<StageId, StageStatus>;
  createdAt: string;
  updatedAt: string;
}

// ---------- 事件（契约 §4） ----------

export type AgentEventType =
  | 'status'
  | 'stage'
  | 'text'
  | 'text_delta'
  | 'thinking'
  | 'tool_use'
  | 'tool_result'
  | 'file'
  | 'subagent'
  | 'usage'
  | 'review'
  | 'error'
  | 'done';

export interface AgentRef {
  id: string;
  name: string;
}

export interface AgentEventInput {
  type: AgentEventType;
  agent: AgentRef;
  data: Record<string, unknown>;
}

export interface AgentEvent extends AgentEventInput {
  seq: number;
  ts: string;
  projectId: string;
}

// ---------- 赛道 / 知识库 ----------

export interface TrackInfo {
  id: TrackId;
  name: string;
  docName: string;
  pageLimit: number;
  wordMin: number;
  scoring: { written: number; pitch: number };
  elements: string[];
  groups: string[];
  assignedCases: { company: string; title: string }[];
  notes?: string;
}

export interface KnowledgeEntry {
  slug: string;
  title: string;
  kind: 'knowledge' | 'track';
  path: string;
  size: number;
  mtime: string;
}

// ---------- 检索 ----------

export interface SearchResult {
  title: string;
  url: string;
  snippet: string;
  date?: string;
  provider: SearchProviderId;
}

// ---------- 错误 ----------

export class ApiError extends Error {
  constructor(
    public readonly statusCode: number,
    public readonly code: string,
    message: string,
  ) {
    super(message);
    this.name = 'ApiError';
  }
}
