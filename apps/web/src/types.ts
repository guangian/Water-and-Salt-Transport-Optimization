// 与 docs/contracts.md 对齐的数据类型（camelCase）

export type TrackId = 'brand' | 'startup' | 'digital' | 'case' | 'mice';

export type ProjectStatus = 'idle' | 'running' | 'paused' | 'completed' | 'failed' | 'stopped';

export type StageId =
  | 'parse'
  | 'research'
  | 'strategy'
  | 'finance'
  | 'charts'
  | 'write'
  | 'pitch'
  | 'review'
  | 'package';

export type StageStatus = 'pending' | 'running' | 'done' | 'failed' | 'skipped';

export interface StageState {
  status: StageStatus;
  startedAt?: string;
  endedAt?: string;
  note?: string;
  round?: number;
}

export type SubjectType = 'brand' | 'startup' | 'case-company' | 'event';

export interface ProjectSubject {
  type: SubjectType;
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

export type ResearchMode = 'synthetic' | 'provided';

export type ModelRole =
  | 'main'
  | 'researcher'
  | 'finance'
  | 'writer'
  | 'pitch'
  | 'judgeA'
  | 'judgeB'
  | 'judgeC'
  | 'inspector'
  | 'fast';

export type ModelMap = Record<ModelRole, string>;

export interface ProjectParams {
  researchMode: ResearchMode;
  pageTarget: number;
  pitchMinutes: number;
  reviewRounds: number;
  targetScore: number;
  language: string;
  models: Partial<ModelMap>;
}

export type InputKind = 'pdf' | 'docx' | 'xlsx' | 'csv' | 'image' | 'text' | 'other';

export interface ProjectInput {
  name: string;
  path: string;
  size: number;
  kind: InputKind;
  preread?: string;
  prereadStatus?: 'pending' | 'done' | 'failed';
  prereadChars?: number;
}

export interface ProjectUsage {
  inputTokens: number;
  outputTokens: number;
  cacheReadTokens: number;
  turns: number;
  costEstimateUsd: number;
  durationMs: number;
}

export type JudgeId = 'investor' | 'professor' | 'executive';

export interface ReviewRoundSummary {
  round: number;
  score: number;
  byJudge: Partial<Record<JudgeId, number>>;
  hardGatePassed: boolean;
  path?: string;
}

export type ArtifactKey = 'book' | 'pdf' | 'deck' | 'model' | 'script' | 'qa' | 'package';

export interface ProjectSession {
  sessionId?: string;
  lastPromptAt?: string;
  mode?: 'run' | 'modify' | 'resume';
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
  currentStage: StageId | null;
  stages: Record<StageId, StageState>;
  session?: ProjectSession;
  usage: ProjectUsage;
  review: { rounds: ReviewRoundSummary[]; latestScore: number | null };
  artifacts: Partial<Record<ArtifactKey, string>>;
  createdAt: string;
  updatedAt: string;
}

export type ProjectCreateInput = Pick<
  Project,
  'name' | 'teamName' | 'track' | 'group' | 'subject' | 'team' | 'params'
>;

export interface Track {
  id: TrackId;
  name: string;
  docName: string;
  pageLimit: number;
  wordMin: number;
  scoring: { written: number; pitch: number };
  elements: string[];
  groups: string[];
  assignedCases: { company: string; title: string }[];
  description?: string;
  deliverables?: string[];
}

export interface Health {
  ok: boolean;
  version: string;
  python: { found: boolean; path: string; version: string; missing: string[] };
  pdfConverter: 'word' | 'wps' | 'libreoffice' | null;
  fonts: { cjk: string[] };
  gitBash: boolean;
  workspaceDir: string;
}

export type SearchProvider = 'bocha' | 'tavily' | 'serper' | 'bing' | 'duckduckgo';

export interface Settings {
  api: {
    baseUrl: string;
    apiKey: string;
    apiKeySource?: 'settings' | 'env' | 'none';
    fallbackBaseUrl: string;
    models: ModelMap;
    contextTokens: number;
    disableThinking: boolean;
  };
  search: {
    providers: SearchProvider[];
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

export interface SettingsTestResult {
  ok: boolean;
  latencyMs: number;
  models: string[];
  error?: string;
  baseUrlUsed?: string;
}

// ---------- 事件 ----------

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

export interface EventDataMap {
  status: { status: ProjectStatus; message?: string };
  stage: { stage: StageId; status: StageStatus; note?: string; round?: number };
  text: { text: string };
  text_delta: { text: string };
  thinking: { text: string; delta: boolean };
  tool_use: { toolUseId: string; name: string; input: unknown };
  tool_result: { toolUseId: string; name: string; isError: boolean; output: string };
  file: { path: string; op: 'write' | 'edit' | 'delete'; size?: number };
  subagent: { toolUseId: string; name: string; status: 'start' | 'end'; summary?: string };
  usage: { inputTokens: number; outputTokens: number; cacheReadTokens: number; turns: number; costEstimateUsd: number };
  review: { round: number; score: number; byJudge: Partial<Record<JudgeId, number>>; hardGatePassed: boolean };
  error: { message: string; recoverable: boolean };
  done: { subtype: 'success' | 'error' | 'stopped'; result?: string; usage?: Partial<ProjectUsage> };
}

export type AgentEvent = {
  [K in AgentEventType]: {
    seq: number;
    ts: string;
    projectId: string;
    type: K;
    agent: AgentRef;
    data: EventDataMap[K];
  };
}[AgentEventType];

export type AgentEventOf<K extends AgentEventType> = Extract<AgentEvent, { type: K }>;

// ---------- 产出物 ----------

export interface ArtifactNode {
  name: string;
  path: string;
  type: 'dir' | 'file';
  size?: number;
  mtime?: string;
  children?: ArtifactNode[];
}

export type PreviewKind = 'markdown' | 'image' | 'slides' | 'table' | 'pdf' | 'text' | 'binary';

export interface SheetData {
  name: string;
  rows: (string | number | null)[][];
}

export interface Preview {
  kind: PreviewKind;
  content?: string;
  url?: string;
  sheets?: SheetData[];
  slides?: SlidesDeck;
  size?: number;
  message?: string;
}

// ---------- slides.json ----------

export type SlideLayout =
  | 'cover'
  | 'agenda'
  | 'section'
  | 'bullets'
  | 'image'
  | 'two_column'
  | 'chart'
  | 'table'
  | 'kpi'
  | 'timeline'
  | 'team'
  | 'quote'
  | 'closing'
  | 'appendix';

export interface SlideChart {
  type: 'bar' | 'line' | 'pie' | 'column_clustered';
  categories: string[];
  series: { name: string; values: number[] }[];
  unit?: string;
}

export interface Slide {
  layout: SlideLayout;
  title?: string;
  subtitle?: string;
  kicker?: string;
  items?: string[];
  bullets?: string[];
  highlight?: string;
  highlightLabel?: string;
  image?: string;
  caption?: string;
  left?: { heading: string; bullets: string[] };
  right?: { heading: string; bullets: string[] };
  chart?: SlideChart;
  insight?: string;
  columns?: string[];
  rows?: string[][];
  kpis?: { value: string; label: string }[];
  steps?: { when: string; what: string }[];
  members?: { role: string; duty: string }[];
  quote?: string;
  by?: string;
  contact?: string;
  note?: string;
  duration?: number;
}

export interface SlidesDeck {
  title: string;
  subtitle?: string;
  aspect?: string;
  theme?: string;
  minutes?: number;
  footer?: string;
  slides: Slide[];
  output?: string;
}

// ---------- 评审 ----------

export interface JudgeScore {
  dimension: string;
  weight: number;
  score: number;
  cap: number;
  evidence?: string;
  issues?: string[];
  fixes?: string[];
}

export interface HardGate {
  id?: string;
  gate?: string;
  passed: boolean;
  detail?: string;
  hits?: string[];
  pages?: number;
  limit?: number;
}

export interface ReviewFix {
  priority: number;
  where: string;
  what: string;
  gain: string;
}

export interface JudgeReview {
  round: number;
  judge: JudgeId;
  model: string;
  track: TrackId;
  rubric?: string;
  scores: JudgeScore[];
  total: number;
  officialEstimate?: { written: number; pitch: number; total: number };
  hardGates: HardGate[];
  topFixes: ReviewFix[];
  verdict: string;
}

export type CheckStatus = 'pass' | 'warn' | 'fail' | 'skipped';

export interface MachineCheck {
  id: string;
  name?: string;
  value: number | string | boolean | null;
  threshold?: { pass: number; warn: number };
  status: CheckStatus;
  detail: string;
}

export interface CrossDocMismatch {
  key: string;
  label?: string;
  book?: number | string | null;
  deck?: number | string | null;
  script?: number | string | null;
  qa?: number | string | null;
}

export interface MachineCheckReport {
  track: TrackId;
  generatedAt: string;
  checks: MachineCheck[];
  hardGates: HardGate[];
  crossDoc: { consistent: boolean; mismatches: CrossDocMismatch[] };
  summary: { pass: number; warn: number; fail: number; skipped?: number };
}

export interface ReviewRound {
  round: number;
  score: number;
  byJudge: Partial<Record<JudgeId, number>>;
  hardGatePassed: boolean;
  converged: boolean;
  fixes: ReviewFix[];
  judges: JudgeReview[];
  machine?: MachineCheckReport;
  path?: string;
  at?: string;
}

export interface ReviewData {
  latest: ReviewRound | null;
  rounds: ReviewRound[];
  machine: MachineCheckReport | null;
}

// ---------- 数字账本 ----------

export interface NumberEntry {
  key: string;
  label: string;
  value: number;
  unit: string;
  source?: string;
  usedIn?: string[];
  synthetic?: boolean;
  consistency?: 'ok' | 'mismatch' | 'unchecked';
}

export interface NumbersData {
  generatedAt: string;
  unit: string;
  entries: NumberEntry[];
  consistency: { consistent: boolean; mismatches: CrossDocMismatch[] };
}

// ---------- 知识库 ----------

export interface KnowledgeDoc {
  slug: string;
  title: string;
  group: string;
  summary?: string;
  updatedAt?: string;
}

export interface KnowledgeContent {
  slug: string;
  title: string;
  content: string;
}

export interface ApiErrorBody {
  error: { code: string; message: string };
}
