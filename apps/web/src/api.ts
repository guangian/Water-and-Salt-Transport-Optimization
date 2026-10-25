import type {
  AgentEvent,
  ApiErrorBody,
  ArtifactNode,
  Health,
  KnowledgeContent,
  KnowledgeDoc,
  NumbersData,
  Preview,
  Project,
  ProjectCreateInput,
  ProjectInput,
  ReviewData,
  Settings,
  SettingsTestResult,
  StageId,
  Track,
} from './types';

import { ApiError, toApiError } from './lib/errors';
import { createMockApi } from './mock';

export { ApiError, toApiError };

export const IS_MOCK = import.meta.env.VITE_MOCK === '1';
export const API_BASE = import.meta.env.VITE_API_BASE ?? '/api';

export type SseStatus = 'connecting' | 'open' | 'reconnecting' | 'closed';

export interface EventHandlers {
  onEvent: (ev: AgentEvent) => void;
  onStatus?: (status: SseStatus, info?: { retryInMs?: number; attempt?: number }) => void;
}

export interface RequestOptions {
  /** 出错时不弹 toast */
  silent?: boolean;
}

export interface BizgoldApi {
  health(opts?: RequestOptions): Promise<Health>;
  getSettings(): Promise<Settings>;
  putSettings(settings: Settings): Promise<Settings>;
  testSettings(body: { baseUrl?: string; apiKey?: string }): Promise<SettingsTestResult>;
  getTracks(): Promise<Track[]>;
  getKnowledge(): Promise<KnowledgeDoc[]>;
  getKnowledgeDoc(slug: string): Promise<KnowledgeContent>;
  getProjects(): Promise<Project[]>;
  createProject(input: ProjectCreateInput): Promise<Project>;
  getProject(id: string, opts?: RequestOptions): Promise<Project>;
  deleteProject(id: string): Promise<void>;
  patchProject(id: string, patch: Partial<ProjectCreateInput>): Promise<Project>;
  uploadFiles(id: string, files: File[]): Promise<ProjectInput[]>;
  deleteFile(id: string, name: string): Promise<ProjectInput[]>;
  run(id: string, body: { mode: 'run' | 'resume'; fromStage?: StageId }): Promise<{ ok: boolean }>;
  stop(id: string): Promise<{ ok: boolean }>;
  chat(id: string, message: string): Promise<{ ok: boolean }>;
  getArtifacts(id: string): Promise<ArtifactNode>;
  getPreview(id: string, path: string): Promise<Preview>;
  fileUrl(id: string, path: string): string;
  getReview(id: string): Promise<ReviewData>;
  getNumbers(id: string): Promise<NumbersData>;
  repackage(id: string): Promise<{ path: string }>;
  openFolder(id: string): Promise<{ ok: boolean }>;
  subscribeEvents(id: string, since: number, handlers: EventHandlers): () => void;
}

// ---------- 统一错误通知 ----------

type ErrorListener = (err: ApiError) => void;
const errorListeners = new Set<ErrorListener>();

export function onApiError(listener: ErrorListener): () => void {
  errorListeners.add(listener);
  return () => errorListeners.delete(listener);
}

function notify(err: ApiError): void {
  for (const l of errorListeners) l(err);
}

// ---------- 真实模式 ----------

async function request<T>(path: string, init: RequestInit & RequestOptions = {}): Promise<T> {
  const { silent, ...rest } = init;
  let res: Response;
  try {
    res = await fetch(`${API_BASE}${path}`, {
      ...rest,
      headers: {
        ...(rest.body && !(rest.body instanceof FormData) ? { 'Content-Type': 'application/json' } : {}),
        ...(rest.headers ?? {}),
      },
    });
  } catch (e) {
    const err = new ApiError('无法连接到服务端（请确认 apps/server 已在 5178 端口启动）', 'network', 0);
    if (!silent) notify(err);
    throw err;
  }
  if (!res.ok) {
    let body: Partial<ApiErrorBody> | null = null;
    try {
      body = (await res.json()) as ApiErrorBody;
    } catch {
      body = null;
    }
    const err = new ApiError(body?.error?.message ?? `${res.status} ${res.statusText}`, body?.error?.code ?? String(res.status), res.status);
    if (!silent) notify(err);
    throw err;
  }
  if (res.status === 204) return undefined as T;
  const text = await res.text();
  if (!text) return undefined as T;
  return JSON.parse(text) as T;
}

function normalizePreview(p: Preview): Preview {
  if (p.kind === 'slides' && !p.slides && typeof p.content === 'string') {
    try {
      return { ...p, slides: JSON.parse(p.content) };
    } catch {
      return { ...p, kind: 'text' };
    }
  }
  return p;
}

function subscribeEventsReal(id: string, since: number, handlers: EventHandlers): () => void {
  let closed = false;
  let lastSeq = since;
  let attempt = 0;
  let es: EventSource | null = null;
  let timer: ReturnType<typeof setTimeout> | null = null;

  const connect = () => {
    if (closed) return;
    handlers.onStatus?.(attempt === 0 ? 'connecting' : 'reconnecting', { attempt });
    es = new EventSource(`${API_BASE}/projects/${encodeURIComponent(id)}/events?since=${lastSeq}`);
    es.addEventListener('agent', (e: MessageEvent<string>) => {
      try {
        const ev = JSON.parse(e.data) as AgentEvent;
        if (typeof ev.seq === 'number') lastSeq = Math.max(lastSeq, ev.seq);
        attempt = 0;
        handlers.onEvent(ev);
      } catch {
        // 忽略无法解析的数据
      }
    });
    es.onopen = () => handlers.onStatus?.('open');
    es.onerror = () => {
      es?.close();
      es = null;
      if (closed) return;
      const delay = Math.min(15_000, 1000 * 2 ** Math.min(attempt, 4));
      attempt += 1;
      handlers.onStatus?.('reconnecting', { retryInMs: delay, attempt });
      timer = setTimeout(connect, delay);
    };
  };
  connect();
  return () => {
    closed = true;
    if (timer) clearTimeout(timer);
    es?.close();
    handlers.onStatus?.('closed');
  };
}

const realApi: BizgoldApi = {
  health: (opts) => request<Health>('/health', { silent: true, ...opts }),
  getSettings: () => request<Settings>('/settings'),
  putSettings: (settings) => request<Settings>('/settings', { method: 'PUT', body: JSON.stringify(settings) }),
  testSettings: (body) => request<SettingsTestResult>('/settings/test', { method: 'POST', body: JSON.stringify(body) }),
  getTracks: () => request<Track[]>('/tracks'),
  getKnowledge: () => request<KnowledgeDoc[]>('/knowledge'),
  getKnowledgeDoc: (slug) => request<KnowledgeContent>(`/knowledge/${encodeURIComponent(slug)}`),
  getProjects: () => request<Project[]>('/projects'),
  createProject: (input) => request<Project>('/projects', { method: 'POST', body: JSON.stringify(input) }),
  getProject: (id, opts) => request<Project>(`/projects/${encodeURIComponent(id)}`, opts),
  deleteProject: (id) => request<void>(`/projects/${encodeURIComponent(id)}?confirm=1`, { method: 'DELETE' }),
  patchProject: (id, patch) => request<Project>(`/projects/${encodeURIComponent(id)}`, { method: 'PATCH', body: JSON.stringify(patch) }),
  uploadFiles: (id, files) => {
    const fd = new FormData();
    for (const f of files) fd.append('files', f, f.name);
    return request<ProjectInput[]>(`/projects/${encodeURIComponent(id)}/files`, { method: 'POST', body: fd });
  },
  deleteFile: (id, name) =>
    request<ProjectInput[]>(`/projects/${encodeURIComponent(id)}/files?name=${encodeURIComponent(name)}`, { method: 'DELETE' }),
  run: (id, body) => request<{ ok: boolean }>(`/projects/${encodeURIComponent(id)}/run`, { method: 'POST', body: JSON.stringify(body) }),
  stop: (id) => request<{ ok: boolean }>(`/projects/${encodeURIComponent(id)}/stop`, { method: 'POST' }),
  chat: (id, message) => request<{ ok: boolean }>(`/projects/${encodeURIComponent(id)}/chat`, { method: 'POST', body: JSON.stringify({ message }) }),
  getArtifacts: async (id) => {
    const tree = await request<ArtifactNode | ArtifactNode[]>(`/projects/${encodeURIComponent(id)}/artifacts`);
    if (Array.isArray(tree)) return { name: '工作区', path: '', type: 'dir', children: tree };
    return tree;
  },
  getPreview: async (id, path) =>
    normalizePreview(await request<Preview>(`/projects/${encodeURIComponent(id)}/preview?path=${encodeURIComponent(path)}`)),
  fileUrl: (id, path) => `${API_BASE}/projects/${encodeURIComponent(id)}/file?path=${encodeURIComponent(path)}`,
  getReview: (id) => request<ReviewData>(`/projects/${encodeURIComponent(id)}/review`),
  getNumbers: (id) => request<NumbersData>(`/projects/${encodeURIComponent(id)}/numbers`),
  repackage: (id) => request<{ path: string }>(`/projects/${encodeURIComponent(id)}/package`, { method: 'POST' }),
  openFolder: (id) => request<{ ok: boolean }>(`/projects/${encodeURIComponent(id)}/open-folder`, { method: 'POST' }),
  subscribeEvents: subscribeEventsReal,
};

// ---------- Mock 模式 ----------

function withErrorNotify(impl: BizgoldApi): BizgoldApi {
  const wrap = <A extends unknown[], R>(fn: (...args: A) => Promise<R>, silentIf?: (args: A) => boolean) =>
    async (...args: A): Promise<R> => {
      try {
        return await fn(...args);
      } catch (e) {
        const err = toApiError(e);
        if (!silentIf?.(args)) notify(err);
        throw err;
      }
    };
  return {
    ...impl,
    health: wrap(impl.health, () => true),
    getSettings: wrap(impl.getSettings),
    putSettings: wrap(impl.putSettings),
    testSettings: wrap(impl.testSettings),
    getTracks: wrap(impl.getTracks),
    getKnowledge: wrap(impl.getKnowledge),
    getKnowledgeDoc: wrap(impl.getKnowledgeDoc),
    getProjects: wrap(impl.getProjects),
    createProject: wrap(impl.createProject),
    getProject: wrap(impl.getProject, (args) => Boolean(args[1]?.silent)),
    deleteProject: wrap(impl.deleteProject),
    patchProject: wrap(impl.patchProject),
    uploadFiles: wrap(impl.uploadFiles),
    deleteFile: wrap(impl.deleteFile),
    run: wrap(impl.run),
    stop: wrap(impl.stop),
    chat: wrap(impl.chat),
    getArtifacts: wrap(impl.getArtifacts),
    getPreview: wrap(impl.getPreview),
    getReview: wrap(impl.getReview),
    getNumbers: wrap(impl.getNumbers),
    repackage: wrap(impl.repackage),
    openFolder: wrap(impl.openFolder),
  };
}

export const api: BizgoldApi = IS_MOCK ? withErrorNotify(createMockApi()) : realApi;
