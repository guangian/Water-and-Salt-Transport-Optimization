import path from 'node:path';
import fsp from 'node:fs/promises';
import { DEFAULT_SETTINGS_DIR, DEFAULT_WORKSPACE_DIR } from '../config.js';
import type { ApiKeySource, PublicSettings, Settings, SearchProviderId } from '../types.js';
import { ensureDir, readJsonOrNull, writeJsonAtomic } from '../util/fs.js';
import { createLogger, registerSecret } from '../util/log.js';

const log = createLogger('settings');

export const DEFAULT_SETTINGS: Settings = {
  api: {
    baseUrl: 'https://llm.goaichat.top',
    apiKey: '',
    fallbackBaseUrl: 'https://model.kimitk.top',
    models: {
      main: 'deepseek-v4-pro',
      researcher: 'kimi-k3',
      finance: 'deepseek-v4-pro',
      writer: 'deepseek-v4-pro',
      pitch: 'deepseek-v4-pro',
      judgeA: 'kimi-k3',
      judgeB: 'glm-5.3',
      judgeC: 'qwen3.8-max',
      inspector: 'deepseek-v4-pro',
      fast: 'deepseek-v4-flash',
    },
    contextTokens: 200000,
    disableThinking: false,
    disableExperimentalBetas: false,
    pricing: {},
  },
  search: {
    providers: ['bocha', 'tavily', 'serper', 'so360', 'bing', 'duckduckgo'],
    bochaKey: '',
    tavilyKey: '',
    serperKey: '',
    maxResults: 8,
  },
  python: { path: '', pipIndex: 'https://pypi.tuna.tsinghua.edu.cn/simple' },
  run: {
    maxTurns: 1500,
    totalTimeoutMin: 360,
    idleTimeoutMin: 45,
    bashTimeoutMs: 600000,
    reviewRounds: 3,
    targetScore: 85,
  },
  paths: { workspaceDir: '' },
  ui: { theme: 'light', language: 'zh-CN' },
};

const MASK_RE = /^(sk-)?\*{4}/;

export function maskApiKey(key: string): string {
  if (!key) return '';
  const prefix = key.startsWith('sk-') ? 'sk-' : '';
  return `${prefix}****${key.slice(-4)}`;
}

export function isMaskedKey(value: unknown): boolean {
  return typeof value === 'string' && MASK_RE.test(value);
}

function isPlainObject(v: unknown): v is Record<string, unknown> {
  return !!v && typeof v === 'object' && !Array.isArray(v);
}

/** 深合并：对象递归，数组与标量整体替换 */
export function deepMerge<T>(base: T, patch: unknown): T {
  if (!isPlainObject(patch)) return (patch === undefined ? base : (patch as T));
  const out: Record<string, unknown> = isPlainObject(base) ? { ...(base as Record<string, unknown>) } : {};
  for (const [k, v] of Object.entries(patch)) {
    const cur = out[k];
    out[k] = isPlainObject(v) && isPlainObject(cur) ? deepMerge(cur, v) : v;
  }
  return out as T;
}

const VALID_PROVIDERS: SearchProviderId[] = ['bocha', 'tavily', 'serper', 'so360', 'bing', 'duckduckgo'];

/**
 * 环境变量优先于设置文件（契约 §9）。返回新对象，不改原值。
 */
export function applyEnvOverrides(s: Settings, env: NodeJS.ProcessEnv = process.env): { settings: Settings; apiKeySource: ApiKeySource } {
  const out: Settings = deepMerge(DEFAULT_SETTINGS, s);
  let apiKeySource: ApiKeySource = out.api.apiKey ? 'settings' : 'none';
  if (env.BIZGOLD_API_BASE_URL) out.api.baseUrl = env.BIZGOLD_API_BASE_URL;
  const envKey = env.BIZGOLD_API_KEY || env.ANTHROPIC_AUTH_TOKEN;
  if (envKey) {
    out.api.apiKey = envKey;
    apiKeySource = 'env';
  }
  if (env.BIZGOLD_MODEL) out.api.models = { ...out.api.models, main: env.BIZGOLD_MODEL };
  if (env.BIZGOLD_WORKSPACE_DIR) out.paths.workspaceDir = env.BIZGOLD_WORKSPACE_DIR;
  if (env.BIZGOLD_PYTHON) out.python.path = env.BIZGOLD_PYTHON;
  if (env.BOCHA_API_KEY) out.search.bochaKey = env.BOCHA_API_KEY;
  if (env.TAVILY_API_KEY) out.search.tavilyKey = env.TAVILY_API_KEY;
  if (env.SERPER_API_KEY) out.search.serperKey = env.SERPER_API_KEY;
  out.search.providers = (out.search.providers ?? []).filter((p): p is SearchProviderId => VALID_PROVIDERS.includes(p));
  if (out.search.providers.length === 0) out.search.providers = [...DEFAULT_SETTINGS.search.providers];
  return { settings: out, apiKeySource };
}

export class SettingsStore {
  private cache: Settings | null = null;
  readonly file: string;

  constructor(readonly dir: string = DEFAULT_SETTINGS_DIR) {
    this.file = path.join(dir, 'settings.json');
  }

  /** 文件原样内容（已与默认值合并，不含环境变量覆盖） */
  async loadRaw(): Promise<Settings> {
    if (this.cache) return this.cache;
    const onDisk = await readJsonOrNull<Partial<Settings>>(this.file);
    this.cache = deepMerge(DEFAULT_SETTINGS, onDisk ?? {});
    return this.cache;
  }

  /** 生效设置：文件 + 环境变量 */
  async get(): Promise<Settings> {
    const { settings } = applyEnvOverrides(await this.loadRaw());
    registerSecret(settings.api.apiKey);
    registerSecret(settings.search.bochaKey);
    registerSecret(settings.search.tavilyKey);
    registerSecret(settings.search.serperKey);
    return settings;
  }

  async apiKey(): Promise<{ key: string; source: ApiKeySource }> {
    const { settings, apiKeySource } = applyEnvOverrides(await this.loadRaw());
    return { key: settings.api.apiKey, source: apiKeySource };
  }

  /** 对外返回：apiKey 脱敏 + apiKeySource；搜索 key 同样脱敏 */
  async getPublic(): Promise<PublicSettings> {
    const { settings, apiKeySource } = applyEnvOverrides(await this.loadRaw());
    return {
      ...settings,
      api: { ...settings.api, apiKey: maskApiKey(settings.api.apiKey), apiKeySource },
      search: {
        ...settings.search,
        bochaKey: maskApiKey(settings.search.bochaKey),
        tavilyKey: maskApiKey(settings.search.tavilyKey),
        serperKey: maskApiKey(settings.search.serperKey),
      },
    };
  }

  /** PUT：合并写入。脱敏串不改动原值。文件权限 0600（Windows 忽略 mode）。 */
  async update(patch: unknown): Promise<PublicSettings> {
    const current = await this.loadRaw();
    const p = isPlainObject(patch) ? structuredClone(patch) : {};
    const api = isPlainObject(p.api) ? p.api : undefined;
    if (api && (isMaskedKey(api.apiKey) || api.apiKey === undefined)) delete api.apiKey;
    if (api && 'apiKeySource' in api) delete api.apiKeySource;
    const search = isPlainObject(p.search) ? p.search : undefined;
    for (const k of ['bochaKey', 'tavilyKey', 'serperKey'] as const) {
      if (search && isMaskedKey(search[k])) delete search[k];
    }
    const merged = deepMerge(current, p);
    await ensureDir(this.dir);
    await writeJsonAtomic(this.file, merged, 0o600);
    try {
      await fsp.chmod(this.file, 0o600);
    } catch {
      /* Windows 无 POSIX 权限 */
    }
    this.cache = merged;
    log.info('settings saved', { file: this.file });
    return this.getPublic();
  }

  /** 工作区目录：环境变量 > 设置 > <repo>/workspace */
  async workspaceDir(): Promise<string> {
    const s = await this.get();
    return s.paths.workspaceDir ? path.resolve(s.paths.workspaceDir) : DEFAULT_WORKSPACE_DIR;
  }

  invalidate(): void {
    this.cache = null;
  }
}
