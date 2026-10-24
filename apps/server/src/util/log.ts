/**
 * 日志与脱敏。密钥只存在设置文件 / 环境变量里，任何写到日志、events.jsonl 的文本都先经过 redact()。
 */

const registered = new Set<string>();

/** 登记需要从所有输出中抹掉的敏感值（API key、搜索 key） */
export function registerSecret(value: string | undefined | null): void {
  if (value && value.length >= 8) registered.add(value);
}

const PATTERNS: RegExp[] = [
  /sk-[A-Za-z0-9_-]{8,}/g,
  /(Authorization\s*[:=]\s*)(Bearer\s+)?[A-Za-z0-9._~+/=-]{8,}/gi,
  /(ANTHROPIC_AUTH_TOKEN|ANTHROPIC_API_KEY|BIZGOLD_API_KEY|BOCHA_API_KEY|TAVILY_API_KEY|SERPER_API_KEY)(\s*[=:]\s*)[^\s"']+/g,
];

export function redact(text: string): string {
  if (!text) return text;
  let out = text;
  for (const secret of registered) {
    if (out.includes(secret)) out = out.split(secret).join('****');
  }
  out = out.replace(PATTERNS[0]!, 'sk-****');
  out = out.replace(PATTERNS[1]!, '$1$2****');
  out = out.replace(PATTERNS[2]!, '$1$2****');
  return out;
}

export function redactDeep<T>(value: T): T {
  if (typeof value === 'string') return redact(value) as T;
  if (Array.isArray(value)) return value.map((v) => redactDeep(v)) as T;
  if (value && typeof value === 'object') {
    const out: Record<string, unknown> = {};
    for (const [k, v] of Object.entries(value as Record<string, unknown>)) out[k] = redactDeep(v);
    return out as T;
  }
  return value;
}

type Level = 'debug' | 'info' | 'warn' | 'error';
const LEVELS: Record<Level, number> = { debug: 10, info: 20, warn: 30, error: 40 };
const threshold = LEVELS[(process.env.BIZGOLD_LOG_LEVEL as Level) || 'info'] ?? 20;

function write(level: Level, scope: string, msg: string, extra?: unknown): void {
  if (LEVELS[level] < threshold) return;
  const ts = new Date().toISOString();
  const line = `${ts} [${level.toUpperCase()}] [${scope}] ${redact(msg)}`;
  const tail = extra === undefined ? '' : ' ' + redact(typeof extra === 'string' ? extra : safeStringify(extra));
  if (level === 'error' || level === 'warn') process.stderr.write(line + tail + '\n');
  else process.stdout.write(line + tail + '\n');
}

function safeStringify(v: unknown): string {
  try {
    return JSON.stringify(v);
  } catch {
    return String(v);
  }
}

export function createLogger(scope: string) {
  return {
    debug: (msg: string, extra?: unknown) => write('debug', scope, msg, extra),
    info: (msg: string, extra?: unknown) => write('info', scope, msg, extra),
    warn: (msg: string, extra?: unknown) => write('warn', scope, msg, extra),
    error: (msg: string, extra?: unknown) => write('error', scope, msg, extra),
  };
}

export type Logger = ReturnType<typeof createLogger>;

export function errorMessage(err: unknown): string {
  if (err instanceof Error) return err.message;
  return String(err);
}
