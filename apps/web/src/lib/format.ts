export function formatDuration(ms: number | null | undefined, opts: { compact?: boolean } = {}): string {
  if (ms == null || !Number.isFinite(ms)) return '—';
  const total = Math.max(0, Math.floor(ms / 1000));
  const h = Math.floor(total / 3600);
  const m = Math.floor((total % 3600) / 60);
  const s = total % 60;
  if (opts.compact) {
    if (h > 0) return `${h}h${m}m`;
    if (m > 0) return `${m}m${s}s`;
    return `${s}s`;
  }
  if (h > 0) return `${h} 小时 ${m} 分`;
  if (m > 0) return `${m} 分 ${s} 秒`;
  return `${s} 秒`;
}

export function formatBytes(n: number | undefined | null): string {
  if (n == null || !Number.isFinite(n)) return '—';
  if (n < 1024) return `${n} B`;
  if (n < 1024 * 1024) return `${(n / 1024).toFixed(1)} KB`;
  if (n < 1024 * 1024 * 1024) return `${(n / 1024 / 1024).toFixed(1)} MB`;
  return `${(n / 1024 / 1024 / 1024).toFixed(2)} GB`;
}

export function formatTokens(n: number | undefined | null): string {
  if (n == null) return '0';
  if (n < 1000) return String(n);
  if (n < 1_000_000) return `${(n / 1000).toFixed(1)}k`;
  return `${(n / 1_000_000).toFixed(2)}M`;
}

export function formatUsd(n: number | undefined | null): string {
  if (n == null) return '$0.00';
  return `$${n.toFixed(n < 10 ? 2 : 1)}`;
}

export function formatRelative(iso: string | undefined | null, now: number = Date.now()): string {
  if (!iso) return '—';
  const t = Date.parse(iso);
  if (Number.isNaN(t)) return '—';
  const diff = now - t;
  const min = Math.round(diff / 60000);
  if (min < 1) return '刚刚';
  if (min < 60) return `${min} 分钟前`;
  const h = Math.round(min / 60);
  if (h < 24) return `${h} 小时前`;
  const d = Math.round(h / 24);
  if (d < 30) return `${d} 天前`;
  return formatDate(iso);
}

export function formatDate(iso: string | undefined | null): string {
  if (!iso) return '—';
  const d = new Date(iso);
  if (Number.isNaN(d.getTime())) return '—';
  const pad = (x: number) => String(x).padStart(2, '0');
  return `${d.getFullYear()}-${pad(d.getMonth() + 1)}-${pad(d.getDate())}`;
}

export function formatDateTime(iso: string | undefined | null): string {
  if (!iso) return '—';
  const d = new Date(iso);
  if (Number.isNaN(d.getTime())) return '—';
  const pad = (x: number) => String(x).padStart(2, '0');
  return `${d.getFullYear()}-${pad(d.getMonth() + 1)}-${pad(d.getDate())} ${pad(d.getHours())}:${pad(d.getMinutes())}`;
}

export function formatTime(iso: string | undefined | null): string {
  if (!iso) return '';
  const d = new Date(iso);
  if (Number.isNaN(d.getTime())) return '';
  const pad = (x: number) => String(x).padStart(2, '0');
  return `${pad(d.getHours())}:${pad(d.getMinutes())}:${pad(d.getSeconds())}`;
}

export function formatNumber(v: number | string | boolean | null | undefined, unit?: string): string {
  if (v == null) return '—';
  if (typeof v === 'boolean') return v ? '是' : '否';
  if (typeof v === 'string') return v;
  if (unit === '%') return `${(v * 100).toFixed(v * 100 >= 10 ? 1 : 2).replace(/\.?0+$/, '')}%`;
  if (Number.isInteger(v)) return v.toLocaleString('zh-CN');
  return v.toLocaleString('zh-CN', { maximumFractionDigits: 2 });
}

export function formatScore(v: number | null | undefined): string {
  if (v == null || !Number.isFinite(v)) return '—';
  return v.toFixed(1);
}

export function truncate(s: string, n: number): string {
  if (s.length <= n) return s;
  return `${s.slice(0, n)}…`;
}

export function summarizeInput(input: unknown, max = 80): string {
  if (input == null) return '';
  if (typeof input === 'string') return truncate(input, max);
  if (typeof input !== 'object') return truncate(String(input), max);
  const obj = input as Record<string, unknown>;
  const keys = ['file_path', 'path', 'command', 'query', 'url', 'pattern', 'description', 'prompt', 'stage', 'message'];
  for (const k of keys) {
    const v = obj[k];
    if (typeof v === 'string' && v.trim()) return truncate(v, max);
  }
  try {
    return truncate(JSON.stringify(input), max);
  } catch {
    return '';
  }
}

export function pretty(v: unknown): string {
  if (v == null) return '';
  if (typeof v === 'string') return v;
  try {
    return JSON.stringify(v, null, 2);
  } catch {
    return String(v);
  }
}

export function fileExt(path: string): string {
  const base = path.split('/').pop() ?? path;
  const i = base.lastIndexOf('.');
  return i >= 0 ? base.slice(i + 1).toLowerCase() : '';
}

export function baseName(path: string): string {
  return path.split('/').pop() ?? path;
}
