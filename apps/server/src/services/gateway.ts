import { createLogger, errorMessage } from '../util/log.js';

const log = createLogger('gateway');

export interface GatewayTest {
  ok: boolean;
  baseUrl: string;
  latencyMs: number;
  models: string[];
  error?: string;
  usedFallback?: boolean;
}

function normalizeBase(url: string): string {
  return url.trim().replace(/\/+$/, '').replace(/\/v1$/, '');
}

/** GET <baseUrl>/v1/models（Authorization: Bearer）。网关若不实现该端点但返回 404，也视为“可连通但无法列模型”。 */
export async function testGateway(baseUrl: string, apiKey: string, timeoutMs = 10000): Promise<GatewayTest> {
  const base = normalizeBase(baseUrl);
  const t0 = Date.now();
  try {
    const res = await fetch(`${base}/v1/models`, {
      headers: { Authorization: `Bearer ${apiKey}`, 'x-api-key': apiKey, 'anthropic-version': '2023-06-01', Accept: 'application/json' },
      signal: AbortSignal.timeout(timeoutMs),
    });
    const latencyMs = Date.now() - t0;
    const text = await res.text();
    if (!res.ok) {
      if (res.status === 404 || res.status === 405) return { ok: true, baseUrl: base, latencyMs, models: [], error: `网关未提供 /v1/models（HTTP ${res.status}），连通性正常` };
      return { ok: false, baseUrl: base, latencyMs, models: [], error: `HTTP ${res.status}: ${text.slice(0, 200)}` };
    }
    let models: string[] = [];
    try {
      const data = JSON.parse(text) as { data?: unknown[]; models?: unknown[] };
      const arr = Array.isArray(data.data) ? data.data : Array.isArray(data.models) ? data.models : [];
      models = arr
        .map((m) => (typeof m === 'string' ? m : m && typeof m === 'object' ? String((m as { id?: unknown; name?: unknown }).id ?? (m as { name?: unknown }).name ?? '') : ''))
        .filter(Boolean)
        .sort();
    } catch {
      /* 非 JSON：连通但无法解析 */
    }
    return { ok: true, baseUrl: base, latencyMs, models };
  } catch (e) {
    return { ok: false, baseUrl: base, latencyMs: Date.now() - t0, models: [], error: errorMessage(e) };
  }
}

/** 主地址失败时尝试备用地址 */
export async function testGatewayWithFallback(baseUrl: string, fallbackBaseUrl: string, apiKey: string): Promise<GatewayTest> {
  const primary = await testGateway(baseUrl, apiKey);
  if (primary.ok || !fallbackBaseUrl || normalizeBase(fallbackBaseUrl) === normalizeBase(baseUrl)) return primary;
  log.warn(`primary gateway failed (${primary.error}), trying fallback`);
  const fb = await testGateway(fallbackBaseUrl, apiKey);
  if (fb.ok) return { ...fb, usedFallback: true, error: `主地址失败：${primary.error}` };
  return { ...primary, error: `主地址：${primary.error}；备用地址：${fb.error}` };
}

/** 运行前选地址：主地址可连通用主地址，否则用备用（都不通仍返回主地址，让 CLI 自己重试报错） */
export async function pickBaseUrl(baseUrl: string, fallbackBaseUrl: string, apiKey: string): Promise<{ baseUrl: string; usedFallback: boolean; note?: string }> {
  const r = await testGatewayWithFallback(baseUrl, fallbackBaseUrl, apiKey);
  if (r.ok) return { baseUrl: r.baseUrl, usedFallback: Boolean(r.usedFallback), note: r.usedFallback ? r.error : undefined };
  return { baseUrl: normalizeBase(baseUrl), usedFallback: false, note: `网关连通性检查失败：${r.error}` };
}
