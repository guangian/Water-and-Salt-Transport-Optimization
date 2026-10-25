import type { FastifyInstance } from 'fastify';
import type { AppContext } from '../context.js';
import { ApiError } from '../types.js';
import { isMaskedKey } from '../services/settings.js';
import { testGatewayWithFallback } from '../services/gateway.js';

export async function settingsRoutes(app: FastifyInstance, ctx: AppContext): Promise<void> {
  app.get('/settings', async () => ctx.settings.getPublic());

  app.put<{ Body: unknown }>('/settings', async (req) => {
    if (!req.body || typeof req.body !== 'object') throw new ApiError(400, 'bad_request', '请求体必须是 JSON 对象');
    const result = await ctx.settings.update(req.body);
    ctx.python.invalidate();
    return result;
  });

  /** 连通性测试：GET <baseUrl>/v1/models；body 可覆盖 baseUrl/apiKey（脱敏串或缺省用当前生效值） */
  app.post<{ Body: { baseUrl?: string; apiKey?: string } }>('/settings/test', async (req) => {
    const body = req.body ?? {};
    const current = await ctx.settings.get();
    const key = body.apiKey && !isMaskedKey(body.apiKey) ? body.apiKey : (await ctx.settings.apiKey()).key;
    if (!key) return { ok: false, latencyMs: 0, models: [], error: '未配置 API 密钥' };
    const base = body.baseUrl?.trim() || current.api.baseUrl;
    const r = await testGatewayWithFallback(base, current.api.fallbackBaseUrl, key);
    return { ok: r.ok, latencyMs: r.latencyMs, models: r.models, baseUrl: r.baseUrl, usedFallback: r.usedFallback ?? false, ...(r.error ? { error: r.error } : {}) };
  });
}
