import path from 'node:path';
import Fastify, { type FastifyInstance } from 'fastify';
import cors from '@fastify/cors';
import multipart from '@fastify/multipart';
import fastifyStatic from '@fastify/static';
import type { AppContext } from './context.js';
import { WEB_DIST_DIR } from './config.js';
import { ApiError } from './types.js';
import { existsSync } from './util/fs.js';
import { createLogger, errorMessage, redact } from './util/log.js';
import { healthRoutes } from './routes/health.js';
import { settingsRoutes } from './routes/settings.js';
import { trackRoutes } from './routes/tracks.js';
import { projectRoutes } from './routes/projects.js';
import { runRoutes } from './routes/run.js';
import { artifactRoutes } from './routes/artifacts.js';

const log = createLogger('http');

/** 开发期前端（Vite）与打包后同源访问；只放行本机来源 */
const DEV_ORIGINS = new Set(['http://localhost:5173', 'http://127.0.0.1:5173']);

export interface BuildAppOptions {
  /** 测试时可关闭静态托管 */
  serveWeb?: boolean;
  /** 单个上传文件上限（默认 200MB） */
  maxUploadBytes?: number;
}

export async function buildApp(ctx: AppContext, opts: BuildAppOptions = {}): Promise<FastifyInstance> {
  const app = Fastify({
    // 自带 logger 会把请求体/headers 打出来，这里用自定义精简日志以免泄露密钥
    logger: false,
    bodyLimit: 10 * 1024 * 1024,
  });

  await app.register(cors, {
    origin: (origin, cb) => {
      // 同源请求（无 Origin）与本机开发端口放行；其余拒绝但不报错（浏览器侧会拦）
      if (!origin || DEV_ORIGINS.has(origin) || /^https?:\/\/(localhost|127\.0\.0\.1)(:\d+)?$/.test(origin)) cb(null, true);
      else cb(null, false);
    },
    methods: ['GET', 'POST', 'PUT', 'PATCH', 'DELETE', 'OPTIONS'],
  });

  await app.register(multipart, {
    limits: { fileSize: opts.maxUploadBytes ?? 200 * 1024 * 1024, files: 20, fields: 20 },
  });

  // 统一错误信封 { error: { code, message } }（契约 §3）
  app.setErrorHandler((err: unknown, req, reply) => {
    if (err instanceof ApiError) {
      void reply.code(err.statusCode).send({ error: { code: err.code, message: err.message } });
      return;
    }
    const e = err as { statusCode?: number; code?: string; message?: string; validation?: unknown };
    const status = typeof e.statusCode === 'number' && e.statusCode >= 400 ? e.statusCode : 500;
    const code = status === 500 ? 'internal' : (e.code ?? 'error');
    const message = redact(e.message ?? String(err));
    if (status >= 500) log.error(`${req.method} ${req.url} → ${status}: ${message}`);
    void reply.code(status).send({ error: { code: String(code).toLowerCase(), message } });
  });

  const serveWeb = (opts.serveWeb ?? true) && existsSync(path.join(WEB_DIST_DIR, 'index.html'));

  // Fastify 同一封装层只能设置一次 notFoundHandler：这里按是否托管前端决定 SPA fallback
  app.setNotFoundHandler((req, reply) => {
    if (req.url.startsWith('/api/')) {
      void reply.code(404).send({ error: { code: 'not_found', message: `未知接口：${req.method} ${req.url}` } });
      return;
    }
    if (serveWeb && (req.method === 'GET' || req.method === 'HEAD')) {
      // SPA fallback：非 /api 的未知路径都回 index.html，由前端路由接管
      void reply.type('text/html; charset=utf-8').sendFile('index.html');
      return;
    }
    void reply.code(404).send({ error: { code: 'not_found', message: 'Not Found' } });
  });

  app.addHook('onResponse', (req, reply, done) => {
    // SSE 连接会长期挂着，响应完成才记一条；只记方法/路径/状态，不记 body
    if (reply.statusCode >= 400 || process.env.BIZGOLD_LOG_LEVEL === 'debug') {
      log.info(`${req.method} ${req.url} ${reply.statusCode} ${reply.elapsedTime.toFixed(0)}ms`);
    }
    done();
  });

  await app.register(
    async (api) => {
      await healthRoutes(api, ctx);
      await settingsRoutes(api, ctx);
      await trackRoutes(api, ctx);
      await projectRoutes(api, ctx);
      await runRoutes(api, ctx);
      await artifactRoutes(api, ctx);
    },
    { prefix: '/api' },
  );

  if (serveWeb) {
    await app.register(fastifyStatic, {
      root: WEB_DIST_DIR,
      prefix: '/',
      wildcard: false,
      index: ['index.html'],
    });
    log.info(`serving web from ${WEB_DIST_DIR}`);
  } else if (opts.serveWeb ?? true) {
    log.info('web dist not found, API only (run `pnpm --filter @bizgold/web build` to serve the UI)');
  }

  app.addHook('onClose', async () => {
    try {
      await ctx.runner.shutdown();
    } catch (e) {
      log.warn('runner shutdown error', errorMessage(e));
    }
  });

  return app;
}
