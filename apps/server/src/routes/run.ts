import type { FastifyInstance } from 'fastify';
import type { AppContext } from '../context.js';
import type { AgentEvent, StageId } from '../types.js';
import { ApiError, STAGE_IDS } from '../types.js';
import { createLogger } from '../util/log.js';

const log = createLogger('sse');

const HEARTBEAT_MS = 15_000;

export async function runRoutes(app: FastifyInstance, ctx: AppContext): Promise<void> {
  app.post<{ Params: { id: string }; Body: { mode?: string; fromStage?: string; smoke?: boolean; maxTurns?: number } }>('/projects/:id/run', async (req) => {
    const body = req.body ?? {};
    const mode = body.mode === 'resume' ? 'resume' : 'run';
    const fromStage = body.fromStage as StageId | undefined;
    if (fromStage && !STAGE_IDS.includes(fromStage)) throw new ApiError(400, 'bad_request', `未知阶段：${fromStage}`);
    await ctx.store.get(req.params.id);
    await ctx.runner.start(req.params.id, {
      mode,
      fromStage,
      smoke: body.smoke === true,
      ...(typeof body.maxTurns === 'number' && body.maxTurns > 0 ? { maxTurns: Math.floor(body.maxTurns) } : {}),
    });
    return { ok: true };
  });

  app.post<{ Params: { id: string } }>('/projects/:id/stop', async (req) => {
    await ctx.store.get(req.params.id);
    const stopped = await ctx.runner.stop(req.params.id);
    return { ok: true, stopped };
  });

  app.post<{ Params: { id: string }; Body: { message?: string } }>('/projects/:id/chat', async (req) => {
    const message = (req.body?.message ?? '').trim();
    if (!message) throw new ApiError(400, 'bad_request', 'message 不能为空');
    await ctx.store.get(req.params.id);
    await ctx.runner.start(req.params.id, { mode: 'modify', message });
    return { ok: true };
  });

  /**
   * SSE：先回放 events.jsonl 中 seq > since 的事件，再实时推送。
   * 先订阅再回放（实时事件进入缓冲），回放完按 seq 过滤掉重复，保证既不丢也不重。
   */
  app.get<{ Params: { id: string }; Querystring: { since?: string } }>('/projects/:id/events', async (req, reply) => {
    const id = req.params.id;
    await ctx.store.get(id);
    const since = Number(req.query.since) || 0;

    reply.raw.writeHead(200, {
      'Content-Type': 'text/event-stream; charset=utf-8',
      'Cache-Control': 'no-cache, no-transform',
      Connection: 'keep-alive',
      'X-Accel-Buffering': 'no',
      ...(corsHeader(req.headers.origin)),
    });
    reply.raw.write(': connected\n\n');

    let closed = false;
    let lastSeq = since;
    let replaying = true;
    const buffer: AgentEvent[] = [];
    const send = (ev: AgentEvent) => {
      if (closed || ev.seq <= lastSeq) return;
      lastSeq = ev.seq;
      reply.raw.write(`id: ${ev.seq}\nevent: agent\ndata: ${JSON.stringify(ev)}\n\n`);
    };
    const unsubscribe = ctx.store.subscribe(id, (ev) => {
      if (replaying) buffer.push(ev);
      else send(ev);
    });
    const heartbeat = setInterval(() => {
      if (!closed) reply.raw.write(': ping\n\n');
    }, HEARTBEAT_MS);
    const cleanup = () => {
      if (closed) return;
      closed = true;
      clearInterval(heartbeat);
      unsubscribe();
      log.debug(`sse closed project=${id} lastSeq=${lastSeq}`);
    };
    req.raw.on('close', cleanup);
    req.raw.on('error', cleanup);

    try {
      const history = await ctx.store.readEvents(id, since);
      for (const ev of history) send(ev);
    } finally {
      replaying = false;
      for (const ev of buffer) send(ev);
      buffer.length = 0;
    }
    if (!ctx.runner.isRunning(id)) {
      // 没在运行也保持连接：前端用同一连接等待后续 run/chat 的事件
      reply.raw.write(`event: idle\ndata: ${JSON.stringify({ projectId: id, lastSeq })}\n\n`);
    }
    return reply;
  });
}

function corsHeader(origin: string | undefined): Record<string, string> {
  if (!origin) return {};
  try {
    const u = new URL(origin);
    if (u.hostname === 'localhost' || u.hostname === '127.0.0.1') return { 'Access-Control-Allow-Origin': origin, Vary: 'Origin' };
  } catch {
    /* ignore */
  }
  return {};
}
