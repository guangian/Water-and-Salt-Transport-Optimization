import path from 'node:path';
import fsp from 'node:fs/promises';
import { pipeline } from 'node:stream/promises';
import { createWriteStream } from 'node:fs';
import type { FastifyInstance } from 'fastify';
import type { AppContext } from '../context.js';
import type { Project, ProjectCreateInput, ProjectInput } from '../types.js';
import { ApiError } from '../types.js';
import { buildProject } from '../services/store.js';
import { ensureDir, safeResolve, sanitizeFileName, toPosixRelative } from '../util/fs.js';
import { createLogger, errorMessage } from '../util/log.js';

const log = createLogger('projects');

const KIND_BY_EXT: Record<string, string> = {
  '.pdf': 'pdf',
  '.docx': 'docx',
  '.doc': 'doc',
  '.xlsx': 'xlsx',
  '.xls': 'xlsx',
  '.csv': 'csv',
  '.txt': 'text',
  '.md': 'text',
  '.json': 'json',
  '.pptx': 'pptx',
  '.png': 'image',
  '.jpg': 'image',
  '.jpeg': 'image',
  '.gif': 'image',
  '.webp': 'image',
  '.bmp': 'image',
};
const PREREAD_KINDS = new Set(['pdf', 'docx', 'xlsx', 'csv', 'text', 'json']);

export async function projectRoutes(app: FastifyInstance, ctx: AppContext): Promise<void> {
  // 前端首页卡片按完整 Project 契约渲染（review.latestScore / params.targetScore / stages[*].status），这里返回完整对象；摘要版见 ?summary=1
  app.get<{ Querystring: { summary?: string } }>('/projects', async (req) => (req.query.summary ? ctx.store.list() : ctx.store.listFull()));

  app.post<{ Body: ProjectCreateInput }>('/projects', async (req, reply) => {
    const settings = await ctx.settings.get();
    const project = buildProject(req.body ?? ({} as ProjectCreateInput), {
      reviewRounds: settings.run.reviewRounds,
      targetScore: settings.run.targetScore,
      models: settings.api.models,
    });
    await ctx.store.create(project);
    reply.code(201);
    return project;
  });

  app.get<{ Params: { id: string } }>('/projects/:id', async (req) => {
    const p = await ctx.store.get(req.params.id);
    return { ...p, running: ctx.runner.isRunning(p.id) };
  });

  app.delete<{ Params: { id: string }; Querystring: { confirm?: string } }>('/projects/:id', async (req) => {
    if (req.query.confirm !== '1') throw new ApiError(400, 'confirm_required', '删除需要 ?confirm=1');
    const id = req.params.id;
    if (ctx.runner.isRunning(id)) {
      await ctx.runner.stop(id);
      await new Promise((r) => setTimeout(r, 500));
    }
    if (!(await ctx.store.exists(id))) throw new ApiError(404, 'not_found', '项目不存在');
    await ctx.store.remove(id);
    return { ok: true };
  });

  /** 修改元数据/参数：运行中仅允许改 name */
  app.patch<{ Params: { id: string }; Body: Partial<Project> }>('/projects/:id', async (req) => {
    const body = req.body ?? {};
    const running = ctx.runner.isRunning(req.params.id);
    return ctx.store.mutate(req.params.id, (p) => {
      if (typeof body.name === 'string' && body.name.trim()) p.name = body.name.trim();
      if (running) {
        const others = Object.keys(body).filter((k) => k !== 'name');
        if (others.length) throw new ApiError(409, 'running', `运行中仅允许修改 name（收到：${others.join(', ')}）`);
        return;
      }
      if (typeof body.teamName === 'string') p.teamName = body.teamName.trim();
      if (typeof body.group === 'string') p.group = body.group;
      if (body.subject && typeof body.subject === 'object') p.subject = { ...p.subject, ...body.subject };
      if (Array.isArray(body.team)) p.team = body.team.map((m) => ({ role: String(m?.role ?? ''), major: String(m?.major ?? ''), strength: String(m?.strength ?? '') }));
      if (body.params && typeof body.params === 'object') p.params = { ...p.params, ...body.params, models: { ...p.params.models, ...(body.params.models ?? {}) } };
    });
  });

  /** multipart 上传到 inputs/，自动预读 */
  app.post<{ Params: { id: string } }>('/projects/:id/files', async (req) => {
    const id = req.params.id;
    const project = await ctx.store.get(id);
    if (ctx.runner.isRunning(id)) throw new ApiError(409, 'running', '运行中不能上传文件');
    const ws = ctx.store.dirOf(id);
    const inputsDir = path.join(ws, 'inputs');
    await ensureDir(inputsDir);
    const saved: string[] = [];
    for await (const part of req.parts()) {
      if (part.type !== 'file') continue;
      const name = sanitizeFileName(part.filename, `upload_${Date.now()}`);
      if (name.startsWith('_')) throw new ApiError(400, 'bad_name', '文件名不能以下划线开头（_preread 为保留目录）');
      const dest = path.join(inputsDir, name);
      await pipeline(part.file, createWriteStream(dest));
      if (part.file.truncated) {
        await fsp.rm(dest, { force: true });
        throw new ApiError(413, 'too_large', `文件过大：${name}`);
      }
      saved.push(name);
    }
    if (saved.length === 0) throw new ApiError(400, 'no_file', '没有收到文件');

    const inputs: ProjectInput[] = [...project.inputs.filter((i) => !saved.includes(i.name))];
    for (const name of saved) {
      const abs = path.join(inputsDir, name);
      const st = await fsp.stat(abs);
      const kind = KIND_BY_EXT[path.extname(name).toLowerCase()] ?? 'other';
      const input: ProjectInput = { name, path: `inputs/${name}`, size: st.size, kind };
      if (PREREAD_KINDS.has(kind)) {
        const preread = await runPreread(ctx, ws, abs);
        if (preread) input.preread = preread;
      }
      inputs.push(input);
    }
    const updated = await ctx.store.mutate(id, (p) => {
      p.inputs = inputs;
    });
    return updated.inputs;
  });

  app.delete<{ Params: { id: string }; Querystring: { name?: string } }>('/projects/:id/files', async (req) => {
    const id = req.params.id;
    const name = sanitizeFileName(req.query.name ?? '');
    if (!req.query.name) throw new ApiError(400, 'bad_request', '缺少 name');
    if (ctx.runner.isRunning(id)) throw new ApiError(409, 'running', '运行中不能删除文件');
    const ws = ctx.store.dirOf(id);
    const abs = safeResolve(ws, `inputs/${name}`);
    await fsp.rm(abs, { force: true });
    await fsp.rm(safeResolve(ws, `inputs/_preread/${name}.txt`), { force: true });
    const updated = await ctx.store.mutate(id, (p) => {
      p.inputs = p.inputs.filter((i) => i.name !== name);
    });
    return updated.inputs;
  });
}

async function runPreread(ctx: AppContext, ws: string, absFile: string): Promise<string | undefined> {
  try {
    const r = await ctx.python.module('preread', [absFile], { workspace: ws, cwd: ws, timeoutMs: 180000 });
    const out = typeof r.json?.out === 'string' ? r.json.out : path.join(path.dirname(absFile), '_preread', `${path.basename(absFile)}.txt`);
    const abs = path.isAbsolute(out) ? out : path.resolve(ws, out);
    await fsp.access(abs);
    return toPosixRelative(ws, abs);
  } catch (e) {
    log.warn(`preread failed for ${path.basename(absFile)}`, errorMessage(e));
    return undefined;
  }
}
