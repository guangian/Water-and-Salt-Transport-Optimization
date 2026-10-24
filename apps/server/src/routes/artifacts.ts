import path from 'node:path';
import fsp from 'node:fs/promises';
import { createReadStream, createWriteStream } from 'node:fs';
import { spawn } from 'node:child_process';
import { ZipArchive } from 'archiver';
import type { FastifyInstance } from 'fastify';
import type { AppContext } from '../context.js';
import { ApiError } from '../types.js';
import { buildTree, ensureDir, isInside, safeResolve, toPosixRelative } from '../util/fs.js';
import { mimeFor, previewFile } from '../services/preview.js';
import { readNumbers, readReview, scanArtifacts } from '../services/artifacts.js';
import { createLogger, errorMessage } from '../util/log.js';

const log = createLogger('artifacts');

/** 契约 §3：文件树排除运行期目录；.tmp 原子写残留由 buildTree 自己过滤 */
const TREE_EXCLUDES = new Set(['.claude-home', 'toolkit', 'node_modules', '.cache', '.git', '__pycache__']);

export async function artifactRoutes(app: FastifyInstance, ctx: AppContext): Promise<void> {
  app.get<{ Params: { id: string } }>('/projects/:id/artifacts', async (req) => {
    await ctx.store.get(req.params.id);
    const ws = ctx.store.dirOf(req.params.id);
    return buildTree(ws, ws, TREE_EXCLUDES);
  });

  /**
   * 原样下载/内联。路径安全：safeResolve 会把客户端路径解析到工作区内，越界抛错 → 400。
   * 这里不跟随符号链接到工作区外（realpath 后再校验一次）。
   */
  app.get<{ Params: { id: string }; Querystring: { path?: string; download?: string } }>('/projects/:id/file', async (req, reply) => {
    await ctx.store.get(req.params.id);
    const ws = ctx.store.dirOf(req.params.id);
    const abs = resolveWorkspacePath(ws, req.query.path);
    const real = await fsp.realpath(abs).catch(() => null);
    if (!real || !isInside(await fsp.realpath(ws), real)) throw new ApiError(404, 'not_found', '文件不存在');
    const st = await fsp.stat(real);
    if (!st.isFile()) throw new ApiError(400, 'not_a_file', '路径不是文件');
    reply.header('Content-Type', mimeFor(real));
    reply.header('Content-Length', String(st.size));
    reply.header('Cache-Control', 'no-cache');
    const name = path.basename(real);
    const disposition = req.query.download === '1' ? 'attachment' : 'inline';
    // RFC 5987：中文文件名用 filename*，同时给 ASCII 回退
    reply.header('Content-Disposition', `${disposition}; filename="${asciiFallback(name)}"; filename*=UTF-8''${encodeURIComponent(name)}`);
    return reply.send(createReadStream(real));
  });

  app.get<{ Params: { id: string }; Querystring: { path?: string } }>('/projects/:id/preview', async (req) => {
    await ctx.store.get(req.params.id);
    const ws = ctx.store.dirOf(req.params.id);
    if (!req.query.path) throw new ApiError(400, 'bad_request', '缺少 path');
    try {
      return await previewFile(ws, req.params.id, req.query.path);
    } catch (e) {
      if (e instanceof ApiError) throw e;
      throw new ApiError(400, 'bad_path', errorMessage(e));
    }
  });

  app.get<{ Params: { id: string } }>('/projects/:id/review', async (req) => {
    await ctx.store.get(req.params.id);
    return readReview(ctx.store.dirOf(req.params.id));
  });

  app.get<{ Params: { id: string } }>('/projects/:id/numbers', async (req) => {
    await ctx.store.get(req.params.id);
    return readNumbers(ctx.store.dirOf(req.params.id));
  });

  /**
   * 重新打包：优先 python -m bizkit.package（带命名规范/缺项检查），
   * Python 不可用时退化为 archiver 直接压缩 产出/ 目录，保证按钮永远有产物。
   */
  app.post<{ Params: { id: string } }>('/projects/:id/package', async (req) => {
    const project = await ctx.store.get(req.params.id);
    if (ctx.runner.isRunning(project.id)) throw new ApiError(409, 'running', '运行中不能打包');
    const ws = ctx.store.dirOf(project.id);
    let zipRel: string | null = null;
    let warning: string | undefined;
    try {
      const r = await ctx.python.module('package', ['--track', project.track, '--workspace', ws], { workspace: ws, cwd: ws, timeoutMs: 10 * 60 * 1000 });
      const zip = r.json?.ok === true && typeof r.json.zip === 'string' ? r.json.zip : null;
      if (zip) {
        const abs = path.isAbsolute(zip) ? zip : path.resolve(ws, zip);
        zipRel = toPosixRelative(ws, abs);
      } else {
        warning = typeof r.json?.error === 'string' ? r.json.error : r.stderr.trim().split(/\r?\n/).slice(-3).join(' ') || `python exited ${r.code}`;
      }
    } catch (e) {
      warning = errorMessage(e);
    }
    if (!zipRel) {
      log.warn('python package failed, falling back to archiver', warning);
      zipRel = await zipFallback(ws, project.name);
    }
    const artifacts = await scanArtifacts(ws);
    artifacts.package = zipRel;
    await ctx.store.mutate(project.id, (p) => {
      p.artifacts = artifacts;
    });
    return { path: zipRel, ...(warning ? { warning } : {}) };
  });

  /** 本地便利功能：在系统文件管理器打开工作区。失败不抛错（无桌面环境很常见）。 */
  app.post<{ Params: { id: string } }>('/projects/:id/open-folder', async (req) => {
    await ctx.store.get(req.params.id);
    const ws = ctx.store.dirOf(req.params.id);
    const opened = openFolder(ws);
    return { ok: opened, path: ws };
  });
}

function resolveWorkspacePath(ws: string, rel: string | undefined): string {
  if (!rel) throw new ApiError(400, 'bad_request', '缺少 path');
  try {
    return safeResolve(ws, rel);
  } catch (e) {
    throw new ApiError(400, 'bad_path', errorMessage(e));
  }
}

function asciiFallback(name: string): string {
  const ext = path.extname(name);
  // eslint-disable-next-line no-control-regex
  const base = name.slice(0, name.length - ext.length).replace(/[^\x20-\x7e]/g, '').replace(/["\\]/g, '_').trim();
  return (base || 'file') + ext.replace(/[^\x20-\x7e]/g, '');
}

/** 退化打包：把 产出/ 下所有文件压成 08_提交包/<项目名>.zip */
async function zipFallback(ws: string, projectName: string): Promise<string> {
  const outDir = path.join(ws, '08_提交包');
  await ensureDir(outDir);
  const safeName = projectName.replace(/[\\/:*?"<>|\u0000-\u001f]/g, '_').trim() || '提交包';
  const zipAbs = path.join(outDir, `${safeName}.zip`);
  const srcDir = path.join(ws, '产出');
  const hasOutput = await fsp.stat(srcDir).then((s) => s.isDirectory()).catch(() => false);
  await new Promise<void>((resolve, reject) => {
    const output = createWriteStream(zipAbs);
    const archive = new ZipArchive({ zlib: { level: 6 } });
    output.on('close', () => resolve());
    output.on('error', reject);
    archive.on('error', reject);
    archive.pipe(output);
    if (hasOutput) archive.directory(srcDir, false);
    archive.file(path.join(ws, 'project.json'), { name: 'project.json' });
    void archive.finalize();
  });
  return toPosixRelative(ws, zipAbs);
}

function openFolder(dir: string): boolean {
  const cmd = process.platform === 'win32' ? 'explorer.exe' : process.platform === 'darwin' ? 'open' : 'xdg-open';
  try {
    const child = spawn(cmd, [dir], { detached: true, stdio: 'ignore', windowsHide: false });
    child.on('error', (e) => log.debug('open-folder failed', errorMessage(e)));
    child.unref();
    return true;
  } catch (e) {
    log.debug('open-folder failed', errorMessage(e));
    return false;
  }
}
