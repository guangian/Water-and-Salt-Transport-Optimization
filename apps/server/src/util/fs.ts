import fs from 'node:fs';
import fsp from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { IS_WINDOWS } from '../config.js';

export async function exists(p: string): Promise<boolean> {
  try {
    await fsp.access(p);
    return true;
  } catch {
    return false;
  }
}

export function existsSync(p: string): boolean {
  try {
    fs.accessSync(p);
    return true;
  } catch {
    return false;
  }
}

export async function ensureDir(p: string): Promise<void> {
  await fsp.mkdir(p, { recursive: true });
}

/**
 * 原子写：先写同目录临时文件再 rename。rename 在 POSIX 与 Windows（MoveFileEx 覆盖）上都是替换语义，
 * 避免进程中途被杀时留下半截 project.json。
 */
export async function writeFileAtomic(file: string, data: string | Buffer, mode?: number): Promise<void> {
  await ensureDir(path.dirname(file));
  const tmp = path.join(path.dirname(file), `.${path.basename(file)}.${process.pid}.${crypto.randomBytes(4).toString('hex')}.tmp`);
  try {
    await fsp.writeFile(tmp, data, mode !== undefined ? { mode } : undefined);
    await fsp.rename(tmp, file);
  } catch (err) {
    await fsp.rm(tmp, { force: true }).catch(() => undefined);
    throw err;
  }
}

export async function readJson<T>(file: string): Promise<T> {
  const text = await fsp.readFile(file, 'utf8');
  return JSON.parse(text.replace(/^\uFEFF/, '')) as T;
}

export async function readJsonOrNull<T>(file: string): Promise<T | null> {
  try {
    return await readJson<T>(file);
  } catch {
    return null;
  }
}

export async function writeJsonAtomic(file: string, value: unknown, mode?: number): Promise<void> {
  await writeFileAtomic(file, JSON.stringify(value, null, 2) + '\n', mode);
}

/** 规范化成比较用的形式：Windows 不区分大小写，统一小写比较 */
function comparable(p: string): string {
  const r = path.resolve(p);
  return IS_WINDOWS ? r.toLowerCase() : r;
}

/** target 是否位于 root 内（含 root 本身）。先 path.resolve 再做前缀校验，防止 ../ 与绝对路径穿越。 */
export function isInside(root: string, target: string): boolean {
  const r = comparable(root);
  const t = comparable(target);
  if (t === r) return true;
  const prefix = r.endsWith(path.sep) ? r : r + path.sep;
  return t.startsWith(prefix);
}

/**
 * 把客户端传来的相对路径解析到 root 内；越界抛错。
 * 接受 / 与 \ 两种分隔符（前端统一用 /）。
 */
export function safeResolve(root: string, rel: string): string {
  const cleaned = String(rel ?? '').replace(/\\/g, '/').replace(/^\/+/, '');
  if (cleaned.includes('\0')) throw new Error('非法路径');
  const abs = path.resolve(root, cleaned);
  if (!isInside(root, abs)) throw new Error('路径越界');
  return abs;
}

/** 转成工作区相对 POSIX 路径（事件与 project.json 中统一用 /） */
export function toPosixRelative(root: string, abs: string): string {
  const rel = path.relative(root, abs);
  return rel.split(path.sep).join('/');
}

export interface CopyDirOptions {
  /** 目录或文件名命中即跳过（如 __pycache__、examples、tests） */
  skipNames?: string[];
}

/** 跨平台递归复制（覆盖同名文件）；源不存在返回 false 以便调用方记录警告。 */
export async function copyDir(src: string, dst: string, opts: CopyDirOptions = {}): Promise<boolean> {
  if (!(await exists(src))) return false;
  const skip = new Set(opts.skipNames ?? []);
  const stat = await fsp.stat(src);
  if (stat.isFile()) {
    await ensureDir(path.dirname(dst));
    await fsp.copyFile(src, dst);
    return true;
  }
  await ensureDir(dst);
  for (const entry of await fsp.readdir(src, { withFileTypes: true })) {
    if (skip.has(entry.name)) continue;
    const s = path.join(src, entry.name);
    const d = path.join(dst, entry.name);
    if (entry.isDirectory()) await copyDir(s, d, opts);
    else if (entry.isFile()) {
      await ensureDir(path.dirname(d));
      await fsp.copyFile(s, d);
    }
  }
  return true;
}

export interface TreeNode {
  name: string;
  path: string;
  type: 'dir' | 'file';
  size: number;
  mtime: string;
  children?: TreeNode[];
}

/** 构建文件树（path 为相对 root 的 POSIX 路径），目录在前、按名称排序 */
export async function buildTree(root: string, dir: string, excludeNames: Set<string>, depth = 0): Promise<TreeNode[]> {
  if (depth > 12) return [];
  let entries: fs.Dirent[];
  try {
    entries = await fsp.readdir(dir, { withFileTypes: true });
  } catch {
    return [];
  }
  const nodes: TreeNode[] = [];
  for (const e of entries) {
    if (excludeNames.has(e.name)) continue;
    if (e.name.endsWith('.tmp') && e.name.startsWith('.')) continue;
    const abs = path.join(dir, e.name);
    let st: fs.Stats;
    try {
      st = await fsp.stat(abs);
    } catch {
      continue;
    }
    const node: TreeNode = {
      name: e.name,
      path: toPosixRelative(root, abs),
      type: st.isDirectory() ? 'dir' : 'file',
      size: st.isDirectory() ? 0 : st.size,
      mtime: st.mtime.toISOString(),
    };
    if (node.type === 'dir') node.children = await buildTree(root, abs, excludeNames, depth + 1);
    nodes.push(node);
  }
  nodes.sort((a, b) => (a.type === b.type ? a.name.localeCompare(b.name, 'zh-CN') : a.type === 'dir' ? -1 : 1));
  return nodes;
}

/** 文件名清洗：去掉路径分隔符与控制字符，保留中文；空则给默认名 */
export function sanitizeFileName(name: string, fallback = 'file'): string {
  const base = path.basename(String(name ?? '').replace(/\\/g, '/'));
  // eslint-disable-next-line no-control-regex
  const cleaned = base.replace(/[<>:"|?*\u0000-\u001f]/g, '_').replace(/^\.+$/, '').trim();
  return cleaned || fallback;
}
