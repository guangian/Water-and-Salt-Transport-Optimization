import path from 'node:path';
import fsp from 'node:fs/promises';
import * as XLSX from 'xlsx';
import { isInside, safeResolve, toPosixRelative } from '../util/fs.js';
import { ApiError } from '../types.js';

export type PreviewKind = 'markdown' | 'image' | 'slides' | 'table' | 'pdf' | 'text' | 'binary';

export interface PreviewResult {
  kind: PreviewKind;
  path: string;
  size: number;
  content?: unknown;
  /** 前端幻灯片预览用的统一结构（slides.json 原样；工具包 .preview.json 页数组会被归一化） */
  slides?: SlidesDeck;
  url?: string;
  sheets?: { name: string; rows: unknown[][]; totalRows: number }[];
  truncated?: boolean;
}

export interface SlidesDeck {
  title: string;
  subtitle?: string;
  slides: Record<string, unknown>[];
  [k: string]: unknown;
}

const CHART_TYPE_MAP: Record<string, string> = { bar: 'bar', line: 'line', pie: 'pie', column_clustered: 'column_clustered', stacked: 'column_clustered', doughnut: 'pie' };

/**
 * 两种输入统一成前端 SlidesPreview 需要的 deck：
 * - 06_路演/slides.json：本来就是 deck，图片路径相对 slides.json 所在目录（前端按 basePath 解析）；
 * - 产出/路演.pptx.preview.json：工具包输出的「页数组」，图片路径相对工作区根，这里改写成相对 preview 文件目录的路径，
 *   同时把 texts/table/chart 映射到 Slide 字段，让同一个渲染器能画。
 */
export function toSlidesDeck(parsed: unknown, relFile: string): SlidesDeck {
  if (parsed && typeof parsed === 'object' && !Array.isArray(parsed) && Array.isArray((parsed as { slides?: unknown }).slides)) {
    return parsed as SlidesDeck;
  }
  const pages = Array.isArray(parsed) ? (parsed as Record<string, unknown>[]) : [];
  const baseDir = path.posix.dirname(relFile.replace(/\\/g, '/'));
  const toRel = (img: unknown): string | undefined => {
    if (typeof img !== 'string' || !img) return undefined;
    if (/^(data:|https?:|blob:)/.test(img)) return img;
    return path.posix.relative(baseDir === '.' ? '' : baseDir, img.replace(/\\/g, '/')) || img;
  };
  const slides = pages.map((pg) => {
    const texts = Array.isArray(pg.texts) ? (pg.texts as unknown[]).filter((x): x is string => typeof x === 'string') : [];
    const table = pg.table && typeof pg.table === 'object' ? (pg.table as { columns?: string[]; rows?: string[][] }) : null;
    const chart = pg.chart && typeof pg.chart === 'object' ? (pg.chart as Record<string, unknown>) : null;
    const layout = typeof pg.layout === 'string' ? pg.layout : 'bullets';
    const slide: Record<string, unknown> = { layout, title: pg.title ?? '', note: pg.note ?? undefined };
    if (layout === 'cover' || layout === 'section' || layout === 'closing') slide.subtitle = texts[0];
    else if (layout === 'agenda') slide.items = texts;
    else if (layout === 'quote') { slide.quote = typeof pg.quote === 'string' ? pg.quote : texts[0]; slide.by = texts[1]; }
    else slide.bullets = texts;
    const img = toRel(pg.image);
    if (img) slide.image = img;
    if (chart) slide.chart = { ...chart, type: CHART_TYPE_MAP[String(chart.type)] ?? 'bar' };
    if (typeof chart?.insight === 'string') slide.insight = chart.insight;
    if (table) { slide.columns = table.columns ?? []; slide.rows = table.rows ?? []; }
    for (const k of ['kpis', 'steps', 'members'] as const) if (Array.isArray(pg[k])) slide[k] = pg[k];
    return slide;
  });
  const first = pages[0] ?? {};
  return { title: typeof first.title === 'string' ? first.title : '路演', slides };
}

const IMAGE_EXT = new Set(['.png', '.jpg', '.jpeg', '.gif', '.webp', '.svg', '.bmp']);
const TEXT_EXT = new Set(['.txt', '.json', '.jsonl', '.csv', '.tsv', '.py', '.ts', '.js', '.html', '.htm', '.css', '.yaml', '.yml', '.toml', '.ini', '.log', '.xml', '.sql', '.sh', '.bat', '.ps1']);
const TABLE_EXT = new Set(['.xlsx', '.xlsm', '.xls', '.ods']);
const MAX_TEXT_BYTES = 400 * 1024;
const TABLE_ROWS = 200;

export const MIME: Record<string, string> = {
  '.md': 'text/markdown; charset=utf-8',
  '.txt': 'text/plain; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.jsonl': 'application/x-ndjson; charset=utf-8',
  '.csv': 'text/csv; charset=utf-8',
  '.html': 'text/html; charset=utf-8',
  '.htm': 'text/html; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.py': 'text/x-python; charset=utf-8',
  '.png': 'image/png',
  '.jpg': 'image/jpeg',
  '.jpeg': 'image/jpeg',
  '.gif': 'image/gif',
  '.webp': 'image/webp',
  '.svg': 'image/svg+xml',
  '.bmp': 'image/bmp',
  '.pdf': 'application/pdf',
  '.docx': 'application/vnd.openxmlformats-officedocument.wordprocessingml.document',
  '.pptx': 'application/vnd.openxmlformats-officedocument.presentationml.presentation',
  '.xlsx': 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
  '.zip': 'application/zip',
};

export function mimeFor(file: string): string {
  return MIME[path.extname(file).toLowerCase()] ?? 'application/octet-stream';
}

export function fileUrl(projectId: string, relPath: string): string {
  return `/api/projects/${encodeURIComponent(projectId)}/file?path=${encodeURIComponent(relPath)}`;
}

/**
 * 把 Markdown 中的相对图片链接改写为 /api/projects/:id/file?path=…；
 * http(s)/data: 链接与绝对路径不动。相对路径按 md 文件所在目录解析，再转为工作区相对路径。
 */
export function rewriteMarkdownImages(md: string, projectId: string, mdRelPath: string): string {
  const baseDir = path.posix.dirname(mdRelPath.replace(/\\/g, '/'));
  const resolveRel = (target: string): string | null => {
    const t = target.trim().replace(/\\/g, '/');
    if (!t || /^(https?:|data:|blob:|\/api\/)/i.test(t) || t.startsWith('/')) return null;
    const joined = path.posix.normalize(path.posix.join(baseDir === '.' ? '' : baseDir, t.split('#')[0]!.split('?')[0]!));
    if (joined.startsWith('..')) return null;
    return joined;
  };
  // ![alt](path "title")
  let out = md.replace(/!\[([^\]]*)\]\(\s*<?([^)\s>]+)>?(\s+"[^"]*")?\s*\)/g, (m, alt: string, target: string, title: string | undefined) => {
    const rel = resolveRel(decodeURIComponentSafe(target));
    return rel ? `![${alt}](${fileUrl(projectId, rel)}${title ?? ''})` : m;
  });
  // <img src="path">
  out = out.replace(/(<img\b[^>]*\bsrc=")([^"]+)(")/gi, (m, pre: string, target: string, post: string) => {
    const rel = resolveRel(decodeURIComponentSafe(target));
    return rel ? `${pre}${fileUrl(projectId, rel)}${post}` : m;
  });
  return out;
}

function decodeURIComponentSafe(s: string): string {
  try {
    return decodeURIComponent(s);
  } catch {
    return s;
  }
}

export async function previewFile(workspace: string, projectId: string, relPath: string): Promise<PreviewResult> {
  const abs = safeResolve(workspace, relPath);
  const st = await fsp.stat(abs).catch(() => null);
  if (!st || !st.isFile()) throw new ApiError(404, 'not_found', '文件不存在');
  // 与 /file 一致：符号链接解析后仍须落在工作区内
  const [realWs, real] = await Promise.all([fsp.realpath(workspace), fsp.realpath(abs)]);
  if (!isInside(realWs, real)) throw new ApiError(404, 'not_found', '文件不存在');
  const rel = toPosixRelative(workspace, abs);
  const ext = path.extname(abs).toLowerCase();
  const base = path.basename(abs).toLowerCase();
  const common = { path: rel, size: st.size };
  const url = fileUrl(projectId, rel);

  if (ext === '.md' || ext === '.markdown') {
    const text = await fsp.readFile(abs, 'utf8');
    return { kind: 'markdown', ...common, content: rewriteMarkdownImages(text, projectId, rel) };
  }
  if (IMAGE_EXT.has(ext)) return { kind: 'image', ...common, url };
  if (ext === '.pdf') return { kind: 'pdf', ...common, url };
  // 产出/路演.pptx 本身是二进制：旁边有工具包写的 路演.pptx.preview.json 就直接按幻灯片预览，没有再退回 binary
  if (ext === '.pptx') {
    const sidecar = `${abs}.preview.json`;
    if (await fsp.stat(sidecar).then((x) => x.isFile()).catch(() => false)) {
      return previewFile(workspace, projectId, `${rel}.preview.json`).then((r) => ({ ...r, path: rel, size: st.size, url }));
    }
  }
  if (base === 'slides.json' || base.endsWith('.preview.json')) {
    const text = await fsp.readFile(abs, 'utf8');
    try {
      const parsed = JSON.parse(text.replace(/^\uFEFF/, '')) as unknown;
      const deck = toSlidesDeck(parsed, rel);
      return { kind: 'slides', ...common, content: parsed, slides: deck };
    } catch {
      return { kind: 'text', ...common, content: text };
    }
  }
  if (TABLE_EXT.has(ext)) {
    const wb = XLSX.read(await fsp.readFile(abs), { type: 'buffer', cellDates: true, sheetRows: TABLE_ROWS + 1 });
    const sheets = wb.SheetNames.map((name) => {
      const ws = wb.Sheets[name];
      const rows = ws ? (XLSX.utils.sheet_to_json(ws, { header: 1, raw: true, defval: null, blankrows: false }) as unknown[][]) : [];
      const total = ws?.['!ref'] ? XLSX.utils.decode_range(ws['!ref']).e.r + 1 : rows.length;
      return { name, rows: rows.slice(0, TABLE_ROWS).map((r) => r.map(normalizeCell)), totalRows: total };
    });
    return { kind: 'table', ...common, sheets };
  }
  if (ext === '.csv' || ext === '.tsv') {
    const text = await fsp.readFile(abs, 'utf8');
    const lines = text.replace(/^\uFEFF/, '').split(/\r?\n/).filter((l) => l.length > 0);
    const sep = ext === '.tsv' ? '\t' : ',';
    const rows = lines.slice(0, TABLE_ROWS).map((l) => parseCsvLine(l, sep));
    return { kind: 'table', ...common, sheets: [{ name: path.basename(abs), rows, totalRows: lines.length }] };
  }
  if (TEXT_EXT.has(ext) || st.size < 64 * 1024) {
    const fh = await fsp.open(abs, 'r');
    try {
      const len = Math.min(st.size, MAX_TEXT_BYTES);
      const buf = Buffer.alloc(len);
      await fh.read(buf, 0, len, 0);
      if (looksBinary(buf)) return { kind: 'binary', ...common, url };
      return { kind: 'text', ...common, content: buf.toString('utf8'), truncated: st.size > MAX_TEXT_BYTES };
    } finally {
      await fh.close();
    }
  }
  return { kind: 'binary', ...common, url };
}

function normalizeCell(v: unknown): unknown {
  if (v instanceof Date) return v.toISOString().slice(0, 10);
  if (typeof v === 'number' && !Number.isInteger(v)) return Math.round(v * 1e6) / 1e6;
  return v;
}

function looksBinary(buf: Buffer): boolean {
  const n = Math.min(buf.length, 2048);
  let suspicious = 0;
  for (let i = 0; i < n; i++) {
    const c = buf[i]!;
    if (c === 0) return true;
    if (c < 7 || (c > 14 && c < 32)) suspicious++;
  }
  return n > 0 && suspicious / n > 0.1;
}

export function parseCsvLine(line: string, sep = ','): string[] {
  const out: string[] = [];
  let cur = '';
  let inQ = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i]!;
    if (inQ) {
      if (ch === '"') {
        if (line[i + 1] === '"') {
          cur += '"';
          i++;
        } else inQ = false;
      } else cur += ch;
    } else if (ch === '"') inQ = true;
    else if (ch === sep) {
      out.push(cur);
      cur = '';
    } else cur += ch;
  }
  out.push(cur);
  return out;
}
