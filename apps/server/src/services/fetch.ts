import path from 'node:path';
import fsp from 'node:fs/promises';
import dns from 'node:dns/promises';
import net from 'node:net';
import crypto from 'node:crypto';
import { ensureDir } from '../util/fs.js';
import { createLogger, errorMessage } from '../util/log.js';
import { extractTitle, htmlToMarkdown } from './html.js';
import { BROWSER_UA } from './search.js';
import type { PythonRunner } from './python.js';

const log = createLogger('fetch');

export const FETCH_TIMEOUT_MS = 20000;
export const FETCH_MAX_BYTES = 2 * 1024 * 1024;

export interface FetchedPage {
  title: string;
  url: string;
  finalUrl: string;
  contentType: string;
  content: string;
  truncated: boolean;
  fetchedAt: string;
}

/** 私网/回环地址判定（IPv4 + IPv6），防止 SSRF */
export function isPrivateAddress(ip: string): boolean {
  const v = net.isIP(ip);
  if (v === 4) {
    const parts = ip.split('.').map(Number);
    const [a, b] = [parts[0] ?? 0, parts[1] ?? 0];
    if (a === 127 || a === 10 || a === 0) return true;
    if (a === 192 && b === 168) return true;
    if (a === 172 && b >= 16 && b <= 31) return true;
    if (a === 169 && b === 254) return true;
    if (a === 100 && b >= 64 && b <= 127) return true; // CGNAT
    return false;
  }
  if (v === 6) {
    const low = ip.toLowerCase();
    if (low === '::1' || low === '::') return true;
    if (low.startsWith('fe80:') || low.startsWith('fc') || low.startsWith('fd')) return true;
    const mapped = /^::ffff:(\d+\.\d+\.\d+\.\d+)$/.exec(low);
    if (mapped) return isPrivateAddress(mapped[1] ?? '');
    return false;
  }
  return false;
}

/** 主机名层面的拒绝：localhost / 无点主机名 / 直接写的私网 IP */
export function isBlockedHostname(hostname: string): boolean {
  const h = hostname.replace(/^\[|\]$/g, '').toLowerCase();
  if (!h || h === 'localhost' || h.endsWith('.localhost') || h.endsWith('.local') || h.endsWith('.internal')) return true;
  if (net.isIP(h)) return isPrivateAddress(h);
  if (!h.includes('.')) return true;
  return false;
}

export async function assertPublicUrl(raw: string): Promise<URL> {
  let url: URL;
  try {
    url = new URL(raw);
  } catch {
    throw new Error(`URL 无法解析：${raw}`);
  }
  if (url.protocol !== 'http:' && url.protocol !== 'https:') throw new Error('只支持 http/https');
  if (isBlockedHostname(url.hostname)) throw new Error('拒绝访问内网/本机地址');
  if (!net.isIP(url.hostname.replace(/^\[|\]$/g, ''))) {
    try {
      const addrs = await dns.lookup(url.hostname, { all: true });
      if (addrs.some((a) => isPrivateAddress(a.address))) throw new Error('拒绝访问解析到内网的地址');
    } catch (e) {
      if (e instanceof Error && e.message.startsWith('拒绝')) throw e;
      // DNS 失败交给 fetch 报错，信息更直观
    }
  }
  return url;
}

async function readBodyLimited(res: Response, limit: number): Promise<{ buf: Buffer; truncated: boolean }> {
  if (!res.body) return { buf: Buffer.alloc(0), truncated: false };
  const reader = res.body.getReader();
  const chunks: Uint8Array[] = [];
  let total = 0;
  let truncated = false;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    if (!value) continue;
    if (total + value.byteLength > limit) {
      chunks.push(value.subarray(0, limit - total));
      total = limit;
      truncated = true;
      await reader.cancel().catch(() => undefined);
      break;
    }
    chunks.push(value);
    total += value.byteLength;
  }
  return { buf: Buffer.concat(chunks.map((c) => Buffer.from(c.buffer, c.byteOffset, c.byteLength))), truncated };
}

/** 从 Content-Type / <meta charset> 推断编码并解码（国内站点仍有 GBK） */
export function decodeBody(buf: Buffer, contentType: string): string {
  let charset = /charset=([\w-]+)/i.exec(contentType)?.[1];
  if (!charset) {
    const head = buf.subarray(0, 4096).toString('latin1');
    charset = /<meta[^>]*charset=["']?([\w-]+)/i.exec(head)?.[1];
  }
  const cs = (charset ?? 'utf-8').toLowerCase();
  try {
    return new TextDecoder(cs === 'gb2312' ? 'gbk' : cs, { fatal: false }).decode(buf);
  } catch {
    return buf.toString('utf8');
  }
}

export interface FetchDeps {
  python?: PythonRunner;
  /** 工作区（PDF 临时文件与 preread 输出放在 <ws>/.cache/fetch） */
  workspace?: string;
  fetcher?: typeof fetch;
}

export async function fetchPage(rawUrl: string, maxChars = 12000, deps: FetchDeps = {}): Promise<FetchedPage> {
  const url = await assertPublicUrl(rawUrl);
  const f = deps.fetcher ?? fetch;
  const res = await f(url, {
    headers: {
      'User-Agent': BROWSER_UA,
      Accept: 'text/html,application/xhtml+xml,application/pdf,application/json,text/plain;q=0.9,*/*;q=0.8',
      'Accept-Language': 'zh-CN,zh;q=0.9,en;q=0.8',
    },
    redirect: 'follow',
    signal: AbortSignal.timeout(FETCH_TIMEOUT_MS),
  });
  if (!res.ok) throw new Error(`HTTP ${res.status} ${res.statusText}`);
  const finalUrl = res.url || url.toString();
  if (finalUrl !== url.toString()) await assertPublicUrl(finalUrl);
  const contentType = (res.headers.get('content-type') ?? '').toLowerCase();
  const { buf, truncated } = await readBodyLimited(res, FETCH_MAX_BYTES);
  const fetchedAt = new Date().toISOString();

  if (contentType.includes('application/pdf') || finalUrl.toLowerCase().split('?')[0]!.endsWith('.pdf')) {
    const content = await pdfToText(buf, deps, maxChars);
    return { title: path.basename(url.pathname) || 'PDF', url: rawUrl, finalUrl, contentType: 'application/pdf', content, truncated, fetchedAt };
  }
  const text = decodeBody(buf, contentType);
  if (contentType.includes('html') || /^\s*<(!doctype|html)/i.test(text.slice(0, 200))) {
    const title = extractTitle(text);
    return { title, url: rawUrl, finalUrl, contentType: 'text/html', content: htmlToMarkdown(text, { maxChars }), truncated, fetchedAt };
  }
  const trimmed = text.length > maxChars ? text.slice(0, maxChars) + `\n\n…（已截断，原文约 ${text.length} 字符）` : text;
  return { title: path.basename(url.pathname) || url.hostname, url: rawUrl, finalUrl, contentType: contentType || 'text/plain', content: trimmed, truncated, fetchedAt };
}

async function pdfToText(buf: Buffer, deps: FetchDeps, maxChars: number): Promise<string> {
  if (!deps.python) return '（PDF 内容：未配置 Python，无法抽取文本）';
  const dir = path.join(deps.workspace ?? process.cwd(), '.cache', 'fetch');
  await ensureDir(dir);
  const name = `${crypto.randomBytes(6).toString('hex')}.pdf`;
  const file = path.join(dir, name);
  await fsp.writeFile(file, buf);
  try {
    const r = await deps.python.module('preread', [file, '--out', dir], { workspace: deps.workspace, timeoutMs: 120000 });
    const out = typeof r.json?.out === 'string' ? r.json.out : path.join(dir, `${name}.txt`);
    const text = await fsp.readFile(out, 'utf8').catch(() => '');
    if (!text) return `（PDF 文本抽取失败：${r.stderr.slice(-300) || '无输出'}）`;
    return text.length > maxChars ? text.slice(0, maxChars) + `\n\n…（已截断，原文约 ${text.length} 字符）` : text;
  } catch (e) {
    log.warn('pdf preread failed', errorMessage(e));
    return `（PDF 文本抽取失败：${errorMessage(e)}）`;
  } finally {
    await fsp.rm(file, { force: true }).catch(() => undefined);
  }
}
