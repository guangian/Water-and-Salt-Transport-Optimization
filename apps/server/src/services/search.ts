import path from 'node:path';
import crypto from 'node:crypto';
import type { SearchProviderId, SearchResult, Settings } from '../types.js';
import { readJsonOrNull, writeJsonAtomic } from '../util/fs.js';
import { createLogger, errorMessage } from '../util/log.js';
import { decodeEntities, stripTags } from './html.js';

const log = createLogger('search');

export const BROWSER_UA =
  'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130.0.0.0 Safari/537.36';

export interface SearchOptions {
  count?: number;
  site?: string;
  /** day | week | month | year | noLimit */
  freshness?: string;
}

export interface SearchOutcome {
  results: SearchResult[];
  provider: SearchProviderId | null;
  errors: { provider: SearchProviderId; error: string }[];
  cached: boolean;
}

const MAX_RESULTS = 10;
const CACHE_TTL_MS = 24 * 60 * 60 * 1000;
const MIN_INTERVAL_MS: Record<SearchProviderId, number> = { bocha: 300, tavily: 300, serper: 300, bing: 1500, duckduckgo: 1500 };

type Fetcher = typeof fetch;

function clampCount(n: number | undefined, fallback: number): number {
  const v = Number.isFinite(n) ? Math.floor(n as number) : fallback;
  return Math.min(MAX_RESULTS, Math.max(1, v));
}

function withSite(query: string, site?: string): string {
  return site ? `${query} site:${site}` : query;
}

async function readJsonResponse(res: Response): Promise<unknown> {
  const text = await res.text();
  if (!res.ok) throw new Error(`HTTP ${res.status}: ${text.slice(0, 200)}`);
  try {
    return JSON.parse(text);
  } catch {
    throw new Error(`响应不是 JSON：${text.slice(0, 120)}`);
  }
}

function asRecord(v: unknown): Record<string, unknown> {
  return v && typeof v === 'object' ? (v as Record<string, unknown>) : {};
}

function asString(v: unknown): string {
  return typeof v === 'string' ? v : '';
}

// ---------- API 提供商 ----------

export async function searchBocha(query: string, key: string, opts: SearchOptions, f: Fetcher): Promise<SearchResult[]> {
  const res = await f('https://api.bochaai.com/v1/web-search', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${key}` },
    body: JSON.stringify({ query: withSite(query, opts.site), freshness: opts.freshness ?? 'noLimit', summary: true, count: clampCount(opts.count, 8) }),
    signal: AbortSignal.timeout(20000),
  });
  const data = asRecord(await readJsonResponse(res));
  const pages = asRecord(asRecord(data.data).webPages);
  const items = Array.isArray(pages.value) ? pages.value : [];
  return items.map((it) => {
    const r = asRecord(it);
    return {
      title: asString(r.name),
      url: asString(r.url),
      snippet: asString(r.summary) || asString(r.snippet),
      date: asString(r.datePublished) || asString(r.dateLastCrawled) || undefined,
      provider: 'bocha' as const,
    };
  });
}

export async function searchTavily(query: string, key: string, opts: SearchOptions, f: Fetcher): Promise<SearchResult[]> {
  const body: Record<string, unknown> = { api_key: key, query, max_results: clampCount(opts.count, 8), search_depth: 'basic' };
  if (opts.site) body.include_domains = [opts.site];
  const res = await f('https://api.tavily.com/search', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json', Authorization: `Bearer ${key}` },
    body: JSON.stringify(body),
    signal: AbortSignal.timeout(20000),
  });
  const data = asRecord(await readJsonResponse(res));
  const items = Array.isArray(data.results) ? data.results : [];
  return items.map((it) => {
    const r = asRecord(it);
    return { title: asString(r.title), url: asString(r.url), snippet: asString(r.content), date: asString(r.published_date) || undefined, provider: 'tavily' as const };
  });
}

export async function searchSerper(query: string, key: string, opts: SearchOptions, f: Fetcher): Promise<SearchResult[]> {
  const res = await f('https://google.serper.dev/search', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json', 'X-API-KEY': key },
    body: JSON.stringify({ q: withSite(query, opts.site), num: clampCount(opts.count, 8), gl: 'cn', hl: 'zh-cn' }),
    signal: AbortSignal.timeout(20000),
  });
  const data = asRecord(await readJsonResponse(res));
  const items = Array.isArray(data.organic) ? data.organic : [];
  return items.map((it) => {
    const r = asRecord(it);
    return { title: asString(r.title), url: asString(r.link), snippet: asString(r.snippet), date: asString(r.date) || undefined, provider: 'serper' as const };
  });
}

// ---------- 网页解析提供商 ----------

/** 解析 Bing 搜索结果页：li.b_algo > h2 > a（标题/链接），p 或 .b_caption p（摘要） */
export function parseBingHtml(html: string): SearchResult[] {
  const out: SearchResult[] = [];
  const parts = html.split(/<li[^>]*class="[^"]*\bb_algo\b[^"]*"[^>]*>/i).slice(1);
  for (const part of parts) {
    const block = part.split(/<\/li>/i)[0] ?? part;
    const h2 = /<h2[^>]*>[\s\S]*?<a[^>]*href="([^"]+)"[^>]*>([\s\S]*?)<\/a>/i.exec(block);
    if (!h2) continue;
    const url = decodeEntities(h2[1] ?? '');
    const title = stripTags(h2[2] ?? '');
    const pm = /<p[^>]*>([\s\S]*?)<\/p>/i.exec(block);
    let snippet = pm ? stripTags(pm[1] ?? '') : '';
    let date: string | undefined;
    const dm = /^(\d{4}年\d{1,2}月\d{1,2}日|\d{4}-\d{2}-\d{2}|[A-Z][a-z]{2} \d{1,2}, \d{4})\s*·\s*/.exec(snippet);
    if (dm) {
      date = dm[1];
      snippet = snippet.slice(dm[0].length);
    }
    if (url.startsWith('http') && title) out.push({ title, url, snippet, date, provider: 'bing' });
  }
  return out;
}

/** 解析 DuckDuckGo lite 页：a.result-link（标题/跳转链接，uddg 参数为真实 URL），td.result-snippet（摘要） */
export function parseDuckDuckGoLiteHtml(html: string): SearchResult[] {
  const out: SearchResult[] = [];
  const linkRe = /<a[^>]*class=['"]result-link['"][^>]*href=['"]([^'"]+)['"][^>]*>([\s\S]*?)<\/a>|<a[^>]*href=['"]([^'"]+)['"][^>]*class=['"]result-link['"][^>]*>([\s\S]*?)<\/a>/gi;
  const snippetRe = /<td[^>]*class=['"]result-snippet['"][^>]*>([\s\S]*?)<\/td>/gi;
  const links: { url: string; title: string }[] = [];
  let m: RegExpExecArray | null;
  while ((m = linkRe.exec(html))) {
    const href = decodeEntities(m[1] ?? m[3] ?? '');
    const title = stripTags(m[2] ?? m[4] ?? '');
    let url = href;
    const u = /[?&]uddg=([^&]+)/.exec(href);
    if (u) {
      try {
        url = decodeURIComponent(u[1] ?? '');
      } catch {
        url = href;
      }
    }
    if (url.startsWith('//')) url = 'https:' + url;
    if (url.startsWith('http') && title) links.push({ url, title });
  }
  const snippets: string[] = [];
  while ((m = snippetRe.exec(html))) snippets.push(stripTags(m[1] ?? ''));
  links.forEach((l, i) => out.push({ title: l.title, url: l.url, snippet: snippets[i] ?? '', provider: 'duckduckgo' }));
  return out;
}

/**
 * 查询词 → 匹配片段：拉丁/数字词（≥2 字符）+ 中文二元组。
 * 用于判断抓取到的结果页是否真的在回答这个问题。
 */
export function queryTokens(query: string): string[] {
  const tokens = new Set<string>();
  for (const w of query.toLowerCase().match(/[a-z0-9][a-z0-9.+-]*/g) ?? []) if (w.length >= 2) tokens.add(w);
  const cjk = query.replace(/[^\u3400-\u9fff]/g, ' ');
  for (const run of cjk.split(/\s+/)) {
    if (run.length === 1) tokens.add(run);
    for (let i = 0; i + 2 <= run.length; i++) tokens.add(run.slice(i, i + 2));
  }
  return [...tokens];
}

/**
 * 网页抓取型提供商（Bing/DDG）在被判定为爬虫时会返回一页“看起来正常”的无关结果，
 * 这里要求至少一条结果的标题或摘要命中查询词片段，否则视为无结果，交给下一个提供商。
 */
export function looksRelevant(query: string, results: SearchResult[]): boolean {
  const tokens = queryTokens(query);
  if (tokens.length === 0 || results.length === 0) return results.length > 0;
  return results.some((r) => {
    const hay = `${r.title} ${r.snippet}`.toLowerCase();
    return tokens.some((t) => hay.includes(t));
  });
}

export async function searchBing(query: string, opts: SearchOptions, f: Fetcher): Promise<SearchResult[]> {
  // cn.bing.com/search 现在 301 到 www.bing.com 首页（丢掉 /search），所以直接请求 www 并带 zh-CN 市场参数
  const url = `https://www.bing.com/search?q=${encodeURIComponent(withSite(query, opts.site))}&mkt=zh-CN&setlang=zh-cn&ensearch=0&count=${clampCount(opts.count, 10)}`;
  const res = await f(url, {
    headers: { 'User-Agent': BROWSER_UA, 'Accept-Language': 'zh-CN,zh;q=0.9', Accept: 'text/html' },
    signal: AbortSignal.timeout(20000),
    redirect: 'follow',
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const results = parseBingHtml(await res.text());
  // 无 Cookie 的请求常被 Bing 投喂整页无关结果（本机实测），过滤掉以免挡住后面的 DuckDuckGo
  return looksRelevant(query, results) ? results : [];
}

export async function searchDuckDuckGo(query: string, opts: SearchOptions, f: Fetcher): Promise<SearchResult[]> {
  const url = `https://lite.duckduckgo.com/lite/?q=${encodeURIComponent(withSite(query, opts.site))}&kl=cn-zh`;
  const res = await f(url, {
    headers: { 'User-Agent': BROWSER_UA, 'Accept-Language': 'zh-CN,zh;q=0.9', Accept: 'text/html' },
    signal: AbortSignal.timeout(20000),
    redirect: 'follow',
  });
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  const results = parseDuckDuckGoLiteHtml(await res.text());
  return looksRelevant(query, results) ? results : [];
}

// ---------- 服务：按顺序尝试 + 缓存 ----------

interface CacheEntry {
  ts: number;
  provider: SearchProviderId;
  results: SearchResult[];
}

export class SearchService {
  private mem = new Map<string, CacheEntry>();
  private fileLoaded = false;
  private lastCall = new Map<SearchProviderId, number>();
  private readonly cacheFile: string | null;

  constructor(
    private readonly getSettings: () => Promise<Settings['search']>,
    workspaceDir?: string,
    private readonly fetcher: Fetcher = fetch,
  ) {
    this.cacheFile = workspaceDir ? path.join(workspaceDir, '.cache', 'search.json') : null;
  }

  private key(query: string, opts: SearchOptions): string {
    return crypto.createHash('sha1').update(JSON.stringify([query.trim().toLowerCase(), opts.site ?? '', opts.freshness ?? '', clampCount(opts.count, 0)])).digest('hex');
  }

  private async loadFileCache(): Promise<void> {
    if (this.fileLoaded || !this.cacheFile) return;
    this.fileLoaded = true;
    const data = await readJsonOrNull<Record<string, CacheEntry>>(this.cacheFile);
    if (data) for (const [k, v] of Object.entries(data)) if (Date.now() - v.ts < CACHE_TTL_MS) this.mem.set(k, v);
  }

  private async persist(): Promise<void> {
    if (!this.cacheFile) return;
    const obj: Record<string, CacheEntry> = {};
    for (const [k, v] of this.mem) obj[k] = v;
    await writeJsonAtomic(this.cacheFile, obj).catch((e) => log.warn('cache write failed', errorMessage(e)));
  }

  private async throttle(p: SearchProviderId): Promise<void> {
    const last = this.lastCall.get(p) ?? 0;
    const wait = MIN_INTERVAL_MS[p] - (Date.now() - last);
    if (wait > 0) await new Promise((r) => setTimeout(r, wait));
    this.lastCall.set(p, Date.now());
  }

  async search(query: string, opts: SearchOptions = {}): Promise<SearchOutcome> {
    const q = query.trim();
    if (!q) return { results: [], provider: null, errors: [], cached: false };
    await this.loadFileCache();
    const cfg = await this.getSettings();
    const count = clampCount(opts.count, cfg.maxResults || 8);
    const normalized: SearchOptions = { ...opts, count };
    const k = this.key(q, normalized);
    const hit = this.mem.get(k);
    if (hit && Date.now() - hit.ts < CACHE_TTL_MS) return { results: hit.results.slice(0, count), provider: hit.provider, errors: [], cached: true };

    const errors: SearchOutcome['errors'] = [];
    for (const p of cfg.providers) {
      try {
        let results: SearchResult[] = [];
        if (p === 'bocha') {
          if (!cfg.bochaKey) continue;
          await this.throttle(p);
          results = await searchBocha(q, cfg.bochaKey, normalized, this.fetcher);
        } else if (p === 'tavily') {
          if (!cfg.tavilyKey) continue;
          await this.throttle(p);
          results = await searchTavily(q, cfg.tavilyKey, normalized, this.fetcher);
        } else if (p === 'serper') {
          if (!cfg.serperKey) continue;
          await this.throttle(p);
          results = await searchSerper(q, cfg.serperKey, normalized, this.fetcher);
        } else if (p === 'bing') {
          await this.throttle(p);
          results = await searchBing(q, normalized, this.fetcher);
        } else if (p === 'duckduckgo') {
          await this.throttle(p);
          results = await searchDuckDuckGo(q, normalized, this.fetcher);
        }
        results = dedupe(results).slice(0, count);
        if (results.length > 0) {
          this.mem.set(k, { ts: Date.now(), provider: p, results });
          await this.persist();
          return { results, provider: p, errors, cached: false };
        }
        errors.push({ provider: p, error: '无结果' });
      } catch (e) {
        errors.push({ provider: p, error: errorMessage(e) });
        log.warn(`provider ${p} failed`, errorMessage(e));
      }
    }
    return { results: [], provider: null, errors, cached: false };
  }
}

export function dedupe(results: SearchResult[]): SearchResult[] {
  const seen = new Set<string>();
  const out: SearchResult[] = [];
  for (const r of results) {
    const key = r.url.replace(/[#?].*$/, '').replace(/\/+$/, '').toLowerCase();
    if (!key || seen.has(key)) continue;
    seen.add(key);
    out.push({ ...r, title: r.title.trim(), snippet: r.snippet.trim() });
  }
  return out;
}
