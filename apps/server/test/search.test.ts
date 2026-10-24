import path from 'node:path';
import os from 'node:os';
import fsp from 'node:fs/promises';
import { afterAll, describe, expect, it } from 'vitest';
import { SearchService, dedupe, looksRelevant, parseBingHtml, parseDuckDuckGoLiteHtml, queryTokens } from '../src/services/search.js';
import { DEFAULT_SETTINGS } from '../src/services/settings.js';

/** 按 cn.bing.com 实际结构裁剪的样例（含日期前缀、实体转义、嵌套标签） */
const BING_HTML = `
<html><body><ol id="b_results">
<li class="b_algo" data-id="1"><div class="b_tpcn"><a href="https://www.example.com/a"></a></div>
  <h2><a href="https://www.jdcm.cn/news/2026/brand-contest?utm=1" h="ID=SERP">全国高校商业精英挑战赛 <strong>2026</strong> 品牌策划赛道 &amp; 赛程</a></h2>
  <div class="b_caption"><p class="b_lineclamp2">2026年3月12日 · 第十四届全国高校商业精英挑战赛品牌策划竞赛正式启动，<b>报名</b>截止 4 月 30 日 …</p></div></li>
<li class="b_algo"><h2><a href="https://brand.bizchallenge.cn/rules">品牌策划竞赛规则</a></h2><p>策划书不超过 80 页，含品牌现状分析、品牌定位、品牌元素设计……</p></li>
<li class="b_ad"><h2><a href="https://ad.example.com">广告</a></h2><p>广告内容</p></li>
<li class="b_algo"><h2><span>没有链接的条目</span></h2><p>应被跳过</p></li>
</ol></body></html>`;

const DDG_HTML = `
<html><body><table>
<tr><td><a rel="nofollow" href="//duckduckgo.com/l/?uddg=https%3A%2F%2Fwww.example.cn%2F2026%2Fbrand&amp;rut=abc" class='result-link'>2026 品牌策划竞赛通知</a></td></tr>
<tr><td class='result-snippet'>示例 <b>摘要</b> 一</td></tr>
<tr><td><a class="result-link" href="https://direct.example.org/page">直接链接标题</a></td></tr>
<tr><td class="result-snippet">摘要二</td></tr>
</table></body></html>`;

describe('HTML search parsers', () => {
  it('parses Bing results, strips nested tags, decodes entities and extracts dates', () => {
    const r = parseBingHtml(BING_HTML);
    expect(r).toHaveLength(2);
    expect(r[0]).toEqual({
      title: '全国高校商业精英挑战赛 2026 品牌策划赛道 & 赛程',
      url: 'https://www.jdcm.cn/news/2026/brand-contest?utm=1',
      snippet: '第十四届全国高校商业精英挑战赛品牌策划竞赛正式启动，报名截止 4 月 30 日 …',
      date: '2026年3月12日',
      provider: 'bing',
    });
    expect(r[1]).toMatchObject({ title: '品牌策划竞赛规则', url: 'https://brand.bizchallenge.cn/rules', date: undefined });
    expect(r[1]!.snippet).toContain('80 页');
  });

  it('parses DuckDuckGo lite results and unwraps uddg redirect URLs', () => {
    const r = parseDuckDuckGoLiteHtml(DDG_HTML);
    expect(r).toEqual([
      { title: '2026 品牌策划竞赛通知', url: 'https://www.example.cn/2026/brand', snippet: '示例 摘要 一', provider: 'duckduckgo' },
      { title: '直接链接标题', url: 'https://direct.example.org/page', snippet: '摘要二', provider: 'duckduckgo' },
    ]);
  });

  it('returns empty arrays for unrelated HTML', () => {
    expect(parseBingHtml('<html><body>验证码</body></html>')).toEqual([]);
    expect(parseDuckDuckGoLiteHtml('<html></html>')).toEqual([]);
  });

  it('looksRelevant rejects anti-bot junk pages that share no token with the query', () => {
    expect(queryTokens('全国高校商业精英挑战赛 2026 品牌策划 site:gov.cn')).toEqual(expect.arrayContaining(['2026', '全国', '品牌', '策划', 'site', 'gov.cn']));
    const junk = [
      { title: '先天性颈椎融合畸形 - 百度百科', url: 'https://a.cn/1', snippet: '医学条目', provider: 'bing' as const },
      { title: '知网查重 - 国际版', url: 'https://a.cn/2', snippet: '论文查重', provider: 'bing' as const },
    ];
    expect(looksRelevant('全国高校商业精英挑战赛 2026 品牌策划', junk)).toBe(false);
    expect(looksRelevant('全国高校商业精英挑战赛 2026 品牌策划', parseBingHtml(BING_HTML))).toBe(true);
    expect(looksRelevant('bizgold', [{ title: 'BizGold 服务端', url: 'https://x', snippet: '', provider: 'bing' }])).toBe(true);
  });

  it('dedupe collapses URLs differing only by query/hash/trailing slash', () => {
    const r = dedupe([
      { title: 'a', url: 'https://x.cn/p/', snippet: '', provider: 'bing' },
      { title: 'b', url: 'https://x.cn/p?utm=1', snippet: '', provider: 'bing' },
      { title: 'c', url: 'https://X.cn/p#frag', snippet: '', provider: 'bing' },
      { title: 'd', url: 'https://x.cn/q', snippet: '', provider: 'bing' },
    ]);
    expect(r.map((x) => x.title)).toEqual(['a', 'd']);
  });
});

describe('SearchService fallback chain', () => {
  const ws = path.join(os.tmpdir(), `bizgold-search-${process.pid}`);
  afterAll(async () => {
    await fsp.rm(ws, { recursive: true, force: true });
  });

  it('skips API providers without keys, falls through failures, caches the first non-empty result', async () => {
    const calls: string[] = [];
    const fetcher: typeof fetch = async (input) => {
      const url = String(input instanceof Request ? input.url : input);
      calls.push(new URL(url).hostname);
      if (url.includes('bing.com')) return new Response('<html>blocked</html>', { status: 200, headers: { 'content-type': 'text/html' } });
      if (url.includes('duckduckgo.com')) return new Response(DDG_HTML, { status: 200, headers: { 'content-type': 'text/html' } });
      throw new Error('unexpected ' + url);
    };
    const svc = new SearchService(async () => ({ ...DEFAULT_SETTINGS.search, providers: ['bocha', 'tavily', 'serper', 'bing', 'duckduckgo'], maxResults: 5 }), ws, fetcher);
    const out = await svc.search('全国高校商业精英挑战赛 2026 品牌策划', { count: 1 });
    expect(out.provider).toBe('duckduckgo');
    expect(out.results).toHaveLength(1);
    expect(out.errors).toEqual([{ provider: 'bing', error: '无结果' }]);
    // 没有 key 的 API 提供商根本不该被请求
    expect(calls).toEqual(['www.bing.com', 'lite.duckduckgo.com']);

    const again = await svc.search('全国高校商业精英挑战赛 2026 品牌策划', { count: 1 });
    expect(again.cached).toBe(true);
    expect(calls).toHaveLength(2);
    // 落盘缓存
    const cached = JSON.parse(await fsp.readFile(path.join(ws, '.cache', 'search.json'), 'utf8')) as Record<string, unknown>;
    expect(Object.keys(cached)).toHaveLength(1);
  });

  it('returns an empty list with per-provider errors when everything fails', async () => {
    const fetcher: typeof fetch = async () => {
      throw new Error('ENETUNREACH');
    };
    const svc = new SearchService(async () => ({ ...DEFAULT_SETTINGS.search, providers: ['bing', 'duckduckgo'] }), undefined, fetcher);
    const out = await svc.search('x');
    expect(out.results).toEqual([]);
    expect(out.provider).toBeNull();
    expect(out.errors.map((e) => e.provider)).toEqual(['bing', 'duckduckgo']);
  });
});
