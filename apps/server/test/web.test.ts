import { describe, expect, it } from 'vitest';
import { extractTitle, htmlToMarkdown } from '../src/services/html.js';
import { assertPublicUrl, isBlockedHostname, isPrivateAddress } from '../src/services/fetch.js';
import { rewriteMarkdownImages } from '../src/services/preview.js';

describe('htmlToMarkdown', () => {
  const html = `<!doctype html><html><head><title>测试 &amp; 页面</title><script>var x = 1;</script><style>.a{}</style></head>
<body><nav><a href="/">首页</a></nav>
<article>
<h1>标题一</h1>
<p>第一段 <strong>加粗</strong> 与 <a href="https://x.cn">链接</a>。</p>
<ul><li>甲</li><li>乙</li></ul>
<table><tr><th>年份</th><th>规模</th></tr><tr><td>2025</td><td>12.3 亿</td></tr></table>
<pre>code &lt;here&gt;</pre>
</article>
<footer>版权所有</footer></body></html>`;

  it('keeps headings/paragraphs/lists/tables and drops script/nav/footer', () => {
    const md = htmlToMarkdown(html);
    expect(md).toContain('# 标题一');
    expect(md).toContain('**加粗**');
    expect(md).toMatch(/- 甲\n- 乙/);
    expect(md).toContain('| 年份 | 规模 |');
    expect(md).toContain('| 2025 | 12.3 亿 |');
    expect(md).toContain('code <here>');
    expect(md).not.toContain('var x');
    expect(md).not.toContain('首页');
    expect(md).not.toContain('版权所有');
    // 代码块里解码出的 <here> 是正文内容；其余位置不应残留 HTML 标签
    expect(md.replace(/```[\s\S]*?```/g, '')).not.toMatch(/<\/?[a-z]+>/i);
  });

  it('extracts title and truncates with a note', () => {
    expect(extractTitle(html)).toBe('测试 & 页面');
    const md = htmlToMarkdown(html, { maxChars: 20 });
    expect(md.length).toBeLessThan(120);
    expect(md).toContain('截断');
  });
});

describe('web_fetch SSRF guard', () => {
  it('flags private, loopback, link-local and CGNAT ranges', () => {
    for (const ip of ['127.0.0.1', '10.1.2.3', '192.168.0.1', '172.16.0.1', '172.31.255.255', '169.254.169.254', '100.64.0.1', '0.0.0.0', '::1', 'fe80::1', 'fd00::1', '::ffff:10.0.0.1']) {
      expect(isPrivateAddress(ip), ip).toBe(true);
    }
    for (const ip of ['8.8.8.8', '172.32.0.1', '172.15.0.1', '2400:3200::1', '114.114.114.114']) {
      expect(isPrivateAddress(ip), ip).toBe(false);
    }
  });

  it('blocks localhost-ish hostnames and literal private IPs', () => {
    expect(isBlockedHostname('localhost')).toBe(true);
    expect(isBlockedHostname('foo.localhost')).toBe(true);
    expect(isBlockedHostname('router')).toBe(true);
    expect(isBlockedHostname('192.168.1.1')).toBe(true);
    expect(isBlockedHostname('[::1]')).toBe(true);
    expect(isBlockedHostname('www.example.com')).toBe(false);
  });

  it('assertPublicUrl rejects bad schemes and private targets', async () => {
    await expect(assertPublicUrl('ftp://example.com/x')).rejects.toThrow('http/https');
    await expect(assertPublicUrl('http://127.0.0.1:5178/api/settings')).rejects.toThrow('内网');
    await expect(assertPublicUrl('http://[::1]/')).rejects.toThrow('内网');
    await expect(assertPublicUrl('not a url')).rejects.toThrow('无法解析');
  });
});

describe('rewriteMarkdownImages', () => {
  it('rewrites relative images to the file endpoint, leaves absolute/data URLs', () => {
    const md = '![图1](../04_图表/市场规模.png "标题")\n![远程](https://cdn.example.com/a.png)\n<img src="charts/b.svg">\n![data](data:image/png;base64,AAA)';
    const out = rewriteMarkdownImages(md, 'p_20260101_000000_abcd', '05_正文/正文.md');
    expect(out).toContain('![图1](/api/projects/p_20260101_000000_abcd/file?path=04_%E5%9B%BE%E8%A1%A8%2F%E5%B8%82%E5%9C%BA%E8%A7%84%E6%A8%A1.png "标题")');
    expect(out).toContain('![远程](https://cdn.example.com/a.png)');
    expect(out).toContain('<img src="/api/projects/p_20260101_000000_abcd/file?path=05_%E6%AD%A3%E6%96%87%2Fcharts%2Fb.svg">');
    expect(out).toContain('![data](data:image/png;base64,AAA)');
  });

  it('never points outside the workspace', () => {
    const out = rewriteMarkdownImages('![x](../../etc/passwd)', 'p_20260101_000000_abcd', 'a.md');
    expect(out).toBe('![x](../../etc/passwd)');
  });
});
