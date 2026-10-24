/**
 * 轻量 HTML 处理：实体解码、去标签、HTML → Markdown（保留标题/段落/列表/表格）。
 * 不引入 DOM 库——正则足以应付搜索结果页与资讯正文，且零依赖。
 */

const NAMED: Record<string, string> = {
  amp: '&',
  lt: '<',
  gt: '>',
  quot: '"',
  apos: "'",
  nbsp: ' ',
  ensp: ' ',
  emsp: ' ',
  hellip: '…',
  mdash: '—',
  ndash: '–',
  lsquo: '‘',
  rsquo: '’',
  ldquo: '“',
  rdquo: '”',
  middot: '·',
  copy: '©',
  reg: '®',
  times: '×',
  laquo: '«',
  raquo: '»',
};

export function decodeEntities(s: string): string {
  return s.replace(/&(#x[0-9a-fA-F]+|#\d+|[a-zA-Z]+);/g, (m, body: string) => {
    if (body.startsWith('#x') || body.startsWith('#X')) {
      const cp = parseInt(body.slice(2), 16);
      return Number.isFinite(cp) ? safeFromCodePoint(cp, m) : m;
    }
    if (body.startsWith('#')) {
      const cp = parseInt(body.slice(1), 10);
      return Number.isFinite(cp) ? safeFromCodePoint(cp, m) : m;
    }
    return NAMED[body] ?? m;
  });
}

function safeFromCodePoint(cp: number, fallback: string): string {
  try {
    return String.fromCodePoint(cp);
  } catch {
    return fallback;
  }
}

export function stripTags(s: string): string {
  return decodeEntities(s.replace(/<[^>]+>/g, ''))
    .replace(/\s+/g, ' ')
    .trim();
}

const DROP_BLOCKS = ['script', 'style', 'noscript', 'svg', 'iframe', 'canvas', 'template', 'nav', 'footer', 'header', 'aside', 'form', 'button', 'select', 'video', 'audio'];

function removeBlocks(html: string): string {
  let out = html.replace(/<!--[\s\S]*?-->/g, '');
  for (const tag of DROP_BLOCKS) {
    out = out.replace(new RegExp(`<${tag}\\b[^>]*>[\\s\\S]*?<\\/${tag}\\s*>`, 'gi'), '\n');
  }
  // 常见的广告/分享/评论容器（按 class/id 关键字）
  out = out.replace(/<(div|section|ul)\b[^>]*(class|id)="[^"]*\b(comment|share|sidebar|advert|ad-|footer|breadcrumb|cookie|related|recommend)[^"]*"[^>]*>[\s\S]*?<\/\1\s*>/gi, '\n');
  return out;
}

function convertTables(html: string): string {
  return html.replace(/<table\b[^>]*>([\s\S]*?)<\/table>/gi, (_m, inner: string) => {
    const rows: string[][] = [];
    const rowRe = /<tr\b[^>]*>([\s\S]*?)<\/tr>/gi;
    let rm: RegExpExecArray | null;
    while ((rm = rowRe.exec(inner))) {
      const cells: string[] = [];
      const cellRe = /<t[hd]\b[^>]*>([\s\S]*?)<\/t[hd]>/gi;
      let cm: RegExpExecArray | null;
      while ((cm = cellRe.exec(rm[1] ?? ''))) cells.push(stripTags(cm[1] ?? '').replace(/\|/g, '\\|'));
      if (cells.length) rows.push(cells);
    }
    if (rows.length === 0) return '\n';
    const width = Math.max(...rows.map((r) => r.length));
    const norm = rows.map((r) => [...r, ...Array(width - r.length).fill('')]);
    const lines = [`| ${norm[0]!.join(' | ')} |`, `|${' --- |'.repeat(width)}`, ...norm.slice(1).map((r) => `| ${r.join(' | ')} |`)];
    return `\n\n${lines.join('\n')}\n\n`;
  });
}

function convertLists(html: string): string {
  let out = html;
  // 有序列表：按出现顺序编号
  out = out.replace(/<ol\b[^>]*>([\s\S]*?)<\/ol>/gi, (_m, inner: string) => {
    let i = 0;
    return '\n' + inner.replace(/<li\b[^>]*>([\s\S]*?)<\/li>/gi, (_x, li: string) => `\n${++i}. ${li.trim()}`) + '\n';
  });
  out = out.replace(/<ul\b[^>]*>([\s\S]*?)<\/ul>/gi, (_m, inner: string) => '\n' + inner.replace(/<li\b[^>]*>([\s\S]*?)<\/li>/gi, (_x, li: string) => `\n- ${li.trim()}`) + '\n');
  out = out.replace(/<li\b[^>]*>([\s\S]*?)<\/li>/gi, (_x, li: string) => `\n- ${li.trim()}`);
  return out;
}

export interface HtmlToMarkdownOptions {
  maxChars?: number;
  /** 保留链接为 [text](url)；默认 false 只留文本，避免正文被 URL 淹没 */
  keepLinks?: boolean;
}

export function extractTitle(html: string): string {
  const m = /<title[^>]*>([\s\S]*?)<\/title>/i.exec(html);
  if (m) return stripTags(m[1] ?? '');
  const og = /<meta[^>]*property="og:title"[^>]*content="([^"]*)"/i.exec(html);
  return og ? decodeEntities(og[1] ?? '') : '';
}

/** 优先取 <article> / <main> / 常见正文容器，命中则只转换该段，减少导航噪声 */
function pickMainContent(html: string): string {
  const candidates = [
    /<article\b[^>]*>([\s\S]*?)<\/article>/i,
    /<main\b[^>]*>([\s\S]*?)<\/main>/i,
    /<div\b[^>]*(?:class|id)="[^"]*\b(?:article-content|articleContent|article_content|rich_media_content|post-content|entry-content|content_area|main-content|news-content|detail-content|artical|RichText)\b[^"]*"[^>]*>([\s\S]*?)<\/div>\s*<\/div>/i,
  ];
  for (const re of candidates) {
    const m = re.exec(html);
    if (m && m[1] && stripTags(m[1]).length > 200) return m[1];
  }
  const body = /<body\b[^>]*>([\s\S]*?)<\/body>/i.exec(html);
  return body?.[1] ?? html;
}

export function htmlToMarkdown(html: string, opts: HtmlToMarkdownOptions = {}): string {
  const maxChars = opts.maxChars ?? 12000;
  let s = removeBlocks(html);
  s = pickMainContent(s);
  s = convertTables(s);
  // 这里不解码实体：&lt; 若提前变成 <，会被末尾的通用去标签正则吃掉；统一留给最后的 decodeEntities
  s = s.replace(/<pre\b[^>]*>([\s\S]*?)<\/pre>/gi, (_m, inner: string) => `\n\n\`\`\`\n${inner.replace(/<[^>]+>/g, '')}\n\`\`\`\n\n`);
  s = s.replace(/<h([1-6])\b[^>]*>([\s\S]*?)<\/h\1>/gi, (_m, lvl: string, inner: string) => `\n\n${'#'.repeat(Number(lvl))} ${stripTags(inner)}\n\n`);
  s = convertLists(s);
  s = s.replace(/<blockquote\b[^>]*>([\s\S]*?)<\/blockquote>/gi, (_m, inner: string) => `\n\n> ${stripTags(inner)}\n\n`);
  s = s.replace(/<(strong|b)\b[^>]*>([\s\S]*?)<\/\1>/gi, (_m, _t, inner: string) => {
    const t = stripTags(inner);
    return t ? `**${t}**` : '';
  });
  s = s.replace(/<(em|i)\b[^>]*>([\s\S]*?)<\/\1>/gi, (_m, _t, inner: string) => {
    const t = stripTags(inner);
    return t ? `*${t}*` : '';
  });
  s = s.replace(/<a\b[^>]*href="([^"]*)"[^>]*>([\s\S]*?)<\/a>/gi, (_m, href: string, inner: string) => {
    const t = stripTags(inner);
    if (!t) return '';
    return opts.keepLinks && /^https?:\/\//.test(href) ? `[${t}](${decodeEntities(href)})` : t;
  });
  s = s.replace(/<img\b[^>]*alt="([^"]*)"[^>]*>/gi, (_m, alt: string) => (alt ? `[图：${decodeEntities(alt)}]` : ''));
  s = s.replace(/<br\s*\/?>/gi, '\n');
  s = s.replace(/<\/(p|div|section|article|tr|dd|dt|figcaption)>/gi, '\n\n');
  s = s.replace(/<(p|div|section|article|dd|dt|figcaption)\b[^>]*>/gi, '\n');
  s = s.replace(/<[^>]+>/g, '');
  s = decodeEntities(s);
  s = s
    .split('\n')
    .map((line) => line.replace(/[ \t\u00a0]+/g, ' ').trim())
    .join('\n')
    .replace(/\n{3,}/g, '\n\n')
    .trim();
  if (s.length > maxChars) s = s.slice(0, maxChars) + `\n\n…（已截断，原文约 ${s.length} 字符）`;
  return s;
}
