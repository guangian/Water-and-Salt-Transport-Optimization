import { memo, useMemo } from 'react';
import ReactMarkdown, { type Components } from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { api, API_BASE, IS_MOCK } from '../../api';
import { cn } from '../../lib/cn';

export interface MarkdownProps {
  content: string;
  className?: string;
  /** 用于把相对图片路径或 /api/file 链接解析为可访问的 URL */
  projectId?: string;
  basePath?: string;
  variant?: 'doc' | 'stream';
  withIds?: boolean;
}

function resolveImage(src: string, projectId?: string, basePath?: string): string {
  if (!src) return src;
  if (/^(data:|https?:|blob:)/.test(src)) return src;
  if (projectId) {
    // 服务端已改写为 /api/projects/:id/file?path=…；Mock 模式映射到内置图片
    const m = /\/projects\/([^/]+)\/file\?path=([^&]+)/.exec(src);
    if (m) return IS_MOCK ? api.fileUrl(decodeURIComponent(m[1]), decodeURIComponent(m[2])) : src;
    if (src.startsWith(API_BASE)) return src;
    // 相对路径：相对于 basePath 所在目录
    const dir = basePath ? basePath.split('/').slice(0, -1) : [];
    const parts = src.split('/');
    for (const p of parts) {
      if (p === '..') dir.pop();
      else if (p !== '.') dir.push(p);
    }
    return api.fileUrl(projectId, dir.join('/'));
  }
  return src;
}

export function slugify(text: string): string {
  return text
    .toLowerCase()
    .replace(/[^\p{L}\p{N}\s-]/gu, '')
    .trim()
    .replace(/\s+/g, '-');
}

function textOf(node: unknown): string {
  if (node == null) return '';
  if (typeof node === 'string' || typeof node === 'number') return String(node);
  if (Array.isArray(node)) return node.map(textOf).join('');
  if (typeof node === 'object' && 'props' in (node as Record<string, unknown>)) {
    const props = (node as { props?: { children?: unknown } }).props;
    return textOf(props?.children);
  }
  return '';
}

export const Markdown = memo(function Markdown({ content, className, projectId, basePath, variant = 'doc', withIds = false }: MarkdownProps) {
  const components = useMemo<Components>(() => {
    const heading =
      (Tag: 'h1' | 'h2' | 'h3' | 'h4') =>
      ({ children, ...rest }: { children?: React.ReactNode }) => {
        const id = withIds ? slugify(textOf(children)) : undefined;
        return (
          <Tag id={id} {...rest}>
            {children}
          </Tag>
        );
      };
    return {
      img: ({ src, alt }) => <img src={resolveImage(src ?? '', projectId, basePath)} alt={alt ?? ''} loading="lazy" className="mx-auto" />,
      a: ({ href, children }) => (
        <a href={href} target="_blank" rel="noreferrer">
          {children}
        </a>
      ),
      table: ({ children }) => (
        <div className="not-prose my-3 overflow-x-auto rounded-lg border border-slate-200 dark:border-surface-dark-border">
          <table className="w-full border-collapse text-sm [&_td]:border-t [&_td]:border-slate-100 [&_td]:px-3 [&_td]:py-1.5 [&_th]:bg-slate-50 [&_th]:px-3 [&_th]:py-1.5 [&_th]:text-left [&_th]:font-semibold [&_th]:text-slate-700 dark:[&_td]:border-surface-dark-border dark:[&_th]:bg-slate-800 dark:[&_th]:text-slate-200">
            {children}
          </table>
        </div>
      ),
      h1: heading('h1'),
      h2: heading('h2'),
      h3: heading('h3'),
      h4: heading('h4'),
    };
  }, [projectId, basePath, withIds]);

  return (
    <div className={cn(variant === 'doc' ? 'prose-biz' : 'prose-stream', className)}>
      <ReactMarkdown remarkPlugins={[remarkGfm]} components={components}>
        {content}
      </ReactMarkdown>
    </div>
  );
});

export interface TocItem {
  id: string;
  text: string;
  level: number;
}

export function extractToc(markdown: string): TocItem[] {
  const out: TocItem[] = [];
  let inCode = false;
  for (const line of markdown.split('\n')) {
    if (line.trim().startsWith('```')) inCode = !inCode;
    if (inCode) continue;
    const m = /^(#{1,4})\s+(.+?)\s*#*\s*$/.exec(line);
    if (!m) continue;
    const text = m[2].replace(/[*_`]/g, '').trim();
    out.push({ id: slugify(text), text, level: m[1].length });
  }
  return out;
}
