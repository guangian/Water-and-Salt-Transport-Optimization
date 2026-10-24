import { useEffect, useMemo, useState } from 'react';
import { useNavigate, useParams } from 'react-router-dom';
import { BookOpen, FileText, List, Search } from 'lucide-react';
import { api } from '../api';
import { useAsync } from '../hooks/useAsync';
import { Page } from '../components/layout/AppShell';
import { Markdown, extractToc } from '../components/domain/Markdown';
import { EmptyState, Skeleton } from '../components/ui/Misc';
import { Input } from '../components/ui/Form';
import { cn } from '../lib/cn';

export default function KnowledgePage() {
  const { slug } = useParams();
  const navigate = useNavigate();
  const { data: docs, loading } = useAsync(() => api.getKnowledge(), []);
  const [query, setQuery] = useState('');
  const active = slug ?? docs?.[0]?.slug ?? null;
  const { data: doc, loading: docLoading } = useAsync(() => api.getKnowledgeDoc(active!), [active], { enabled: Boolean(active) });
  const [activeHeading, setActiveHeading] = useState<string | null>(null);

  useEffect(() => {
    if (!slug && docs?.length) navigate(`/knowledge/${docs[0].slug}`, { replace: true });
  }, [slug, docs, navigate]);

  const groups = useMemo(() => {
    const list = (docs ?? []).filter((d) => !query || d.title.includes(query) || d.summary?.includes(query) || d.group.includes(query));
    const map = new Map<string, typeof list>();
    for (const d of list) map.set(d.group, [...(map.get(d.group) ?? []), d]);
    return [...map.entries()];
  }, [docs, query]);

  const toc = useMemo(() => (doc ? extractToc(doc.content).filter((t) => t.level >= 2 && t.level <= 3) : []), [doc]);

  // 目录锁定：滚动时高亮当前章节
  useEffect(() => {
    if (!toc.length) return;
    const headings = toc.map((t) => document.getElementById(t.id)).filter((x): x is HTMLElement => Boolean(x));
    if (!headings.length) return;
    const observer = new IntersectionObserver(
      (entries) => {
        const visible = entries.filter((e) => e.isIntersecting).sort((a, b) => a.boundingClientRect.top - b.boundingClientRect.top);
        if (visible[0]) setActiveHeading(visible[0].target.id);
      },
      { rootMargin: '-72px 0px -70% 0px', threshold: [0, 1] },
    );
    headings.forEach((h) => observer.observe(h));
    setActiveHeading(headings[0].id);
    return () => observer.disconnect();
  }, [toc, doc?.slug]);

  return (
    <Page wide>
      <div className="grid grid-cols-1 gap-6 lg:grid-cols-[280px_minmax(0,1fr)] 2xl:grid-cols-[280px_minmax(0,1fr)_220px]">
        <aside className="card self-start lg:sticky lg:top-20">
          <div className="border-b border-slate-100 p-3 dark:border-surface-dark-border">
            <div className="relative">
              <Search className="pointer-events-none absolute left-2.5 top-1/2 h-4 w-4 -translate-y-1/2 text-slate-400" />
              <Input className="pl-8" placeholder="搜索知识库…" value={query} onChange={(e) => setQuery(e.target.value)} />
            </div>
          </div>
          <div className="max-h-[calc(100vh-11rem)] overflow-y-auto p-2">
            {loading ? (
              <div className="space-y-2 p-2">
                {Array.from({ length: 8 }).map((_, i) => (
                  <Skeleton key={i} className="h-9" />
                ))}
              </div>
            ) : (
              groups.map(([group, items]) => (
                <div key={group} className="mb-2">
                  <div className="px-2 pb-1 pt-2 text-[11px] font-semibold uppercase tracking-wide text-steel-500">{group}</div>
                  {items.map((d) => (
                    <button
                      key={d.slug}
                      onClick={() => navigate(`/knowledge/${d.slug}`)}
                      className={cn(
                        'flex w-full items-start gap-2.5 rounded-lg px-2.5 py-2 text-left text-sm transition',
                        d.slug === active ? 'bg-brand-50 text-brand-700 dark:bg-brand-900/50 dark:text-brand-100' : 'text-slate-700 hover:bg-slate-100 dark:text-slate-300 dark:hover:bg-slate-800',
                      )}
                    >
                      <FileText className={cn('mt-0.5 h-4 w-4 shrink-0', d.slug === active ? 'text-gold-500' : 'text-slate-400')} />
                      <span className="min-w-0">
                        <span className="block truncate font-medium">{d.title}</span>
                        {d.summary ? <span className="mt-0.5 line-clamp-1 block text-[11px] text-steel-500">{d.summary}</span> : null}
                      </span>
                    </button>
                  ))}
                </div>
              ))
            )}
          </div>
        </aside>

        <article className="card min-w-0">
          {docLoading && !doc ? (
            <div className="space-y-3 p-8">
              <Skeleton className="h-8 w-1/2" />
              <Skeleton className="h-4 w-full" />
              <Skeleton className="h-4 w-11/12" />
              <Skeleton className="h-40 w-full" />
            </div>
          ) : doc ? (
            <div className="px-8 py-7">
              <div className="mb-6 flex items-center gap-2 text-xs text-steel-500">
                <BookOpen className="h-3.5 w-3.5" />
                知识库 / {docs?.find((d) => d.slug === doc.slug)?.group ?? ''}
              </div>
              <Markdown content={doc.content} withIds />
            </div>
          ) : (
            <EmptyState icon={<BookOpen />} title="选择一篇文档" description="左侧列出赛道规范、官方评分表、获奖作品特征、评委 50 问、数据源清单与使用手册。" />
          )}
        </article>

        <aside className="hidden self-start 2xl:block 2xl:sticky 2xl:top-20">
          {toc.length ? (
            <div className="card p-4">
              <div className="mb-2 flex items-center gap-1.5 text-xs font-semibold uppercase tracking-wide text-steel-500">
                <List className="h-3.5 w-3.5" />
                目录
              </div>
              <ul className="space-y-0.5 border-l border-slate-200 dark:border-surface-dark-border">
                {toc.map((t) => (
                  <li key={t.id}>
                    <a
                      href={`#${t.id}`}
                      onClick={(e) => {
                        e.preventDefault();
                        document.getElementById(t.id)?.scrollIntoView({ behavior: 'smooth', block: 'start' });
                        setActiveHeading(t.id);
                      }}
                      className={cn(
                        '-ml-px block border-l-2 py-1 text-xs transition',
                        t.level === 3 ? 'pl-6' : 'pl-3',
                        activeHeading === t.id ? 'border-gold-500 font-medium text-brand-700 dark:text-gold-300' : 'border-transparent text-steel-500 hover:text-slate-800 dark:hover:text-slate-200',
                      )}
                    >
                      {t.text}
                    </a>
                  </li>
                ))}
              </ul>
            </div>
          ) : null}
        </aside>
      </div>
    </Page>
  );
}
