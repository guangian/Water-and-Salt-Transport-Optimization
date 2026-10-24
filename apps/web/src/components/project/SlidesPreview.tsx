import { useCallback, useEffect, useMemo, useState } from 'react';
import { ChevronLeft, ChevronRight, Clock, Grid2x2, MessageSquareText, X } from 'lucide-react';
import type { SlidesDeck } from '../../types';
import { api } from '../../api';
import { resolveRelative } from '../../lib/paths';
import { cn } from '../../lib/cn';
import { LAYOUT_LABELS, ScaledSlide, SlideCanvas } from './SlideCanvas';

export function SlidesPreview({ deck, projectId, basePath }: { deck: SlidesDeck; projectId: string; basePath: string }) {
  const [active, setActive] = useState<number | null>(null);
  const total = deck.slides.length;
  const resolveImage = useCallback((src: string) => (/^(data:|https?:|blob:)/.test(src) ? src : api.fileUrl(projectId, resolveRelative(basePath, src))), [projectId, basePath]);
  const totalSeconds = useMemo(() => deck.slides.reduce((s, x) => s + (x.duration ?? 0), 0), [deck.slides]);

  useEffect(() => {
    if (active == null) return;
    const onKey = (e: KeyboardEvent) => {
      if (e.key === 'Escape') setActive(null);
      if (e.key === 'ArrowRight') setActive((a) => (a == null ? a : Math.min(total - 1, a + 1)));
      if (e.key === 'ArrowLeft') setActive((a) => (a == null ? a : Math.max(0, a - 1)));
    };
    window.addEventListener('keydown', onKey);
    return () => window.removeEventListener('keydown', onKey);
  }, [active, total]);

  const activeSlide = active != null ? deck.slides[active] : null;

  return (
    <div>
      <div className="mb-4 flex flex-wrap items-end justify-between gap-3">
        <div className="min-w-0">
          <h3 className="truncate text-base font-semibold text-slate-800 dark:text-slate-100">{deck.title}</h3>
          {deck.subtitle ? <p className="mt-0.5 truncate text-xs text-steel-500">{deck.subtitle}</p> : null}
        </div>
        <div className="flex flex-wrap items-center gap-1.5 text-xs">
          <span className="chip">
            <Grid2x2 className="h-3 w-3" />
            {total} 页
          </span>
          <span className="chip">{deck.aspect ?? '16:9'}</span>
          {deck.theme ? <span className="chip">主题 {deck.theme}</span> : null}
          {deck.minutes ? (
            <span className="chip">
              <Clock className="h-3 w-3" />
              目标 {deck.minutes} 分钟{totalSeconds ? ` · 配时 ${Math.round(totalSeconds / 60)} 分 ${totalSeconds % 60} 秒` : ''}
            </span>
          ) : null}
        </div>
      </div>

      <div className="grid grid-cols-2 gap-4 2xl:grid-cols-3">
        {deck.slides.map((s, i) => (
          <div key={i} className="group">
            <ScaledSlide className="rounded-lg border border-slate-200 shadow-sm ring-gold-400 transition group-hover:ring-2 dark:border-surface-dark-border" onClick={() => setActive(i)} title="点击放大">
              <SlideCanvas slide={s} index={i} total={total} deck={deck} resolveImage={resolveImage} />
            </ScaledSlide>
            <div className="mt-1.5 flex items-center gap-2 px-0.5 text-xs">
              <span className="w-5 shrink-0 font-semibold tabular text-steel-500">{i + 1}</span>
              <span className="min-w-0 flex-1 truncate text-slate-700 dark:text-slate-200" title={s.title}>
                {s.title ?? LAYOUT_LABELS[s.layout]}
              </span>
              <span className="shrink-0 rounded bg-slate-100 px-1.5 py-0.5 text-[10px] text-steel-500 dark:bg-slate-800">{LAYOUT_LABELS[s.layout] ?? s.layout}</span>
              {s.duration ? <span className="shrink-0 tabular text-steel-400">{s.duration}s</span> : null}
            </div>
          </div>
        ))}
      </div>

      {active != null && activeSlide ? (
        <div className="fixed inset-0 z-[95] flex flex-col bg-brand-900/90 backdrop-blur-sm" role="dialog" aria-modal="true" onClick={() => setActive(null)}>
          <div className="flex items-center justify-between px-6 py-3 text-white/80" onClick={(e) => e.stopPropagation()}>
            <div className="min-w-0 truncate text-sm">
              {deck.title} · 第 {active + 1} / {total} 页
              <span className="ml-2 rounded bg-white/10 px-1.5 py-0.5 text-[11px]">{LAYOUT_LABELS[activeSlide.layout] ?? activeSlide.layout}</span>
            </div>
            <button onClick={() => setActive(null)} className="rounded-md p-1.5 hover:bg-white/10" aria-label="关闭">
              <X className="h-5 w-5" />
            </button>
          </div>
          <div className="flex min-h-0 flex-1 items-center justify-center gap-4 px-6" onClick={(e) => e.stopPropagation()}>
            <button onClick={() => setActive(Math.max(0, active - 1))} disabled={active === 0} className="rounded-full bg-white/10 p-3 text-white hover:bg-white/20 disabled:opacity-30" aria-label="上一页">
              <ChevronLeft className="h-6 w-6" />
            </button>
            <div className="w-full max-w-[1200px]" style={{ maxHeight: 'calc(100vh - 220px)', aspectRatio: '16 / 9' }}>
              <ScaledSlide className="rounded-lg shadow-2xl">
                <SlideCanvas slide={activeSlide} index={active} total={total} deck={deck} resolveImage={resolveImage} />
              </ScaledSlide>
            </div>
            <button onClick={() => setActive(Math.min(total - 1, active + 1))} disabled={active === total - 1} className="rounded-full bg-white/10 p-3 text-white hover:bg-white/20 disabled:opacity-30" aria-label="下一页">
              <ChevronRight className="h-6 w-6" />
            </button>
          </div>
          <div className="mx-auto w-full max-w-[1200px] px-6 py-4" onClick={(e) => e.stopPropagation()}>
            <div className={cn('rounded-xl border border-white/10 bg-white/5 px-4 py-3 text-sm text-white/85')}>
              <div className="mb-1 flex items-center gap-2 text-xs text-[#EBD88E]">
                <MessageSquareText className="h-3.5 w-3.5" />
                演讲者备注
                {activeSlide.duration ? <span className="ml-auto tabular text-white/60">建议 {activeSlide.duration} 秒</span> : null}
              </div>
              <p className="leading-relaxed">{activeSlide.note ?? '（本页无备注）'}</p>
            </div>
          </div>
        </div>
      ) : null}
    </div>
  );
}
