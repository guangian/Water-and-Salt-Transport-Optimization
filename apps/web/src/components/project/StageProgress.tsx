import type { ReactNode } from 'react';
import { BarChart3, Calculator, Check, Compass, FileSearch, Gavel, Globe, Loader2, Package, PenLine, Presentation, SkipForward, X } from 'lucide-react';
import type { StageId, StageState } from '../../types';
import { STAGES, stageDurationMs } from '../../lib/stages';
import { formatDuration } from '../../lib/format';
import { cn } from '../../lib/cn';

const STAGE_ICONS: Record<StageId, ReactNode> = {
  parse: <FileSearch />,
  research: <Globe />,
  strategy: <Compass />,
  finance: <Calculator />,
  charts: <BarChart3 />,
  write: <PenLine />,
  pitch: <Presentation />,
  review: <Gavel />,
  package: <Package />,
};

export function StageProgress({ stages, now, className }: { stages: Record<StageId, StageState>; now: number; className?: string }) {
  return (
    <ol className={cn('grid grid-cols-9 gap-0', className)} aria-label="阶段进度">
      {STAGES.map((def, i) => {
        const st = stages[def.id] ?? { status: 'pending' };
        const status = st.status;
        const prevDone = i > 0 && ['done', 'skipped'].includes(stages[STAGES[i - 1].id]?.status ?? 'pending');
        const duration = stageDurationMs(st, now);
        const ring =
          status === 'done'
            ? 'border-brand-600 bg-brand-600 text-white dark:border-brand-300 dark:bg-brand-300 dark:text-brand-900'
            : status === 'running'
              ? 'border-emerald-500 bg-emerald-50 text-emerald-600 shadow-[0_0_0_4px_rgba(16,185,129,0.15)] dark:bg-emerald-900/40 dark:text-emerald-300'
              : status === 'failed'
                ? 'border-red-500 bg-red-50 text-red-600 dark:bg-red-900/40 dark:text-red-300'
                : status === 'skipped'
                  ? 'border-slate-300 bg-slate-100 text-slate-400 dark:border-slate-600 dark:bg-slate-800'
                  : 'border-slate-200 bg-white text-slate-300 dark:border-slate-700 dark:bg-surface-dark-card dark:text-slate-600';
        const label =
          status === 'running' ? (
            <span className="text-emerald-600 dark:text-emerald-300">进行中 · {formatDuration(duration, { compact: true })}</span>
          ) : status === 'done' ? (
            <span>{formatDuration(duration, { compact: true })}</span>
          ) : status === 'failed' ? (
            <span className="text-red-600 dark:text-red-300">失败</span>
          ) : status === 'skipped' ? (
            <span>跳过</span>
          ) : (
            <span className="text-slate-400 dark:text-slate-600">等待</span>
          );
        return (
          <li key={def.id} className="relative flex flex-col items-center px-1 text-center" title={`${def.name} · ${def.description}${st.note ? `\n${st.note}` : ''}`}>
            {i > 0 ? <span className={cn('absolute left-0 right-1/2 top-[18px] h-0.5 -translate-y-1/2', prevDone ? 'bg-brand-500 dark:bg-brand-300' : 'bg-slate-200 dark:bg-slate-700')} style={{ right: 'calc(50% + 18px)' }} /> : null}
            {i < STAGES.length - 1 ? (
              <span className={cn('absolute right-0 top-[18px] h-0.5 -translate-y-1/2', status === 'done' || status === 'skipped' ? 'bg-brand-500 dark:bg-brand-300' : 'bg-slate-200 dark:bg-slate-700')} style={{ left: 'calc(50% + 18px)' }} />
            ) : null}
            <span className={cn('relative z-10 inline-flex h-9 w-9 items-center justify-center rounded-full border-2 transition-colors [&>svg]:h-4 [&>svg]:w-4', ring)}>
              {status === 'done' ? <Check /> : status === 'running' ? <Loader2 className="animate-spin" /> : status === 'failed' ? <X /> : status === 'skipped' ? <SkipForward /> : STAGE_ICONS[def.id]}
              {st.round && def.id === 'review' ? (
                <span className="absolute -right-1.5 -top-1.5 inline-flex h-4 min-w-4 items-center justify-center rounded-full bg-gold-500 px-1 text-[10px] font-semibold text-brand-900">{st.round}</span>
              ) : null}
            </span>
            <span className={cn('mt-2 text-xs font-medium', status === 'pending' ? 'text-slate-400 dark:text-slate-500' : 'text-slate-700 dark:text-slate-200')}>
              {def.index}. {def.name}
            </span>
            <span className="mt-0.5 text-[11px] tabular text-steel-500 dark:text-steel-400">{label}</span>
            {st.note && status !== 'pending' ? <span className="mt-0.5 line-clamp-1 max-w-full px-1 text-[11px] text-steel-500/90 dark:text-steel-400/90">{st.note}</span> : null}
          </li>
        );
      })}
    </ol>
  );
}
