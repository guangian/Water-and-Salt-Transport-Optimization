import { Check } from 'lucide-react';
import { cn } from '../../lib/cn';

export function Stepper({ steps, current, onSelect, maxReachable }: { steps: { title: string; description: string }[]; current: number; onSelect?: (i: number) => void; maxReachable: number }) {
  return (
    <ol className="flex items-center gap-2">
      {steps.map((s, i) => {
        const done = i < current;
        const active = i === current;
        const reachable = i <= maxReachable;
        return (
          <li key={s.title} className="flex flex-1 items-center gap-2 last:flex-none">
            <button
              type="button"
              disabled={!reachable || !onSelect}
              onClick={() => onSelect?.(i)}
              className={cn('flex items-center gap-3 rounded-lg px-2 py-1.5 text-left transition', reachable && onSelect ? 'hover:bg-slate-100 dark:hover:bg-slate-800' : 'cursor-default')}
            >
              <span
                className={cn(
                  'inline-flex h-8 w-8 shrink-0 items-center justify-center rounded-full text-sm font-semibold ring-2 transition',
                  done ? 'bg-brand-600 text-white ring-brand-600' : active ? 'bg-gold-500 text-brand-900 ring-gold-500' : 'bg-white text-slate-400 ring-slate-300 dark:bg-surface-dark-card dark:ring-slate-600',
                )}
              >
                {done ? <Check className="h-4 w-4" /> : i + 1}
              </span>
              <span className="hidden sm:block">
                <span className={cn('block text-sm font-medium', active ? 'text-brand-700 dark:text-gold-300' : done ? 'text-slate-700 dark:text-slate-200' : 'text-slate-400')}>{s.title}</span>
                <span className="block text-[11px] text-steel-500 dark:text-steel-400">{s.description}</span>
              </span>
            </button>
            {i < steps.length - 1 ? <span className={cn('h-px flex-1 rounded', i < current ? 'bg-brand-400' : 'bg-slate-200 dark:bg-slate-700')} /> : null}
          </li>
        );
      })}
    </ol>
  );
}
