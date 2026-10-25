import type { HTMLAttributes, ReactNode } from 'react';
import { cn } from '../../lib/cn';

export type BadgeTone = 'neutral' | 'brand' | 'gold' | 'green' | 'red' | 'amber' | 'blue' | 'purple' | 'slate';

const TONES: Record<BadgeTone, string> = {
  neutral: 'bg-slate-100 text-slate-700 ring-slate-200 dark:bg-slate-800 dark:text-slate-200 dark:ring-slate-700',
  slate: 'bg-steel-100 text-steel-700 ring-steel-200 dark:bg-steel-800/60 dark:text-steel-200 dark:ring-steel-700',
  brand: 'bg-brand-50 text-brand-700 ring-brand-100 dark:bg-brand-900/50 dark:text-brand-100 dark:ring-brand-800',
  gold: 'bg-gold-50 text-gold-800 ring-gold-200 dark:bg-gold-900/30 dark:text-gold-200 dark:ring-gold-800',
  green: 'bg-emerald-50 text-emerald-700 ring-emerald-200 dark:bg-emerald-900/30 dark:text-emerald-200 dark:ring-emerald-800',
  red: 'bg-red-50 text-red-700 ring-red-200 dark:bg-red-900/30 dark:text-red-200 dark:ring-red-800',
  amber: 'bg-amber-50 text-amber-800 ring-amber-200 dark:bg-amber-900/30 dark:text-amber-200 dark:ring-amber-800',
  blue: 'bg-sky-50 text-sky-700 ring-sky-200 dark:bg-sky-900/30 dark:text-sky-200 dark:ring-sky-800',
  purple: 'bg-violet-50 text-violet-700 ring-violet-200 dark:bg-violet-900/30 dark:text-violet-200 dark:ring-violet-800',
};

export function Badge({
  tone = 'neutral',
  className,
  children,
  dot,
  ...rest
}: HTMLAttributes<HTMLSpanElement> & { tone?: BadgeTone; dot?: boolean | 'pulse'; children: ReactNode }) {
  return (
    <span className={cn('inline-flex items-center gap-1.5 whitespace-nowrap rounded-full px-2 py-0.5 text-xs font-medium ring-1 ring-inset', TONES[tone], className)} {...rest}>
      {dot ? <span className={cn('h-1.5 w-1.5 rounded-full bg-current', dot === 'pulse' && 'animate-pulse-dot')} /> : null}
      {children}
    </span>
  );
}
