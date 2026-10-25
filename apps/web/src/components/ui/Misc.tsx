import type { ReactNode } from 'react';
import { Loader2 } from 'lucide-react';
import { cn } from '../../lib/cn';

export function ProgressRing({
  value,
  size = 56,
  stroke = 6,
  label,
  tone = 'auto',
  className,
}: {
  value: number | null | undefined;
  size?: number;
  stroke?: number;
  label?: ReactNode;
  tone?: 'auto' | 'brand' | 'gold' | 'green' | 'red' | 'amber';
  className?: string;
}) {
  const v = value == null || !Number.isFinite(value) ? 0 : Math.max(0, Math.min(100, value));
  const r = (size - stroke) / 2;
  const c = 2 * Math.PI * r;
  const offset = c - (v / 100) * c;
  const color =
    tone === 'auto' ? (value == null ? '#94A3B8' : v >= 85 ? '#2E8B6E' : v >= 75 ? '#C9A227' : v >= 60 ? '#D97706' : '#B5533C') : { brand: '#14305C', gold: '#C9A227', green: '#2E8B6E', red: '#B5533C', amber: '#D97706' }[tone];
  return (
    <div className={cn('relative inline-flex items-center justify-center', className)} style={{ width: size, height: size }}>
      <svg width={size} height={size} className="-rotate-90">
        <circle cx={size / 2} cy={size / 2} r={r} stroke="currentColor" strokeWidth={stroke} fill="none" className="text-slate-200 dark:text-slate-700" />
        <circle cx={size / 2} cy={size / 2} r={r} stroke={color} strokeWidth={stroke} fill="none" strokeLinecap="round" strokeDasharray={c} strokeDashoffset={offset} className="transition-all duration-700" />
      </svg>
      <div className="absolute inset-0 flex items-center justify-center text-center leading-none">{label}</div>
    </div>
  );
}

export function Spinner({ className, label }: { className?: string; label?: string }) {
  return (
    <div className={cn('flex items-center justify-center gap-2 text-sm text-steel-500', className)}>
      <Loader2 className="h-4 w-4 animate-spin" />
      {label}
    </div>
  );
}

export function Skeleton({ className }: { className?: string }) {
  return <div className={cn('animate-pulse rounded-md bg-slate-200/80 dark:bg-slate-700/60', className)} />;
}

export function EmptyState({
  icon,
  title,
  description,
  action,
  className,
  compact,
}: {
  icon?: ReactNode;
  title: ReactNode;
  description?: ReactNode;
  action?: ReactNode;
  className?: string;
  compact?: boolean;
}) {
  return (
    <div className={cn('flex flex-col items-center justify-center text-center', compact ? 'px-4 py-8' : 'px-6 py-16', className)}>
      {icon ? (
        <div className={cn('mb-4 inline-flex items-center justify-center rounded-2xl bg-gradient-to-br from-brand-50 to-gold-50 text-brand-500 ring-1 ring-brand-100 dark:from-brand-900/60 dark:to-gold-900/20 dark:text-brand-200 dark:ring-brand-800', compact ? 'h-12 w-12 [&>svg]:h-6 [&>svg]:w-6' : 'h-16 w-16 [&>svg]:h-8 [&>svg]:w-8')}>
          {icon}
        </div>
      ) : null}
      <h3 className={cn('font-semibold text-slate-800 dark:text-slate-100', compact ? 'text-sm' : 'text-base')}>{title}</h3>
      {description ? <p className={cn('mt-1.5 max-w-md text-steel-500 dark:text-steel-400', compact ? 'text-xs' : 'text-sm')}>{description}</p> : null}
      {action ? <div className="mt-5">{action}</div> : null}
    </div>
  );
}

export function Stat({ label, value, sub, icon, tone = 'brand' }: { label: ReactNode; value: ReactNode; sub?: ReactNode; icon?: ReactNode; tone?: 'brand' | 'gold' | 'green' | 'steel' }) {
  const tones = {
    brand: 'bg-brand-50 text-brand-600 dark:bg-brand-900/50 dark:text-brand-200',
    gold: 'bg-gold-50 text-gold-700 dark:bg-gold-900/30 dark:text-gold-300',
    green: 'bg-emerald-50 text-emerald-600 dark:bg-emerald-900/30 dark:text-emerald-300',
    steel: 'bg-steel-100 text-steel-600 dark:bg-steel-800/60 dark:text-steel-200',
  };
  return (
    <div className="card flex items-center gap-4 px-5 py-4">
      {icon ? <div className={cn('inline-flex h-11 w-11 shrink-0 items-center justify-center rounded-xl [&>svg]:h-5 [&>svg]:w-5', tones[tone])}>{icon}</div> : null}
      <div className="min-w-0">
        <p className="text-xs font-medium text-steel-500 dark:text-steel-400">{label}</p>
        <div className="mt-0.5 truncate text-2xl font-semibold tabular text-slate-800 dark:text-slate-100">{value}</div>
        {sub ? <div className="mt-0.5 truncate text-xs text-steel-500 dark:text-steel-400">{sub}</div> : null}
      </div>
    </div>
  );
}

export function Tabs<T extends string>({
  value,
  onChange,
  tabs,
  className,
}: {
  value: T;
  onChange: (v: T) => void;
  tabs: { value: T; label: ReactNode; icon?: ReactNode; badge?: ReactNode }[];
  className?: string;
}) {
  return (
    <div className={cn('flex gap-1 overflow-x-auto border-b border-slate-200 scrollbar-none dark:border-surface-dark-border', className)} role="tablist">
      {tabs.map((t) => {
        const active = t.value === value;
        return (
          <button
            key={t.value}
            role="tab"
            aria-selected={active}
            onClick={() => onChange(t.value)}
            className={cn(
              '-mb-px inline-flex shrink-0 items-center gap-2 border-b-2 px-3.5 py-2.5 text-sm font-medium transition',
              active ? 'border-gold-500 text-brand-700 dark:text-gold-300' : 'border-transparent text-slate-500 hover:border-slate-300 hover:text-slate-800 dark:text-slate-400 dark:hover:text-slate-100',
            )}
          >
            {t.icon ? <span className="[&>svg]:h-4 [&>svg]:w-4">{t.icon}</span> : null}
            {t.label}
            {t.badge}
          </button>
        );
      })}
    </div>
  );
}

export function KeyValue({ items, className, cols = 2 }: { items: { k: ReactNode; v: ReactNode }[]; className?: string; cols?: 1 | 2 | 3 }) {
  return (
    <dl className={cn('grid gap-x-6 gap-y-3 text-sm', cols === 1 ? 'grid-cols-1' : cols === 2 ? 'grid-cols-1 sm:grid-cols-2' : 'grid-cols-1 sm:grid-cols-3', className)}>
      {items.map((it, i) => (
        <div key={i} className="min-w-0">
          <dt className="text-xs text-steel-500 dark:text-steel-400">{it.k}</dt>
          <dd className="mt-0.5 break-words text-slate-800 dark:text-slate-100">{it.v}</dd>
        </div>
      ))}
    </dl>
  );
}
