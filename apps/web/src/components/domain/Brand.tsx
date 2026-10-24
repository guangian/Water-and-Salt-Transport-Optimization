import type { ProjectStatus, StageStatus, TrackId } from '../../types';
import { Badge, type BadgeTone } from '../ui/Badge';
import { cn } from '../../lib/cn';

export const TRACK_META: Record<TrackId, { short: string; tone: BadgeTone; color: string; bg: string; text: string }> = {
  brand: { short: '品牌', tone: 'brand', color: '#14305C', bg: 'bg-brand-600', text: 'text-white' },
  startup: { short: '创业', tone: 'gold', color: '#C9A227', bg: 'bg-gold-500', text: 'text-brand-900' },
  digital: { short: '数营', tone: 'purple', color: '#6D4FC2', bg: 'bg-violet-600', text: 'text-white' },
  case: { short: '案例', tone: 'green', color: '#2E8B6E', bg: 'bg-emerald-600', text: 'text-white' },
  mice: { short: '会奖', tone: 'blue', color: '#2C7FB8', bg: 'bg-sky-600', text: 'text-white' },
};

export const TRACK_NAMES: Record<TrackId, string> = {
  brand: '品牌策划赛道',
  startup: '创业计划赛道',
  digital: '数字营销赛道',
  case: '会计与商业管理案例',
  mice: '商务会奖旅游策划',
};

export function TrackBadge({ track, size = 'md', withName = false, className }: { track: TrackId; size?: 'sm' | 'md' | 'lg'; withName?: boolean; className?: string }) {
  const m = TRACK_META[track] ?? TRACK_META.brand;
  return (
    <span className={cn('inline-flex items-center gap-2', className)}>
      <span
        className={cn(
          'inline-flex shrink-0 items-center justify-center rounded-lg font-semibold tracking-wide',
          m.bg,
          m.text,
          size === 'sm' ? 'h-6 min-w-6 px-1.5 text-[11px]' : size === 'lg' ? 'h-10 min-w-10 px-2 text-sm' : 'h-8 min-w-8 px-2 text-xs',
        )}
      >
        {m.short}
      </span>
      {withName ? <span className="text-sm text-slate-700 dark:text-slate-200">{TRACK_NAMES[track]}</span> : null}
    </span>
  );
}

export const STATUS_META: Record<ProjectStatus, { label: string; tone: BadgeTone; dot?: boolean | 'pulse' }> = {
  idle: { label: '待运行', tone: 'neutral', dot: true },
  running: { label: '运行中', tone: 'green', dot: 'pulse' },
  paused: { label: '已暂停', tone: 'amber', dot: true },
  completed: { label: '已完成', tone: 'brand', dot: true },
  failed: { label: '失败', tone: 'red', dot: true },
  stopped: { label: '已停止', tone: 'slate', dot: true },
};

export function StatusPill({ status, className }: { status: ProjectStatus; className?: string }) {
  const m = STATUS_META[status] ?? STATUS_META.idle;
  return (
    <Badge tone={m.tone} dot={m.dot} className={className}>
      {m.label}
    </Badge>
  );
}

export const STAGE_STATUS_META: Record<StageStatus, { label: string; tone: BadgeTone }> = {
  pending: { label: '等待', tone: 'neutral' },
  running: { label: '进行中', tone: 'green' },
  done: { label: '完成', tone: 'brand' },
  failed: { label: '失败', tone: 'red' },
  skipped: { label: '跳过', tone: 'slate' },
};

export function Logo({ size = 36, className }: { size?: number; className?: string }) {
  return (
    <svg width={size} height={size} viewBox="0 0 64 64" className={className} aria-hidden="true">
      <defs>
        <linearGradient id="bg-logo-gold" x1="0" y1="0" x2="1" y2="1">
          <stop offset="0" stopColor="#E4C65A" />
          <stop offset="1" stopColor="#A8861F" />
        </linearGradient>
      </defs>
      <rect x="4" y="4" width="56" height="56" rx="14" fill="#14305C" />
      <path d="M32 12 L50 24 L50 42 L32 54 L14 42 L14 24 Z" fill="none" stroke="url(#bg-logo-gold)" strokeWidth="3.5" strokeLinejoin="round" />
      <path d="M22 40 L30 30 L36 36 L44 24" fill="none" stroke="url(#bg-logo-gold)" strokeWidth="4" strokeLinecap="round" strokeLinejoin="round" />
      <circle cx="44" cy="24" r="3.2" fill="#E4C65A" />
    </svg>
  );
}

export function ScoreText({ score, target, className }: { score: number | null | undefined; target?: number; className?: string }) {
  if (score == null) return <span className={cn('text-steel-400', className)}>—</span>;
  const color = target != null ? (score >= target ? 'text-emerald-600 dark:text-emerald-400' : score >= target - 10 ? 'text-gold-600 dark:text-gold-300' : 'text-red-600 dark:text-red-400') : 'text-slate-800 dark:text-slate-100';
  return <span className={cn('tabular font-semibold', color, className)}>{score.toFixed(1)}</span>;
}
