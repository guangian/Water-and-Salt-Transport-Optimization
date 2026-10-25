import { Check, FileText } from 'lucide-react';
import type { Track } from '../../types';
import type { Draft, Errors } from './draft';
import { TRACK_META } from '../domain/Brand';
import { Field, Select } from '../ui/Form';
import { Skeleton } from '../ui/Misc';
import { cn } from '../../lib/cn';

export function StepTrack({ tracks, draft, onChange, errors }: { tracks: Track[] | null; draft: Draft; onChange: (patch: Partial<Draft>) => void; errors: Errors }) {
  const selected = tracks?.find((t) => t.id === draft.track) ?? null;
  return (
    <div className="space-y-6">
      <div>
        <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">选择赛道</h2>
        <p className="mt-1 text-sm text-steel-500 dark:text-steel-400">每个赛道有独立的章节大纲、评分量规与机器检查阈值。卡片展示官方评分权重与必备要素。</p>
        {errors.track ? <p className="mt-1 text-xs text-red-600">{errors.track}</p> : null}
      </div>
      {!tracks ? (
        <div className="grid grid-cols-1 gap-4 lg:grid-cols-2 2xl:grid-cols-3">
          {Array.from({ length: 5 }).map((_, i) => (
            <Skeleton key={i} className="h-56" />
          ))}
        </div>
      ) : (
        <div className="grid grid-cols-1 gap-4 lg:grid-cols-2 2xl:grid-cols-3">
          {tracks.map((t) => {
            const meta = TRACK_META[t.id];
            const active = draft.track === t.id;
            const total = t.scoring.written + t.scoring.pitch;
            return (
              <button
                type="button"
                key={t.id}
                onClick={() => {
                  const subjectType = t.id === 'startup' ? 'startup' : t.id === 'case' ? 'case-company' : t.id === 'mice' ? 'event' : 'brand';
                  onChange({
                    track: t.id,
                    group: t.groups.includes(draft.group) ? draft.group : t.groups[0] ?? '',
                    subject: { ...draft.subject, type: subjectType, isAssignedCase: t.id === 'case' ? true : draft.subject.isAssignedCase, caseTitle: '' },
                    params: { ...draft.params, pageTarget: Math.min(draft.params.pageTarget, t.pageLimit) },
                  });
                }}
                className={cn(
                  'card card-hover relative flex flex-col p-5 text-left transition',
                  active && 'border-gold-400 ring-2 ring-gold-300/60 dark:border-gold-500 dark:ring-gold-500/40',
                )}
              >
                {active ? (
                  <span className="absolute right-4 top-4 inline-flex h-6 w-6 items-center justify-center rounded-full bg-gold-500 text-brand-900">
                    <Check className="h-4 w-4" />
                  </span>
                ) : null}
                <div className="flex items-center gap-3">
                  <span className={cn('inline-flex h-10 w-10 items-center justify-center rounded-xl text-sm font-bold', meta.bg, meta.text)}>{meta.short}</span>
                  <div className="min-w-0 pr-8">
                    <h3 className="truncate text-[15px] font-semibold text-slate-800 dark:text-slate-100">{t.name}</h3>
                    <p className="mt-0.5 flex items-center gap-1 text-xs text-steel-500 dark:text-steel-400">
                      <FileText className="h-3.5 w-3.5" />
                      {t.docName} · ≤{t.pageLimit} 页{t.wordMin ? ` · ≥${t.wordMin} 字` : ''}
                    </p>
                  </div>
                </div>
                <p className="mt-3 line-clamp-2 text-xs leading-relaxed text-slate-600 dark:text-slate-300">{t.description}</p>

                <div className="mt-4">
                  <div className="mb-1 flex items-center justify-between text-[11px] text-steel-500">
                    <span>评分权重</span>
                    <span className="tabular">
                      书面 {t.scoring.written} / 展示答辩 {t.scoring.pitch}
                    </span>
                  </div>
                  <div className="flex h-2 overflow-hidden rounded-full bg-slate-100 dark:bg-slate-800">
                    <span className="h-full bg-brand-600 dark:bg-brand-400" style={{ width: `${(t.scoring.written / total) * 100}%` }} />
                    <span className="h-full bg-gold-500" style={{ width: `${(t.scoring.pitch / total) * 100}%` }} />
                  </div>
                </div>

                <div className="mt-4">
                  <div className="mb-1.5 text-[11px] text-steel-500">必备要素（{t.elements.length}）</div>
                  <div className="flex flex-wrap gap-1">
                    {t.elements.slice(0, 7).map((el) => (
                      <span key={el} className="chip">
                        {el}
                      </span>
                    ))}
                    {t.elements.length > 7 ? <span className="chip">+{t.elements.length - 7}</span> : null}
                  </div>
                </div>
              </button>
            );
          })}
        </div>
      )}

      {selected ? (
        <div className="card grid grid-cols-1 gap-5 p-5 lg:grid-cols-3">
          <Field label="组别" required error={errors.group ?? null} hint="组别将写入封面与任务书；命题组须选择命题方向。">
            <Select value={draft.group} onChange={(e) => onChange({ group: e.target.value })}>
              {selected.groups.map((g) => (
                <option key={g} value={g}>
                  {g}
                </option>
              ))}
            </Select>
          </Field>
          <div className="lg:col-span-2">
            <div className="label">交付物清单</div>
            <ul className="grid grid-cols-1 gap-1.5 text-sm text-slate-700 sm:grid-cols-2 dark:text-slate-300">
              {(selected.deliverables ?? []).map((d) => (
                <li key={d} className="flex items-start gap-2">
                  <Check className="mt-0.5 h-4 w-4 shrink-0 text-emerald-500" />
                  {d}
                </li>
              ))}
            </ul>
          </div>
        </div>
      ) : null}
    </div>
  );
}
