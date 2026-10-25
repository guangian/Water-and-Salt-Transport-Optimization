import { useEffect, useState } from 'react';
import { Bot, Database, FlaskConical, RefreshCw, ShieldAlert } from 'lucide-react';
import type { ModelRole, Settings, Track } from '../../types';
import type { Draft, Errors } from './draft';
import { Field, Segmented, Select, Slider } from '../ui/Form';
import { Button } from '../ui/Button';
import { api } from '../../api';
import { cn } from '../../lib/cn';

export const ROLE_LABELS: { role: ModelRole; label: string; desc: string }[] = [
  { role: 'main', label: '主代理 / 主笔 / 路演', desc: '长文档写作与流程调度' },
  { role: 'researcher', label: '行业研究员', desc: '检索归纳与长上下文' },
  { role: 'finance', label: '财务建模师', desc: '数值与 JSON 严谨' },
  { role: 'judgeA', label: '评委 A（投资人）', desc: '多模型评委团' },
  { role: 'judgeB', label: '评委 B（学者）', desc: '多模型评委团' },
  { role: 'judgeC', label: '评委 C（企业高管）', desc: '多模型评委团' },
];

export function useModelOptions(settings: Settings | null) {
  const [models, setModels] = useState<string[]>([]);
  const [loading, setLoading] = useState(false);
  const [source, setSource] = useState<'gateway' | 'settings' | null>(null);

  const refresh = async () => {
    setLoading(true);
    try {
      const r = await api.testSettings({});
      if (r.ok && r.models.length) {
        setModels(r.models);
        setSource('gateway');
        return;
      }
      throw new Error(r.error ?? '网关未返回模型');
    } catch {
      const fromSettings = settings ? Array.from(new Set(Object.values(settings.api.models).filter(Boolean))) : [];
      setModels(fromSettings);
      setSource('settings');
    } finally {
      setLoading(false);
    }
  };

  useEffect(() => {
    void refresh();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [settings?.api.baseUrl]);

  return { models, loading, source, refresh };
}

export function StepParams({ track, draft, onChange, errors, settings }: { track: Track | null; draft: Draft; onChange: (patch: Partial<Draft>) => void; errors: Errors; settings: Settings | null }) {
  const { models, loading, source, refresh } = useModelOptions(settings);
  const setParams = (patch: Partial<Draft['params']>) => onChange({ params: { ...draft.params, ...patch } });
  const pageLimit = track?.pageLimit ?? 80;
  const options = Array.from(new Set([...models, ...Object.values(draft.params.models).filter((x): x is string => Boolean(x))]));

  return (
    <div className="grid grid-cols-1 gap-6 xl:grid-cols-3">
      <section className="card p-5 xl:col-span-2">
        <div className="flex flex-wrap items-start justify-between gap-3">
          <div>
            <h2 className="flex items-center gap-2 text-base font-semibold text-slate-800 dark:text-slate-100">
              <Bot className="h-4 w-4 text-brand-500" />
              角色模型
            </h2>
            <p className="mt-1 text-xs text-steel-500 dark:text-steel-400">
              {loading ? '正在从网关获取模型列表…' : source === 'gateway' ? `模型列表来自网关 /v1/models（${models.length} 个）` : source === 'settings' ? '网关不可达，使用设置中的模型名' : ''}
            </p>
          </div>
          <Button variant="secondary" size="sm" icon={<RefreshCw className={cn(loading && 'animate-spin')} />} onClick={() => void refresh()} loading={loading}>
            刷新列表
          </Button>
        </div>
        {errors.models ? <p className="mt-2 text-xs text-red-600">{errors.models}</p> : null}
        <div className="mt-4 grid grid-cols-1 gap-4 md:grid-cols-2">
          {ROLE_LABELS.map(({ role, label, desc }) => (
            <Field key={role} label={label} hint={desc}>
              <Select value={draft.params.models[role] ?? ''} onChange={(e) => setParams({ models: { ...draft.params.models, [role]: e.target.value } })}>
                <option value="">（使用设置默认）</option>
                {options.map((m) => (
                  <option key={m} value={m}>
                    {m}
                  </option>
                ))}
              </Select>
            </Field>
          ))}
        </div>
      </section>

      <section className="card p-5">
        <h2 className="flex items-center gap-2 text-base font-semibold text-slate-800 dark:text-slate-100">
          <Database className="h-4 w-4 text-brand-500" />
          调研模式
        </h2>
        <div className="mt-4 space-y-3">
          {(
            [
              { value: 'synthetic', title: '示例数据', icon: <FlaskConical className="h-4 w-4" />, desc: '生成问卷/访谈提纲/执行计划，并用明确标注的模拟结果占位；提交前须替换为真实调研结果（系统生成《数据真实性清单》）。' },
              { value: 'provided', title: '仅用我的资料', icon: <ShieldAlert className="h-4 w-4" />, desc: '只引用上传资料与公开检索数据，不生成任何模拟调研结果；缺少一手数据的章节将标注【待补充】。' },
            ] as const
          ).map((o) => (
            <button
              type="button"
              key={o.value}
              onClick={() => setParams({ researchMode: o.value })}
              className={cn(
                'flex w-full items-start gap-3 rounded-xl border p-3.5 text-left transition',
                draft.params.researchMode === o.value ? 'border-gold-400 bg-gold-50/60 ring-2 ring-gold-200 dark:border-gold-500 dark:bg-gold-900/10 dark:ring-gold-700/40' : 'border-slate-200 hover:border-brand-300 dark:border-surface-dark-border',
              )}
            >
              <span className={cn('mt-0.5 inline-flex h-7 w-7 shrink-0 items-center justify-center rounded-lg', draft.params.researchMode === o.value ? 'bg-gold-500 text-brand-900' : 'bg-slate-100 text-slate-500 dark:bg-slate-800')}>{o.icon}</span>
              <span>
                <span className="block text-sm font-medium text-slate-800 dark:text-slate-100">{o.title}</span>
                <span className="mt-0.5 block text-xs leading-relaxed text-steel-500 dark:text-steel-400">{o.desc}</span>
              </span>
            </button>
          ))}
        </div>
      </section>

      <section className="card p-5 xl:col-span-3">
        <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">生成参数</h2>
        <div className="mt-5 grid grid-cols-1 gap-x-10 gap-y-6 lg:grid-cols-2">
          <Field label="正文目标页数" error={errors.pageTarget ?? null} hint={`赛道上限 ${pageLimit} 页（不含封面目录）；目标页数影响各章篇幅配额。`}>
            <Slider value={draft.params.pageTarget} min={30} max={80} step={5} onChange={(v) => setParams({ pageTarget: v })} format={(v) => `${v} 页`} marks={[30, 40, 50, 60, 70, 80]} />
          </Field>
          <Field label="讲稿时长" hint="决定路演页数与讲稿字数（约 260 字/分钟）。">
            <Segmented value={draft.params.pitchMinutes} onChange={(v) => setParams({ pitchMinutes: v })} options={[{ value: 8, label: '8 分钟' }, { value: 10, label: '10 分钟' }, { value: 15, label: '15 分钟' }]} />
          </Field>
          <Field label="评审轮数" hint="每轮三位评委独立打分 + 机器检查；未达目标分自动修改并重评，提升 <2 分视为收敛。">
            <Slider value={draft.params.reviewRounds} min={1} max={5} onChange={(v) => setParams({ reviewRounds: v })} format={(v) => `${v} 轮`} marks={[1, 2, 3, 4, 5]} />
          </Field>
          <Field label="目标分" hint="国金线默认 85 分；每个维度不低于 3/5。">
            <Slider value={draft.params.targetScore} min={70} max={95} onChange={(v) => setParams({ targetScore: v })} format={(v) => `${v} 分`} marks={[70, 75, 80, 85, 90, 95]} />
          </Field>
        </div>
      </section>
    </div>
  );
}
