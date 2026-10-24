import { useEffect, useMemo, useState, type ReactNode } from 'react';
import { AlertTriangle, CheckCircle2, Gavel, MinusCircle, ShieldAlert, ShieldCheck, Target, TrendingUp, XCircle } from 'lucide-react';
import {
  Bar,
  BarChart,
  CartesianGrid,
  Cell,
  Legend,
  Line,
  LineChart,
  PolarAngleAxis,
  PolarGrid,
  PolarRadiusAxis,
  Radar,
  RadarChart,
  ReferenceLine,
  ResponsiveContainer,
  Tooltip,
  XAxis,
  YAxis,
} from 'recharts';
import type { CheckStatus, CrossDocMismatch, HardGate, JudgeId, JudgeReview, MachineCheckReport, Project, ReviewData, ReviewRound } from '../../types';
import { api } from '../../api';
import { useAsync } from '../../hooks/useAsync';
import { formatDateTime, formatNumber, formatScore } from '../../lib/format';
import { cn } from '../../lib/cn';
import { Card, CardBody, CardHeader } from '../ui/Card';
import { Badge } from '../ui/Badge';
import { EmptyState, ProgressRing, Spinner } from '../ui/Misc';
import { Segmented } from '../ui/Form';
import { useTheme, type Theme } from '../../hooks/useTheme';

const JUDGES: Record<JudgeId, { name: string; short: string; persona: string; color: string; colorDark: string }> = {
  investor: { name: '评委 A · 投资人', short: '投资人', persona: '商业逻辑、财务可信度与增长路径', color: '#14305C', colorDark: '#8FB3E8' },
  professor: { name: '评委 B · 学者', short: '学者', persona: '理论框架、调研方法与论证严谨', color: '#C9A227', colorDark: '#E3C25A' },
  executive: { name: '评委 C · 企业高管', short: '企业高管', persona: '落地可行性、预算与执行细节', color: '#2E8B6E', colorDark: '#5FC3A0' },
};
const JUDGE_ORDER: JudgeId[] = ['investor', 'professor', 'executive'];
// 图表里的品牌藏青在深色底上不可见，按主题切换
const NAVY: Record<Theme, string> = { light: '#14305C', dark: '#8FB3E8' };
const judgeColor = (id: JudgeId, theme: Theme): string => (theme === 'dark' ? JUDGES[id].colorDark : JUDGES[id].color);

const GATE_LABELS: Record<string, string> = {
  anonymity: '匿名要求',
  anonymity_hits: '匿名要求',
  page_limit: '页数上限',
  required_elements: '必备要素',
  forbidden_terms: '禁语检查',
  cover_fields: '封面字段',
  file_naming: '文件命名',
};

const CHECK_META: Record<CheckStatus, { label: string; icon: ReactNode; tone: 'green' | 'amber' | 'red' | 'slate' }> = {
  pass: { label: '通过', icon: <CheckCircle2 className="h-4 w-4 text-emerald-500" />, tone: 'green' },
  warn: { label: '警告', icon: <AlertTriangle className="h-4 w-4 text-amber-500" />, tone: 'amber' },
  fail: { label: '失败', icon: <XCircle className="h-4 w-4 text-red-500" />, tone: 'red' },
  skipped: { label: '跳过', icon: <MinusCircle className="h-4 w-4 text-slate-400" />, tone: 'slate' },
};

function ChartTooltip({ active, payload, label, suffix = '' }: { active?: boolean; payload?: { name: string; value: number | string; color?: string }[]; label?: string; suffix?: string }) {
  if (!active || !payload?.length) return null;
  return (
    <div className="rounded-lg border border-slate-200 bg-white/95 px-3 py-2 text-xs shadow-card dark:border-surface-dark-border dark:bg-surface-dark-card/95">
      {label ? <div className="mb-1 font-medium text-slate-800 dark:text-slate-100">{label}</div> : null}
      {payload.map((p, i) => (
        <div key={i} className="flex items-center gap-2 text-steel-600 dark:text-steel-300">
          <span className="h-2 w-2 rounded-full" style={{ background: p.color ?? '#5B6B82' }} />
          {p.name}
          <span className="ml-auto pl-3 font-semibold tabular text-slate-800 dark:text-slate-100">
            {typeof p.value === 'number' ? p.value.toFixed(1) : p.value}
            {suffix}
          </span>
        </div>
      ))}
    </div>
  );
}

function gateLabel(g: HardGate): string {
  const id = g.id ?? g.gate ?? '';
  return GATE_LABELS[id] ?? id;
}

function gateDetail(g: HardGate): string {
  if (g.detail) return g.detail;
  if (g.pages != null && g.limit != null) return `${g.pages} / ${g.limit} 页`;
  if (g.hits) return g.hits.length ? `命中 ${g.hits.length} 处：${g.hits.join('、')}` : '0 处命中';
  return g.passed ? '通过' : '未通过';
}

function JudgeCard({ review }: { review: JudgeReview }) {
  const meta = JUDGES[review.judge];
  const failedGates = review.hardGates.filter((g) => !g.passed);
  return (
    <Card className="flex flex-col">
      <div className="flex items-start gap-3 border-b border-slate-100 px-4 py-3 dark:border-surface-dark-border">
        <span className="mt-0.5 inline-flex h-9 w-9 shrink-0 items-center justify-center rounded-lg text-white" style={{ background: meta.color }}>
          <Gavel className="h-4 w-4" />
        </span>
        <div className="min-w-0 flex-1">
          <h4 className="truncate text-sm font-semibold text-slate-800 dark:text-slate-100">{meta.name}</h4>
          <p className="mt-1 flex items-center gap-1.5 text-[11px] text-steel-500">
            <span className="chip shrink-0 font-mono text-[10px]">{review.model}</span>
            <span className="truncate">{meta.persona}</span>
          </p>
        </div>
        <div className="text-right">
          <div className="text-2xl font-semibold tabular text-brand-700 dark:text-gold-200">{formatScore(review.total)}</div>
          {review.officialEstimate ? <div className="text-[10px] text-steel-500">书面 {review.officialEstimate.written} · 路演 {review.officialEstimate.pitch}</div> : null}
        </div>
      </div>
      <CardBody className="flex flex-1 flex-col gap-3 px-4 py-3">
        <p className="text-[13px] leading-relaxed text-slate-700 dark:text-slate-200">{review.verdict}</p>
        {failedGates.length ? (
          <div className="flex flex-wrap gap-1">
            {failedGates.map((g, i) => (
              <Badge key={i} tone="red">
                {gateLabel(g)}未过
              </Badge>
            ))}
          </div>
        ) : (
          <Badge tone="green" className="w-fit">
            硬门槛 {review.hardGates.length}/{review.hardGates.length} 通过
          </Badge>
        )}
        <div>
          <div className="mb-1.5 text-xs font-medium text-steel-500">优先修改（Top 3）</div>
          <ol className="space-y-1.5">
            {review.topFixes.slice(0, 3).map((f) => (
              <li key={f.priority} className="flex items-start gap-2 text-xs">
                <span className="mt-px inline-flex h-4 w-4 shrink-0 items-center justify-center rounded bg-gold-100 text-[10px] font-semibold text-gold-800 dark:bg-gold-900/40 dark:text-gold-200">{f.priority}</span>
                <span className="min-w-0 flex-1">
                  <span className="text-slate-700 dark:text-slate-200">{f.what}</span>
                  <span className="ml-1 text-steel-500">@ {f.where}</span>
                </span>
                <span className="shrink-0 font-semibold tabular text-emerald-600 dark:text-emerald-300">{f.gain}</span>
              </li>
            ))}
            {review.topFixes.length === 0 ? <li className="text-xs text-steel-500">无需修改</li> : null}
          </ol>
        </div>
      </CardBody>
    </Card>
  );
}

function MachineChecks({ report }: { report: MachineCheckReport }) {
  const [filter, setFilter] = useState<'all' | CheckStatus>('all');
  const list = report.checks.filter((c) => filter === 'all' || c.status === filter);
  const summary = report.summary;
  return (
    <Card>
      <CardHeader
        title="机器检查"
        subtitle={`${report.checks.length} 项 · 生成于 ${formatDateTime(report.generatedAt)}`}
        actions={
          <Segmented
            size="sm"
            value={filter}
            onChange={setFilter}
            options={[
              { value: 'all', label: `全部 ${report.checks.length}` },
              { value: 'pass', label: `通过 ${summary.pass}` },
              { value: 'warn', label: `警告 ${summary.warn}` },
              { value: 'fail', label: `失败 ${summary.fail}` },
              { value: 'skipped', label: `跳过 ${summary.skipped ?? report.checks.filter((c) => c.status === 'skipped').length}` },
            ]}
          />
        }
      />
      <div className="max-h-[480px] overflow-auto">
        <table className="table-base">
          <thead>
            <tr>
              <th className="w-20">状态</th>
              <th>检查项</th>
              <th className="w-28 text-right">值</th>
              <th className="w-32">阈值</th>
              <th>说明</th>
            </tr>
          </thead>
          <tbody>
            {list.map((c) => {
              const m = CHECK_META[c.status];
              const isPct = typeof c.value === 'number' && c.value > 0 && c.value <= 1 && (c.threshold?.pass ?? 0) <= 1 && c.id !== 'forbidden_terms' && c.id !== 'anonymity_hits' && c.threshold != null;
              const fmt = (v: number | string | boolean | null | undefined) => (v == null ? '—' : isPct && typeof v === 'number' ? `${Math.round(v * 100)}%` : formatNumber(v));
              return (
                <tr key={c.id}>
                  <td>
                    <span className="inline-flex items-center gap-1.5 text-xs">
                      {m.icon}
                      {m.label}
                    </span>
                  </td>
                  <td>
                    <div className="text-[13px] text-slate-800 dark:text-slate-100">{c.name ?? c.id}</div>
                    <div className="font-mono text-[11px] text-steel-400">{c.id}</div>
                  </td>
                  <td className="text-right tabular">{fmt(c.value)}</td>
                  <td className="text-xs text-steel-500">{c.threshold ? `通过 ${fmt(c.threshold.pass)} · 警告 ${fmt(c.threshold.warn)}` : '—'}</td>
                  <td className="text-xs">{c.detail}</td>
                </tr>
              );
            })}
            {list.length === 0 ? (
              <tr>
                <td colSpan={5} className="py-8 text-center text-steel-500">
                  没有该状态的检查项
                </td>
              </tr>
            ) : null}
          </tbody>
        </table>
      </div>
    </Card>
  );
}

function CrossDocTable({ mismatches, consistent }: { mismatches: CrossDocMismatch[]; consistent: boolean }) {
  const cols: { key: keyof CrossDocMismatch; label: string }[] = [
    { key: 'book', label: '策划书' },
    { key: 'deck', label: '路演 PPT' },
    { key: 'script', label: '讲稿' },
    { key: 'qa', label: '题库' },
  ];
  return (
    <Card>
      <CardHeader title="跨材料一致性" subtitle="同一关键数字在策划书 / 路演 / 讲稿 / 题库中的取值对照" actions={<Badge tone={consistent ? 'green' : 'red'}>{consistent ? '一致' : `${mismatches.length} 处不一致`}</Badge>} />
      {mismatches.length === 0 ? (
        <EmptyState compact icon={<ShieldCheck />} title="全部关键数字一致" />
      ) : (
        <table className="table-base">
          <thead>
            <tr>
              <th>关键数字</th>
              {cols.map((c) => (
                <th key={c.key}>{c.label}</th>
              ))}
            </tr>
          </thead>
          <tbody>
            {mismatches.map((m) => {
              const values = cols.map((c) => m[c.key]).filter((v) => v != null && v !== '');
              const mode = values.reduce<Record<string, number>>((acc, v) => {
                acc[String(v)] = (acc[String(v)] ?? 0) + 1;
                return acc;
              }, {});
              const majority = Object.entries(mode).sort((a, b) => b[1] - a[1])[0]?.[0];
              return (
                <tr key={m.key}>
                  <td>
                    <div className="text-[13px] text-slate-800 dark:text-slate-100">{m.label ?? m.key}</div>
                    <div className="font-mono text-[11px] text-steel-400">{m.key}</div>
                  </td>
                  {cols.map((c) => {
                    const v = m[c.key];
                    const bad = v != null && v !== '' && String(v) !== majority;
                    return (
                      <td key={c.key} className={cn('tabular', bad && 'font-semibold text-red-600 dark:text-red-300')}>
                        {v == null || v === '' ? <span className="text-steel-400">—</span> : formatNumber(v)}
                      </td>
                    );
                  })}
                </tr>
              );
            })}
          </tbody>
        </table>
      )}
    </Card>
  );
}

function RoundView({ round, target }: { round: ReviewRound; target: number }) {
  const [theme] = useTheme();
  const judges = JUDGE_ORDER.map((id) => round.judges.find((j) => j.judge === id)).filter((j): j is JudgeReview => Boolean(j));
  const radarData = useMemo(() => {
    const dims = judges[0]?.scores.map((s) => s.dimension) ?? [];
    return dims.map((d) => {
      const row: Record<string, string | number> = { dimension: d };
      for (const j of judges) {
        const s = j.scores.find((x) => x.dimension === d);
        row[JUDGES[j.judge].short] = s ? Math.round((s.score / s.cap) * 1000) / 10 : 0;
      }
      return row;
    });
  }, [judges]);
  const barData = useMemo(() => {
    const dims = judges[0]?.scores ?? [];
    return dims.map((d) => {
      const got = judges.reduce((s, j) => {
        const x = j.scores.find((y) => y.dimension === d.dimension);
        return s + (x ? (x.score / x.cap) * x.weight : 0);
      }, 0);
      return { dimension: d.dimension, 得分: Math.round((got / Math.max(1, judges.length)) * 10) / 10, 满分: d.weight };
    });
  }, [judges]);
  const machine = round.machine;
  const gates = machine?.hardGates ?? judges[0]?.hardGates ?? [];

  return (
    <div className="space-y-4">
      <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_minmax(0,1fr)]">
        <Card>
          <CardHeader title="维度雷达" subtitle="三位评委逐维得分（折算为百分制）" />
          <CardBody className="text-slate-600 dark:text-slate-300">
            <div className="h-[340px]">
              <ResponsiveContainer width="100%" height="100%">
                <RadarChart data={radarData} outerRadius="72%">
                  <PolarGrid stroke="currentColor" strokeOpacity={0.18} />
                  <PolarAngleAxis dataKey="dimension" tick={{ fill: 'currentColor', fontSize: 11 }} />
                  <PolarRadiusAxis domain={[0, 100]} tick={{ fill: 'currentColor', fontSize: 10 }} tickCount={5} axisLine={false} />
                  {judges.map((j) => (
                    <Radar key={j.judge} name={JUDGES[j.judge].short} dataKey={JUDGES[j.judge].short} stroke={judgeColor(j.judge, theme)} fill={judgeColor(j.judge, theme)} fillOpacity={0.12} strokeWidth={2} />
                  ))}
                  <Legend wrapperStyle={{ fontSize: 12 }} />
                  <Tooltip content={<ChartTooltip />} />
                </RadarChart>
              </ResponsiveContainer>
            </div>
          </CardBody>
        </Card>
        <Card>
          <CardHeader title="维度得分与权重" subtitle="三评委平均加权得分 / 该维度满分" />
          <CardBody className="text-slate-600 dark:text-slate-300">
            <div className="h-[340px]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={barData} layout="vertical" margin={{ left: 8, right: 24, top: 4, bottom: 4 }} barCategoryGap="30%">
                  <CartesianGrid horizontal={false} stroke="currentColor" strokeOpacity={0.12} />
                  <XAxis type="number" tick={{ fill: 'currentColor', fontSize: 11 }} axisLine={false} tickLine={false} />
                  <YAxis type="category" dataKey="dimension" width={112} tick={{ fill: 'currentColor', fontSize: 11 }} axisLine={false} tickLine={false} />
                  <Tooltip content={<ChartTooltip />} cursor={{ fill: 'currentColor', fillOpacity: 0.05 }} />
                  <Bar dataKey="满分" fill="currentColor" fillOpacity={0.12} radius={[0, 4, 4, 0]} barSize={14} />
                  <Bar dataKey="得分" radius={[0, 4, 4, 0]} barSize={14}>
                    {barData.map((d, i) => (
                      <Cell key={i} fill={d.得分 / d.满分 >= 0.9 ? judgeColor('executive', theme) : d.得分 / d.满分 >= 0.8 ? NAVY[theme] : judgeColor('professor', theme)} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>
          </CardBody>
        </Card>
      </div>

      <div className="grid gap-4 lg:grid-cols-3">
        {judges.map((j) => (
          <JudgeCard key={j.judge} review={j} />
        ))}
      </div>

      <div className="grid gap-4 xl:grid-cols-[360px_minmax(0,1fr)]">
        <Card>
          <CardHeader title="硬门槛" subtitle="任一项未通过都视为不合格" actions={<Badge tone={gates.every((g) => g.passed) ? 'green' : 'red'}>{gates.filter((g) => g.passed).length}/{gates.length} 通过</Badge>} />
          <ul className="divide-y divide-slate-100 dark:divide-surface-dark-border">
            {gates.map((g, i) => (
              <li key={i} className="flex items-start gap-3 px-5 py-3">
                {g.passed ? <ShieldCheck className="mt-0.5 h-4 w-4 shrink-0 text-emerald-500" /> : <ShieldAlert className="mt-0.5 h-4 w-4 shrink-0 text-red-500" />}
                <div className="min-w-0 flex-1">
                  <div className="text-sm font-medium text-slate-800 dark:text-slate-100">{gateLabel(g)}</div>
                  <div className={cn('mt-0.5 break-words text-xs', g.passed ? 'text-steel-500' : 'text-red-600 dark:text-red-300')}>{gateDetail(g)}</div>
                </div>
              </li>
            ))}
            {gates.length === 0 ? <li className="px-5 py-6 text-center text-sm text-steel-500">无硬门槛数据</li> : null}
          </ul>
        </Card>
        <div className="space-y-4">
          {machine ? <MachineChecks report={machine} /> : null}
          {machine ? <CrossDocTable mismatches={machine.crossDoc.mismatches} consistent={machine.crossDoc.consistent} /> : null}
        </div>
      </div>

      {round.fixes.length ? (
        <Card>
          <CardHeader title={`修改清单 · 第 ${round.round} 轮`} subtitle="评审汇总后的合并清单，按优先级执行" actions={<Badge tone="gold">目标 {target} 分</Badge>} />
          <table className="table-base">
            <thead>
              <tr>
                <th className="w-16">优先级</th>
                <th>位置</th>
                <th>修改内容</th>
                <th className="w-24 text-right">预计提升</th>
              </tr>
            </thead>
            <tbody>
              {round.fixes.map((f) => (
                <tr key={f.priority}>
                  <td className="tabular">{f.priority}</td>
                  <td className="font-mono text-xs">{f.where}</td>
                  <td>{f.what}</td>
                  <td className="text-right font-semibold tabular text-emerald-600 dark:text-emerald-300">{f.gain}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </Card>
      ) : null}
    </div>
  );
}

export function ReviewTab({ project, refreshKey }: { project: Project; refreshKey: number }) {
  const { data, loading, error } = useAsync<ReviewData>(() => api.getReview(project.id), [project.id, refreshKey]);
  const [theme] = useTheme();
  const [roundNo, setRoundNo] = useState<number | null>(null);
  const rounds = data?.rounds ?? [];
  const latest = data?.latest ?? rounds[rounds.length - 1] ?? null;
  useEffect(() => {
    if (latest) setRoundNo((r) => (r != null && rounds.some((x) => x.round === r) ? r : latest.round));
  }, [latest, rounds]);
  const round = rounds.find((r) => r.round === roundNo) ?? latest;
  const target = project.params.targetScore;

  if (loading && !data) return <Spinner className="py-20" label="加载评审数据…" />;
  if (error) return <EmptyState icon={<Gavel />} title="评审数据加载失败" description={error.message} />;
  if (!round) {
    return (
      <Card>
        <EmptyState
          icon={<Gavel />}
          title="尚未进行评审"
          description={project.status === 'running' ? '评审阶段开始后，三位评委（投资人 / 学者 / 企业高管）与质检员的结果将出现在这里。' : '运行到「评审」阶段后，这里会显示三评委打分、机器检查、硬门槛与修改清单。'}
        />
      </Card>
    );
  }

  const trend = rounds.map((r) => ({ name: `第 ${r.round} 轮`, 均分: r.score, 投资人: r.byJudge.investor ?? null, 学者: r.byJudge.professor ?? null, 企业高管: r.byJudge.executive ?? null }));
  const delta = rounds.length >= 2 ? round.score - (rounds.find((r) => r.round === round.round - 1)?.score ?? round.score) : null;

  return (
    <div className="space-y-4">
      <div className="grid gap-4 xl:grid-cols-[minmax(0,1.1fr)_minmax(0,1fr)]">
        <Card className="relative overflow-hidden">
          <div className="pointer-events-none absolute -right-16 -top-16 h-56 w-56 rounded-full bg-gold-100/60 blur-2xl dark:bg-gold-900/20" />
          <CardBody className="flex flex-wrap items-center gap-6 py-5">
            <ProgressRing
              value={round.score}
              size={132}
              stroke={10}
              label={
                <div>
                  <div className="text-3xl font-semibold tabular text-slate-800 dark:text-slate-100">{formatScore(round.score)}</div>
                  <div className="mt-0.5 text-[11px] text-steel-500">第 {round.round} 轮均分</div>
                </div>
              }
            />
            <div className="min-w-0 flex-1">
              <div className="flex flex-wrap items-center gap-2">
                <h3 className="text-lg font-semibold text-slate-800 dark:text-slate-100">{round.score >= target ? '已达到目标分' : `距目标还差 ${(target - round.score).toFixed(1)} 分`}</h3>
                <Badge tone={round.hardGatePassed ? 'green' : 'red'}>{round.hardGatePassed ? '硬门槛通过' : '硬门槛未通过'}</Badge>
                {round.converged ? <Badge tone="brand">已收敛</Badge> : null}
                {delta != null ? (
                  <Badge tone={delta >= 0 ? 'green' : 'red'}>
                    <TrendingUp className="h-3 w-3" />
                    较上轮 {delta >= 0 ? '+' : ''}
                    {delta.toFixed(1)}
                  </Badge>
                ) : null}
              </div>
              <div className="mt-2 flex items-center gap-2 text-xs text-steel-500">
                <Target className="h-3.5 w-3.5" />
                目标 {target} 分 · 国金线 85 分 · 共 {rounds.length} 轮
                {round.at ? ` · ${formatDateTime(round.at)}` : ''}
              </div>
              <div className="mt-4 grid grid-cols-3 gap-3">
                {JUDGE_ORDER.map((id) => (
                  <div key={id} className="rounded-lg border border-slate-200 px-3 py-2 dark:border-surface-dark-border">
                    <div className="flex items-center gap-1.5 text-[11px] text-steel-500">
                      <span className="h-2 w-2 rounded-full" style={{ background: judgeColor(id, theme) }} />
                      {JUDGES[id].short}
                    </div>
                    <div className="mt-0.5 text-xl font-semibold tabular text-slate-800 dark:text-slate-100">{formatScore(round.byJudge[id])}</div>
                  </div>
                ))}
              </div>
            </div>
          </CardBody>
          {rounds.length > 1 ? (
            <div className="flex items-center gap-2 border-t border-slate-100 px-5 py-2.5 dark:border-surface-dark-border">
              <span className="text-xs text-steel-500">查看轮次</span>
              <Segmented size="sm" value={round.round} onChange={(v) => setRoundNo(v)} options={rounds.map((r) => ({ value: r.round, label: `第 ${r.round} 轮 · ${formatScore(r.score)}` }))} />
            </div>
          ) : null}
        </Card>

        <Card>
          <CardHeader title="轮次得分趋势" subtitle="均分与三位评委分数随修改轮次的变化" />
          <CardBody className="text-slate-600 dark:text-slate-300">
            <div className="h-[220px]">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={trend} margin={{ left: -8, right: 16, top: 8, bottom: 0 }}>
                  <CartesianGrid vertical={false} stroke="currentColor" strokeOpacity={0.12} />
                  <XAxis dataKey="name" tick={{ fill: 'currentColor', fontSize: 11 }} axisLine={false} tickLine={false} />
                  <YAxis domain={[Math.max(0, Math.floor(Math.min(...rounds.map((r) => Math.min(r.score, ...Object.values(r.byJudge).map((x) => x ?? r.score)))) / 5) * 5 - 5), 100]} tick={{ fill: 'currentColor', fontSize: 11 }} axisLine={false} tickLine={false} />
                  <Tooltip content={<ChartTooltip />} />
                  <Legend wrapperStyle={{ fontSize: 12 }} />
                  <ReferenceLine y={target} stroke="#C9A227" strokeDasharray="4 4" label={{ value: `目标 ${target}`, fill: '#C9A227', fontSize: 11, position: 'insideTopRight' }} />
                  <Line type="monotone" dataKey="均分" stroke={NAVY[theme]} strokeWidth={3} dot={{ r: 4 }} />
                  <Line type="monotone" dataKey="投资人" stroke="#7F97BF" strokeWidth={1.5} strokeDasharray="3 3" dot={false} />
                  <Line type="monotone" dataKey="学者" stroke={judgeColor('professor', theme)} strokeWidth={1.5} strokeDasharray="3 3" dot={false} />
                  <Line type="monotone" dataKey="企业高管" stroke={judgeColor('executive', theme)} strokeWidth={1.5} strokeDasharray="3 3" dot={false} />
                </LineChart>
              </ResponsiveContainer>
            </div>
          </CardBody>
        </Card>
      </div>

      <RoundView round={round} target={target} />
    </div>
  );
}
