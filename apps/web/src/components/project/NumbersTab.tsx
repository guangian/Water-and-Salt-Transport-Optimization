import { useMemo, useState } from 'react';
import { AlertTriangle, CheckCircle2, Database, FlaskConical, HelpCircle, Search } from 'lucide-react';
import type { NumberEntry, Project } from '../../types';
import { api } from '../../api';
import { useAsync } from '../../hooks/useAsync';
import { formatDateTime, formatNumber } from '../../lib/format';
import { cn } from '../../lib/cn';
import { Card, CardHeader } from '../ui/Card';
import { Badge } from '../ui/Badge';
import { EmptyState, Spinner, Stat } from '../ui/Misc';
import { Input, Segmented } from '../ui/Form';

type Filter = 'all' | 'synthetic' | 'mismatch';

function Consistency({ v }: { v: NumberEntry['consistency'] }) {
  if (v === 'ok')
    return (
      <span className="inline-flex items-center gap-1 text-xs text-emerald-600 dark:text-emerald-300">
        <CheckCircle2 className="h-3.5 w-3.5" />
        一致
      </span>
    );
  if (v === 'mismatch')
    return (
      <span className="inline-flex items-center gap-1 text-xs font-medium text-red-600 dark:text-red-300">
        <AlertTriangle className="h-3.5 w-3.5" />
        不一致
      </span>
    );
  return (
    <span className="inline-flex items-center gap-1 text-xs text-steel-400">
      <HelpCircle className="h-3.5 w-3.5" />
      未校验
    </span>
  );
}

export function NumbersTab({ project, refreshKey }: { project: Project; refreshKey: number }) {
  const { data, loading, error } = useAsync(() => api.getNumbers(project.id), [project.id, refreshKey]);
  const [query, setQuery] = useState('');
  const [filter, setFilter] = useState<Filter>('all');

  const entries = data?.entries ?? [];
  const synthetic = entries.filter((e) => e.synthetic).length;
  const mismatches = entries.filter((e) => e.consistency === 'mismatch').length;
  const list = useMemo(() => {
    const q = query.trim().toLowerCase();
    return entries.filter((e) => {
      if (filter === 'synthetic' && !e.synthetic) return false;
      if (filter === 'mismatch' && e.consistency !== 'mismatch') return false;
      if (!q) return true;
      return [e.key, e.label, e.source ?? '', ...(e.usedIn ?? [])].some((s) => s.toLowerCase().includes(q));
    });
  }, [entries, query, filter]);

  if (loading && !data) return <Spinner className="py-20" label="加载数字账本…" />;
  if (error) return <EmptyState icon={<Database />} title="数字账本加载失败" description={error.message} />;
  if (!entries.length) {
    return (
      <Card>
        <EmptyState icon={<Database />} title="数字账本为空" description="财务阶段完成后，财务建模师会把关键数字（预算、KPI、市场规模、样本量…）登记到 03_财务/numbers.json；正文、路演与讲稿必须统一引用这里的数值。" />
      </Card>
    );
  }

  return (
    <div className="space-y-4">
      <div className="grid grid-cols-2 gap-4 lg:grid-cols-4">
        <Stat label="登记数字" value={entries.length} sub={data?.generatedAt ? `更新于 ${formatDateTime(data.generatedAt)}` : undefined} icon={<Database />} tone="brand" />
        <Stat label="示例数据" value={synthetic} sub="提交前须替换为真实调研结果" icon={<FlaskConical />} tone="gold" />
        <Stat label="跨材料不一致" value={mismatches} sub={data?.consistency.consistent ? '策划书 / 路演 / 讲稿一致' : `${data?.consistency.mismatches.length ?? 0} 处需修正`} icon={<AlertTriangle />} tone={mismatches ? 'steel' : 'green'} />
        <Stat label="已被引用" value={entries.filter((e) => e.usedIn?.length).length} sub={`默认单位 ${data?.unit ?? '万元'}`} icon={<CheckCircle2 />} tone="green" />
      </div>

      <Card>
        <CardHeader
          title="数字账本"
          subtitle="key 是全文唯一引用标识；来源精确到 文件!工作表!单元格 或章节"
          actions={
            <>
              <Segmented
                size="sm"
                value={filter}
                onChange={setFilter}
                options={[
                  { value: 'all', label: `全部 ${entries.length}` },
                  { value: 'synthetic', label: `示例 ${synthetic}` },
                  { value: 'mismatch', label: `不一致 ${mismatches}` },
                ]}
              />
              <div className="relative w-48">
                <Search className="pointer-events-none absolute left-2.5 top-1/2 h-3.5 w-3.5 -translate-y-1/2 text-slate-400" />
                <Input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="搜索 key / 名称 / 来源" className="h-8 py-1 pl-8 text-xs" />
              </div>
            </>
          }
        />
        <div className="max-h-[calc(100vh-24rem)] overflow-auto">
          <table className="table-base">
            <thead>
              <tr className="[&>th]:whitespace-nowrap">
                <th>Key</th>
                <th>名称</th>
                <th className="text-right">值</th>
                <th>单位</th>
                <th>来源</th>
                <th>使用位置</th>
                <th>示例数据</th>
                <th>一致性</th>
              </tr>
            </thead>
            <tbody>
              {list.map((e) => (
                <tr key={e.key} className={cn(e.consistency === 'mismatch' && 'bg-red-50/40 dark:bg-red-900/10')}>
                  <td className="font-mono text-xs text-brand-700 dark:text-brand-200">{e.key}</td>
                  <td className="text-[13px] text-slate-800 dark:text-slate-100">{e.label}</td>
                  <td className="text-right text-[13px] font-semibold tabular">{formatNumber(e.value, e.unit)}</td>
                  <td className="text-xs text-steel-500">{e.unit}</td>
                  <td className="max-w-[260px] break-all font-mono text-[11px] text-steel-600 dark:text-steel-300">{e.source ?? '—'}</td>
                  <td>
                    <div className="flex max-w-[280px] flex-wrap gap-1">
                      {(e.usedIn ?? []).map((u) => (
                        <span key={u} className="chip font-mono text-[10px]" title={u}>
                          {u.length > 28 ? `…${u.slice(-26)}` : u}
                        </span>
                      ))}
                      {!e.usedIn?.length ? <span className="text-xs text-steel-400">尚未引用</span> : null}
                    </div>
                  </td>
                  <td>{e.synthetic ? <Badge tone="gold">示例</Badge> : <span className="text-xs text-steel-400">—</span>}</td>
                  <td className="whitespace-nowrap">
                    <Consistency v={e.consistency} />
                  </td>
                </tr>
              ))}
              {list.length === 0 ? (
                <tr>
                  <td colSpan={8} className="py-10 text-center text-steel-500">
                    没有匹配的数字
                  </td>
                </tr>
              ) : null}
            </tbody>
          </table>
        </div>
      </Card>

      {data?.consistency.mismatches.length ? (
        <Card>
          <CardHeader title="不一致明细" subtitle="同一 key 在不同材料中的取值" actions={<Badge tone="red">{data.consistency.mismatches.length} 处</Badge>} />
          <table className="table-base">
            <thead>
              <tr>
                <th>关键数字</th>
                <th>策划书</th>
                <th>路演 PPT</th>
                <th>讲稿</th>
                <th>题库</th>
              </tr>
            </thead>
            <tbody>
              {data.consistency.mismatches.map((m) => {
                const vals = [m.book, m.deck, m.script, m.qa];
                const counts = vals.reduce<Record<string, number>>((acc, v) => {
                  if (v != null && v !== '') acc[String(v)] = (acc[String(v)] ?? 0) + 1;
                  return acc;
                }, {});
                const majority = Object.entries(counts).sort((a, b) => b[1] - a[1])[0]?.[0];
                return (
                  <tr key={m.key}>
                    <td>
                      <div className="text-[13px] text-slate-800 dark:text-slate-100">{m.label ?? m.key}</div>
                      <div className="font-mono text-[11px] text-steel-400">{m.key}</div>
                    </td>
                    {vals.map((v, i) => (
                      <td key={i} className={cn('tabular', v != null && v !== '' && String(v) !== majority && 'font-semibold text-red-600 dark:text-red-300')}>
                        {v == null || v === '' ? <span className="text-steel-400">—</span> : formatNumber(v)}
                      </td>
                    ))}
                  </tr>
                );
              })}
            </tbody>
          </table>
        </Card>
      ) : null}
    </div>
  );
}
