import { useMemo, useState, type ReactNode } from 'react';
import { Activity, Bot, CheckCircle2, Coins, Filter, Loader2, Pause, Play, Repeat, Timer, Wifi, WifiOff, Zap } from 'lucide-react';
import type { Project } from '../../types';
import type { SseStatus } from '../../api';
import { filterItems, type AgentLane, type ConsoleState, type StreamFilter } from '../../lib/events';
import { formatDuration, formatTokens, formatUsd } from '../../lib/format';
import { cn } from '../../lib/cn';
import { useTheme } from '../../hooks/useTheme';
import { Card, CardBody, CardHeader } from '../ui/Card';
import { Segmented, Select } from '../ui/Form';
import { StageProgress } from './StageProgress';
import { EventStream, type AgentColorFn } from './EventStream';

const AGENT_PALETTE = ['#2E8B6E', '#6D4FC2', '#2C7FB8', '#D97706', '#B5533C', '#0E7490', '#9333EA', '#65A30D'];
const AGENT_PALETTE_DARK = ['#5FC3A0', '#A994F0', '#6FB3E8', '#F2B35A', '#E58A74', '#4FC2D9', '#C08BFA', '#A3D65A'];
// 主代理用品牌藏青；深色背景下藏青几乎不可见，改用浅品牌蓝
const MAIN_COLOR = { light: '#14305C', dark: '#9DBCF0' } as const;
const FALLBACK_COLOR = { light: '#5B6B82', dark: '#9AA8BC' } as const;

export function useAgentColors(agents: AgentLane[]): AgentColorFn {
  const [theme] = useTheme();
  return useMemo(() => {
    const palette = theme === 'dark' ? AGENT_PALETTE_DARK : AGENT_PALETTE;
    const map = new Map<string, string>();
    let i = 0;
    for (const a of agents) {
      if (a.id === 'main') map.set(a.id, MAIN_COLOR[theme]);
      else map.set(a.id, palette[i++ % palette.length]);
    }
    return (id: string) => map.get(id) ?? (id === 'main' ? MAIN_COLOR[theme] : FALLBACK_COLOR[theme]);
  }, [agents, theme]);
}

export function SseIndicator({ status, retryInMs }: { status: SseStatus; retryInMs: number | null }) {
  const meta =
    status === 'open'
      ? { cls: 'text-emerald-600 dark:text-emerald-300', icon: <Wifi />, label: '实时连接' }
      : status === 'connecting'
        ? { cls: 'text-steel-500', icon: <Loader2 className="animate-spin" />, label: '连接中' }
        : status === 'reconnecting'
          ? { cls: 'text-amber-600 dark:text-amber-300', icon: <Loader2 className="animate-spin" />, label: `重连中${retryInMs ? `（${Math.round(retryInMs / 1000)}s）` : ''}` }
          : { cls: 'text-steel-400', icon: <WifiOff />, label: '已断开' };
  return (
    <span className={cn('inline-flex items-center gap-1.5 text-xs font-medium [&>svg]:h-3.5 [&>svg]:w-3.5', meta.cls)}>
      {meta.icon}
      {meta.label}
    </span>
  );
}

function LaneCard({ lane, color, active, onClick, now }: { lane: AgentLane; color: string; active: boolean; onClick: () => void; now: number }) {
  const start = Date.parse(lane.startedAt);
  const end = lane.endedAt ? Date.parse(lane.endedAt) : now;
  const dur = Number.isNaN(start) ? null : Math.max(0, end - start);
  return (
    <button
      onClick={onClick}
      className={cn(
        'w-full rounded-lg border px-3 py-2.5 text-left transition',
        active ? 'border-brand-300 bg-brand-50/60 dark:border-brand-400/60 dark:bg-brand-900/30' : 'border-slate-200 bg-white hover:border-slate-300 dark:border-surface-dark-border dark:bg-surface-dark-card dark:hover:border-slate-600',
      )}
    >
      <div className="flex items-center gap-2">
        <span className="inline-flex h-6 w-6 shrink-0 items-center justify-center rounded-full" style={{ backgroundColor: `${color}22`, color }}>
          <Bot className="h-3.5 w-3.5" />
        </span>
        <span className="min-w-0 flex-1 truncate text-sm font-medium text-slate-800 dark:text-slate-100">{lane.name}</span>
        {lane.status === 'running' ? (
          <span className="inline-flex items-center gap-1 text-[11px] font-medium text-emerald-600 dark:text-emerald-300">
            <span className="h-1.5 w-1.5 animate-pulse-dot rounded-full bg-emerald-500" />
            运行中
          </span>
        ) : lane.status === 'done' ? (
          <span className="inline-flex items-center gap-1 text-[11px] text-steel-500">
            <CheckCircle2 className="h-3.5 w-3.5 text-brand-500 dark:text-brand-300" />
            完成
          </span>
        ) : (
          <span className="text-[11px] text-steel-400">空闲</span>
        )}
      </div>
      <p className="mt-1.5 line-clamp-2 min-h-[2rem] text-xs leading-4 text-steel-600 dark:text-steel-300">{lane.lastText || '（尚无输出）'}</p>
      <div className="mt-1.5 flex items-center gap-3 text-[11px] tabular text-steel-500">
        <span>工具 {lane.toolCalls}</span>
        <span>文本 {formatTokens(lane.textChars)} 字</span>
        <span className="ml-auto">{formatDuration(dur, { compact: true })}</span>
      </div>
    </button>
  );
}

function Counter({ icon, label, value, sub, wide }: { icon: ReactNode; label: string; value: string; sub?: string; wide?: boolean }) {
  return (
    <div className={cn('min-w-0', wide && 'col-span-2')}>
      <div className="flex items-center gap-1.5 text-[11px] text-steel-500 [&>svg]:h-3 [&>svg]:w-3">
        {icon}
        {label}
      </div>
      <div className="truncate text-sm font-semibold tabular leading-5 text-slate-800 dark:text-slate-100">
        {value}
        {sub ? <span className="ml-1.5 text-[11px] font-normal text-steel-400">{sub}</span> : null}
      </div>
    </div>
  );
}

export function ConsoleTab({
  project,
  state,
  sse,
  retryInMs,
  now,
  fetchedAt,
  onOpenFile,
}: {
  project: Project;
  state: ConsoleState;
  sse: SseStatus;
  retryInMs: number | null;
  now: number;
  fetchedAt: number;
  onOpenFile: (path: string) => void;
}) {
  const [filter, setFilter] = useState<StreamFilter>('all');
  const [agent, setAgent] = useState<string>('all');
  const [autoScroll, setAutoScroll] = useState(true);
  const colorOf = useAgentColors(state.agents);
  const items = useMemo(() => filterItems(state.items, filter, agent), [state.items, filter, agent]);

  const usage = state.usage.inputTokens || state.usage.outputTokens ? state.usage : project.usage;
  const running = project.status === 'running';
  const durationMs = project.usage.durationMs + (running ? Math.max(0, now - fetchedAt) : 0);
  const subagents = state.agents.filter((a) => a.id !== 'main');
  const runningLanes = subagents.filter((a) => a.status === 'running').length;

  return (
    <div className="space-y-4">
      <Card>
        <CardBody className="py-4">
          <StageProgress stages={project.stages} now={now} />
        </CardBody>
      </Card>

      <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_340px] 2xl:grid-cols-[minmax(0,1fr)_380px]">
        <Card className="flex min-h-[360px] flex-col overflow-hidden xl:h-[calc(100vh-30rem)]">
          <div className="flex flex-wrap items-center gap-x-3 gap-y-2 border-b border-slate-100 px-4 py-2.5 dark:border-surface-dark-border">
            <h3 className="flex items-center gap-2 text-sm font-semibold text-slate-800 dark:text-slate-100">
              <Activity className="h-4 w-4 text-brand-500 dark:text-brand-200" />
              实时事件流
            </h3>
            <SseIndicator status={sse} retryInMs={retryInMs} />
            <span className="text-xs tabular text-steel-400">{state.items.length} 项 · seq {state.lastSeq}</span>
            <div className="ml-auto flex flex-wrap items-center gap-2">
              <Segmented
                size="sm"
                value={filter}
                onChange={setFilter}
                options={[
                  { value: 'all', label: '全部' },
                  { value: 'text', label: '只看文本' },
                  { value: 'tools', label: '只看工具' },
                ]}
              />
              <div className="w-40">
                <Select value={agent} onChange={(e) => setAgent(e.target.value)} className="h-8 py-1 text-xs">
                  <option value="all">全部代理</option>
                  {state.agents.map((a) => (
                    <option key={a.id} value={a.id}>
                      {a.name}
                    </option>
                  ))}
                </Select>
              </div>
              <button
                onClick={() => setAutoScroll((v) => !v)}
                className={cn(
                  'inline-flex h-8 items-center gap-1.5 rounded-lg border px-2.5 text-xs font-medium transition',
                  autoScroll ? 'border-emerald-200 bg-emerald-50 text-emerald-700 dark:border-emerald-800 dark:bg-emerald-900/30 dark:text-emerald-200' : 'border-slate-200 bg-white text-steel-600 dark:border-surface-dark-border dark:bg-surface-dark-card dark:text-steel-300',
                )}
                title="切换自动滚动"
              >
                {autoScroll ? <Play className="h-3.5 w-3.5" /> : <Pause className="h-3.5 w-3.5" />}
                自动滚动
              </button>
            </div>
          </div>
          <EventStream
            className="min-h-0 flex-1 bg-slate-50/50 dark:bg-surface-dark/40"
            items={items}
            colorOf={colorOf}
            projectId={project.id}
            onOpenFile={onOpenFile}
            autoScroll={autoScroll}
            onAutoScrollChange={setAutoScroll}
            emptyHint={
              state.items.length === 0
                ? project.status === 'idle'
                  ? '项目尚未运行。点击右上角「开始」后，这里将实时显示主代理与子代理的工作过程。'
                  : sse === 'open'
                    ? '暂无事件'
                    : '正在连接事件流…'
                : (
                  <span className="inline-flex items-center gap-1.5">
                    <Filter className="h-4 w-4" />
                    当前筛选条件下没有事件
                  </span>
                )
            }
          />
        </Card>

        <div className="flex min-h-0 flex-col gap-4 xl:h-[calc(100vh-30rem)] xl:min-h-[360px]">
          <Card className="flex min-h-0 flex-1 flex-col">
            <CardHeader
              title="子代理泳道"
              subtitle={subagents.length ? `${subagents.length} 个子代理 · ${runningLanes} 个运行中` : '主代理按需派出研究员、建模师、主笔、评委…'}
              className="px-4 py-3"
              actions={
                agent !== 'all' ? (
                  <button onClick={() => setAgent('all')} className="text-xs text-brand-600 hover:underline dark:text-brand-200">
                    清除筛选
                  </button>
                ) : null
              }
            />
            <div className="max-h-[420px] min-h-[120px] flex-1 space-y-2 overflow-y-auto px-4 py-3 xl:max-h-none">
              {state.agents.length === 0 ? (
                <div className="py-6 text-center text-xs text-steel-500">尚无代理活动</div>
              ) : (
                state.agents.map((lane) => <LaneCard key={lane.id} lane={lane} color={colorOf(lane.id)} active={agent === lane.id} onClick={() => setAgent(agent === lane.id ? 'all' : lane.id)} now={now} />)
              )}
            </div>
          </Card>

          <Card className="shrink-0">
            <CardHeader title="用量统计" subtitle={running ? '运行中实时更新' : '最近一次运行'} className="px-4 py-2.5" />
            <CardBody className="grid grid-cols-2 gap-x-4 gap-y-2.5 px-4 py-3">
              <Counter icon={<Zap />} label="输入 Token" value={formatTokens(usage.inputTokens)} sub={`缓存 ${formatTokens(usage.cacheReadTokens)}`} />
              <Counter icon={<Zap />} label="输出 Token" value={formatTokens(usage.outputTokens)} />
              <Counter icon={<Coins />} label="费用估算" value={formatUsd(usage.costEstimateUsd)} />
              <Counter icon={<Repeat />} label="轮次" value={String(usage.turns)} />
              <Counter wide icon={<Timer />} label="用时" value={formatDuration(durationMs)} sub={state.firstTs ? `首个事件 ${new Date(state.firstTs).toLocaleTimeString('zh-CN', { hour12: false })}` : undefined} />
            </CardBody>
          </Card>
        </div>
      </div>
    </div>
  );
}
