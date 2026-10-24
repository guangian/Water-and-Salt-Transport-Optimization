import { memo, useEffect, useMemo, useRef, useState, type ReactNode } from 'react';
import {
  AlertTriangle,
  ArrowDown,
  Bot,
  Brain,
  CheckCircle2,
  ChevronDown,
  ChevronRight,
  FileEdit,
  FilePlus2,
  FileX2,
  Flag,
  Gavel,
  Globe,
  Loader2,
  Pencil,
  Search,
  StopCircle,
  Terminal,
  Users,
  Wrench,
  XCircle,
  BookOpen,
  FolderSearch,
  ListChecks,
  Database,
} from 'lucide-react';
import type { StreamItem } from '../../lib/events';
import { STAGE_BY_ID } from '../../lib/stages';
import { formatBytes, formatTime, pretty, summarizeInput } from '../../lib/format';
import { cn } from '../../lib/cn';
import { Markdown } from '../domain/Markdown';
import { Badge } from '../ui/Badge';
import { STATUS_META } from '../domain/Brand';

export type AgentColorFn = (agentId: string) => string;

const TOOL_META: Record<string, { label: string; icon: ReactNode }> = {
  Read: { label: '读取文件', icon: <BookOpen /> },
  Write: { label: '写入文件', icon: <FilePlus2 /> },
  Edit: { label: '编辑文件', icon: <Pencil /> },
  Glob: { label: '查找文件', icon: <FolderSearch /> },
  Grep: { label: '搜索内容', icon: <Search /> },
  Bash: { label: '执行命令', icon: <Terminal /> },
  Agent: { label: '派出子代理', icon: <Users /> },
  mcp__bizgold__web_search: { label: '联网检索', icon: <Globe /> },
  mcp__bizgold__web_fetch: { label: '抓取网页', icon: <Globe /> },
  mcp__bizgold__update_progress: { label: '更新进度', icon: <ListChecks /> },
  mcp__bizgold__record_numbers: { label: '登记数字', icon: <Database /> },
};

function toolMeta(name: string): { label: string; icon: ReactNode } {
  return TOOL_META[name] ?? { label: name.replace(/^mcp__bizgold__/, ''), icon: <Wrench /> };
}

function AgentTag({ name, color, time, right }: { name: string; color: string; time: string; right?: ReactNode }) {
  return (
    <div className="mb-1 flex items-center gap-2 text-[11px] text-steel-500 dark:text-steel-400">
      <span className="inline-flex items-center gap-1.5 font-medium" style={{ color }}>
        <Bot className="h-3.5 w-3.5" />
        {name}
      </span>
      <span className="tabular">{formatTime(time)}</span>
      {right}
    </div>
  );
}

const TextItem = memo(function TextItem({ item, color, projectId }: { item: Extract<StreamItem, { kind: 'text' }>; color: string; projectId: string }) {
  const isMain = item.agent.id === 'main';
  return (
    <div className={cn('rounded-lg px-3 py-2', isMain ? 'bg-white dark:bg-surface-dark-card' : 'bg-slate-50/80 dark:bg-slate-800/40')} style={{ borderLeft: `3px solid ${color}` }}>
      <AgentTag name={item.agent.name} color={color} time={item.ts} right={item.streaming ? <Loader2 className="h-3 w-3 animate-spin text-emerald-500" /> : null} />
      <Markdown content={item.text} variant="stream" projectId={projectId} />
      {item.streaming ? <span className="ml-0.5 inline-block h-3.5 w-1.5 animate-pulse bg-gold-500 align-text-bottom" /> : null}
    </div>
  );
});

const ThinkingItem = memo(function ThinkingItem({ item, color }: { item: Extract<StreamItem, { kind: 'thinking' }>; color: string }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="rounded-lg border border-dashed border-slate-200 bg-slate-50/60 px-3 py-2 dark:border-slate-700 dark:bg-slate-800/30">
      <button onClick={() => setOpen((o) => !o)} className="flex w-full items-center gap-2 text-left text-xs text-steel-600 dark:text-steel-300">
        {open ? <ChevronDown className="h-3.5 w-3.5" /> : <ChevronRight className="h-3.5 w-3.5" />}
        <Brain className="h-3.5 w-3.5" style={{ color }} />
        <span className="font-medium">思考过程</span>
        <span className="text-steel-400">· {item.agent.name} · {formatTime(item.ts)}</span>
        {item.streaming ? <Loader2 className="h-3 w-3 animate-spin text-emerald-500" /> : null}
        {!open ? <span className="ml-2 min-w-0 flex-1 truncate text-steel-400">{item.text.slice(0, 80)}</span> : null}
      </button>
      {open ? <p className="mt-2 whitespace-pre-wrap text-xs leading-relaxed text-steel-600 dark:text-steel-300">{item.text}</p> : null}
    </div>
  );
});

const ToolItem = memo(function ToolItem({ item, color }: { item: Extract<StreamItem, { kind: 'tool' }>; color: string }) {
  const [open, setOpen] = useState(false);
  const meta = toolMeta(item.name);
  const summary = summarizeInput(item.input, 96);
  const pending = !item.result;
  const isError = item.result?.isError;
  const isAgent = item.name === 'Agent';
  const input = item.input as Record<string, unknown> | undefined;
  return (
    <div className={cn('rounded-lg border bg-white dark:bg-surface-dark-card', isError ? 'border-red-200 dark:border-red-900/60' : 'border-slate-200 dark:border-surface-dark-border')}>
      <button onClick={() => setOpen((o) => !o)} className="flex w-full items-center gap-2 px-3 py-2 text-left">
        <span className={cn('inline-flex h-6 w-6 shrink-0 items-center justify-center rounded-md [&>svg]:h-3.5 [&>svg]:w-3.5', isAgent ? 'bg-gold-50 text-gold-700 dark:bg-gold-900/30 dark:text-gold-300' : 'bg-slate-100 text-slate-600 dark:bg-slate-800 dark:text-slate-300')}>{meta.icon}</span>
        <span className="shrink-0 text-xs font-medium text-slate-700 dark:text-slate-200">{meta.label}</span>
        <span className="hidden shrink-0 font-mono text-[11px] text-steel-400 sm:inline">{item.name}</span>
        <span className="min-w-0 flex-1 truncate font-mono text-xs text-steel-600 dark:text-steel-300">{isAgent && input ? String(input.description ?? summary) : summary}</span>
        <span className="hidden shrink-0 text-[11px] tabular text-steel-400 md:inline" style={{ color }}>
          {item.agent.name}
        </span>
        <span className="shrink-0 text-[11px] tabular text-steel-400">{formatTime(item.ts)}</span>
        {pending ? <Loader2 className="h-3.5 w-3.5 shrink-0 animate-spin text-emerald-500" /> : isError ? <XCircle className="h-3.5 w-3.5 shrink-0 text-red-500" /> : <CheckCircle2 className="h-3.5 w-3.5 shrink-0 text-emerald-500" />}
        {open ? <ChevronDown className="h-3.5 w-3.5 shrink-0 text-slate-400" /> : <ChevronRight className="h-3.5 w-3.5 shrink-0 text-slate-400" />}
      </button>
      {open ? (
        <div className="grid gap-3 border-t border-slate-100 px-3 py-3 text-xs dark:border-surface-dark-border md:grid-cols-2">
          <div className="min-w-0">
            <div className="mb-1 font-medium text-steel-500">输入</div>
            <pre className="max-h-64 overflow-auto whitespace-pre-wrap break-words rounded-md bg-slate-50 p-2.5 font-mono text-[11px] leading-relaxed text-slate-700 dark:bg-slate-900 dark:text-slate-300">{item.input === undefined ? '（无）' : pretty(item.input)}</pre>
          </div>
          <div className="min-w-0">
            <div className="mb-1 flex items-center gap-2 font-medium text-steel-500">
              结果
              {item.result ? <span className="text-steel-400">{formatTime(item.result.ts)}</span> : null}
            </div>
            <pre className={cn('max-h-64 overflow-auto whitespace-pre-wrap break-words rounded-md p-2.5 font-mono text-[11px] leading-relaxed', isError ? 'bg-red-50 text-red-700 dark:bg-red-900/20 dark:text-red-200' : 'bg-slate-50 text-slate-700 dark:bg-slate-900 dark:text-slate-300')}>
              {item.result ? item.result.output || '（空）' : '等待结果…'}
            </pre>
          </div>
        </div>
      ) : null}
    </div>
  );
});

function FileItem({ item, color, onOpenFile }: { item: Extract<StreamItem, { kind: 'file' }>; color: string; onOpenFile?: (path: string) => void }) {
  const icon = item.op === 'write' ? <FilePlus2 /> : item.op === 'edit' ? <FileEdit /> : <FileX2 />;
  const verb = item.op === 'write' ? '写入' : item.op === 'edit' ? '编辑' : '删除';
  return (
    <button
      onClick={() => item.op !== 'delete' && onOpenFile?.(item.path)}
      disabled={item.op === 'delete' || !onOpenFile}
      className="flex w-full items-center gap-2 rounded-md border border-brand-100 bg-brand-50/60 px-3 py-1.5 text-left text-xs transition hover:bg-brand-50 disabled:cursor-default dark:border-brand-900/60 dark:bg-brand-900/20 dark:hover:bg-brand-900/40"
      title={item.op === 'delete' ? undefined : '在产出物中预览'}
    >
      <span className="inline-flex text-brand-500 dark:text-brand-200 [&>svg]:h-3.5 [&>svg]:w-3.5">{icon}</span>
      <span className="text-brand-700 dark:text-brand-100">{verb}</span>
      <span className="min-w-0 flex-1 truncate font-mono text-slate-700 dark:text-slate-200">{item.path}</span>
      {item.size != null ? <span className="tabular text-steel-500">{formatBytes(item.size)}</span> : null}
      <span className="hidden text-[11px] text-steel-400 md:inline" style={{ color }}>
        {item.agent.name}
      </span>
      <span className="text-[11px] tabular text-steel-400">{formatTime(item.ts)}</span>
    </button>
  );
}

function StageItem({ item }: { item: Extract<StreamItem, { kind: 'stage' }> }) {
  const def = STAGE_BY_ID[item.stage];
  const tone = item.status === 'done' ? 'text-brand-600 dark:text-brand-200' : item.status === 'running' ? 'text-emerald-600 dark:text-emerald-300' : item.status === 'failed' ? 'text-red-600 dark:text-red-300' : 'text-steel-500';
  const text = item.status === 'running' ? '开始' : item.status === 'done' ? '完成' : item.status === 'failed' ? '失败' : item.status === 'skipped' ? '跳过' : '等待';
  return (
    <div className="my-1 flex items-center gap-3 text-xs">
      <span className="h-px flex-1 bg-slate-200 dark:bg-slate-700" />
      <span className={cn('inline-flex items-center gap-1.5 font-medium', tone)}>
        <Flag className="h-3.5 w-3.5" />
        阶段 {def?.index ?? ''} {def?.name ?? item.stage}
        {item.round ? `（第 ${item.round} 轮）` : ''} {text}
        {item.note ? <span className="font-normal text-steel-500">· {item.note}</span> : null}
        <span className="font-normal tabular text-steel-400">{formatTime(item.ts)}</span>
      </span>
      <span className="h-px flex-1 bg-slate-200 dark:bg-slate-700" />
    </div>
  );
}

function StatusItem({ item }: { item: Extract<StreamItem, { kind: 'status' }> }) {
  const meta = STATUS_META[item.status];
  return (
    <div className="flex justify-center py-0.5">
      <Badge tone={meta.tone} dot={meta.dot}>
        {meta.label}
        {item.message ? ` · ${item.message}` : ''}
        <span className="font-normal tabular opacity-70">{formatTime(item.ts)}</span>
      </Badge>
    </div>
  );
}

function SubagentItem({ item, color }: { item: Extract<StreamItem, { kind: 'subagent' }>; color: string }) {
  return (
    <div className="flex items-center gap-2 pl-1 text-xs text-steel-600 dark:text-steel-300">
      <span className="inline-flex h-5 w-5 items-center justify-center rounded-full" style={{ backgroundColor: `${color}22`, color }}>
        {item.status === 'start' ? <Bot className="h-3 w-3" /> : <CheckCircle2 className="h-3 w-3" />}
      </span>
      <span className="font-medium" style={{ color }}>
        {item.name}
      </span>
      <span>{item.status === 'start' ? '已启动' : '已完成'}</span>
      {item.summary ? <span className="min-w-0 flex-1 truncate text-steel-500">· {item.summary}</span> : <span className="flex-1" />}
      <span className="tabular text-steel-400">{formatTime(item.ts)}</span>
    </div>
  );
}

function ReviewItem({ item }: { item: Extract<StreamItem, { kind: 'review' }> }) {
  return (
    <div className="flex items-center gap-3 rounded-lg border border-gold-200 bg-gold-50/70 px-3 py-2 text-sm dark:border-gold-800 dark:bg-gold-900/20">
      <Gavel className="h-4 w-4 text-gold-700 dark:text-gold-300" />
      <span className="font-medium text-slate-800 dark:text-slate-100">第 {item.round} 轮评审</span>
      <span className="text-2xl font-semibold tabular text-brand-700 dark:text-gold-200">{item.score.toFixed(1)}</span>
      <Badge tone={item.hardGatePassed ? 'green' : 'red'}>{item.hardGatePassed ? '硬门槛通过' : '硬门槛未通过'}</Badge>
      <span className="ml-auto text-xs tabular text-steel-400">{formatTime(item.ts)}</span>
    </div>
  );
}

function ErrorItem({ item }: { item: Extract<StreamItem, { kind: 'error' }> }) {
  return (
    <div className="flex items-start gap-3 rounded-lg border border-red-200 bg-red-50 px-3 py-2.5 text-sm dark:border-red-900/60 dark:bg-red-900/20">
      <AlertTriangle className="mt-0.5 h-4 w-4 shrink-0 text-red-600 dark:text-red-300" />
      <div className="min-w-0 flex-1">
        <div className="font-medium text-red-800 dark:text-red-200">{item.recoverable ? '可恢复错误' : '错误'}</div>
        <div className="mt-0.5 break-words text-xs text-red-700 dark:text-red-200/90">{item.message}</div>
      </div>
      <span className="text-xs tabular text-red-400">{formatTime(item.ts)}</span>
    </div>
  );
}

function DoneItem({ item }: { item: Extract<StreamItem, { kind: 'done' }> }) {
  const meta =
    item.subtype === 'success'
      ? { cls: 'border-emerald-200 bg-emerald-50 text-emerald-800 dark:border-emerald-900/60 dark:bg-emerald-900/20 dark:text-emerald-200', icon: <CheckCircle2 />, label: '运行完成' }
      : item.subtype === 'stopped'
        ? { cls: 'border-slate-200 bg-slate-50 text-slate-700 dark:border-slate-700 dark:bg-slate-800/60 dark:text-slate-200', icon: <StopCircle />, label: '已停止' }
        : { cls: 'border-red-200 bg-red-50 text-red-800 dark:border-red-900/60 dark:bg-red-900/20 dark:text-red-200', icon: <XCircle />, label: '运行失败' };
  return (
    <div className={cn('flex items-center gap-3 rounded-lg border px-3 py-2.5 text-sm', meta.cls)}>
      <span className="[&>svg]:h-4 [&>svg]:w-4">{meta.icon}</span>
      <span className="font-medium">{meta.label}</span>
      {item.result ? <span className="min-w-0 flex-1 truncate text-xs opacity-80">{item.result}</span> : <span className="flex-1" />}
      <span className="text-xs tabular opacity-60">{formatTime(item.ts)}</span>
    </div>
  );
}

export const EventItem = memo(function EventItem({ item, colorOf, projectId, onOpenFile }: { item: StreamItem; colorOf: AgentColorFn; projectId: string; onOpenFile?: (path: string) => void }) {
  const color = colorOf(item.kind === 'subagent' ? item.toolUseId : item.agent.id);
  switch (item.kind) {
    case 'text':
      return <TextItem item={item} color={color} projectId={projectId} />;
    case 'thinking':
      return <ThinkingItem item={item} color={color} />;
    case 'tool':
      return <ToolItem item={item} color={color} />;
    case 'file':
      return <FileItem item={item} color={color} onOpenFile={onOpenFile} />;
    case 'stage':
      return <StageItem item={item} />;
    case 'status':
      return <StatusItem item={item} />;
    case 'subagent':
      return <SubagentItem item={item} color={color} />;
    case 'review':
      return <ReviewItem item={item} />;
    case 'error':
      return <ErrorItem item={item} />;
    case 'done':
      return <DoneItem item={item} />;
    default:
      return null;
  }
});

const PAGE = 300;

export function EventStream({
  items,
  colorOf,
  projectId,
  onOpenFile,
  autoScroll,
  onAutoScrollChange,
  className,
  emptyHint,
}: {
  items: StreamItem[];
  colorOf: AgentColorFn;
  projectId: string;
  onOpenFile?: (path: string) => void;
  autoScroll: boolean;
  onAutoScrollChange: (v: boolean) => void;
  className?: string;
  emptyHint?: ReactNode;
}) {
  const ref = useRef<HTMLDivElement>(null);
  const [atBottom, setAtBottom] = useState(true);
  const [limit, setLimit] = useState(PAGE);
  const visible = useMemo(() => (items.length > limit ? items.slice(items.length - limit) : items), [items, limit]);
  const hidden = items.length - visible.length;

  useEffect(() => {
    const el = ref.current;
    if (!el || !autoScroll || !atBottom) return;
    el.scrollTop = el.scrollHeight;
  }, [items, autoScroll, atBottom]);

  const onScroll = () => {
    const el = ref.current;
    if (!el) return;
    const dist = el.scrollHeight - el.scrollTop - el.clientHeight;
    const nowAtBottom = dist < 32;
    if (nowAtBottom !== atBottom) setAtBottom(nowAtBottom);
    // 用户向上滚动 → 暂停自动滚动；回到底部 → 恢复
    if (!nowAtBottom && autoScroll) onAutoScrollChange(false);
  };

  const jumpToBottom = () => {
    const el = ref.current;
    if (!el) return;
    el.scrollTop = el.scrollHeight;
    setAtBottom(true);
    onAutoScrollChange(true);
  };

  return (
    <div className={cn('relative', className)}>
      <div ref={ref} onScroll={onScroll} className="h-full overflow-y-auto px-4 py-3">
        {hidden > 0 ? (
          <div className="mb-2 flex justify-center">
            <button onClick={() => setLimit((l) => l + PAGE)} className="rounded-full border border-slate-200 bg-white px-3 py-1 text-xs text-steel-600 hover:bg-slate-50 dark:border-surface-dark-border dark:bg-surface-dark-card dark:text-steel-300">
              显示更早的 {Math.min(hidden, PAGE)} 条（共隐藏 {hidden} 条）
            </button>
          </div>
        ) : null}
        {visible.length === 0 ? (
          <div className="flex h-full min-h-[200px] flex-col items-center justify-center gap-2 text-center text-sm text-steel-500">
            <Bot className="h-8 w-8 text-slate-300 dark:text-slate-600" />
            {emptyHint ?? '暂无事件'}
          </div>
        ) : (
          <div className="space-y-2">
            {visible.map((it) => (
              <EventItem key={it.id} item={it} colorOf={colorOf} projectId={projectId} onOpenFile={onOpenFile} />
            ))}
          </div>
        )}
      </div>
      {!atBottom || !autoScroll ? (
        <button
          onClick={jumpToBottom}
          className="absolute bottom-4 left-1/2 inline-flex -translate-x-1/2 items-center gap-1.5 rounded-full bg-brand-600 px-3 py-1.5 text-xs font-medium text-white shadow-card-hover transition hover:bg-brand-500"
        >
          <ArrowDown className="h-3.5 w-3.5" />
          {autoScroll ? '回到底部' : '已暂停自动滚动 · 回到底部'}
        </button>
      ) : null}
    </div>
  );
}
