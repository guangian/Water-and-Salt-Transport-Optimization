import { useEffect, useMemo, useRef, useState } from 'react';
import { Bot, CheckCircle2, ChevronDown, ChevronRight, FileEdit, Loader2, MessageSquareText, Send, Sparkles, Trash2, User, Wrench, XCircle } from 'lucide-react';
import type { Project } from '../../types';
import { api } from '../../api';
import type { ConsoleState, StreamItem } from '../../lib/events';
import { formatBytes, formatTime, summarizeInput } from '../../lib/format';
import { loadJson, saveJson } from '../../lib/storage';
import { cn } from '../../lib/cn';
import { Card } from '../ui/Card';
import { Button } from '../ui/Button';
import { Textarea } from '../ui/Form';
import { useToast } from '../ui/Toast';
import { Markdown } from '../domain/Markdown';

interface ChatMessage {
  id: string;
  text: string;
  at: string;
  /** 发送时的事件 seq，之后的事件视为本轮回复 */
  afterSeq: number;
}

const SAMPLES = ['把第三章竞品分析再加一个与东道汝窑、朱氏汝瓷的对比表', '预算总额口径统一为 292 万元，并同步修改讲稿与路演', '讲稿压缩到 9 分钟以内，保留三个核心数据', '第二章补充问卷样本结构表（年龄 / 城市 / 收入）'];

type ReplyPart = { kind: 'text'; item: Extract<StreamItem, { kind: 'text' }> } | { kind: 'tools'; items: Extract<StreamItem, { kind: 'tool' | 'file' }>[] } | { kind: 'error'; item: Extract<StreamItem, { kind: 'error' }> } | { kind: 'done'; item: Extract<StreamItem, { kind: 'done' }> };

function groupReply(items: StreamItem[]): ReplyPart[] {
  const parts: ReplyPart[] = [];
  for (const it of items) {
    if (it.kind === 'text') parts.push({ kind: 'text', item: it });
    else if (it.kind === 'tool' || it.kind === 'file') {
      const last = parts[parts.length - 1];
      if (last?.kind === 'tools') last.items.push(it);
      else parts.push({ kind: 'tools', items: [it] });
    } else if (it.kind === 'error') parts.push({ kind: 'error', item: it });
    else if (it.kind === 'done') parts.push({ kind: 'done', item: it });
  }
  return parts;
}

function ToolsGroup({ items }: { items: Extract<StreamItem, { kind: 'tool' | 'file' }>[] }) {
  const [open, setOpen] = useState(false);
  const tools = items.filter((i) => i.kind === 'tool').length;
  const files = items.filter((i) => i.kind === 'file').length;
  const pending = items.some((i) => i.kind === 'tool' && !i.result);
  return (
    <div className="rounded-lg border border-slate-200 bg-slate-50/70 text-xs dark:border-surface-dark-border dark:bg-slate-800/40">
      <button onClick={() => setOpen((o) => !o)} className="flex w-full items-center gap-2 px-3 py-1.5 text-left text-steel-600 dark:text-steel-300">
        {open ? <ChevronDown className="h-3.5 w-3.5" /> : <ChevronRight className="h-3.5 w-3.5" />}
        <Wrench className="h-3.5 w-3.5" />
        {tools ? `${tools} 次工具调用` : ''}
        {tools && files ? ' · ' : ''}
        {files ? `${files} 个文件更新` : ''}
        {pending ? <Loader2 className="ml-1 h-3 w-3 animate-spin text-emerald-500" /> : null}
      </button>
      {open ? (
        <ul className="space-y-1 border-t border-slate-200 px-3 py-2 dark:border-surface-dark-border">
          {items.map((it) => (
            <li key={it.id} className="flex items-center gap-2 font-mono text-[11px] text-slate-700 dark:text-slate-300">
              {it.kind === 'tool' ? (
                <>
                  {!it.result ? <Loader2 className="h-3 w-3 animate-spin text-emerald-500" /> : it.result.isError ? <XCircle className="h-3 w-3 text-red-500" /> : <CheckCircle2 className="h-3 w-3 text-emerald-500" />}
                  <span className="font-semibold">{it.name.replace(/^mcp__bizgold__/, '')}</span>
                  <span className="min-w-0 flex-1 truncate text-steel-500">{summarizeInput(it.input, 90)}</span>
                </>
              ) : (
                <>
                  <FileEdit className="h-3 w-3 text-brand-500" />
                  <span className="font-semibold">{it.op === 'write' ? '写入' : it.op === 'edit' ? '编辑' : '删除'}</span>
                  <span className="min-w-0 flex-1 truncate text-steel-500">{it.path}</span>
                  {it.size != null ? <span className="text-steel-400">{formatBytes(it.size)}</span> : null}
                </>
              )}
              <span className="tabular text-steel-400">{formatTime(it.ts)}</span>
            </li>
          ))}
        </ul>
      ) : null}
    </div>
  );
}

export function ChatTab({ project, state, onProjectChanged }: { project: Project; state: ConsoleState; onProjectChanged: () => void }) {
  const key = `bizgold.chat.${project.id}`;
  const [messages, setMessages] = useState<ChatMessage[]>(() => loadJson<ChatMessage[]>(key, []));
  const [draft, setDraft] = useState('');
  const [sending, setSending] = useState(false);
  const toast = useToast();
  const listRef = useRef<HTMLDivElement>(null);

  useEffect(() => setMessages(loadJson<ChatMessage[]>(key, [])), [key]);
  useEffect(() => saveJson(key, messages), [key, messages]);

  const running = project.status === 'running';
  const modifying = running && project.session?.mode === 'modify';
  const disabled = running || sending;

  const threads = useMemo(() => {
    return messages.map((m, i) => {
      const next = messages[i + 1];
      const items = state.items.filter((it) => it.seq > m.afterSeq && (!next || it.seq <= next.afterSeq) && it.kind !== 'status' && it.kind !== 'stage' && it.kind !== 'subagent' && it.kind !== 'thinking' && it.kind !== 'review');
      const done = items.some((it) => it.kind === 'done');
      return { message: m, parts: groupReply(items), done, pending: !done && i === messages.length - 1 && (modifying || sending || items.length === 0) };
    });
  }, [messages, state.items, modifying, sending]);

  useEffect(() => {
    const el = listRef.current;
    if (el) el.scrollTop = el.scrollHeight;
  }, [threads]);

  const send = async (text: string) => {
    const t = text.trim();
    if (!t || disabled) return;
    setSending(true);
    try {
      await api.chat(project.id, t);
      setMessages((prev) => [...prev, { id: `${Date.now()}`, text: t, at: new Date().toISOString(), afterSeq: state.lastSeq }]);
      setDraft('');
      onProjectChanged();
    } catch {
      // 错误已由 api 统一 toast
    } finally {
      setSending(false);
    }
  };

  const clear = () => {
    setMessages([]);
    toast.info('已清空本地对话记录', '工作区内的修改不受影响');
  };

  return (
    <div className="grid gap-4 xl:grid-cols-[minmax(0,1fr)_300px]">
      <Card className="flex h-[calc(100vh-22rem)] min-h-[460px] flex-col overflow-hidden">
        <div className="flex items-center gap-3 border-b border-slate-100 px-4 py-2.5 dark:border-surface-dark-border">
          <MessageSquareText className="h-4 w-4 text-brand-500 dark:text-brand-200" />
          <h3 className="text-sm font-semibold text-slate-800 dark:text-slate-100">修改对话</h3>
          <span className="text-xs text-steel-500">用自然语言指挥主代理在同一会话中修改产出，自动重新渲染并自检</span>
          {messages.length ? (
            <Button size="xs" variant="ghost" icon={<Trash2 />} className="ml-auto" onClick={clear}>
              清空记录
            </Button>
          ) : null}
        </div>

        <div ref={listRef} className="min-h-0 flex-1 overflow-y-auto bg-slate-50/50 px-4 py-4 dark:bg-surface-dark/40">
          {threads.length === 0 ? (
            <div className="flex h-full flex-col items-center justify-center text-center">
              <div className="mb-3 inline-flex h-14 w-14 items-center justify-center rounded-2xl bg-gradient-to-br from-brand-50 to-gold-50 text-brand-500 ring-1 ring-brand-100 dark:from-brand-900/60 dark:to-gold-900/20 dark:text-brand-200 dark:ring-brand-800">
                <Sparkles className="h-7 w-7" />
              </div>
              <h4 className="text-sm font-semibold text-slate-800 dark:text-slate-100">还没有修改指令</h4>
              <p className="mt-1 max-w-md text-xs leading-relaxed text-steel-500">
                项目完成（或停止）后，可以在这里下达修改指令，例如补充章节、统一口径、压缩讲稿。主代理会沿用同一会话上下文，只改需要改的地方，并重新渲染 Word / PPT。
              </p>
              <div className="mt-4 flex flex-wrap justify-center gap-2">
                {SAMPLES.map((s) => (
                  <button key={s} onClick={() => setDraft(s)} className="chip cursor-pointer py-1 hover:border-brand-300 hover:text-brand-700 dark:hover:border-brand-400">
                    {s}
                  </button>
                ))}
              </div>
            </div>
          ) : (
            <div className="space-y-5">
              {threads.map(({ message, parts, pending, done }) => (
                <div key={message.id} className="space-y-3">
                  <div className="flex justify-end">
                    <div className="max-w-[80%]">
                      <div className="rounded-2xl rounded-tr-sm bg-brand-600 px-4 py-2.5 text-sm leading-relaxed text-white shadow-sm dark:bg-brand-500">{message.text}</div>
                      <div className="mt-1 flex items-center justify-end gap-1 text-[11px] text-steel-400">
                        <User className="h-3 w-3" />
                        {formatTime(message.at)}
                      </div>
                    </div>
                  </div>
                  <div className="flex gap-3">
                    <span className="mt-1 inline-flex h-7 w-7 shrink-0 items-center justify-center rounded-full bg-gold-100 text-gold-800 dark:bg-gold-900/40 dark:text-gold-200">
                      <Bot className="h-4 w-4" />
                    </span>
                    <div className="min-w-0 max-w-[85%] flex-1 space-y-2">
                      {parts.map((p, i) => {
                        if (p.kind === 'text') {
                          return (
                            <div key={p.item.id} className="rounded-2xl rounded-tl-sm border border-slate-200 bg-white px-4 py-2.5 shadow-sm dark:border-surface-dark-border dark:bg-surface-dark-card">
                              <Markdown content={p.item.text} variant="stream" projectId={project.id} />
                              {p.item.streaming ? <span className="ml-0.5 inline-block h-3.5 w-1.5 animate-pulse bg-gold-500 align-text-bottom" /> : null}
                            </div>
                          );
                        }
                        if (p.kind === 'tools') return <ToolsGroup key={`tools-${i}`} items={p.items} />;
                        if (p.kind === 'error') {
                          return (
                            <div key={p.item.id} className="rounded-lg border border-red-200 bg-red-50 px-3 py-2 text-xs text-red-700 dark:border-red-900/60 dark:bg-red-900/20 dark:text-red-200">
                              {p.item.message}
                            </div>
                          );
                        }
                        return (
                          <div key={p.item.id} className={cn('inline-flex items-center gap-1.5 text-[11px]', p.item.subtype === 'success' ? 'text-emerald-600 dark:text-emerald-300' : 'text-red-600 dark:text-red-300')}>
                            {p.item.subtype === 'success' ? <CheckCircle2 className="h-3.5 w-3.5" /> : <XCircle className="h-3.5 w-3.5" />}
                            {p.item.subtype === 'success' ? '修改完成' : p.item.subtype === 'stopped' ? '已停止' : '修改失败'} · {formatTime(p.item.ts)}
                          </div>
                        );
                      })}
                      {pending && !done ? (
                        <div className="inline-flex items-center gap-2 rounded-2xl rounded-tl-sm border border-slate-200 bg-white px-4 py-2.5 text-xs text-steel-500 dark:border-surface-dark-border dark:bg-surface-dark-card">
                          <span className="flex gap-1">
                            <span className="h-1.5 w-1.5 animate-bounce rounded-full bg-gold-500 [animation-delay:-0.3s]" />
                            <span className="h-1.5 w-1.5 animate-bounce rounded-full bg-gold-500 [animation-delay:-0.15s]" />
                            <span className="h-1.5 w-1.5 animate-bounce rounded-full bg-gold-500" />
                          </span>
                          主代理正在处理…
                        </div>
                      ) : null}
                    </div>
                  </div>
                </div>
              ))}
            </div>
          )}
        </div>

        <div className="border-t border-slate-100 px-4 py-3 dark:border-surface-dark-border">
          {running ? (
            <div className="mb-2 rounded-md bg-amber-50 px-3 py-1.5 text-xs text-amber-800 dark:bg-amber-900/20 dark:text-amber-200">
              {modifying ? '正在执行上一条修改指令，请稍候…' : '项目正在运行中，完成或停止后才能发送修改指令。'}
            </div>
          ) : null}
          <div className="flex items-end gap-2">
            <Textarea
              value={draft}
              onChange={(e) => setDraft(e.target.value)}
              disabled={disabled}
              rows={2}
              placeholder="输入修改指令，Ctrl / ⌘ + Enter 发送"
              className="min-h-[56px]"
              onKeyDown={(e) => {
                if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') void send(draft);
              }}
            />
            <Button variant="gold" icon={<Send />} loading={sending} disabled={disabled || !draft.trim()} onClick={() => void send(draft)} className="h-[56px]">
              发送
            </Button>
          </div>
        </div>
      </Card>

      <div className="space-y-4">
        <Card className="p-4">
          <h4 className="text-sm font-semibold text-slate-800 dark:text-slate-100">示例指令</h4>
          <div className="mt-3 space-y-2">
            {SAMPLES.map((s) => (
              <button key={s} onClick={() => setDraft(s)} disabled={disabled} className="w-full rounded-lg border border-slate-200 px-3 py-2 text-left text-xs text-slate-700 transition hover:border-brand-300 hover:bg-brand-50/50 disabled:opacity-50 dark:border-surface-dark-border dark:text-slate-200 dark:hover:border-brand-400">
                {s}
              </button>
            ))}
          </div>
        </Card>
        <Card className="p-4 text-xs leading-relaxed text-steel-600 dark:text-steel-300">
          <h4 className="mb-2 text-sm font-semibold text-slate-800 dark:text-slate-100">修改模式说明</h4>
          <ul className="list-disc space-y-1 pl-4">
            <li>沿用上次运行的会话上下文（session），主代理知道全部产出物。</li>
            <li>只修改指令涉及的章节 / 页面，数字与数字账本保持一致。</li>
            <li>修改后自动重新渲染 Word / PPT，并运行章节自检。</li>
            <li>每条指令都会在「运行控制台」留下完整事件记录。</li>
          </ul>
        </Card>
      </div>
    </div>
  );
}
