import type { AgentEvent, AgentRef, ProjectStatus, StageId, StageState, StageStatus } from '../types';

export interface StreamItemBase {
  id: string;
  seq: number;
  ts: string;
  agent: AgentRef;
}

export type StreamItem =
  | (StreamItemBase & { kind: 'text'; text: string; streaming: boolean })
  | (StreamItemBase & { kind: 'thinking'; text: string; streaming: boolean })
  | (StreamItemBase & {
      kind: 'tool';
      toolUseId: string;
      name: string;
      input: unknown;
      result?: { isError: boolean; output: string; ts: string };
    })
  | (StreamItemBase & { kind: 'file'; path: string; op: 'write' | 'edit' | 'delete'; size?: number })
  | (StreamItemBase & { kind: 'stage'; stage: StageId; status: StageStatus; note?: string; round?: number })
  | (StreamItemBase & { kind: 'status'; status: ProjectStatus; message?: string })
  | (StreamItemBase & { kind: 'subagent'; name: string; status: 'start' | 'end'; summary?: string; toolUseId: string })
  | (StreamItemBase & { kind: 'error'; message: string; recoverable: boolean })
  | (StreamItemBase & { kind: 'review'; round: number; score: number; hardGatePassed: boolean })
  | (StreamItemBase & { kind: 'done'; subtype: 'success' | 'error' | 'stopped'; result?: string });

export interface AgentLane {
  id: string;
  name: string;
  status: 'running' | 'done' | 'idle';
  lastText: string;
  toolCalls: number;
  textChars: number;
  startedAt: string;
  endedAt?: string;
  lastSeq: number;
  parentToolUseId?: string;
}

export interface UsageCounters {
  inputTokens: number;
  outputTokens: number;
  cacheReadTokens: number;
  turns: number;
  costEstimateUsd: number;
}

export interface ConsoleState {
  items: StreamItem[];
  agents: AgentLane[];
  usage: UsageCounters;
  lastSeq: number;
  firstTs: string | null;
  lastTs: string | null;
  stageUpdates: Partial<Record<StageId, StageState>>;
  currentStage: StageId | null;
  status: ProjectStatus | null;
  errors: number;
  done: { subtype: 'success' | 'error' | 'stopped'; result?: string } | null;
  /** 每个代理当前“开放中”的文本/思考项 id（用于合并增量） */
  open: Record<string, { textId?: string; thinkingId?: string }>;
  toolIndex: Record<string, string>;
}

export const MAX_ITEMS = 2500;

export function initialConsoleState(): ConsoleState {
  return {
    items: [],
    agents: [],
    usage: { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, turns: 0, costEstimateUsd: 0 },
    lastSeq: 0,
    firstTs: null,
    lastTs: null,
    stageUpdates: {},
    currentStage: null,
    status: null,
    errors: 0,
    done: null,
    open: {},
    toolIndex: {},
  };
}

function itemId(ev: AgentEvent, suffix = ''): string {
  return `${ev.seq}${suffix}`;
}

function ensureLane(agents: AgentLane[], ref: AgentRef, ev: AgentEvent): AgentLane[] {
  const idx = agents.findIndex((a) => a.id === ref.id);
  if (idx >= 0) {
    const lane = agents[idx];
    const next = { ...lane, lastSeq: ev.seq };
    if (ref.name && ref.name !== lane.name && lane.name === lane.id) next.name = ref.name;
    const copy = agents.slice();
    copy[idx] = next;
    return copy;
  }
  const lane: AgentLane = {
    id: ref.id,
    name: ref.name || ref.id,
    status: 'running',
    lastText: '',
    toolCalls: 0,
    textChars: 0,
    startedAt: ev.ts,
    lastSeq: ev.seq,
  };
  // 主代理永远排第一
  if (ref.id === 'main') return [lane, ...agents];
  return [...agents, lane];
}

function updateLane(agents: AgentLane[], id: string, patch: (l: AgentLane) => AgentLane): AgentLane[] {
  const idx = agents.findIndex((a) => a.id === id);
  if (idx < 0) return agents;
  const copy = agents.slice();
  copy[idx] = patch(copy[idx]);
  return copy;
}

function lastSentence(text: string): string {
  const clean = text.replace(/[#*`>|_-]+/g, ' ').replace(/\s+/g, ' ').trim();
  if (!clean) return '';
  const parts = clean.split(/(?<=[。！？!?；;])/).map((s) => s.trim()).filter(Boolean);
  const last = parts.length ? parts[parts.length - 1] : clean;
  return last.length > 60 ? `${last.slice(-60)}` : last;
}

function pushItem(items: StreamItem[], item: StreamItem): StreamItem[] {
  const next = items.length >= MAX_ITEMS ? items.slice(items.length - MAX_ITEMS + 1) : items.slice();
  next.push(item);
  return next;
}

function replaceItem(items: StreamItem[], id: string, patch: (it: StreamItem) => StreamItem): StreamItem[] {
  // 从尾部查找（最近的项大概率在末尾）
  for (let i = items.length - 1; i >= 0; i--) {
    if (items[i].id === id) {
      const copy = items.slice();
      copy[i] = patch(copy[i]);
      return copy;
    }
  }
  return items;
}

function closeOpen(state: ConsoleState, agentId: string): ConsoleState {
  const open = state.open[agentId];
  if (!open) return state;
  let items = state.items;
  if (open.textId) items = replaceItem(items, open.textId, (it) => (it.kind === 'text' ? { ...it, streaming: false } : it));
  if (open.thinkingId)
    items = replaceItem(items, open.thinkingId, (it) => (it.kind === 'thinking' ? { ...it, streaming: false } : it));
  const nextOpen = { ...state.open };
  delete nextOpen[agentId];
  return { ...state, items, open: nextOpen };
}

/** 事件归约：把 SSE 事件流折叠为可渲染的流项、代理泳道与计数器。 */
export function reduceEvent(prev: ConsoleState, ev: AgentEvent): ConsoleState {
  if (ev.seq <= prev.lastSeq && prev.lastSeq !== 0 && ev.seq !== 0) {
    // 重复/乱序事件（断线重连回放）直接忽略
    return prev;
  }
  const agentId = ev.agent?.id ?? 'main';
  const agentRef: AgentRef = ev.agent ?? { id: 'main', name: '主代理' };
  let state: ConsoleState = {
    ...prev,
    lastSeq: Math.max(prev.lastSeq, ev.seq),
    firstTs: prev.firstTs ?? ev.ts,
    lastTs: ev.ts,
    agents: ensureLane(prev.agents, agentRef, ev),
  };

  switch (ev.type) {
    case 'text_delta': {
      const open = state.open[agentId];
      if (open?.textId) {
        state = {
          ...state,
          items: replaceItem(state.items, open.textId, (it) =>
            it.kind === 'text' ? { ...it, text: it.text + ev.data.text, ts: ev.ts } : it,
          ),
        };
      } else {
        const id = itemId(ev);
        state = {
          ...state,
          items: pushItem(state.items, { kind: 'text', id, seq: ev.seq, ts: ev.ts, agent: agentRef, text: ev.data.text, streaming: true }),
          open: { ...state.open, [agentId]: { ...state.open[agentId], textId: id } },
        };
      }
      state = {
        ...state,
        agents: updateLane(state.agents, agentId, (l) => ({
          ...l,
          status: 'running',
          textChars: l.textChars + ev.data.text.length,
          lastText: lastSentence((l.lastText + ev.data.text).slice(-200)),
        })),
      };
      return state;
    }
    case 'text': {
      const open = state.open[agentId];
      if (open?.textId) {
        // 完整文本块到达：覆盖增量拼接结果并收尾
        state = {
          ...state,
          items: replaceItem(state.items, open.textId, (it) =>
            it.kind === 'text' ? { ...it, text: ev.data.text || it.text, streaming: false } : it,
          ),
          open: { ...state.open, [agentId]: { ...state.open[agentId], textId: undefined } },
        };
      } else {
        state = {
          ...state,
          items: pushItem(state.items, { kind: 'text', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, text: ev.data.text, streaming: false }),
        };
      }
      state = {
        ...state,
        agents: updateLane(state.agents, agentId, (l) => ({
          ...l,
          status: 'running',
          textChars: open?.textId ? l.textChars : l.textChars + ev.data.text.length,
          lastText: lastSentence(ev.data.text) || l.lastText,
        })),
      };
      return state;
    }
    case 'thinking': {
      const open = state.open[agentId];
      if (ev.data.delta && open?.thinkingId) {
        return {
          ...state,
          items: replaceItem(state.items, open.thinkingId, (it) =>
            it.kind === 'thinking' ? { ...it, text: it.text + ev.data.text, ts: ev.ts } : it,
          ),
        };
      }
      // 文本流被思考打断：先收尾文本
      if (open?.textId) state = closeOpen(state, agentId);
      const id = itemId(ev);
      return {
        ...state,
        items: pushItem(state.items, { kind: 'thinking', id, seq: ev.seq, ts: ev.ts, agent: agentRef, text: ev.data.text, streaming: ev.data.delta }),
        open: { ...state.open, [agentId]: { textId: undefined, thinkingId: ev.data.delta ? id : undefined } },
      };
    }
    case 'tool_use': {
      state = closeOpen(state, agentId);
      const id = itemId(ev);
      return {
        ...state,
        items: pushItem(state.items, {
          kind: 'tool',
          id,
          seq: ev.seq,
          ts: ev.ts,
          agent: agentRef,
          toolUseId: ev.data.toolUseId,
          name: ev.data.name,
          input: ev.data.input,
        }),
        toolIndex: { ...state.toolIndex, [ev.data.toolUseId]: id },
        agents: updateLane(state.agents, agentId, (l) => ({ ...l, status: 'running', toolCalls: l.toolCalls + 1 })),
      };
    }
    case 'tool_result': {
      const id = state.toolIndex[ev.data.toolUseId];
      if (!id) {
        return {
          ...state,
          items: pushItem(state.items, {
            kind: 'tool',
            id: itemId(ev),
            seq: ev.seq,
            ts: ev.ts,
            agent: agentRef,
            toolUseId: ev.data.toolUseId,
            name: ev.data.name,
            input: undefined,
            result: { isError: ev.data.isError, output: ev.data.output, ts: ev.ts },
          }),
        };
      }
      return {
        ...state,
        items: replaceItem(state.items, id, (it) =>
          it.kind === 'tool' ? { ...it, result: { isError: ev.data.isError, output: ev.data.output, ts: ev.ts } } : it,
        ),
      };
    }
    case 'file': {
      return {
        ...state,
        items: pushItem(state.items, { kind: 'file', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, path: ev.data.path, op: ev.data.op, size: ev.data.size }),
      };
    }
    case 'stage': {
      state = closeOpen(state, agentId);
      const prevStage: Partial<StageState> = state.stageUpdates[ev.data.stage] ?? {};
      const patch: StageState = {
        ...prevStage,
        status: ev.data.status,
        note: ev.data.note ?? prevStage.note,
        round: ev.data.round ?? prevStage.round,
      };
      if (ev.data.status === 'running') patch.startedAt = prevStage.startedAt ?? ev.ts;
      if (ev.data.status === 'done' || ev.data.status === 'failed' || ev.data.status === 'skipped') patch.endedAt = ev.ts;
      return {
        ...state,
        items: pushItem(state.items, { kind: 'stage', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, stage: ev.data.stage, status: ev.data.status, note: ev.data.note, round: ev.data.round }),
        stageUpdates: { ...state.stageUpdates, [ev.data.stage]: patch },
        currentStage: ev.data.status === 'running' ? ev.data.stage : state.currentStage,
      };
    }
    case 'status': {
      return {
        ...state,
        status: ev.data.status,
        items: pushItem(state.items, { kind: 'status', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, status: ev.data.status, message: ev.data.message }),
      };
    }
    case 'subagent': {
      const laneId = ev.data.toolUseId;
      let agents = state.agents;
      if (ev.data.status === 'start') {
        const existing = agents.find((a) => a.id === laneId);
        if (!existing) {
          agents = [...agents, { id: laneId, name: ev.data.name, status: 'running', lastText: '', toolCalls: 0, textChars: 0, startedAt: ev.ts, lastSeq: ev.seq, parentToolUseId: ev.data.toolUseId }];
        } else {
          agents = updateLane(agents, laneId, (l) => ({ ...l, name: ev.data.name || l.name, status: 'running' }));
        }
      } else {
        agents = updateLane(agents, laneId, (l) => ({ ...l, status: 'done', endedAt: ev.ts, lastText: ev.data.summary ? lastSentence(ev.data.summary) : l.lastText }));
      }
      return {
        ...state,
        agents,
        items: pushItem(state.items, { kind: 'subagent', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, name: ev.data.name, status: ev.data.status, summary: ev.data.summary, toolUseId: ev.data.toolUseId }),
      };
    }
    case 'usage': {
      return {
        ...state,
        usage: {
          inputTokens: ev.data.inputTokens,
          outputTokens: ev.data.outputTokens,
          cacheReadTokens: ev.data.cacheReadTokens,
          turns: ev.data.turns,
          costEstimateUsd: ev.data.costEstimateUsd,
        },
      };
    }
    case 'review': {
      return {
        ...state,
        items: pushItem(state.items, { kind: 'review', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, round: ev.data.round, score: ev.data.score, hardGatePassed: ev.data.hardGatePassed }),
      };
    }
    case 'error': {
      return {
        ...state,
        errors: state.errors + 1,
        items: pushItem(state.items, { kind: 'error', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, message: ev.data.message, recoverable: ev.data.recoverable }),
      };
    }
    case 'done': {
      // 关闭所有开放文本
      let s = state;
      for (const id of Object.keys(s.open)) s = closeOpen(s, id);
      const agents = s.agents.map((a) => (a.status === 'running' ? { ...a, status: 'done' as const, endedAt: ev.ts } : a));
      return {
        ...s,
        agents,
        done: { subtype: ev.data.subtype, result: ev.data.result },
        items: pushItem(s.items, { kind: 'done', id: itemId(ev), seq: ev.seq, ts: ev.ts, agent: agentRef, subtype: ev.data.subtype, result: ev.data.result }),
      };
    }
    default:
      return state;
  }
}

export function reduceEvents(prev: ConsoleState, events: AgentEvent[]): ConsoleState {
  let s = prev;
  for (const ev of events) s = reduceEvent(s, ev);
  return s;
}

export type StreamFilter = 'all' | 'text' | 'tools';

export function filterItems(items: StreamItem[], filter: StreamFilter, agentId: string | 'all'): StreamItem[] {
  return items.filter((it) => {
    if (agentId !== 'all' && it.agent.id !== agentId && it.kind !== 'stage' && it.kind !== 'error') return false;
    if (filter === 'text') return it.kind === 'text' || it.kind === 'thinking' || it.kind === 'stage' || it.kind === 'error';
    if (filter === 'tools') return it.kind === 'tool' || it.kind === 'file' || it.kind === 'stage' || it.kind === 'error';
    return true;
  });
}
