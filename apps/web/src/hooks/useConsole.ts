import { useCallback, useEffect, useReducer, useRef, useState } from 'react';
import { api, type SseStatus } from '../api';
import { initialConsoleState, reduceEvent, type ConsoleState } from '../lib/events';
import type { AgentEvent } from '../types';

interface ConsoleHook {
  state: ConsoleState;
  sse: SseStatus;
  retryInMs: number | null;
  reset: () => void;
  /** 最近一次收到事件的时间戳（ms） */
  lastEventAt: number | null;
}

type Action = { type: 'events'; events: AgentEvent[] } | { type: 'reset' };

function reducer(state: ConsoleState, action: Action): ConsoleState {
  if (action.type === 'reset') return initialConsoleState();
  let s = state;
  for (const ev of action.events) s = reduceEvent(s, ev);
  return s;
}

/**
 * 订阅项目事件流（SSE / Mock），归约为控制台状态。
 * - 事件以 ~60ms 批量合并进 reducer，避免高频渲染；
 * - 断线后自动以 since=lastSeq 续传（由 api.subscribeEvents 处理重连）。
 */
export function useConsole(projectId: string | undefined, onEvent?: (ev: AgentEvent) => void): ConsoleHook {
  const [state, dispatch] = useReducer(reducer, undefined, initialConsoleState);
  const [sse, setSse] = useState<SseStatus>('connecting');
  const [retryInMs, setRetryInMs] = useState<number | null>(null);
  const [lastEventAt, setLastEventAt] = useState<number | null>(null);
  const queue = useRef<AgentEvent[]>([]);
  const flushTimer = useRef<number | null>(null);
  const onEventRef = useRef(onEvent);
  onEventRef.current = onEvent;
  const [epoch, setEpoch] = useState(0);

  const flush = useCallback(() => {
    flushTimer.current = null;
    if (!queue.current.length) return;
    const batch = queue.current;
    queue.current = [];
    dispatch({ type: 'events', events: batch });
    setLastEventAt(Date.now());
  }, []);

  useEffect(() => {
    if (!projectId) return;
    dispatch({ type: 'reset' });
    queue.current = [];
    const unsubscribe = api.subscribeEvents(projectId, 0, {
      onEvent: (ev) => {
        queue.current.push(ev);
        onEventRef.current?.(ev);
        if (flushTimer.current == null) flushTimer.current = window.setTimeout(flush, 60);
      },
      onStatus: (s, info) => {
        setSse(s);
        setRetryInMs(s === 'reconnecting' ? info?.retryInMs ?? null : null);
      },
    });
    return () => {
      unsubscribe();
      if (flushTimer.current != null) window.clearTimeout(flushTimer.current);
      flushTimer.current = null;
    };
  }, [projectId, epoch, flush]);

  const reset = useCallback(() => setEpoch((e) => e + 1), []);

  return { state, sse, retryInMs, reset, lastEventAt };
}
