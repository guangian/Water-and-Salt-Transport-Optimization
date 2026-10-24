import { useCallback, useEffect, useRef, useState } from 'react';
import { toApiError, type ApiError } from '../lib/errors';

export interface AsyncState<T> {
  data: T | null;
  loading: boolean;
  error: ApiError | null;
  reload: () => Promise<void>;
  setData: (updater: T | ((prev: T | null) => T | null)) => void;
}

/** 简单的数据加载 hook：首次加载 + 手动刷新；依赖变化时重新加载。 */
export function useAsync<T>(fn: () => Promise<T>, deps: unknown[], opts: { enabled?: boolean; pollMs?: number } = {}): AsyncState<T> {
  const { enabled = true, pollMs } = opts;
  const [data, setDataState] = useState<T | null>(null);
  const [loading, setLoading] = useState<boolean>(enabled);
  const [error, setError] = useState<ApiError | null>(null);
  const fnRef = useRef(fn);
  fnRef.current = fn;
  const seq = useRef(0);

  const load = useCallback(async (silent = false) => {
    const id = ++seq.current;
    if (!silent) setLoading(true);
    try {
      const d = await fnRef.current();
      if (id === seq.current) {
        setDataState(d);
        setError(null);
      }
    } catch (e) {
      if (id === seq.current) setError(toApiError(e));
    } finally {
      if (id === seq.current) setLoading(false);
    }
  }, []);

  useEffect(() => {
    if (!enabled) {
      setLoading(false);
      return;
    }
    void load();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [enabled, ...deps]);

  useEffect(() => {
    if (!enabled || !pollMs) return;
    const t = window.setInterval(() => void load(true), pollMs);
    return () => window.clearInterval(t);
  }, [enabled, pollMs, load]);

  const setData = useCallback((updater: T | ((prev: T | null) => T | null)) => {
    setDataState((prev) => (typeof updater === 'function' ? (updater as (p: T | null) => T | null)(prev) : updater));
  }, []);

  return { data, loading, error, reload: () => load(true), setData };
}

export function useInterval(cb: () => void, ms: number | null): void {
  const ref = useRef(cb);
  ref.current = cb;
  useEffect(() => {
    if (ms == null) return;
    const t = window.setInterval(() => ref.current(), ms);
    return () => window.clearInterval(t);
  }, [ms]);
}

export function useNow(ms = 1000, enabled = true): number {
  const [now, setNow] = useState(Date.now());
  useEffect(() => {
    if (!enabled) return;
    const t = window.setInterval(() => setNow(Date.now()), ms);
    return () => window.clearInterval(t);
  }, [ms, enabled]);
  return now;
}
