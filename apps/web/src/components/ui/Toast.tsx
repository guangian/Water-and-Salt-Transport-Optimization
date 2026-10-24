import { createContext, useCallback, useContext, useEffect, useMemo, useRef, useState, type ReactNode } from 'react';
import { AlertTriangle, CheckCircle2, Info, X, XCircle } from 'lucide-react';
import { onApiError } from '../../api';
import { cn } from '../../lib/cn';

export type ToastKind = 'success' | 'error' | 'info' | 'warning';

export interface ToastItem {
  id: number;
  kind: ToastKind;
  title: string;
  description?: string;
}

interface ToastApi {
  push: (t: Omit<ToastItem, 'id'>) => void;
  success: (title: string, description?: string) => void;
  error: (title: string, description?: string) => void;
  info: (title: string, description?: string) => void;
  warning: (title: string, description?: string) => void;
}

const ToastContext = createContext<ToastApi | null>(null);

const ICONS: Record<ToastKind, ReactNode> = {
  success: <CheckCircle2 className="h-5 w-5 text-emerald-500" />,
  error: <XCircle className="h-5 w-5 text-red-500" />,
  info: <Info className="h-5 w-5 text-brand-400" />,
  warning: <AlertTriangle className="h-5 w-5 text-amber-500" />,
};

export function ToastProvider({ children }: { children: ReactNode }) {
  const [items, setItems] = useState<ToastItem[]>([]);
  const counter = useRef(0);

  const dismiss = useCallback((id: number) => setItems((xs) => xs.filter((x) => x.id !== id)), []);

  const push = useCallback(
    (t: Omit<ToastItem, 'id'>) => {
      counter.current += 1;
      const id = counter.current;
      setItems((xs) => [...xs.slice(-4), { ...t, id }]);
      window.setTimeout(() => dismiss(id), t.kind === 'error' ? 7000 : 4000);
    },
    [dismiss],
  );

  const api = useMemo<ToastApi>(
    () => ({
      push,
      success: (title, description) => push({ kind: 'success', title, description }),
      error: (title, description) => push({ kind: 'error', title, description }),
      info: (title, description) => push({ kind: 'info', title, description }),
      warning: (title, description) => push({ kind: 'warning', title, description }),
    }),
    [push],
  );

  useEffect(() => {
    let last = 0;
    let lastMsg = '';
    return onApiError((err) => {
      const now = Date.now();
      if (err.message === lastMsg && now - last < 2500) return;
      last = now;
      lastMsg = err.message;
      api.error('请求失败', err.message);
    });
  }, [api]);

  return (
    <ToastContext.Provider value={api}>
      {children}
      <div className="pointer-events-none fixed bottom-5 right-5 z-[100] flex w-[360px] max-w-[calc(100vw-2rem)] flex-col gap-2">
        {items.map((t) => (
          <div
            key={t.id}
            className={cn(
              'pointer-events-auto flex animate-fade-in items-start gap-3 rounded-xl border bg-white/95 p-3.5 shadow-card-hover backdrop-blur dark:bg-surface-dark-card/95',
              t.kind === 'error' ? 'border-red-200 dark:border-red-900/60' : t.kind === 'success' ? 'border-emerald-200 dark:border-emerald-900/60' : 'border-slate-200 dark:border-surface-dark-border',
            )}
            role="status"
          >
            <span className="mt-0.5 shrink-0">{ICONS[t.kind]}</span>
            <div className="min-w-0 flex-1">
              <p className="text-sm font-medium text-slate-800 dark:text-slate-100">{t.title}</p>
              {t.description ? <p className="mt-0.5 break-words text-xs leading-relaxed text-steel-500 dark:text-steel-300">{t.description}</p> : null}
            </div>
            <button onClick={() => dismiss(t.id)} className="rounded-md p-1 text-slate-400 hover:bg-slate-100 hover:text-slate-600 dark:hover:bg-slate-800" aria-label="关闭">
              <X className="h-4 w-4" />
            </button>
          </div>
        ))}
      </div>
    </ToastContext.Provider>
  );
}

export function useToast(): ToastApi {
  const ctx = useContext(ToastContext);
  if (!ctx) throw new Error('useToast 必须在 ToastProvider 内使用');
  return ctx;
}
