import { forwardRef, type InputHTMLAttributes, type ReactNode, type SelectHTMLAttributes, type TextareaHTMLAttributes } from 'react';
import { ChevronDown } from 'lucide-react';
import { cn } from '../../lib/cn';

export function Field({
  label,
  hint,
  error,
  required,
  children,
  className,
  right,
}: {
  label?: ReactNode;
  hint?: ReactNode;
  error?: string | null;
  required?: boolean;
  children: ReactNode;
  className?: string;
  right?: ReactNode;
}) {
  return (
    <div className={className}>
      {label ? (
        <div className="mb-1.5 flex items-center justify-between">
          <label className="block text-sm font-medium text-slate-700 dark:text-slate-300">
            {label}
            {required ? <span className="ml-0.5 text-red-500">*</span> : null}
          </label>
          {right}
        </div>
      ) : null}
      {children}
      {error ? <p className="mt-1.5 text-xs text-red-600 dark:text-red-400">{error}</p> : hint ? <p className="hint">{hint}</p> : null}
    </div>
  );
}

export const Input = forwardRef<HTMLInputElement, InputHTMLAttributes<HTMLInputElement> & { invalid?: boolean }>(function Input({ className, invalid, ...rest }, ref) {
  return <input ref={ref} className={cn('input', invalid && 'border-red-400 focus:border-red-400 focus:ring-red-100', className)} {...rest} />;
});

export const Textarea = forwardRef<HTMLTextAreaElement, TextareaHTMLAttributes<HTMLTextAreaElement> & { invalid?: boolean }>(function Textarea({ className, invalid, ...rest }, ref) {
  return <textarea ref={ref} className={cn('input min-h-[96px] resize-y leading-relaxed', invalid && 'border-red-400 focus:border-red-400 focus:ring-red-100', className)} {...rest} />;
});

export function Select({ className, children, invalid, ...rest }: SelectHTMLAttributes<HTMLSelectElement> & { invalid?: boolean }) {
  return (
    <div className="relative">
      <select className={cn('input appearance-none pr-9', invalid && 'border-red-400', className)} {...rest}>
        {children}
      </select>
      <ChevronDown className="pointer-events-none absolute right-3 top-1/2 h-4 w-4 -translate-y-1/2 text-slate-400" />
    </div>
  );
}

export function Switch({ checked, onChange, label, disabled }: { checked: boolean; onChange: (v: boolean) => void; label?: ReactNode; disabled?: boolean }) {
  return (
    <button
      type="button"
      role="switch"
      aria-checked={checked}
      disabled={disabled}
      onClick={() => onChange(!checked)}
      className={cn('inline-flex items-center gap-2.5 text-sm text-slate-700 disabled:opacity-50 dark:text-slate-300')}
    >
      <span className={cn('relative inline-flex h-5 w-9 shrink-0 items-center rounded-full transition', checked ? 'bg-brand-600' : 'bg-slate-300 dark:bg-slate-600')}>
        <span className={cn('inline-block h-4 w-4 transform rounded-full bg-white shadow transition', checked ? 'translate-x-4' : 'translate-x-0.5')} />
      </span>
      {label}
    </button>
  );
}

export function Slider({
  value,
  min,
  max,
  step = 1,
  onChange,
  format = (v: number) => String(v),
  marks,
}: {
  value: number;
  min: number;
  max: number;
  step?: number;
  onChange: (v: number) => void;
  format?: (v: number) => string;
  marks?: number[];
}) {
  const pct = ((value - min) / (max - min)) * 100;
  return (
    <div>
      <div className="flex items-center gap-3">
        <input
          type="range"
          min={min}
          max={max}
          step={step}
          value={value}
          onChange={(e) => onChange(Number(e.target.value))}
          className="h-2 w-full cursor-pointer appearance-none rounded-full accent-brand-600"
          style={{ background: `linear-gradient(to right, var(--range-fill) ${pct}%, var(--range-track) ${pct}%)` }}
        />
        <span className="w-16 shrink-0 text-right text-sm font-semibold tabular text-brand-700 dark:text-brand-100">{format(value)}</span>
      </div>
      {marks ? (
        <div className="mt-1 flex justify-between pr-[76px] text-[11px] text-steel-400">
          {marks.map((m) => (
            <span key={m}>{format(m)}</span>
          ))}
        </div>
      ) : null}
    </div>
  );
}

export function Segmented<T extends string | number>({
  value,
  options,
  onChange,
  size = 'md',
  className,
}: {
  value: T;
  options: { value: T; label: ReactNode; icon?: ReactNode; disabled?: boolean }[];
  onChange: (v: T) => void;
  size?: 'sm' | 'md';
  className?: string;
}) {
  return (
    <div className={cn('inline-flex rounded-lg bg-slate-100 p-0.5 dark:bg-slate-800', className)} role="tablist">
      {options.map((o) => (
        <button
          key={String(o.value)}
          type="button"
          disabled={o.disabled}
          onClick={() => onChange(o.value)}
          className={cn(
            'inline-flex items-center gap-1.5 whitespace-nowrap rounded-md font-medium transition disabled:opacity-40',
            size === 'sm' ? 'px-2.5 py-1 text-xs' : 'px-3 py-1.5 text-sm',
            value === o.value ? 'bg-white text-brand-700 shadow-sm dark:bg-surface-dark-card dark:text-brand-100' : 'text-slate-600 hover:text-slate-900 dark:text-slate-300 dark:hover:text-white',
          )}
        >
          {o.icon ? <span className="[&>svg]:h-3.5 [&>svg]:w-3.5">{o.icon}</span> : null}
          {o.label}
        </button>
      ))}
    </div>
  );
}
