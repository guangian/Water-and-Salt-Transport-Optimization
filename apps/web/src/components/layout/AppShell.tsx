import { type ReactNode } from 'react';
import { NavLink, useLocation } from 'react-router-dom';
import { BookOpen, FlaskConical, LayoutDashboard, Moon, PlusCircle, Settings, Sun } from 'lucide-react';
import { Logo } from '../domain/Brand';
import { useTheme } from '../../hooks/useTheme';
import { useAsync } from '../../hooks/useAsync';
import { api, IS_MOCK } from '../../api';
import { cn } from '../../lib/cn';
import type { Health } from '../../types';

const NAV = [
  { to: '/', label: '工作台', icon: <LayoutDashboard />, end: true },
  { to: '/new', label: '新建项目', icon: <PlusCircle /> },
  { to: '/knowledge', label: '知识库', icon: <BookOpen /> },
  { to: '/settings', label: '设置', icon: <Settings /> },
];

function HealthDot({ ok, label, detail }: { ok: 'ok' | 'warn' | 'bad' | 'unknown'; label: string; detail: string }) {
  const color = ok === 'ok' ? 'bg-emerald-500' : ok === 'warn' ? 'bg-amber-500' : ok === 'bad' ? 'bg-red-500' : 'bg-slate-300 dark:bg-slate-600';
  return (
    <span className="group relative inline-flex items-center gap-1.5 text-xs text-slate-600 dark:text-slate-300" title={detail}>
      <span className={cn('h-2 w-2 rounded-full', color, ok === 'ok' && 'shadow-[0_0_0_3px_rgba(16,185,129,0.15)]')} />
      {label}
    </span>
  );
}

export function HealthDots({ health, loading }: { health: Health | null; loading: boolean }) {
  if (!health) {
    return (
      <div className="flex items-center gap-4">
        <HealthDot ok={loading ? 'unknown' : 'bad'} label="API" detail={loading ? '检测中…' : '无法连接服务端'} />
        <HealthDot ok="unknown" label="Python" detail="未知" />
        <HealthDot ok="unknown" label="PDF 转换器" detail="未知" />
      </div>
    );
  }
  const py = !health.python.found ? 'bad' : health.python.missing.length ? 'warn' : 'ok';
  const pdf = health.pdfConverter ? 'ok' : 'warn';
  return (
    <div className="flex items-center gap-4">
      <HealthDot ok={health.ok ? 'ok' : 'bad'} label="API" detail={`服务端 v${health.version}${IS_MOCK ? '（Mock）' : ''}`} />
      <HealthDot ok={py} label="Python" detail={health.python.found ? `${health.python.version} · ${health.python.path}${health.python.missing.length ? ` · 缺少 ${health.python.missing.join(', ')}` : ''}` : '未找到 Python'} />
      <HealthDot ok={pdf} label="PDF 转换器" detail={health.pdfConverter ? `使用 ${health.pdfConverter}` : '未安装 Word/WPS/LibreOffice，无法自动生成 PDF'} />
    </div>
  );
}

export function AppShell({ children }: { children: ReactNode }) {
  const [theme, , toggle] = useTheme();
  const { data: health, loading } = useAsync(() => api.health({ silent: true }), [], { pollMs: 30_000 });
  const location = useLocation();

  return (
    <div className="flex h-full min-h-screen">
      <aside className="fixed inset-y-0 left-0 z-40 flex w-16 flex-col bg-gradient-to-b from-brand-700 via-brand-700 to-brand-900 text-white xl:w-60">
        <div className="flex h-16 items-center gap-3 px-3 xl:px-5">
          <Logo size={36} className="shrink-0 drop-shadow" />
          <div className="hidden min-w-0 xl:block">
            <div className="truncate text-[15px] font-semibold tracking-wide">金策 BizGold</div>
            <div className="truncate text-[11px] text-gold-300/90">商赛国金级智能体工作台</div>
          </div>
        </div>
        <nav className="mt-2 flex flex-1 flex-col gap-1 px-2 xl:px-3">
          {NAV.map((n) => (
            <NavLink key={n.to} to={n.to} end={n.end} className={({ isActive }) => cn('nav-item justify-center xl:justify-start', isActive && 'nav-item-active')} title={n.label}>
              <span className="[&>svg]:h-5 [&>svg]:w-5">{n.icon}</span>
              <span className="hidden xl:inline">{n.label}</span>
            </NavLink>
          ))}
        </nav>
        <div className="hidden px-5 pb-5 xl:block">
          <div className="rounded-xl border border-white/10 bg-white/5 p-3 text-[11px] leading-relaxed text-brand-100/80">
            <div className="mb-1 flex items-center gap-1.5 font-medium text-gold-300">
              <FlaskConical className="h-3.5 w-3.5" />
              {IS_MOCK ? 'Mock 演示模式' : '本地运行'}
            </div>
            {IS_MOCK ? '使用内置假数据与模拟事件流，不调用模型 API。' : '事件与产出保存在 workspace/ 目录。'}
          </div>
        </div>
      </aside>

      <div className="flex min-w-0 flex-1 flex-col pl-16 xl:pl-60">
        <header className="sticky top-0 z-30 flex h-14 items-center justify-between gap-4 border-b border-slate-200/80 bg-white/80 px-6 backdrop-blur dark:border-surface-dark-border dark:bg-surface-dark/80">
          <div className="flex min-w-0 items-center gap-3 text-sm text-steel-500 dark:text-steel-400">
            <span className="hidden sm:inline">环境健康</span>
            <HealthDots health={health} loading={loading} />
          </div>
          <div className="flex items-center gap-2">
            {IS_MOCK ? <span className="chip border-gold-200 bg-gold-50 text-gold-800 dark:border-gold-800 dark:bg-gold-900/30 dark:text-gold-200">MOCK</span> : null}
            <button
              onClick={toggle}
              className="inline-flex h-9 w-9 items-center justify-center rounded-lg border border-slate-200 bg-white text-slate-600 transition hover:bg-slate-50 dark:border-surface-dark-border dark:bg-surface-dark-card dark:text-slate-300 dark:hover:bg-slate-800"
              title={theme === 'dark' ? '切换到浅色' : '切换到深色'}
              aria-label="切换主题"
            >
              {theme === 'dark' ? <Sun className="h-4 w-4" /> : <Moon className="h-4 w-4" />}
            </button>
          </div>
        </header>
        <main key={location.pathname} className="min-w-0 flex-1 animate-fade-in">
          {children}
        </main>
      </div>
    </div>
  );
}

export function PageHeader({ title, description, actions, breadcrumb }: { title: ReactNode; description?: ReactNode; actions?: ReactNode; breadcrumb?: ReactNode }) {
  return (
    <div className="mb-6 flex flex-wrap items-end justify-between gap-4">
      <div className="min-w-0">
        {breadcrumb ? <div className="mb-1 text-xs text-steel-500">{breadcrumb}</div> : null}
        <h1 className="text-2xl font-semibold tracking-tight text-brand-700 dark:text-slate-50">{title}</h1>
        {description ? <p className="mt-1 text-sm text-steel-500 dark:text-steel-400">{description}</p> : null}
      </div>
      {actions ? <div className="flex flex-wrap items-center gap-2">{actions}</div> : null}
    </div>
  );
}

export function Page({ children, className, wide }: { children: ReactNode; className?: string; wide?: boolean }) {
  return <div className={cn('mx-auto w-full px-6 py-6', wide ? 'max-w-[1680px]' : 'max-w-[1440px]', className)}>{children}</div>;
}
