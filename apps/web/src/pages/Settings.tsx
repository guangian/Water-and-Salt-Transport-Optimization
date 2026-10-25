import { useEffect, useMemo, useState } from 'react';
import { CheckCircle2, Cpu, Eye, EyeOff, FolderOpen, Globe, KeyRound, Palette, PlugZap, Save, Search, SlidersHorizontal, XCircle, Wifi } from 'lucide-react';
import { api } from '../api';
import { useAsync } from '../hooks/useAsync';
import { useTheme } from '../hooks/useTheme';
import { Page, PageHeader } from '../components/layout/AppShell';
import { Card, CardBody, CardHeader } from '../components/ui/Card';
import { Button } from '../components/ui/Button';
import { Badge } from '../components/ui/Badge';
import { Field, Input, Select, Switch } from '../components/ui/Form';
import { Skeleton } from '../components/ui/Misc';
import { useToast } from '../components/ui/Toast';
import type { ModelRole, SearchProvider, Settings, SettingsTestResult } from '../types';
import { cn } from '../lib/cn';

const ROLES: { role: ModelRole; label: string; desc: string }[] = [
  { role: 'main', label: '主代理', desc: '流程调度与读题' },
  { role: 'writer', label: '主笔', desc: '各章撰写' },
  { role: 'pitch', label: '路演设计师', desc: 'slides / 讲稿 / 题库' },
  { role: 'researcher', label: '行业研究员', desc: '检索归纳' },
  { role: 'finance', label: '财务建模师', desc: '模型与账本' },
  { role: 'judgeA', label: '评委 A（投资人）', desc: '多模型评委' },
  { role: 'judgeB', label: '评委 B（学者）', desc: '多模型评委' },
  { role: 'judgeC', label: '评委 C（企业高管）', desc: '多模型评委' },
  { role: 'inspector', label: '质检员', desc: '机器检查解释' },
  { role: 'fast', label: '快速模型', desc: '后台辅助调用' },
];

const PROVIDERS: { id: SearchProvider; label: string; keyField?: 'bochaKey' | 'tavilyKey' | 'serperKey'; note: string }[] = [
  { id: 'bocha', label: '博查 Bocha', keyField: 'bochaKey', note: '中文检索质量高，推荐' },
  { id: 'tavily', label: 'Tavily', keyField: 'tavilyKey', note: '英文资料与学术' },
  { id: 'serper', label: 'Serper（Google）', keyField: 'serperKey', note: 'Google 结果' },
  { id: 'bing', label: 'Bing 网页解析', note: '无需 key，稳定性一般' },
  { id: 'duckduckgo', label: 'DuckDuckGo lite', note: '无需 key，兜底' },
];

function SecretInput({ value, onChange, placeholder }: { value: string; onChange: (v: string) => void; placeholder?: string }) {
  const [show, setShow] = useState(false);
  return (
    <div className="relative">
      <Input type={show ? 'text' : 'password'} value={value} onChange={(e) => onChange(e.target.value)} placeholder={placeholder} className="pr-10 font-mono" autoComplete="off" />
      <button type="button" onClick={() => setShow((s) => !s)} className="absolute right-2 top-1/2 -translate-y-1/2 rounded p-1 text-slate-400 hover:text-slate-600" aria-label="显示/隐藏">
        {show ? <EyeOff className="h-4 w-4" /> : <Eye className="h-4 w-4" />}
      </button>
    </div>
  );
}

export default function SettingsPage() {
  const toast = useToast();
  const [theme, setTheme] = useTheme();
  const { data: loaded, loading } = useAsync(() => api.getSettings(), []);
  const { data: health, reload: reloadHealth } = useAsync(() => api.health({ silent: true }), []);
  const [form, setForm] = useState<Settings | null>(null);
  const [saving, setSaving] = useState(false);
  const [testing, setTesting] = useState(false);
  const [test, setTest] = useState<SettingsTestResult | null>(null);

  useEffect(() => {
    if (loaded) setForm(loaded);
  }, [loaded]);

  const dirty = useMemo(() => JSON.stringify(form) !== JSON.stringify(loaded), [form, loaded]);

  if (loading || !form) {
    return (
      <Page>
        <PageHeader title="设置" />
        <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
          {Array.from({ length: 4 }).map((_, i) => (
            <Skeleton key={i} className="h-64" />
          ))}
        </div>
      </Page>
    );
  }

  const set = <K extends keyof Settings>(key: K, patch: Partial<Settings[K]>) => setForm((f) => (f ? { ...f, [key]: { ...f[key], ...patch } } : f));

  const runTest = async () => {
    setTesting(true);
    setTest(null);
    try {
      const r = await api.testSettings({ baseUrl: form.api.baseUrl, apiKey: form.api.apiKey.startsWith('sk-****') ? undefined : form.api.apiKey });
      setTest(r);
      if (r.ok) toast.success('连接成功', `${r.latencyMs} ms · ${r.models.length} 个模型`);
      else toast.error('连接失败', r.error);
    } finally {
      setTesting(false);
    }
  };

  const save = async () => {
    setSaving(true);
    try {
      const saved = await api.putSettings({ ...form, ui: { ...form.ui, theme } });
      setForm(saved);
      toast.success('设置已保存');
      void reloadHealth();
    } finally {
      setSaving(false);
    }
  };

  const modelOptions = Array.from(new Set([...(test?.models ?? []), ...Object.values(form.api.models).filter(Boolean)]));

  return (
    <Page>
      <PageHeader
        title="设置"
        description="API 网关、模型映射、搜索提供商、Python 环境与运行参数。密钥只保存在本机 ~/.bizgold/settings.json。"
        actions={
          <>
            {dirty ? <Badge tone="amber" dot>有未保存的修改</Badge> : null}
            <Button icon={<Save />} onClick={() => void save()} loading={saving} disabled={!dirty && !saving}>
              保存设置
            </Button>
          </>
        }
      />

      <div className="grid grid-cols-1 gap-6 xl:grid-cols-2">
        <Card>
          <CardHeader icon={<Globe />} title="API 网关" subtitle="Anthropic 兼容格式根地址（不带 /v1）" actions={<Button variant="secondary" size="sm" icon={<PlugZap />} onClick={() => void runTest()} loading={testing}>测试连接</Button>} />
          <CardBody className="space-y-4">
            <Field label="API 地址" required>
              <Input value={form.api.baseUrl} onChange={(e) => set('api', { baseUrl: e.target.value })} placeholder="https://llm.goaichat.top" className="font-mono" />
            </Field>
            <Field
              label="API 密钥"
              hint={
                <span className="inline-flex items-center gap-1.5">
                  <KeyRound className="h-3 w-3" />
                  来源：{form.api.apiKeySource === 'env' ? '环境变量 BIZGOLD_API_KEY / ANTHROPIC_AUTH_TOKEN（优先级高于设置）' : form.api.apiKeySource === 'settings' ? '设置文件（已脱敏显示）' : '未配置'}；留空或保持脱敏值则不修改。
                </span>
              }
            >
              <SecretInput value={form.api.apiKey} onChange={(v) => set('api', { apiKey: v })} placeholder="sk-…" />
            </Field>
            <Field label="备用地址" hint="主地址不可达时自动切换。">
              <Input value={form.api.fallbackBaseUrl} onChange={(e) => set('api', { fallbackBaseUrl: e.target.value })} placeholder="https://model.kimitk.top" className="font-mono" />
            </Field>
            <div className="grid grid-cols-2 gap-4">
              <Field label="上下文窗口（tokens）">
                <Input type="number" value={form.api.contextTokens} onChange={(e) => set('api', { contextTokens: Number(e.target.value) })} />
              </Field>
              <Field label="思考模式">
                <div className="pt-2">
                  <Switch checked={!form.api.disableThinking} onChange={(v) => set('api', { disableThinking: !v })} label={form.api.disableThinking ? '已关闭' : '已开启'} />
                </div>
              </Field>
            </div>
            {test ? (
              <div className={cn('rounded-lg border p-3 text-sm', test.ok ? 'border-emerald-200 bg-emerald-50/60 dark:border-emerald-900 dark:bg-emerald-900/10' : 'border-red-200 bg-red-50/60 dark:border-red-900 dark:bg-red-900/10')}>
                <div className="flex items-center gap-2 font-medium">
                  {test.ok ? <CheckCircle2 className="h-4 w-4 text-emerald-600" /> : <XCircle className="h-4 w-4 text-red-600" />}
                  {test.ok ? `连接成功 · 延迟 ${test.latencyMs} ms` : '连接失败'}
                  {test.baseUrlUsed ? <span className="ml-auto font-mono text-xs text-steel-500">{test.baseUrlUsed}</span> : null}
                </div>
                {test.ok ? (
                  <div className="mt-2 flex flex-wrap gap-1.5">
                    {test.models.map((m) => (
                      <span key={m} className="chip font-mono">
                        {m}
                      </span>
                    ))}
                  </div>
                ) : (
                  <p className="mt-1 text-xs text-red-700 dark:text-red-300">{test.error}</p>
                )}
              </div>
            ) : null}
          </CardBody>
        </Card>

        <Card>
          <CardHeader icon={<Cpu />} title="角色模型映射" subtitle={modelOptions.length && test?.ok ? `可选模型来自网关（${test.models.length} 个），点击「测试连接」刷新` : '先「测试连接」拉取网关模型列表，也可直接输入'} />
          <CardBody>
            <div className="overflow-hidden rounded-lg border border-slate-200 dark:border-surface-dark-border">
              <table className="table-base">
                <thead>
                  <tr>
                    <th>角色</th>
                    <th>模型</th>
                  </tr>
                </thead>
                <tbody>
                  {ROLES.map(({ role, label, desc }) => (
                    <tr key={role}>
                      <td className="w-[42%]">
                        <div className="font-medium text-slate-800 dark:text-slate-100">{label}</div>
                        <div className="text-xs text-steel-500">{desc}</div>
                      </td>
                      <td>
                        <Input list="model-list" value={form.api.models[role]} onChange={(e) => set('api', { models: { ...form.api.models, [role]: e.target.value } })} className="font-mono text-xs" />
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
              <datalist id="model-list">
                {modelOptions.map((m) => (
                  <option key={m} value={m} />
                ))}
              </datalist>
            </div>
          </CardBody>
        </Card>

        <Card>
          <CardHeader icon={<Search />} title="搜索提供商" subtitle="按顺序尝试；没有 key 的 API 自动跳过" />
          <CardBody className="space-y-3">
            {PROVIDERS.map((p) => {
              const enabled = form.search.providers.includes(p.id);
              const order = form.search.providers.indexOf(p.id);
              return (
                <div key={p.id} className="flex flex-col gap-2 rounded-lg border border-slate-200 p-3 dark:border-surface-dark-border md:flex-row md:items-center">
                  <div className="flex min-w-0 flex-1 items-center gap-3">
                    <Switch checked={enabled} onChange={(v) => set('search', { providers: v ? [...form.search.providers, p.id] : form.search.providers.filter((x) => x !== p.id) })} />
                    <div className="min-w-0">
                      <div className="flex items-center gap-2 text-sm font-medium text-slate-800 dark:text-slate-100">
                        {p.label}
                        {enabled ? <Badge tone="brand">#{order + 1}</Badge> : null}
                      </div>
                      <div className="text-xs text-steel-500">{p.note}</div>
                    </div>
                  </div>
                  {p.keyField ? (
                    <div className="md:w-72">
                      <SecretInput value={form.search[p.keyField]} onChange={(v) => set('search', { [p.keyField as string]: v } as Partial<Settings['search']>)} placeholder="API key" />
                    </div>
                  ) : null}
                </div>
              );
            })}
            <Field label="每次检索结果数">
              <Input type="number" min={3} max={20} value={form.search.maxResults} onChange={(e) => set('search', { maxResults: Number(e.target.value) })} className="w-32" />
            </Field>
          </CardBody>
        </Card>

        <Card>
          <CardHeader icon={<Wifi />} title="Python 与环境" subtitle="工具包渲染 Word / PPT / Excel / 图表依赖本机 Python ≥ 3.10" />
          <CardBody className="space-y-4">
            <Field label="Python 路径" hint="留空则按 python → python3 自动探测。">
              <Input value={form.python.path} onChange={(e) => set('python', { path: e.target.value })} placeholder={health?.python.path || 'python'} className="font-mono" />
            </Field>
            <Field label="pip 镜像">
              <Input value={form.python.pipIndex} onChange={(e) => set('python', { pipIndex: e.target.value })} className="font-mono" />
            </Field>
            <div className="rounded-lg border border-slate-200 p-3 text-sm dark:border-surface-dark-border">
              <div className="mb-2 text-xs font-medium uppercase tracking-wide text-steel-500">环境状态（来自 /health）</div>
              {health ? (
                <ul className="space-y-1.5">
                  <li className="flex items-center justify-between gap-2">
                    <span>Python</span>
                    <span className="flex items-center gap-2">
                      {health.python.found ? <Badge tone="green">{health.python.version}</Badge> : <Badge tone="red">未找到</Badge>}
                      <span className="max-w-[260px] truncate font-mono text-xs text-steel-500">{health.python.path}</span>
                    </span>
                  </li>
                  <li className="flex items-center justify-between gap-2">
                    <span>依赖包</span>
                    {health.python.missing.length ? <Badge tone="amber">缺少 {health.python.missing.join(', ')}</Badge> : <Badge tone="green">齐全</Badge>}
                  </li>
                  <li className="flex items-center justify-between gap-2">
                    <span>PDF 转换器</span>
                    {health.pdfConverter ? <Badge tone="green">{health.pdfConverter}</Badge> : <Badge tone="amber">未安装（Word / WPS / LibreOffice）</Badge>}
                  </li>
                  <li className="flex items-center justify-between gap-2">
                    <span>中文字体</span>
                    <span className="truncate text-xs text-steel-500">{health.fonts.cjk.join('、') || '—'}</span>
                  </li>
                  <li className="flex items-center justify-between gap-2">
                    <span>Git Bash（Windows）</span>
                    {health.gitBash ? <Badge tone="green">已检测</Badge> : <Badge tone="neutral">—</Badge>}
                  </li>
                </ul>
              ) : (
                <p className="text-xs text-steel-500">无法获取环境信息（服务端未连接）。</p>
              )}
            </div>
          </CardBody>
        </Card>

        <Card>
          <CardHeader icon={<SlidersHorizontal />} title="运行参数" subtitle="超时与保护、默认评审参数" />
          <CardBody>
            <div className="grid grid-cols-2 gap-4 md:grid-cols-3">
              <Field label="最大轮次 maxTurns">
                <Input type="number" value={form.run.maxTurns} onChange={(e) => set('run', { maxTurns: Number(e.target.value) })} />
              </Field>
              <Field label="总时长上限（分钟）">
                <Input type="number" value={form.run.totalTimeoutMin} onChange={(e) => set('run', { totalTimeoutMin: Number(e.target.value) })} />
              </Field>
              <Field label="无输出超时（分钟）">
                <Input type="number" value={form.run.idleTimeoutMin} onChange={(e) => set('run', { idleTimeoutMin: Number(e.target.value) })} />
              </Field>
              <Field label="Bash 超时（毫秒）">
                <Input type="number" value={form.run.bashTimeoutMs} onChange={(e) => set('run', { bashTimeoutMs: Number(e.target.value) })} />
              </Field>
              <Field label="默认评审轮数">
                <Input type="number" min={1} max={5} value={form.run.reviewRounds} onChange={(e) => set('run', { reviewRounds: Number(e.target.value) })} />
              </Field>
              <Field label="默认目标分">
                <Input type="number" min={60} max={100} value={form.run.targetScore} onChange={(e) => set('run', { targetScore: Number(e.target.value) })} />
              </Field>
            </div>
          </CardBody>
        </Card>

        <Card>
          <CardHeader icon={<Palette />} title="工作区与外观" />
          <CardBody className="space-y-4">
            <Field label="工作区目录" hint={`留空 = <repo>/workspace${health?.workspaceDir ? `；当前：${health.workspaceDir}` : ''}`}>
              <div className="flex gap-2">
                <Input value={form.paths.workspaceDir} onChange={(e) => set('paths', { workspaceDir: e.target.value })} placeholder={health?.workspaceDir || 'D:\\bizgold\\workspace'} className="font-mono" />
                <Button variant="secondary" icon={<FolderOpen />} onClick={() => toast.info('请在系统文件管理器中打开工作区目录', form.paths.workspaceDir || health?.workspaceDir)}>
                  打开
                </Button>
              </div>
            </Field>
            <div className="grid grid-cols-2 gap-4">
              <Field label="主题">
                <Select
                  value={theme}
                  onChange={(e) => {
                    const t = e.target.value as 'light' | 'dark';
                    setTheme(t);
                    set('ui', { theme: t });
                  }}
                >
                  <option value="light">浅色</option>
                  <option value="dark">深色</option>
                </Select>
              </Field>
              <Field label="语言">
                <Select value={form.ui.language} onChange={(e) => set('ui', { language: e.target.value })}>
                  <option value="zh-CN">简体中文</option>
                </Select>
              </Field>
            </div>
          </CardBody>
        </Card>
      </div>
    </Page>
  );
}
