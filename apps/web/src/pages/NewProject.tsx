import { useCallback, useEffect, useMemo, useState } from 'react';
import { useNavigate } from 'react-router-dom';
import { ArrowLeft, ArrowRight, CheckCircle2, Circle, Loader2, Rocket, RotateCcw, Save, XCircle } from 'lucide-react';
import { api } from '../api';
import { useAsync } from '../hooks/useAsync';
import { Page, PageHeader } from '../components/layout/AppShell';
import { Button } from '../components/ui/Button';
import { KeyValue } from '../components/ui/Misc';
import { Modal } from '../components/ui/Modal';
import { useToast } from '../components/ui/Toast';
import { Stepper } from '../components/wizard/Stepper';
import { StepTrack } from '../components/wizard/StepTrack';
import { StepInfo, type PendingFile } from '../components/wizard/StepInfo';
import { StepParams, ROLE_LABELS } from '../components/wizard/StepParams';
import { clearDraft, DRAFT_STEP_KEY, loadDraft, saveDraft, toCreateInput, validateStep1, validateStep2, validateStep3, type Draft, type Errors } from '../components/wizard/draft';
import { TRACK_NAMES } from '../components/domain/Brand';
import { loadJson, saveJson } from '../lib/storage';
import type { ProjectInput } from '../types';
import { cn } from '../lib/cn';

const STEPS = [
  { title: '赛道与组别', description: '选择赛道，查看评分权重' },
  { title: '项目信息', description: '作品、团队、资料' },
  { title: '生成参数', description: '模型与运行参数' },
];

type LaunchStep = { key: string; label: string; status: 'pending' | 'running' | 'done' | 'failed'; detail?: string };

export default function NewProjectPage() {
  const navigate = useNavigate();
  const toast = useToast();
  const { data: tracks } = useAsync(() => api.getTracks(), []);
  const { data: settings } = useAsync(() => api.getSettings(), []);
  const [draft, setDraft] = useState<Draft>(loadDraft);
  const [step, setStep] = useState<number>(() => Math.min(2, Math.max(0, loadJson<number>(DRAFT_STEP_KEY, 0))));
  const [maxReached, setMaxReached] = useState(step);
  const [errors, setErrors] = useState<Errors>({});
  const [touched, setTouched] = useState(false);
  const [files, setFiles] = useState<PendingFile[]>([]);
  const [launch, setLaunch] = useState<{ open: boolean; steps: LaunchStep[]; projectId?: string; uploaded?: ProjectInput[]; done: boolean; failed: boolean }>({ open: false, steps: [], done: false, failed: false });
  const [savedAt, setSavedAt] = useState<number | null>(null);

  const track = useMemo(() => tracks?.find((t) => t.id === draft.track) ?? null, [tracks, draft.track]);

  useEffect(() => {
    const t = window.setTimeout(() => {
      saveDraft(draft);
      setSavedAt(Date.now());
    }, 400);
    return () => window.clearTimeout(t);
  }, [draft]);

  useEffect(() => {
    saveJson(DRAFT_STEP_KEY, step);
    setMaxReached((m) => Math.max(m, step));
  }, [step]);

  // 用设置中的默认模型填充未指定的角色
  useEffect(() => {
    if (!settings) return;
    setDraft((d) => {
      const models = { ...d.params.models };
      let changed = false;
      for (const { role } of ROLE_LABELS) {
        if (!models[role] && settings.api.models[role]) {
          models[role] = settings.api.models[role];
          changed = true;
        }
      }
      return changed ? { ...d, params: { ...d.params, models } } : d;
    });
  }, [settings]);

  const update = useCallback((patch: Partial<Draft>) => {
    setDraft((d) => ({ ...d, ...patch }));
  }, []);

  const validate = useCallback(
    (s: number): Errors => {
      if (s === 0) return validateStep1(draft);
      if (s === 1) return validateStep2(draft);
      return validateStep3(draft, track?.pageLimit ?? 80);
    },
    [draft, track],
  );

  useEffect(() => {
    if (touched) setErrors(validate(step));
  }, [draft, step, touched, validate]);

  const next = () => {
    const e = validate(step);
    setErrors(e);
    setTouched(true);
    if (Object.keys(e).length) {
      toast.warning('请先完善必填项', Object.values(e)[0]);
      return;
    }
    setTouched(false);
    setStep((s) => Math.min(2, s + 1));
    window.scrollTo({ top: 0, behavior: 'smooth' });
  };

  const back = () => {
    setTouched(false);
    setStep((s) => Math.max(0, s - 1));
  };

  const reset = () => {
    clearDraft();
    setDraft(loadDraft());
    setFiles([]);
    setStep(0);
    setMaxReached(0);
    setTouched(false);
    setErrors({});
    toast.info('已清空草稿');
  };

  const addFiles = (fs: File[]) => {
    setFiles((prev) => {
      const names = new Set(prev.map((p) => p.file.name));
      const add = fs.filter((f) => !names.has(f.name)).map((f) => ({ id: `${f.name}-${f.size}-${Math.random().toString(36).slice(2, 7)}`, file: f }));
      return [...prev, ...add];
    });
  };

  const setLaunchStep = (key: string, patch: Partial<LaunchStep>) =>
    setLaunch((l) => ({ ...l, steps: l.steps.map((s) => (s.key === key ? { ...s, ...patch } : s)) }));

  const submit = async (run: boolean) => {
    const all = { ...validateStep1(draft), ...validateStep2(draft), ...validateStep3(draft, track?.pageLimit ?? 80) };
    if (Object.keys(all).length) {
      setErrors(all);
      setTouched(true);
      const firstStep = validateStep1(draft).track || validateStep1(draft).group ? 0 : Object.keys(validateStep2(draft)).length ? 1 : 2;
      setStep(firstStep);
      toast.warning('请先完善必填项', Object.values(all)[0]);
      return;
    }
    const steps: LaunchStep[] = [
      { key: 'create', label: '创建项目工作区', status: 'running' },
      ...(files.length ? [{ key: 'upload', label: `上传并预读资料（${files.length} 个）`, status: 'pending' as const }] : []),
      ...(run ? [{ key: 'run', label: '启动智能体运行', status: 'pending' as const }] : []),
    ];
    setLaunch({ open: true, steps, done: false, failed: false });
    try {
      const project = await api.createProject(toCreateInput(draft));
      setLaunchStep('create', { status: 'done', detail: project.id });
      setLaunch((l) => ({ ...l, projectId: project.id }));
      if (files.length) {
        setLaunchStep('upload', { status: 'running' });
        const inputs = await api.uploadFiles(
          project.id,
          files.map((f) => f.file),
        );
        setLaunch((l) => ({ ...l, uploaded: inputs }));
        setLaunchStep('upload', { status: 'done', detail: `${inputs.filter((i) => i.prereadStatus === 'done').length}/${inputs.length} 已预读` });
      }
      if (run) {
        setLaunchStep('run', { status: 'running' });
        await api.run(project.id, { mode: 'run' });
        setLaunchStep('run', { status: 'done' });
      }
      clearDraft();
      setLaunch((l) => ({ ...l, done: true }));
      toast.success(run ? '已开始生成' : '项目已创建', project.name);
      window.setTimeout(() => navigate(`/projects/${project.id}`), 700);
    } catch (e) {
      setLaunch((l) => ({ ...l, failed: true, steps: l.steps.map((s) => (s.status === 'running' ? { ...s, status: 'failed', detail: (e as Error).message } : s)) }));
    }
  };

  const summary = (
    <KeyValue
      cols={3}
      items={[
        { k: '赛道 / 组别', v: `${draft.track ? TRACK_NAMES[draft.track] : '—'} · ${draft.group || '—'}` },
        { k: '作品名', v: draft.name || '—' },
        { k: '团队名', v: draft.teamName || '—' },
        { k: '策划对象', v: `${draft.subject.name || '—'}${draft.subject.company ? `（${draft.subject.company}）` : ''}` },
        { k: '命题', v: draft.subject.isAssignedCase ? draft.subject.caseTitle || '—' : '自选' },
        { k: '团队成员', v: `${draft.team.length} 人（匿名）` },
        { k: '调研模式', v: draft.params.researchMode === 'synthetic' ? '示例数据（标注）' : '仅用我的资料' },
        { k: '页数 / 讲稿', v: `${draft.params.pageTarget} 页 · ${draft.params.pitchMinutes} 分钟` },
        { k: '评审 / 目标分', v: `${draft.params.reviewRounds} 轮 · ${draft.params.targetScore} 分` },
        { k: '主代理模型', v: draft.params.models.main || '设置默认' },
        { k: '评委团', v: [draft.params.models.judgeA, draft.params.models.judgeB, draft.params.models.judgeC].filter(Boolean).join(' / ') || '设置默认' },
        { k: '资料', v: files.length ? `${files.length} 个文件待上传` : '无' },
      ]}
    />
  );

  return (
    <Page>
      <PageHeader
        title="新建项目"
        description="三步完成：赛道与组别 → 项目信息 → 生成参数。草稿自动保存到本地。"
        actions={
          <>
            {savedAt ? (
              <span className="inline-flex items-center gap-1 text-xs text-steel-500">
                <Save className="h-3.5 w-3.5" />
                草稿已保存
              </span>
            ) : null}
            <Button variant="ghost" size="sm" icon={<RotateCcw />} onClick={reset}>
              清空草稿
            </Button>
          </>
        }
      />

      <div className="card mb-6 px-5 py-4">
        <Stepper steps={STEPS} current={step} maxReachable={maxReached} onSelect={(i) => setStep(i)} />
      </div>

      {step === 0 ? <StepTrack tracks={tracks} draft={draft} onChange={update} errors={errors} /> : null}
      {step === 1 ? <StepInfo track={track} draft={draft} onChange={update} errors={errors} files={files} onAddFiles={addFiles} onRemoveFile={(id) => setFiles((fs) => fs.filter((f) => f.id !== id))} /> : null}
      {step === 2 ? (
        <div className="space-y-6">
          <StepParams track={track} draft={draft} onChange={update} errors={errors} settings={settings} />
          <section className="card overflow-hidden">
            <div className="border-b border-slate-100 bg-gradient-to-r from-brand-50 to-gold-50/60 px-5 py-4 dark:border-surface-dark-border dark:from-brand-900/40 dark:to-gold-900/10">
              <h2 className="text-base font-semibold text-brand-700 dark:text-slate-100">摘要确认</h2>
              <p className="mt-0.5 text-xs text-steel-500 dark:text-steel-400">确认无误后开始生成；智能体将按 9 个阶段自动完成全部材料并迭代到目标分。</p>
            </div>
            <div className="px-5 py-4">{summary}</div>
          </section>
        </div>
      ) : null}

      <div className="sticky bottom-0 z-20 mt-6 -mx-6 border-t border-slate-200/80 bg-white/85 px-6 py-3 backdrop-blur dark:border-surface-dark-border dark:bg-surface-dark/85">
        <div className="mx-auto flex max-w-[1440px] items-center justify-between gap-3">
          <Button variant="secondary" icon={<ArrowLeft />} onClick={back} disabled={step === 0}>
            上一步
          </Button>
          <div className="flex items-center gap-2">
            {step < 2 ? (
              <Button onClick={next} iconRight={<ArrowRight />}>
                下一步
              </Button>
            ) : (
              <>
                <Button variant="secondary" icon={<Save />} onClick={() => void submit(false)}>
                  仅保存
                </Button>
                <Button variant="gold" size="lg" icon={<Rocket />} onClick={() => void submit(true)}>
                  开始生成
                </Button>
              </>
            )}
          </div>
        </div>
      </div>

      <Modal open={launch.open} onClose={() => (launch.done || launch.failed ? setLaunch((l) => ({ ...l, open: false })) : undefined)} title={launch.failed ? '启动失败' : launch.done ? '已就绪' : '正在创建项目…'} size="sm"
        footer={
          launch.failed ? (
            <>
              <Button variant="secondary" onClick={() => setLaunch((l) => ({ ...l, open: false }))}>
                关闭
              </Button>
              {launch.projectId ? <Button onClick={() => navigate(`/projects/${launch.projectId}`)}>进入项目</Button> : null}
            </>
          ) : launch.done && launch.projectId ? (
            <Button onClick={() => navigate(`/projects/${launch.projectId}`)}>进入项目</Button>
          ) : null
        }
      >
        <ul className="space-y-3">
          {launch.steps.map((s) => (
            <li key={s.key} className="flex items-start gap-3 text-sm">
              {s.status === 'done' ? <CheckCircle2 className="mt-0.5 h-4 w-4 text-emerald-500" /> : s.status === 'running' ? <Loader2 className="mt-0.5 h-4 w-4 animate-spin text-brand-500" /> : s.status === 'failed' ? <XCircle className="mt-0.5 h-4 w-4 text-red-500" /> : <Circle className="mt-0.5 h-4 w-4 text-slate-300" />}
              <div>
                <div className={cn('font-medium', s.status === 'pending' ? 'text-slate-400' : 'text-slate-800 dark:text-slate-100')}>{s.label}</div>
                {s.detail ? <div className="mt-0.5 break-all text-xs text-steel-500">{s.detail}</div> : null}
              </div>
            </li>
          ))}
        </ul>
        {launch.uploaded?.length ? (
          <ul className="mt-4 space-y-1 rounded-lg bg-slate-50 p-3 text-xs dark:bg-slate-800/60">
            {launch.uploaded.map((u) => (
              <li key={u.name} className="flex items-center justify-between gap-2">
                <span className="truncate text-slate-700 dark:text-slate-200">{u.name}</span>
                <span className="shrink-0 text-emerald-600 dark:text-emerald-300">{u.prereadStatus === 'done' ? (u.prereadChars ? `已预读 ${(u.prereadChars / 10000).toFixed(1)} 万字` : '已上传') : u.prereadStatus === 'failed' ? '预读失败' : '预读中'}</span>
              </li>
            ))}
          </ul>
        ) : null}
      </Modal>
    </Page>
  );
}
