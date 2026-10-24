import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { Link, useParams, useSearchParams } from 'react-router-dom';
import { Activity, ArrowLeft, Building2, Database, FolderOpen, FolderTree, Gavel, MapPin, MessageSquareText, Package, Play, RotateCcw, Square, Target, Users } from 'lucide-react';
import type { AgentEvent, Project, StageId } from '../types';
import { api, IS_MOCK } from '../api';
import { useAsync, useNow } from '../hooks/useAsync';
import { useConsole } from '../hooks/useConsole';
import { STAGE_IDS, STAGES, stageName, stageProgressPercent } from '../lib/stages';
import { formatDuration, formatRelative } from '../lib/format';
import { Page } from '../components/layout/AppShell';
import { Button } from '../components/ui/Button';
import { Badge } from '../components/ui/Badge';
import { EmptyState, Skeleton, Tabs } from '../components/ui/Misc';
import { ConfirmDialog, Modal } from '../components/ui/Modal';
import { Select, Field } from '../components/ui/Form';
import { useToast } from '../components/ui/Toast';
import { ScoreText, StatusPill, TrackBadge, TRACK_NAMES } from '../components/domain/Brand';
import { ConsoleTab, SseIndicator } from '../components/project/ConsoleTab';
import { ArtifactsTab } from '../components/project/ArtifactsTab';
import { ReviewTab } from '../components/project/ReviewTab';
import { ChatTab } from '../components/project/ChatTab';
import { NumbersTab } from '../components/project/NumbersTab';

type Tab = 'console' | 'artifacts' | 'review' | 'chat' | 'numbers';
const TABS: Tab[] = ['console', 'artifacts', 'review', 'chat', 'numbers'];

const SUBJECT_TYPE: Record<Project['subject']['type'], string> = { brand: '品牌', startup: '创业项目', 'case-company': '案例企业', event: '会奖活动' };

export default function ProjectPage() {
  const { id = '' } = useParams();
  const toast = useToast();
  const [search, setSearch] = useSearchParams();
  const tabParam = search.get('tab');
  const tab: Tab = TABS.includes(tabParam as Tab) ? (tabParam as Tab) : 'console';
  const setTab = useCallback(
    (t: Tab) => {
      const next = new URLSearchParams(search);
      if (t === 'console') next.delete('tab');
      else next.set('tab', t);
      setSearch(next, { replace: true });
    },
    [search, setSearch],
  );

  const fetchedAt = useRef(Date.now());
  const { data: project, error, reload } = useAsync(
    async () => {
      const p = await api.getProject(id, { silent: true });
      fetchedAt.current = Date.now();
      return p;
    },
    [id],
    { pollMs: 4000 },
  );

  const [artifactsKey, setArtifactsKey] = useState(0);
  const [reviewKey, setReviewKey] = useState(0);
  const [numbersKey, setNumbersKey] = useState(0);
  const [selectedPath, setSelectedPath] = useState<string | null>(null);
  const timers = useRef<{ project?: number; artifacts?: number }>({});

  const onEvent = useCallback(
    (ev: AgentEvent) => {
      if (ev.type === 'stage' || ev.type === 'status' || ev.type === 'review' || ev.type === 'done' || ev.type === 'usage') {
        if (timers.current.project == null) {
          timers.current.project = window.setTimeout(() => {
            timers.current.project = undefined;
            void reload();
          }, 250);
        }
        if (ev.type === 'review' || ev.type === 'done') setReviewKey((k) => k + 1);
      }
      if (ev.type === 'file') {
        if (timers.current.artifacts == null) {
          timers.current.artifacts = window.setTimeout(() => {
            timers.current.artifacts = undefined;
            setArtifactsKey((k) => k + 1);
          }, 800);
        }
        if (ev.data.path.endsWith('numbers.json')) setNumbersKey((k) => k + 1);
      }
    },
    [reload],
  );
  const { state, sse, retryInMs } = useConsole(id || undefined, onEvent);
  useEffect(() => {
    const t = timers.current;
    return () => {
      if (t.project) window.clearTimeout(t.project);
      if (t.artifacts) window.clearTimeout(t.artifacts);
    };
  }, []);

  const running = project?.status === 'running';
  const now = useNow(1000, Boolean(running));

  const [busy, setBusy] = useState<string | null>(null);
  const [confirmRerun, setConfirmRerun] = useState(false);
  const [resumeOpen, setResumeOpen] = useState(false);
  const [resumeFrom, setResumeFrom] = useState<StageId>('parse');

  const firstPending = useMemo(() => (project ? STAGE_IDS.find((s) => project.stages[s].status !== 'done' && project.stages[s].status !== 'skipped') ?? 'package' : 'parse'), [project]);

  const act = async (name: string, fn: () => Promise<unknown>, success?: string) => {
    setBusy(name);
    try {
      await fn();
      if (success) toast.success(success);
      await reload();
    } catch {
      // api 已统一 toast
    } finally {
      setBusy(null);
    }
  };

  const start = () => {
    if (!project) return;
    if (project.status !== 'idle') {
      setConfirmRerun(true);
      return;
    }
    void act('run', () => api.run(project.id, { mode: 'run' }), '已开始运行');
  };
  const confirmStart = () => {
    if (!project) return;
    setConfirmRerun(false);
    void act('run', () => api.run(project.id, { mode: 'run' }), '已重新开始运行');
    setTab('console');
  };
  const stop = () => project && void act('stop', () => api.stop(project.id), '已发送停止指令');
  const openResume = () => {
    setResumeFrom(firstPending);
    setResumeOpen(true);
  };
  const resume = () => {
    if (!project) return;
    setResumeOpen(false);
    void act('resume', () => api.run(project.id, { mode: 'resume', fromStage: resumeFrom }), `从「${stageName(resumeFrom)}」续跑`);
    setTab('console');
  };
  const repackage = () =>
    project &&
    void act(
      'package',
      async () => {
        const r = await api.repackage(project.id);
        setArtifactsKey((k) => k + 1);
        setSelectedPath(r.path);
        setTab('artifacts');
      },
      '提交包已重新生成',
    );
  const openFolder = () => project && void act('folder', () => api.openFolder(project.id), IS_MOCK ? 'Mock 模式：已模拟打开工作区文件夹' : '已在资源管理器中打开工作区');

  const openFile = useCallback(
    (path: string) => {
      setSelectedPath(path);
      setTab('artifacts');
    },
    [setTab],
  );

  if (error && !project) {
    return (
      <Page>
        <EmptyState
          icon={<FolderTree />}
          title="项目不存在或已被删除"
          description={error.message}
          action={
            <Link to="/">
              <Button variant="secondary" icon={<ArrowLeft />}>
                返回工作台
              </Button>
            </Link>
          }
        />
      </Page>
    );
  }

  if (!project) {
    return (
      <Page wide>
        <Skeleton className="h-6 w-40" />
        <Skeleton className="mt-4 h-10 w-2/3" />
        <Skeleton className="mt-3 h-4 w-1/2" />
        <Skeleton className="mt-8 h-24 w-full" />
        <Skeleton className="mt-4 h-96 w-full" />
      </Page>
    );
  }

  const score = project.review.latestScore;
  const progress = stageProgressPercent(project.stages);
  const canStop = running;
  const canResume = !running && (project.status === 'paused' || project.status === 'failed' || project.status === 'stopped');
  const canPackage = !running && Object.keys(project.artifacts).length > 0;

  return (
    <Page wide className="pb-10">
      <div className="mb-1 flex items-center gap-2 text-xs text-steel-500">
        <Link to="/" className="hover:text-brand-600 hover:underline dark:hover:text-brand-200">
          工作台
        </Link>
        <span>/</span>
        <span className="truncate">{project.name}</span>
      </div>

      <div className="flex flex-wrap items-start justify-between gap-4">
        <div className="flex min-w-0 items-start gap-4">
          <TrackBadge track={project.track} size="lg" />
          <div className="min-w-0">
            <div className="flex flex-wrap items-center gap-2">
              <h1 className="text-balance text-2xl font-semibold tracking-tight text-brand-700 dark:text-slate-50">{project.name}</h1>
              <StatusPill status={project.status} />
              {running ? (
                <Badge tone="neutral">
                  <Activity className="h-3 w-3" />
                  {project.session?.mode === 'modify' ? '修改模式' : `${stageName(project.currentStage)} · ${progress}%`}
                </Badge>
              ) : null}
            </div>
            <div className="mt-1.5 flex flex-wrap items-center gap-x-3 gap-y-1 text-sm text-steel-500 dark:text-steel-400">
              <span>{TRACK_NAMES[project.track]}</span>
              <span>·</span>
              <span>{project.group}</span>
              <span>·</span>
              <span className="inline-flex items-center gap-1">
                <Users className="h-3.5 w-3.5" />
                {project.teamName}
              </span>
              <span>·</span>
              <span className="inline-flex items-center gap-1">
                <Building2 className="h-3.5 w-3.5" />
                {SUBJECT_TYPE[project.subject.type]} {project.subject.name}
                {project.subject.company && project.subject.company !== project.subject.name ? `（${project.subject.company}）` : ''}
              </span>
              {project.subject.region ? (
                <span className="inline-flex items-center gap-1">
                  <MapPin className="h-3.5 w-3.5" />
                  {project.subject.region}
                </span>
              ) : null}
              {project.subject.isAssignedCase ? <Badge tone="gold">命题企业</Badge> : null}
            </div>
            <div className="mt-2 flex flex-wrap gap-1.5 text-[11px]">
              <span className="chip">正文目标 {project.params.pageTarget} 页</span>
              <span className="chip">路演 {project.params.pitchMinutes} 分钟</span>
              <span className="chip">评审 ≤ {project.params.reviewRounds} 轮</span>
              <span className="chip">
                <Target className="h-3 w-3" />
                目标 {project.params.targetScore} 分
              </span>
              <span className="chip">{project.params.researchMode === 'synthetic' ? '示例数据模式' : '仅用我的资料'}</span>
              <span className="chip">{project.inputs.length} 份资料</span>
              <span className="chip">更新于 {formatRelative(project.updatedAt, now)}</span>
              <span className="chip">
                <SseIndicator status={sse} retryInMs={retryInMs} />
              </span>
            </div>
          </div>
        </div>

        <div className="flex flex-wrap items-center gap-2">
          <div className="mr-2 text-right">
            <div className="text-[11px] text-steel-500">最新得分 / 目标 {project.params.targetScore}</div>
            <div className="text-2xl leading-tight">
              <ScoreText score={score} target={project.params.targetScore} />
            </div>
          </div>
          {canStop ? (
            <Button variant="danger" icon={<Square />} onClick={stop} loading={busy === 'stop'}>
              停止
            </Button>
          ) : (
            <Button variant="gold" icon={<Play />} onClick={start} loading={busy === 'run'} disabled={busy != null}>
              {project.status === 'idle' ? '开始生成' : '重新开始'}
            </Button>
          )}
          <Button variant="secondary" icon={<RotateCcw />} onClick={openResume} disabled={!canResume || busy != null} loading={busy === 'resume'} title={canResume ? '从未完成的阶段继续' : '仅暂停 / 失败 / 已停止的项目可续跑'}>
            续跑
          </Button>
          <Button variant="secondary" icon={<Package />} onClick={repackage} disabled={!canPackage || busy != null} loading={busy === 'package'} title={canPackage ? '重新生成 08_提交包' : '尚无产出物'}>
            重新打包
          </Button>
          <Button variant="secondary" icon={<FolderOpen />} onClick={openFolder} loading={busy === 'folder'} disabled={busy != null}>
            打开文件夹
          </Button>
        </div>
      </div>

      <Tabs
        className="mt-5"
        value={tab}
        onChange={setTab}
        tabs={[
          { value: 'console', label: '运行控制台', icon: <Activity />, badge: running ? <span className="h-1.5 w-1.5 animate-pulse-dot rounded-full bg-emerald-500" /> : null },
          { value: 'artifacts', label: '产出物', icon: <FolderTree />, badge: Object.keys(project.artifacts).length ? <span className="rounded-full bg-slate-100 px-1.5 text-[10px] tabular text-steel-600 dark:bg-slate-800 dark:text-steel-300">{Object.keys(project.artifacts).length}</span> : null },
          { value: 'review', label: '评审', icon: <Gavel />, badge: score != null ? <span className="rounded-full bg-gold-100 px-1.5 text-[10px] font-semibold tabular text-gold-800 dark:bg-gold-900/40 dark:text-gold-200">{score.toFixed(1)}</span> : null },
          { value: 'chat', label: '修改对话', icon: <MessageSquareText /> },
          { value: 'numbers', label: '数字账本', icon: <Database /> },
        ]}
      />

      <div className="mt-4">
        {tab === 'console' ? <ConsoleTab project={project} state={state} sse={sse} retryInMs={retryInMs} now={now} fetchedAt={fetchedAt.current} onOpenFile={openFile} /> : null}
        {tab === 'artifacts' ? <ArtifactsTab project={project} refreshKey={artifactsKey} selectedPath={selectedPath} onSelectPath={setSelectedPath} now={now} /> : null}
        {tab === 'review' ? <ReviewTab project={project} refreshKey={reviewKey} /> : null}
        {tab === 'chat' ? <ChatTab project={project} state={state} onProjectChanged={() => void reload()} /> : null}
        {tab === 'numbers' ? <NumbersTab project={project} refreshKey={numbersKey} /> : null}
      </div>

      <ConfirmDialog
        open={confirmRerun}
        onClose={() => setConfirmRerun(false)}
        onConfirm={confirmStart}
        title="重新开始运行"
        description={
          <>
            将从「读题解析」重新执行全部 9 个阶段，现有产出物、评审记录与用量将被重置（工作区文件会被覆盖）。如只想继续未完成的阶段，请使用「续跑」。
            {project.usage.durationMs ? <span className="mt-1 block text-xs text-steel-400">上次运行用时 {formatDuration(project.usage.durationMs)}。</span> : null}
          </>
        }
        confirmText="重新开始"
        danger
        loading={busy === 'run'}
      />

      <Modal
        open={resumeOpen}
        onClose={() => setResumeOpen(false)}
        title="续跑"
        description="沿用上次会话上下文，从所选阶段继续；之前阶段的产出保留。"
        size="sm"
        footer={
          <>
            <Button variant="secondary" onClick={() => setResumeOpen(false)}>
              取消
            </Button>
            <Button variant="primary" icon={<RotateCcw />} onClick={resume}>
              开始续跑
            </Button>
          </>
        }
      >
        <Field label="从哪个阶段开始" hint={`系统建议：${stageName(firstPending)}（第一个未完成的阶段）`}>
          <Select value={resumeFrom} onChange={(e) => setResumeFrom(e.target.value as StageId)}>
            {STAGES.map((s) => (
              <option key={s.id} value={s.id}>
                {s.index}. {s.name} · {project.stages[s.id].status === 'done' ? '已完成' : project.stages[s.id].status === 'failed' ? '失败' : '未完成'}
              </option>
            ))}
          </Select>
        </Field>
      </Modal>
    </Page>
  );
}
