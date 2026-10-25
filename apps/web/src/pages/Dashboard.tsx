import { useMemo, useState } from 'react';
import { Link, useNavigate } from 'react-router-dom';
import { Activity, Clock, FileStack, FolderKanban, Plus, Sparkles, Trash2, Trophy } from 'lucide-react';
import { api } from '../api';
import { useAsync, useNow } from '../hooks/useAsync';
import { Page, PageHeader } from '../components/layout/AppShell';
import { Button } from '../components/ui/Button';
import { EmptyState, ProgressRing, Skeleton, Stat } from '../components/ui/Misc';
import { ConfirmDialog } from '../components/ui/Modal';
import { useToast } from '../components/ui/Toast';
import { ScoreText, StatusPill, TrackBadge, TRACK_NAMES } from '../components/domain/Brand';
import { STAGES, stageName, stageProgressPercent } from '../lib/stages';
import { formatRelative } from '../lib/format';
import { cn } from '../lib/cn';
import type { Project } from '../types';

function StageMini({ project }: { project: Project }) {
  return (
    <div className="mt-3">
      <div className="mb-1.5 flex items-center justify-between text-[11px] text-steel-500">
        <span className="inline-flex items-center gap-1">
          <span className="h-1.5 w-1.5 animate-pulse-dot rounded-full bg-emerald-500" />
          当前阶段：{stageName(project.currentStage)}
        </span>
        <span className="tabular">{stageProgressPercent(project.stages)}%</span>
      </div>
      <div className="flex gap-0.5">
        {STAGES.map((s) => {
          const st = project.stages[s.id]?.status;
          return (
            <span
              key={s.id}
              title={`${s.name}：${st}`}
              className={cn(
                'h-1.5 flex-1 rounded-sm',
                st === 'done' ? 'bg-brand-500 dark:bg-brand-300' : st === 'running' ? 'animate-pulse bg-emerald-500' : st === 'failed' ? 'bg-red-500' : 'bg-slate-200 dark:bg-slate-700',
              )}
            />
          );
        })}
      </div>
    </div>
  );
}

function ProjectCard({ project, onDelete, now }: { project: Project; onDelete: (p: Project) => void; now: number }) {
  const navigate = useNavigate();
  const score = project.review.latestScore;
  const target = project.params.targetScore;
  return (
    <div
      className="card card-hover group relative flex cursor-pointer flex-col p-5"
      onClick={() => navigate(`/projects/${project.id}`)}
      role="link"
      tabIndex={0}
      onKeyDown={(e) => {
        if (e.key === 'Enter') navigate(`/projects/${project.id}`);
      }}
    >
      <div className="flex items-start justify-between gap-3">
        <div className="flex min-w-0 items-start gap-3">
          <TrackBadge track={project.track} />
          <div className="min-w-0">
            <h3 className="line-clamp-2 text-[15px] font-semibold leading-snug text-slate-800 dark:text-slate-100" title={project.name}>
              {project.name}
            </h3>
            <p className="mt-1 truncate text-xs text-steel-500 dark:text-steel-400">
              {TRACK_NAMES[project.track]} · {project.group} · {project.teamName}
            </p>
          </div>
        </div>
        <ProgressRing
          value={score}
          size={54}
          stroke={5}
          label={
            <div>
              <div className="text-[13px] leading-none">
                <ScoreText score={score} target={target} />
              </div>
              <div className="mt-0.5 text-[9px] text-steel-400">{score == null ? '未评' : `目标 ${target}`}</div>
            </div>
          }
        />
      </div>

      <div className="mt-4 flex items-center justify-between gap-2">
        <StatusPill status={project.status} />
        <span className="inline-flex items-center gap-1 text-xs text-steel-500 dark:text-steel-400">
          <Clock className="h-3.5 w-3.5" />
          {formatRelative(project.updatedAt, now)}
        </span>
      </div>

      {project.status === 'running' ? (
        <StageMini project={project} />
      ) : (
        <div className="mt-3 flex flex-wrap gap-1.5 text-[11px] text-steel-500">
          {project.status === 'completed' ? (
            <>
              <span className="chip">{project.review.rounds.length} 轮评审</span>
              <span className="chip">{Object.keys(project.artifacts).length} 项产出</span>
            </>
          ) : project.status === 'failed' || project.status === 'stopped' ? (
            <span className="chip text-red-600 dark:text-red-300">停在「{stageName(project.currentStage)}」· {project.stages[project.currentStage ?? 'parse']?.note ?? '可续跑'}</span>
          ) : (
            <span className="chip">尚未开始运行</span>
          )}
        </div>
      )}

      <button
        className="absolute right-3 top-3 hidden rounded-md p-1.5 text-slate-400 hover:bg-red-50 hover:text-red-600 group-hover:inline-flex dark:hover:bg-red-900/30"
        onClick={(e) => {
          e.stopPropagation();
          onDelete(project);
        }}
        title="删除项目"
        aria-label="删除项目"
      >
        <Trash2 className="h-4 w-4" />
      </button>
    </div>
  );
}

export default function DashboardPage() {
  const { data: projects, loading, reload, setData } = useAsync(() => api.getProjects(), [], { pollMs: 5000 });
  const [toDelete, setToDelete] = useState<Project | null>(null);
  const [deleting, setDeleting] = useState(false);
  const toast = useToast();
  const now = useNow(30_000);

  const stats = useMemo(() => {
    const list = projects ?? [];
    const running = list.filter((p) => p.status === 'running').length;
    const scored = list.filter((p) => p.review.latestScore != null);
    const avg = scored.length ? scored.reduce((s, p) => s + (p.review.latestScore ?? 0), 0) / scored.length : null;
    const monthStart = new Date();
    monthStart.setDate(1);
    monthStart.setHours(0, 0, 0, 0);
    const materials = list.filter((p) => Date.parse(p.updatedAt) >= monthStart.getTime()).reduce((s, p) => s + Object.keys(p.artifacts).length, 0);
    return { total: list.length, running, avg, materials, scored: scored.length };
  }, [projects]);

  const confirmDelete = async () => {
    if (!toDelete) return;
    setDeleting(true);
    try {
      await api.deleteProject(toDelete.id);
      setData((prev) => (prev ?? []).filter((p) => p.id !== toDelete.id));
      toast.success('项目已删除', toDelete.name);
      setToDelete(null);
    } catch {
      // toast 已由 api 统一处理
    } finally {
      setDeleting(false);
      void reload();
    }
  };

  return (
    <Page>
      <PageHeader
        title="工作台"
        description="管理参赛项目，查看运行进度与评审得分。"
        actions={
          <Link to="/new">
            <Button variant="gold" icon={<Plus />}>
              新建项目
            </Button>
          </Link>
        }
      />

      <div className="grid grid-cols-2 gap-4 lg:grid-cols-4">
        <Stat label="项目数" value={loading && !projects ? <Skeleton className="h-7 w-12" /> : stats.total} sub={`${stats.scored} 个已评审`} icon={<FolderKanban />} tone="brand" />
        <Stat label="运行中" value={loading && !projects ? <Skeleton className="h-7 w-12" /> : stats.running} sub={stats.running ? '实时流推送中' : '当前无运行任务'} icon={<Activity />} tone="green" />
        <Stat label="平均得分" value={loading && !projects ? <Skeleton className="h-7 w-16" /> : stats.avg == null ? '—' : stats.avg.toFixed(1)} sub="国金线 85 分" icon={<Trophy />} tone="gold" />
        <Stat label="本月生成材料" value={loading && !projects ? <Skeleton className="h-7 w-12" /> : stats.materials} sub="策划书 / 模型 / 路演 / 讲稿 / 题库 / 提交包" icon={<FileStack />} tone="steel" />
      </div>

      <div className="mt-8">
        <div className="mb-3 flex items-center justify-between">
          <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">项目</h2>
          <span className="text-xs text-steel-500">按更新时间排序</span>
        </div>
        {loading && !projects ? (
          <div className="grid grid-cols-1 gap-4 md:grid-cols-2 2xl:grid-cols-3">
            {Array.from({ length: 6 }).map((_, i) => (
              <div key={i} className="card p-5">
                <div className="flex items-start gap-3">
                  <Skeleton className="h-8 w-8" />
                  <div className="flex-1 space-y-2">
                    <Skeleton className="h-4 w-3/4" />
                    <Skeleton className="h-3 w-1/2" />
                  </div>
                  <Skeleton className="h-12 w-12 rounded-full" />
                </div>
                <Skeleton className="mt-5 h-4 w-20" />
              </div>
            ))}
          </div>
        ) : !projects?.length ? (
          <div className="card">
            <EmptyState
              icon={<Sparkles />}
              title="还没有项目"
              description="从新建项目开始：选择赛道与组别、填写项目信息与团队、上传企业资料，系统将自动完成读题、调研、财务、撰写、路演、评审与打包。"
              action={
                <Link to="/new">
                  <Button variant="gold" icon={<Plus />}>
                    新建第一个项目
                  </Button>
                </Link>
              }
            />
          </div>
        ) : (
          <div className="grid grid-cols-1 gap-4 md:grid-cols-2 2xl:grid-cols-3">
            {projects.map((p) => (
              <ProjectCard key={p.id} project={p} onDelete={setToDelete} now={now} />
            ))}
          </div>
        )}
      </div>

      <ConfirmDialog
        open={Boolean(toDelete)}
        onClose={() => setToDelete(null)}
        onConfirm={() => void confirmDelete()}
        title="删除项目"
        description={
          toDelete ? (
            <>
              将删除「{toDelete.name}」及其工作区内全部产出物（{toDelete.status === 'running' ? '当前运行将被中止，' : ''}此操作不可恢复）。
            </>
          ) : null
        }
        confirmText="确认删除"
        danger
        loading={deleting}
      />
    </Page>
  );
}
