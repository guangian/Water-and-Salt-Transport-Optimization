import { useEffect, useMemo, useState, type ReactNode } from 'react';
import {
  Archive,
  BookText,
  ChevronDown,
  ChevronRight,
  File,
  FileCode2,
  FileImage,
  FileSpreadsheet,
  FileText,
  FileType2,
  Folder,
  FolderOpen,
  HelpCircle,
  Loader2,
  Mic,
  Presentation,
  RefreshCw,
  Search,
} from 'lucide-react';
import type { ArtifactKey, ArtifactNode, Project } from '../../types';
import { api } from '../../api';
import { useAsync } from '../../hooks/useAsync';
import { baseName, fileExt, formatBytes, formatRelative } from '../../lib/format';
import { parentDir } from '../../lib/paths';
import { cn } from '../../lib/cn';
import { Card } from '../ui/Card';
import { Button } from '../ui/Button';
import { EmptyState } from '../ui/Misc';
import { Input } from '../ui/Form';
import { ArtifactPreview } from './ArtifactPreview';

const QUICK: { key: ArtifactKey; label: string; icon: ReactNode; hint: string }[] = [
  { key: 'book', label: '策划书 · Word', icon: <BookText />, hint: 'docx' },
  { key: 'pdf', label: '策划书 · PDF', icon: <FileType2 />, hint: 'pdf' },
  { key: 'deck', label: '路演 PPT', icon: <Presentation />, hint: 'pptx' },
  { key: 'model', label: '财务模型', icon: <FileSpreadsheet />, hint: 'xlsx' },
  { key: 'script', label: '路演讲稿', icon: <Mic />, hint: 'md' },
  { key: 'qa', label: '答辩题库', icon: <HelpCircle />, hint: 'md' },
  { key: 'package', label: '提交包', icon: <Archive />, hint: 'zip' },
];

function fileIcon(name: string): ReactNode {
  const ext = fileExt(name);
  if (['png', 'jpg', 'jpeg', 'gif', 'webp', 'svg'].includes(ext)) return <FileImage className="text-violet-500" />;
  if (['xlsx', 'xls', 'csv'].includes(ext)) return <FileSpreadsheet className="text-emerald-600" />;
  if (ext === 'pptx') return <Presentation className="text-orange-500" />;
  if (ext === 'docx' || ext === 'doc') return <BookText className="text-brand-500" />;
  if (ext === 'pdf') return <FileType2 className="text-red-500" />;
  if (ext === 'zip') return <Archive className="text-gold-600" />;
  if (ext === 'json') return <FileCode2 className="text-sky-600" />;
  if (ext === 'md' || ext === 'txt') return <FileText className="text-steel-500" />;
  return <File className="text-steel-400" />;
}

function flatten(node: ArtifactNode, out: ArtifactNode[] = []): ArtifactNode[] {
  if (node.type === 'file') out.push(node);
  node.children?.forEach((c) => flatten(c, out));
  return out;
}

function TreeNode({ node, depth, expanded, toggle, selected, onSelect, now }: { node: ArtifactNode; depth: number; expanded: Set<string>; toggle: (p: string) => void; selected: string | null; onSelect: (n: ArtifactNode) => void; now: number }) {
  const isDir = node.type === 'dir';
  const open = expanded.has(node.path);
  const fresh = node.type === 'file' && node.mtime && now - Date.parse(node.mtime) < 2 * 60_000;
  return (
    <div>
      <button
        onClick={() => (isDir ? toggle(node.path) : onSelect(node))}
        className={cn(
          'flex w-full items-center gap-1.5 rounded-md py-1 pr-2 text-left text-[13px] transition',
          selected === node.path ? 'bg-brand-50 text-brand-700 dark:bg-brand-900/40 dark:text-brand-100' : 'text-slate-700 hover:bg-slate-100 dark:text-slate-300 dark:hover:bg-slate-800',
        )}
        style={{ paddingLeft: 8 + depth * 14 }}
        title={node.path}
      >
        {isDir ? (
          <span className="inline-flex text-slate-400 [&>svg]:h-3.5 [&>svg]:w-3.5">{open ? <ChevronDown /> : <ChevronRight />}</span>
        ) : (
          <span className="inline-block w-3.5" />
        )}
        <span className="inline-flex shrink-0 [&>svg]:h-4 [&>svg]:w-4">{isDir ? open ? <FolderOpen className="text-gold-500" /> : <Folder className="text-gold-500" /> : fileIcon(node.name)}</span>
        <span className="min-w-0 flex-1 truncate">{node.name}</span>
        {fresh ? <span className="h-1.5 w-1.5 shrink-0 rounded-full bg-emerald-500" title="刚刚更新" /> : null}
        {isDir ? <span className="shrink-0 text-[11px] tabular text-steel-400">{node.children?.length ?? 0}</span> : node.size != null ? <span className="shrink-0 text-[11px] tabular text-steel-400">{formatBytes(node.size)}</span> : null}
      </button>
      {isDir && open ? node.children?.map((c) => <TreeNode key={c.path} node={c} depth={depth + 1} expanded={expanded} toggle={toggle} selected={selected} onSelect={onSelect} now={now} />) : null}
    </div>
  );
}

export function ArtifactsTab({
  project,
  refreshKey,
  selectedPath,
  onSelectPath,
  now,
}: {
  project: Project;
  refreshKey: number;
  selectedPath: string | null;
  onSelectPath: (path: string | null) => void;
  now: number;
}) {
  const { data: tree, loading, reload } = useAsync(() => api.getArtifacts(project.id), [project.id, refreshKey]);
  const files = useMemo(() => (tree ? flatten(tree) : []), [tree]);
  const byPath = useMemo(() => new Map(files.map((f) => [f.path, f])), [files]);
  const [query, setQuery] = useState('');
  const [expanded, setExpanded] = useState<Set<string>>(new Set());

  // 默认展开一级目录与所选文件所在目录
  useEffect(() => {
    if (!tree) return;
    setExpanded((prev) => {
      const next = new Set(prev);
      if (prev.size === 0) tree.children?.forEach((c) => c.type === 'dir' && next.add(c.path));
      if (selectedPath) {
        let dir = parentDir(selectedPath);
        while (dir) {
          next.add(dir);
          dir = parentDir(dir);
        }
      }
      return next;
    });
  }, [tree, selectedPath]);

  // 默认选中：路演稿 → PDF → 讲稿 → 第一个 Markdown
  useEffect(() => {
    if (!tree || selectedPath) return;
    const a = project.artifacts;
    const candidates = [a.deck, a.pdf, a.script, a.model, files.find((f) => f.path.startsWith('05_') && fileExt(f.name) === 'md')?.path, files.find((f) => fileExt(f.name) === 'md')?.path, files[0]?.path];
    const first = candidates.find((p): p is string => Boolean(p && byPath.has(p)));
    if (first) onSelectPath(first);
  }, [tree, selectedPath, project.artifacts, files, byPath, onSelectPath]);

  const selectedNode = selectedPath ? byPath.get(selectedPath) ?? (selectedPath ? { name: baseName(selectedPath), path: selectedPath, type: 'file' as const } : null) : null;

  const toggle = (p: string) =>
    setExpanded((prev) => {
      const next = new Set(prev);
      if (next.has(p)) next.delete(p);
      else next.add(p);
      return next;
    });

  const matches = useMemo(() => {
    const q = query.trim().toLowerCase();
    if (!q) return null;
    return files.filter((f) => f.path.toLowerCase().includes(q)).slice(0, 60);
  }, [files, query]);

  const latest = useMemo(() => files.reduce<string | null>((acc, f) => (f.mtime && (!acc || Date.parse(f.mtime) > Date.parse(acc)) ? f.mtime : acc), null), [files]);

  return (
    <div className="grid gap-4 xl:grid-cols-[360px_minmax(0,1fr)]">
      <div className="flex min-h-0 flex-col gap-4 xl:h-[calc(100vh-22rem)] xl:min-h-[460px]">
        <Card className="shrink-0 p-3">
          <div className="mb-2 flex items-center justify-between px-1">
            <h3 className="text-sm font-semibold text-slate-800 dark:text-slate-100">关键产出</h3>
            <span className="text-[11px] text-steel-500">{Object.keys(project.artifacts).length} / {QUICK.length} 已生成</span>
          </div>
          <div className="grid grid-cols-2 gap-2">
            {QUICK.map((q) => {
              const path = project.artifacts[q.key];
              const active = path != null && selectedPath === path;
              return (
                <button
                  key={q.key}
                  disabled={!path}
                  onClick={() => path && onSelectPath(path)}
                  className={cn(
                    'flex items-center gap-2.5 rounded-lg border px-2.5 py-2 text-left transition',
                    !path
                      ? 'cursor-not-allowed border-dashed border-slate-200 text-slate-400 dark:border-slate-700 dark:text-slate-600'
                      : active
                        ? 'border-brand-300 bg-brand-50 dark:border-brand-400/60 dark:bg-brand-900/40'
                        : 'border-slate-200 bg-white hover:border-brand-200 hover:bg-brand-50/40 dark:border-surface-dark-border dark:bg-surface-dark-card dark:hover:border-brand-500/50',
                  )}
                  title={path ?? '尚未生成'}
                >
                  <span className={cn('inline-flex h-8 w-8 shrink-0 items-center justify-center rounded-md [&>svg]:h-4 [&>svg]:w-4', path ? 'bg-gold-50 text-gold-700 dark:bg-gold-900/30 dark:text-gold-300' : 'bg-slate-100 text-slate-400 dark:bg-slate-800')}>{q.icon}</span>
                  <span className="min-w-0">
                    <span className={cn('block truncate text-xs font-medium', path ? 'text-slate-800 dark:text-slate-100' : '')}>{q.label}</span>
                    <span className="block truncate text-[11px] text-steel-500">{path ? baseName(path) : '未生成'}</span>
                  </span>
                </button>
              );
            })}
          </div>
        </Card>

        <Card className="flex max-h-[480px] min-h-[260px] flex-1 flex-col xl:max-h-none">
          <div className="flex items-center gap-2 border-b border-slate-100 px-3 py-2 dark:border-surface-dark-border">
            <div className="relative min-w-0 flex-1">
              <Search className="pointer-events-none absolute left-2.5 top-1/2 h-3.5 w-3.5 -translate-y-1/2 text-slate-400" />
              <Input value={query} onChange={(e) => setQuery(e.target.value)} placeholder="搜索文件…" className="h-8 py-1 pl-8 text-xs" />
            </div>
            <span className="text-[11px] tabular text-steel-500">{files.length} 个文件</span>
            <Button size="xs" variant="ghost" icon={<RefreshCw className={cn(loading && 'animate-spin')} />} onClick={() => void reload()} title="刷新目录" aria-label="刷新目录" />
          </div>
          <div className="min-h-0 flex-1 overflow-y-auto p-2">
            {loading && !tree ? (
              <div className="flex h-full items-center justify-center text-sm text-steel-500">
                <Loader2 className="mr-2 h-4 w-4 animate-spin" />
                加载目录…
              </div>
            ) : !tree || files.length === 0 ? (
              <EmptyState compact icon={<Folder />} title="工作区还是空的" description="运行后，任务书、调研笔记、财务模型、图表、正文、路演与提交包会出现在这里。" />
            ) : matches ? (
              matches.length === 0 ? (
                <EmptyState compact title="没有匹配的文件" />
              ) : (
                <div className="space-y-0.5">
                  {matches.map((f) => (
                    <button key={f.path} onClick={() => onSelectPath(f.path)} className={cn('flex w-full items-center gap-2 rounded-md px-2 py-1 text-left text-[13px] hover:bg-slate-100 dark:hover:bg-slate-800', selectedPath === f.path && 'bg-brand-50 dark:bg-brand-900/40')}>
                      <span className="inline-flex shrink-0 [&>svg]:h-4 [&>svg]:w-4">{fileIcon(f.name)}</span>
                      <span className="min-w-0 flex-1 truncate">{f.path}</span>
                    </button>
                  ))}
                </div>
              )
            ) : (
              tree.children?.map((c) => <TreeNode key={c.path} node={c} depth={0} expanded={expanded} toggle={toggle} selected={selectedPath} onSelect={(n) => onSelectPath(n.path)} now={now} />)
            )}
          </div>
          {latest ? <div className="border-t border-slate-100 px-3 py-1.5 text-[11px] text-steel-500 dark:border-surface-dark-border">最近更新 {formatRelative(latest, now)}</div> : null}
        </Card>
      </div>

      <Card className="h-[calc(100vh-22rem)] min-h-[460px] overflow-hidden">
        <ArtifactPreview projectId={project.id} node={selectedNode} refreshKey={refreshKey} />
      </Card>
    </div>
  );
}
