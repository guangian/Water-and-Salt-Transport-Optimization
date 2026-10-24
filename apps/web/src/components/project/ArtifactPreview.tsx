import { useEffect, useState } from 'react';
import { Download, ExternalLink, FileQuestion, Loader2, RefreshCw } from 'lucide-react';
import type { ArtifactNode, Preview, SheetData } from '../../types';
import { api } from '../../api';
import { useAsync } from '../../hooks/useAsync';
import { baseName, fileExt, formatBytes, formatDateTime } from '../../lib/format';
import { cn } from '../../lib/cn';
import { Button } from '../ui/Button';
import { EmptyState } from '../ui/Misc';
import { Markdown } from '../domain/Markdown';
import { SlidesPreview } from './SlidesPreview';

function cellText(v: string | number | null): string {
  if (v == null) return '';
  if (typeof v === 'number') {
    if (Number.isInteger(v)) return v.toLocaleString('zh-CN');
    return v.toLocaleString('zh-CN', { maximumFractionDigits: 4 });
  }
  return v;
}

export function SheetTable({ sheets }: { sheets: SheetData[] }) {
  const [idx, setIdx] = useState(0);
  useEffect(() => setIdx(0), [sheets]);
  const sheet = sheets[idx] ?? sheets[0];
  if (!sheet) return <EmptyState compact title="空表格" />;
  const [header, ...rows] = sheet.rows;
  return (
    <div className="flex h-full flex-col">
      {sheets.length > 1 ? (
        <div className="mb-3 flex flex-wrap gap-1">
          {sheets.map((s, i) => (
            <button
              key={s.name}
              onClick={() => setIdx(i)}
              className={cn('rounded-md px-3 py-1 text-xs font-medium transition', i === idx ? 'bg-brand-600 text-white dark:bg-brand-300 dark:text-brand-900' : 'bg-slate-100 text-slate-600 hover:bg-slate-200 dark:bg-slate-800 dark:text-slate-300 dark:hover:bg-slate-700')}
            >
              {s.name}
              <span className="ml-1 opacity-60">{Math.max(0, s.rows.length - 1)}</span>
            </button>
          ))}
        </div>
      ) : null}
      <div className="min-h-0 flex-1 overflow-auto rounded-lg border border-slate-200 dark:border-surface-dark-border">
        <table className="table-base min-w-full">
          <thead>
            <tr>
              <th className="w-10 text-right">#</th>
              {(header ?? []).map((h, i) => (
                <th key={i} className="whitespace-nowrap">
                  {cellText(h)}
                </th>
              ))}
            </tr>
          </thead>
          <tbody>
            {rows.map((r, i) => (
              <tr key={i}>
                <td className="text-right text-xs tabular text-steel-400">{i + 1}</td>
                {(header ?? r).map((_, j) => {
                  const v = r[j] ?? null;
                  return (
                    <td key={j} className={cn('whitespace-nowrap', typeof v === 'number' && 'text-right tabular')}>
                      {cellText(v)}
                    </td>
                  );
                })}
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
}

function PreviewBody({ preview, projectId, path, downloadUrl }: { preview: Preview; projectId: string; path: string; downloadUrl: string }) {
  switch (preview.kind) {
    case 'markdown':
      return (
        <div className="h-full overflow-auto px-6 py-5">
          <Markdown content={preview.content ?? ''} projectId={projectId} basePath={path} />
        </div>
      );
    case 'image':
      return (
        <div className="flex h-full items-center justify-center overflow-auto bg-[linear-gradient(45deg,#f1f5f9_25%,transparent_25%,transparent_75%,#f1f5f9_75%),linear-gradient(45deg,#f1f5f9_25%,transparent_25%,transparent_75%,#f1f5f9_75%)] bg-[length:20px_20px] bg-[position:0_0,10px_10px] p-6 dark:bg-none dark:bg-slate-900">
          <img src={preview.url} alt={baseName(path)} className="max-h-full max-w-full rounded-lg border border-slate-200 bg-white shadow-sm dark:border-surface-dark-border" />
        </div>
      );
    case 'slides':
      return (
        <div className="h-full overflow-auto px-6 py-5">
          {preview.slides ? <SlidesPreview deck={preview.slides} projectId={projectId} basePath={path} /> : <EmptyState compact title="slides.json 无法解析" />}
        </div>
      );
    case 'table':
      return <div className="h-full px-5 py-4">{preview.sheets?.length ? <SheetTable sheets={preview.sheets} /> : <EmptyState compact title="空表格" />}</div>;
    case 'pdf':
      return <iframe src={preview.url} title={baseName(path)} className="h-full w-full bg-slate-100 dark:bg-slate-900" />;
    case 'text':
      return (
        <pre className="h-full overflow-auto whitespace-pre-wrap break-words px-6 py-5 font-mono text-xs leading-relaxed text-slate-700 dark:text-slate-300">{preview.content ?? ''}</pre>
      );
    default:
      return (
        <EmptyState
          icon={<FileQuestion />}
          title="此文件不支持在线预览"
          description={preview.message ?? `${fileExt(path).toUpperCase() || '该'} 文件请下载后用本地软件打开。`}
          action={
            <a href={downloadUrl} download={baseName(path)}>
              <Button variant="secondary" icon={<Download />}>
                下载文件{preview.size ? `（${formatBytes(preview.size)}）` : ''}
              </Button>
            </a>
          }
        />
      );
  }
}

export function ArtifactPreview({ projectId, node, refreshKey }: { projectId: string; node: ArtifactNode | null; refreshKey: number }) {
  const path = node?.path ?? '';
  const { data, loading, error, reload } = useAsync(() => api.getPreview(projectId, path), [projectId, path, refreshKey], { enabled: Boolean(node) });
  const downloadUrl = node ? api.fileUrl(projectId, node.path) : '#';

  if (!node) {
    return (
      <div className="flex h-full items-center justify-center">
        <EmptyState icon={<FileQuestion />} title="选择一个文件开始预览" description="支持 Markdown（含图表）、图片、slides.json 路演稿、Excel/CSV 表格、PDF 与文本文件。" />
      </div>
    );
  }
  return (
    <div className="flex h-full flex-col">
      <div className="flex items-center gap-3 border-b border-slate-100 px-4 py-2.5 dark:border-surface-dark-border">
        <div className="min-w-0 flex-1">
          <div className="truncate text-sm font-semibold text-slate-800 dark:text-slate-100" title={node.path}>
            {node.name}
          </div>
          <div className="truncate text-[11px] text-steel-500">
            {node.path}
            {node.size != null ? ` · ${formatBytes(node.size)}` : ''}
            {node.mtime ? ` · ${formatDateTime(node.mtime)}` : ''}
          </div>
        </div>
        <Button size="xs" variant="ghost" icon={<RefreshCw className={cn(loading && 'animate-spin')} />} onClick={() => void reload()} title="重新加载">
          刷新
        </Button>
        <a href={downloadUrl} target="_blank" rel="noreferrer" className="hidden sm:block">
          <Button size="xs" variant="ghost" icon={<ExternalLink />}>
            新标签打开
          </Button>
        </a>
        <a href={downloadUrl} download={node.name}>
          <Button size="xs" variant="secondary" icon={<Download />}>
            下载
          </Button>
        </a>
      </div>
      <div className="relative min-h-0 flex-1">
        {loading && !data ? (
          <div className="flex h-full items-center justify-center text-sm text-steel-500">
            <Loader2 className="mr-2 h-4 w-4 animate-spin" />
            加载预览…
          </div>
        ) : error ? (
          <EmptyState compact icon={<FileQuestion />} title="预览加载失败" description={error.message} action={<Button variant="secondary" size="sm" onClick={() => void reload()}>重试</Button>} />
        ) : data ? (
          // 按路径重建子树，切换文件时滚动位置归零
          <div key={node.path} className="h-full">
            <PreviewBody preview={data} projectId={projectId} path={node.path} downloadUrl={downloadUrl} />
          </div>
        ) : null}
      </div>
    </div>
  );
}
