import { useRef, useState, type DragEvent } from 'react';
import { FileText, Image as ImageIcon, Plus, Sheet, Trash2, UploadCloud, UserRound, ShieldCheck } from 'lucide-react';
import type { ProjectInput, SubjectType, Track } from '../../types';
import type { Draft, Errors } from './draft';
import { Field, Input, Select, Switch, Textarea } from '../ui/Form';
import { Button } from '../ui/Button';
import { formatBytes } from '../../lib/format';
import { cn } from '../../lib/cn';

const SUBJECT_TYPES: { value: SubjectType; label: string }[] = [
  { value: 'brand', label: '品牌' },
  { value: 'startup', label: '创业项目' },
  { value: 'case-company', label: '案例企业' },
  { value: 'event', label: '活动 / 会奖项目' },
];

export interface PendingFile {
  id: string;
  file: File;
}

function kindIcon(name: string) {
  const ext = name.split('.').pop()?.toLowerCase() ?? '';
  if (['png', 'jpg', 'jpeg', 'webp', 'gif'].includes(ext)) return <ImageIcon className="h-4 w-4 text-violet-500" />;
  if (['xlsx', 'xls', 'csv'].includes(ext)) return <Sheet className="h-4 w-4 text-emerald-600" />;
  return <FileText className="h-4 w-4 text-brand-500" />;
}

export function FileList({ files, uploaded, onRemove }: { files: PendingFile[]; uploaded?: ProjectInput[]; onRemove?: (id: string) => void }) {
  if (!files.length && !uploaded?.length) return null;
  return (
    <ul className="mt-3 divide-y divide-slate-100 rounded-lg border border-slate-200 dark:divide-surface-dark-border dark:border-surface-dark-border">
      {uploaded?.map((u) => (
        <li key={`u-${u.name}`} className="flex items-center gap-3 px-3 py-2 text-sm">
          {kindIcon(u.name)}
          <span className="min-w-0 flex-1 truncate text-slate-700 dark:text-slate-200">{u.name}</span>
          <span className="text-xs text-steel-500">{formatBytes(u.size)}</span>
          <span className={cn('chip', u.prereadStatus === 'done' ? 'text-emerald-700 dark:text-emerald-300' : u.prereadStatus === 'failed' ? 'text-red-600' : 'text-amber-700')}>
            {u.prereadStatus === 'done' ? (u.prereadChars ? `已预读 ${(u.prereadChars / 10000).toFixed(1)} 万字` : '已上传') : u.prereadStatus === 'failed' ? '预读失败' : '预读中…'}
          </span>
        </li>
      ))}
      {files.map((f) => (
        <li key={f.id} className="flex items-center gap-3 px-3 py-2 text-sm">
          {kindIcon(f.file.name)}
          <span className="min-w-0 flex-1 truncate text-slate-700 dark:text-slate-200">{f.file.name}</span>
          <span className="text-xs text-steel-500">{formatBytes(f.file.size)}</span>
          <span className="chip">待上传</span>
          {onRemove ? (
            <button type="button" onClick={() => onRemove(f.id)} className="rounded p-1 text-slate-400 hover:bg-red-50 hover:text-red-600 dark:hover:bg-red-900/30" aria-label="移除">
              <Trash2 className="h-3.5 w-3.5" />
            </button>
          ) : null}
        </li>
      ))}
    </ul>
  );
}

export function Dropzone({ onFiles, compact }: { onFiles: (files: File[]) => void; compact?: boolean }) {
  const [over, setOver] = useState(false);
  const inputRef = useRef<HTMLInputElement>(null);
  const onDrop = (e: DragEvent) => {
    e.preventDefault();
    setOver(false);
    const files = Array.from(e.dataTransfer.files ?? []);
    if (files.length) onFiles(files);
  };
  return (
    <div
      onDragOver={(e) => {
        e.preventDefault();
        setOver(true);
      }}
      onDragLeave={() => setOver(false)}
      onDrop={onDrop}
      onClick={() => inputRef.current?.click()}
      className={cn(
        'flex cursor-pointer flex-col items-center justify-center rounded-xl border-2 border-dashed text-center transition',
        compact ? 'px-4 py-5' : 'px-6 py-8',
        over ? 'border-gold-400 bg-gold-50/60 dark:bg-gold-900/10' : 'border-slate-300 bg-slate-50/60 hover:border-brand-300 hover:bg-brand-50/40 dark:border-slate-600 dark:bg-slate-900/40 dark:hover:bg-slate-800/60',
      )}
    >
      <input
        ref={inputRef}
        type="file"
        multiple
        className="hidden"
        accept=".pdf,.doc,.docx,.xlsx,.xls,.csv,.png,.jpg,.jpeg,.webp,.txt,.md,.pptx"
        onChange={(e) => {
          const files = Array.from(e.target.files ?? []);
          if (files.length) onFiles(files);
          e.target.value = '';
        }}
      />
      <UploadCloud className={cn('text-brand-400', compact ? 'h-6 w-6' : 'h-8 w-8')} />
      <p className="mt-2 text-sm font-medium text-slate-700 dark:text-slate-200">拖拽文件到此处，或点击选择</p>
      <p className="mt-1 text-xs text-steel-500 dark:text-steel-400">支持 PDF / Word / Excel / CSV / 图片；上传后自动预读为文本供智能体引用</p>
    </div>
  );
}

export function StepInfo({
  track,
  draft,
  onChange,
  errors,
  files,
  onAddFiles,
  onRemoveFile,
}: {
  track: Track | null;
  draft: Draft;
  onChange: (patch: Partial<Draft>) => void;
  errors: Errors;
  files: PendingFile[];
  onAddFiles: (files: File[]) => void;
  onRemoveFile: (id: string) => void;
}) {
  const setSubject = (patch: Partial<Draft['subject']>) => onChange({ subject: { ...draft.subject, ...patch } });
  const cases = track?.assignedCases ?? [];
  const briefLen = draft.subject.brief.length;

  return (
    <div className="grid grid-cols-1 gap-6 xl:grid-cols-3">
      <div className="space-y-6 xl:col-span-2">
        <section className="card p-5">
          <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">作品与团队</h2>
          <div className="mt-4 grid grid-cols-1 gap-4 md:grid-cols-2">
            <Field label="作品名" required error={errors.name ?? null} hint="建议“主标题（意象）——品牌/项目 + 文体”，将写入封面与文件名。" className="md:col-span-2">
              <Input value={draft.name} onChange={(e) => onChange({ name: e.target.value })} placeholder="例：宋韵天青——弘宝汝瓷品牌策划书" maxLength={80} invalid={Boolean(errors.name)} />
            </Field>
            <Field label="团队名" required error={errors.teamName ?? null} hint="匿名制：团队名代替院校出现在封面。">
              <Input value={draft.teamName} onChange={(e) => onChange({ teamName: e.target.value })} placeholder="例：青瓷工作室" maxLength={40} invalid={Boolean(errors.teamName)} />
            </Field>
            <Field label="策划对象类型" required>
              <Select value={draft.subject.type} onChange={(e) => setSubject({ type: e.target.value as SubjectType })}>
                {SUBJECT_TYPES.map((t) => (
                  <option key={t.value} value={t.value}>
                    {t.label}
                  </option>
                ))}
              </Select>
            </Field>
          </div>
        </section>

        <section className="card p-5">
          <div className="flex flex-wrap items-center justify-between gap-3">
            <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">策划对象 / 项目</h2>
            <Switch
              checked={draft.subject.isAssignedCase}
              onChange={(v) => setSubject({ isAssignedCase: v, caseTitle: v ? draft.subject.caseTitle : '' })}
              label={<span className="text-sm">命题企业 / 命题组</span>}
              disabled={!cases.length}
            />
          </div>
          {draft.subject.isAssignedCase && cases.length ? (
            <Field label="选择命题" required error={errors['subject.caseTitle'] ?? null} className="mt-4" hint="来自组委会公布的命题企业与方向；选择后自动填充企业名。">
              <Select
                value={draft.subject.caseTitle}
                invalid={Boolean(errors['subject.caseTitle'])}
                onChange={(e) => {
                  const c = cases.find((x) => x.title === e.target.value);
                  setSubject({ caseTitle: e.target.value, company: c?.company ?? draft.subject.company, name: draft.subject.name || (c?.title.replace(/品牌.*$|策划.*$|方案.*$/, '') ?? '') });
                }}
              >
                <option value="">请选择命题…</option>
                {cases.map((c) => (
                  <option key={c.title} value={c.title}>
                    {c.company} — {c.title}
                  </option>
                ))}
              </Select>
            </Field>
          ) : null}
          <div className="mt-4 grid grid-cols-1 gap-4 md:grid-cols-2">
            <Field label="品牌 / 项目名" required error={errors['subject.name'] ?? null}>
              <Input value={draft.subject.name} onChange={(e) => setSubject({ name: e.target.value })} placeholder="例：弘宝汝瓷" invalid={Boolean(errors['subject.name'])} />
            </Field>
            <Field label="企业 / 公司">
              <Input value={draft.subject.company} onChange={(e) => setSubject({ company: e.target.value })} placeholder="例：汝州市弘宝汝瓷坊有限公司" />
            </Field>
            <Field label="地区">
              <Input value={draft.subject.region} onChange={(e) => setSubject({ region: e.target.value })} placeholder="例：河南汝州" />
            </Field>
            <Field label="行业">
              <Input value={draft.subject.industry} onChange={(e) => setSubject({ industry: e.target.value })} placeholder="例：非遗工艺品 / 文创" />
            </Field>
            <Field
              label="项目简介"
              required
              className="md:col-span-2"
              error={errors['subject.brief'] ?? null}
              hint="写清楚现状、目标与约束（预算、人群、渠道等），智能体据此读题；越具体产出越贴合。"
              right={<span className={cn('text-xs tabular', briefLen < 20 ? 'text-amber-600' : 'text-steel-500')}>{briefLen} / 800 字</span>}
            >
              <Textarea
                value={draft.subject.brief}
                onChange={(e) => setSubject({ brief: e.target.value.slice(0, 800) })}
                rows={5}
                placeholder="例：弘宝汝瓷是汝州本地非遗传承品牌，以天青釉茶器与文创礼品为主，线下渠道为主、线上份额不足 3%。希望面向 25–40 岁新中产重塑品牌定位……"
                invalid={Boolean(errors['subject.brief'])}
              />
            </Field>
          </div>
        </section>

        <section className="card p-5">
          <div className="flex items-start justify-between gap-3">
            <div>
              <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">团队成员（匿名）</h2>
              <p className="mt-1 flex items-center gap-1.5 text-xs text-steel-500 dark:text-steel-400">
                <ShieldCheck className="h-3.5 w-3.5 text-emerald-500" />
                匿名制：只填角色 / 专业 / 特长，不填姓名与院校，系统据此撰写团队章节与分工。
              </p>
            </div>
            <Button variant="secondary" size="sm" icon={<Plus />} onClick={() => onChange({ team: [...draft.team, { role: '', major: '', strength: '' }] })} disabled={draft.team.length >= 8}>
              添加成员
            </Button>
          </div>
          {errors.team ? <p className="mt-2 text-xs text-red-600">{errors.team}</p> : null}
          <div className="mt-4 space-y-2">
            <div className="hidden grid-cols-[1.2fr_1fr_1.4fr_32px] gap-2 px-1 text-[11px] font-medium uppercase tracking-wide text-steel-500 md:grid">
              <span>角色</span>
              <span>专业</span>
              <span>特长</span>
              <span />
            </div>
            {draft.team.map((m, i) => (
              <div key={i} className="grid grid-cols-1 gap-2 md:grid-cols-[1.2fr_1fr_1.4fr_32px]">
                <div className="relative">
                  <UserRound className="pointer-events-none absolute left-2.5 top-1/2 h-4 w-4 -translate-y-1/2 text-slate-400" />
                  <Input className="pl-8" placeholder="角色，如 队长 / 总策划" value={m.role} onChange={(e) => onChange({ team: draft.team.map((x, j) => (j === i ? { ...x, role: e.target.value } : x)) })} />
                </div>
                <Input placeholder="专业，如 市场营销" value={m.major} onChange={(e) => onChange({ team: draft.team.map((x, j) => (j === i ? { ...x, major: e.target.value } : x)) })} />
                <Input placeholder="特长，如 品牌传播与策划" value={m.strength} onChange={(e) => onChange({ team: draft.team.map((x, j) => (j === i ? { ...x, strength: e.target.value } : x)) })} />
                <button
                  type="button"
                  className="inline-flex h-9 w-8 items-center justify-center rounded-md text-slate-400 hover:bg-red-50 hover:text-red-600 disabled:opacity-30 dark:hover:bg-red-900/30"
                  onClick={() => onChange({ team: draft.team.filter((_, j) => j !== i) })}
                  disabled={draft.team.length <= 1}
                  aria-label="删除成员"
                >
                  <Trash2 className="h-4 w-4" />
                </button>
              </div>
            ))}
          </div>
        </section>
      </div>

      <div className="space-y-6">
        <section className="card p-5">
          <h2 className="text-base font-semibold text-slate-800 dark:text-slate-100">资料上传</h2>
          <p className="mt-1 text-xs text-steel-500 dark:text-steel-400">企业资料、年报、经营数据、产品图片等。创建项目后将上传到 inputs/ 并自动预读。</p>
          <div className="mt-4">
            <Dropzone onFiles={onAddFiles} />
            <FileList files={files} onRemove={onRemoveFile} />
            {files.length ? <p className="mt-2 text-[11px] text-steel-500">共 {files.length} 个文件，{formatBytes(files.reduce((s, f) => s + f.file.size, 0))}。刷新页面后需重新添加文件。</p> : null}
          </div>
        </section>
        <section className="rounded-xl border border-gold-200 bg-gold-50/60 p-4 text-xs leading-relaxed text-gold-900 dark:border-gold-800 dark:bg-gold-900/10 dark:text-gold-100">
          <p className="font-medium">写好简介的三个要点</p>
          <ul className="mt-1.5 list-disc space-y-1 pl-4">
            <li>现状：规模、渠道结构、已有品牌资产或业务基础。</li>
            <li>目标：面向谁、要达到什么效果、预算范围。</li>
            <li>约束：命题要求、不能改动的事项、必须使用的资料。</li>
          </ul>
        </section>
      </div>
    </div>
  );
}
