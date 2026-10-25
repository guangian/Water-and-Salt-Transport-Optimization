import { useLayoutEffect, useRef, useState, type ReactNode } from 'react';
import { Bar, BarChart, CartesianGrid, Cell, Legend, Line, LineChart, Pie, PieChart, XAxis, YAxis } from 'recharts';
import type { Slide, SlideChart, SlidesDeck } from '../../types';
import { cn } from '../../lib/cn';

export const SLIDE_W = 960;
export const SLIDE_H = 540;

const NAVY = '#14305C';
const GOLD = '#C9A227';
const STEEL = '#5B6B82';
const PALETTE = [NAVY, GOLD, STEEL, '#2E8B6E', '#B5533C', '#7F97BF'];

export const LAYOUT_LABELS: Record<Slide['layout'], string> = {
  cover: '封面',
  agenda: '目录',
  section: '章节页',
  bullets: '要点',
  image: '图片',
  two_column: '双栏',
  chart: '图表',
  table: '表格',
  kpi: '指标',
  timeline: '时间线',
  team: '团队',
  quote: '引言',
  closing: '结尾',
  appendix: '附录',
};

/** 把 960×540 的逻辑画布按容器宽度等比缩放 */
export function ScaledSlide({ children, className, onClick, title }: { children: ReactNode; className?: string; onClick?: () => void; title?: string }) {
  const ref = useRef<HTMLDivElement>(null);
  const [scale, setScale] = useState(0);
  useLayoutEffect(() => {
    const el = ref.current;
    if (!el) return;
    const update = () => setScale(el.clientWidth / SLIDE_W);
    update();
    const ro = new ResizeObserver(update);
    ro.observe(el);
    return () => ro.disconnect();
  }, []);
  return (
    <div ref={ref} onClick={onClick} title={title} className={cn('relative aspect-video w-full overflow-hidden bg-white', onClick && 'cursor-zoom-in', className)}>
      <div style={{ width: SLIDE_W, height: SLIDE_H, transform: `scale(${scale})`, transformOrigin: 'top left', visibility: scale ? 'visible' : 'hidden' }}>{children}</div>
    </div>
  );
}

function Title({ children, kicker }: { children: ReactNode; kicker?: string }) {
  return (
    <div className="mb-5">
      {kicker ? <div className="mb-1 text-[14px] font-semibold tracking-[0.2em] text-[#C9A227]">{kicker}</div> : null}
      <div className="flex items-start gap-3">
        <span className="mt-1.5 h-7 w-1.5 shrink-0 rounded-sm bg-[#C9A227]" />
        <h2 className="text-balance text-[28px] font-bold leading-tight text-[#14305C]">{children}</h2>
      </div>
    </div>
  );
}

function Bullets({ items, size = 20 }: { items: string[]; size?: number }) {
  return (
    <ul className="space-y-3">
      {items.map((b, i) => (
        <li key={i} className="flex items-start gap-3 leading-snug text-slate-700" style={{ fontSize: size }}>
          <span className="mt-[0.45em] h-2 w-2 shrink-0 rounded-full bg-[#C9A227]" />
          <span>{b}</span>
        </li>
      ))}
    </ul>
  );
}

function Footer({ deck, index, total, dark }: { deck: SlidesDeck; index: number; total: number; dark?: boolean }) {
  return (
    <div className={cn('absolute inset-x-0 bottom-0 flex items-center justify-between px-10 pb-4 text-[12px]', dark ? 'text-white/60' : 'text-slate-400')}>
      <span className="truncate">{deck.footer ?? deck.title}</span>
      <span className="tabular">
        {index + 1} / {total}
      </span>
    </div>
  );
}

function SlideChartView({ chart }: { chart: SlideChart }) {
  const data = chart.categories.map((c, i) => {
    const row: Record<string, string | number> = { name: c };
    for (const s of chart.series) row[s.name] = s.values[i] ?? 0;
    return row;
  });
  const common = { width: 560, height: 330, margin: { top: 10, right: 10, left: 0, bottom: 0 } };
  if (chart.type === 'pie') {
    const s = chart.series[0];
    const pieData = chart.categories.map((c, i) => ({ name: c, value: s?.values[i] ?? 0 }));
    return (
      <PieChart {...common}>
        <Pie data={pieData} dataKey="value" nameKey="name" cx="50%" cy="50%" outerRadius={120} innerRadius={60} isAnimationActive={false} label={({ name, percent }) => `${name} ${Math.round((percent ?? 0) * 100)}%`}>
          {pieData.map((_, i) => (
            <Cell key={i} fill={PALETTE[i % PALETTE.length]} />
          ))}
        </Pie>
      </PieChart>
    );
  }
  if (chart.type === 'line') {
    return (
      <LineChart {...common} data={data}>
        <CartesianGrid vertical={false} stroke="#E2E8F0" />
        <XAxis dataKey="name" tick={{ fontSize: 13, fill: '#475569' }} axisLine={{ stroke: '#CBD5E1' }} tickLine={false} />
        <YAxis tick={{ fontSize: 12, fill: '#64748B' }} axisLine={false} tickLine={false} width={48} />
        <Legend wrapperStyle={{ fontSize: 13 }} />
        {chart.series.map((s, i) => (
          <Line key={s.name} type="monotone" dataKey={s.name} stroke={PALETTE[i % PALETTE.length]} strokeWidth={3} dot={{ r: 4 }} isAnimationActive={false} />
        ))}
      </LineChart>
    );
  }
  return (
    <BarChart {...common} data={data} barCategoryGap="28%">
      <CartesianGrid vertical={false} stroke="#E2E8F0" />
      <XAxis dataKey="name" tick={{ fontSize: 13, fill: '#475569' }} axisLine={{ stroke: '#CBD5E1' }} tickLine={false} />
      <YAxis tick={{ fontSize: 12, fill: '#64748B' }} axisLine={false} tickLine={false} width={48} />
      <Legend wrapperStyle={{ fontSize: 13 }} />
      {chart.series.map((s, i) => (
        <Bar key={s.name} dataKey={s.name} fill={PALETTE[i % PALETTE.length]} radius={[3, 3, 0, 0]} isAnimationActive={false} />
      ))}
    </BarChart>
  );
}

export function SlideCanvas({ slide, index, total, deck, resolveImage }: { slide: Slide; index: number; total: number; deck: SlidesDeck; resolveImage: (src: string) => string }) {
  const items = slide.bullets ?? slide.items ?? [];
  const body = (() => {
    switch (slide.layout) {
      case 'cover':
        return (
          <div className="relative h-full w-full overflow-hidden bg-gradient-to-br from-[#10274B] via-[#14305C] to-[#0C1E3A] px-16 pt-24 text-white">
            <div className="absolute -right-24 -top-24 h-80 w-80 rounded-full border-[18px] border-[#C9A227]/20" />
            <div className="absolute -bottom-28 right-24 h-72 w-72 rounded-full border-[14px] border-[#C9A227]/10" />
            <div className="mb-6 h-1.5 w-20 rounded bg-[#C9A227]" />
            <h1 className="max-w-[720px] text-balance text-[46px] font-bold leading-tight">{slide.title}</h1>
            {slide.subtitle ? <p className="mt-6 text-[22px] text-[#EBD88E]">{slide.subtitle}</p> : null}
            <div className="absolute bottom-10 left-16 text-[14px] text-white/60">{deck.title}</div>
          </div>
        );
      case 'section':
        return (
          <div className="relative flex h-full w-full flex-col justify-center bg-[#14305C] px-16 text-white">
            <div className="absolute left-0 top-0 h-full w-3 bg-[#C9A227]" />
            {slide.kicker ? <div className="mb-3 text-[16px] font-semibold tracking-[0.3em] text-[#EBD88E]">{slide.kicker}</div> : null}
            <h2 className="text-balance text-[48px] font-bold leading-tight">{slide.title}</h2>
            {slide.subtitle ? <p className="mt-4 text-[20px] text-white/70">{slide.subtitle}</p> : null}
          </div>
        );
      case 'closing':
        return (
          <div className="relative flex h-full w-full flex-col items-center justify-center bg-gradient-to-br from-[#14305C] to-[#081528] text-center text-white">
            <div className="mb-6 h-1.5 w-16 rounded bg-[#C9A227]" />
            <h2 className="text-[52px] font-bold">{slide.title ?? '谢谢'}</h2>
            {slide.subtitle ? <p className="mt-4 text-[22px] text-[#EBD88E]">{slide.subtitle}</p> : null}
            {slide.contact ? <p className="mt-6 text-[16px] text-white/60">{slide.contact}</p> : null}
            <div className="absolute bottom-8 text-[13px] text-white/50">{deck.footer ?? deck.title}</div>
          </div>
        );
      case 'agenda':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title ?? '目录'}</Title>
            <ol className="grid grid-cols-2 gap-x-10 gap-y-4 pt-2">
              {items.map((it, i) => (
                <li key={i} className="flex items-center gap-4 border-b border-slate-200 pb-3 text-[22px] text-slate-700">
                  <span className="w-10 text-[26px] font-bold tabular text-[#C9A227]">{String(i + 1).padStart(2, '0')}</span>
                  {it}
                </li>
              ))}
            </ol>
          </div>
        );
      case 'bullets':
      case 'appendix':
        return (
          <div className="h-full px-12 pt-10">
            <Title kicker={slide.kicker}>{slide.title}</Title>
            <div className={cn('grid gap-8', slide.highlight ? 'grid-cols-[1fr_280px]' : 'grid-cols-1')}>
              <Bullets items={items} />
              {slide.highlight ? (
                <div className="flex h-fit flex-col items-center justify-center rounded-2xl border border-[#EBD88E] bg-[#FBF7E9] px-6 py-8 text-center">
                  <div className="text-[44px] font-bold leading-none tabular text-[#14305C]">{slide.highlight}</div>
                  {slide.highlightLabel ? <div className="mt-3 text-[15px] text-[#846818]">{slide.highlightLabel}</div> : null}
                </div>
              ) : null}
            </div>
          </div>
        );
      case 'image':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <div className="flex gap-8">
              <div className="flex-1">
                {slide.image ? <img src={resolveImage(slide.image)} alt={slide.caption ?? ''} className="mx-auto max-h-[300px] rounded-lg border border-slate-200 object-contain" /> : null}
                {slide.caption ? <p className="mt-2 text-center text-[12px] text-slate-500">{slide.caption}</p> : null}
              </div>
              {items.length ? (
                <div className="w-[260px] pt-2">
                  <Bullets items={items} size={17} />
                </div>
              ) : null}
            </div>
          </div>
        );
      case 'two_column':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <div className="grid grid-cols-2 gap-8">
              {[slide.left, slide.right].map((col, i) =>
                col ? (
                  <div key={i} className={cn('rounded-xl border p-6', i === 0 ? 'border-slate-200 bg-slate-50' : 'border-[#EBD88E] bg-[#FBF7E9]')}>
                    <h3 className={cn('mb-4 border-b-2 pb-2 text-[22px] font-bold', i === 0 ? 'border-slate-300 text-slate-700' : 'border-[#C9A227] text-[#14305C]')}>{col.heading}</h3>
                    <Bullets items={col.bullets} size={18} />
                  </div>
                ) : null,
              )}
            </div>
          </div>
        );
      case 'chart':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <div className="flex items-start gap-6">
              <div className="shrink-0">{slide.chart ? <SlideChartView chart={slide.chart} /> : null}</div>
              <div className="flex-1 pt-6">
                {slide.chart?.unit ? <div className="mb-3 text-[13px] text-slate-500">单位：{slide.chart.unit}</div> : null}
                {slide.insight ? (
                  <div className="rounded-xl border-l-4 border-[#C9A227] bg-[#FBF7E9] px-5 py-4 text-[19px] font-medium leading-snug text-[#14305C]">{slide.insight}</div>
                ) : null}
              </div>
            </div>
          </div>
        );
      case 'table':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <table className="w-full border-collapse text-[17px]">
              <thead>
                <tr>
                  {(slide.columns ?? []).map((c, i) => (
                    <th key={i} className="bg-[#14305C] px-4 py-2.5 text-left font-semibold text-white first:rounded-tl-lg last:rounded-tr-lg">
                      {c}
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {(slide.rows ?? []).map((r, i) => (
                  <tr key={i} className={i % 2 ? 'bg-slate-50' : 'bg-white'}>
                    {r.map((c, j) => (
                      <td key={j} className={cn('border-b border-slate-200 px-4 py-2.5 text-slate-700', j === r.length - 1 && 'font-semibold text-[#14305C]')}>
                        {c}
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
            {slide.note ? <p className="mt-3 text-[12px] text-slate-400">{slide.note}</p> : null}
          </div>
        );
      case 'kpi':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <div className={cn('grid gap-5 pt-4', (slide.kpis?.length ?? 0) <= 3 ? 'grid-cols-3' : 'grid-cols-4')}>
              {(slide.kpis ?? []).map((k, i) => (
                <div key={i} className="rounded-2xl border border-slate-200 bg-gradient-to-b from-white to-slate-50 px-5 py-8 text-center shadow-sm">
                  <div className="text-[40px] font-bold leading-none tabular text-[#14305C]">{k.value}</div>
                  <div className="mt-4 text-[15px] leading-snug text-slate-500">{k.label}</div>
                </div>
              ))}
            </div>
          </div>
        );
      case 'timeline':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <div className="relative mt-10 flex items-start">
              <div className="absolute left-0 right-0 top-[18px] h-1 bg-slate-200" />
              {(slide.steps ?? []).map((s, i) => (
                <div key={i} className="relative flex-1 px-3 text-center">
                  <div className="mx-auto flex h-10 w-10 items-center justify-center rounded-full border-4 border-white bg-[#C9A227] text-[15px] font-bold text-[#14305C] shadow">{i + 1}</div>
                  <div className="mt-4 text-[17px] font-semibold text-[#14305C]">{s.when}</div>
                  <div className="mt-2 text-[15px] leading-snug text-slate-600">{s.what}</div>
                </div>
              ))}
            </div>
          </div>
        );
      case 'team':
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <div className={cn('grid gap-5 pt-2', (slide.members?.length ?? 0) <= 3 ? 'grid-cols-3' : 'grid-cols-4')}>
              {(slide.members ?? []).map((m, i) => (
                <div key={i} className="rounded-2xl border border-slate-200 bg-slate-50 px-5 py-6 text-center">
                  <div className="mx-auto mb-4 flex h-16 w-16 items-center justify-center rounded-full bg-[#14305C] text-[22px] font-bold text-[#EBD88E]">{m.role.slice(0, 1)}</div>
                  <div className="text-[18px] font-semibold text-[#14305C]">{m.role}</div>
                  <div className="mt-2 text-[14px] leading-snug text-slate-500">{m.duty}</div>
                </div>
              ))}
            </div>
          </div>
        );
      case 'quote':
        return (
          <div className="flex h-full flex-col items-center justify-center px-20 text-center">
            {slide.title ? <div className="mb-6 text-[15px] font-semibold tracking-[0.3em] text-[#C9A227]">{slide.title}</div> : null}
            <div className="text-[80px] leading-none text-[#EBD88E]">“</div>
            <p className="-mt-6 text-balance text-[40px] font-bold leading-snug text-[#14305C]">{slide.quote}</p>
            {slide.by ? <p className="mt-8 text-[18px] text-slate-500">—— {slide.by}</p> : null}
          </div>
        );
      default:
        return (
          <div className="h-full px-12 pt-10">
            <Title>{slide.title}</Title>
            <Bullets items={items} />
          </div>
        );
    }
  })();
  const dark = slide.layout === 'cover' || slide.layout === 'section' || slide.layout === 'closing';
  return (
    <div className="relative h-full w-full overflow-hidden bg-white font-sans text-slate-800" style={{ width: SLIDE_W, height: SLIDE_H }}>
      {body}
      {slide.layout !== 'cover' && slide.layout !== 'closing' ? <Footer deck={deck} index={index} total={total} dark={dark} /> : null}
    </div>
  );
}
