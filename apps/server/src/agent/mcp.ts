import path from 'node:path';
import { createSdkMcpServer, tool } from '@anthropic-ai/claude-agent-sdk';
import { z } from 'zod';
import type { StageId, StageStatus } from '../types.js';
import { STAGE_IDS } from '../types.js';
import { readJsonOrNull, writeJsonAtomic } from '../util/fs.js';
import { createLogger, errorMessage } from '../util/log.js';
import { fetchPage } from '../services/fetch.js';
import type { PythonRunner } from '../services/python.js';
import type { SearchService } from '../services/search.js';

const log = createLogger('mcp');

export const NUMBERS_FILE = path.join('03_财务', 'numbers.json');

export interface NumberEntry {
  key: string;
  label: string;
  value: number | string;
  unit?: string;
  source?: string;
  usedIn?: string[];
  synthetic?: boolean;
  note?: string;
}

export interface NumbersLedger {
  generatedAt: string;
  unit?: string;
  entries: NumberEntry[];
}

/** 每个项目运行一个绑定该项目上下文的 MCP server 实例 */
export interface McpContext {
  projectId: string;
  workspaceDir: string;
  search: SearchService;
  python: PythonRunner;
  updateStage: (stage: StageId, status: StageStatus, note?: string, round?: number) => Promise<void>;
  onFileWritten?: (relPath: string, size: number) => void;
}

function text(s: string, isError = false) {
  return { content: [{ type: 'text' as const, text: s }], ...(isError ? { isError: true } : {}) };
}

/** 合并写入数字账本：按 key 去重（后写覆盖），保持原顺序 */
export async function mergeNumbers(workspaceDir: string, entries: NumberEntry[], unit?: string): Promise<{ file: string; total: number; added: number; replaced: number }> {
  const file = path.join(workspaceDir, NUMBERS_FILE);
  const existing = (await readJsonOrNull<NumbersLedger>(file)) ?? { generatedAt: '', entries: [] };
  const list = Array.isArray(existing.entries) ? [...existing.entries] : [];
  let added = 0;
  let replaced = 0;
  for (const e of entries) {
    const idx = list.findIndex((x) => x.key === e.key);
    const merged: NumberEntry = idx >= 0 ? { ...list[idx]!, ...e } : e;
    if (idx >= 0) {
      list[idx] = merged;
      replaced++;
    } else {
      list.push(merged);
      added++;
    }
  }
  const ledger: NumbersLedger = { ...existing, generatedAt: new Date().toISOString(), unit: unit ?? existing.unit ?? '万元', entries: list };
  await writeJsonAtomic(file, ledger);
  return { file: NUMBERS_FILE.split(path.sep).join('/'), total: list.length, added, replaced };
}

export function createBizgoldMcpServer(ctx: McpContext) {
  const webSearch = tool(
    'web_search',
    '联网检索（中文优先）。按设置顺序尝试博查 / Tavily / Serper / Bing 网页 / DuckDuckGo，返回 JSON 数组 [{title,url,snippet,date?,provider}]，单次最多 10 条，结果会缓存 24 小时。' +
      '参数：query 检索词（必填，建议含年份与地区）；count 条数 1–10（默认取设置值）；site 限定站点域名（如 stats.gov.cn）；freshness 时效 day|week|month|year|noLimit。' +
      '所有提供商都失败时返回空数组并说明原因，请换关键词重试或改用 web_fetch 直接抓取已知来源。',
    {
      query: z.string().min(1).max(200).describe('检索词'),
      count: z.number().int().min(1).max(10).optional().describe('返回条数，1–10'),
      site: z.string().max(100).optional().describe('限定站点域名，例如 stats.gov.cn'),
      freshness: z.enum(['day', 'week', 'month', 'year', 'noLimit']).optional().describe('时效范围'),
    },
    async (args) => {
      try {
        const out = await ctx.search.search(args.query, { count: args.count, site: args.site, freshness: args.freshness });
        if (out.results.length === 0) {
          const why = out.errors.length ? out.errors.map((e) => `${e.provider}: ${e.error}`).join('；') : '没有可用的搜索提供商（请在设置中配置搜索 API key）';
          return text(JSON.stringify({ results: [], provider: null, note: `检索无结果。${why}` }, null, 2));
        }
        return text(JSON.stringify({ provider: out.provider, cached: out.cached, results: out.results }, null, 2));
      } catch (e) {
        log.warn('web_search failed', errorMessage(e));
        return text(JSON.stringify({ results: [], note: `检索失败：${errorMessage(e)}` }), true);
      }
    },
  );

  const webFetch = tool(
    'web_fetch',
    '抓取网页或 PDF 并抽取正文为 Markdown（保留标题、段落、列表、表格；去掉脚本/导航/页脚）。返回 {title,url,finalUrl,contentType,content,fetchedAt}。' +
      '参数：url（必填，http/https）；max_chars 正文最大字符数（默认 12000，最大 60000）。超时 20 秒，正文上限 2MB；拒绝内网/本机地址；PDF 自动用工具包抽文本。' +
      '抓到的数据引用时必须保留 URL 与访问日期。',
    {
      url: z.string().url().describe('要抓取的网址'),
      max_chars: z.number().int().min(500).max(60000).optional().describe('正文最大字符数，默认 12000'),
    },
    async (args) => {
      try {
        const page = await fetchPage(args.url, args.max_chars ?? 12000, { python: ctx.python, workspace: ctx.workspaceDir });
        return text(JSON.stringify(page, null, 2));
      } catch (e) {
        log.warn('web_fetch failed', errorMessage(e));
        return text(JSON.stringify({ url: args.url, error: `抓取失败：${errorMessage(e)}` }), true);
      }
    },
  );

  const updateProgress = tool(
    'update_progress',
    '上报阶段进度（每个阶段开始时 status=running、结束时 done/failed；可选 note 一句话说明、round 评审轮次）。' +
      '阶段 id：parse 读题解析 | research 调研 | strategy 策略 | finance 财务 | charts 图表 | write 撰写 | pitch 路演 | review 评审 | package 打包。服务端据此更新 project.json 并推送给前端。',
    {
      stage: z.enum(STAGE_IDS).describe('阶段 id'),
      status: z.enum(['running', 'done', 'failed', 'skipped', 'pending']).describe('阶段状态'),
      note: z.string().max(500).optional().describe('一句话说明（可选）'),
      round: z.number().int().min(0).max(20).optional().describe('评审轮次（review 阶段可选）'),
    },
    async (args) => {
      try {
        await ctx.updateStage(args.stage, args.status as StageStatus, args.note, args.round);
        return text(`已记录：${args.stage} → ${args.status}${args.round !== undefined ? `（第 ${args.round} 轮）` : ''}`);
      } catch (e) {
        return text(`进度记录失败：${errorMessage(e)}`, true);
      }
    },
  );

  const recordNumbers = tool(
    'record_numbers',
    '把关键数字登记到数字账本 03_财务/numbers.json（按 key 合并去重，同 key 后写覆盖）。正文、PPT、讲稿必须引用账本中的数字，质检员据此比对跨材料一致性。' +
      '固定 key 集合：revenue_y1..y5, net_profit_y1..y5, gross_margin_y1..y5, breakeven_month, min_cash, financing_amount, pre_money, post_money, equity_offered, cac, ltv, ltv_cac, payback_months, tam, sam, som_y1, sample_size, budget_total, kpi_*；其他数字可自定义 key（小写下划线）。' +
      '每条：key、label 中文名、value 数值、unit 单位（万元/%/月/份…）、source 来源（文件!表!单元格 或 URL）、synthetic 是否示例数据、usedIn 被引用位置。',
    {
      entries: z
        .array(
          z.object({
            key: z.string().min(1).max(60).regex(/^[a-z0-9_]+$/, 'key 只能用小写字母、数字、下划线'),
            label: z.string().min(1).max(100),
            value: z.union([z.number(), z.string().max(100)]),
            unit: z.string().max(20).optional(),
            source: z.string().max(300).optional(),
            usedIn: z.array(z.string().max(200)).max(50).optional(),
            synthetic: z.boolean().optional(),
            note: z.string().max(300).optional(),
          }),
        )
        .min(1)
        .max(200)
        .describe('要登记的数字条目'),
      unit: z.string().max(20).optional().describe('账本默认展示单位（可选，默认 万元）'),
    },
    async (args) => {
      try {
        const r = await mergeNumbers(ctx.workspaceDir, args.entries as NumberEntry[], args.unit);
        const size = JSON.stringify(r).length;
        ctx.onFileWritten?.(r.file, size);
        return text(`已写入 ${r.file}：新增 ${r.added}，更新 ${r.replaced}，账本共 ${r.total} 条。`);
      } catch (e) {
        return text(`写入数字账本失败：${errorMessage(e)}`, true);
      }
    },
  );

  return createSdkMcpServer({
    name: 'bizgold',
    version: '1.0.0',
    tools: [webSearch, webFetch, updateProgress, recordNumbers],
  });
}
