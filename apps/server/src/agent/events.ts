import fs from 'node:fs';
import path from 'node:path';
import type { SDKMessage } from '@anthropic-ai/claude-agent-sdk';
import type { AgentEventInput, AgentRef, ModelPrice } from '../types.js';
import { isInside, toPosixRelative } from '../util/fs.js';
import { estimateCost, type ModelUsageLike } from './pricing.js';
import { parseProgressChecklist } from './progress.js';

export const MAIN_AGENT: AgentRef = { id: 'main', name: '主代理' };

/** 子代理类型 → 中文显示名（与 agent/.claude/agents/*.md 文件名对应） */
export const SUBAGENT_NAMES: Record<string, string> = {
  'industry-researcher': '行业研究员',
  'consumer-insight': '消费者洞察',
  'financial-modeler': '财务建模师',
  'chief-writer': '主笔',
  'pitch-designer': '路演设计师',
  'judge-investor': '评委·投资人',
  'judge-professor': '评委·学者',
  'judge-executive': '评委·企业高管',
  'qa-inspector': '质检员',
  'general-purpose': '通用子代理',
  Explore: '探索子代理',
  Plan: '规划子代理',
};

export const TOOL_INPUT_LIMIT = 4096;
export const TOOL_OUTPUT_LIMIT = 4096;
const FILE_TOOL_PREVIEW = 300;
const WRITE_TOOLS = new Set(['Write', 'Edit', 'MultiEdit', 'NotebookEdit']);
const AGENT_TOOLS = new Set(['Agent', 'Task']);

export interface NormalizerOptions {
  workspaceDir: string;
  pricing?: Record<string, ModelPrice>;
  /** 文本增量合并窗口（契约 §4：16ms） */
  deltaWindowMs?: number;
  /** 思考增量合并窗口（默认 200ms）：思考只作折叠展示，不需要逐字刷新，合并后事件量约为 1/10 */
  thinkingWindowMs?: number;
  /** 定时器注入（测试用） */
  schedule?: (fn: () => void, ms: number) => void;
  /** 文件大小查询（file 事件）；默认 fs.statSync */
  fileSize?: (absPath: string) => number | undefined;
  onSession?: (sessionId: string) => void;
}

interface ToolUseInfo {
  name: string;
  input: Record<string, unknown>;
  agent: AgentRef;
}

export function truncate(s: string, limit: number): string {
  if (s.length <= limit) return s;
  return s.slice(0, limit) + `…（已截断，共 ${s.length} 字符）`;
}

/** 工具输入裁剪：Write/Edit 只留路径与前 300 字；其余整体 ≤4KB */
export function compactToolInput(name: string, input: unknown): Record<string, unknown> {
  const obj: Record<string, unknown> = input && typeof input === 'object' ? { ...(input as Record<string, unknown>) } : { value: input };
  if (WRITE_TOOLS.has(name)) {
    const out: Record<string, unknown> = {};
    for (const key of ['file_path', 'notebook_path', 'path']) if (typeof obj[key] === 'string') out[key] = obj[key];
    for (const key of ['content', 'new_string', 'old_string', 'new_source']) {
      if (typeof obj[key] === 'string') out[key] = truncate(obj[key] as string, FILE_TOOL_PREVIEW);
    }
    if (Array.isArray(obj.edits)) out.edits = `${obj.edits.length} 处修改`;
    return out;
  }
  const json = JSON.stringify(obj);
  if (json.length <= TOOL_INPUT_LIMIT) return obj;
  const out: Record<string, unknown> = {};
  let budget = TOOL_INPUT_LIMIT;
  for (const [k, v] of Object.entries(obj)) {
    const s = typeof v === 'string' ? v : JSON.stringify(v) ?? '';
    const keep = Math.max(0, Math.min(s.length, budget));
    out[k] = typeof v === 'string' ? truncate(v, keep) : keep >= s.length ? v : truncate(s, keep);
    budget -= Math.min(s.length, keep) + k.length + 8;
    if (budget <= 0) {
      out._truncated = true;
      break;
    }
  }
  return out;
}

export function toolResultText(content: unknown): string {
  if (typeof content === 'string') return content;
  if (Array.isArray(content)) {
    return content
      .map((b) => {
        if (b && typeof b === 'object') {
          const r = b as Record<string, unknown>;
          if (r.type === 'text' && typeof r.text === 'string') return r.text;
          if (r.type === 'image') return '[图片]';
        }
        return '';
      })
      .filter(Boolean)
      .join('\n');
  }
  if (content === undefined || content === null) return '';
  return JSON.stringify(content);
}

function subagentName(subagentType: unknown, description: unknown): string {
  if (typeof subagentType === 'string' && subagentType) return SUBAGENT_NAMES[subagentType] ?? subagentType;
  if (typeof description === 'string' && description) return description.slice(0, 30);
  return '子代理';
}

/**
 * SDKMessage → 契约 §4 AgentEvent（不含 seq/ts/projectId，由存储层补）。
 * 子代理归属按 parent_tool_use_id；名称来自 Agent 工具 input.subagent_type / description。
 */
export class EventNormalizer {
  private subagents = new Map<string, AgentRef>();
  private subagentEnded = new Set<string>();
  private toolUses = new Map<string, ToolUseInfo>();
  private deltaBuf = new Map<string, { agent: AgentRef; type: 'text_delta' | 'thinking'; text: string }>();
  private flushScheduled: { text: boolean; thinking: boolean } = { text: false, thinking: false };
  private aborted = false;
  sessionId: string | undefined;
  lastResult: { subtype: string; text: string } | undefined;

  constructor(
    private readonly emit: (ev: AgentEventInput) => void,
    private readonly opts: NormalizerOptions,
  ) {}

  markAborted(): void {
    this.aborted = true;
  }

  private agentFor(parentToolUseId: string | null | undefined): AgentRef {
    if (!parentToolUseId) return MAIN_AGENT;
    return this.subagents.get(parentToolUseId) ?? { id: parentToolUseId, name: '子代理' };
  }

  private registerSubagent(toolUseId: string, name: string): boolean {
    if (this.subagents.has(toolUseId)) return false;
    this.subagents.set(toolUseId, { id: toolUseId, name });
    return true;
  }

  private endSubagent(toolUseId: string, summary: string | undefined): void {
    const ref = this.subagents.get(toolUseId);
    if (!ref || this.subagentEnded.has(toolUseId)) return;
    this.subagentEnded.add(toolUseId);
    this.emit({ type: 'subagent', agent: MAIN_AGENT, data: { toolUseId, name: ref.name, status: 'end', ...(summary ? { summary: truncate(summary, 500) } : {}) } });
  }

  // ---------- 增量合并 ----------

  private pushDelta(agent: AgentRef, type: 'text_delta' | 'thinking', text: string): void {
    if (!text) return;
    const key = `${agent.id}:${type}`;
    const cur = this.deltaBuf.get(key);
    if (cur) cur.text += text;
    else this.deltaBuf.set(key, { agent, type, text });
    const kind = type === 'text_delta' ? 'text' : 'thinking';
    if (!this.flushScheduled[kind]) {
      this.flushScheduled[kind] = true;
      const schedule = this.opts.schedule ?? ((fn, ms) => void setTimeout(fn, ms));
      const ms = kind === 'text' ? (this.opts.deltaWindowMs ?? 16) : (this.opts.thinkingWindowMs ?? 200);
      schedule(() => this.flush(type), ms);
    }
  }

  /** 冲出缓冲的增量；不传 only 则全部冲出（完整块/工具调用到达前必须先冲出，保证顺序） */
  flush(only?: 'text_delta' | 'thinking'): void {
    if (!only || only === 'text_delta') this.flushScheduled.text = false;
    if (!only || only === 'thinking') this.flushScheduled.thinking = false;
    if (this.deltaBuf.size === 0) return;
    for (const [key, it] of [...this.deltaBuf]) {
      if (only && it.type !== only) continue;
      this.deltaBuf.delete(key);
      if (it.type === 'text_delta') this.emit({ type: 'text_delta', agent: it.agent, data: { text: it.text } });
      else this.emit({ type: 'thinking', agent: it.agent, data: { text: it.text, delta: true } });
    }
  }

  // ---------- 主入口 ----------

  handle(msg: SDKMessage): void {
    switch (msg.type) {
      case 'system':
        this.handleSystem(msg);
        return;
      case 'stream_event':
        this.handleStream(msg);
        return;
      case 'assistant':
        this.handleAssistant(msg);
        return;
      case 'user':
        this.handleUser(msg);
        return;
      case 'result':
        this.handleResult(msg);
        return;
      default:
        return; // tool_progress / rate_limit / 未知类型：忽略
    }
  }

  private handleSystem(msg: Extract<SDKMessage, { type: 'system' }>): void {
    switch (msg.subtype) {
      case 'init': {
        this.sessionId = msg.session_id;
        this.opts.onSession?.(msg.session_id);
        const mcp = (msg.mcp_servers ?? []).map((s) => `${s.name}:${s.status}`).join(', ');
        this.emit({
          type: 'status',
          agent: MAIN_AGENT,
          data: { status: 'running', message: `会话已启动：模型 ${msg.model}，工具 ${msg.tools?.length ?? 0} 个${mcp ? `，MCP ${mcp}` : ''}`, sessionId: msg.session_id, model: msg.model },
        });
        return;
      }
      case 'compact_boundary':
        this.emit({
          type: 'status',
          agent: MAIN_AGENT,
          data: { status: 'running', message: `上下文已压缩（${msg.compact_metadata.trigger}，${msg.compact_metadata.pre_tokens} tokens）` },
        });
        return;
      case 'api_retry':
        this.emit({
          type: 'status',
          agent: MAIN_AGENT,
          data: { status: 'running', message: `API 重试 ${msg.attempt}/${msg.max_retries}（${msg.error_status ?? '无响应'}），${Math.round(msg.retry_delay_ms / 1000)}s 后重试` },
        });
        return;
      case 'permission_denied':
        this.emit({ type: 'error', agent: MAIN_AGENT, data: { message: `工具调用被拒绝：${msg.tool_name}${msg.decision_reason ? `（${msg.decision_reason}）` : ''}`, recoverable: true } });
        return;
      case 'task_started': {
        if (!msg.tool_use_id) return;
        const name = subagentName(msg.subagent_type, msg.description);
        if (this.registerSubagent(msg.tool_use_id, name)) {
          this.emit({ type: 'subagent', agent: MAIN_AGENT, data: { toolUseId: msg.tool_use_id, name, status: 'start', description: msg.description } });
        }
        return;
      }
      case 'task_notification':
        if (msg.tool_use_id) this.endSubagent(msg.tool_use_id, msg.summary);
        return;
      default:
        return;
    }
  }

  private handleStream(msg: Extract<SDKMessage, { type: 'stream_event' }>): void {
    const ev = msg.event;
    if (ev.type !== 'content_block_delta') return;
    const agent = this.agentFor(msg.parent_tool_use_id);
    if (ev.delta.type === 'text_delta') this.pushDelta(agent, 'text_delta', ev.delta.text);
    else if (ev.delta.type === 'thinking_delta') this.pushDelta(agent, 'thinking', ev.delta.thinking);
  }

  private handleAssistant(msg: Extract<SDKMessage, { type: 'assistant' }>): void {
    const agent = this.agentFor(msg.parent_tool_use_id);
    const content = Array.isArray(msg.message?.content) ? msg.message.content : [];
    for (const block of content) {
      if (block.type === 'text') {
        this.flush();
        this.emit({ type: 'text', agent, data: { text: block.text } });
        if (agent.id === 'main') {
          const items = parseProgressChecklist(block.text);
          for (const it of items) {
            this.emit({ type: 'stage', agent, data: { stage: it.stage, status: it.status, source: 'fallback', ...(it.round !== undefined ? { round: it.round } : {}) } });
          }
        }
      } else if (block.type === 'thinking') {
        this.flush();
        if (block.thinking) this.emit({ type: 'thinking', agent, data: { text: block.thinking, delta: false } });
      } else if (block.type === 'tool_use') {
        this.flush();
        const input = (block.input ?? {}) as Record<string, unknown>;
        this.toolUses.set(block.id, { name: block.name, input, agent });
        this.emit({ type: 'tool_use', agent, data: { toolUseId: block.id, name: block.name, input: compactToolInput(block.name, input) } });
        if (AGENT_TOOLS.has(block.name)) {
          const name = subagentName(input.subagent_type, input.description);
          if (this.registerSubagent(block.id, name)) {
            this.emit({ type: 'subagent', agent, data: { toolUseId: block.id, name, status: 'start', description: typeof input.description === 'string' ? input.description : undefined } });
          }
        }
      }
    }
    if (msg.error) {
      const e = msg.error as unknown;
      const text = typeof e === 'string' ? e : JSON.stringify(e);
      this.emit({ type: 'error', agent, data: { message: `模型请求错误：${truncate(text, 500)}`, recoverable: true } });
    }
  }

  private handleUser(msg: Extract<SDKMessage, { type: 'user' }>): void {
    const content = msg.message?.content;
    if (!Array.isArray(content)) return;
    const agent = this.agentFor(msg.parent_tool_use_id);
    for (const block of content) {
      if (!block || typeof block !== 'object' || (block as { type?: string }).type !== 'tool_result') continue;
      const tr = block as { tool_use_id: string; content?: unknown; is_error?: boolean };
      const info = this.toolUses.get(tr.tool_use_id);
      const output = toolResultText(tr.content);
      const isError = Boolean(tr.is_error);
      this.flush();
      this.emit({
        type: 'tool_result',
        agent: info?.agent ?? agent,
        data: { toolUseId: tr.tool_use_id, name: info?.name ?? 'unknown', isError, output: truncate(output, TOOL_OUTPUT_LIMIT) },
      });
      if (info && !isError && WRITE_TOOLS.has(info.name)) this.emitFileEvent(info);
      if (this.subagents.has(tr.tool_use_id)) {
        const result = msg.tool_use_result as Record<string, unknown> | undefined;
        const summary = typeof result?.content === 'string' ? result.content : output;
        if (result?.status !== 'async_launched') this.endSubagent(tr.tool_use_id, summary);
      }
    }
  }

  private emitFileEvent(info: ToolUseInfo): void {
    const raw = ['file_path', 'notebook_path', 'path'].map((k) => info.input[k]).find((v): v is string => typeof v === 'string');
    if (!raw) return;
    const abs = path.isAbsolute(raw) ? raw : path.resolve(this.opts.workspaceDir, raw);
    const rel = isInside(this.opts.workspaceDir, abs) ? toPosixRelative(this.opts.workspaceDir, abs) : abs.split(path.sep).join('/');
    const sizeFn = this.opts.fileSize ?? defaultFileSize;
    const size = sizeFn(abs);
    this.emit({ type: 'file', agent: info.agent, data: { path: rel, op: info.name === 'Write' ? 'write' : 'edit', size: size ?? null } });
  }

  private handleResult(msg: Extract<SDKMessage, { type: 'result' }>): void {
    this.flush();
    this.sessionId = msg.session_id;
    const modelUsage = (msg as { modelUsage?: Record<string, ModelUsageLike> }).modelUsage;
    let totals = estimateCost(modelUsage, this.opts.pricing ?? {});
    if (!modelUsage || Object.keys(modelUsage).length === 0) {
      const u = msg.usage ?? ({} as Record<string, number>);
      totals = {
        inputTokens: (u.input_tokens ?? 0) + (u.cache_creation_input_tokens ?? 0),
        outputTokens: u.output_tokens ?? 0,
        cacheReadTokens: u.cache_read_input_tokens ?? 0,
        costEstimateUsd: 0,
      };
    }
    const usage = { ...totals, turns: msg.num_turns, durationMs: msg.duration_ms, sdkCostUsd: msg.total_cost_usd };
    this.emit({ type: 'usage', agent: MAIN_AGENT, data: usage });

    const isSuccess = msg.subtype === 'success' && !msg.is_error;
    const text = msg.subtype === 'success' ? msg.result : (msg.errors ?? []).join('\n');
    const subtype = this.aborted ? 'stopped' : isSuccess ? 'success' : 'error';
    this.lastResult = { subtype, text };
    this.emit({
      type: 'done',
      agent: MAIN_AGENT,
      data: {
        subtype,
        result: truncate(text ?? '', 8000),
        usage,
        sdkSubtype: msg.subtype,
        terminalReason: (msg as { terminal_reason?: string }).terminal_reason,
        sessionId: msg.session_id,
        permissionDenials: msg.permission_denials?.length ?? 0,
      },
    });
  }
}

function defaultFileSize(abs: string): number | undefined {
  try {
    return fs.statSync(abs).size;
  } catch {
    return undefined;
  }
}
