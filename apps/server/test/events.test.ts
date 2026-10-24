import { describe, expect, it } from 'vitest';
import type { SDKMessage } from '@anthropic-ai/claude-agent-sdk';
import { EventNormalizer, compactToolInput, toolResultText, truncate } from '../src/agent/events.js';
import type { AgentEventInput } from '../src/types.js';

const WS = '/ws/p_20260101_000000_abcd';

/** 事件形状来自 /tmp/sdktest/smoke_*.log 与 sdk.d.ts；只保留归一化用到的字段 */
function msg(m: Record<string, unknown>): SDKMessage {
  return { uuid: 'u', session_id: 'sess_1', ...m } as unknown as SDKMessage;
}

function setup(opts: { pricing?: Record<string, { input: number; output: number }> } = {}) {
  const events: AgentEventInput[] = [];
  const pending: { fn: () => void; ms: number }[] = [];
  const n = new EventNormalizer((e) => events.push(e), {
    workspaceDir: WS,
    pricing: opts.pricing,
    schedule: (fn, ms) => pending.push({ fn, ms }),
    fileSize: (abs) => (abs.endsWith('冒烟.md') ? 1234 : undefined),
  });
  const tick = () => {
    while (pending.length) pending.shift()!.fn();
  };
  return { n, events, tick, pending };
}

describe('EventNormalizer', () => {
  it('init → status with session id', () => {
    const { n, events } = setup();
    n.handle(msg({ type: 'system', subtype: 'init', model: 'deepseek-v4-pro', tools: ['Bash', 'Write', 'Task'], mcp_servers: [{ name: 'bizgold', status: 'connected' }], cwd: WS, apiKeySource: 'none' }));
    expect(n.sessionId).toBe('sess_1');
    expect(events).toHaveLength(1);
    expect(events[0]).toMatchObject({ type: 'status', agent: { id: 'main' }, data: { status: 'running', sessionId: 'sess_1', model: 'deepseek-v4-pro' } });
    expect(String(events[0]!.data.message)).toContain('bizgold:connected');
  });

  it('merges text deltas within the window and flushes before a tool_use', () => {
    const { n, events, tick, pending } = setup();
    const delta = (text: string) => msg({ type: 'stream_event', parent_tool_use_id: null, event: { type: 'content_block_delta', index: 0, delta: { type: 'text_delta', text } } });
    n.handle(delta('你'));
    n.handle(delta('好'));
    n.handle(delta('，世界'));
    expect(events).toHaveLength(0);
    // 文本 16ms（契约），思考用更宽的窗口
    expect(pending.map((p) => p.ms)).toEqual([16]);
    n.handle(msg({ type: 'stream_event', parent_tool_use_id: null, event: { type: 'content_block_delta', index: 1, delta: { type: 'thinking_delta', thinking: '嗯' } } }));
    expect(pending.map((p) => p.ms)).toEqual([16, 200]);
    tick();
    expect(events).toEqual([
      { type: 'text_delta', agent: { id: 'main', name: '主代理' }, data: { text: '你好，世界' } },
      { type: 'thinking', agent: { id: 'main', name: '主代理' }, data: { text: '嗯', delta: true } },
    ]);
    events.length = 0;

    n.handle(delta('第二段'));
    n.handle(msg({ type: 'stream_event', parent_tool_use_id: null, event: { type: 'content_block_delta', index: 1, delta: { type: 'thinking_delta', thinking: '思考中' } } }));
    n.handle(msg({ type: 'assistant', parent_tool_use_id: null, message: { role: 'assistant', content: [{ type: 'tool_use', id: 'toolu_1', name: 'Bash', input: { command: 'ls' } }] } }));
    // 缓冲的增量必须在 tool_use 之前被冲出，保证顺序
    expect(events.map((e) => e.type)).toEqual(['text_delta', 'thinking', 'tool_use']);
    expect(events[1]!.data).toEqual({ text: '思考中', delta: true });
  });

  it('Write tool_use + tool_result → tool_use/tool_result/file events with workspace-relative path', () => {
    const { n, events } = setup();
    const content = 'x'.repeat(1000);
    n.handle(
      msg({
        type: 'assistant',
        parent_tool_use_id: null,
        message: { role: 'assistant', content: [{ type: 'tool_use', id: 'toolu_w', name: 'Write', input: { file_path: `${WS}/00_任务解析/冒烟.md`, content } }] },
      }),
    );
    n.handle(
      msg({
        type: 'user',
        parent_tool_use_id: null,
        message: { role: 'user', content: [{ type: 'tool_result', tool_use_id: 'toolu_w', content: 'The file has been updated successfully.' }] },
      }),
    );
    const types = events.map((e) => e.type);
    expect(types).toEqual(['tool_use', 'tool_result', 'file']);
    const use = events[0]!.data as { name: string; input: Record<string, unknown> };
    expect(use.name).toBe('Write');
    // Write 的 content 只保留前 300 字
    expect(String(use.input.content).length).toBeLessThan(400);
    expect(use.input.file_path).toBe(`${WS}/00_任务解析/冒烟.md`);
    expect(events[1]!.data).toMatchObject({ toolUseId: 'toolu_w', name: 'Write', isError: false, output: 'The file has been updated successfully.' });
    expect(events[2]!.data).toEqual({ path: '00_任务解析/冒烟.md', op: 'write', size: 1234 });
  });

  it('does not emit file event when the write failed', () => {
    const { n, events } = setup();
    n.handle(msg({ type: 'assistant', parent_tool_use_id: null, message: { role: 'assistant', content: [{ type: 'tool_use', id: 'toolu_e', name: 'Edit', input: { file_path: 'a.md', old_string: 'a', new_string: 'b' } }] } }));
    n.handle(msg({ type: 'user', parent_tool_use_id: null, message: { role: 'user', content: [{ type: 'tool_result', tool_use_id: 'toolu_e', is_error: true, content: [{ type: 'text', text: 'String not found' }] }] } }));
    expect(events.map((e) => e.type)).toEqual(['tool_use', 'tool_result']);
    expect(events[1]!.data).toMatchObject({ isError: true, output: 'String not found' });
  });

  it('attributes subagent messages via parent_tool_use_id and names them from subagent_type', () => {
    const { n, events } = setup();
    n.handle(
      msg({
        type: 'assistant',
        parent_tool_use_id: null,
        message: { role: 'assistant', content: [{ type: 'tool_use', id: 'toolu_agent', name: 'Agent', input: { subagent_type: 'industry-researcher', description: '调研行业', prompt: '…' } }] },
      }),
    );
    expect(events.map((e) => e.type)).toEqual(['tool_use', 'subagent']);
    expect(events[1]!.data).toMatchObject({ toolUseId: 'toolu_agent', name: '行业研究员', status: 'start' });

    // 子代理内部的文本与工具调用：parent_tool_use_id 指向 Agent 的 tool_use id
    n.handle(msg({ type: 'assistant', parent_tool_use_id: 'toolu_agent', message: { role: 'assistant', content: [{ type: 'text', text: '开始检索' }, { type: 'tool_use', id: 'toolu_sub1', name: 'mcp__bizgold__web_search', input: { query: '米酒 市场规模' } }] } }));
    n.handle(msg({ type: 'user', parent_tool_use_id: 'toolu_agent', message: { role: 'user', content: [{ type: 'tool_result', tool_use_id: 'toolu_sub1', content: '[]' }] } }));
    const sub = events.slice(2);
    expect(sub.map((e) => e.type)).toEqual(['text', 'tool_use', 'tool_result']);
    for (const e of sub) expect(e.agent).toEqual({ id: 'toolu_agent', name: '行业研究员' });

    // 子代理 task_started 重复通知不会再发 start；task_notification 结束它
    n.handle(msg({ type: 'system', subtype: 'task_started', task_id: 't1', tool_use_id: 'toolu_agent', description: '调研行业', subagent_type: 'industry-researcher' }));
    expect(events.filter((e) => e.type === 'subagent')).toHaveLength(1);
    n.handle(msg({ type: 'system', subtype: 'task_notification', task_id: 't1', tool_use_id: 'toolu_agent', status: 'completed', summary: '完成行业调研' }));
    const ends = events.filter((e) => e.type === 'subagent' && e.data.status === 'end');
    expect(ends).toHaveLength(1);
    expect(ends[0]!.data).toMatchObject({ name: '行业研究员', summary: '完成行业调研' });

    // Agent 工具结果到达时不再重复 end
    n.handle(msg({ type: 'user', parent_tool_use_id: null, message: { role: 'user', content: [{ type: 'tool_result', tool_use_id: 'toolu_agent', content: 'done' }] } }));
    expect(events.filter((e) => e.type === 'subagent' && e.data.status === 'end')).toHaveLength(1);
  });

  it('Task tool with unknown type falls back to description', () => {
    const { n, events } = setup();
    n.handle(msg({ type: 'assistant', parent_tool_use_id: null, message: { role: 'assistant', content: [{ type: 'tool_use', id: 'toolu_t', name: 'Task', input: { description: '整理参考资料', prompt: '…' } }] } }));
    expect(events[1]!.data).toMatchObject({ name: '整理参考资料', status: 'start' });
  });

  it('main agent text containing the §6 checklist emits fallback stage events', () => {
    const { n, events } = setup();
    const text = '当前进度：\n1. ✅ 读题解析\n2. 🔄 调研\n3. ⏳ 策略\n8. ⏳ 评审（第 0/3 轮）';
    n.handle(msg({ type: 'assistant', parent_tool_use_id: null, message: { role: 'assistant', content: [{ type: 'text', text }] } }));
    const stages = events.filter((e) => e.type === 'stage');
    expect(stages.map((e) => e.data)).toEqual([
      { stage: 'parse', status: 'done', source: 'fallback' },
      { stage: 'research', status: 'running', source: 'fallback' },
      { stage: 'strategy', status: 'pending', source: 'fallback' },
      { stage: 'review', status: 'pending', source: 'fallback', round: 0 },
    ]);
    // 子代理文本不触发兜底解析
    n.handle(msg({ type: 'assistant', parent_tool_use_id: 'toolu_x', message: { role: 'assistant', content: [{ type: 'text', text }] } }));
    expect(events.filter((e) => e.type === 'stage')).toHaveLength(4);
  });

  it('result → usage (priced from modelUsage) + done', () => {
    const { n, events } = setup({ pricing: { 'deepseek-v4': { input: 1, output: 2 } } });
    n.handle(
      msg({
        type: 'result',
        subtype: 'success',
        is_error: false,
        duration_ms: 13600,
        duration_api_ms: 12000,
        num_turns: 3,
        result: '完成。',
        stop_reason: 'end_turn',
        total_cost_usd: 0.0897,
        usage: { input_tokens: 19997, cache_creation_input_tokens: 0, cache_read_input_tokens: 0, output_tokens: 490 },
        modelUsage: { 'deepseek-v4-pro': { inputTokens: 1_000_000, outputTokens: 500_000, cacheReadInputTokens: 0, cacheCreationInputTokens: 0, webSearchRequests: 0, costUSD: 0, contextWindow: 128000 } },
        permission_denials: [],
      }),
    );
    expect(events.map((e) => e.type)).toEqual(['usage', 'done']);
    expect(events[0]!.data).toMatchObject({ inputTokens: 1_000_000, outputTokens: 500_000, turns: 3, durationMs: 13600, costEstimateUsd: 2 });
    expect(events[1]!.data).toMatchObject({ subtype: 'success', result: '完成。', sdkSubtype: 'success' });
    expect(n.lastResult).toEqual({ subtype: 'success', text: '完成。' });
  });

  it('result without modelUsage falls back to usage; error subtype and abort are reported', () => {
    const { n, events } = setup();
    n.markAborted();
    n.handle(msg({ type: 'result', subtype: 'error_during_execution', is_error: true, duration_ms: 10, duration_api_ms: 5, num_turns: 1, stop_reason: null, total_cost_usd: 0, usage: { input_tokens: 10, output_tokens: 2, cache_read_input_tokens: 3, cache_creation_input_tokens: 1 }, errors: ['aborted'], permission_denials: [] }));
    expect(events[0]!.data).toMatchObject({ inputTokens: 11, outputTokens: 2, cacheReadTokens: 3 });
    expect(events[1]!.data).toMatchObject({ subtype: 'stopped', sdkSubtype: 'error_during_execution', result: 'aborted' });
  });

  it('api_retry / compact_boundary / permission_denied map to status/error', () => {
    const { n, events } = setup();
    n.handle(msg({ type: 'system', subtype: 'api_retry', attempt: 2, max_retries: 10, retry_delay_ms: 4000, error_status: 529 }));
    n.handle(msg({ type: 'system', subtype: 'compact_boundary', compact_metadata: { trigger: 'auto', pre_tokens: 150000 } }));
    n.handle(msg({ type: 'system', subtype: 'permission_denied', tool_name: 'Bash', tool_use_id: 't', input: {}, decision_reason: 'policy' }));
    expect(events.map((e) => e.type)).toEqual(['status', 'status', 'error']);
    expect(String(events[0]!.data.message)).toContain('2/10');
    expect(String(events[1]!.data.message)).toContain('150000');
  });
});

describe('helpers', () => {
  it('truncate appends a note only when exceeding the limit', () => {
    expect(truncate('abc', 5)).toBe('abc');
    expect(truncate('abcdefgh', 5)).toMatch(/^abcde…/);
  });

  it('compactToolInput bounds large non-file inputs', () => {
    const big = compactToolInput('Bash', { command: 'x'.repeat(10000), description: 'd' });
    expect(JSON.stringify(big).length).toBeLessThan(5000);
    expect(big._truncated).toBe(true);
    expect(compactToolInput('Bash', { command: 'ls' })).toEqual({ command: 'ls' });
  });

  it('toolResultText joins text blocks and marks images', () => {
    expect(toolResultText([{ type: 'text', text: 'a' }, { type: 'image' }, { type: 'text', text: 'b' }])).toBe('a\n[图片]\nb');
    expect(toolResultText('plain')).toBe('plain');
    expect(toolResultText(undefined)).toBe('');
  });
});
