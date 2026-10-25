import { describe, expect, it } from 'vitest';
import { filterItems, initialConsoleState, reduceEvents } from './events';
import type { AgentEvent } from '../types';

const MAIN = { id: 'main', name: '主代理' };
const SUB = { id: 'toolu_01', name: '行业研究员' };

function ev<T extends AgentEvent['type']>(
  seq: number,
  type: T,
  data: Extract<AgentEvent, { type: T }>['data'],
  agent = MAIN,
): AgentEvent {
  return { seq, ts: new Date(1_700_000_000_000 + seq * 1000).toISOString(), projectId: 'p1', type, agent, data } as AgentEvent;
}

describe('reduceEvents', () => {
  it('合并同一代理的连续 text_delta 为一个流式文本项，并在 text 到达时收尾', () => {
    const state = reduceEvents(initialConsoleState(), [
      ev(1, 'text_delta', { text: '汝瓷市场' }),
      ev(2, 'text_delta', { text: '6 年增长 2.3 倍。' }),
      ev(3, 'text', { text: '汝瓷市场 6 年增长 2.3 倍。' }),
    ]);
    const texts = state.items.filter((i) => i.kind === 'text');
    expect(texts).toHaveLength(1);
    expect(texts[0].kind === 'text' && texts[0].text).toBe('汝瓷市场 6 年增长 2.3 倍。');
    expect(texts[0].kind === 'text' && texts[0].streaming).toBe(false);
    expect(state.lastSeq).toBe(3);
  });

  it('不同代理的增量互不干扰', () => {
    const state = reduceEvents(initialConsoleState(), [
      ev(1, 'text_delta', { text: '主代理说' }),
      ev(2, 'text_delta', { text: '研究员说' }, SUB),
      ev(3, 'text_delta', { text: '话' }),
    ]);
    const texts = state.items.filter((i) => i.kind === 'text');
    expect(texts).toHaveLength(2);
    expect(texts[0].kind === 'text' && texts[0].text).toBe('主代理说话');
    expect(texts[1].kind === 'text' && texts[1].text).toBe('研究员说');
    expect(state.agents.map((a) => a.id)).toEqual(['main', 'toolu_01']);
  });

  it('tool_use 与 tool_result 按 toolUseId 配对，并统计泳道调用次数', () => {
    const state = reduceEvents(initialConsoleState(), [
      ev(1, 'tool_use', { toolUseId: 't1', name: 'Read', input: { file_path: '任务书.md' } }),
      ev(2, 'tool_result', { toolUseId: 't1', name: 'Read', isError: false, output: '# 任务书' }),
      ev(3, 'tool_use', { toolUseId: 't2', name: 'web_search', input: { query: '汝瓷 市场规模' } }, SUB),
    ]);
    const tools = state.items.filter((i) => i.kind === 'tool');
    expect(tools).toHaveLength(2);
    expect(tools[0].kind === 'tool' && tools[0].result?.output).toBe('# 任务书');
    expect(state.agents.find((a) => a.id === 'main')?.toolCalls).toBe(1);
    expect(state.agents.find((a) => a.id === 'toolu_01')?.toolCalls).toBe(1);
  });

  it('stage 事件更新阶段状态并记录起止时间', () => {
    const state = reduceEvents(initialConsoleState(), [
      ev(1, 'stage', { stage: 'parse', status: 'running' }),
      ev(2, 'stage', { stage: 'parse', status: 'done', note: '完成' }),
      ev(3, 'stage', { stage: 'research', status: 'running' }),
    ]);
    expect(state.stageUpdates.parse?.status).toBe('done');
    expect(state.stageUpdates.parse?.startedAt).toBeDefined();
    expect(state.stageUpdates.parse?.endedAt).toBeDefined();
    expect(state.currentStage).toBe('research');
  });

  it('subagent start/end 维护泳道状态，done 关闭所有运行中泳道', () => {
    const state = reduceEvents(initialConsoleState(), [
      ev(1, 'subagent', { toolUseId: 'toolu_01', name: '行业研究员', status: 'start' }),
      ev(2, 'text_delta', { text: '市场规模 18.6 亿元。' }, SUB),
      ev(3, 'subagent', { toolUseId: 'toolu_01', name: '行业研究员', status: 'end', summary: '完成市场规模测算' }),
      ev(4, 'usage', { inputTokens: 1000, outputTokens: 200, cacheReadTokens: 50, turns: 3, costEstimateUsd: 0.12 }),
      ev(5, 'done', { subtype: 'success', result: 'ok' }),
    ]);
    const lane = state.agents.find((a) => a.id === 'toolu_01');
    expect(lane?.status).toBe('done');
    expect(lane?.lastText).toContain('完成市场规模测算');
    expect(state.usage.turns).toBe(3);
    expect(state.done?.subtype).toBe('success');
    expect(state.agents.every((a) => a.status !== 'running')).toBe(true);
  });

  it('忽略重复回放的旧事件（断线重连 since=seq）', () => {
    const base = reduceEvents(initialConsoleState(), [ev(1, 'text', { text: 'A' }), ev(2, 'text', { text: 'B' })]);
    const again = reduceEvents(base, [ev(1, 'text', { text: 'A' }), ev(2, 'text', { text: 'B' }), ev(3, 'text', { text: 'C' })]);
    expect(again.items.filter((i) => i.kind === 'text')).toHaveLength(3);
    expect(again.lastSeq).toBe(3);
  });
});

describe('filterItems', () => {
  const state = reduceEvents(initialConsoleState(), [
    ev(1, 'stage', { stage: 'parse', status: 'running' }),
    ev(2, 'text', { text: 'hello' }),
    ev(3, 'tool_use', { toolUseId: 't1', name: 'Read', input: {} }),
    ev(4, 'text', { text: 'sub' }, SUB),
  ]);
  it('只看文本', () => {
    const kinds = filterItems(state.items, 'text', 'all').map((i) => i.kind);
    expect(kinds).toEqual(['stage', 'text', 'text']);
  });
  it('只看工具', () => {
    const kinds = filterItems(state.items, 'tools', 'all').map((i) => i.kind);
    expect(kinds).toEqual(['stage', 'tool']);
  });
  it('按代理过滤', () => {
    const items = filterItems(state.items, 'all', 'toolu_01');
    expect(items.map((i) => i.kind)).toEqual(['stage', 'text']);
  });
});
