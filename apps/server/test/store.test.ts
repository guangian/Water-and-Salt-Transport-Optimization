import path from 'node:path';
import os from 'node:os';
import fsp from 'node:fs/promises';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { PROJECT_ID_RE, ProjectStore, buildProject, newProjectId, summarize } from '../src/services/store.js';
import { registerSecret } from '../src/util/log.js';
import type { AgentEvent } from '../src/types.js';

const ws = path.join(os.tmpdir(), `bizgold-store-${process.pid}`);

describe('ProjectStore', () => {
  beforeAll(async () => {
    await fsp.rm(ws, { recursive: true, force: true });
  });
  afterAll(async () => {
    await fsp.rm(ws, { recursive: true, force: true });
  });

  it('generates ids in the p_YYYYMMDD_HHmmss_xxxx format', () => {
    const id = newProjectId(new Date('2026-03-05T14:07:09Z'));
    expect(id).toMatch(PROJECT_ID_RE);
    expect(id.startsWith('p_2026')).toBe(true);
    expect(() => new ProjectStore(ws).dirOf('../etc')).toThrow('项目 id 格式不正确');
  });

  it('buildProject validates input and fills contract defaults', () => {
    expect(() => buildProject({ name: '', teamName: 't', track: 'brand' }, { reviewRounds: 3, targetScore: 85, models: {} })).toThrow();
    expect(() => buildProject({ name: 'x', teamName: 't', track: 'nope' as 'brand' }, { reviewRounds: 3, targetScore: 85, models: {} })).toThrow();
    const p = buildProject(
      { name: '弘宝汝瓷品牌策划', teamName: '瓷韵队', track: 'brand', subject: { name: '弘宝汝瓷' } as never, params: { pageTarget: 60 } as never },
      { reviewRounds: 3, targetScore: 85, models: { main: 'deepseek-v4-pro' } },
    );
    expect(p.status).toBe('idle');
    expect(p.stages.parse.status).toBe('pending');
    expect(p.stages.review.round).toBe(0);
    expect(p.params.reviewRounds).toBe(3);
    expect(p.params.pageTarget).toBe(60);
    expect(p.params.models.main).toBe('deepseek-v4-pro');
    expect(p.subject.name).toBe('弘宝汝瓷');
    expect(p.subject.type).toBe('brand');
    expect(summarize(p)).toMatchObject({ id: p.id, name: p.name, status: 'idle', latestScore: null, stages: expect.objectContaining({ parse: 'pending' }) });
  });

  it('create/get/list/mutate/remove round-trip through project.json', async () => {
    const store = new ProjectStore(ws);
    const p = buildProject({ name: 'A', teamName: 'T', track: 'startup' }, { reviewRounds: 2, targetScore: 80, models: {} });
    await store.create(p);
    expect(await store.exists(p.id)).toBe(true);
    expect((await fsp.readdir(path.join(ws, p.id))).sort()).toEqual(['events.jsonl', 'inputs', 'project.json']);
    const got = await store.get(p.id);
    expect(got).toEqual(p);
    expect((await store.list()).map((s) => s.id)).toEqual([p.id]);

    // 并发 mutate 串行化：两次自增都生效
    await Promise.all([store.mutate(p.id, (x) => void (x.usage.turns += 1)), store.mutate(p.id, (x) => void (x.usage.turns += 1))]);
    expect((await store.get(p.id)).usage.turns).toBe(2);
    await store.remove(p.id);
    expect(await store.exists(p.id)).toBe(false);
    await expect(store.get(p.id)).rejects.toMatchObject({ statusCode: 404 });
  });

  it('appendEvent assigns increasing seq/ts, redacts secrets, broadcasts, and replays since N', async () => {
    const store = new ProjectStore(ws);
    const p = buildProject({ name: 'B', teamName: 'T', track: 'brand' }, { reviewRounds: 3, targetScore: 85, models: {} });
    await store.create(p);
    const secret = 'sk-live-ABCDEFGHIJKLMNOP';
    registerSecret(secret);
    const received: AgentEvent[] = [];
    const unsub = store.subscribe(p.id, (ev) => received.push(ev));

    const all = await Promise.all(
      [1, 2, 3, 4, 5].map((i) => store.appendEvent(p.id, { type: 'text', agent: { id: 'main', name: '主代理' }, data: { text: `第${i}条 ${i === 3 ? secret : ''}` } })),
    );
    expect(all.map((e) => e.seq)).toEqual([1, 2, 3, 4, 5]);
    expect(all.every((e) => e.projectId === p.id && typeof e.ts === 'string')).toBe(true);
    expect(received.map((e) => e.seq)).toEqual([1, 2, 3, 4, 5]);
    unsub();

    const raw = await fsp.readFile(path.join(ws, p.id, 'events.jsonl'), 'utf8');
    expect(raw.trim().split('\n')).toHaveLength(5);
    expect(raw).not.toContain(secret);

    const since3 = await store.readEvents(p.id, 3);
    expect(since3.map((e) => e.seq)).toEqual([4, 5]);

    // persist=false：分配 seq 并广播，但不落盘（增量事件）
    const live: AgentEvent[] = [];
    const unsub2 = store.subscribe(p.id, (ev) => live.push(ev));
    const delta = await store.appendEvent(p.id, { type: 'text_delta', agent: { id: 'main', name: '主代理' }, data: { text: '增' } }, { persist: false });
    unsub2();
    expect(delta.seq).toBe(6);
    expect(live.map((e) => e.seq)).toEqual([6]);
    expect((await store.readEvents(p.id)).map((e) => e.seq)).toEqual([1, 2, 3, 4, 5]);
    const after = await store.appendEvent(p.id, { type: 'text', agent: { id: 'main', name: '主代理' }, data: { text: '完整' } });
    expect(after.seq).toBe(7);

    // 新实例从文件尾部恢复 seq，继续递增
    const fresh = new ProjectStore(ws);
    const next = await fresh.appendEvent(p.id, { type: 'done', agent: { id: 'main', name: '主代理' }, data: {} });
    expect(next.seq).toBe(8);
    expect(received).toHaveLength(5);
  });
});
