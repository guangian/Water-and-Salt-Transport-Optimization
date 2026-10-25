import path from 'node:path';
import os from 'node:os';
import fsp from 'node:fs/promises';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import type { FastifyInstance } from 'fastify';
import { buildApp } from '../src/app.js';
import { createContext } from '../src/context.js';

const base = path.join(os.tmpdir(), `bizgold-app-${process.pid}`);
const settingsDir = path.join(base, 'settings');
const workspaceDir = path.join(base, 'workspace');

let app: FastifyInstance;
let projectId = '';

beforeAll(async () => {
  await fsp.rm(base, { recursive: true, force: true });
  await fsp.mkdir(settingsDir, { recursive: true });
  // 测试环境不读取真实密钥
  delete process.env.BIZGOLD_API_KEY;
  delete process.env.ANTHROPIC_AUTH_TOKEN;
  const ctx = await createContext({ settingsDir, workspaceDir });
  app = await buildApp(ctx, { serveWeb: false });
  await app.ready();
});

afterAll(async () => {
  await app.close();
  await fsp.rm(base, { recursive: true, force: true });
});

describe('HTTP API', () => {
  it('returns the error envelope for unknown /api routes', async () => {
    const res = await app.inject({ method: 'GET', url: '/api/does-not-exist' });
    expect(res.statusCode).toBe(404);
    expect(res.json()).toEqual({ error: { code: 'not_found', message: expect.stringContaining('/api/does-not-exist') } });
  });

  it('GET /api/settings masks the key and reports apiKeySource', async () => {
    const res = await app.inject({ method: 'GET', url: '/api/settings' });
    expect(res.statusCode).toBe(200);
    const body = res.json();
    expect(body.api.apiKeySource).toBe('none');
    expect(body.api.apiKey).toBe('');
    const put = await app.inject({ method: 'PUT', url: '/api/settings', payload: { api: { apiKey: 'sk-unit-test-key-ABCD1234' } } });
    expect(put.json().api.apiKey).toBe('sk-****1234');
    expect(put.json().api.apiKeySource).toBe('settings');
    const bad = await app.inject({ method: 'PUT', url: '/api/settings', payload: '"string"', headers: { 'content-type': 'application/json' } });
    expect(bad.statusCode).toBe(400);
    expect(bad.json().error.code).toBe('bad_request');
  });

  it('GET /api/tracks lists the five tracks with the brand assigned cases', async () => {
    const res = await app.inject({ method: 'GET', url: '/api/tracks' });
    const tracks = res.json() as { id: string; assignedCases: unknown[]; elements: string[] }[];
    expect(tracks.map((t) => t.id)).toEqual(['brand', 'startup', 'digital', 'case', 'mice']);
    // agent/tracks/brand/track.json 存在时以其为准（15 家），否则内置常量（10 家）
    expect(tracks[0]!.assignedCases.length).toBeGreaterThanOrEqual(10);
    expect(tracks[0]!.elements.length).toBeGreaterThanOrEqual(9);
  });

  it('creates, reads, patches and lists a project', async () => {
    const create = await app.inject({
      method: 'POST',
      url: '/api/projects',
      payload: { name: '伟伟冰酒酿品牌策划', teamName: '酿造者', track: 'brand', group: '本科组', subject: { name: '伟伟冰酒酿', company: '蚌埠市伟伟冰米酒餐饮管理有限公司', isAssignedCase: true } },
    });
    expect(create.statusCode).toBe(201);
    const p = create.json();
    projectId = p.id;
    expect(p.id).toMatch(/^p_\d{8}_\d{6}_[a-z0-9]{4}$/);
    expect(p.status).toBe('idle');
    expect(p.params.reviewRounds).toBe(3);
    expect(p.subject.isAssignedCase).toBe(true);

    const invalid = await app.inject({ method: 'POST', url: '/api/projects', payload: { name: 'x', teamName: 'y', track: 'unknown' } });
    expect(invalid.statusCode).toBe(400);

    const get = await app.inject({ method: 'GET', url: `/api/projects/${projectId}` });
    expect(get.json()).toMatchObject({ id: projectId, running: false });

    const patch = await app.inject({ method: 'PATCH', url: `/api/projects/${projectId}`, payload: { name: '新名字', params: { pageTarget: 70 } } });
    expect(patch.json()).toMatchObject({ name: '新名字', params: { pageTarget: 70, reviewRounds: 3 } });

    const list = await app.inject({ method: 'GET', url: '/api/projects' });
    expect(list.json().map((s: { id: string }) => s.id)).toEqual([projectId]);

    const missing = await app.inject({ method: 'GET', url: '/api/projects/p_20000101_000000_zzzz' });
    expect(missing.statusCode).toBe(404);
    const badId = await app.inject({ method: 'GET', url: '/api/projects/..%2F..%2Fetc' });
    expect([400, 404]).toContain(badId.statusCode);
  });

  it('serves files only from inside the workspace', async () => {
    const ws = path.join(workspaceDir, projectId);
    await fsp.mkdir(path.join(ws, '00_任务解析'), { recursive: true });
    await fsp.writeFile(path.join(ws, '00_任务解析', '冒烟.md'), '# 冒烟\n\n![图](../04_图表/a.png)\n');
    await fsp.writeFile(path.join(base, 'outside.txt'), 'secret');

    const ok = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/file?path=${encodeURIComponent('00_任务解析/冒烟.md')}` });
    expect(ok.statusCode).toBe(200);
    expect(ok.headers['content-type']).toContain('text/markdown');
    expect(ok.body).toContain('# 冒烟');

    for (const bad of ['../outside.txt', '../../outside.txt', '..%5Coutside.txt', '/etc/passwd/../../../../outside.txt']) {
      const res = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/file?path=${encodeURIComponent(bad)}` });
      expect(res.statusCode, bad).toBeGreaterThanOrEqual(400);
      expect(res.body, bad).not.toContain('secret');
    }
    const noPath = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/file` });
    expect(noPath.statusCode).toBe(400);
    const dir = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/file?path=inputs` });
    expect(dir.statusCode).toBe(400);
  });

  it('preview rewrites markdown images and classifies json/text', async () => {
    const md = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/preview?path=${encodeURIComponent('00_任务解析/冒烟.md')}` });
    expect(md.json().kind).toBe('markdown');
    expect(md.json().content).toContain(`/api/projects/${projectId}/file?path=04_`);
    const json = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/preview?path=project.json` });
    expect(json.json().kind).toBe('text');
    const traversal = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/preview?path=../outside.txt` });
    expect(traversal.statusCode).toBe(400);
    const nope = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/preview?path=missing.md` });
    expect(nope.statusCode).toBe(404);
  });

  it('refuses to serve or preview a symlink that escapes the workspace', async () => {
    if (process.platform === 'win32') return; // 创建符号链接需要特权
    const ws = path.join(workspaceDir, projectId);
    await fsp.symlink(path.join(base, 'outside.txt'), path.join(ws, 'link.txt'));
    const file = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/file?path=link.txt` });
    expect(file.statusCode).toBe(404);
    expect(file.body).not.toContain('secret');
    const preview = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/preview?path=link.txt` });
    expect(preview.statusCode).toBe(404);
    expect(preview.body).not.toContain('secret');
  });

  it('artifacts tree excludes runtime dirs; review/numbers are empty when missing', async () => {
    const ws = path.join(workspaceDir, projectId);
    await fsp.mkdir(path.join(ws, '.claude-home'), { recursive: true });
    await fsp.mkdir(path.join(ws, 'toolkit'), { recursive: true });
    const tree = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/artifacts` });
    const names = (tree.json() as { name: string }[]).map((n) => n.name);
    expect(names).toContain('00_任务解析');
    expect(names).toContain('project.json');
    expect(names).not.toContain('.claude-home');
    expect(names).not.toContain('toolkit');

    const review = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/review` });
    expect(review.json()).toEqual({ latest: null, rounds: [], machine: null, judges: [] });
    const numbers = await app.inject({ method: 'GET', url: `/api/projects/${projectId}/numbers` });
    expect(numbers.json()).toEqual({ ledger: null, consistency: { checked: 0, consistent: true, issues: [] } });
  });

  it('run without an API key is rejected with 400; chat requires a message', async () => {
    await app.inject({ method: 'PUT', url: '/api/settings', payload: { api: { apiKey: '' } } });
    const run = await app.inject({ method: 'POST', url: `/api/projects/${projectId}/run`, payload: { mode: 'run' } });
    expect(run.statusCode).toBe(400);
    expect(run.json().error.message).toContain('密钥');
    const chat = await app.inject({ method: 'POST', url: `/api/projects/${projectId}/chat`, payload: { message: '  ' } });
    expect(chat.statusCode).toBe(400);
    const stop = await app.inject({ method: 'POST', url: `/api/projects/${projectId}/stop` });
    expect(stop.json()).toEqual({ ok: true, stopped: false });
  });

  it('SSE replays stored events then goes idle', async () => {
    const ctxEvents = await import('../src/services/store.js');
    const store = new ctxEvents.ProjectStore(workspaceDir);
    await store.appendEvent(projectId, { type: 'status', agent: { id: 'main', name: '主代理' }, data: { status: 'idle' } });
    await store.appendEvent(projectId, { type: 'text', agent: { id: 'main', name: '主代理' }, data: { text: 'hello' } });

    const address = await app.listen({ port: 0, host: '127.0.0.1' });
    const ac = new AbortController();
    const res = await fetch(`${address}/api/projects/${projectId}/events?since=1`, { signal: ac.signal, headers: { Origin: 'http://localhost:5173' } });
    expect(res.headers.get('content-type')).toContain('text/event-stream');
    expect(res.headers.get('access-control-allow-origin')).toBe('http://localhost:5173');
    const reader = res.body!.getReader();
    let text = '';
    while (!text.includes('event: idle')) {
      const { value, done } = await reader.read();
      if (done) break;
      text += new TextDecoder().decode(value);
    }
    ac.abort();
    expect(text).toContain(': connected');
    expect(text).not.toContain('"seq":1,');
    expect(text).toContain('id: 2\nevent: agent\ndata: {"seq":2,');
    expect(text).toContain('event: idle');
  });

  it('deletes a project only with ?confirm=1', async () => {
    const noConfirm = await app.inject({ method: 'DELETE', url: `/api/projects/${projectId}` });
    expect(noConfirm.statusCode).toBe(400);
    const del = await app.inject({ method: 'DELETE', url: `/api/projects/${projectId}?confirm=1` });
    expect(del.json()).toEqual({ ok: true });
    await expect(fsp.access(path.join(workspaceDir, projectId))).rejects.toThrow();
  });
});
