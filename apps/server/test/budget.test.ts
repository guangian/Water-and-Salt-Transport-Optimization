import { describe, expect, it } from 'vitest';
import { overBudgetHint, type McpContext } from '../src/agent/mcp.js';
import { STAGE_BUDGET_MIN, STAGE_IDS } from '../src/types.js';

function ctxWith(budget: McpContext['budget']): McpContext {
  return {
    projectId: 'p_test',
    workspaceDir: '/tmp/none',
    search: {} as McpContext['search'],
    python: {} as McpContext['python'],
    updateStage: async () => '',
    budget,
  };
}

describe('stage budget hints', () => {
  it('has a budget for every stage', () => {
    for (const s of STAGE_IDS) expect(STAGE_BUDGET_MIN[s]).toBeGreaterThan(0);
  });

  it('stays silent inside the budget or without a snapshot', async () => {
    expect(await overBudgetHint(ctxWith(undefined))).toBeUndefined();
    expect(await overBudgetHint(ctxWith(async () => null))).toBeUndefined();
    expect(await overBudgetHint(ctxWith(async () => ({ stage: 'research', elapsedMin: 40, budgetMin: 40, totalMin: 50 })))).toBeUndefined();
  });

  it('asks the agent to wrap up once the stage is over budget', async () => {
    const hint = await overBudgetHint(ctxWith(async () => ({ stage: 'research', elapsedMin: 48, budgetMin: 40, totalMin: 62 })));
    expect(hint).toContain('research');
    expect(hint).toContain('48');
    expect(hint).toContain('信息空白');
  });

  it('never throws when the snapshot fails', async () => {
    expect(await overBudgetHint(ctxWith(async () => { throw new Error('boom'); }))).toBeUndefined();
  });
});

describe('python shims', () => {
  it('writes python/python3 shims that forward to the configured interpreter', async () => {
    const { writePythonShims } = await import('../src/agent/workspace.js');
    const os = await import('node:os');
    const fsp = await import('node:fs/promises');
    const path = await import('node:path');
    const ws = await fsp.mkdtemp(path.join(os.tmpdir(), 'bizgold-shim-'));
    const dir = await writePythonShims(ws, process.execPath);
    for (const n of ['python', 'python3', 'python.cmd', 'python3.cmd']) await fsp.access(path.join(dir, n));
    const sh = await fsp.readFile(path.join(dir, 'python3'), 'utf-8');
    expect(sh.startsWith('#!/bin/sh')).toBe(true);
    expect(sh).toContain(process.execPath.replace(/\\/g, '/'));
    if (process.platform !== 'win32') {
      const { execFile } = await import('node:child_process');
      const out = await new Promise<string>((res, rej) => execFile(path.join(dir, 'python3'), ['-e', 'process.stdout.write("shim-ok")'], (e, so) => (e ? rej(e) : res(so))));
      expect(out).toBe('shim-ok');
    }
    await fsp.rm(ws, { recursive: true, force: true });
  });
});
