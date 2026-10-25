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
