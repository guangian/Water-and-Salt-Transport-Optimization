import { describe, expect, it } from 'vitest';
import { PROGRESS_LINE_RE, parseProgressChecklist, shouldApplyFallback } from '../src/agent/progress.js';

describe('progress checklist (contract §6)', () => {
  it('parses the canonical 9-line checklist', () => {
    const text = [
      '1. ✅ 读题解析',
      '2. ✅ 调研',
      '3. 🔄 策略',
      '4. ⏳ 财务',
      '5. ⏳ 图表',
      '6. ⏳ 撰写',
      '7. ⏳ 路演',
      '8. ⏳ 评审（第 0/3 轮）',
      '9. ❌ 打包',
    ].join('\n');
    const items = parseProgressChecklist(text);
    expect(items.map((i) => `${i.stage}:${i.status}`)).toEqual([
      'parse:done',
      'research:done',
      'strategy:running',
      'finance:pending',
      'charts:pending',
      'write:pending',
      'pitch:pending',
      'review:pending',
      'package:failed',
    ]);
    expect(items[7]).toEqual({ stage: 'review', status: 'pending', round: 0, totalRounds: 3 });
  });

  it('tolerates markdown bullets, half-width parens, extra spaces and duplicate lines', () => {
    const text = '  3.🔄策略\n8. 🔄 评审 (第 2 / 3 轮)\n8. ✅ 评审\n随便一句话 ✅ 没有编号';
    const items = parseProgressChecklist(text);
    expect(items).toEqual([
      { stage: 'strategy', status: 'running' },
      { stage: 'review', status: 'running', round: 2, totalRounds: 3 },
    ]);
  });

  it('ignores prose that merely mentions stage names', () => {
    expect(parseProgressChecklist('我们现在进入调研阶段，稍后更新 ✅')).toEqual([]);
    expect(parseProgressChecklist('')).toEqual([]);
    expect(PROGRESS_LINE_RE.test('1) ✅ 读题解析')).toBe(false);
  });

  it('fallback only upgrades stage status', () => {
    expect(shouldApplyFallback(undefined, 'running')).toBe(true);
    expect(shouldApplyFallback('pending', 'running')).toBe(true);
    expect(shouldApplyFallback('running', 'done')).toBe(true);
    expect(shouldApplyFallback('done', 'running')).toBe(false);
    expect(shouldApplyFallback('done', 'pending')).toBe(false);
    expect(shouldApplyFallback('running', 'running')).toBe(false);
    expect(shouldApplyFallback('failed', 'done')).toBe(true);
  });
});
