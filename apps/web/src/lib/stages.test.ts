import { describe, expect, it } from 'vitest';
import { currentStageId, emptyStages, parseChecklist, stageDurationMs, stageProgressPercent } from './stages';

describe('stages', () => {
  it('emptyStages 生成 9 个 pending 阶段', () => {
    const s = emptyStages();
    expect(Object.keys(s)).toHaveLength(9);
    expect(Object.values(s).every((x) => x.status === 'pending')).toBe(true);
  });

  it('stageProgressPercent：done 计 1、running 计 0.5', () => {
    const s = emptyStages();
    s.parse.status = 'done';
    s.research.status = 'done';
    s.strategy.status = 'running';
    expect(stageProgressPercent(s)).toBe(Math.round((2.5 / 9) * 100));
  });

  it('currentStageId 取 running，否则取最后 done 的下一阶段', () => {
    const s = emptyStages();
    s.parse.status = 'done';
    s.research.status = 'done';
    expect(currentStageId(s)).toBe('strategy');
    s.strategy.status = 'running';
    expect(currentStageId(s)).toBe('strategy');
    for (const k of Object.keys(s) as (keyof typeof s)[]) s[k].status = 'done';
    expect(currentStageId(s)).toBeNull();
  });

  it('stageDurationMs 计算起止用时，运行中按 now', () => {
    const start = '2026-10-24T10:00:00.000Z';
    expect(stageDurationMs({ status: 'done', startedAt: start, endedAt: '2026-10-24T10:05:00.000Z' })).toBe(5 * 60_000);
    expect(stageDurationMs({ status: 'running', startedAt: start }, Date.parse(start) + 90_000)).toBe(90_000);
    expect(stageDurationMs({ status: 'pending' })).toBeNull();
  });

  it('parseChecklist 解析 ✅/🔄/⏳ 清单', () => {
    const text = ['1. ✅ 读题解析', '2. ✅ 调研', '3. 🔄 策略', '4. ⏳ 财务', '8. ⏳ 评审（第 0/3 轮）', '9. ❌ 打包'].join('\n');
    const parsed = parseChecklist(text);
    expect(parsed.parse).toBe('done');
    expect(parsed.strategy).toBe('running');
    expect(parsed.finance).toBe('pending');
    expect(parsed.review).toBe('pending');
    expect(parsed.package).toBe('failed');
  });
});
