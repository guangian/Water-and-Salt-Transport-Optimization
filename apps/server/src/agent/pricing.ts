import type { ModelPrice } from '../types.js';

export interface ModelUsageLike {
  inputTokens: number;
  outputTokens: number;
  cacheReadInputTokens?: number;
  cacheCreationInputTokens?: number;
}

export interface UsageTotals {
  inputTokens: number;
  outputTokens: number;
  cacheReadTokens: number;
  costEstimateUsd: number;
}

/**
 * 第三方网关模型在 SDK 里 costBasis='unknown'，total_cost_usd 不可信；
 * 按设置里的价格表（USD / 百万 token）自行估算，没配价格的模型按 0 计。
 * 模型名匹配：精确 → 去掉版本后缀的前缀匹配（如 "deepseek-v4-pro" 匹配 "deepseek-v4"）。
 */
export function estimateCost(modelUsage: Record<string, ModelUsageLike> | undefined, pricing: Record<string, ModelPrice>): UsageTotals {
  const totals: UsageTotals = { inputTokens: 0, outputTokens: 0, cacheReadTokens: 0, costEstimateUsd: 0 };
  if (!modelUsage) return totals;
  for (const [model, u] of Object.entries(modelUsage)) {
    const input = u.inputTokens ?? 0;
    const output = u.outputTokens ?? 0;
    const cacheRead = u.cacheReadInputTokens ?? 0;
    const cacheWrite = u.cacheCreationInputTokens ?? 0;
    totals.inputTokens += input + cacheWrite;
    totals.outputTokens += output;
    totals.cacheReadTokens += cacheRead;
    const price = findPrice(model, pricing);
    if (price) {
      const cacheReadPrice = price.cacheRead ?? price.input * 0.1;
      totals.costEstimateUsd += ((input + cacheWrite) * price.input + output * price.output + cacheRead * cacheReadPrice) / 1e6;
    }
  }
  totals.costEstimateUsd = Math.round(totals.costEstimateUsd * 1e6) / 1e6;
  return totals;
}

export function findPrice(model: string, pricing: Record<string, ModelPrice>): ModelPrice | undefined {
  if (pricing[model]) return pricing[model];
  const lower = model.toLowerCase();
  let best: { key: string; price: ModelPrice } | undefined;
  for (const [key, price] of Object.entries(pricing)) {
    const k = key.toLowerCase();
    if ((lower.startsWith(k) || k.startsWith(lower)) && (!best || k.length > best.key.length)) best = { key: k, price };
  }
  return best?.price;
}
