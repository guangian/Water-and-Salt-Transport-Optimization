import type { Health, Settings } from '../types';

export const MOCK_MODELS = [
  'deepseek-v4-pro',
  'deepseek-v4-flash',
  'kimi-k3',
  'glm-5.3',
  'qwen3.8-max',
  'qwen3.8-plus',
  'claude-sonnet-4-5',
  'claude-opus-4-1',
  'gpt-5.2',
  'minimax-m3',
];

export const MOCK_SETTINGS: Settings = {
  api: {
    baseUrl: 'https://llm.goaichat.top',
    apiKey: 'sk-****9f3a',
    apiKeySource: 'settings',
    fallbackBaseUrl: 'https://model.kimitk.top',
    models: {
      main: 'deepseek-v4-pro',
      researcher: 'kimi-k3',
      finance: 'deepseek-v4-pro',
      writer: 'deepseek-v4-pro',
      pitch: 'deepseek-v4-pro',
      judgeA: 'kimi-k3',
      judgeB: 'glm-5.3',
      judgeC: 'qwen3.8-max',
      inspector: 'deepseek-v4-pro',
      fast: 'deepseek-v4-flash',
    },
    contextTokens: 200000,
    disableThinking: false,
  },
  search: {
    providers: ['bocha', 'tavily', 'serper', 'bing', 'duckduckgo'],
    bochaKey: 'sk-****c2d1',
    tavilyKey: '',
    serperKey: '',
    maxResults: 8,
  },
  python: { path: '', pipIndex: 'https://pypi.tuna.tsinghua.edu.cn/simple' },
  run: { maxTurns: 1500, totalTimeoutMin: 360, idleTimeoutMin: 45, bashTimeoutMs: 600000, reviewRounds: 3, targetScore: 85 },
  paths: { workspaceDir: '' },
  ui: { theme: 'light', language: 'zh-CN' },
};

export const MOCK_HEALTH: Health = {
  ok: true,
  version: '1.0.0-mock',
  python: {
    found: true,
    path: 'C:\\Users\\bizgold\\anaconda3\\python.exe',
    version: '3.12.7',
    missing: ['pdfplumber'],
  },
  pdfConverter: null,
  fonts: { cjk: ['Microsoft YaHei', 'SimSun', 'SimHei', 'Noto Sans CJK SC'] },
  gitBash: true,
  workspaceDir: 'D:\\bizgold\\workspace',
};
