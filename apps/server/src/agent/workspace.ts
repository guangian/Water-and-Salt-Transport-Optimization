import path from 'node:path';
import fsp from 'node:fs/promises';
import { AGENT_DIR } from '../config.js';
import type { ModelMap, Project, TrackInfo } from '../types.js';
import { copyDir, ensureDir, exists, writeFileAtomic } from '../util/fs.js';
import { createLogger } from '../util/log.js';
import { prereadWordCount, renderTaskSheet } from './prompt.js';

const log = createLogger('workspace');

export const STAGE_DIRS = ['00_任务解析', '01_调研', '02_策略', '03_财务', '04_图表', '05_正文', '06_路演', '07_评审', '08_提交包', '产出', '任务', 'inputs'];

/** agents 模板里的 `model: {{MODEL_XXX}}` 占位符 → 设置里的角色模型 */
const PLACEHOLDER_TO_ROLE: Record<string, keyof ModelMap> = {
  MAIN: 'main',
  RESEARCHER: 'researcher',
  FINANCE: 'finance',
  WRITER: 'writer',
  PITCH: 'pitch',
  JUDGE_A: 'judgeA',
  JUDGE_B: 'judgeB',
  JUDGE_C: 'judgeC',
  INSPECTOR: 'inspector',
  FAST: 'fast',
};

export function replaceModelPlaceholders(text: string, models: Partial<ModelMap>, fallback: string): string {
  return text.replace(/\{\{\s*MODEL_([A-Z_]+)\s*\}\}/g, (m, key: string) => {
    const role = PLACEHOLDER_TO_ROLE[key];
    const model = role ? models[role] : undefined;
    if (!model) log.warn(`unknown/empty model placeholder ${m}, using ${fallback}`);
    return model || fallback;
  });
}

export interface PrepareResult {
  warnings: string[];
  copied: string[];
}

/**
 * 在工作区 .bin 下写 python / python3 垫片，全部转发到配置的解释器（runner 把 .bin 放到 PATH 最前面）。
 * Git Bash / Linux 走无扩展名的 sh 脚本，cmd / PowerShell 走 .cmd；两套都写，谁调用都能命中。
 */
export async function writePythonShims(ws: string, pythonPath: string): Promise<string> {
  const dir = path.join(ws, '.bin');
  await ensureDir(dir);
  const sh = `#!/bin/sh\nexec "${pythonPath.replace(/\\/g, '/')}" "$@"\n`;
  const cmd = `@echo off\r\n"${pythonPath}" %*\r\n`;
  for (const name of ['python', 'python3']) {
    await fsp.writeFile(path.join(dir, name), sh, { encoding: 'utf-8', mode: 0o755 });
    await fsp.chmod(path.join(dir, name), 0o755).catch(() => undefined);
    await fsp.writeFile(path.join(dir, `${name}.cmd`), cmd, 'utf-8');
  }
  return dir;
}

/**
 * 准备工作区：复制 agent/ 大脑（CLAUDE.md、.claude/、toolkit/、knowledge/、tracks/<track>/），
 * 渲染 任务书.md，替换 agents 模型占位符，创建阶段目录。源缺失时只记警告（规范文件由他人并行撰写）。
 */
export async function prepareWorkspace(ws: string, project: Project, models: ModelMap, track: TrackInfo | undefined, agentDir = AGENT_DIR): Promise<PrepareResult> {
  const warnings: string[] = [];
  const copied: string[] = [];
  await ensureDir(ws);
  for (const d of STAGE_DIRS) await ensureDir(path.join(ws, d));
  await ensureDir(path.join(ws, '.claude-home'));
  await ensureDir(path.join(ws, '.cache'));
  // 预先接受“工作区信任”，否则 CLI 会忽略 .claude/settings.json 里的 permissions/env（日志：this workspace has not been trusted）
  const claudeJson = path.join(ws, '.claude-home', '.claude.json');
  if (!(await exists(claudeJson))) {
    await writeFileAtomic(
      claudeJson,
      JSON.stringify({ hasCompletedOnboarding: true, projects: { [ws]: { hasTrustDialogAccepted: true } } }, null, 2),
    );
  }

  const jobs: { src: string; dst: string; skip?: string[] }[] = [
    { src: path.join(agentDir, 'CLAUDE.md'), dst: path.join(ws, 'CLAUDE.md') },
    { src: path.join(agentDir, '.claude'), dst: path.join(ws, '.claude') },
    { src: path.join(agentDir, 'toolkit'), dst: path.join(ws, 'toolkit'), skip: ['__pycache__', 'examples', 'tests', '.pytest_cache'] },
    { src: path.join(agentDir, 'knowledge'), dst: path.join(ws, 'knowledge') },
    { src: path.join(agentDir, 'tracks', project.track), dst: path.join(ws, 'tracks', project.track) },
  ];
  for (const j of jobs) {
    const ok = await copyDir(j.src, j.dst, { skipNames: j.skip ?? ['__pycache__'] });
    if (ok) copied.push(path.relative(agentDir, j.src) || path.basename(j.src));
    else {
      const w = `源不存在，跳过：${path.relative(agentDir, j.src)}`;
      warnings.push(w);
      log.warn(w);
    }
  }

  // 子代理模板：model: {{MODEL_RESEARCHER}} → 设置中的模型
  const agentsDir = path.join(ws, '.claude', 'agents');
  const effective: ModelMap = { ...models, ...(project.params.models as Partial<ModelMap>) } as ModelMap;
  if (await exists(agentsDir)) {
    for (const f of await fsp.readdir(agentsDir)) {
      if (!f.endsWith('.md')) continue;
      const file = path.join(agentsDir, f);
      const text = await fsp.readFile(file, 'utf8');
      if (text.includes('{{')) await fsp.writeFile(file, replaceModelPlaceholders(text, effective, effective.main), 'utf8');
    }
  }

  // 任务书（预读字数用于提示信息量）
  const words: Record<string, number> = {};
  for (const inp of project.inputs) {
    if (!inp.preread) continue;
    const text = await fsp.readFile(path.join(ws, inp.preread), 'utf8').catch(() => '');
    if (text) words[inp.name] = prereadWordCount(text);
  }
  await writeFileAtomic(path.join(ws, '任务书.md'), renderTaskSheet(project, track, words));
  // 工具包与 Python 脚本读 project.json 判断赛道，保持一份最新快照即可（store 已写在同一位置）
  return { warnings, copied };
}
