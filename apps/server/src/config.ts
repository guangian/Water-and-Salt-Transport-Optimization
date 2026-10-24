import path from 'node:path';
import os from 'node:os';
import { fileURLToPath } from 'node:url';

/**
 * 路径常量。src/ 与 dist/ 到仓库根的层级相同（apps/server/<src|dist>/config.*），
 * 所以用 import.meta.url 向上三级即可，开发（tsx）与构建产物行为一致。
 */
const here = path.dirname(fileURLToPath(import.meta.url));

export const REPO_ROOT = path.resolve(here, '..', '..', '..');
export const SERVER_ROOT = path.resolve(here, '..');
export const AGENT_DIR = path.join(REPO_ROOT, 'agent');
export const TOOLKIT_DIR = path.join(AGENT_DIR, 'toolkit');
export const WEB_DIST_DIR = path.join(REPO_ROOT, 'apps', 'web', 'dist');
export const DEFAULT_WORKSPACE_DIR = path.join(REPO_ROOT, 'workspace');
export const DEFAULT_SETTINGS_DIR = process.env.BIZGOLD_SETTINGS_DIR || path.join(os.homedir(), '.bizgold');

export const PORT = Number(process.env.BIZGOLD_PORT) || 5178;
/** 本地应用默认只监听回环地址，避免把带密钥的服务暴露到局域网 */
export const HOST = process.env.BIZGOLD_HOST || '127.0.0.1';

export const VERSION = '1.0.0';
export const IS_WINDOWS = process.platform === 'win32';
