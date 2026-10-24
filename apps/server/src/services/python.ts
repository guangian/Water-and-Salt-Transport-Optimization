import path from 'node:path';
import fsp from 'node:fs/promises';
import { spawn } from 'node:child_process';
import { IS_WINDOWS, TOOLKIT_DIR } from '../config.js';
import { existsSync } from '../util/fs.js';
import { createLogger, errorMessage } from '../util/log.js';

const log = createLogger('python');

export interface RunResult {
  code: number | null;
  stdout: string;
  stderr: string;
  /** stdout 中最后一个顶层 JSON 对象（工具包 CLI 约定输出一段 JSON） */
  json: Record<string, unknown> | null;
  timedOut: boolean;
}

export interface RunOptions {
  cwd?: string;
  /** 工作区根：存在 <ws>/toolkit 时优先用它作 PYTHONPATH，否则用 agent/toolkit */
  workspace?: string;
  timeoutMs?: number;
  env?: Record<string, string | undefined>;
  input?: string;
}

/** 在 PATH 中查找可执行文件（不走 shell；Windows 按 PATHEXT 补扩展名） */
export function findInPath(name: string, env: NodeJS.ProcessEnv = process.env): string | null {
  const dirs = (env.PATH ?? env.Path ?? '').split(path.delimiter).filter(Boolean);
  const exts = IS_WINDOWS ? (env.PATHEXT ?? '.EXE;.CMD;.BAT;.COM').split(';').map((e) => e.toLowerCase()) : [''];
  for (const dir of dirs) {
    for (const ext of exts) {
      const candidate = path.join(dir, IS_WINDOWS && !name.toLowerCase().endsWith(ext) ? name + ext : name);
      if (existsSync(candidate)) {
        // Microsoft Store 的 python.exe 占位符会弹商店窗口而不是运行，必须跳过
        if (IS_WINDOWS && /\\WindowsApps\\/i.test(candidate)) continue;
        return candidate;
      }
    }
  }
  return null;
}

export function parseLastJson(stdout: string): Record<string, unknown> | null {
  const text = stdout.trim();
  if (!text) return null;
  // 先试整体，再从最后一个 '{' 往前找可解析的片段（CLI 可能在 JSON 前打印日志）
  try {
    const v = JSON.parse(text);
    if (v && typeof v === 'object' && !Array.isArray(v)) return v as Record<string, unknown>;
  } catch {
    /* fallthrough */
  }
  let idx = text.lastIndexOf('\n{');
  while (idx >= 0) {
    try {
      const v = JSON.parse(text.slice(idx + 1));
      if (v && typeof v === 'object') return v as Record<string, unknown>;
    } catch {
      /* 继续往前 */
    }
    idx = text.lastIndexOf('\n{', idx - 1);
  }
  if (text.startsWith('{')) {
    const end = text.indexOf('\n}');
    if (end > 0) {
      try {
        return JSON.parse(text.slice(0, end + 2)) as Record<string, unknown>;
      } catch {
        /* ignore */
      }
    }
  }
  return null;
}

/**
 * Python 解释器定位与调用：设置 > BIZGOLD_PYTHON > PATH 中的 python > python3。
 * 所有调用不经 shell；统一 UTF-8；PYTHONPATH 指向工具包目录。
 */
export class PythonRunner {
  private resolved: string | null | undefined;
  private depsChecked = false;
  private installing: Promise<void> | null = null;

  constructor(private readonly getConfig: () => Promise<{ path: string; pipIndex: string }>) {}

  invalidate(): void {
    this.resolved = undefined;
    this.depsChecked = false;
  }

  async interpreter(): Promise<string | null> {
    if (this.resolved !== undefined) return this.resolved;
    const cfg = await this.getConfig();
    const candidates: string[] = [];
    if (cfg.path) candidates.push(cfg.path);
    if (process.env.BIZGOLD_PYTHON) candidates.push(process.env.BIZGOLD_PYTHON);
    for (const name of IS_WINDOWS ? ['python', 'python3', 'py'] : ['python3', 'python']) {
      const found = findInPath(name);
      if (found) candidates.push(found);
    }
    for (const c of candidates) {
      const abs = path.isAbsolute(c) ? c : (findInPath(c) ?? c);
      if (IS_WINDOWS && /\\WindowsApps\\/i.test(abs)) continue;
      const ok = await this.probe(abs);
      if (ok) {
        this.resolved = abs;
        log.info('python interpreter', { path: abs, version: ok });
        return abs;
      }
    }
    this.resolved = null;
    log.warn('no usable python interpreter found');
    return null;
  }

  private async probe(exe: string): Promise<string | null> {
    try {
      const r = await this.spawnRaw(exe, ['-c', 'import sys; print(sys.version.split()[0])'], { timeoutMs: 15000 });
      if (r.code === 0 && /^\d+\.\d+/.test(r.stdout.trim())) return r.stdout.trim();
    } catch {
      /* not usable */
    }
    return null;
  }

  async version(): Promise<string | null> {
    const exe = await this.interpreter();
    return exe ? this.probe(exe) : null;
  }

  /** 工具包所在目录：工作区内的 toolkit/ 优先（随项目分发），否则服务级 agent/toolkit */
  toolkitDir(workspace?: string): string {
    if (workspace) {
      const local = path.join(workspace, 'toolkit');
      if (existsSync(path.join(local, 'bizkit'))) return local;
    }
    return TOOLKIT_DIR;
  }

  async run(args: string[], opts: RunOptions = {}): Promise<RunResult> {
    const exe = await this.interpreter();
    if (!exe) throw new Error('未找到可用的 Python 解释器（请在设置中指定 python.path 或设置 BIZGOLD_PYTHON）');
    const toolkit = this.toolkitDir(opts.workspace);
    const pythonPath = [toolkit, process.env.PYTHONPATH].filter(Boolean).join(path.delimiter);
    const r = await this.spawnRaw(exe, args, {
      cwd: opts.cwd ?? opts.workspace ?? toolkit,
      timeoutMs: opts.timeoutMs ?? 10 * 60 * 1000,
      input: opts.input,
      env: { PYTHONPATH: pythonPath, ...(opts.env ?? {}) },
    });
    if (r.code !== 0) log.warn(`python exited ${r.code}: ${args.slice(0, 3).join(' ')}`, r.stderr.slice(-800));
    return r;
  }

  /** 运行工具包模块：python -m bizkit.<module> args */
  async module(mod: string, args: string[], opts: RunOptions = {}): Promise<RunResult> {
    return this.run(['-X', 'utf8', '-m', `bizkit.${mod}`, ...args], opts);
  }

  /** 检查依赖是否可导入；返回缺失包名 */
  async missingDeps(packages = ['docx', 'pptx', 'openpyxl', 'matplotlib', 'pdfplumber', 'pypdf', 'pandas']): Promise<string[]> {
    const exe = await this.interpreter();
    if (!exe) return [...packages];
    const script = [
      'import importlib, json, sys',
      'missing = []',
      `for m in ${JSON.stringify(packages)}:`,
      '    try:',
      '        importlib.import_module(m)',
      '    except Exception:',
      '        missing.append(m)',
      'print(json.dumps(missing))',
    ].join('\n');
    const r = await this.spawnRaw(exe, ['-X', 'utf8', '-c', script], { timeoutMs: 60000 });
    try {
      return JSON.parse(r.stdout.trim()) as string[];
    } catch {
      return [...packages];
    }
  }

  /**
   * 首次调用时检查核心依赖，缺失则后台 pip install -r requirements.txt（不阻塞启动）。
   */
  ensureDeps(): void {
    if (this.depsChecked) return;
    this.depsChecked = true;
    void (async () => {
      const exe = await this.interpreter();
      if (!exe) return;
      const r = await this.spawnRaw(exe, ['-c', 'import docx,pptx,openpyxl,matplotlib'], { timeoutMs: 60000 }).catch(() => null);
      if (r && r.code === 0) {
        log.info('toolkit python deps ok');
        return;
      }
      const req = path.join(TOOLKIT_DIR, 'requirements.txt');
      if (!existsSync(req)) {
        log.warn('requirements.txt not found, skip pip install');
        return;
      }
      const cfg = await this.getConfig();
      log.warn('toolkit deps missing, running pip install in background');
      this.installing = this.spawnRaw(exe, ['-m', 'pip', 'install', '-r', req, ...(cfg.pipIndex ? ['-i', cfg.pipIndex] : [])], {
        timeoutMs: 20 * 60 * 1000,
      })
        .then((res) => {
          if (res.code === 0) log.info('pip install finished');
          else log.error('pip install failed', res.stderr.slice(-1500));
        })
        .catch((e) => log.error('pip install error', errorMessage(e)))
        .finally(() => {
          this.installing = null;
        });
    })();
  }

  get isInstalling(): boolean {
    return this.installing !== null;
  }

  spawnRaw(
    exe: string,
    args: string[],
    opts: { cwd?: string; timeoutMs?: number; env?: Record<string, string | undefined>; input?: string } = {},
  ): Promise<RunResult> {
    return new Promise((resolve, reject) => {
      const env: NodeJS.ProcessEnv = {
        ...process.env,
        PYTHONUTF8: '1',
        PYTHONIOENCODING: 'utf-8',
        PYTHONDONTWRITEBYTECODE: '1',
        ...(opts.env ?? {}),
      };
      let child;
      try {
        child = spawn(exe, args, { cwd: opts.cwd, env, windowsHide: true, stdio: ['pipe', 'pipe', 'pipe'] });
      } catch (e) {
        reject(e);
        return;
      }
      const out: Buffer[] = [];
      const err: Buffer[] = [];
      let timedOut = false;
      const timer = setTimeout(() => {
        timedOut = true;
        child.kill();
      }, opts.timeoutMs ?? 600000);
      child.stdout.on('data', (d: Buffer) => out.push(d));
      child.stderr.on('data', (d: Buffer) => err.push(d));
      child.on('error', (e) => {
        clearTimeout(timer);
        reject(e);
      });
      child.on('close', (code) => {
        clearTimeout(timer);
        const stdout = Buffer.concat(out).toString('utf8');
        const stderr = Buffer.concat(err).toString('utf8');
        resolve({ code, stdout, stderr, json: parseLastJson(stdout), timedOut });
      });
      if (opts.input !== undefined) child.stdin.end(opts.input, 'utf8');
      else child.stdin.end();
    });
  }
}

/** 探测 Windows Git Bash（Claude Code 的 Bash 工具依赖），顺序同技术手册 §4.7 */
export async function detectGitBash(): Promise<string | null> {
  if (!IS_WINDOWS) return null;
  const candidates: string[] = [];
  if (process.env.CLAUDE_CODE_GIT_BASH_PATH) candidates.push(process.env.CLAUDE_CODE_GIT_BASH_PATH);
  const pf = process.env.ProgramFiles ?? 'C:\\Program Files';
  const pf86 = process.env['ProgramFiles(x86)'] ?? 'C:\\Program Files (x86)';
  candidates.push(path.join(pf, 'Git', 'bin', 'bash.exe'), path.join(pf86, 'Git', 'bin', 'bash.exe'));
  const git = findInPath('git');
  if (git) candidates.push(path.resolve(path.dirname(git), '..', 'bin', 'bash.exe'), path.resolve(path.dirname(git), '..', '..', 'bin', 'bash.exe'));
  for (const c of candidates) {
    if (c.toLowerCase().endsWith('bash.exe') && existsSync(c)) return c;
  }
  return null;
}

/** 读取文本文件（容错） */
export async function readTextOrEmpty(file: string): Promise<string> {
  try {
    return await fsp.readFile(file, 'utf8');
  } catch {
    return '';
  }
}
