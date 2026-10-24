import path from 'node:path';
import { spawn } from 'node:child_process';
import { IS_WINDOWS, VERSION } from '../config.js';
import { existsSync } from '../util/fs.js';
import { createLogger, errorMessage } from '../util/log.js';
import { detectGitBash, findInPath, type PythonRunner } from './python.js';

const log = createLogger('health');

export type PdfConverter = 'word' | 'wps' | 'libreoffice' | null;

export interface HealthReport {
  ok: boolean;
  version: string;
  python: { found: boolean; path: string | null; version: string | null; missing: string[]; installing: boolean };
  pdfConverter: PdfConverter;
  fonts: { cjk: string[] };
  gitBash: boolean;
  gitBashPath: string | null;
  workspaceDir: string;
  apiKeySource: 'settings' | 'env' | 'none';
  activeRuns: number;
  checkedAt: string;
}

const CACHE_MS = 20000;

export class HealthService {
  private cache: { ts: number; report: HealthReport } | null = null;

  constructor(
    private readonly python: PythonRunner,
    private readonly ctx: () => Promise<{ workspaceDir: string; apiKeySource: 'settings' | 'env' | 'none'; activeRuns: number }>,
  ) {}

  async report(force = false): Promise<HealthReport> {
    const dyn = await this.ctx();
    if (!force && this.cache && Date.now() - this.cache.ts < CACHE_MS) {
      return { ...this.cache.report, ...dyn, workspaceDir: dyn.workspaceDir };
    }
    const [py, converter, fonts, gitBashPath] = await Promise.all([this.checkPython(), this.detectPdfConverter(), this.detectFonts(), detectGitBash()]);
    const report: HealthReport = {
      ok: py.found && py.missing.length === 0,
      version: VERSION,
      python: py,
      pdfConverter: converter,
      fonts,
      gitBash: IS_WINDOWS ? gitBashPath !== null : true,
      gitBashPath,
      ...dyn,
      checkedAt: new Date().toISOString(),
    };
    this.cache = { ts: Date.now(), report };
    return report;
  }

  private async checkPython(): Promise<HealthReport['python']> {
    const p = await this.python.interpreter();
    if (!p) return { found: false, path: null, version: null, missing: [], installing: false };
    const [version, missing] = await Promise.all([this.python.version(), this.python.missingDeps()]);
    return { found: true, path: p, version, missing, installing: this.python.isInstalling };
  }

  /** 优先用工具包探测（与渲染时一致），不可用时自行找 soffice / WINWORD / WPS */
  private async detectPdfConverter(): Promise<PdfConverter> {
    try {
      if (await this.python.interpreter()) {
        const r = await this.python.module('pdf', ['--detect'], { timeoutMs: 60000 });
        const converters = Array.isArray(r.json?.converters) ? (r.json!.converters as string[]) : [];
        for (const c of converters) {
          if (c === 'word' || c === 'docx2pdf') return 'word';
          if (c === 'wps') return 'wps';
          if (c === 'libreoffice') return 'libreoffice';
        }
        if (r.code === 0 && converters.length === 0) return null;
      }
    } catch (e) {
      log.debug('bizkit.pdf --detect unavailable', errorMessage(e));
    }
    if (IS_WINDOWS) {
      const pf = process.env.ProgramFiles ?? 'C:\\Program Files';
      const pf86 = process.env['ProgramFiles(x86)'] ?? 'C:\\Program Files (x86)';
      const wordPaths = ['root\\Office16', 'Office16', 'root\\Office15', 'Office15'].flatMap((d) => [path.join(pf, 'Microsoft Office', d, 'WINWORD.EXE'), path.join(pf86, 'Microsoft Office', d, 'WINWORD.EXE')]);
      if (wordPaths.some(existsSync) || findInPath('winword')) return 'word';
      const local = process.env.LOCALAPPDATA ?? '';
      if ((local && existsSync(path.join(local, 'Kingsoft', 'WPS Office'))) || findInPath('wps')) return 'wps';
      if ([path.join(pf, 'LibreOffice', 'program', 'soffice.exe'), path.join(pf86, 'LibreOffice', 'program', 'soffice.exe')].some(existsSync) || findInPath('soffice')) return 'libreoffice';
      return null;
    }
    if (findInPath('soffice') || findInPath('libreoffice') || existsSync('/Applications/LibreOffice.app/Contents/MacOS/soffice')) return 'libreoffice';
    return null;
  }

  private async detectFonts(): Promise<{ cjk: string[] }> {
    try {
      if (await this.python.interpreter()) {
        const r = await this.python.module('fonts', [], { timeoutMs: 60000 });
        const list = Array.isArray(r.json?.cjkInstalled) ? (r.json!.cjkInstalled as string[]) : null;
        if (list) return { cjk: list };
      }
    } catch (e) {
      log.debug('bizkit.fonts unavailable', errorMessage(e));
    }
    if (IS_WINDOWS) {
      const fontsDir = path.join(process.env.WINDIR ?? 'C:\\Windows', 'Fonts');
      const known: [string, string][] = [['msyh.ttc', '微软雅黑'], ['simsun.ttc', '宋体'], ['simhei.ttf', '黑体'], ['simkai.ttf', '楷体'], ['simfang.ttf', '仿宋']];
      return { cjk: known.filter(([f]) => existsSync(path.join(fontsDir, f))).map(([, n]) => n) };
    }
    const fcList = findInPath('fc-list');
    if (!fcList) return { cjk: [] };
    const out = await new Promise<string>((resolve) => {
      const child = spawn(fcList, [':lang=zh', 'family'], { windowsHide: true });
      const chunks: Buffer[] = [];
      child.stdout.on('data', (d: Buffer) => chunks.push(d));
      child.on('close', () => resolve(Buffer.concat(chunks).toString('utf8')));
      child.on('error', () => resolve(''));
    });
    const names = new Set<string>();
    for (const line of out.split('\n')) {
      const first = line.split(',')[0]?.trim();
      if (first) names.add(first);
    }
    return { cjk: [...names].sort().slice(0, 40) };
  }
}
