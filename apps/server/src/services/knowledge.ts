import path from 'node:path';
import fsp from 'node:fs/promises';
import { AGENT_DIR } from '../config.js';
import type { KnowledgeEntry } from '../types.js';
import { ApiError } from '../types.js';
import { isInside } from '../util/fs.js';

/**
 * 知识库：agent/knowledge/*.md（slug = 文件名去 .md）与 agent/tracks/<id>/spec.md（slug = track-<id>）。
 * slug 只允许 [\w\u4e00-\u9fff.-]，再做前缀校验，避免路径穿越。
 */
export class KnowledgeService {
  constructor(private readonly agentDir = AGENT_DIR) {}

  private get knowledgeDir(): string {
    return path.join(this.agentDir, 'knowledge');
  }
  private get tracksDir(): string {
    return path.join(this.agentDir, 'tracks');
  }

  async list(): Promise<KnowledgeEntry[]> {
    const out: KnowledgeEntry[] = [];
    const files = await fsp.readdir(this.knowledgeDir).catch(() => [] as string[]);
    for (const f of files.sort()) {
      if (!f.toLowerCase().endsWith('.md')) continue;
      const abs = path.join(this.knowledgeDir, f);
      const st = await fsp.stat(abs).catch(() => null);
      if (!st?.isFile()) continue;
      out.push({ slug: f.slice(0, -3), title: await this.titleOf(abs, f.slice(0, -3)), kind: 'knowledge', path: `agent/knowledge/${f}`, size: st.size, mtime: st.mtime.toISOString() });
    }
    const tracks = await fsp.readdir(this.tracksDir).catch(() => [] as string[]);
    for (const t of tracks.sort()) {
      const abs = path.join(this.tracksDir, t, 'spec.md');
      const st = await fsp.stat(abs).catch(() => null);
      if (!st?.isFile()) continue;
      out.push({ slug: `track-${t}`, title: await this.titleOf(abs, `${t} 赛道规范`), kind: 'track', path: `agent/tracks/${t}/spec.md`, size: st.size, mtime: st.mtime.toISOString() });
    }
    return out;
  }

  async read(slug: string): Promise<{ slug: string; title: string; content: string; path: string }> {
    if (!/^[\w\u4e00-\u9fff.\- ·（）()]+$/u.test(slug) || slug.includes('..')) throw new ApiError(400, 'bad_slug', 'slug 不合法');
    const tm = /^track-(.+)$/.exec(slug);
    const abs = tm ? path.join(this.tracksDir, tm[1]!, 'spec.md') : path.join(this.knowledgeDir, `${slug}.md`);
    const root = tm ? this.tracksDir : this.knowledgeDir;
    if (!isInside(root, abs)) throw new ApiError(400, 'bad_slug', 'slug 不合法');
    let content: string;
    try {
      content = await fsp.readFile(abs, 'utf8');
    } catch {
      throw new ApiError(404, 'not_found', '知识条目不存在');
    }
    const rel = path.relative(this.agentDir, abs).split(path.sep).join('/');
    return { slug, title: await this.titleOf(abs, slug, content), content, path: `agent/${rel}` };
  }

  private async titleOf(abs: string, fallback: string, content?: string): Promise<string> {
    const text = content ?? (await fsp.readFile(abs, 'utf8').catch(() => ''));
    const m = /^#\s+(.+)$/m.exec(text.slice(0, 4000));
    return m?.[1]?.trim() || fallback;
  }
}
