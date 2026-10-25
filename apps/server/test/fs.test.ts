import path from 'node:path';
import os from 'node:os';
import fsp from 'node:fs/promises';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { buildTree, copyDir, isInside, safeResolve, sanitizeFileName, toPosixRelative, writeFileAtomic } from '../src/util/fs.js';

const ROOT = path.resolve(os.tmpdir(), 'bizgold-fs-test', 'p_20260101_000000_abcd');

describe('path safety', () => {
  it('isInside accepts root itself and descendants only', () => {
    expect(isInside(ROOT, ROOT)).toBe(true);
    expect(isInside(ROOT, path.join(ROOT, 'a', 'b.md'))).toBe(true);
    expect(isInside(ROOT, path.join(ROOT, '..'))).toBe(false);
    // 前缀相同但不是子目录（…/p_x 与 …/p_x_evil）
    expect(isInside(ROOT, ROOT + '_evil')).toBe(false);
    expect(isInside(ROOT, path.resolve(ROOT, '..', 'other'))).toBe(false);
  });

  it('safeResolve resolves relative paths inside the workspace', () => {
    expect(safeResolve(ROOT, '00_任务解析/冒烟.md')).toBe(path.join(ROOT, '00_任务解析', '冒烟.md'));
    expect(safeResolve(ROOT, '/leading/slash.md')).toBe(path.join(ROOT, 'leading', 'slash.md'));
    expect(safeResolve(ROOT, 'a\\b\\c.txt')).toBe(path.join(ROOT, 'a', 'b', 'c.txt'));
    expect(safeResolve(ROOT, 'a/../b.md')).toBe(path.join(ROOT, 'b.md'));
    expect(safeResolve(ROOT, '')).toBe(ROOT);
  });

  it('safeResolve rejects traversal, absolute paths and NUL bytes', () => {
    expect(() => safeResolve(ROOT, '../../etc/passwd')).toThrow('路径越界');
    expect(() => safeResolve(ROOT, 'a/../../x')).toThrow('路径越界');
    expect(() => safeResolve(ROOT, '..\\..\\windows\\win.ini')).toThrow('路径越界');
    expect(() => safeResolve(ROOT, '/etc/passwd')).not.toThrow(); // 前导 / 被剥掉后落在工作区内
    expect(() => safeResolve(ROOT, 'a\0b')).toThrow('非法路径');
    expect(() => safeResolve(ROOT, '%2e%2e/%2e%2e/etc')).not.toThrow(); // 未解码的 %2e 只是普通文件名
  });

  it('toPosixRelative always uses forward slashes', () => {
    expect(toPosixRelative(ROOT, path.join(ROOT, 'a', 'b', 'c.md'))).toBe('a/b/c.md');
  });

  it('sanitizeFileName strips directories and forbidden characters', () => {
    expect(sanitizeFileName('../../x/评审 报告:v1?.docx')).toBe('评审 报告_v1_.docx');
    expect(sanitizeFileName('..\\..\\evil.exe')).toBe('evil.exe');
    expect(sanitizeFileName('')).toBe('file');
    expect(sanitizeFileName('..', 'fb')).toBe('fb');
  });
});

describe('fs helpers', () => {
  const base = path.join(os.tmpdir(), `bizgold-fs-${process.pid}`);
  beforeAll(async () => {
    await fsp.rm(base, { recursive: true, force: true });
    await fsp.mkdir(path.join(base, 'src', 'keep', '__pycache__'), { recursive: true });
    await fsp.writeFile(path.join(base, 'src', 'a.txt'), 'a');
    await fsp.writeFile(path.join(base, 'src', 'keep', 'b.txt'), 'b');
    await fsp.writeFile(path.join(base, 'src', 'keep', '__pycache__', 'c.pyc'), 'c');
  });
  afterAll(async () => {
    await fsp.rm(base, { recursive: true, force: true });
  });

  it('writeFileAtomic leaves no temp file behind', async () => {
    const f = path.join(base, 'out', 'project.json');
    await writeFileAtomic(f, '{"a":1}');
    await writeFileAtomic(f, '{"a":2}');
    expect(await fsp.readFile(f, 'utf8')).toBe('{"a":2}');
    expect(await fsp.readdir(path.join(base, 'out'))).toEqual(['project.json']);
  });

  it('copyDir skips named entries and reports missing source', async () => {
    expect(await copyDir(path.join(base, 'missing'), path.join(base, 'dst0'))).toBe(false);
    expect(await copyDir(path.join(base, 'src'), path.join(base, 'dst'), { skipNames: ['__pycache__'] })).toBe(true);
    expect(await fsp.readFile(path.join(base, 'dst', 'keep', 'b.txt'), 'utf8')).toBe('b');
    await expect(fsp.access(path.join(base, 'dst', 'keep', '__pycache__'))).rejects.toThrow();
  });

  it('buildTree excludes names, sorts dirs first and uses posix paths', async () => {
    const tree = await buildTree(base, base, new Set(['dst', 'dst0', 'out']));
    expect(tree.map((n) => n.name)).toEqual(['src']);
    const src = tree[0]!;
    expect(src.type).toBe('dir');
    expect(src.children!.map((n) => `${n.type}:${n.path}`)).toEqual(['dir:src/keep', 'file:src/a.txt']);
  });
});
