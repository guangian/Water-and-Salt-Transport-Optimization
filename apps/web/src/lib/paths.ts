/** 以 basePath 所在目录为基准解析相对路径（`../04_图表/x.png` → `04_图表/x.png`）。 */
export function resolveRelative(basePath: string | undefined, rel: string): string {
  if (/^(data:|https?:|blob:)/.test(rel)) return rel;
  const dir = basePath ? basePath.split('/').slice(0, -1) : [];
  for (const part of rel.split('/')) {
    if (part === '..') dir.pop();
    else if (part && part !== '.') dir.push(part);
  }
  return dir.join('/');
}

export function parentDir(path: string): string {
  const i = path.lastIndexOf('/');
  return i >= 0 ? path.slice(0, i) : '';
}
