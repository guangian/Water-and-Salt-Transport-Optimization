"""用户上传资料预读：PDF / docx / xlsx / csv / txt / md / 图片 → UTF-8 文本摘要。

用法：
    python -m bizkit.preread <file> [--out dir] [--max-pages 200] [--max-rows 20]

输出 `<dir>/<文件名>.txt`（默认 dir = 文件所在目录下的 _preread/，文件名保留原扩展名，如 企业资料.pdf.txt），
并打印 JSON：{ok, kind, out, words, pages, tables, sheets?, warnings}。

- PDF：pdfplumber 逐页抽文本，页内表格转为管道表（| a | b |）；
- docx：按文档顺序输出段落（标题按样式加 #）与表格；
- xlsx / csv：每个 sheet 的列名、类型、行数、前 N 行与数值列统计；
- txt / md：原文；
- 图片：只记录尺寸，不做 OCR。
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from ._util import force_utf8_stdout, write_text

_CJK_RE = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff\uf900-\ufaff]")
IMAGE_EXT = {".png", ".jpg", ".jpeg", ".gif", ".webp", ".bmp", ".tif", ".tiff"}


def count_words(text: str) -> int:
    return len(_CJK_RE.findall(text)) + len(re.findall(r"[A-Za-z0-9][A-Za-z0-9.%,]*", text))


def _pipe_table(rows: List[List[Any]]) -> str:
    rows = [[("" if c is None else str(c)).replace("\n", " ").replace("|", "\\|").strip() for c in r] for r in rows if r]
    if not rows:
        return ""
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    out = ["| " + " | ".join(rows[0]) + " |", "|" + "---|" * width]
    out += ["| " + " | ".join(r) + " |" for r in rows[1:]]
    return "\n".join(out)


# ---- 各格式 ----

def read_pdf(path: Path, max_pages: int) -> Dict[str, Any]:
    import pdfplumber
    parts: List[str] = []
    tables = 0
    warnings: List[str] = []
    with pdfplumber.open(str(path)) as pdf:
        total = len(pdf.pages)
        for i, page in enumerate(pdf.pages[:max_pages], 1):
            text = page.extract_text() or ""
            parts.append(f"\n## 第 {i} 页\n{text.strip()}")
            try:
                for tb in page.extract_tables() or []:
                    pt = _pipe_table(tb)
                    if pt:
                        tables += 1
                        parts.append(f"\n[表格 {tables}]\n{pt}")
            except Exception as e:  # 单页表格抽取失败不影响全文
                warnings.append(f"第 {i} 页表格抽取失败：{e}")
        if total > max_pages:
            warnings.append(f"PDF 共 {total} 页，只预读了前 {max_pages} 页")
    body = "\n".join(parts).strip()
    if count_words(body) < 20 and total > 0:
        warnings.append("几乎没有抽到文本，可能是扫描件（需要 OCR）")
    return {"kind": "pdf", "text": body, "pages": total, "tables": tables, "warnings": warnings}


def read_docx(path: Path) -> Dict[str, Any]:
    from docx import Document
    from docx.oxml.ns import qn
    from docx.table import Table
    from docx.text.paragraph import Paragraph
    doc = Document(str(path))
    parts: List[str] = []
    tables = 0
    images = 0
    for child in doc.element.body.iterchildren():
        if child.tag == qn("w:p"):
            p = Paragraph(child, doc)
            if child.findall(".//" + qn("w:drawing")) or child.findall(".//" + qn("w:pict")):
                images += 1
            t = p.text.strip()
            if not t:
                continue
            style = (p.style.name or "") if p.style is not None else ""
            m = re.match(r"(?i)heading\s*(\d)", style)
            if m:
                parts.append("#" * int(m.group(1)) + " " + t)
            elif style.lower().startswith("title"):
                parts.append("# " + t)
            else:
                parts.append(t)
        elif child.tag == qn("w:tbl"):
            tb = Table(child, doc)
            rows = [[c.text for c in r.cells] for r in tb.rows]
            tables += 1
            parts.append(f"\n[表格 {tables}]\n{_pipe_table(rows)}\n")
    text = "\n".join(parts)
    return {"kind": "docx", "text": text, "pages": None, "tables": tables, "images": images, "warnings": []}


def _describe_frame(name: str, df, max_rows: int) -> str:
    import pandas as pd
    lines = [f"\n## 工作表：{name}", f"行数：{len(df)}，列数：{len(df.columns)}", "", "### 列信息", "| 列名 | 类型 | 非空数 | 示例 |", "|---|---|---|---|"]
    for col in df.columns:
        s = df[col]
        dtype = str(s.dtype)
        kind = "数值" if pd.api.types.is_numeric_dtype(s) else ("日期" if pd.api.types.is_datetime64_any_dtype(s) else "文本")
        example = next((str(v) for v in s.dropna().head(3).tolist()), "")
        lines.append(f"| {col} | {kind}（{dtype}） | {int(s.notna().sum())} | {example[:30]} |")
    num = df.select_dtypes("number")
    if not num.empty:
        lines += ["", "### 数值列统计", "| 列名 | 计数 | 均值 | 最小 | 最大 | 合计 |", "|---|---|---|---|---|---|"]
        for col in num.columns:
            s = num[col].dropna()
            if s.empty:
                continue
            lines.append(f"| {col} | {len(s)} | {s.mean():,.2f} | {s.min():,.2f} | {s.max():,.2f} | {s.sum():,.2f} |")
    lines += ["", f"### 前 {min(max_rows, len(df))} 行"]
    head = df.head(max_rows)
    rows = [[str(c) for c in head.columns]] + [[("" if (isinstance(v, float) and v != v) else v) for v in r] for r in head.itertuples(index=False)]
    lines.append(_pipe_table(rows))
    return "\n".join(lines)


def read_table_file(path: Path, max_rows: int) -> Dict[str, Any]:
    import pandas as pd
    warnings: List[str] = []
    sheets: List[Dict[str, Any]] = []
    parts: List[str] = []
    if path.suffix.lower() == ".csv":
        df = None
        for enc in ("utf-8-sig", "gbk", "gb18030", "latin-1"):
            try:
                df = pd.read_csv(path, encoding=enc)
                break
            except UnicodeDecodeError:
                continue
            except Exception as e:
                warnings.append(f"CSV 读取失败：{e}")
                break
        if df is None:
            return {"kind": "csv", "text": "", "pages": None, "tables": 0, "sheets": [], "warnings": warnings or ["无法识别 CSV 编码"]}
        frames = {path.stem: df}
    else:
        frames = pd.read_excel(path, sheet_name=None)
    for name, df in frames.items():
        sheets.append({"name": str(name), "rows": int(len(df)), "cols": int(len(df.columns)), "columns": [str(c) for c in df.columns]})
        parts.append(_describe_frame(str(name), df, max_rows))
    return {"kind": "csv" if path.suffix.lower() == ".csv" else "xlsx", "text": "\n".join(parts), "pages": None,
            "tables": len(sheets), "sheets": sheets, "warnings": warnings}


def read_text_file(path: Path) -> Dict[str, Any]:
    raw = path.read_bytes()
    text = None
    for enc in ("utf-8-sig", "gbk", "gb18030", "latin-1"):
        try:
            text = raw.decode(enc)
            break
        except UnicodeDecodeError:
            continue
    return {"kind": path.suffix.lower().lstrip(".") or "txt", "text": text or "", "pages": None, "tables": 0, "warnings": []}


def read_image(path: Path) -> Dict[str, Any]:
    from PIL import Image
    with Image.open(path) as im:
        w, h = im.size
        mode = im.mode
    text = f"[图片] {path.name}：{w}×{h} 像素，模式 {mode}。未做 OCR，请在需要时用 Read 工具直接查看图片。"
    return {"kind": "image", "text": text, "pages": None, "tables": 0, "width": w, "height": h, "warnings": []}


def preread(path: Path, out_dir: Optional[Path] = None, max_pages: int = 200, max_rows: int = 20) -> Dict[str, Any]:
    path = Path(path).resolve()
    if not path.exists():
        return {"ok": False, "error": f"文件不存在：{path}"}
    out_dir = Path(out_dir).resolve() if out_dir else path.parent / "_preread"
    ext = path.suffix.lower()
    try:
        if ext == ".pdf":
            info = read_pdf(path, max_pages)
        elif ext == ".docx":
            info = read_docx(path)
        elif ext in (".xlsx", ".xlsm", ".xls", ".csv"):
            info = read_table_file(path, max_rows)
        elif ext in IMAGE_EXT:
            info = read_image(path)
        elif ext in (".txt", ".md", ".markdown", ".json", ".log"):
            info = read_text_file(path)
        else:
            return {"ok": False, "error": f"不支持的文件类型：{ext}（支持 pdf/docx/xlsx/csv/txt/md/图片）"}
    except Exception as e:
        return {"ok": False, "error": f"{type(e).__name__}: {e}", "file": str(path)}
    header = f"# 预读摘要：{path.name}\n类型：{info['kind']}" + (f"｜页数：{info['pages']}" if info.get("pages") else "") + (f"｜表格：{info['tables']}" if info.get("tables") else "") + "\n"
    out_path = out_dir / f"{path.name}.txt"
    write_text(out_path, header + "\n" + (info.get("text") or ""))
    result: Dict[str, Any] = {
        "ok": True, "file": str(path), "kind": info["kind"], "out": str(out_path), "words": count_words(info.get("text") or ""),
        "pages": info.get("pages"), "tables": info.get("tables", 0), "warnings": info.get("warnings", []),
    }
    for k in ("sheets", "images", "width", "height"):
        if k in info:
            result[k] = info[k]
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    force_utf8_stdout()
    ap = argparse.ArgumentParser(description="资料预读为文本摘要")
    ap.add_argument("file")
    ap.add_argument("--out", default=None, help="输出目录（默认 <文件目录>/_preread）")
    ap.add_argument("--max-pages", type=int, default=200)
    ap.add_argument("--max-rows", type=int, default=20)
    args = ap.parse_args(argv)
    result = preread(Path(args.file), Path(args.out) if args.out else None, args.max_pages, args.max_rows)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
