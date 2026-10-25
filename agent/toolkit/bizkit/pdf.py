"""docx → PDF 转换与 PDF 页数统计。

用法：
    python -m bizkit.pdf in.docx [--out out.pdf] [--prefer libreoffice]
    python -m bizkit.pdf --pages file.pdf
    python -m bizkit.pdf --detect

转换器探测顺序：Windows Word COM（win32com）→ WPS COM → docx2pdf → LibreOffice soffice --headless。
都不可用时返回 {"ok": false, "reason": "..."}，不抛异常（上层只警告）。
`detect()` 返回可用转换器名称列表，供服务端 /health 使用。
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from ._util import force_utf8_stdout, log

WD_FORMAT_PDF = 17


def _is_windows() -> bool:
    return platform.system() == "Windows"


def _com_available(progid: str) -> bool:
    if not _is_windows():
        return False
    try:
        import pythoncom  # type: ignore
        import win32com.client  # type: ignore
    except Exception:
        return False
    try:
        pythoncom.CoInitialize()
        app = win32com.client.DispatchEx(progid)
        try:
            app.Quit()
        except Exception:
            pass
        return True
    except Exception:
        return False
    finally:
        try:
            pythoncom.CoUninitialize()
        except Exception:
            pass


def soffice_path() -> Optional[str]:
    """LibreOffice 可执行文件：PATH → 各平台常见安装位置。"""
    for name in ("soffice", "libreoffice", "soffice.exe"):
        p = shutil.which(name)
        if p:
            return p
    candidates: List[str] = []
    if _is_windows():
        for base in (os.environ.get("ProgramFiles", r"C:\Program Files"), os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")):
            candidates += glob.glob(os.path.join(base, "LibreOffice*", "program", "soffice.exe"))
    elif platform.system() == "Darwin":
        candidates += ["/Applications/LibreOffice.app/Contents/MacOS/soffice"]
    else:
        candidates += ["/usr/bin/soffice", "/usr/lib/libreoffice/program/soffice", "/usr/local/bin/soffice", "/snap/bin/libreoffice"]
        candidates += glob.glob("/opt/libreoffice*/program/soffice")
    for c in candidates:
        if Path(c).exists():
            return c
    return None


@lru_cache(maxsize=1)
def detect() -> List[str]:
    """按优先级返回可用转换器：word / wps / docx2pdf / libreoffice。"""
    found: List[str] = []
    if _com_available("Word.Application"):
        found.append("word")
    if _com_available("KWPS.Application"):
        found.append("wps")
    if _is_windows() or platform.system() == "Darwin":
        try:
            import docx2pdf  # type: ignore  # noqa: F401
            found.append("docx2pdf")
        except Exception:
            pass
    if soffice_path():
        found.append("libreoffice")
    return found


def page_count(pdf_path: Path) -> Optional[int]:
    try:
        from pypdf import PdfReader
        return len(PdfReader(str(pdf_path)).pages)
    except Exception:
        return None


# ---- 各转换器 ----

def _convert_com(progid: str, src: Path, dst: Path) -> None:
    import pythoncom  # type: ignore
    import win32com.client  # type: ignore
    pythoncom.CoInitialize()
    app = None
    doc = None
    try:
        app = win32com.client.DispatchEx(progid)
        app.Visible = False
        try:
            app.DisplayAlerts = 0
        except Exception:
            pass
        doc = app.Documents.Open(str(src), ReadOnly=True, AddToRecentFiles=False)
        # 导出前刷新目录与页码域，否则 PDF 里的目录仍是渲染时的静态占位
        try:
            doc.Fields.Update()
            for toc in doc.TablesOfContents:
                toc.Update()
        except Exception:
            pass
        doc.ExportAsFixedFormat(str(dst), WD_FORMAT_PDF)
    finally:
        try:
            if doc is not None:
                doc.Close(False)
        finally:
            try:
                if app is not None:
                    app.Quit()
            except Exception:
                pass
            pythoncom.CoUninitialize()


def _convert_docx2pdf(src: Path, dst: Path) -> None:
    from docx2pdf import convert as _c  # type: ignore
    _c(str(src), str(dst))


def _convert_libreoffice(src: Path, dst: Path, timeout: int = 240) -> None:
    exe = soffice_path()
    if not exe:
        raise RuntimeError("未找到 soffice")
    # 独立的用户配置目录：避免与正在运行的 LibreOffice 实例争抢锁导致静默失败
    with tempfile.TemporaryDirectory(prefix="bizkit_lo_") as tmp:
        profile = Path(tmp) / "profile"
        outdir = Path(tmp) / "out"
        outdir.mkdir()
        cmd = [exe, "--headless", "--norestore", "--nolockcheck", f"-env:UserInstallation={profile.as_uri()}",
               "--convert-to", "pdf", "--outdir", str(outdir), str(src)]
        env = dict(os.environ)
        env.setdefault("SAL_USE_VCLPLUGIN", "gen")
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, env=env)
        produced = outdir / (src.stem + ".pdf")
        if not produced.exists():
            raise RuntimeError(f"soffice 退出码 {proc.returncode}: {(proc.stderr or proc.stdout).strip()[:400]}")
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(produced), str(dst))


def convert(docx_path: Path | str, out: Optional[Path | str] = None, prefer: Optional[str] = None) -> Dict[str, Any]:
    """docx → pdf。返回 {ok, pdf, converter, pages, tried, reason?}。"""
    src = Path(docx_path).resolve()
    if not src.exists():
        return {"ok": False, "reason": f"文件不存在：{src}", "tried": []}
    dst = Path(out).resolve() if out else src.with_suffix(".pdf")
    order = list(detect())
    if prefer and prefer in order:
        order.remove(prefer)
        order.insert(0, prefer)
    if not order:
        return {"ok": False, "reason": "未找到可用的 PDF 转换器（Word / WPS / docx2pdf / LibreOffice）", "tried": []}
    tried: List[Dict[str, str]] = []
    for name in order:
        t0 = time.time()
        try:
            if name == "word":
                _convert_com("Word.Application", src, dst)
            elif name == "wps":
                _convert_com("KWPS.Application", src, dst)
            elif name == "docx2pdf":
                _convert_docx2pdf(src, dst)
            elif name == "libreoffice":
                _convert_libreoffice(src, dst)
            if dst.exists() and dst.stat().st_size > 0:
                return {"ok": True, "pdf": str(dst), "converter": name, "pages": page_count(dst), "tried": tried,
                        "seconds": round(time.time() - t0, 1)}
            tried.append({"converter": name, "error": "未生成文件"})
        except Exception as e:
            tried.append({"converter": name, "error": f"{type(e).__name__}: {e}"[:400]})
            log(f"{name} 转换失败：{e}")
    return {"ok": False, "reason": "所有转换器均失败", "tried": tried}


def main(argv: Optional[Sequence[str]] = None) -> int:
    force_utf8_stdout()
    ap = argparse.ArgumentParser(description="docx → pdf / PDF 页数")
    ap.add_argument("docx", nargs="?", help="输入 docx")
    ap.add_argument("--out", default=None, help="输出 pdf 路径（默认同名）")
    ap.add_argument("--prefer", default=None, help="优先使用的转换器：word|wps|docx2pdf|libreoffice")
    ap.add_argument("--pages", default=None, metavar="PDF", help="统计 PDF 页数")
    ap.add_argument("--detect", action="store_true", help="列出可用转换器")
    args = ap.parse_args(argv)
    if args.detect:
        print(json.dumps({"ok": True, "converters": detect(), "soffice": soffice_path()}, ensure_ascii=False))
        return 0
    if args.pages:
        n = page_count(Path(args.pages))
        print(json.dumps({"ok": n is not None, "pdf": args.pages, "pages": n}, ensure_ascii=False))
        return 0 if n is not None else 1
    if not args.docx:
        ap.print_help()
        return 1
    result = convert(args.docx, args.out, args.prefer)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("ok") else 1


if __name__ == "__main__":
    sys.exit(main())
