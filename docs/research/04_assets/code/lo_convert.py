#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""用 LibreOffice UNO 在服务端完成：更新目录/域 → 另存 docx → 导出 PDF；xlsx 强制重算后另存。

必须用系统 python3（apt install python3-uno），不是 venv 的 python。
用法：
  python3 lo_convert.py in.docx  out_dir   # 生成 out_dir/in.updated.docx + out_dir/in.pdf
  python3 lo_convert.py in.xlsx  out_dir   # 重算公式后生成 out_dir/in.recalc.xlsx + out_dir/in.pdf
  python3 lo_convert.py in.pptx  out_dir   # 生成 out_dir/in.pdf
"""
import os
import subprocess
import sys
import time

import uno  # noqa: E402  (python3-uno)
from com.sun.star.beans import PropertyValue  # noqa: E402

PROFILE = "file:///tmp/officetest/lo_profile"  # 独立用户配置，避免与桌面版冲突、可并发多实例
PIPE = "officetest_pipe"


def _prop(name, value):
    p = PropertyValue()
    p.Name, p.Value = name, value
    return p


def start_office():
    cmd = ["soffice", "--headless", "--invisible", "--nologo", "--norestore",
           f"-env:UserInstallation={PROFILE}", f"--accept=pipe,name={PIPE};urp;StarOffice.ComponentContext"]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    local_ctx = uno.getComponentContext()
    resolver = local_ctx.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local_ctx)
    for _ in range(60):
        try:
            ctx = resolver.resolve(f"uno:pipe,name={PIPE};urp;StarOffice.ComponentContext")
            return proc, ctx
        except Exception:
            time.sleep(0.5)
    proc.kill()
    raise RuntimeError("LibreOffice 启动超时")


def convert(path, out_dir):
    path = os.path.abspath(path)
    out_dir = os.path.abspath(out_dir)  # systemPathToFileUrl 需要绝对路径，否则 IOException
    os.makedirs(out_dir, exist_ok=True)
    stem, ext = os.path.splitext(os.path.basename(path))
    ext = ext.lower()
    proc, ctx = start_office()
    try:
        desktop = ctx.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", ctx)
        doc = desktop.loadComponentFromURL(uno.systemPathToFileUrl(path), "_blank", 0,
                                           (_prop("Hidden", True), _prop("ReadOnly", False)))
        pdf_filter = {"docx": "writer_pdf_Export", "xlsx": "calc_pdf_Export", "pptx": "impress_pdf_Export"}[ext[1:]]
        if ext == ".docx":
            idx = doc.getDocumentIndexes()
            for i in range(idx.getCount()):
                idx.getByIndex(i).update()          # 目录/图表目录
            doc.getTextFields().refresh()           # PAGE/NUMPAGES/SEQ 等域
            # 目录更新后页码可能变化，再更新一次
            for i in range(idx.getCount()):
                idx.getByIndex(i).update()
            out_docx = os.path.join(out_dir, f"{stem}.updated.docx")
            doc.storeToURL(uno.systemPathToFileUrl(out_docx), (_prop("FilterName", "MS Word 2007 XML"),))
            print("docx (TOC updated):", out_docx, "| indexes:", idx.getCount())
        elif ext == ".xlsx":
            doc.calculateAll()
            out_xlsx = os.path.join(out_dir, f"{stem}.recalc.xlsx")
            doc.storeToURL(uno.systemPathToFileUrl(out_xlsx), (_prop("FilterName", "Calc MS Excel 2007 XML"),))
            print("xlsx (recalculated):", out_xlsx)
        out_pdf = os.path.join(out_dir, f"{stem}.pdf")
        doc.storeToURL(uno.systemPathToFileUrl(out_pdf), (_prop("FilterName", pdf_filter),))
        print("pdf:", out_pdf)
        doc.close(True)
    finally:
        try:
            desktop.terminate()
        except Exception:
            pass
        proc.wait(timeout=20)


if __name__ == "__main__":
    convert(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(sys.argv[1])))
