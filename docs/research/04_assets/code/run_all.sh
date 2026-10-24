#!/usr/bin/env bash
# 一键复现：生成全部示例并转 PDF/PNG。前置：apt 装好 libreoffice-{writer,calc,impress} python3-uno fonts-noto-cjk poppler-utils
set -e
cd "$(dirname "$0")"
export FONTCONFIG_FILE="$PWD/fonts.conf"
PY=./venv/bin/python
echo "== versions =="; $PY --version; $PY -m pip list 2>/dev/null | grep -iE "^(python-docx|python-pptx|openpyxl|matplotlib|pandas|pillow|markdown-it-py|formulas|PyMuPDF|docxtpl|lxml|numpy|fonttools) " ; soffice --version | head -1; pandoc --version | head -1
echo "== 4 matplotlib =="; $PY 04_mpl_demo.py
echo "== 1 docx =="; $PY 01_docx_demo.py; /usr/bin/python3 lo_convert.py "out/01_商业计划书示例.docx" out/pdf
echo "== 2 pptx =="; $PY 02_pptx_demo.py; soffice --headless --convert-to pdf --outdir out/pptx_render "out/02_路演示例.pptx" >/dev/null; pdftoppm -png -r 50 "out/pptx_render/02_路演示例.pdf" out/pptx_render/slide; ls out/pptx_render/slide-*.png | wc -l
echo "== 3 xlsx =="; $PY 03_xlsx_demo.py; $PY 03b_xlsx_values.py 2>&1 | grep -E "^(A|B|C)\)"
echo "== 5 md2docx =="; $PY 05_md2docx_demo.py; /usr/bin/python3 lo_convert.py "out/05_md渲染示例.docx" out/md_render
echo "ALL DONE"
