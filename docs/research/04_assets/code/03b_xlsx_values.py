# -*- coding: utf-8 -*-
"""公式缓存值为空的三种对策，全部实际运行：
A) LibreOffice 重算（命令行 / UNO lo_convert.py）
B) formulas 库纯 Python 计算
C) Python 自算后写『数值版』：复制公式版工作簿，把每个公式单元格替换为计算值（样式不变）
"""
import os
import shutil
import subprocess
import sys
import time

from openpyxl import load_workbook

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
SRC = os.path.join(OUT, "03_财务模型.xlsx")
ENV = dict(os.environ, FONTCONFIG_FILE=os.path.join(HERE, "fonts.conf"))

# ---------- A) LibreOffice 重算 ----------
t0 = time.time()
recalc_dir = os.path.join(OUT, "xlsx_recalc")
subprocess.run(["soffice", "--headless", "--convert-to", "xlsx", "--outdir", recalc_dir, SRC], check=True,
               capture_output=True, env=ENV)
recalc_path = os.path.join(recalc_dir, os.path.basename(SRC))
wb_vals = load_workbook(recalc_path, data_only=True)
print(f"A) LibreOffice 重算 {time.time() - t0:.1f}s  NPV={wb_vals['估值']['C10'].value:.2f}  IRR={wb_vals['估值']['C11'].value:.4f}"
      f"  检验行={[wb_vals['资产负债表'][f'{c}13'].value for c in 'CDEFG']}")

# ---------- B) formulas 库 ----------
try:
    import formulas
    t0 = time.time()
    xl = formulas.ExcelModel().loads(SRC).finish()
    sol = xl.calculate()
    def get(sheet, cell):
        key = f"'[{os.path.basename(SRC)}]{sheet.upper()}'!{cell}"
        return sol[key].value[0][0] if key in sol else None
    print(f"B) formulas {formulas.__version__} 计算 {time.time() - t0:.1f}s  NPV={get('估值','C10')}  IRR={get('估值','C11')}"
          f"  净利润2026={get('利润表','E14')}  情景选择={get('假设','C23')} 情景系数={get('假设','C24')}  加权NPV={get('估值','C19')}")
except Exception as e:  # noqa: BLE001
    print("B) formulas 库失败:", type(e).__name__, str(e)[:300])

# ---------- C) 数值版 ----------
# 做法：以公式版为模板，用 LibreOffice(或 formulas/自算) 得到的值覆盖公式单元格 → 样式/图表/命名区域全部保留，
# 适合给评委看的“只读数值版”，也适合 pandas 读取。
t0 = time.time()
wb_f = load_workbook(SRC)              # 公式版（含样式）
n = 0
for ws in wb_f.worksheets:
    ws_v = wb_vals[ws.title]
    for row in ws.iter_rows():
        for c in row:
            if isinstance(c.value, str) and c.value.startswith("="):
                c.value = ws_v[c.coordinate].value
                n += 1
wb_f.calculation.fullCalcOnLoad = False
dst = os.path.join(OUT, "03_财务模型_数值版.xlsx")
wb_f.save(dst)
chk = load_workbook(dst, data_only=True)
print(f"C) 数值版写入 {n} 个单元格 {time.time() - t0:.1f}s  NPV={chk['估值']['C10'].value:.2f}  (pandas 可直接读)")
import pandas as pd
df = pd.read_excel(dst, sheet_name="利润表", header=3, index_col=0).dropna(how="all")
print(df.iloc[:, 1:].head(6).round(1).to_string())
