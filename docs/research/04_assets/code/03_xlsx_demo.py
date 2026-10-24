# -*- coding: utf-8 -*-
"""openpyxl 3.1.5 财务模型演示：生成 out/03_财务模型.xlsx（公式版）与 out/03_财务模型_数值版.xlsx。

颜色约定（投行/咨询建模惯例）：蓝色字=手工输入假设，黑色字=本表公式，绿色字=跨表引用，黄色底=关键假设。
"""
import os

import numpy_financial as npf
from openpyxl import Workbook
from openpyxl.chart import BarChart, LineChart, Reference
from openpyxl.comments import Comment
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side
from openpyxl.utils import FORMULAE, get_column_letter
from openpyxl.workbook.defined_name import DefinedName
from openpyxl.workbook.properties import CalcProperties
from openpyxl.worksheet.datavalidation import DataValidation

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)

YEARS = [2024, 2025, 2026, 2027, 2028]
N = len(YEARS)
COLS = [get_column_letter(3 + i) for i in range(N)]  # C..G 为年份列

# ---------------------------------------------------------------- 样式 ----
F_INPUT = Font(name="Arial", size=10, color="0000FF")          # 蓝：输入
F_FORMULA = Font(name="Arial", size=10, color="000000")        # 黑：公式
F_LINK = Font(name="Arial", size=10, color="008000")           # 绿：跨表引用
F_HEADER = Font(name="微软雅黑", size=10, bold=True, color="FFFFFF")
F_TITLE = Font(name="微软雅黑", size=14, bold=True, color="1F3864")
F_LABEL = Font(name="微软雅黑", size=10)
F_BOLD = Font(name="微软雅黑", size=10, bold=True)
FILL_HEADER = PatternFill("solid", fgColor="1F3864")
FILL_INPUT = PatternFill("solid", fgColor="FFF2CC")
FILL_TOTAL = PatternFill("solid", fgColor="D9E2F3")
THIN = Side(style="thin", color="BFBFBF")
BORDER = Border(top=THIN, bottom=THIN, left=THIN, right=THIN)
TOP_BOTTOM = Border(top=Side(style="thin", color="000000"), bottom=Side(style="double", color="000000"))
NUM = '#,##0;[Red](#,##0);"-"'            # 负数红色加括号，零显示 -
NUM1 = '#,##0.0;[Red](#,##0.0);"-"'
PCT = "0.0%"
MULT = '0.0"x"'
WAN = '#,##0"万元"'
RIGHT = Alignment(horizontal="right", vertical="center")
CENTER = Alignment(horizontal="center", vertical="center")


def title(ws, text, sub=None):
    ws["A1"] = text
    ws["A1"].font = F_TITLE
    if sub:
        ws["A2"] = sub
        ws["A2"].font = Font(name="微软雅黑", size=9, color="7F7F7F")
    ws.sheet_view.showGridLines = False


def header_row(ws, row, labels, start_col=1):
    for j, lab in enumerate(labels):
        c = ws.cell(row=row, column=start_col + j, value=lab)
        c.font = F_HEADER
        c.fill = FILL_HEADER
        c.alignment = CENTER
        c.border = BORDER


def label(ws, row, text, bold=False, indent=0):
    c = ws.cell(row=row, column=1, value=text)
    c.font = F_BOLD if bold else F_LABEL
    c.alignment = Alignment(indent=indent)
    return c


def put(ws, row, col, value, kind="formula", fmt=NUM, bold=False, fill=None):
    """kind: input / formula / link"""
    c = ws.cell(row=row, column=col, value=value)
    c.font = Font(name="Arial", size=10, bold=bold,
                  color={"input": "0000FF", "formula": "000000", "link": "008000"}[kind])
    c.number_format = fmt
    c.alignment = RIGHT
    c.border = BORDER
    if fill:
        c.fill = fill
    return c


def year_row(ws, row, text, cells, kind="formula", fmt=NUM, bold=False, fill=None, indent=0):
    label(ws, row, text, bold, indent)
    for i, v in enumerate(cells):
        put(ws, row, 3 + i, v, kind, fmt, bold, fill)


def set_widths(ws, widths):
    for i, w in enumerate(widths, start=1):
        ws.column_dimensions[get_column_letter(i)].width = w


wb = Workbook()

# ================================================================ 假设页 ====
ws = wb.active
ws.title = "假设"
title(ws, "关键假设（蓝字 = 可修改输入；黄底 = 核心假设）", "单位：万元 / 户 / %")
header_row(ws, 4, ["项目", "单位"] + [f"{y}E" if y > 2024 else f"{y}A" for y in YEARS])
assump = [
    # (标签, 单位, 值列表, 格式, 名称)
    ("付费家庭数", "户", [2500, 8000, 30000, 65000, 110000], "#,##0", None),
    ("客单价（设备）", "元/户", [1200, 1150, 1100, 1050, 1000], "#,##0", None),
    ("年订阅费", "元/户/年", [360, 360, 400, 400, 420], "#,##0", None),
    ("订阅转化率", "%", [0.6, 0.65, 0.7, 0.72, 0.75], PCT, None),
    ("硬件成本率", "%", [0.62, 0.58, 0.55, 0.52, 0.50], PCT, None),
    ("订阅服务成本率", "%", [0.35, 0.32, 0.30, 0.28, 0.26], PCT, None),
    ("销售费用率", "%", [0.25, 0.20, 0.16, 0.14, 0.12], PCT, None),
    ("研发费用", "万元", [150, 180, 220, 260, 300], "#,##0", None),
    ("管理费用", "万元", [80, 100, 130, 160, 200], "#,##0", None),
    ("资本开支", "万元", [200, 150, 250, 300, 300], "#,##0", None),
    ("折旧年限", "年", [5, 5, 5, 5, 5], "0", None),
]
ROW_A = {}
r = 5
for lab, unit, vals, fmt, _ in assump:
    label(ws, r, lab)
    ws.cell(row=r, column=2, value=unit).font = F_LABEL
    for i, v in enumerate(vals):
        put(ws, r, 3 + i, v, "input", fmt, fill=FILL_INPUT if lab in ("付费家庭数", "客单价（设备）") else None)
    ROW_A[lab] = r
    r += 1
r += 1
label(ws, r, "全局参数", bold=True)
r += 1
globals_ = [("所得税率", 0.25, PCT, "税率"), ("折现率（WACC）", 0.12, PCT, "折现率"), ("永续增长率", 0.03, PCT, "永续增长率"),
            ("初始投资（2023 末）", 800, "#,##0", "初始投资"), ("营运资金占收入比", 0.10, PCT, "营运资金比")]
for lab, v, fmt, name in globals_:
    label(ws, r, lab)
    put(ws, r, 3, v, "input", fmt, fill=FILL_INPUT)
    # 命名区域：公式里可写 税率 而不是 假设!$C$20
    wb.defined_names[name] = DefinedName(name, attr_text=f"'假设'!$C${r}")
    ROW_A[lab] = r
    r += 1
label(ws, r, "情景选择")
put(ws, r, 3, "基准", "input", "@", fill=FILL_INPUT)
dv = DataValidation(type="list", formula1='"基准,乐观,悲观"', allow_blank=False, showDropDown=False)
dv.error, dv.errorTitle = "请选择 基准/乐观/悲观", "无效输入"
ws.add_data_validation(dv)
dv.add(ws.cell(row=r, column=3))
ROW_A["情景选择"] = r
r += 1
label(ws, r, "情景收入系数")
put(ws, r, 3, f'=IF(C{r - 1}="乐观",1.2,IF(C{r - 1}="悲观",0.8,1))', "formula", "0.00")
wb.defined_names["情景系数"] = DefinedName("情景系数", attr_text=f"'假设'!$C${r}")
ROW_A["情景系数"] = r
dv2 = DataValidation(type="decimal", operator="between", formula1="0", formula2="1", allow_blank=False)
dv2.error, dv2.errorTitle = "比率应在 0~1 之间", "范围错误"
ws.add_data_validation(dv2)
dv2.add(f"C{ROW_A['订阅转化率']}:G{ROW_A['订阅服务成本率']}")
ws.cell(row=ROW_A["付费家庭数"], column=3).comment = Comment("2024 年实际数，来源：公司后台。", "财务模型")
ws.freeze_panes = "C5"
set_widths(ws, [24, 12, 12, 12, 12, 12, 12])

# ================================================================ 收入页 ====
wr = wb.create_sheet("收入")
title(wr, "收入预测", "单位：万元；绿色 = 引用假设页")
header_row(wr, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wr, 5, "付费家庭数（户）", [f"='假设'!{c}{ROW_A['付费家庭数']}*情景系数" for c in COLS], "link", "#,##0")
year_row(wr, 6, "设备收入", [f"={c}5*'假设'!{c}{ROW_A['客单价（设备）']}/10000" for c in COLS], "link")
year_row(wr, 7, "订阅收入", [f"={c}5*'假设'!{c}{ROW_A['订阅转化率']}*'假设'!{c}{ROW_A['年订阅费']}/10000" for c in COLS], "link")
year_row(wr, 8, "营业收入合计", [f"=SUM({c}6:{c}7)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wr, 9, "收入增速", ["" if i == 0 else f"=IF({COLS[i - 1]}8=0,\"\",{c}8/{COLS[i - 1]}8-1)" for i, c in enumerate(COLS)],
         "formula", PCT)
year_row(wr, 10, "订阅收入占比", [f"=IF({c}8=0,0,{c}7/{c}8)" for c in COLS], "formula", PCT)
wr.freeze_panes = "C5"
set_widths(wr, [24, 4, 12, 12, 12, 12, 12])

# ================================================================ 成本页 ====
wc = wb.create_sheet("成本")
title(wc, "成本与费用", "单位：万元")
header_row(wc, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wc, 5, "硬件成本", [f"='收入'!{c}6*'假设'!{c}{ROW_A['硬件成本率']}" for c in COLS], "link")
year_row(wc, 6, "订阅服务成本", [f"='收入'!{c}7*'假设'!{c}{ROW_A['订阅服务成本率']}" for c in COLS], "link")
year_row(wc, 7, "营业成本合计", [f"=SUM({c}5:{c}6)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wc, 8, "销售费用", [f"='收入'!{c}8*'假设'!{c}{ROW_A['销售费用率']}" for c in COLS], "link")
year_row(wc, 9, "研发费用", [f"='假设'!{c}{ROW_A['研发费用']}" for c in COLS], "link")
year_row(wc, 10, "管理费用", [f"='假设'!{c}{ROW_A['管理费用']}" for c in COLS], "link")
year_row(wc, 11, "期间费用合计", [f"=SUM({c}8:{c}10)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
# 折旧：SUMPRODUCT 对历年资本开支按折旧年限分摊（简化：直线法，当年起计）
dep = []
for i, c in enumerate(COLS):
    rng = f"'假设'!$C${ROW_A['资本开支']}:{c}${ROW_A['资本开支']}"
    life = f"'假设'!$C${ROW_A['折旧年限']}:{c}${ROW_A['折旧年限']}"
    dep.append(f"=SUMPRODUCT({rng}/{life})")
year_row(wc, 12, "折旧摊销（SUMPRODUCT）", dep, "link")
wc.freeze_panes = "C5"
set_widths(wc, [24, 4, 12, 12, 12, 12, 12])

# ================================================================ 利润表 ====
wi = wb.create_sheet("利润表")
title(wi, "利润表", "单位：万元")
header_row(wi, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wi, 5, "营业收入", [f"='收入'!{c}8" for c in COLS], "link", bold=True)
year_row(wi, 6, "减：营业成本", [f"=-'成本'!{c}7" for c in COLS], "link")
year_row(wi, 7, "毛利", [f"={c}5+{c}6" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wi, 8, "毛利率", [f"=IF({c}5=0,0,{c}7/{c}5)" for c in COLS], "formula", PCT)
year_row(wi, 9, "减：期间费用", [f"=-'成本'!{c}11" for c in COLS], "link")
year_row(wi, 10, "EBITDA", [f"={c}7+{c}9" for c in COLS], "formula", bold=True)
year_row(wi, 11, "减：折旧摊销", [f"=-'成本'!{c}12" for c in COLS], "link")
year_row(wi, 12, "EBIT", [f"={c}10+{c}11" for c in COLS], "formula", bold=True)
year_row(wi, 13, "减：所得税（IF 亏损不缴）", [f"=-IF({c}12>0,{c}12*税率,0)" for c in COLS], "formula")
year_row(wi, 14, "净利润", [f"={c}12+{c}13" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wi, 15, "净利率", [f"=IF({c}5=0,0,{c}14/{c}5)" for c in COLS], "formula", PCT)
for c in COLS:
    wi[f"{c}14"].border = TOP_BOTTOM
wi.freeze_panes = "C5"
set_widths(wi, [26, 4, 12, 12, 12, 12, 12])

# ================================================================ 现金流量表 ====
wf = wb.create_sheet("现金流量表")
title(wf, "现金流量表（自由现金流口径）", "单位：万元")
header_row(wf, 4, ["项目", "2023"] + [str(y) for y in YEARS])
year_row(wf, 5, "净利润", [f"='利润表'!{c}14" for c in COLS], "link")
year_row(wf, 6, "加：折旧摊销", [f"='成本'!{c}12" for c in COLS], "link")
year_row(wf, 7, "减：资本开支", [f"=-'假设'!{c}{ROW_A['资本开支']}" for c in COLS], "link")
year_row(wf, 8, "减：营运资金增加", [f"=-('收入'!{c}8-{'0' if i == 0 else f'收入!{COLS[i - 1]}8'})*营运资金比"
                                 if i else f"=-'收入'!{c}8*营运资金比" for i, c in enumerate(COLS)], "link")
year_row(wf, 9, "自由现金流（FCF）", [f"=SUM({c}5:{c}8)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
put(wf, 9, 2, "=-初始投资", "link", NUM, bold=True, fill=FILL_TOTAL)  # B9 = 期初投资
year_row(wf, 10, "累计现金流", [f"=SUM($B$9:{c}9)" for c in COLS], "formula", bold=True)
put(wf, 10, 2, "=B9", "formula", NUM, bold=True)
wf.freeze_panes = "C5"
set_widths(wf, [26, 12, 12, 12, 12, 12, 12])

# ================================================================ 资产负债表 ====
wbsh = wb.create_sheet("资产负债表")
title(wbsh, "资产负债表（简化）", "单位：万元；检验行应为 0")
header_row(wbsh, 4, ["项目", ""] + [str(y) for y in YEARS])
year_row(wbsh, 5, "货币资金", [f"='现金流量表'!{c}10+初始投资" for c in COLS], "link")
year_row(wbsh, 6, "营运资金", [f"='收入'!{c}8*营运资金比" for c in COLS], "link")
year_row(wbsh, 7, "固定资产净值", [f"=SUM('假设'!$C${ROW_A['资本开支']}:{c}${ROW_A['资本开支']})-SUM('成本'!$C$12:{c}$12)" for c in COLS], "link")
year_row(wbsh, 8, "初始投资形成的长期资产", ["=初始投资"] * N, "link")
year_row(wbsh, 9, "资产合计", [f"=SUM({c}5:{c}8)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wbsh, 10, "实收资本", ["=初始投资"] * N, "link")
year_row(wbsh, 11, "留存收益", [f"=SUM('利润表'!$C$14:{c}$14)" for c in COLS], "link")
year_row(wbsh, 12, "所有者权益合计", [f"=SUM({c}10:{c}11)" for c in COLS], "formula", bold=True, fill=FILL_TOTAL)
year_row(wbsh, 13, "检验：资产-权益（应为 0）", [f"=ROUND({c}9-{c}12,2)" for c in COLS], "formula", NUM)
wbsh.freeze_panes = "C5"
set_widths(wbsh, [26, 4, 12, 12, 12, 12, 12])

# ================================================================ 估值页 ====
wv = wb.create_sheet("估值")
title(wv, "DCF 估值", "单位：万元")
label(wv, 4, "折现率（WACC）")
put(wv, 4, 3, "=折现率", "link", PCT)
label(wv, 5, "永续增长率")
put(wv, 5, 3, "=永续增长率", "link", PCT)
label(wv, 6, "预测期 FCF 现值（NPV）")
put(wv, 6, 3, f"=NPV(C4,'现金流量表'!C9:G9)", "link")
label(wv, 7, "终值（Gordon）")
put(wv, 7, 3, "='现金流量表'!G9*(1+C5)/(C4-C5)", "link")
label(wv, 8, "终值现值")
put(wv, 8, 3, f"=C7/(1+C4)^{N}", "formula")
label(wv, 9, "企业价值（EV）", bold=True)
put(wv, 9, 3, "=C6+C8", "formula", bold=True, fill=FILL_TOTAL)
label(wv, 10, "项目 NPV（含初始投资）", bold=True)
put(wv, 10, 3, "=C6-初始投资", "formula", bold=True, fill=FILL_TOTAL)
label(wv, 11, "项目 IRR", bold=True)
put(wv, 11, 3, "=IRR('现金流量表'!B9:G9)", "link", PCT, bold=True, fill=FILL_TOTAL)
label(wv, 12, "投资回收期（年，静态）")
put(wv, 12, 3, "=COUNTIF('现金流量表'!C10:G10,\"<0\")+1", "link", "0")
label(wv, 14, "情景加权 NPV（SUMPRODUCT）", bold=True)
header_row(wv, 15, ["情景", "概率", "NPV"])
for i, (name, p, mult) in enumerate([("悲观", 0.25, 0.7), ("基准", 0.5, 1.0), ("乐观", 0.25, 1.3)]):
    wv.cell(row=16 + i, column=1, value=name).font = F_LABEL
    put(wv, 16 + i, 2, p, "input", PCT)
    put(wv, 16 + i, 3, f"=$C$10*{mult}", "formula")
label(wv, 19, "加权 NPV", bold=True)
put(wv, 19, 3, "=SUMPRODUCT(B16:B18,C16:C18)", "formula", bold=True, fill=FILL_TOTAL)
set_widths(wv, [28, 12, 14])

# ================================================================ 敏感性 ====
wsn = wb.create_sheet("敏感性")
title(wsn, "敏感性分析", "上表：公式实时计算（折现率一维）；下表：Python 预计算数值（二维）")
header_row(wsn, 4, ["折现率", "NPV（含初始投资）", "相对基准"])
for i, rate in enumerate([0.08, 0.10, 0.12, 0.14, 0.16]):
    put(wsn, 5 + i, 1, rate, "input", PCT)
    put(wsn, 5 + i, 2, f"=NPV(A{5 + i},'现金流量表'!$C$9:$G$9)-初始投资", "link")
    put(wsn, 5 + i, 3, f"=IF('估值'!$C$10=0,0,B{5 + i}/'估值'!$C$10-1)", "link", PCT)
set_widths(wsn, [16, 20, 14, 14, 14, 14])

# ================================================================ 图表页 ====
wch = wb.create_sheet("图表")
title(wch, "图表（原生 Excel 图表，数据引用利润表/现金流量表）")
header_row(wch, 3, ["年份", "营业收入", "净利润", "自由现金流"])
for i, (y, c) in enumerate(zip(YEARS, COLS)):
    wch.cell(row=4 + i, column=1, value=str(y)).font = F_LABEL
    put(wch, 4 + i, 2, f"='利润表'!{c}5", "link")
    put(wch, 4 + i, 3, f"='利润表'!{c}14", "link")
    put(wch, 4 + i, 4, f"='现金流量表'!{c}9", "link")
line = LineChart()
line.title = "营业收入与净利润（万元）"
line.y_axis.title = "万元"
line.x_axis.title = "年份"
line.height, line.width = 8, 16
line.add_data(Reference(wch, min_col=2, min_row=3, max_col=3, max_row=3 + N), titles_from_data=True)
line.set_categories(Reference(wch, min_col=1, min_row=4, max_row=3 + N))
for s, color in zip(line.series, ["1F3864", "C9A227"]):
    s.graphicalProperties.line.solidFill = color
    s.graphicalProperties.line.width = 28000  # EMU
    s.marker.symbol = "circle"
    s.marker.graphicalProperties.solidFill = color
    s.smooth = False  # 不写 <smooth> 时 LibreOffice 会按平滑曲线渲染，显式关闭
line.y_axis.majorGridlines.spPr = None
line.y_axis.delete = False  # openpyxl 3.1 默认轴可能被隐藏，显式打开
line.x_axis.delete = False
wch.add_chart(line, "F3")
bar = BarChart()
bar.type = "col"
bar.title = "自由现金流（万元）"
bar.height, bar.width = 8, 16
bar.add_data(Reference(wch, min_col=4, min_row=3, max_row=3 + N), titles_from_data=True)
bar.set_categories(Reference(wch, min_col=1, min_row=4, max_row=3 + N))
bar.series[0].graphicalProperties.solidFill = "2E75B6"
bar.series[0].invertIfNegative = False
bar.dataLabels = None
bar.y_axis.delete = False
bar.x_axis.delete = False
bar.legend = None
wch.add_chart(bar, "F21")
set_widths(wch, [10, 14, 14, 14])

# 打印/PDF 导出设置：横向、缩放到一页宽（LibreOffice 转 PDF 时也生效）
for sh in wb.worksheets:
    sh.page_setup.orientation = "landscape"
    sh.page_setup.paperSize = sh.PAPERSIZE_A4
    sh.sheet_properties.pageSetUpPr.fitToPage = True
    sh.page_setup.fitToWidth = 1
    sh.page_setup.fitToHeight = 0
    sh.print_options.horizontalCentered = True

# 让 Excel 打开时强制全量重算（python 写的公式没有缓存值）
wb.calculation = CalcProperties(fullCalcOnLoad=True)
wb.properties.creator = "bizgold-office-engine"

path = os.path.join(OUT, "03_财务模型.xlsx")
wb.save(path)
print("saved", path)
print("FORMULAE check:", {f: (f in FORMULAE) for f in ["NPV", "IRR", "SUMPRODUCT", "IF", "XLOOKUP", "IFS", "TEXTJOIN", "COUNTIF"]})
