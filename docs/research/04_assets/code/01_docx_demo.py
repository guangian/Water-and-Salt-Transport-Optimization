# -*- coding: utf-8 -*-
"""python-docx 1.2.0 全功能演示：生成 out/01_商业计划书示例.docx"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from docx import Document  # noqa: E402
from docx.enum.text import WD_ALIGN_PARAGRAPH  # noqa: E402
from docx.shared import Pt  # noqa: E402

import docx_helpers as H  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
CHART1 = os.path.join(OUT, "charts", "01_市场规模堆叠柱.png")
CHART2 = os.path.join(OUT, "charts", "08_现金流组合.png")

doc = Document()
H.setup_chinese_document(doc, body_east_asia="宋体", heading_east_asia="黑体", latin="Times New Roman", body_size=12)
H.add_heading_numbering(doc, levels=3)
fn = H.Footnotes(doc)
cap = H.CaptionCounter()

# ---------- 第 1 节：封面 + 目录（无页码） ----------
sec_cover = doc.sections[0]
H.add_cover_page(doc, "智守康联：居家养老 AI 健康监测平台", "商业计划书",
                 meta=[("参赛团队", "智守康联团队"), ("指导教师", "张三 教授"), ("日期", "2026 年 10 月")],
                 new_section_after=False)
doc.add_page_break()
H.add_toc(doc, levels="1-3")

# ---------- 第 2 节：正文（页码从 1 开始，页眉页脚） ----------
sec_body = doc.add_section()
H.set_header_text(sec_body, "智守康联商业计划书  |  2026 全国高校商业精英挑战赛")
H.add_page_number_footer(sec_body, "第 {PAGE} 页 / 共 {NUMPAGES} 页")
H.restart_page_numbering(sec_body, start=1)

h = doc.add_heading("执行摘要", level=1)
H.add_bookmark(h, "exec_summary")
H.add_body_paragraph(doc, "智守康联聚焦独居老人的居家健康风险，提供“硬件 + 订阅服务”的一体化解决方案。"
                     "2026 年目标覆盖 3 万户家庭，实现营业收入 1,800 万元，净利润 260 万元。",
                     bold_spans=["3 万户家庭", "1,800 万元"])
p = H.add_body_paragraph(doc, "本方案中的市场数据主要来自艾瑞咨询与国家统计局公开报告")
fn.add(p, "艾瑞咨询《2025 年中国智慧养老行业研究报告》，2025 年 6 月。")
p.add_run("，行业规模测算采用自下而上法。")
H.set_run_font(p.runs[-1])

doc.add_heading("项目概述", level=1)
doc.add_heading("行业痛点", level=2)
H.add_bullets(doc, ["独居老人意外跌倒后平均 4.5 小时才被发现", "现有穿戴设备佩戴依从性低于 40%",
                    "社区养老服务人手不足，人均服务 200 户以上"])
doc.add_heading("解决方案", level=2)
H.add_numbered(doc, ["毫米波雷达无感监测，无需佩戴", "边缘 AI 模型 3 分钟内识别跌倒/长时间静止",
                     "家属 App + 社区大屏双端联动"])
H.add_quote_block(doc, "“我们不是在卖设备，而是在卖一份让子女安心的确定性。” —— 团队创始人")

doc.add_heading("市场分析", level=1)
doc.add_heading("市场规模", level=2)
H.add_body_paragraph(doc, "如图 1 所示，智慧养老市场规模自 2021 年的 160 亿元增长至 2026 年的 790 亿元，"
                     "年复合增长率达 37.6%。")
H.add_figure(doc, CHART1, "智慧养老市场规模（2021—2026E，亿元）", cap, width_cm=14)
doc.add_heading("竞争格局", level=2)
H.add_body_paragraph(doc, "主要竞争对手对比如表 1 所示（三线表）。")
H.add_caption(doc, "表", "主要竞争对手对比", cap)
H.add_three_line_table(doc, ["公司", "核心产品", "2025 年营收（万元）", "市场份额"],
                       [["A 公司", "穿戴手环", "12,500", "28.5%"], ["B 公司", "摄像头监护", "9,200", "21.0%"],
                        ["本项目", "毫米波雷达", "1,350", "6.5%"]], widths_cm=[3, 4, 4.5, 3], header_fill="D9E2F3")
doc.add_paragraph()
H.add_body_paragraph(doc, "带底纹与合并单元格的全框线表格示例见表 2。")
H.add_caption(doc, "表", "三年财务摘要（万元）", cap)
t = H.add_grid_table(doc, ["指标", "2024A", "2025E", "2026E"],
                     [["营业收入", "480", "980", "1,800"], ["毛利率", "33%", "38%", "42%"], ["净利润", "-120", "60", "260"],
                      ["合并行示例", "", "", ""]], widths_cm=[4, 3.5, 3.5, 3.5])
merged = t.cell(4, 1).merge(t.cell(4, 3))
H._format_cell_text(merged, "此单元格由 3 列合并而成", size=10.5)
H.set_repeat_table_header(t.rows[0])

doc.add_heading("财务预测", level=1)
doc.add_heading("现金流", level=2)
H.add_body_paragraph(doc, "现金流预测见图 2，项目预计在 2027 年实现累计现金流转正。")
H.add_figure(doc, CHART2, "现金流预测（万元）", cap, width_cm=14)
doc.add_heading("段落格式示例", level=2)
H.add_body_paragraph(doc, "本段使用固定行距 20 磅、段前 6 磅、段后 6 磅、首行缩进 2 字符。" * 3,
                     line_spacing=Pt(20), space_before=6, space_after=6)
H.add_body_paragraph(doc, "本段不缩进且两端对齐，1.5 倍行距。" * 4, first_line_indent=False,
                     align=WD_ALIGN_PARAGRAPH.JUSTIFY)
p = H.add_body_paragraph(doc, "更多信息见团队官网：")
H.add_hyperlink(p, "https://example.com/zhishou", url="https://example.com/zhishou")
p.add_run("；返回 ")
H.set_run_font(p.runs[-1])
H.add_hyperlink(p, "执行摘要", anchor="exec_summary")
p.add_run("（文档内部链接）。")
H.set_run_font(p.runs[-1])
# python-docx 1.2.0 新增：批注（评审意见）
doc.add_comment(runs=p.runs[0], text="请补充官网上线时间。", author="评审老师", initials="PS")

# ---------- 第 3 节：横向页 ----------
sec_land = H.add_landscape_section(doc)
doc.add_heading("附录：横向页宽表", level=1)
H.add_caption(doc, "表", "五年利润表预测（横向页）", cap)
H.add_three_line_table(doc, ["项目", "2024A", "2025E", "2026E", "2027E", "2028E", "CAGR"],
                       [["营业收入", "480", "980", "1,800", "3,000", "4,500", "75%"],
                        ["营业成本", "320", "610", "1,040", "1,650", "2,400", "66%"],
                        ["净利润", "-120", "60", "260", "520", "900", "—"]],
                       widths_cm=[4, 3, 3, 3, 3, 3, 3])

# ---------- 第 4 节：恢复纵向 ----------
H.add_portrait_section(doc)
doc.add_heading("参考文献", level=1)
H.add_numbered(doc, ["艾瑞咨询. 2025 年中国智慧养老行业研究报告[R]. 2025.",
                     "国家统计局. 2024 年国民经济和社会发展统计公报[R]. 2025."])

path = os.path.join(OUT, "01_商业计划书示例.docx")
doc.save(path)
print("saved", path)
