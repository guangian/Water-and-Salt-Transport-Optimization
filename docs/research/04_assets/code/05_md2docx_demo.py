# -*- coding: utf-8 -*-
"""Markdown→docx 流水线演示 + pandoc / docxtpl 对比。"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from md2docx import build_document  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

CH1 = """# 市场分析

## 市场规模与增长

智慧养老市场正处于**高速成长期**。根据艾瑞咨询数据，2021—2026 年市场规模的年复合增长率达到 *37.6%*，其中软件订阅占比持续提升（参见 `图 1`）。

![智慧养老市场规模（亿元）](charts/01_市场规模堆叠柱.png "width=13cm")

核心驱动因素包括：

- 人口老龄化加速：60 岁以上人口占比 **21.1%**
- 政策支持：民政部“十四五”养老服务规划明确智慧养老方向
  - 地方配套补贴覆盖 28 个省份
  - 社区适老化改造专项资金
- 技术成熟：毫米波雷达成本三年下降 60%

## 竞争格局

表：主要竞争对手对比（2025 年）

| 公司 | 核心产品 | 营收（万元） | 市场份额 |
|:-----|:--------|------------:|--------:|
| A 公司 | 穿戴手环 | 12,500 | 28.5% |
| B 公司 | 摄像头监护 | 9,200 | 21.0% |
| 本项目 | 毫米波雷达 | 1,350 | 6.5% |

> 评审提示：对比表应突出本项目差异化，而非罗列参数。

实施步骤如下：

1. 完成 3 个社区试点（2026 Q1）
2. 取得医疗器械二类备案
3. 启动 A 轮融资，详见 [融资计划](https://example.com/financing)

<!-- pagebreak -->

# 财务预测

## 收入预测

| 年份 | 2024A | 2025E | 2026E |
|------|------:|------:|------:|
| 营业收入 | 480 | 980 | 1,800 |
| 净利润 | -120 | 60 | 260 |

![现金流预测（万元）](charts/08_现金流组合.png)

```
假设：客单价 1,200 元 → 1,000 元（年降 4%）
      订阅转化率 60% → 75%
```
"""

path = build_document([CH1], os.path.join(OUT, "05_md渲染示例.docx"), base_dir=OUT, title="智守康联商业计划书",
                      subtitle="（Markdown 渲染流水线示例）", meta=[("团队", "智守康联"), ("日期", "2026-10")],
                      header_text="智守康联商业计划书")
print("saved", path)

# ---- pandoc 对比：同一份 Markdown ----
md_path = os.path.join(OUT, "05_chapter.md")
open(md_path, "w", encoding="utf-8").write(CH1)
r = subprocess.run(["pandoc", md_path, "-o", os.path.join(OUT, "05_pandoc.docx"), "--toc", "--resource-path", OUT],
                   capture_output=True, text=True)
print("pandoc rc:", r.returncode, r.stderr.strip()[:200])

# ---- docxtpl 对比：Jinja2 模板填充 ----
from docx import Document  # noqa: E402
from docxtpl import DocxTemplate  # noqa: E402

tpl_doc = Document()
tpl_doc.add_heading("{{ project }} 报名信息表", 1)
tpl_doc.add_paragraph("团队：{{ team }}；指导教师：{{ advisor }}")
t = tpl_doc.add_table(rows=4, cols=3)  # {%tr for %} 与 {%tr endfor %} 必须各占独立的一行
t.style = "Table Grid"
t.cell(0, 0).text, t.cell(0, 1).text, t.cell(0, 2).text = "姓名", "角色", "学院"
t.cell(1, 0).text = "{%tr for m in members %}"
t.cell(2, 0).text, t.cell(2, 1).text, t.cell(2, 2).text = "{{ m.name }}", "{{ m.role }}", "{{ m.college }}"
t.cell(3, 0).text = "{%tr endfor %}"
tpl_path = os.path.join(OUT, "05_docxtpl_模板.docx")
tpl_doc.save(tpl_path)
dt = DocxTemplate(tpl_path)
dt.render({"project": "智守康联", "team": "智守康联团队", "advisor": "张三",
           "members": [{"name": "李华", "role": "CEO", "college": "经管学院"}, {"name": "王敏", "role": "CTO", "college": "电子学院"}]})
dt.save(os.path.join(OUT, "05_docxtpl_输出.docx"))
d = Document(os.path.join(OUT, "05_docxtpl_输出.docx"))
print("docxtpl:", d.paragraphs[0].text, "| 表格行数:", len(d.tables[0].rows), [c.text for c in d.tables[0].rows[1].cells])
