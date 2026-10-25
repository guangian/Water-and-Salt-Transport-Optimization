---
name: docx-report
description: 用 bizkit.docx_render 把章节 Markdown 渲染成国金级 Word 策划书 / 商业计划书（封面、自动目录、章节编号、图注表注、三线表或斑马表、GB/T 7714 参考文献、匿名清洗、可选 PDF）。撰写正文或需要生成 / 重新生成 Word、控制页数时使用。
---

# docx-report：章节 Markdown → Word

## 1. 章节 Markdown 子集（渲染器只认这些）

| 语法 | 渲染结果 |
|---|---|
| `# 标题` / `## 小节` / `### 条目` | 自动编号 1 / 1.1 / 1.1.1（作者写的“第一章”“1.2”会被剥掉重编）；一级标题另起一页 |
| 普通段落 | 宋体小四、1.5 倍行距、首行缩进 2 字符 |
| `**粗体**`、`*斜体*` | 对应字体样式（每页加粗 ≤ 2 处） |
| `- 项` / `1. 项` | 列表 |
| `表：表题` + 下一行开始的 GFM 表格 | “表 N-M 表题”在表上；表后可跟 `{资料来源：…}` |
| `![图题](../04_图表/fig_x.png)` | 居中图 + “图 N-M 图题”在图下；路径相对章节 md 所在目录 |
| `[3]`、`[2-4]` | 上标文献引用，编号 = `01_调研/参考文献.json` 的 id |
| `{资料来源：…}` 单独一行 | 小号灰色来源行（图内已有来源时不必重复） |
| `> 引用` | 引用块（用于洞察句 / 定位陈述句） |
| `---` 或 `<!-- pagebreak -->` | 分页 |
| `[^1]` + 文末 `[^1]: 内容` | 章末注释 |

不要使用 HTML、嵌套列表、代码块、多级引用；表格列数 ≤ 7，单元格不换行；每个表 ≤ 25 行（更长的拆分或移附录）。

## 2. book.json（契约 §7.2）

```jsonc
{ "track": "brand", "docType": "品牌策划书",
  "cover": { "title": "宋韵落天青，一釉生万象", "subtitle": "弘宝汝瓷品牌策划书", "teamName": "青瓷工作室",
             "company": "指导企业：汝州市弘宝汝瓷坊有限公司", "group": "本科组", "date": "2026 年 3 月" },
  "abstract": "05_正文/00_概要.md",
  "chapters": [ { "file": "05_正文/01_引言.md" }, { "file": "05_正文/02_品牌现状分析.md" } ],
  "appendices": [ { "file": "05_正文/A1_调研问卷.md" } ],
  "references": "01_调研/参考文献.json",
  "style": { "bodyFont": "宋体", "bodySize": 12, "headingFont": "黑体", "latinFont": "Times New Roman", "lineSpacing": 1.5,
             "firstLineIndent": 2, "toc": true, "tocDepth": 2, "header": "弘宝汝瓷品牌策划书", "tableStyle": "banded", "h1Align": "center" },
  "anonymize": { "forbidden": ["某某大学", "某某学院"], "whitelist": ["大学生"] },
  "output": "产出/弘宝汝瓷品牌策划书.docx", "pdf": true }
```
- 封面字段按 `tracks/<赛道>/naming.json.coverFields`：品牌策划 = 作品名 / 团队名 / 指导企业名；创业计划 = 组别 / 团队名 / 项目名（`company` 留空）。
- `anonymize.forbidden` 填任务书与资料里出现的院校词；渲染结果 `anonymizeHits` 必须为 0（命中会被替换为“XX”，但你应回到 md 源头改掉）。
- 创业计划 `style.bodyFont` 必须 `宋体`、`bodySize` 12（细则硬性要求）。

## 3. 命令与结果

```bash
export PYTHONPATH=toolkit PYTHONUTF8=1
"${BIZGOLD_PYTHON:-python}" -m bizkit.docx_render 05_正文/book.json          # 输出 docx（pdf:true 时尝试转 PDF）
"${BIZGOLD_PYTHON:-python}" -m bizkit.pdf --pages 产出/xxx.pdf              # 真实页数
```
结果 JSON 关注：`pages`（真实）或 `pagesEstimate`（= 字数/650 + 图数×0.45 + 表行数/40）、`figures`、`tables`、`anonymizeHits`、`missingImages`、`warnings`。

## 4. 页数控制策略

目标 = 任务书 pageTarget（±10%），上限见 `checks.json.pageLimit`（品牌 / 创业 80，会计案例 50）。
- 超页：删除重复或弱论证的图（保留“结论句决定的图”）→ 合并同类表 → 问卷附录缩小为“题目 + 选项”紧凑表 → 删除过渡段。
- 不足：补“读图两句”、基线对照、敏感性、竞品案例卡片、一手调研交叉表；**不**用放大图片、空行凑页。
- 一级标题另起一页会产生留白，章节数不要超过 outline。

## 5. 常见问题
- “【缺图：…】”：路径写错或图未生成，先 `ls 04_图表`。
- 目录页码为估算：本机没有 PDF 转换器；提交前用 Word 打开按 F9 更新域即可。
- 中文变方块：仅图表会出现，见 biz-charts 技能。
