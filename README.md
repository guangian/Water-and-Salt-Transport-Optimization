# 金策 BizGold —— 全国高校商业精英挑战赛 国金级智能体

> 把一支顶级商赛团队（行业研究员、消费者洞察、财务建模师、主笔、路演设计师、三位评委、质检员）装进一个本地 Web 工作台。输入赛道、命题 / 品牌 / 项目信息与资料，端到端产出全国总决赛一等奖（“国金”）水准的全套参赛材料，并按官方评分表反复自评迭代。

## 能做什么

| 产出 | 说明 |
|---|---|
| 策划书 / 商业计划书 | Word + PDF，章节与官方评分项一一对应，图表、三线表、GB/T 7714 参考文献、自动目录、匿名 |
| 财务模型 | Excel（带公式）：收入驱动、三表、盈亏平衡、单位经济、三情景、敏感性、估值、融资；策划类为预算表 + KPI + ROI |
| 路演 PPT | 16:9，18–24 页 + 备用页，一页一结论句，演讲者备注 |
| 讲稿 / 答辩题库 / 视频分镜 | 10 分钟讲稿（字数与段落时间受控）、≥ 40 题六维度题库（结论–证据–行动）、陈述视频分镜 |
| 调研工具包 | 问卷（≥ 25 题）、访谈提纲、执行计划（样本量 / 配额 / 质控）、示例结果（明确标注，供替换） |
| 评审与合规 | 三位不同模型的评委按量规打分 + 机器检查（匿名 / 页数 / 来源 / 图表 / 禁语 / 财务勾稽 / 跨材料一致性），不达标自动返工 |
| 提交包 | 按赛道命名规则重命名、声明与授权书模板、数据真实性清单、README |

支持赛道（v1）：品牌策划（含留学生）、创业计划（命题 / 创业 / 创客）、数字营销、案例分析（会计案例 / 流通业）、商务会奖旅游。新增赛道只需在 `agent/tracks/` 加一个目录。

## 架构一句话

Claude Agent SDK（Claude Code 同源内核）+ 超详细领域规范（`agent/CLAUDE.md` + 赛道规范 + 子代理 + 技能 + 知识库）+ Python 确定性渲染工具包（Word / PPT / Excel / 图表）+ Fastify 服务 + React 工作台。详见 `docs/architecture.md`、`docs/contracts.md`。

## 快速开始（Windows）

```powershell
# 0. 前置：Node 20+、pnpm、Python 3.10+、Git for Windows（Bash 工具）；任意电脑先跑一键体检，缺什么它会直接给修复命令
powershell -ExecutionPolicy Bypass -File scripts\doctor.ps1   # Linux/macOS：bash scripts/doctor.sh
powershell -ExecutionPolicy Bypass -File scripts\setup.ps1     # pnpm install + pip install + 自检
# 1. 密钥：系统环境变量 BIZGOLD_API_BASE_URL / BIZGOLD_API_KEY，或启动后在「设置」页填写
pnpm dev        # 开发模式：http://localhost:5173（后端 5178）
# 或
pnpm build; pnpm start   # 生产模式：http://localhost:5178
# 或（Windows 常驻，推荐长时间运行）
pnpm build; powershell -ExecutionPolicy Bypass -File scripts\service.ps1 -Install   # 注册为计划任务 BizGoldServer
```

Linux / macOS：`bash scripts/setup.sh`。

## 目录

```
agent/        智能体大脑：CLAUDE.md、tracks/、.claude/{agents,skills}、toolkit/bizkit、knowledge/
apps/server   Fastify + Claude Agent SDK 运行器 + 进程内 MCP 工具（web_search / web_fetch / update_progress / record_numbers）
apps/web      React + Vite 工作台
docs/         架构、契约、调研报告、评审记录、用户手册
scripts/      setup / dev / start / 离线分发打包
workspace/    运行期项目（自动创建，已 gitignore）
```

## 重要说明

- 本系统是**辅助工作台**：赛事细则禁止“完全使用 AI 工具创作”。示例模式下的一手调研数据全部带【示例数据，须替换为真实调研结果】标记并列入《数据真实性清单》，团队须用真实调研替换并自行修改确认后再提交。
- 匿名制：产出不含院校 / 姓名 / 联系方式；文件名中的“学校名称”由团队提交时填写。
- 密钥只保存在 `%USERPROFILE%\.bizgold\settings.json` 或环境变量，仓库内不含任何密钥。

## 许可证

MIT
