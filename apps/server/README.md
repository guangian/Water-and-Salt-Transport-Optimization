# @bizgold/server — 金策 BizGold 服务端

Fastify 5 API + Claude Agent SDK 运行器。负责项目/设置持久化（JSON 文件）、把 Agent 消息流归一化为契约事件（SSE），并给 Agent 提供 `bizgold` MCP 工具（联网检索、网页抓取、进度上报、数字账本）。

依赖契约：`docs/contracts.md`（§1 设置、§2 项目、§3 接口、§4 事件、§5 任务书、§6 进度协议）、`docs/architecture.md` §3.1。

## 启动

```bash
pnpm install                                  # 仓库根目录
pnpm --filter @bizgold/server dev             # tsx watch，默认 http://127.0.0.1:5178
pnpm --filter @bizgold/server typecheck
pnpm --filter @bizgold/server test            # vitest（59 个用例，不访问网络/网关）
pnpm --filter @bizgold/server build && pnpm --filter @bizgold/server start
```

Node ≥ 20（开发环境 22）。若 `apps/web/dist/index.html` 存在则同源托管前端（SPA fallback 到 `index.html`）；否则只提供 `/api`。

## 环境变量

| 变量 | 作用 |
|---|---|
| `BIZGOLD_PORT` / `BIZGOLD_HOST` | 监听端口（5178）/ 地址（默认 `127.0.0.1`，不暴露到局域网） |
| `BIZGOLD_SETTINGS_DIR` | 设置目录，默认 `~/.bizgold`（`settings.json`，0600） |
| `BIZGOLD_WORKSPACE_DIR` | 工作区根，默认仓库 `workspace/`（已 gitignore） |
| `BIZGOLD_API_BASE_URL` / `BIZGOLD_API_KEY` / `BIZGOLD_MODEL` | 网关地址 / 密钥 / 主模型；**环境变量优先于设置文件**（契约 §9），`GET /api/settings` 以 `apiKeySource: env` 标识 |
| `ANTHROPIC_AUTH_TOKEN` | 等价于 `BIZGOLD_API_KEY`（备用名） |
| `BIZGOLD_PYTHON` | Python 解释器全路径（否则按 设置 → PATH `python`/`python3` 探测，跳过 Windows Store 占位符） |
| `BOCHA_API_KEY` / `TAVILY_API_KEY` / `SERPER_API_KEY` | 搜索提供商密钥（可选；没有则用 Bing 网页 / DuckDuckGo lite） |
| `BIZGOLD_LOG_LEVEL` | `debug` 时输出 CLI stderr 与每个请求日志 |

密钥只存在于设置文件/环境变量中：所有写入日志与 `events.jsonl` 的文本先经过 `redact()`（登记的密钥、`sk-…`、`Authorization/Bearer`、`*_API_KEY=` 模式）。

## 接口（契约 §3，前缀 `/api`）

| 方法 | 路径 | 说明 |
|---|---|---|
| GET | `/health` | Python / PDF 转换器（`bizkit.pdf --detect`）/ 中文字体（`bizkit.fonts`）/ Git Bash（Windows）/ 密钥来源；结果缓存 20s，`?force=1` 重新探测 |
| GET / PUT | `/settings` | 读取（apiKey 脱敏为 `sk-****后四位` + `apiKeySource`）/ 合并写入（脱敏串原样回传不会覆盖真值） |
| POST | `/settings/test` | `GET <baseUrl>/v1/models` 连通性测试，失败回退 `fallbackBaseUrl` |
| GET | `/tracks` | 五赛道参数：`agent/tracks/<id>/track.json` 优先，缺省用内置常量 |
| GET | `/knowledge`, `/knowledge/:slug` | `agent/knowledge/*.md` 与各赛道 `spec.md` |
| GET / POST | `/projects` | 列表 / 创建（201，id `p_YYYYMMDD_HHmmss_xxxx`） |
| GET / PATCH / DELETE | `/projects/:id` | 详情（含 `running`）/ 修改（运行中仅允许 `name`）/ 删除（需 `?confirm=1`，运行中先停止） |
| POST / DELETE | `/projects/:id/files` | multipart 上传到 `inputs/` 并自动 `bizkit.preread` → `inputs/_preread/<name>.txt` / 删除 |
| POST | `/projects/:id/run` | `{ mode: "run"|"resume", fromStage?, smoke?, maxTurns? }`；运行中 409 |
| POST | `/projects/:id/stop` | AbortController 终止子进程 → `done{subtype:"stopped"}` |
| POST | `/projects/:id/chat` | `{ message }` → modify 模式（resume 会话 + 「不得修改 inputs/」约束） |
| GET | `/projects/:id/events?since=` | SSE：回放 `events.jsonl` 中 `seq > since` 的事件再实时推送；`: ping` 每 15s；未运行时额外发 `event: idle` |
| GET | `/projects/:id/artifacts` | 文件树（排除 `.claude-home` `toolkit` `node_modules` `.cache` `.git` `__pycache__`） |
| GET | `/projects/:id/file?path=&download=` | 原样下载/内联；路径解析到工作区内后再做 realpath 前缀校验（符号链接逃逸 → 404） |
| GET | `/projects/:id/preview?path=` | markdown（相对图片改写为 `/file?path=`）/ image / pdf / slides（`.preview.json`）/ table（xlsx 每表前 200 行、csv）/ text / binary；同样做 realpath 校验 |
| GET | `/projects/:id/review` | `07_评审/第N轮/{评审汇总,机器检查,评审_*}.json`，缺省为空 |
| GET | `/projects/:id/numbers` | `03_财务/numbers.json` + 账本自洽性检查（重复 key、缺来源、非数值） |
| POST | `/projects/:id/package` | `python -m bizkit.package --track … --workspace …`；Python 不可用时退化为 archiver 压缩 `产出/` |
| POST | `/projects/:id/open-folder` | explorer / open / xdg-open 打开工作区（失败不报错） |

错误统一为 `{ error: { code, message } }`。

## 设计要点

- **工作区准备**（`agent/workspace.ts`）：每次运行把 `agent/{CLAUDE.md,.claude,toolkit,knowledge,tracks/<track>}` 复制进 `workspace/<id>/`（跳过 `__pycache__/examples/tests`），渲染 `任务书.md`（§5），并把 `.claude/agents/*.md` 里的 `model: {{MODEL_xxx}}` 占位符替换为设置中的模型。
- **query() 选项**（`agent/runner.ts`）：`claude_code` 预设 + 环境说明、`settingSources: ['project']`、`bypassPermissions`、`includePartialMessages`、`persistSession`、`strictMcpConfig`、`disallowedTools: ['WebSearch']`；子进程 env 只用 `ANTHROPIC_AUTH_TOKEN`（Bearer），`CLAUDE_CONFIG_DIR=<ws>/.claude-home` 隔离会话，关闭遥测/自动更新/自动记忆，`ANTHROPIC_DEFAULT_HAIKU_MODEL=models.fast`、`CLAUDE_CODE_SUBAGENT_MODEL=models.main`。网关 400 时在设置里打开 `api.disableExperimentalBetas`（→ `CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1`）；本机对 llm.goaichat.top 冒烟无需该开关。
- **事件归一化**（`agent/events.ts`）：子代理按 `parent_tool_use_id` 归属，名称来自 `Agent/Task` 的 `input.subagent_type`（映射为中文）或 `description`；`text_delta` 按 16ms 合并、思考增量按 200ms 合并，并在完整块/工具调用前冲出以保证顺序；`Write/Edit` 成功后的 `tool_result` 产生 `file` 事件；`result.modelUsage` × 设置价格表估算费用（第三方网关的 `total_cost_usd` 不可信，另存为 `sdkCostUsd`）。
- **进度**：`update_progress` MCP 工具是主通道（写 `project.json.stages` + `stage` 事件）；主代理文本里的 §6 清单作兜底，只允许「推进」状态（done 不会被 pending/running 覆盖）。
- **持久化**：`project.json` 原子写（临时文件 + rename）；`events.jsonl` 每项目串行追加、`seq` 单调递增（重启后从文件尾恢复）。`text_delta` 与 `thinking{delta:true}` 只走 SSE 不落盘——它们随后总以完整块到达，落盘只会让长跑的回放膨胀几十倍（实测 212 条实时事件对应 33 条落盘）。
- **检索**（`services/search.ts`）：博查 → Tavily → Serper → Bing 网页 → DuckDuckGo lite，没有 key 的 API 提供商直接跳过；网页型结果需至少一条命中查询词片段（Bing 对无 Cookie 请求会返回整页无关结果），否则交给下一家；内存 + `<ws>/.cache/search.json` 缓存 24h。
- **抓取**（`services/fetch.ts`）：浏览器 UA、20s 超时、2MB 上限、拒绝内网/回环/链路本地/CGNAT（含 DNS 解析后校验）；HTML → Markdown（保留标题/段落/列表/表格/代码块，去掉脚本/导航/页脚），PDF 交给 `bizkit.preread`。
- **Python**（`services/python.ts`）：不经 shell 的 `spawn`，强制 UTF-8，`PYTHONPATH` 指向 `<ws>/toolkit` 或 `agent/toolkit`；首次启动检查 `docx/pptx/openpyxl/matplotlib`，缺失则后台 `pip install -r requirements.txt -i <pipIndex>`。
- **重启恢复**：启动时把上次异常退出留下的 `running` 项目标为 `paused`（进行中的阶段记 `failed` 并留注释），并追加一条 `status` 事件。
- **Windows**：全部路径用 `path.join/resolve`，`safeResolve` 在 Windows 下不区分大小写比较；探测 Git Bash 并注入 `CLAUDE_CODE_GIT_BASH_PATH`；`open-folder` 用 `explorer.exe`。

## 已知局限

- 子代理归属依赖 `parent_tool_use_id`；SDK 若对后台子代理（`async_launched`）不转发消息，则只能看到其 start/end。
- `file` 事件来自 `Write/Edit` 工具结果，不监听文件系统：Bash/Python 脚本直接写出的文件只会在阶段 `done` 时被 `scanArtifacts` 收入 `project.artifacts`。
- 费用估算依赖设置里的 `api.pricing`（USD / 百万 token），默认空表 → `costEstimateUsd = 0`。
- Bing/DuckDuckGo 网页抓取属于尽力而为：被限流或页面结构变化时会落到「无结果」，需配置 API 提供商 key 才稳定。
- `/health` 首次探测要调用 Python（字体/PDF 转换器检测），耗时数秒；之后 20s 内走缓存。
- 事件回放会把 `events.jsonl` 整个读入内存（无分页）；超长运行建议前端带 `since` 增量拉取。
