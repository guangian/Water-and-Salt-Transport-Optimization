# Claude Agent SDK（TypeScript）技术手册 —— 面向实现者

| 项 | 内容 |
| --- | --- |
| 适用版本 | `@anthropic-ai/claude-agent-sdk@0.3.310`，内置 Claude Code CLI 2.1.310（`SDK/package.json` 的 `version` / `claudeCodeVersion`；`SDK/manifest.json` 的 `version: 2.1.310`） |
| 事实优先级 | 本地安装文件（`sdk.d.ts`、`sdk.mjs`、`README.md`、`manifest.json`、平台二进制包）＞ 官方文档（`https://code.claude.com/docs/en/...`，即用户给出的 `docs.claude.com/en/api/agent-sdk/...` 与 `docs.claude.com/en/docs/claude-code/...` 的重定向目标）＞ `/tmp/sdktest/smoke_*.log` 与会话 JSONL |
| 方法 | 纯静态分析 + 文档阅读；未调用任何模型 API |
| 标记约定 | 『本版本不存在』= 在 0.3.310 的 `.d.ts` / `sdk.mjs` / 二进制里 grep 不到；『已弃用』= `.d.ts` 标注 `@deprecated` |
| 引用缩写 | `SDK/` = `/tmp/sdktest/node_modules/@anthropic-ai/claude-agent-sdk/`；`BIN/` = `/tmp/sdktest/node_modules/@anthropic-ai/claude-agent-sdk-linux-x64/`；`T/` = `/tmp/sdktest/`；`JSONL` = `T/home/.claude/projects/-tmp-sdktest-ws/42e923ab-a86f-4c9e-9491-65c4b8a6fa63.jsonl`（glm-5.3 会话）；`DOC/` = `/tmp/docs/`（官方文档的 Markdown 快照，URL 另注） |

目录

0. 架构速览：query() 到底启动了什么
1. `query()` Options 完整字段表
2. 消息流：SDKMessage 每个成员的字段与用法
3. 项目级配置：CLAUDE.md / agents / skills / commands / settings.json
4. 第三方 Anthropic 兼容网关、环境变量、Windows 与二进制定位
5. 最小可运行示例（Node 22 ESM TypeScript）+ Fastify/SSE 骨架
6. 版本说明（0.3.x 相对早期版本的关键变化）

---

## 0. 架构速览：query() 到底启动了什么

来源：`SDK/sdk.d.ts:2915-3325`、`SDK/sdk.mjs:157`（ProcessTransport）、`SDK/sdk.mjs:163`（`query()` 实现）、`SDK/package.json`、`SDK/README.md`、`JSONL:3-17`。

1. `query({ prompt, options })` 返回 `Query`，它是 `AsyncGenerator<SDKMessage, void>` 并附带控制方法（`SDK/sdk.d.ts:2915`、`3322-3325`）。
2. SDK **不在本进程内运行模型循环**，而是 spawn 一个原生二进制（bun 单文件，约 291 MB）：
   - 定位逻辑 `SDK/sdk.mjs:157` 函数 `N9`：候选包名 = `@anthropic-ai/claude-agent-sdk-${platform}-${arch}`（Linux 还按 musl/glibc 排序尝试 `-musl` 变体，Android 用 `-linux-${arch}-android`），文件名 = `claude`，Windows 追加 `.exe`；逐个 `createRequire(import.meta.url).resolve()` 后 `existsSync` 命中即返回。
   - musl 判定 `SDK/sdk.mjs:157`：`process.report.getReport().header?.glibcVersionRuntime === void 0`。
   - 找不到时抛错（`SDK/sdk.mjs:163`）：`Native CLI binary for ${process.platform}-${process.arch} not found. Reinstall @anthropic-ai/claude-agent-sdk without --omit=optional, or set options.pathToClaudeCodeExecutable.`
   - `pathToClaudeCodeExecutable` 可指向原生二进制或 JS 文件；是否 JS 由后缀判断（`SDK/sdk.mjs:163` 的 `FYe`：`.js/.mjs/.tsx/.ts/.jsx`），只有 JS 文件才会用 `executable`/`executableArgs` 指定的运行时启动。
3. 子进程以 `--output-format stream-json --verbose --input-format stream-json` 方式运行（`SDK/sdk.mjs:157`），Options 中一部分变成命令行 flag（`--model`、`--max-turns`、`--setting-sources`、`--allowedTools`、`--permission-mode`、`--resume` 等），另一部分通过 stdin 的 `initialize` 控制请求发送：`systemPrompt / appendSystemPrompt / systemPromptSnapshot / excludeDynamicSections / agents / hooks / sdkMcpServers / plugins(initialize 模式)` 等（`SDK/sdk.mjs:163`）。
4. 子进程环境：当 `options.env` 省略时继承 `process.env`；无论如何都会追加 `CLAUDE_CODE_ENTRYPOINT="sdk-ts"`、`CLAUDE_AGENT_SDK_VERSION="0.3.310"`，开启文件检查点时追加 `CLAUDE_CODE_ENABLE_SDK_FILE_CHECKPOINTING="true"`，并 `delete NODE_OPTIONS`（`SDK/sdk.mjs:157`）。会话 JSONL 每条记录的 `entrypoint:"sdk-ts"`、`version:"2.1.310"` 印证了这一点（`JSONL:3-13`）。
5. 省略 `systemPrompt` 时发往模型的 system 实际是一句话：`You are a Claude agent, built on Anthropic's Claude Agent SDK.`（`JSONL:14` 的 `api-request-shape.shape.system`），并非 Claude Code 完整提示词。要 Claude Code 行为必须显式 `systemPrompt: { type: 'preset', preset: 'claude_code' }`（第 1 节）。
6. 依赖：`peerDependencies` = `@anthropic-ai/sdk >=0.93.0`、`@modelcontextprotocol/sdk ^1.29.0`、`zod ^4.0.0`；`optionalDependencies` 为 8 个平台二进制包；`engines.node >= 18`；`exports` 含 `.`、`./core`、`./extract`、`./browser`、`./bridge`、`./sdk-tools`（`SDK/package.json`）。`SDK/README.md`：`/core` 入口从你的安装解析 `zod/v4`（兼容 `zod/v3`），根入口与 `/core` 不要混用；需要 TypeScript 5.0+；Windows 下二进制子路径是 `/claude.exe`。

---

## 1. `query()` Options 完整字段表

来源：`SDK/sdk.d.ts:1558-2509`（`Options` 类型本体，每行号即字段声明位置）；flag 映射见 `SDK/sdk.mjs:157`；官方参考 `https://code.claude.com/docs/en/agent-sdk/typescript`。
核验方法：`awk 'NR>=1558&&NR<=2509&&/^    [A-Za-z_]+\??: /' SDK/sdk.d.ts` 共得 **66 个**顶层字段，下表逐一列出（按用途分组，不遗漏）。「推荐值」针对"服务端常驻 + 第三方网关 + 前端 SSE"场景。

### 1.1 进程、运行时、调试

| 字段 | 类型 | 含义（d.ts 原意） | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `abortController` | `AbortController` | 取消查询并清理资源 | 每个 HTTP 请求一个；前端断连即 `abort()` | 1563 |
| `cwd` | `string` | 会话工作目录，默认 `process.cwd()` | 每个任务一个独立绝对路径的工作区 | 1635 |
| `additionalDirectories` | `string[]` | 除 cwd 外可访问的目录（绝对路径） | 需要跨目录读写时才设 | 1568 |
| `projectConfigRoot` | `string` | worktree 场景：hooks/permissions/.mcp.json/.claude 配置树与 `CLAUDE_PROJECT_DIR` 从此目录取而非 cwd | 一般不设 | 1578 |
| `env` | `{[k]: string \| undefined}` | 子进程环境变量；**设置后整体替换**，不与 `process.env` 合并；可设 `CLAUDE_AGENT_SDK_CLIENT_APP` 标识 UA | `{ ...process.env, ANTHROPIC_BASE_URL, ANTHROPIC_AUTH_TOKEN, ... }`（见第 4 节） | 1698 |
| `executable` | `'bun' \| 'deno' \| 'node'` | 仅当 `pathToClaudeCodeExecutable` 指向 JS 文件时使用的运行时 | 不设 | 1705 |
| `executableArgs` | `string[]` | 传给上述运行时的参数 | 不设 | 1709 |
| `extraArgs` | `Record<string, string \| null>` | 额外 CLI 参数（键不带 `--`，布尔 flag 用 `null`） | 仅用于 d.ts 未封装的 flag（见 `/tmp/claude_help.txt`） | 1715 |
| `pathToClaudeCodeExecutable` | `string` | 覆盖内置二进制路径 | 默认不设；离线打包/自带 `claude.exe` 时设绝对路径 | 2017 |
| `spawnClaudeCodeProcess` | `(SpawnOptions) => SpawnedProcess` | 自定义 spawn（VM/容器/远程） | 不设 | 2508 |
| `stderr` | `(data: string) => void` | 子进程 stderr 回调 | 必设，接入日志 | 2357 |
| `debug` | `boolean` | 等价 `--debug` | 排障时 `true` | 2347 |
| `debugFile` | `string` | 等价 `--debug-file <path>`，隐式开启 debug | 排障时设到 `/tmp` | 2352 |
| `loadTimeoutMs` | `number` | `sessionStore.load()` 超时，默认 60000（@alpha） | 不设 | 1885 |

### 1.2 提示词、模型、推理

| 字段 | 类型 | 含义 | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `systemPrompt` | `string \| string[] \| {type:'custom', prompt, snapshot?} \| {type:'preset', preset:'claude_code', append?, excludeDynamicSections?, snapshot?}` | 省略 → 极简默认提示词（`JSONL:14`）；字符串/数组 → 完全自定义（数组中可放 `SYSTEM_PROMPT_DYNAMIC_BOUNDARY` 划分可缓存前缀，`sdk.d.ts:10002`）；preset → Claude Code 完整提示词，`append` 追加；`excludeDynamicSections` 去掉 cwd/git/memory 等动态段并改放首条 user 消息；`snapshot` 把提示词记录进会话并在 resume 时原样复用（默认 true 行为） | 编码类任务：`{type:'preset', preset:'claude_code', append: 业务规则, snapshot: true}`；非编码业务：自定义字符串 | 2467 |
| `model` | `string` | 模型 ID，默认 CLI 默认模型 | 网关场景与 `env.ANTHROPIC_MODEL` 一致地显式传 | 2000 |
| `fallbackModel` | `string`（逗号分隔多个） | 主模型过载时的后备；每个用户轮次重试主模型。与 `model` 相同会抛错（`sdk.mjs:157`） | 网关通常不设 | 1722 |
| `thinking` | `{type:'adaptive', display?} \| {type:'enabled', budgetTokens?, display?} \| {type:'disabled'}` | 推理控制；优先于已弃用的 `maxThinkingTokens`（`sdk.d.ts:10053-10077`） | 第三方模型若 400/不支持 → `{type:'disabled'}`；支持者 `{type:'enabled', budgetTokens: N}` | 1944 |
| `effort` | `'low' \| 'medium' \| 'high' \| 'xhigh' \| 'max'` | 努力程度，配合 adaptive thinking（`sdk.d.ts:709`） | 网关下不设（默认请求带 `output_config.effort:"high"`，`JSONL:17`） | 1951 |
| `maxThinkingTokens` | `number` | 『已弃用』→ 用 `thinking` | 不用 | 1960 |
| `betas` | `SdkBeta[]`，当前仅 `'context-1m-2025-08-07'` | 启用 beta 头（`sdk.d.ts:3823`） | 网关不设 | 1754 |
| `outputFormat` | `{type:'json_schema', schema}` | 结构化输出；结果在 `result.structured_output`（`sdk.d.ts:6078`）；失败 subtype `error_max_structured_output_retries` | 需要机器可读结果时用 | 2013 |
| `taskBudget` | `{total: number}` | API 侧 token 预算（beta 头 `task-budgets-2026-03-13`，@alpha） | 不设 | 1978 |
| `maxTurns` | `number` | 最大代理轮次（user+assistant 为一轮） | 服务端必设（如 30），超限 result.subtype=`error_max_turns` | 1965 |
| `maxBudgetUsd` | `number` | 预算上限，超限 `error_max_budget_usd`；只计本次 `query()` 起（或上次 `/clear` 后）的花费 | 第三方模型成本为 `unknown` 基准时仅作兜底 | 1970 |
| `title` | `string` | 新会话标题；resume 时沿用已保存标题 | 可设业务 ID | 2486 |
| `verbatimPrompts` | `boolean` | 所有 user 消息按原文发送：不展开 `@path`、不分发 `/命令`（需 CLI ≥ 2.1.248） | 当 prompt 由非用户输入拼接时 `true` | 1931 |
| `planModeInstructions` | `string` | `permissionMode:'plan'` 时替换计划工作流正文 | 不设 | 2039 |

### 1.3 工具与权限

| 字段 | 类型 | 含义 | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `tools` | `string[] \| [] \| {type:'preset', preset:'claude_code'}` | 内置工具基础集合；`[]` 禁用全部内置工具。**注意**：原生版默认可能不提供 Grep/Glob，而用 Bash `find/grep` 替代——要它们必须在 `tools` 或 `allowedTools` 中列出（d.ts 1675-1676；`JSONL:14` 首轮工具列表确无 Grep/Glob） | 显式白名单，如 `['Read','Write','Edit','Bash','Grep','Glob','Agent','Skill']` | 1678 |
| `allowedTools` | `string[]` | 免提示自动放行的工具（可带规则，如 `Bash(npm test:*)`、`mcp__srv__*`、`Skill(name)`）。传 `'Skill'` 已弃用 → 用 `skills` | 与 `permissionMode:'dontAsk'` 组合形成硬白名单 | 1621 |
| `disallowedTools` | `string[]` | 从模型上下文移除的工具 | 服务端禁 `WebSearch`/`AskUserQuestion` 等 | 1641 |
| `toolAliases` | `Record<string,string>` | 模型发出的工具名重定向（单跳），如 `{Bash:'mcp__workspace__bash'}` | 远程沙箱时用 | 1667 |
| `toolConfig` | `ToolConfig` | 内置工具的逐项配置，如 `askUserQuestion.previewFormat` | 不设 | 1742 |
| `permissionMode` | `'default' \| 'acceptEdits' \| 'bypassPermissions' \| 'plan' \| 'dontAsk' \| 'auto'`（`sdk.d.ts:2548`） | 省略时按 `claude -p` 规则：settings 的 `permissions.defaultMode`，否则 `auto`（不可用时 `default`）。`dontAsk` = 未预批准即拒；`auto` = 分类器模型裁决；`bypassPermissions` 需同时 `allowDangerouslySkipPermissions:true` | 服务端：`'dontAsk'` + `allowedTools`；需要人工审批：`'default'` + `canUseTool`；`bypassPermissions` 仅限容器沙箱 | 2032 |
| `allowDangerouslySkipPermissions` | `boolean` | 使用 `bypassPermissions` 的安全开关 | 仅沙箱内 `true` | 2044 |
| `canUseTool` | `CanUseTool`（`sdk.d.ts:228-320`） | 工具权限回调 `(toolName, input, {signal, suggestions, blockedPath, mcpServer, decisionReason, title, toolUseID, agentID, requestId, ...}) => Promise<PermissionResult \| null>`；返回 `{behavior:'allow', updatedInput?, updatedPermissions?}` 或 `{behavior:'deny', message, interrupt?}`（`sdk.d.ts:2571`）。与 `permissionPromptToolName` 互斥（`sdk.mjs:157` 抛错，SDK 改传 `--permission-prompt-tool=stdio`）。字符串 prompt 模式也可用：设置了 canUseTool/hooks/进程内 MCP 时 SDK 会保持 stdin 直到本轮结束（`sdk.mjs:157` 的 `hasBidirectionalNeeds()`；字符串 prompt 在 `sdk.mjs:163` 被包装成单条 user 消息） | 人审场景用，并通过 SSE 推给前端 | 1626 |
| `permissionPromptToolName` | `string` | 用某个 MCP 工具接管权限提示 | 不设 | 2049 |
| `permissionPrompts` | `'host' \| 'none'` | `'none'`：任何会提示的调用立即拒绝，`canUseTool` 永不被调 | 无人值守服务 `'none'` | 2058 |
| `sandbox` | `SandboxSettings` | 命令沙箱（Linux 需 bubblewrap）；`enabled:true` 时 `failIfUnavailable` 默认 true | 可用时开启 | 2246 |
| `hooks` | `Partial<Record<HookEvent, HookCallbackMatcher[]>>` | 进程内钩子回调（事件名 33 个，见 1.10） | 至少 `PreToolUse`/`PostToolUse` 记审计日志 | 1768 |
| `includeHookEvents` | `boolean` | 把 `hook_started/hook_progress/hook_response` 系统消息发到流中（SessionStart/Setup 总是发） | 调试时 `true` | 1895 |
| `onElicitation` | `OnElicitation` | MCP 服务器请求用户输入时的回调；未处理则自动拒绝 | 不设 | 1789 |
| `onUserDialog` / `supportedDialogKinds` | 回调 / `string[]` | CLI 请求宿主渲染阻塞对话框；不声明 kind 则不会发 | 不设 | 1803 / 1825 |
| `perTaskStopAffordance` | `boolean` | 声明宿主有"按任务停止"控件，interrupt 时不杀后台任务 | 不设 | 1845 |

### 1.4 子代理、技能、插件、设置来源

| 字段 | 类型 | 含义 | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `agents` | `Record<string, AgentDefinition>` | 程序化子代理；同名覆盖文件系统 agents（`AgentDefinition` 见 1.9） | 服务端用它而不是 `.claude/agents` | 1613 |
| `agent` | `string` | 主线程使用某个 agent 的提示词/工具/模型（等价 `--agent`） | 需要"主代理即某角色"时用 | 1597 |
| `skills` | `string[] \| 'all'` | 开启技能并自动加 Skill 工具到 allowedTools（`sdk.mjs:157` 生成 `Skill(name)`）；省略 ≠ 关闭；只是上下文过滤，文件仍可被 Read | 显式列出技能名 | 2339 |
| `plugins` | `SdkPluginConfig[]`，仅 `{type:'local', path}`（`sdk.d.ts:5852`） | 加载本地插件 | 需要时用 | 2073 |
| `pluginDelivery` | `'argv' \| 'initialize'` | 插件路径走命令行（默认）还是 initialize 请求；Windows 命令行超 32,767 字符会拒绝启动 → 用 `'initialize'`（需 CLI ≥ 2.1.261） | Windows 多插件时 `'initialize'` | 2088 |
| `forwardSubagentText` | `boolean` | 转发子代理的 text/thinking（默认只转 tool_use/tool_result） | 需渲染子代理对话时 `true` | 1907 |
| `agentProgressSummaries` | `boolean` | 每 ~30s 为运行中的子代理生成一句进度摘要（`task_progress.summary`） | 可选 | 2120 |
| `settingSources` | `('user' \| 'project' \| 'local')[]` | 控制加载哪些文件系统设置；**省略 = 全部加载**；`[]` = 隔离模式；含 `'project'` 才加载 CLAUDE.md（d.ts 2312-2315；`sdk.mjs:157` 默认 `lYe=["user","project","local"]`） | 服务端必显式：`[]`（完全隔离）或 `['project']`（只认工作区配置） | 2316 |
| `settings` | `string \| Settings` | 文件路径或内联对象，进入最高优先级的 flag settings 层（等价 `--settings`） | 用内联对象注入 `permissions`/`env`/`skipWebFetchPreflight` | 2264 |
| `managedSettings` | `Settings` | 宿主注入的策略层设置（仅限制性键生效） | 不设 | 2305 |

### 1.5 MCP

| 字段 | 类型 | 含义 | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `mcpServers` | `Record<string, McpServerConfig>`，变体：stdio `{type?:'stdio', command, args?, env?, timeout?}`（d.ts 1370）、sse `{type:'sse', url, headers?, tools?, timeout?}`（1353）、http `{type:'http', url, headers?, tools?, timeout?, alwaysLoad?}`、sdk `createSdkMcpServer()` 返回的 `{type:'sdk', name, instance, tools?, timeout?}`（1222） | 工具名形如 `mcp__<server>__<tool>`，放行用 `allowedTools:['mcp__srv__*']` | 进程内工具一律用 `createSdkMcpServer` | 1995 |
| `strictMcpConfig` | `boolean` | 只用 `mcpServers`（和 `agents` 里声明的）MCP，忽略 `.mcp.json`/用户设置/插件/磁盘 agent frontmatter | 服务端 `true` | 2365 |

### 1.6 会话

| 字段 | 类型 | 含义 | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `continue` | `boolean` | 继续 cwd 下最近一次会话；与 `resume` 互斥 | 服务端不用（多租户不确定） | 1631 |
| `resume` | `string` | 按 session_id 恢复历史 | 多轮对话用，ID 来自 `init`/`result` 的 `session_id` | 2125 |
| `forkSession` | `boolean` | 恢复时分叉为新 session_id | 需保留原会话时 `true` | 1747 |
| `sessionId` | `string`（UUID） | 指定新会话 ID；与 `continue`/`resume` 同用须配 `forkSession` | 用业务生成的 UUID 便于关联 | 2131 |
| `resumeSessionAt` | `string`（消息 uuid） | 只恢复到该消息为止 | 回滚场景 | 2139 |
| `resumeDropsTurn` | `string` | 与 `resumeSessionAt` 配合校验被丢弃范围；仅 headless 路径有效 | 不设 | 2190 |
| `persistSession` | `boolean`，默认 true | false 则不写 `~/.claude/projects/`，不可 resume | 一次性任务 `false`；需多轮则 `true` | 1855 |
| `sessionStore` / `sessionStoreFlush` | `SessionStore` / 策略（@alpha） | 把会话记录镜像到外部存储（仍需本地写入） | 不设 | 1867 / 1875 |
| `enableFileCheckpointing` | `boolean` | 文件检查点，配合 `Query.rewindFiles()` | 需要"撤销模型改动"时 `true` | 1731 |

### 1.7 输出流

| 字段 | 类型 | 含义 | 推荐值 | 行号 |
| --- | --- | --- | --- | --- |
| `includePartialMessages` | `boolean` | 发出 `stream_event`（`SDKPartialAssistantMessage`） | 前端打字机效果 `true` | 1900 |
| `promptSuggestions` | `boolean` | 每轮 result 后多发一条 `prompt_suggestion` | 不设 | 2110 |

### 1.8 用户清单中出现、但本版本不存在/已改名的名称

| 名称 | 结论 | 替代 | 依据 |
| --- | --- | --- | --- |
| `appendSystemPrompt`（顶层 Option） | 『本版本不存在』（只是 initialize 请求的内部字段） | `systemPrompt: {type:'preset', preset:'claude_code', append}` | `rg 'appendSystemPrompt' SDK/sdk.d.ts` 无顶层字段；`SDK/sdk.mjs:163` 内部映射 |
| `continueSession` / `continue_conversation` | 『本版本不存在』（Python 命名） | `continue: true` | `SDK/sdk.d.ts:1631` |
| `apiKey` / `baseUrl` / `authToken` Option | 『本版本不存在』 | 通过 `env` 传 `ANTHROPIC_API_KEY` / `ANTHROPIC_AUTH_TOKEN` / `ANTHROPIC_BASE_URL` | `SDK/sdk.d.ts:1698`；第 4 节 |
| `maxThinkingTokens` | 存在但『已弃用』 | `thinking` | `SDK/sdk.d.ts:1955-1960` |
| `allowedTools: ['Skill']` / `AgentDefinition.tools: ['Skill']` | 存在但『已弃用』 | `skills` 字段 | `SDK/sdk.d.ts:1619`、`:45` |
| `Query.setMaxThinkingTokens` | 存在但『已弃用』 | `thinking` 选项 | `SDK/sdk.d.ts:2992` |
| `structuredOutput` / `jsonSchema` Option | 『本版本不存在』 | `outputFormat` | `SDK/sdk.d.ts:2013` |
| `telemetry` / `disableTelemetry` Option | 『本版本不存在』 | 环境变量 `DISABLE_TELEMETRY` 等 | 第 4.3 节 |
| 环境变量 `DISABLE_NON_ESSENTIAL_MODEL_CALLS` | 『本版本不存在』（二进制与 `DOC/env-vars.md` 均 0 命中） | `CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=1` | `rg -a 'DISABLE_NON_ESSENTIAL_MODEL_CALLS' BIN/claude` 无结果 |
| 环境变量 `ANTHROPIC_SMALL_FAST_MODEL` | 存在（`/tmp/claude_envlist2.txt:3`）但文档标为弃用 | `ANTHROPIC_DEFAULT_HAIKU_MODEL` | `https://code.claude.com/docs/en/model-config` |

### 1.9 `AgentDefinition`（`agents` 的值）字段表

来源：`SDK/sdk.d.ts:38-108`；文档 `https://code.claude.com/docs/en/agent-sdk/subagents`。

| 字段 | 类型 | 含义 | 建议 |
| --- | --- | --- | --- |
| `description` | `string`（必填） | 何时委派给它——主代理据此自动选择 | 写清触发条件 |
| `prompt` | `string`（必填） | 子代理系统提示词 | 子代理是全新上下文，只看到 Agent 工具的委派文本 + 本提示词（+ CLAUDE.md，除非 `omitClaudeMd`） |
| `tools` | `string[]` | 允许工具；省略则继承父级全部；`'Skill'` 已弃用 | 最小集合，如 `['Read','Grep','Glob']` |
| `disallowedTools` | `string[]` | 显式禁用；`mcp__server`/`mcp__server__*`/`mcp__*` 可整服务器移除 | |
| `model` | `string` | 别名 `'sonnet'/'opus'/'haiku'/'fable'`、完整 ID 或 `'inherit'`；省略 → 默认子代理模型（`CLAUDE_CODE_SUBAGENT_MODEL`）否则主模型 | 网关下填网关模型 ID 或 `'inherit'` |
| `mcpServers` | `AgentMcpServerSpec[]` | 子代理专属 MCP | |
| `skills` | `string[]` | 预加载技能内容到子代理上下文 | |
| `initialPrompt` | `string` | 作为主线程 agent 时自动提交的首条 user 消息 | |
| `maxTurns` | `number` | 子代理轮次上限 | 必设 |
| `autoCompactWindow` | `number` | 子代理自动压缩阈值（仅能调低） | |
| `background` | `boolean` | 后台运行 | |
| `omitClaudeMd` | `boolean` | 不加载 user/project/local CLAUDE.md | 委派文本自足时 `true` |
| `memory` | `'user' \| 'project' \| 'local'` | 记忆目录 `~/.claude/agent-memory/<type>/`、`.claude/agent-memory/<type>/`、`.claude/agent-memory-local/<type>/` | |
| `effort` | 等级或整数 | 子代理努力程度 | |
| `permissionMode` | `PermissionMode` | 子代理权限模式（文档：对子代理从不应用 `bypassPermissions`） | |
| `observer` / `observerMessage` | `string` | 后台观察者代理 | 不设 |
| `criticalSystemReminder_EXPERIMENTAL` | `string` | 实验性 | 不设 |

### 1.10 `hooks` 的事件名与回调签名

来源：`SDK/sdk.d.ts:977`（`HOOK_EVENTS`）、`:982`（`HookCallback`）、`:989`（`HookCallbackMatcher`）、`:998`（`HookInput` 联合）、`:1000-1002`（输出与决策）、`:175-186`（`BaseHookInput`）、`:2804-2819`（PreToolUse）、`:2721-2740`（PostToolUse）、`:137-140`（异步输出）、`:9977-9990`（同步输出）；文档 `https://code.claude.com/docs/en/agent-sdk/hooks`、`https://code.claude.com/docs/en/hooks`。

全部 33 个事件名（`HOOK_EVENTS`）：`PreToolUse, PostToolUse, PostToolUseFailure, PostToolBatch, Notification, UserPromptSubmit, UserPromptExpansion, SessionStart, SessionEnd, Stop, StopFailure, SubagentStart, SubagentStop, PreCompact, PostCompact, PreModelSwitch, PostModelSwitch, PermissionRequest, PermissionDenied, Setup, TeammateIdle, TaskCreated, TaskCompleted, Elicitation, ElicitationResult, ConfigChange, WorktreeCreate, WorktreeRemove, InstructionsLoaded, CwdChanged, FileChanged, DirectoryAdded, MessageDisplay`。

```ts
// SDK/sdk.d.ts:982-993
type HookCallback = (input: HookInput, toolUseID: string | undefined,
                     options: { signal: AbortSignal }) => Promise<HookJSONOutput>;
interface HookCallbackMatcher { matcher?: string; hooks: HookCallback[]; timeout?: number }
// Options.hooks: Partial<Record<HookEvent, HookCallbackMatcher[]>>   (sdk.d.ts:1768)

// 输入公共字段 (sdk.d.ts:175-186)
type BaseHookInput = { session_id; transcript_path; cwd; prompt_id?; permission_mode?;
                       agent_id?; agent_type?; effort?: { level: string } };
// PreToolUse 输入 (2804): hook_event_name:'PreToolUse', tool_name, tool_input, tool_use_id, mcp_server?
// PostToolUse 输入 (2721): + tool_response, duration_ms?
// 输出 (9977-9990 / 137-140)
type SyncHookJSONOutput = { continue?; suppressOutput?; stopReason?; decision?: 'approve'|'block';
  systemMessage?; reason?; hookSpecificOutput?: PreToolUseHookSpecificOutput | ... };
type PreToolUseHookSpecificOutput = { hookEventName:'PreToolUse'; permissionDecision?: 'allow'|'deny'|'ask'|'defer';
  permissionDecisionReason?; updatedInput?; additionalContext? };                        // 2812
type PostToolUseHookSpecificOutput = { hookEventName:'PostToolUse'; additionalContext?; classifierContext?;
  updatedToolOutput?; updatedMCPToolOutput? };                                            // 2734
type AsyncHookJSONOutput = { async: true; asyncTimeout?: number };                        // 137
```

要点：
- `matcher` 是对工具名的正则（PreToolUse/PostToolUse），对 SubagentStart/Stop 匹配 agent 类型名（`DOC/sub-agents.md:771-804`）。Windows 上 shell 工具可能是 `PowerShell`，matcher 写 `Bash|PowerShell`（第 4.7 节）。
- 决策优先级 deny > defer > ask > allow（文档 hooks 页）。程序化 hooks 与 settings.json 的 hooks 都会在子代理内触发。
- 必须返回对象（至少 `{}` 或 `{continue:true}`），否则 CLI 等待直到 `timeout`。

### 1.11 推荐配置模板（服务端 + 第三方网关）

```ts
const options: Options = {
  cwd: workspaceDir,
  model: process.env.AGENT_MODEL!,                   // 与网关模型 ID 一致
  systemPrompt: { type: 'preset', preset: 'claude_code', append: BUSINESS_RULES, snapshot: true },
  settingSources: ['project'],                        // 或 [] 完全隔离
  tools: ['Read', 'Write', 'Edit', 'Bash', 'Grep', 'Glob', 'Agent', 'Skill'],
  allowedTools: ['Read', 'Write', 'Edit', 'Grep', 'Glob', 'Agent', 'Bash(npm test:*)', 'mcp__biz__*'],
  disallowedTools: ['WebSearch', 'AskUserQuestion'],
  permissionMode: 'dontAsk',
  permissionPrompts: 'none',
  maxTurns: 30,
  thinking: { type: 'disabled' },                     // 网关不支持 thinking 时
  includePartialMessages: true,
  persistSession: true,
  strictMcpConfig: true,
  env: { ...process.env, ANTHROPIC_BASE_URL: GATEWAY, ANTHROPIC_AUTH_TOKEN: TOKEN,
         ANTHROPIC_MODEL: process.env.AGENT_MODEL!, CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC: '1',
         DISABLE_TELEMETRY: '1', DISABLE_ERROR_REPORTING: '1', CLAUDE_CONFIG_DIR: agentHomeDir },
  stderr: (d) => log.warn(d),
};
```

---

## 2. 消息流：`SDKMessage` 每个成员的字段与用法

来源：`SDK/sdk.d.ts:5579`（`SDKMessage` 联合，39 个成员）及各成员声明行（下表），`T/smoke_glm.log`、`T/smoke_deepseek.log`、`T/smoke_kimi.log`，`JSONL`；文档 `https://code.claude.com/docs/en/agent-sdk/typescript`、`https://code.claude.com/docs/en/agent-sdk/streaming-output`、`https://code.claude.com/docs/en/agent-sdk/streaming-vs-single-mode`。

### 2.1 联合成员速查（判别键 `type` / `subtype`）

| TS 类型（声明行） | `type` | `subtype` | 何时出现 / 关键字段 |
| --- | --- | --- | --- |
| `SDKSystemMessage`（6192） | `system` | `init` | 会话开始第一条；见 2.2 |
| `SDKAssistantMessage`（3688） | `assistant` | — | 模型一条完整消息（含 text / thinking / tool_use 块）；见 2.3 |
| `SDKUserMessage`（6533） | `user` | — | 工具结果（`tool_result` 块）或回显的用户消息；见 2.3 |
| `SDKUserMessageReplay`（6630） | `user` | — | resume 时回放的历史用户消息（带 `uuid`） |
| `SDKPartialAssistantMessage`（5739） | `stream_event` | — | `includePartialMessages:true` 时的原始流事件；见 2.4 |
| `SDKResultSuccess` / `SDKResultError`（5997 / 5933） | `result` | `success` / `error_during_execution` / `error_max_turns` / `error_max_budget_usd` / `error_max_structured_output_retries` | 每轮恰好一条，轮次结束信号；见 2.5 |
| `SDKCompactBoundaryMessage`（3836） | `system` | `compact_boundary` | 上下文压缩发生时；`compact_metadata{trigger:'manual'\|'auto', pre_tokens, post_tokens?, duration_ms?, preserved_segment?, preserved_messages?}` |
| `SDKStatusMessage`（6175） | `system` | `status` | `status: SDKStatus`、`permissionMode?`、`compact_result?: 'success'\|'failed'`、`compact_error?` |
| `SDKAPIRetryMessage`（3666） | `system` | `api_retry` | API 重试通知 |
| `SDKThinkingTokensMessage`（6415） | `system` | `thinking_tokens` | `estimated_tokens`、`estimated_tokens_delta`：思考进度估计（smoke 日志中大量出现） |
| `SDKToolProgressMessage`（6428） | `tool_progress` | — | `tool_use_id, tool_name, parent_tool_use_id, elapsed_time_seconds, task_id?, heartbeat?, subagent_type?, subagent_retry?` |
| `SDKTaskStartedMessage` / `SDKTaskUpdatedMessage` / `SDKTaskProgressMessage` / `SDKTaskNotificationMessage`（6344 / 6389 / 6315 / 6279） | `system` | `task_started` / `task_updated` / `task_progress` / `task_notification` | 子代理与后台任务生命周期；`task_progress` 含 `task_id, run_id?, tool_use_id?, description, subagent_type?, usage{total_tokens, tool_uses, duration_ms}, last_tool_name?, summary?` |
| `SDKBackgroundTasksChangedMessage`（3791） | `system` | `background_tasks_changed` | 后台任务列表变化 |
| `SDKPermissionDeniedMessage`（5776） | `system` | `permission_denied` | 无提示即被拒的工具调用：`tool_name, tool_use_id, agent_id?, decision_reason_type?, decision_reason?`；权威记录仍是 `result.permission_denials` |
| `SDKHookStartedMessage` / `SDKHookProgressMessage` / `SDKHookResponseMessage`（5471 / 5443 / 5456） | `system` | `hook_started` / `hook_progress` / `hook_response` | 需 `includeHookEvents:true`（SessionStart/Setup 除外） |
| `SDKControlRequestProgressMessage`（5239） | `system` | `control_request_progress` | 控制请求进度 |
| `SDKModelRefusalFallbackMessage` / `SDKModelRefusalNoFallbackMessage`（5671 / 5708） | `system` | `model_refusal_fallback` / `model_refusal_no_fallback` | 模型拒答回退 |
| `SDKLocalCommandOutputMessage`（5519） | `system` | `local_command_output` | 本地斜杠命令（如 `/compact`）的文本输出 |
| `SDKPluginInstallMessage`（5870） | `system` | `plugin_install` | 插件安装 |
| `SDKAuthStatusMessage`（3779） | `auth_status` | — | 认证状态 |
| `SDKSessionStateChangedMessage`（6142） | `system` | `session_state_changed` | 会话状态变化 |
| `SDKWorkerShuttingDownMessage`（6711） | `system` | `worker_shutting_down` | 远程 worker 关闭 |
| `SDKCommandsChangedMessage`（3828） | `system` | `commands_changed` | 斜杠命令集合变化（如 reloadSkills 后） |
| `SDKNotificationMessage`（5724） | `system` | `notification` | 通知 |
| `SDKFilesPersistedEvent`（5403） | `system` | `files_persisted` | 文件持久化事件 |
| `SDKToolUseSummaryMessage`（6449） | `tool_use_summary` | — | `summary, preceding_tool_use_ids` |
| `SDKMemoryRecallMessage`（5554） | `system` | `memory_recall` | 自动记忆召回 |
| `SDKRateLimitEvent`（5893） | `rate_limit_event` | — | 速率限制 |
| `SDKElicitationCompleteMessage`（5394） | `system` | `elicitation_complete` | MCP elicitation 完成 |
| `SDKPromptSuggestionMessage`（5883） | `prompt_suggestion` | — | 需 `promptSuggestions:true`，在 `result` 之后到达 |
| `SDKMirrorErrorMessage`（5655） | `system` | `mirror_error` | `sessionStore` 镜像写入错误 |
| `SDKInformationalMessage`（5484） | `system` | `informational` | 信息性提示 |
| `SDKConversationResetMessage`（5366） | `conversation_reset` | — | `/clear` 后 |

处理建议：以 `type` 为第一层 switch，`system` 再按 `subtype`；对未知 `subtype` 一律忽略而不是报错（联合在版本间持续增长，见第 6 节）。

### 2.2 `system/init`（`SDKSystemMessage`，`SDK/sdk.d.ts:6192-6277`）

| 字段 | 含义 |
| --- | --- |
| `agents?: string[]` | 可用子代理名（内置 + 文件系统 + 程序化）。smoke 实测：`['claude','Explore','general-purpose','Plan','statusline-setup']`（`T/smoke_glm.log:1`） |
| `apiKeySource` | `'ANTHROPIC_API_KEY'` / `'apiKeyHelper'` / `'/login managed key'` / `'none'`（Bearer token、OAuth、第三方云均为 `none`） |
| `betas?` | 生效的 beta 列表 |
| `claude_code_version` | 内置 CLI 版本，实测 `2.1.310` |
| `cwd` | 工作目录 |
| `tools: string[]` | 工具名列表。实测 22 个；**列表里子代理工具显示为 `Task`**，但模型实际发出的 `tool_use.name` 是 `Agent`（`SDK/sdk.d.ts:1103` + `SDK/sdk.mjs:126` 的 `LEGACY_TOOL_NAME_ALIASES = {Task:"Agent", KillShell:"TaskStop", KillBash:"TaskStop", ListPeers:"ListAgents", Brief:"SendUserMessage", ListMcpResources:"ListMcpResourcesTool", ReadMcpResource:"ReadMcpResourceTool", ReadMcpResourceDir:"ReadMcpResourceDirTool"}`）。写判断时两个名字都要认 |
| `mcp_servers: {name, status, source?}[]` | MCP 服务器及连接状态（`.mcp.json`/插件来源的常为 `pending`） |
| `model` | 实际模型 ID（网关下即你传的字符串，如 `glm-5.3`） |
| `permissionMode` | 生效的权限模式 |
| `slash_commands: string[]` | 可用斜杠命令（内置 + 捆绑技能 + 用户技能 + `.claude/commands`）；实测 44 个 |
| `terminal_slash_commands?` | 仅终端可用的命令 |
| `output_style` | 输出风格名 |
| `skills: string[]` | 用户可调用技能；实测 16 个**捆绑技能**：`deep-research, dataviz, update-config, verify, debug, code-review, simplify, batch, fewer-permission-prompts, doctor, loop, claude-api, workflow-authoring, run, run-skill-generator, plugin-authoring`（`T/smoke_glm.log:1`） |
| `plugins[]` / `plugin_errors?` | 插件及错误 |
| `fast_mode_state?` / `fast_mode_disabled_reason?` / `effort?` / `view_mode?` / `capabilities?` | 其他运行态 |
| `uuid` / `session_id` | 消息 ID / **会话 ID（保存它用于 resume）** |

### 2.3 `assistant` 与 `user`

`SDKAssistantMessage`（`SDK/sdk.d.ts:3688-3732`）：

| 字段 | 含义 |
| --- | --- |
| `message: BetaMessage` | Anthropic Messages API 的 assistant 消息：`id, model, content[], stop_reason, usage`。`content` 块类型：`text`、`thinking`（含 `signature`）、`tool_use`（`id, name, input`）、`redacted_thinking`。流式时一条 assistant 消息通常只含一个块，`stop_reason` 为 null |
| `parent_tool_use_id: string \| null` | null = 主线程；非 null = 产生它的子代理所对应的 `Agent` 工具调用 id（2.7） |
| `agent_id?` / `subagent_type?` | 子代理 ID（= 其 task 事件的 `task_id`）与类型；主线程消息没有 |
| `error?` | API 错误时的结构化错误 |
| `aborted?: true` | 被中断截断 |
| `uuid`、`session_id`、`request_id?` | 标识 |

`SDKUserMessage`（`SDK/sdk.d.ts:6533-6608`）：

| 字段 | 含义 |
| --- | --- |
| `message: MessageParam` | `role:'user'`，`content` 为字符串或块数组；工具结果为 `{type:'tool_result', tool_use_id, content, is_error?}` |
| `parent_tool_use_id` | 同上 |
| `tool_use_result?: unknown` | 工具的结构化完整输出（非发给模型的文本）；Agent 工具完成时是子代理最终报告 + 用量 |
| `isSynthetic?` | CLI 合成的用户消息 |
| `uuid?` / `session_id?` | 作为**输入**（流式输入）时可省略，SDK 会补 |
| `priority?`、`origin?`、`shouldQuery?`、`timestamp?`、`client_composed?: true`、`agent_id?`、`subagent_type?` | 高级字段；`client_composed` 单条禁用 `@path` 展开与斜杠分发 |

实测顺序（`T/smoke_glm.log:133-136`）：`assistant(tool_use Write)` → `user(tool_result)` → `assistant(tool_use Bash {"command":"python3 -c \"print(1+1)\""})` → `user(tool_result "2")`。

### 2.4 `stream_event`（`SDKPartialAssistantMessage`，`SDK/sdk.d.ts:5739-5765`）

`event: BetaRawMessageStreamEvent`，即 Messages API 流事件：`message_start` → `content_block_start` → `content_block_delta`* → `content_block_stop` → `message_delta` → `message_stop`。`content_block_delta.delta.type`：

| delta.type | 字段 | 用途 |
| --- | --- | --- |
| `text_delta` | `text` | 正文增量 → 打字机输出 |
| `thinking_delta` | `thinking` | 推理增量（可选展示） |
| `input_json_delta` | `partial_json` | `tool_use.input` 的 JSON 片段，需按 `content_block_start` 的 `index` 累加后 `JSON.parse` |
| `signature_delta` | `signature` | thinking 签名，忽略 |

还带 `parent_tool_use_id`（子代理流事件同样会来）、`uuid`、`session_id`、`ttft_ms?`。文档的处理模式（`https://code.claude.com/docs/en/agent-sdk/streaming-output`）：`if (m.type==='stream_event' && m.event.type==='content_block_delta' && m.event.delta.type==='text_delta') write(m.event.delta.text)`。注意：开启后 `assistant` 整条消息仍会到达，不要把 `text_delta` 与 `assistant.message.content[].text` 重复渲染。

### 2.5 `result`（`SDK/sdk.d.ts:5933-6090`）

| 字段 | 含义与注意 |
| --- | --- |
| `subtype` | `success` 或四种错误；`terminal_reason?`（`SDK/sdk.d.ts:10048`）给出更细原因：`max_turns, budget_exhausted, aborted_streaming, aborted_tools, stop_hook_prevented, hook_stopped, api_error, prompt_too_long, completed, ...` |
| `result: string`（success） | 最终文本；`is_error:true` 时为错误文本（API 错误仍可能是 `success` 子类型） |
| `errors: string[]`（error） | 错误文本；`startup_failure_reason?: SDKStartupFailureReason`（`SDK/sdk.d.ts:6171`，含 `shell_tool_missing`、`cwd_unavailable`、`proxy_invalid`、`cli_version_too_old` 等） |
| `num_turns` | 本轮内模型往返次数（smoke: 3） |
| `duration_ms` / `duration_api_ms` / `ttft_ms?` | 耗时 |
| `total_cost_usd` | **累计**（流式输入会话中跨轮累加；每条 result 是到目前为止的总额，读最新一条而不是求和）；resume 时从会话保存的总额继续；`/clear` 清零；估算值 |
| `usage` | **仅主循环、每轮**；不含子代理/侧链；字段同 Messages API（`input_tokens, output_tokens, cache_read_input_tokens, cache_creation_input_tokens, server_tool_use{web_search_requests, web_fetch_requests}` 等，`T/smoke_glm.log:215-224`）。d.ts 明确"Prefer modelUsage" |
| `modelUsage: Record<modelId, ModelUsage>` | 按模型累计：`inputTokens, outputTokens, thinkingTokens?, cacheReadInputTokens, cacheCreationInputTokens, webSearchRequests, costUSD, contextWindow, maxOutputTokens, canonicalModel?, provider?, costBasis?: 'list'\|'managed'\|'unknown'`（`SDK/sdk.d.ts:1469`）。第三方模型 ID → `costBasis:'unknown'`，会话记录里 `hasUnknownModelCost:true`（`JSONL:34`）——金额不可信，自己按 token 计价 |
| `permission_denials: SDKPermissionDenial[]` | `{tool_name, tool_use_id, tool_input}`（`SDK/sdk.d.ts:5767`） |
| `structured_output?` | `outputFormat` 的校验后结果（`SDK/sdk.d.ts:6078`） |
| `stop_reason`、`api_error_status?`、`queued_turn_count?`、`result_index?`、`user_message_uuid?` | 其他 |
| `uuid` / `session_id` | 会话 ID |

实测：`T/smoke_glm.log:215` `RESULT success {turns: 3, cost: 0.0555632, ...}`；`T/smoke_deepseek.log:199`（0.0575224）；`T/smoke_kimi.log:313`（0.089788）。文档 `https://code.claude.com/docs/en/agent-sdk/cost-tracking`：多步工具调用时每步 assistant 消息的 `usage.output_tokens` 可能是占位值，最终以 result 为准。

### 2.6 子代理输出识别：`parent_tool_use_id` 与 `agent_id`

来源：`SDK/sdk.d.ts:3694-3732`、`:1901-1907`、`:6428-6447`；文档 `https://code.claude.com/docs/en/agent-sdk/subagents`。

1. 主代理发出 `assistant.message.content[]` 中 `type:'tool_use'`、`name:'Agent'`（旧版/init 列表为 `Task`）、`input:{subagent_type, description, prompt, ...}` 的块 → 记下 `tool_use.id`。
2. 其后凡 `parent_tool_use_id === 该 id` 的 `assistant` / `user` / `stream_event` / `tool_progress` 都是该子代理内部活动；默认只转发子代理的 tool_use/tool_result，要正文需 `forwardSubagentText:true`。
3. `assistant.agent_id` 与 `task_started/task_progress.task_id` 相同，可把子代理消息与任务进度行关联；`task_progress.subagent_type` 给类型名。
4. 子代理结束：主线程收到 `user` 消息，其 `tool_result.tool_use_id === 该 id`，`tool_use_result` 是结构化最终报告。
5. 子代理默认后台化：Agent 工具可能立刻返回 `status:"async_launched"`，再由 `task_notification` 送达结果（文档 subagents 页）。深度上限 `CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH`（默认 3），并发上限 `CLAUDE_CODE_MAX_CONCURRENT_SUBAGENTS`（默认 20）。
6. 分叉子代理 `subagent_type:"fork"`（继承父上下文），可用 `disallowedTools:['Agent(fork)']` 禁止。

### 2.7 `uuid` 与 `session_id`

- 每条消息都有 `uuid`（消息级）与 `session_id`（会话级）。`session_id` 在 `init` 与每条 `result` 上都有，保存后用于 `resume`。
- `SDKAssistantMessage.uuid` 可作为 `resumeSessionAt` 的锚点（`SDK/sdk.d.ts:2133-2139`）；`enableFileCheckpointing` 时用用户消息 uuid 调 `rewindFiles()`。
- 会话文件：`<CLAUDE_CONFIG_DIR 或 ~/.claude>/projects/<cwd 路径把 / 替换为 ->/<session_id>.jsonl`（实测 `T/home/.claude/projects/-tmp-sdktest-ws/42e923ab-....jsonl`，因 `HOME` 被指到 `T/home`）。`persistSession:false` 不写文件。

### 2.8 中断与多轮

- 单次 prompt 为字符串时，一轮结束进程退出，`interrupt()`/`setModel()` 等控制方法**不可用**（`SDK/sdk.d.ts:2915` 注释："only supported when streaming input/output is used"）。
- 流式输入：`prompt` 传 `AsyncIterable<SDKUserMessage>`（消息至少 `{type:'user', message:{role:'user', content:'...'}, parent_tool_use_id:null}`）。生成器不结束则会话保持；每个用户消息对应一条 `result`，以 `result` 作为"该轮完成"的信号再发下一条。也可用 `q.streamInput(iterable)` 追加（`SDK/sdk.d.ts:3290`）。
- `q.interrupt()` 中断当前轮（`SDK/sdk.d.ts:2929`）；`abortController.abort()` 终止整个进程。
- `/compact`、`/clear` 作为 prompt 文本发送即可（`https://code.claude.com/docs/en/agent-sdk/slash-commands`）；`/clear` 后会收到 `conversation_reset`，成本累计清零。

### 2.9 `Query` 控制方法（`SDK/sdk.d.ts:2915-3320`，流式输入模式）

| 方法 | 行号 | 用途 |
| --- | --- | --- |
| `interrupt()` | 2929 | 中断当前轮 |
| `setPermissionMode(mode)` | 2936 | 运行中切换权限模式 |
| `setMcpPermissionModeOverride(server, mode)` | 2953 | 某 MCP 服务器的模式覆盖 |
| `setModel(model?)` | 2965 | 切换模型 |
| `setMaxThinkingTokens(...)` | 2992 | 『已弃用』 |
| `applyFlagSettings(settings)` / `updateSettings(source, settings)` | 3020 / 3035 | 运行中改设置 |
| `initializationResult()` / `reinitialize()` | 3044 / 3070 | 初始化结果（含 `plugins_applied`） |
| `supportedCommands()` / `supportedModels()` / `supportedAgents()` | 3076 / 3082 / 3088 | 列表查询（`SlashCommand`、`ModelInfo`、`AgentInfo`） |
| `mcpServerStatus()` | 3094 | MCP 状态 |
| `getContextUsage()` | 3105 | 上下文占用 |
| `readFile(path)` | 3147 | 让子进程读文件 |
| `reloadPlugins()` / `reloadSkills()` / `reloadOutputStyles()` | 3170 / 3179 / 3189 | 热重载 |
| `accountInfo()` | 3195 | 账户信息 |
| `rewindFiles(userMessageId)` | 3204 | 文件回滚（需 `enableFileCheckpointing`） |
| `seedReadState(path, mtime)` | 3217 | 预置已读状态 |
| `reconnectMcpServer()` / `toggleMcpServer()` / `readMcpResource()` / `setMcpServers()` | 3233 / 3241 / 3254 / 3283 | MCP 管理（`setMcpServers` 可运行中增删服务器） |
| `streamInput(iterable)` | 3290 | 追加用户消息流 |
| `stopTask(taskId)` / `backgroundTasks(toolUseId?)` | 3295 / 3310 | 后台任务控制 |
| `close()` | 3319 | 关闭 |

### 2.10 一次真实消息流（`T/smoke_glm.log`）

```
[0.2s]  system init {model:'glm-5.3', tools:22, agents:[claude,Explore,general-purpose,Plan,statusline-setup], slash:44, skills:16}   # 行 1
        system thinking_tokens × N ；assistant(thinking)
[15.4s] assistant tool_use: Write {"file_path":"/tmp/sdktest/ws/hello.md",...}                                                        # 行 133
        user tool_result                                                                                                                # 行 134
[15.9s] assistant tool_use: Bash {"command":"python3 -c \"print(1+1)\""}                                                               # 行 135
        user tool_result: 2                                                                                                             # 行 136
        assistant text
[20.5s] result success {turns:3, cost:0.0555632, usage:{input_tokens:10735, cache_read_input_tokens:9216, output_tokens:539,
        server_tool_use:{web_search_requests:0, web_fetch_requests:0}}}                                                                  # 行 215-224
```

对应的 API 请求形状（`JSONL:14,17`）：system 为极简一句话；首轮工具 `Agent, Bash, Edit, ListAgents, Read, ReportFindings, ScheduleWakeup, Skill, ToolSearch, Workflow, DeferredToolPlaceholder, Write`；`WebFetch/WebSearch/NotebookEdit/Cron*/EnterWorktree/ExitWorktree/Monitor/SendMessage/TaskStop` 为"延迟工具"，只有模型先调用 `ToolSearch` 才会加载其 schema（`JSONL:6` 的 `deferred_tools_delta`）；请求带 `betas:[claude-code-20250219, interleaved-thinking-2025-05-14, ...]`、`thinking:{type:"adaptive"}`、`output_config:{effort:"high"}`、`max_tokens:32000`——这些在第三方网关上可能被拒绝，见第 4 节。

---

## 3. 项目级配置：CLAUDE.md / agents / skills / commands / settings.json

来源：`SDK/sdk.d.ts:2306-2316`（`settingSources`）、`SDK/sdk.mjs:157`（默认值 `lYe=["user","project","local"]` 与 `--setting-sources` 映射）、`SDK/sdk.d.ts:6995-9665`（`Settings` 接口）；文档 `https://code.claude.com/docs/en/agent-sdk/claude-code-features`（`DOC/agent-sdk_claude-code-features.md:74-114, 268-291`）、`https://code.claude.com/docs/en/memory`（`DOC/memory.md:57-64, 102-176`）、`https://code.claude.com/docs/en/sub-agents`（`DOC/sub-agents.md:298-316, 600-605, 652-654`）、`https://code.claude.com/docs/en/skills`（`DOC/skills.md:17, 158-211, 399-418`）、`https://code.claude.com/docs/en/agent-sdk/skills`（`DOC/agent-sdk_skills.md:17-48, 108`）、`https://code.claude.com/docs/en/settings`（`DOC/settings.md:25-30, 206-210`）、`https://code.claude.com/docs/en/settings-reference`、`https://code.claude.com/docs/en/agent-sdk/slash-commands`、`https://code.claude.com/docs/en/claude-directory`。

### 3.1 `settingSources` 决定一切：加载规则总表

| 值 | 加载内容 | 位置 |
| --- | --- | --- |
| `'project'` | 项目 `settings.json` 与其中 hooks/permissions；项目 CLAUDE.md 与 `.claude/rules/*.md`；项目 skills、commands、subagents；`.mcp.json` | `<cwd>/.claude/`；CLAUDE.md 取 `<cwd>` 及所有父目录 |
| `'user'` | 用户 `settings.json`；`~/.claude/CLAUDE.md` 与 `~/.claude/rules/*.md`；`~/.claude/skills/`、`~/.claude/commands/`、`~/.claude/agents/` | `~/.claude/`（Windows 为 `%USERPROFILE%\.claude`；可用 `CLAUDE_CONFIG_DIR` 整体改位） |
| `'local'` | `CLAUDE.local.md`、`.claude/settings.local.json` | `<cwd>/.claude/`；`CLAUDE.local.md` 取 `<cwd>` 及父目录 |
| 省略 | **三者全部加载**（`SDK/sdk.d.ts:2312`、`SDK/sdk.mjs:157`、文档："Omitting settingSources is equivalent to ["user","project","local"]"） | — |
| `[]` | 隔离模式：不读任何文件系统设置；CLAUDE.md、项目 skills/agents/commands、settings hooks 与权限规则全部不加载 | — |

**不受 `settingSources` 控制、始终生效的输入**（`DOC/agent-sdk_claude-code-features.md:88-94`）：

| 输入 | 行为 | 关闭方法 |
| --- | --- | --- |
| 托管策略设置（MDM / `managed-settings.json` / 服务器托管） | 总是加载 | 无（管理员控制） |
| `~/.claude.json` 全局配置 | 总是读取 | `env.CLAUDE_CONFIG_DIR` 改位 |
| 自动记忆 `~/.claude/projects/<project>/memory/` | 会话开始注入系统提示 | settings `autoMemoryEnabled:false` 或 `env.CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` |
| claude.ai MCP 连接器 | claude.ai 登录时加载；`mcpServers:{}` 不能屏蔽 | `strictMcpConfig:true` 或 settings `disableClaudeAiConnectors:true` |
| `~/.claude/settings.json` 里 `sandbox.credentials` 的 deny/mask | 沙箱运行时仍应用 | 删条目 |

**服务端建议**：显式写 `settingSources`。只认工作区配置用 `['project']`；完全由代码控制用 `[]` 并把需要的规则放进 `options.settings`（内联对象，等价 `--settings`，优先级最高，`SDK/sdk.d.ts:2247-2264`）。

### 3.2 CLAUDE.md

加载条件：`settingSources` 含 `'project'`（项目/父目录/子目录）、`'user'`（`~/.claude/CLAUDE.md`）、`'local'`（`CLAUDE.local.md`）；`SDK/sdk.d.ts:2315` 明确 "Must include 'project' to load CLAUDE.md files"。

层级与顺序（`DOC/memory.md:57-64, 142-146`；`DOC/agent-sdk_claude-code-features.md:106-114`）：

| 级别 | 位置 | 何时加载 |
| --- | --- | --- |
| 托管策略 | Linux `/etc/claude-code/CLAUDE.md`、macOS `/Library/Application Support/ClaudeCode/CLAUDE.md`、Windows `C:\Program Files\ClaudeCode\CLAUDE.md`；或 settings `claudeMd` 字符串（仅策略层，`SDK/sdk.d.ts:9531`） | 总是 |
| 用户 | `~/.claude/CLAUDE.md`、`~/.claude/rules/*.md` | `'user'` |
| 项目根 | `<cwd>/CLAUDE.md` 或 `<cwd>/.claude/CLAUDE.md`（`AGENTS.md` 可替代/并存） | `'project'` |
| 项目规则 | `<cwd>/.claude/rules/*.md` 及每个父目录的 `.claude/rules/*.md`（可用 frontmatter `paths` 限定文件类型） | `'project'` |
| 父目录 | cwd 之上每一级目录的 `CLAUDE.md` | `'project'`，启动时加载 |
| 子目录 | cwd 之下的 `CLAUDE.md` | `'project'`，模型首次读写该目录文件时按需加载 |
| 本地 | `<cwd>/CLAUDE.local.md` 及父目录 | `'local'` |

规则：所有发现的文件**拼接**而非覆盖，顺序从文件系统根到 cwd；`@path/to/file` 可导入其他文件（`DOC/memory.md:102`）；`settings.claudeMdExcludes` 可按 glob 排除 User/Project/Local 文件（`SDK/sdk.d.ts:9535`）；`additionalDirectories` 中的目录也会加载其 CLAUDE.md/.claude/rules（`DOC/memory.md:164`）。子代理默认也读 CLAUDE.md，`AgentDefinition.omitClaudeMd:true` 可关（`SDK/sdk.d.ts:88`）。

### 3.3 `.claude/agents/*.md`（文件系统子代理）

位置：项目 `.claude/agents/<name>.md`（需 `'project'`）、用户 `~/.claude/agents/<name>.md`（需 `'user'`）；程序化 `agents` 同名者优先（文档 subagents 页）。格式 = YAML frontmatter + 正文（正文即系统提示词）。frontmatter 字段（`DOC/sub-agents.md:298-316`）：

| 字段 | 必填 | 说明 |
| --- | --- | --- |
| `name` | 是 | ≤256 字符的唯一标识 |
| `description` | 是 | 何时委派 |
| `tools` | 否 | 逗号分隔或 YAML 列表；省略继承全部 |
| `disallowedTools` | 否 | 从继承集合中移除；支持 `Bash(git push:*)` 这类带限定符的条目 |
| `model` | 否 | `sonnet`/`opus`/`haiku`/`fable`/完整 ID/`inherit` |
| `permissionMode` | 否 | `default`/`acceptEdits`/`auto`/`dontAsk`/`bypassPermissions`/`plan` |
| `maxTurns`、`autoCompactWindow` | 否 | 轮次上限；自动压缩窗口（100000 起） |
| `skills` | 否 | 启动时预载入技能全文 |
| `mcpServers` | 否 | 子代理专属 MCP（名称引用或内联定义） |
| `hooks` | 否 | 子代理范围的 hooks（`PreToolUse`/`PostToolUse`/`Stop`→运行时转为 `SubagentStop`，`DOC/sub-agents.md:771-773`） |
| `memory` | 否 | `user`/`project`/`local`（目录见 1.9） |
| `background`、`omitClaudeMd`、`effort`、`isolation: worktree`、`color`、`initialPrompt`、`experimental.cacheTtl` | 否 | 见文档 |

与程序化 `AgentDefinition` 的差异：程序化版本没有 `name`（键即名字）、`isolation`、`color`、`hooks`、`experimental`；文件版没有 `observer*`、`criticalSystemReminder_EXPERIMENTAL`（对照 `SDK/sdk.d.ts:38-108`）。`strictMcpConfig:true` 会忽略磁盘 agent frontmatter 里的 MCP（`SDK/sdk.d.ts:2359-2364`）。内置代理（`general-purpose`、`Explore`、`Plan`、`claude`、`statusline-setup`）可用 `CLAUDE_AGENT_SDK_DISABLE_BUILTIN_AGENTS=1` 关闭（文档 subagents 页）。

### 3.4 `.claude/skills/<name>/SKILL.md`

**只能是文件系统工件，SDK 没有程序化定义技能的 API**（`DOC/agent-sdk_skills.md:17-24`）。位置与优先级（`DOC/skills.md:158-162, 211`）：

| 范围 | 位置 | 需要的 settingSource |
| --- | --- | --- |
| 企业 | 托管设置目录下 `.claude/skills/<name>/SKILL.md` | 总是 |
| 个人 | `~/.claude/skills/<name>/SKILL.md` | `'user'` |
| 项目 | `.claude/skills/<name>/SKILL.md`（启动目录到仓库根的每一级都会扫描） | `'project'` |
| 嵌套 | `<subdir>/.claude/skills/...` | `'project'`，模型首次读写该子目录时加载 |
| 附加目录 | `additionalDirectories` 中的 `.claude/skills/` | — |

同名时企业 > 个人 > 项目。保留名 `synced` 不可用。

文件格式：`---` 包裹的 YAML frontmatter + Markdown 正文；目录内可放脚本/模板等支撑文件，正文里用相对路径引用；支持 `$ARGUMENTS`、`$0`/`$1` 位置参数替换与 `` !`command` `` 前置命令（`DOC/skills.md:102-116, 399-418`）。frontmatter 字段：

| 字段 | 说明 |
| --- | --- |
| `name` | `/` 菜单命令名，默认目录名 |
| `description` | **自动触发依据**：模型据此决定何时应用；省略则用正文首段 |
| `argument-hint`、`arguments` | 自动补全提示；命名位置参数 |
| `disable-model-invocation: true` | 禁止模型自动调用，只能用户 `/name` 触发 |
| `user-invocable: false` | 只允许模型调用，不出现在 `/` 菜单 |
| `allowed-tools` / `disallowed-tools` | 技能生效的那一轮免提示工具 / 移除工具 |
| `model`、`effort` | 技能生效期间的模型/努力度覆盖 |
| `context: fork` + `agent` + `background` | 在子代理上下文运行 |
| `hooks` | 技能被调用时注册、会话期间持续的 hooks |
| `paths` | glob，限定只在处理匹配文件时激活 |
| `shell: bash\|powershell` | `` !`command` `` 用的 shell（Windows 相关） |
| `metadata`、`license`、`compatibility` | 元数据（Agent Skills 规范） |

运行机制：启动时只把技能的 **name + description** 放进上下文（渐进加载），模型通过 **`Skill` 工具**按名调用才读入全文；用户也可在 prompt 中直接写 `/name 参数` 分发（不依赖 `skills` 选项，`DOC/agent-sdk_slash-commands.md:164`）。SDK 侧：`skills` 选项省略时"发现的技能全部可用且 Skill 工具可用"；设置后 SDK 自动把 Skill 加入 `allowedTools`，但若你同时传了显式 `tools` 列表，须把 `'Skill'` 写进去（`DOC/agent-sdk_skills.md:32-46`）。`init.skills` 列出用户可调用技能，其中 16 个是 CLI 捆绑技能（2.2 节）。

### 3.5 `.claude/commands/*.md`（旧格式斜杠命令）

"Custom commands have been merged into skills"（`DOC/skills.md:17`）：`.claude/commands/deploy.md` 与 `.claude/skills/deploy/SKILL.md` 都创建 `/deploy`，行为相同；commands 文件支持同一套 frontmatter，**除 `name` 和 `paths`**（`DOC/skills.md:172`），命令名即文件名。项目级需 `'project'`，用户级 `~/.claude/commands/` 需 `'user'`。SDK 中：`init.slash_commands` 混合列出内置命令、捆绑技能、用户技能与 commands 文件；发送方式就是把 `/name args` 作为 prompt 文本（`DOC/agent-sdk_slash-commands.md:116-164`）。常用内置：`/compact`（触发 `compact_boundary`）、`/clear`（触发 `conversation_reset`，成本累计清零）。

### 3.6 `.claude/settings.json`

来源：`SDK/sdk.d.ts:6995-9665`（`Settings` 接口，所有键均有 JSDoc）、`DOC/settings.md:25-30, 206-210`、`https://code.claude.com/docs/en/settings-reference`。

优先级（高→低）：托管设置 > 命令行/`options.settings`（flag settings）> `.claude/settings.local.json` > `.claude/settings.json` > `~/.claude/settings.json`。

与 SDK 实现者最相关的键：

| 键 | 类型/取值 | 说明 | d.ts 行 |
| --- | --- | --- | --- |
| `env` | `Record<string,string>` | 注入会话环境变量（如 `ANTHROPIC_BASE_URL`）；与 `options.env` 叠加 | 7080 |
| `permissions.allow` / `deny` / `ask` | `string[]` 规则，如 `"Bash(npm run test:*)"`、`"Read(./.env)"`、`"mcp__srv__*"`、`"WebFetch(domain:example.com)"` | 权限规则；评估顺序 hooks → deny → ask → mode → allow → canUseTool（文档 permissions 页） | 7125 |
| `permissions.defaultMode` | 6 种模式之一 | 省略 `options.permissionMode` 时的起始模式 | 7129 |
| `permissions.disableBypassPermissionsMode` | `'disable'` | 禁止 bypass | 7130 |
| `permissions.additionalDirectories` / `blockReadsOutsideWorkingDirectories` | | 目录范围 | 7131-7132 |
| `model` / `fallbackModel` | `string` / `string[]` | 默认模型 | 7159 / 7163 |
| `hooks` | `{ [event]: [{matcher?, hooks:[{type:'command'\|'http'\|'mcp_tool'\|'prompt'\|'agent', command?/url?/..., timeout?, async?}]}] }` | 文件系统 hooks（与程序化 hooks 并存，均在子代理内触发） | 7305-7501 |
| `skipWebFetchPreflight` | `boolean` | 跳过 WebFetch 对 `api.anthropic.com/api/web/domain_info` 的预检（第 4.6 节） | 8968 |
| `sandbox` | 对象 | 命令沙箱 | 8969 |
| `effortLevel` | `'low'\|'medium'\|'high'\|'xhigh'` | 默认努力度 | 9297 |
| `agent` | `string` | 主线程默认 agent | 9366 |
| `claudeMd` / `claudeMdExcludes` | `string` / `string[]` | 托管记忆 / 排除 CLAUDE.md | 9531 / 9535 |
| `enableAllProjectMcpServers`、`enabledMcpjsonServers`、`disabledMcpjsonServers`、`disableClaudeAiConnectors` | | `.mcp.json` 与连接器控制 | （grep `SDK/sdk.d.ts`） |
| `autoMemoryEnabled`、`promptSuggestionEnabled`、`alwaysThinkingEnabled`、`disableAllHooks`、`attribution`（替代已弃用 `includeCoAuthoredBy`）、`apiKeyHelper`、`cleanupPeriodDays`、`outputStyle`、`enabledPlugins` | | 其他常用键 | （grep `SDK/sdk.d.ts`） |

SDK 模式下是否默认加载：**是**（省略 `settingSources` 时三层都读）；`[]` 时全部不读，但 `options.settings` 内联对象仍生效。

### 3.7 `.mcp.json`

项目根的 `.mcp.json` 随 `'project'` 加载（`DOC/agent-sdk_mcp.md:130`）；`strictMcpConfig:true` 忽略它。工具命名 `mcp__<server>__<tool>`；放行 `allowedTools:['mcp__<server>__*']`。`init.mcp_servers[].status` 对文件/插件来源常为 `pending`，用 `q.mcpServerStatus()` 复查。

### 3.8 判定表：SDK 模式下"默认"到底加载什么

| 工件 | `settingSources` 省略 | `['project']` | `[]` |
| --- | --- | --- | --- |
| `~/.claude/settings.json`、`~/.claude/CLAUDE.md`、`~/.claude/{skills,commands,agents}` | 加载 | 不加载 | 不加载 |
| `<cwd>/CLAUDE.md`、`.claude/rules`、父目录 CLAUDE.md | 加载 | 加载 | 不加载 |
| `.claude/settings.json`（hooks、permissions、env） | 加载 | 加载 | 不加载 |
| `.claude/{skills,commands,agents}`、`.mcp.json` | 加载 | 加载 | 不加载 |
| `.claude/settings.local.json`、`CLAUDE.local.md` | 加载 | 不加载 | 不加载 |
| 托管策略、`~/.claude.json`、自动记忆、claude.ai 连接器 | 加载 | 加载 | 加载（需用上表方法单独关） |
| `options.agents` / `options.hooks` / `options.mcpServers` / `options.settings` | 生效 | 生效 | 生效 |
| 16 个捆绑技能、内置代理 | 可用 | 可用 | 可用（`skills:[...]` 过滤；`CLAUDE_AGENT_SDK_DISABLE_BUILTIN_AGENTS=1`） |

---

## 4. 第三方 Anthropic 兼容网关、环境变量、Windows 与二进制定位

来源：`SDK/sdk.d.ts:1682-1700`（`env` 语义）、`SDK/sdk.mjs:157,163`、`BIN/claude`（二进制字符串搜索，命令形如 `rg -a -o '<串>' BIN/claude`，无行号）、`T/smoke.mjs:13-30`、`JSONL:14,17,34`、`/tmp/claude_envlist2.txt`（从二进制提取的 192 个配置变量名）；文档 `https://code.claude.com/docs/en/llm-gateway-connect`（`DOC/llm-gateway-connect.md:60-64, 141-159, 360-386, 505-535`）、`https://code.claude.com/docs/en/llm-gateway-protocol`（`DOC/llm-gateway-protocol.md:43, 101-114, 220-270`）、`https://code.claude.com/docs/en/env-vars`（`DOC/env-vars.md`，行号见各表）、`https://code.claude.com/docs/en/model-config`（`DOC/model-config.md:235-237, 331, 733`）、`https://code.claude.com/docs/en/tools-reference`（`DOC/tools-reference.md:618-652, 670-684`）、`https://code.claude.com/docs/en/setup`（`DOC/setup.md:66, 118-130`）。

### 4.1 凭据与地址

| 变量 | 作用 | 依据 |
| --- | --- | --- |
| `ANTHROPIC_BASE_URL` | 把 `/v1/messages`（及可选 `/v1/messages/count_tokens`）请求指向网关；网关需接受 `anthropic-beta`、`anthropic-version` 头 | `DOC/env-vars.md:138`；`DOC/llm-gateway-protocol.md:43` |
| `ANTHROPIC_AUTH_TOKEN` | 放入 `Authorization: Bearer <token>` | `DOC/env-vars.md:134`；`DOC/llm-gateway-connect.md:60,155` |
| `ANTHROPIC_API_KEY` | 放入 `x-api-key` 头；与 Bearer 二选一（网关没说明就先用 AUTH_TOKEN） | `DOC/env-vars.md:133`；`DOC/llm-gateway-connect.md:61-64` |
| `apiKeyHelper`（settings） | 两个头都带；凭据轮换/从 vault 取 | `DOC/llm-gateway-connect.md:62,155` |
| `ANTHROPIC_CUSTOM_HEADERS` | `Name: Value`（多行分隔）附加头，如租户路由 | `DOC/env-vars.md:144`；`DOC/llm-gateway-connect.md:275-286` |
| `CLAUDE_CODE_OAUTH_TOKEN` | claude.ai 长效 token（非网关场景） | `DOC/env-vars.md:339` |

两个凭据同时存在会有 "auth may not work as expected" 的启动警告（`DOC/llm-gateway-connect.md:509`）。`init.apiKeySource` 在 Bearer 场景为 `'none'`（`SDK/sdk.d.ts:6196-6198`）。

**TS 侧必须注意**：`options.env` 一旦设置就**整体替换**子进程环境（`SDK/sdk.d.ts:1684-1688`），必须 `{ ...process.env, ... }`，否则连 `PATH` 都没有；`T/smoke.mjs:20-29` 正是这样写的。

### 4.2 模型相关环境变量

| 变量 | 作用 | 状态/依据 |
| --- | --- | --- |
| `ANTHROPIC_MODEL` | 默认模型名（与 `options.model` 作用相同；`options.model` 走 `--model` flag，优先） | `DOC/env-vars.md:171`；`/tmp/claude_envlist2.txt:1` |
| `ANTHROPIC_DEFAULT_HAIKU_MODEL` | `haiku` 别名解析目标，**也是后台/辅助调用（标题生成、摘要、WebFetch 提取等）用的模型**；网关下必须设成网关可用的 ID，否则后台调用会打到网关不认识的 `claude-haiku-*` | `DOC/env-vars.md:153`；`DOC/llm-gateway-protocol.md:105` |
| `ANTHROPIC_SMALL_FAST_MODEL` | 旧名，文档标 `[DEPRECATED]`，二进制仍识别 | `DOC/env-vars.md:174`；`/tmp/claude_envlist2.txt:3` |
| `ANTHROPIC_DEFAULT_SONNET_MODEL` / `ANTHROPIC_DEFAULT_OPUS_MODEL` / `ANTHROPIC_DEFAULT_FABLE_MODEL` | 别名解析；**不能**重定向完整模型 ID | `DOC/env-vars.md:158,162`；`DOC/model-config.md:235,331` |
| `CLAUDE_CODE_SUBAGENT_MODEL` | 子代理默认模型（优先级：Agent 调用参数 > frontmatter/`AgentDefinition.model` > 本变量 > 主模型） | `DOC/env-vars.md:396`；`SDK/sdk.d.ts:56` |
| `CLAUDE_CODE_SUBAGENT_MODEL_FORCE=1` | 强制所有子代理用同一模型 | `DOC/env-vars.md:397` |
| `CLAUDE_CODE_MAX_CONTEXT_TOKENS` | 覆盖未知模型 ID 的上下文窗口假设（默认 200K，ID 带 `[1m]` 为 1M） | `DOC/env-vars.md:318`；`DOC/llm-gateway-protocol.md:113` |
| `CLAUDE_CODE_DISABLE_THINKING=1` | 请求中完全不带 `thinking` 字段（网关兼容选项） | `DOC/env-vars.md:272`；`/tmp/claude_envlist2.txt:115` |
| `MAX_THINKING_TOKENS` | 固定思考预算 | `DOC/env-vars.md:477` |
| `CLAUDE_CODE_EFFORT_LEVEL` | 努力度 | `DOC/env-vars.md:279` |
| `CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1` | 去掉预发布 `anthropic-beta` 头及配对字段（`context_management`、`defer_loading` 等）；**网关 400 的首选开关** | `DOC/env-vars.md:246`；`DOC/llm-gateway-connect.md:517`；`DOC/llm-gateway-protocol.md:257-270` |
| `DISABLE_PROMPT_CACHING=1` | 关闭 `cache_control` | `DOC/env-vars.md:455` |
| `API_TIMEOUT_MS` | 单次请求等待响应头的超时（默认 600000） | `DOC/env-vars.md:180` |
| settings `modelOverrides` | 把网关别名映射为某 Anthropic 模型 ID 以获得其能力认定 | `DOC/llm-gateway-protocol.md:114` |

### 4.3 关闭非必要流量与遥测

| 变量 | 效果 | 依据 |
| --- | --- | --- |
| `CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=1` | 关闭自动更新、遥测、错误上报、`/feedback`、版本/发行说明检查、fast mode 可用性探测等所有非模型请求；**不**影响 WebFetch 的域名预检（见 4.6） | `DOC/env-vars.md:261`；`DOC/llm-gateway-connect.md:360-386` |
| `DISABLE_TELEMETRY` | 关遥测；注意**任何非空值（含 `0`/`false`）都算关闭** | `DOC/env-vars.md:460` |
| `DISABLE_ERROR_REPORTING` | 关错误上报（同上语义） | `DOC/env-vars.md:446` |
| `DISABLE_AUTOUPDATER=1` | 关后台自动更新 | `DOC/env-vars.md:441` |
| `DISABLE_BUG_COMMAND=1` | 关 `/bug` | `/tmp/claude_envlist2.txt:39` |
| `DISABLE_COST_WARNINGS=1` | 关成本告警 | `DOC/env-vars.md:444` |
| `DISABLE_GROWTHBOOK=1` | 关特性开关拉取，全部用代码默认值 | `DOC/env-vars.md:449` |
| `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` | 关自动记忆（避免把用户机器上的记忆注入服务端会话） | `DOC/env-vars.md:235` |
| `CLAUDE_CODE_ENABLE_TELEMETRY` / `OTEL_*` | 只有想要 OTel 时才设 | `DOC/env-vars.md:292` |
| `DISABLE_NON_ESSENTIAL_MODEL_CALLS` | 『本版本不存在』（二进制 0 命中、`DOC/env-vars.md` 无此行） | 见 1.8 |

遥测请求只在目标主机与凭据归属一致时才附凭据；`ANTHROPIC_BASE_URL` 指向网关时，遥测发往 Anthropic 且不带网关凭据（`DOC/llm-gateway-connect.md:364`）。服务端建议全部关闭。

### 4.4 与用户 `~/.claude` 隔离

| 手段 | 效果 | 依据 |
| --- | --- | --- |
| `env.CLAUDE_CONFIG_DIR=<独立目录>` | 设置、会话历史（`projects/`）、插件、`~/.claude.json` 全部改到该目录；**推荐做法** | `DOC/env-vars.md:428`；`DOC/agent-sdk_claude-code-features.md:91` |
| `env.HOME=<独立目录>`（`T/smoke.mjs:29` 的做法） | 等效但更粗暴，会影响 git、npm、python 等子命令的 HOME | 实测会话落在 `T/home/.claude/projects/` |
| `settingSources: []` 或 `['project']` | 不读 `~/.claude/settings.json`、`~/.claude/CLAUDE.md`、用户 skills/agents | `SDK/sdk.d.ts:2312-2314` |
| `strictMcpConfig: true` | 不读 `.mcp.json`/用户 MCP/claude.ai 连接器 | `SDK/sdk.d.ts:2359-2364` |
| `persistSession: false` | 不写会话文件（代价：不可 resume） | `SDK/sdk.d.ts:1849-1855` |
| `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` | 不加载/写入自动记忆 | `DOC/env-vars.md:235` |

### 4.5 未知模型 ID 走 `ANTHROPIC_BASE_URL` 时的真实请求形状与应对

实测（`JSONL:17`，模型 `glm-5.3`；另两份 JSONL 的 `deepseek-v4-pro`、`kimi-k3` 相同）：请求带 `betas:[claude-code-20250219, interleaved-thinking-2025-05-14, thinking-token-count-2026-05-13, context-management-2025-06-27, mid-conversation-tool-changes-2026-07-01, advanced-tool-use-2025-11-20, effort-2025-11-24, thinking-display-updates-2026-08-18]`、`max_tokens:32000`、`thinking:{type:"adaptive",display:"updates"}`、`context_management:{edits:[{type:"clear_thinking_20251015",keep:"all"}]}`、`output_config:{effort:"high"}`。这与 `DOC/llm-gateway-protocol.md:103` 一致：Anthropic Messages 格式网关会收到"当前 Claude 模型接受的全部字段"，而非 Bedrock 格式下对未知 ID 的保守请求。

| 现象 | 原因 | 处理 |
| --- | --- | --- |
| 400 提到 `context_management` / `Extra inputs are not permitted` / 未知字段 | 上游不接受 beta 字段 | `CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1`（`DOC/llm-gateway-connect.md:517`） |
| 400 提到 `thinking` / `adaptive`（`Input tag 'adaptive' found`） | 上游不支持 adaptive；CLI 会自动用 `enabled` 重试一次 | `thinking:{type:'enabled', budgetTokens}` 或 `CLAUDE_CODE_DISABLE_THINKING=1`（`DOC/llm-gateway-protocol.md:236`；`DOC/env-vars.md:272`） |
| 400 拒绝 `role:"system"` 消息项 | CLI 对未知模型会把提醒文本作为 `messages` 中的 `role:"system"` 发送；被拒会自动重试一次改放 user 消息 | 一般无需处理（`DOC/llm-gateway-protocol.md:225,237`） |
| 400 拒绝 `output_config.effort` | 自动重试去掉 effort | 无需处理（`DOC/llm-gateway-protocol.md:241`） |
| 上下文溢出被网关改写成自家错误文本 | CLI 不识别为 prompt-too-long，不触发自动压缩 | `CLAUDE_CODE_MAX_CONTEXT_TOKENS` 设为真实窗口（`DOC/llm-gateway-connect.md:519`） |
| `result.total_cost_usd` 不可信 | 未知模型无价格表 → `modelUsage[].costBasis:'unknown'`、`hasUnknownModelCost:true`（`JSONL:34`） | 用 `modelUsage` 的 token 数按网关价目自行计费 |
| 后台调用打到 `claude-haiku-*` | 辅助模型默认 Haiku 别名 | `ANTHROPIC_DEFAULT_HAIKU_MODEL=<网关模型>` |
| 工具 `defer_loading`/ToolSearch 被网关拒 | beta 工具字段 | `CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1` 或 `ENABLE_TOOL_SEARCH=false`（`DOC/env-vars.md:468`） |
| 签名思考块被拒（`bound to a different conversation`） | 上游不校验/不保留签名 | CLI 自动重试去签名（`DOC/llm-gateway-protocol.md:239`） |
| `Files API is unavailable on third-party providers` | 二进制内置限制 | 不在网关下依赖 Files API（`BIN/claude` 字符串） |
| `count_tokens` 404 | 可选端点 | 自动回退字符估算（`DOC/llm-gateway-protocol.md:228`） |

### 4.6 WebSearch / WebFetch 在第三方网关下

| 工具 | 实现位置 | 在第三方网关下 | 依据 |
| --- | --- | --- | --- |
| `WebSearch` | **服务端工具**：请求 `tools` 中加入 `{type:"web_search_20250305", name:"web_search", allowed_domains, blocked_domains, max_uses:8}`，由 Anthropic API 执行；不可配置替代实现 | 需要网关/上游支持该服务端工具；不支持则报 "server-side web search is not available in this environment"（二进制内置正则）。Bedrock 不提供；Claude API、AWS Claude Platform、Foundry 可用 | `BIN/claude` 字符串 `web_search_20250305`（2 处）；`DOC/tools-reference.md:670` |
| `WebFetch` | **客户端工具**：CLI 自己 HTTP GET → HTML 转 Markdown → 用**单独一次模型调用**按 prompt 提取（走网关）；拒绝 `localhost`/无点主机名；HTTP 自动升 HTTPS；15 分钟缓存 | 可用，但有**硬编码预检** `https://api.anthropic.com/api/web/domain_info?domain=...`，不走 `ANTHROPIC_BASE_URL`，且 `CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC` 不关它；内网需 settings `skipWebFetchPreflight: true`（`SDK/sdk.d.ts:8968`）；要彻底关闭用 `CLAUDE_CODE_DISABLE_WEB_FETCH=1` | `BIN/claude` 字符串 `api.anthropic.com/api/web/domain_info`；`DOC/tools-reference.md:618-650`；`DOC/llm-gateway-connect.md:384` |

两者在 0.3.310 都是**延迟工具**（`JSONL:6`）：首轮不在 `tools` 里，模型需先 `ToolSearch`。`result.usage.server_tool_use.web_search_requests` 可统计服务端搜索次数（`T/smoke_glm.log:224`）。服务端建议：`disallowedTools:['WebSearch']`，WebFetch 视网络策略决定，或改为自己的 MCP 搜索工具。

### 4.7 Windows 注意事项

| 事项 | 事实 | 依据 |
| --- | --- | --- |
| 二进制位置 | `node_modules/@anthropic-ai/claude-agent-sdk-win32-x64/claude.exe`（arm64 为 `-win32-arm64/claude.exe`）；`manifest.json` 中 `win32-x64.binary = "claude.exe"`，约 297 MB | `SDK/manifest.json`；`SDK/README.md:42`；`SDK/sdk.mjs:157`（`a=n==="win32"?".exe":""`） |
| 安装 | 必须允许可选依赖；`npm i --omit=optional` 会导致上面的 "Native CLI binary ... not found" 错误 | `SDK/sdk.mjs:163` |
| Bash 工具依赖 Git Bash | 查找顺序（二进制逻辑）：`CLAUDE_CODE_GIT_BASH_PATH`（须存在且文件名 `bash.exe`）→ `C:\Program Files\Git\bin\bash.exe` → `C:\Program Files (x86)\Git\bin\bash.exe` → `where git` 所在目录的 `..\..\bin\bash.exe`；都找不到时提示 "Or set CLAUDE_CODE_GIT_BASH_PATH to your bash.exe location." | `BIN/claude` 字符串；`DOC/env-vars.md:306`；`DOC/setup.md:66,118-130` |
| 没有 Git Bash | 自动改用 **PowerShell 工具**（工具名 `PowerShell`，在 `BUILTIN_TOOL_NAMES` 中）；`CLAUDE_CODE_USE_POWERSHELL_TOOL=0` 可禁；两者都没有 → `startup_failure_reason:'shell_tool_missing'` | `SDK/sdk.d.ts:216,6171`；`DOC/env-vars.md:417` |
| hooks/allowedTools 要同时覆盖两个 shell 工具 | matcher 写 `Bash\|PowerShell`；`allowedTools` 里 `Bash(...)` 规则不会放行 `PowerShell` | 推论自上两行 |
| 命令行长度 | Windows 进程命令行上限 32,767 字符；插件多时用 `pluginDelivery:'initialize'` | `SDK/sdk.d.ts:2075-2088` |
| 配置目录 | `~/.claude` = `%USERPROFILE%\.claude`；`CLAUDE_CONFIG_DIR` 同样生效 | `DOC/env-vars.md:428` |
| 路径编码 | 会话文件名用 cwd 路径替换分隔符生成（Linux 实测 `-tmp-sdktest-ws`），Windows 盘符与反斜杠同样会被替换；传 `cwd` 时用 `path.resolve()` 的绝对路径，避免混用 `/` 与 `\`；`additionalDirectories` 要绝对路径；给 `env.CLAUDE_CODE_GIT_BASH_PATH` 的路径在 JSON/TS 里要双反斜杠 | 实测目录名 + `DOC/setup.md:128-130` |
| 技能脚本 | SKILL.md 的 `` !`command` `` 默认 bash，可 `shell: powershell` | `DOC/skills.md:415` |
| WSL | WSL 2 内运行即 Linux 流程，不需要 Git for Windows；WSL 1 不支持沙箱 | `DOC/setup.md:113-114` |
| 托管 CLAUDE.md | `C:\Program Files\ClaudeCode\CLAUDE.md` | `DOC/memory.md:59` |

### 4.8 平台二进制包与定位代码（可直接引用）

- 8 个可选依赖包（`SDK/package.json` `optionalDependencies`）：`@anthropic-ai/claude-agent-sdk-{linux-x64, linux-arm64, linux-x64-musl, linux-arm64-musl, darwin-x64, darwin-arm64, win32-x64, win32-arm64}`，版本号与主包一致 `0.3.310`。
- 平台包内容（`BIN/package.json`）：`{"name":"@anthropic-ai/claude-agent-sdk-linux-x64","version":"0.3.310","os":["linux"],"cpu":["x64"],"libc":["glibc"],"files":["claude","README.md","LICENSE.md"]}`；文件名 `claude`，Windows 包为 `claude.exe`（`SDK/manifest.json` 各平台 `binary` 字段）。
- 定位代码（`SDK/sdk.mjs:157`，变量名为混淆后名称，`Kc="@anthropic-ai/claude-agent-sdk"`）：

```js
function N9(e,t={}){let n=t.platform??process.platform,r=t.arch??process.arch,o=t.exists??existsSync,
 s=t.preferMusl??BYe(),i="claude",a=n==="win32"?".exe":"",
 c=(n==="android"?[`${Kc}-linux-${r}-android`]
   :n==="linux"?s?[`${Kc}-linux-${r}-musl`,`${Kc}-linux-${r}`]:[`${Kc}-linux-${r}`,`${Kc}-linux-${r}-musl`]
   :[`${Kc}-${n}-${r}`]).map((d)=>`${d}/claude${a}`);
 for(let d of c)try{let p=e(d);if(o(p))return p}catch{} return null}
// BYe(): musl 判定 = process.report.getReport().header?.glibcVersionRuntime === void 0
```

```js
// SDK/sdk.mjs:163（query() 内）
Fy=p.pathToClaudeCodeExecutable;
if(!Fy){let Ut=fileURLToPath(import.meta.url),In=createRequire(Ut),Hr=N9((xa)=>In.resolve(xa));
 if(!Hr)throw Error(`Native CLI binary for ${process.platform}-${process.arch} not found. Reinstall @anthropic-ai/claude-agent-sdk without --omit=optional, or set options.pathToClaudeCodeExecutable.`);
 Fy=Hr}
```

- 单文件打包（`bun build --compile`）时 `require.resolve` 失效，需用 `@anthropic-ai/claude-agent-sdk/extract` 的 `extractFromBunfs()` 解出二进制再传 `pathToClaudeCodeExecutable`（`SDK/README.md:25-42`）。
- zod：`sdk.d.ts:13-15` 从 `zod/v4` 与 `zod/v3` 导入类型，TypeScript 项目必须安装 zod（peer `^4.0.0`；3.25+ 实际可用）；根入口 `sdk.mjs` 内联了 zod 运行时，`/core` 入口则从你的安装导入 `zod/v4`（`SDK/README.md:21`）。`tool()` 的 schema 参数是**原始 shape 对象** `{ a: z.string() }` 而非 `z.object()`（`SDK/sdk.d.ts:130,10079`）。`T/` 下实际安装的是 zod 4.6.7。

---

## 5. 最小可运行示例（Node 22 ESM TypeScript）+ Fastify/SSE 骨架

来源：示例文件 `T/example/agent-demo.ts`、`T/example/server.ts`、`T/example/tsconfig.json`（本次编写）；API 依据为第 1、2 节引用的 `SDK/sdk.d.ts` 行号；文档 `https://code.claude.com/docs/en/agent-sdk/custom-tools`、`https://code.claude.com/docs/en/agent-sdk/hooks`、`https://code.claude.com/docs/en/agent-sdk/subagents`、`https://code.claude.com/docs/en/agent-sdk/sessions`、`https://code.claude.com/docs/en/agent-sdk/streaming-output`、`https://code.claude.com/docs/en/agent-sdk/streaming-vs-single-mode`。

验证方式（未调用任何模型 API）：
- `tsc -p tsconfig.json`（TypeScript 5.9.3，`strict`，`module: NodeNext`）对两份文件类型检查通过，解析到的正是 `T/node_modules/@anthropic-ai/claude-agent-sdk/sdk.d.ts` 与 `T/node_modules/zod/v4`；反向对照：把 `appendSystemPrompt` 或 `permissionMode:'yolo'` 写进 `Options` 会被 tsc 报错（证明检查有效，也证明 `appendSystemPrompt` 不是 Options 字段）。
- 编译后 `import('./dist/agent-demo.js')`：`createSdkMcpServer` 返回 `{type:'sdk', name:'biz', instance}`，`agents`/`hooks`/`env`（含展开的 `PATH`）构造正常；用按 `sdk.d.ts` 形状手工构造的 `init / stream_event(text_delta, input_json_delta) / assistant(tool_use Agent) / 子代理 assistant+user / task_progress / result / 未知 subtype` 消息序列驱动 `handleMessage()`，输出事件与 2.6 节的识别规则一致（子代理消息 `scope:'subagent'`，session_id 在 init 与 result 返回，未知 subtype 被忽略）。
- `query()` 本身未执行（无可用密钥）。

目录结构与依赖：

```
example/
├─ package.json      {"type":"module"}            # ESM
├─ tsconfig.json     target ES2022, module NodeNext, strict
├─ agent-demo.ts     子代理 + 进程内 MCP 工具 + hooks + 流式处理 + resume + 流式输入多轮
└─ server.ts         Fastify + SSE 骨架（未运行）
npm i @anthropic-ai/claude-agent-sdk zod fastify      # zod ^4（sdk.d.ts 从 zod/v4、zod/v3 导入类型）
npm i -D typescript @types/node                        # 运行：npx tsx agent-demo.ts
```

### 5.1 `agent-demo.ts`

```ts
// agent-demo.ts —— Node 22 ESM + TypeScript 5.x，@anthropic-ai/claude-agent-sdk@0.3.310
// 运行：npx tsx agent-demo.ts   或   tsc -p . && node dist/agent-demo.js
// 依赖：npm i @anthropic-ai/claude-agent-sdk zod   （zod ^4，类型定义从 zod/v4、zod/v3 导入）
import path from 'node:path';
import {
  query,
  tool,
  createSdkMcpServer,
  type Options,
  type SDKMessage,
  type SDKUserMessage,
  type HookCallback,
} from '@anthropic-ai/claude-agent-sdk';
import { z } from 'zod';

// ---------- 1. 进程内 MCP 工具（createSdkMcpServer + tool，schema 是"原始 shape"而非 z.object） ----------
async function fakeDb(name: string, year?: number) {
  return [{ name, year: year ?? 2026, organizer: '示例主办方', level: '国家级' }];
}

const bizServer = createSdkMcpServer({
  name: 'biz', // 工具名将是 mcp__biz__lookup_competition
  version: '1.0.0',
  tools: [
    tool(
      'lookup_competition',
      '按名称查询赛事基本信息（本地数据库）',
      {
        name: z.string().describe('赛事名称关键字'),
        year: z.number().int().optional().describe('年份，缺省为当年'),
      },
      async (args) => {
        const rows = await fakeDb(args.name, args.year);
        return { content: [{ type: 'text', text: JSON.stringify(rows, null, 2) }] };
      },
    ),
  ],
});

// ---------- 2. 钩子：记录每次工具调用，并拦截危险命令 ----------
const logPre: HookCallback = async (input, toolUseID) => {
  if (input.hook_event_name === 'PreToolUse') {
    const brief = JSON.stringify(input.tool_input).slice(0, 200);
    console.error(`[hook:pre] ${input.tool_name} id=${toolUseID} agent=${input.agent_id ?? 'main'} ${brief}`);
    const cmd = String((input.tool_input as { command?: unknown })?.command ?? '');
    if ((input.tool_name === 'Bash' || input.tool_name === 'PowerShell') && /rm\s+-rf|Remove-Item.*-Recurse/.test(cmd)) {
      return {
        hookSpecificOutput: {
          hookEventName: 'PreToolUse',
          permissionDecision: 'deny',
          permissionDecisionReason: '破坏性命令已被策略拦截',
        },
      };
    }
  }
  return {}; // 必须返回对象，否则 CLI 等到 timeout
};

const logPost: HookCallback = async (input) => {
  if (input.hook_event_name === 'PostToolUse') {
    console.error(`[hook:post] ${input.tool_name} ${input.duration_ms ?? '?'}ms`);
  }
  return {};
};

// ---------- 3. 基础 Options（第三方网关 + 服务端无人值守） ----------
const MODEL = process.env.AGENT_MODEL ?? 'glm-5.3';
const WORKSPACE = path.resolve('./ws');

export const baseOptions: Options = {
  cwd: WORKSPACE,
  model: MODEL,
  systemPrompt: { type: 'preset', preset: 'claude_code', append: '回答一律使用中文。', snapshot: true },
  settingSources: ['project'], // 只读工作区 .claude/ 与 CLAUDE.md；[] 为完全隔离
  tools: ['Read', 'Write', 'Edit', 'Bash', 'Grep', 'Glob', 'Agent', 'Skill'],
  allowedTools: ['Read', 'Write', 'Edit', 'Grep', 'Glob', 'Agent', 'Bash(python3:*)', 'mcp__biz__*'],
  disallowedTools: ['WebSearch', 'AskUserQuestion'],
  permissionMode: 'dontAsk', // 未在 allowedTools 的调用直接拒绝
  permissionPrompts: 'none',
  maxTurns: 30,
  thinking: { type: 'disabled' }, // 网关不支持 thinking 时；支持则 { type: 'enabled', budgetTokens: 4096 }
  includePartialMessages: true,
  persistSession: true,
  strictMcpConfig: true,
  mcpServers: { biz: bizServer },
  agents: {
    researcher: {
      description: '只读调研：在工作区内检索文件并归纳要点，不修改任何文件。需要查资料时委派给它。',
      prompt: '你是只读调研员。用 Read/Grep/Glob 查找资料，输出要点列表并附文件路径。不要写文件。',
      tools: ['Read', 'Grep', 'Glob'],
      model: 'inherit',
      maxTurns: 10,
      omitClaudeMd: true,
    },
  },
  hooks: {
    PreToolUse: [{ matcher: 'Bash|PowerShell|Write|Edit|Agent|mcp__biz__.*', hooks: [logPre] }],
    PostToolUse: [{ hooks: [logPost] }],
  },
  env: {
    ...process.env, // 必须展开：env 会整体替换子进程环境
    ANTHROPIC_BASE_URL: process.env.GATEWAY_URL ?? 'https://gateway.example.com',
    ANTHROPIC_AUTH_TOKEN: process.env.GATEWAY_TOKEN ?? '',
    ANTHROPIC_MODEL: MODEL,
    ANTHROPIC_DEFAULT_HAIKU_MODEL: MODEL, // 后台/辅助调用也走网关可用的模型
    CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC: '1',
    DISABLE_TELEMETRY: '1',
    DISABLE_ERROR_REPORTING: '1',
    DISABLE_AUTOUPDATER: '1',
    CLAUDE_CODE_DISABLE_AUTO_MEMORY: '1',
    CLAUDE_CONFIG_DIR: path.resolve('./agent-home'), // 与用户 ~/.claude 隔离
  },
  stderr: (d) => console.error('[cli]', d.trimEnd()),
};

// ---------- 4. 流式输出处理：把 SDKMessage 归一成前端事件 ----------
export type AgentEvent = { kind: string; [k: string]: unknown };
export type Sink = (ev: AgentEvent) => void;

const partialJson = new Map<string, string>(); // `${messageUuid}:${blockIndex}` -> 累积的 partial_json

/** 返回 session_id（init / result 上携带），其余情况返回 undefined */
export function handleMessage(m: SDKMessage, emit: Sink): string | undefined {
  switch (m.type) {
    case 'system': {
      if (m.subtype === 'init') {
        emit({ kind: 'init', sessionId: m.session_id, model: m.model, tools: m.tools, agents: m.agents, skills: m.skills });
        return m.session_id;
      }
      if (m.subtype === 'compact_boundary') emit({ kind: 'compact', trigger: m.compact_metadata.trigger, preTokens: m.compact_metadata.pre_tokens });
      if (m.subtype === 'task_progress') emit({ kind: 'subagent_progress', taskId: m.task_id, desc: m.description, tokens: m.usage.total_tokens, summary: m.summary });
      if (m.subtype === 'permission_denied') emit({ kind: 'denied', tool: m.tool_name, reason: m.decision_reason });
      return undefined;
    }
    case 'stream_event': {
      const ev = m.event;
      const scope = m.parent_tool_use_id ? 'subagent' : 'main';
      if (ev.type === 'content_block_delta') {
        if (ev.delta.type === 'text_delta') emit({ kind: 'text', scope, text: ev.delta.text });
        else if (ev.delta.type === 'thinking_delta') emit({ kind: 'thinking', scope, text: ev.delta.thinking });
        else if (ev.delta.type === 'input_json_delta') {
          const key = `${m.uuid}:${ev.index}`;
          partialJson.set(key, (partialJson.get(key) ?? '') + ev.delta.partial_json);
          emit({ kind: 'tool_input_partial', scope, index: ev.index, json: partialJson.get(key) });
        }
      } else if (ev.type === 'content_block_stop') {
        partialJson.delete(`${m.uuid}:${ev.index}`);
      }
      return undefined;
    }
    case 'assistant': {
      const scope = m.parent_tool_use_id ? 'subagent' : 'main';
      for (const block of m.message.content) {
        if (block.type === 'tool_use') {
          // name 为 'Agent'（旧版/init 列表为 'Task'）即子代理调用；记下 block.id 以匹配后续 parent_tool_use_id
          emit({ kind: 'tool_use', scope, id: block.id, name: block.name, input: block.input, agentId: m.agent_id, isSubagentCall: block.name === 'Agent' || block.name === 'Task' });
        } else if (block.type === 'text' && !m.parent_tool_use_id) {
          emit({ kind: 'text_final', text: block.text }); // 与 text_delta 二选一渲染，避免重复
        }
      }
      return undefined;
    }
    case 'user': {
      const content = m.message.content;
      if (Array.isArray(content)) {
        for (const block of content) {
          if (block.type === 'tool_result') emit({ kind: 'tool_result', id: block.tool_use_id, isError: block.is_error ?? false, scope: m.parent_tool_use_id ? 'subagent' : 'main' });
        }
      }
      return undefined;
    }
    case 'tool_progress':
      emit({ kind: 'tool_progress', name: m.tool_name, seconds: m.elapsed_time_seconds, scope: m.parent_tool_use_id ? 'subagent' : 'main' });
      return undefined;
    case 'result': {
      emit({
        kind: 'result',
        subtype: m.subtype,
        sessionId: m.session_id,
        turns: m.num_turns,
        costUsd: m.total_cost_usd, // 累计值；第三方模型为 unknown 计价，按 modelUsage 自行计费
        modelUsage: m.modelUsage,
        denials: m.permission_denials,
        terminalReason: m.terminal_reason,
        text: m.subtype === 'success' ? m.result : m.errors.join('\n'),
      });
      return m.session_id;
    }
    default:
      return undefined; // 其余 system 子类型 / 新增类型一律忽略
  }
}

// ---------- 5. 单次任务（字符串 prompt），返回 session_id ----------
export async function runOnce(prompt: string, emit: Sink, extra: Partial<Options> = {}): Promise<string | undefined> {
  const abort = new AbortController();
  let sessionId: string | undefined;
  const q = query({ prompt, options: { ...baseOptions, ...extra, abortController: abort } });
  try {
    for await (const m of q) {
      const sid = handleMessage(m, emit);
      if (sid) sessionId = sid;
    }
  } catch (e) {
    // 单消息模式下，error_* 结果在 yield 之后还会再抛一次异常（文档 streaming-vs-single-mode）；result 已处理，这里只记录
    if (!sessionId) throw e;
    console.error('[query] ended with error:', e instanceof Error ? e.message : String(e));
  }
  return sessionId;
}

// ---------- 6. 会话恢复：第二次 query 用 resume ----------
export async function demoResume(): Promise<void> {
  const emit: Sink = (ev) => console.log(JSON.stringify(ev));
  const sid = await runOnce(
    '先用 researcher 子代理调研工作区里的赛事资料并汇总要点，再调用 lookup_competition 工具核对主办方，最后把结果写入 summary.md。',
    emit,
  );
  if (!sid) throw new Error('未拿到 session_id');
  await runOnce('把 summary.md 精简到 200 字以内。', emit, { resume: sid });
  // 需要保留原会话时：{ resume: sid, forkSession: true }
}

// ---------- 7. 多轮（流式输入）：同一子进程内连续对话，可用 q.interrupt()/setModel() ----------
export async function demoStreamingInput(): Promise<void> {
  const inbox: SDKUserMessage[] = [];
  let wake: (() => void) | undefined;
  let closed = false;
  const push = (text: string) => {
    inbox.push({ type: 'user', message: { role: 'user', content: text }, parent_tool_use_id: null });
    wake?.();
  };
  async function* feed(): AsyncGenerator<SDKUserMessage> {
    while (!closed) {
      const next = inbox.shift();
      if (next) { yield next; continue; }
      await new Promise<void>((resolve) => { wake = resolve; });
    }
  }
  const q = query({ prompt: feed(), options: baseOptions });
  push('列出工作区里的文件。');
  let results = 0;
  for await (const m of q) {
    handleMessage(m, (ev) => console.log(JSON.stringify(ev)));
    if (m.type === 'result') {
      results += 1;
      if (results === 1) push('把这些文件按类型统计数量。');
      else { closed = true; wake?.(); q.close(); }
    }
  }
}

const isMain = process.argv[1] && path.resolve(process.argv[1]) === path.resolve(new URL(import.meta.url).pathname);
if (isMain) {
  demoResume().catch((e) => { console.error(e); process.exit(1); });
}
```

要点对照：
- 子代理：`agents.researcher`（`SDK/sdk.d.ts:1613`、`38-108`）；`tools` 列表里含 `Agent`，`allowedTools` 里也放 `Agent` 以免 `dontAsk` 拒绝委派；子代理消息靠 `parent_tool_use_id` / `agent_id` 区分（2.6 节）。
- 进程内工具：`createSdkMcpServer({name:'biz', tools:[tool(...)]})`（`SDK/sdk.d.ts:633-663`、`10079-10083`），工具名 `mcp__biz__lookup_competition`，放行用 `allowedTools:['mcp__biz__*']`。
- 钩子：`PreToolUse` 的 `matcher` 用正则同时覆盖 `Bash|PowerShell`（Windows 无 Git Bash 时工具名为 `PowerShell`），返回 `hookSpecificOutput.permissionDecision:'deny'` 可拦截（`SDK/sdk.d.ts:2812-2818`）。
- 流式：`includePartialMessages:true` 后处理 `stream_event.content_block_delta` 的 `text_delta / thinking_delta / input_json_delta`（`SDK/sdk.d.ts:5739-5765`）；最终文本与增量二选一渲染。
- 会话恢复：从 `init`/`result` 拿 `session_id`，第二次 `query` 传 `resume`（`SDK/sdk.d.ts:2125`）；要分叉加 `forkSession:true`。
- 多轮：`prompt` 传 `AsyncGenerator<SDKUserMessage>`，每收到一条 `result` 再推下一条消息；结束时 `q.close()`（`SDK/sdk.d.ts:3319`）。
- 错误：字符串 prompt 模式下 `error_*` 结果在 yield 后还会抛异常，`runOnce` 用 try/catch 包住（`DOC/agent-sdk_streaming-vs-single-mode.md:189`）。

### 5.2 `server.ts`（Fastify + SSE 骨架，未运行）

```ts
// server.ts —— Fastify + SSE 骨架：把 Agent SDK 消息流推给前端（骨架，未运行）
// 依赖：npm i fastify @anthropic-ai/claude-agent-sdk zod
import Fastify from 'fastify';
import { query, type Options } from '@anthropic-ai/claude-agent-sdk';
import { baseOptions, handleMessage } from './agent-demo.js';

const app = Fastify({ logger: true });

/** 正在运行的任务，便于 /api/agent/abort 取消 */
const running = new Map<string, AbortController>();

interface RunBody {
  prompt: string;
  sessionId?: string; // 传入则 resume 该会话（多轮）
  fork?: boolean;     // 与 sessionId 同用：分叉而不是续写
}

app.post<{ Body: RunBody }>('/api/agent/run', async (req, reply) => {
  const { prompt, sessionId, fork } = req.body;
  const runId = crypto.randomUUID();
  const abort = new AbortController();
  running.set(runId, abort);

  // SSE 头；X-Accel-Buffering 防止 nginx 缓冲
  reply.raw.writeHead(200, {
    'Content-Type': 'text/event-stream; charset=utf-8',
    'Cache-Control': 'no-cache, no-transform',
    Connection: 'keep-alive',
    'X-Accel-Buffering': 'no',
  });
  const send = (event: string, data: unknown) => {
    reply.raw.write(`event: ${event}\ndata: ${JSON.stringify(data)}\n\n`);
  };
  const heartbeat = setInterval(() => reply.raw.write(': ping\n\n'), 15_000);
  req.raw.on('close', () => abort.abort()); // 前端断开即终止子进程

  send('run', { runId });

  const options: Options = {
    ...baseOptions,
    abortController: abort,
    ...(sessionId ? { resume: sessionId, forkSession: !!fork } : {}),
  };

  try {
    for await (const m of query({ prompt, options })) {
      const sid = handleMessage(m, (ev) => send(ev.kind, ev));
      if (sid && m.type === 'system') send('session', { sessionId: sid }); // init 时尽早把 session_id 给前端
    }
  } catch (e) {
    send('error', { message: e instanceof Error ? e.message : String(e) });
  } finally {
    clearInterval(heartbeat);
    running.delete(runId);
    send('done', { runId });
    reply.raw.end();
  }
  return reply; // 已手动接管响应
});

app.post<{ Body: { runId: string } }>('/api/agent/abort', async (req) => {
  running.get(req.body.runId)?.abort();
  return { ok: true };
});

app.get('/healthz', async () => ({ ok: true }));

app.listen({ port: Number(process.env.PORT ?? 3000), host: '0.0.0.0' }).catch((e) => {
  app.log.error(e);
  process.exit(1);
});

/* 前端（浏览器）示意：
const es = new EventSource('/api/agent/run') // 实际需 POST，可用 fetch + ReadableStream 解析 SSE，或改成 GET+query 参数
es.addEventListener('text', (e) => appendDelta(JSON.parse(e.data).text));
es.addEventListener('tool_use', (e) => showToolCard(JSON.parse(e.data)));
es.addEventListener('subagent_progress', (e) => updateTaskRow(JSON.parse(e.data)));
es.addEventListener('result', (e) => { const r = JSON.parse(e.data); saveSessionId(r.sessionId); showCost(r.modelUsage); });
*/
```

SSE 事件约定（由 `handleMessage` 的 `kind` 直接作为 `event:` 名）：`run`、`session`、`init`、`text`、`thinking`、`tool_input_partial`、`tool_use`、`tool_result`、`tool_progress`、`subagent_progress`、`compact`、`denied`、`result`、`error`、`done`。前端以 `result.sessionId` 作为下一次请求的 `sessionId` 实现多轮；一条 SSE 连接对应一次 `query()`；断连触发 `AbortController.abort()` 终止子进程。生产环境还需：按租户隔离 `cwd` 与 `CLAUDE_CONFIG_DIR`、限制并发（每个 `query()` 都是一个约 300 MB 的子进程）、把 `modelUsage` 落库计费。

---

## 6. 版本说明：0.3.x 相对早期版本的关键变化

来源：`https://code.claude.com/docs/en/agent-sdk/migration-guide`（`DOC/agent-sdk_migration-guide.md:11-20, 93-190`）、`https://code.claude.com/docs/en/agent-sdk/typescript`（`DOC/agent-sdk_typescript.md:17, 27, 53, 246, 525, 533, 549-554, 561, 654, 875, 1514, 3093, 3273`）、`https://code.claude.com/docs/en/agent-sdk/sessions`（`DOC/agent-sdk_sessions.md:144`）、`https://code.claude.com/docs/en/agent-sdk/subagents`（`DOC/agent-sdk_subagents.md:155, 157, 326-330, 632, 723`）、`https://code.claude.com/docs/en/agent-sdk/hooks`（`DOC/agent-sdk_hooks.md:250`）、`https://code.claude.com/docs/en/agent-sdk/permissions`（`DOC/agent-sdk_permissions.md:120`）、`https://code.claude.com/docs/en/agent-sdk/modifying-system-prompts`（`DOC/agent-sdk_modifying-system-prompts.md:377-383`）、`https://code.claude.com/docs/en/agent-sdk/streaming-output`（`DOC/agent-sdk_streaming-output.md:123`）、`SDK/sdk.d.ts`、`SDK/sdk.mjs`、`SDK/manifest.json`、`SDK/package.json`。

### 6.1 包与架构

| 变化 | 说明 | 依据 |
| --- | --- | --- |
| 改名 | "Claude Code SDK" → "Claude Agent SDK"：npm `@anthropic-ai/claude-agent-sdk`，Python `claude-code-sdk` → `claude-agent-sdk`；v0.1.0 为破坏性版本 | `DOC/agent-sdk_migration-guide.md:11-20, 97` |
| 原生平台二进制 | 0.3.x 通过 8 个 `optionalDependencies` 平台包携带 bun 编译的 `claude`/`claude.exe`，**不再需要单独安装 Claude Code**；SDK 版本号与内置 CLI 版本号同步（0.3.191 ↔ 2.1.191，本机 0.3.310 ↔ 2.1.310）。`executable`/`executableArgs` 只在 `pathToClaudeCodeExecutable` 指向 `.js/.mjs/.ts` 文件时生效（沿用自早期以 JS `cli.js` 运行 CLI 的方式，推断自 `FYe` 判断与 `executable` 注释） | `SDK/package.json`；`DOC/agent-sdk_typescript.md:17`；`SDK/sdk.mjs:157,163`；`SDK/sdk.d.ts:1702-1709` |
| `manifest.json` | 新增：记录 CLI 版本、各平台二进制校验和、`sdkCompat.testedWrapperVersions`（0.3.277–0.3.306） | `SDK/manifest.json` |
| `/core` 入口 | 0.3.282 起提供精简入口（不内联 zod / MCP SDK，供自带打包器的应用使用）；需 TS 5.0+ | `DOC/agent-sdk_typescript.md:53`；`SDK/README.md:19-21` |
| `extractFromBunfs()` | 0.3.144 起，用于 `bun build --compile` 场景 | `DOC/agent-sdk_typescript.md:27` |
| V2 会话 API 移除 | 实验性 `createSession()`（`send/stream` 模式）在 0.3.142 删除，统一用 `query()` + `resume/continue/forkSession` | `DOC/agent-sdk_sessions.md:144` |

### 6.2 默认值变化（最容易踩坑）

| 项 | 早期 | 现在（0.3.310） | 依据 |
| --- | --- | --- | --- |
| `systemPrompt` 省略 | v0.0.x 默认使用 Claude Code 完整系统提示词 | v0.1.0 起为**极简默认提示词**（实测只有一句 "You are a Claude agent, built on Anthropic's Claude Agent SDK."）；要 Claude Code 行为需 `{type:'preset', preset:'claude_code'}` | `DOC/agent-sdk_migration-guide.md:117-142`；`JSONL:14`；`DOC/agent-sdk_typescript.md:554` |
| `settingSources` 省略 | v0.1.0 曾短暂改为"不加载任何文件系统设置"，随后回退 | 省略 = `['user','project','local']`（与 CLI 一致）；`[]` 才是隔离 | `DOC/agent-sdk_migration-guide.md:178-190`；`SDK/sdk.d.ts:2312`；`SDK/sdk.mjs:157` |
| `permissionMode` 省略 | 早期固定 `default` | d.ts：settings `permissions.defaultMode`，否则 `auto`（不可用则 `default`）；文档：CLI ≥ 2.1.302 无其他设置时为 `default`。两处表述不一致 → **一律显式设置** | `SDK/sdk.d.ts:2019-2024`；`DOC/agent-sdk_permissions.md:120`；`DOC/agent-sdk_typescript.md:533` |
| 内置工具集合 | 含 Grep/Glob | 原生版可能不带 Grep/Glob（改用 Bash find/grep），需在 `tools`/`allowedTools` 显式列出；WebSearch/WebFetch/NotebookEdit 等变为 ToolSearch 延迟加载 | `SDK/sdk.d.ts:1675-1676`；`JSONL:6,14` |
| 子代理默认 | 前台阻塞 | 默认后台化（`async_launched`），结果经 `task_notification` 送达；并发 20 / 深度 3 上限（0.3.219 起） | `DOC/agent-sdk_subagents.md:632` 及同页 |
| 系统提示词快照 | 每次请求重新渲染 | CLI ≥ 2.1.265 默认记录（`snapshot: true` 语义），resume 时沿用旧 append | `SDK/sdk.d.ts:2423-2450`；`DOC/agent-sdk_modifying-system-prompts.md:377-383` |

### 6.3 弃用与改名

| 旧 | 新 | 依据 |
| --- | --- | --- |
| `maxThinkingTokens` / `Query.setMaxThinkingTokens()` | `thinking: {type:'adaptive'\|'enabled'\|'disabled'}` | `SDK/sdk.d.ts:1955-1960, 2992`；`DOC/agent-sdk_typescript.md:525, 654` |
| `allowedTools: ['Skill']`、`AgentDefinition.tools: ['Skill']` | `skills` 选项 / 字段 | `SDK/sdk.d.ts:1619, 45` |
| 工具名 `Task` | `Agent`（CLI ≥ 2.1.63 的 `tool_use` 用 `Agent`，init 列表仍显示 `Task`）；其他别名 `KillShell`/`KillBash`→`TaskStop`、`ListPeers`→`ListAgents`、`Brief`→`SendUserMessage`、`ListMcpResources`→`ListMcpResourcesTool`、`ReadMcpResource`→`ReadMcpResourceTool`、`ReadMcpResourceDir`→`ReadMcpResourceDirTool`（`LEGACY_TOOL_NAME_ALIASES`） | `SDK/sdk.d.ts:1103`；`SDK/sdk.mjs:126`；`DOC/agent-sdk_subagents.md:330` |
| Agent 工具输入的 `mode`、`team_name` | 忽略（CLI ≥ 2.1.212），子代理权限模式由继承规则/`permissionMode` 决定 | `DOC/agent-sdk_typescript.md:3093-3105` |
| `shell_id` | `task_id` | `DOC/agent-sdk_typescript.md:3273` |
| settings `includeCoAuthoredBy` | `attribution` | `SDK/sdk.d.ts:7114-7117` |
| 环境变量 `ANTHROPIC_SMALL_FAST_MODEL` | `ANTHROPIC_DEFAULT_HAIKU_MODEL` | `DOC/env-vars.md:153,174` |
| PostToolUse 的 `updatedMCPToolOutput` | `updatedToolOutput`（对所有工具有效） | `DOC/agent-sdk_typescript.md:2964` |

### 6.4 功能的版本门槛（升级/兼容时查）

| 功能 | 需要 | 依据 |
| --- | --- | --- |
| `SDKAssistantMessage.agent_id` | SDK ≥ 0.3.292 | `DOC/agent-sdk_typescript.md:1514` |
| `AgentDefinition.autoCompactWindow` | SDK ≥ 0.3.296（Python 无） | `DOC/agent-sdk_subagents.md:155` |
| `AgentDefinition.omitClaudeMd` | SDK ≥ 0.3.271（Python 无） | `DOC/agent-sdk_subagents.md:157` |
| `skills` 名称校验（拒绝通配/畸形名） | SDK ≥ 0.3.221 | `DOC/agent-sdk_typescript.md:550` |
| PostToolUse `classifierContext` | SDK ≥ 0.3.236 | `DOC/agent-sdk_hooks.md:250` |
| 子代理并发/深度/预算限制 | SDK ≥ 0.3.219 | `DOC/agent-sdk_subagents.md:632` |
| `Workflow` 工具 | SDK ≥ 0.3.149 | `DOC/agent-sdk_subagents.md:723` |
| `createSdkMcpServer({timeout})` | SDK ≥ 0.3.248 | `DOC/agent-sdk_typescript.md:246` |
| `verbatimPrompts` / `client_composed` | CLI ≥ 2.1.248 | `SDK/sdk.d.ts:1927` |
| `pluginDelivery: 'initialize'` | CLI ≥ 2.1.261 | `SDK/sdk.d.ts:2083` |
| 系统提示词记录（snapshot） | CLI ≥ 2.1.265 | `DOC/agent-sdk_modifying-system-prompts.md:383` |
| 流被切断仍保证 `message_stop` | CLI ≥ 2.1.290（SDK ≥ 0.3.290） | `DOC/agent-sdk_streaming-output.md:123` |
| 默认权限模式 `default` | CLI ≥ 2.1.302 | `DOC/agent-sdk_permissions.md:120` |
| `getContextUsage()` 新字段 | SDK ≥ 0.3.304 | `DOC/agent-sdk_typescript.md:875` |
| 数据驻留计价（`inference_geo:"us"` ×1.1） | SDK ≥ 0.3.239 | `DOC/agent-sdk_cost-tracking.md:21` |

### 6.5 对实现者的兼容策略

1. 锁定 SDK 精确版本（`"@anthropic-ai/claude-agent-sdk": "0.3.310"`），升级时重跑第 5 节的 `tsc` 类型检查——`Options`/`SDKMessage` 每个小版本都在增删字段。
2. 消息消费端对未知 `type`/`subtype` 静默忽略；工具名判断同时接受 `Agent` 与 `Task`。
3. 不依赖任何默认值：显式写 `systemPrompt`、`settingSources`、`permissionMode`、`tools`、`model`、`thinking`。
4. 成本以 `modelUsage` 的 token 数自行计价，不信任第三方模型的 `total_cost_usd`。
5. 第三方网关先用 `CLAUDE_CODE_DISABLE_EXPERIMENTAL_BETAS=1` + `thinking:{type:'disabled'}` 跑通，再逐项放开。

---

## 附录 A：引用索引

| 缩写 | 实际路径 / URL |
| --- | --- |
| `SDK/sdk.d.ts` | `/tmp/sdktest/node_modules/@anthropic-ai/claude-agent-sdk/sdk.d.ts`（10302 行） |
| `SDK/sdk.mjs` | 同目录 `sdk.mjs`（压缩代码，169 行；第 126 行别名表，157 行 ProcessTransport/flag 映射/二进制定位，163 行 `query()`） |
| `SDK/README.md`、`SDK/package.json`、`SDK/manifest.json` | 同目录 |
| `BIN/` | `/tmp/sdktest/node_modules/@anthropic-ai/claude-agent-sdk-linux-x64/`（`claude` 二进制 291 MB，`package.json`） |
| `T/smoke.mjs`、`T/smoke_{glm,deepseek,kimi}.log` | `/tmp/sdktest/` |
| `T/example/*` | `/tmp/sdktest/example/`（本手册第 5 节示例） |
| `JSONL` | `/tmp/sdktest/home/.claude/projects/-tmp-sdktest-ws/42e923ab-a86f-4c9e-9491-65c4b8a6fa63.jsonl`（另两份 `9cb53a77-*.jsonl`、`c156f34d-*.jsonl` 对应 deepseek/kimi） |
| `/tmp/claude_help.txt` | `claude --help` 输出（flag 清单） |
| `/tmp/claude_envlist.txt`、`/tmp/claude_envlist2.txt` | 从二进制提取的环境变量名清单（95 / 192 个） |
| `DOC/<name>.md` | `https://code.claude.com/docs/en/agent-sdk/<name去掉 agent-sdk_ 前缀>` 或 `https://code.claude.com/docs/en/<name>`，即用户给出的 `docs.claude.com/en/api/agent-sdk/...`、`docs.claude.com/en/docs/claude-code/...` 的重定向目标；行号为 Markdown 快照行号 |
