import { buildApp } from './app.js';
import { createContext } from './context.js';
import { HOST, PORT, VERSION } from './config.js';
import { createLogger, errorMessage } from './util/log.js';

const log = createLogger('server');

async function main(): Promise<void> {
  const ctx = await createContext();
  const app = await buildApp(ctx);

  // 首次启动检查工具包 Python 依赖，缺失则后台安装（不阻塞监听）
  ctx.python.ensureDeps();
  // 上次进程异常退出时留下的 running 项目 → paused
  await ctx.runner.recoverStale().catch((e) => log.warn('recoverStale failed', errorMessage(e)));

  let shuttingDown = false;
  const shutdown = async (signal: string) => {
    if (shuttingDown) return;
    shuttingDown = true;
    log.info(`${signal} received, shutting down`);
    // 先停 Agent 子进程（runner 自己也挂了信号钩子，这里保证 HTTP 层也正常收尾）
    const timer = setTimeout(() => process.exit(1), 10_000);
    try {
      await app.close();
    } catch (e) {
      log.warn('close error', errorMessage(e));
    } finally {
      clearTimeout(timer);
      process.exit(0);
    }
  };
  process.on('SIGINT', () => void shutdown('SIGINT'));
  process.on('SIGTERM', () => void shutdown('SIGTERM'));

  try {
    await app.listen({ port: PORT, host: HOST });
  } catch (e) {
    log.error(`listen failed on ${HOST}:${PORT}`, errorMessage(e));
    process.exit(1);
  }
  const apiKey = await ctx.settings.apiKey();
  log.info(`BizGold server v${VERSION} listening on http://${HOST}:${PORT}  workspace=${ctx.workspaceDir}  apiKey=${apiKey.source}`);
}

main().catch((e) => {
  log.error('fatal', errorMessage(e));
  process.exit(1);
});
