import { AgentRunner } from './agent/runner.js';
import { HealthService } from './services/health.js';
import { KnowledgeService } from './services/knowledge.js';
import { PythonRunner } from './services/python.js';
import { SettingsStore } from './services/settings.js';
import { ProjectStore } from './services/store.js';
import { TrackService } from './services/tracks.js';

export interface AppContext {
  settings: SettingsStore;
  store: ProjectStore;
  python: PythonRunner;
  tracks: TrackService;
  knowledge: KnowledgeService;
  health: HealthService;
  runner: AgentRunner;
  workspaceDir: string;
}

export async function createContext(opts: { settingsDir?: string; workspaceDir?: string } = {}): Promise<AppContext> {
  const settings = new SettingsStore(opts.settingsDir);
  const workspaceDir = opts.workspaceDir ?? (await settings.workspaceDir());
  const store = new ProjectStore(workspaceDir);
  const python = new PythonRunner(async () => (await settings.get()).python);
  const tracks = new TrackService();
  const knowledge = new KnowledgeService();
  const runner = new AgentRunner({ store, settings, python, tracks });
  const health = new HealthService(python, async () => ({
    workspaceDir,
    apiKeySource: (await settings.apiKey()).source,
    activeRuns: runner.activeCount,
  }));
  return { settings, store, python, tracks, knowledge, health, runner, workspaceDir };
}
