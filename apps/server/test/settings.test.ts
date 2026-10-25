import path from 'node:path';
import os from 'node:os';
import fsp from 'node:fs/promises';
import { afterAll, beforeAll, describe, expect, it } from 'vitest';
import { DEFAULT_SETTINGS, SettingsStore, applyEnvOverrides, deepMerge, isMaskedKey, maskApiKey } from '../src/services/settings.js';
import { redact, registerSecret } from '../src/util/log.js';

const FAKE_KEY = 'sk-testkey-0123456789abcdefXYZ';

describe('settings helpers', () => {
  it('maskApiKey keeps only the prefix and the last 4 chars', () => {
    expect(maskApiKey(FAKE_KEY)).toBe('sk-****fXYZ');
    expect(maskApiKey('plainsecret')).toBe('****cret');
    expect(maskApiKey('')).toBe('');
    expect(isMaskedKey('sk-****fXYZ')).toBe(true);
    expect(isMaskedKey('****cret')).toBe(true);
    expect(isMaskedKey(FAKE_KEY)).toBe(false);
  });

  it('deepMerge merges objects recursively and replaces arrays/scalars', () => {
    const out = deepMerge({ a: { b: 1, c: [1, 2] }, d: 'x' }, { a: { c: [3] }, d: 'y', e: 1 });
    expect(out).toEqual({ a: { b: 1, c: [3] }, d: 'y', e: 1 });
  });

  it('env vars override the file (contract §9) and report apiKeySource', () => {
    const base = deepMerge(DEFAULT_SETTINGS, { api: { apiKey: 'from-file-key-1234' } });
    const none = applyEnvOverrides(DEFAULT_SETTINGS, {});
    expect(none.apiKeySource).toBe('none');
    const file = applyEnvOverrides(base, {});
    expect(file.apiKeySource).toBe('settings');
    expect(file.settings.api.apiKey).toBe('from-file-key-1234');
    const env = applyEnvOverrides(base, { BIZGOLD_API_KEY: FAKE_KEY, BIZGOLD_API_BASE_URL: 'https://gw.example.com', BIZGOLD_MODEL: 'kimi-k3', BIZGOLD_PYTHON: '/opt/py/bin/python' });
    expect(env.apiKeySource).toBe('env');
    expect(env.settings.api.apiKey).toBe(FAKE_KEY);
    expect(env.settings.api.baseUrl).toBe('https://gw.example.com');
    expect(env.settings.api.models.main).toBe('kimi-k3');
    expect(env.settings.api.models.fast).toBe(DEFAULT_SETTINGS.api.models.fast);
    expect(env.settings.python.path).toBe('/opt/py/bin/python');
    // 原对象不被修改
    expect(base.api.apiKey).toBe('from-file-key-1234');
  });

  it('redact scrubs registered secrets and sk- tokens from logs/events', () => {
    registerSecret(FAKE_KEY);
    expect(redact(`Authorization: Bearer ${FAKE_KEY}`)).not.toContain(FAKE_KEY);
    expect(redact('token sk-abcdefghijklmnop end')).toBe('token sk-**** end');
    expect(redact('ANTHROPIC_AUTH_TOKEN=secret-value-123')).toBe('ANTHROPIC_AUTH_TOKEN=****');
  });
});

describe('SettingsStore', () => {
  const dir = path.join(os.tmpdir(), `bizgold-settings-${process.pid}`);
  beforeAll(async () => {
    await fsp.rm(dir, { recursive: true, force: true });
  });
  afterAll(async () => {
    await fsp.rm(dir, { recursive: true, force: true });
  });

  it('returns defaults when no file exists, masks on read, merges on update, never persists masked values', async () => {
    const store = new SettingsStore(dir);
    const pub = await store.getPublic();
    expect(pub.api.apiKey).toBe('');
    expect(pub.api.models.main).toBe(DEFAULT_SETTINGS.api.models.main);
    expect(pub.run.maxTurns).toBe(1500);

    const updated = await store.update({ api: { apiKey: FAKE_KEY, models: { main: 'glm-5.3' } }, run: { reviewRounds: 2 }, search: { bochaKey: 'bocha-key-abcdef' } });
    expect(updated.api.apiKey).toBe('sk-****fXYZ');
    expect(updated.api.apiKeySource).toBe('settings');
    expect(updated.api.models.main).toBe('glm-5.3');
    expect(updated.api.models.fast).toBe(DEFAULT_SETTINGS.api.models.fast);
    expect(updated.run.reviewRounds).toBe(2);
    expect(updated.run.maxTurns).toBe(1500);
    expect(updated.search.bochaKey).toBe('****cdef');

    // 把脱敏串原样 PUT 回来不能覆盖真实密钥
    const again = await store.update({ api: { apiKey: 'sk-****fXYZ', apiKeySource: 'env', baseUrl: 'https://b.example.com' }, search: { bochaKey: '****cdef' } });
    expect(again.api.baseUrl).toBe('https://b.example.com');
    expect((await store.apiKey()).key).toBe(FAKE_KEY);
    expect((await store.get()).search.bochaKey).toBe('bocha-key-abcdef');

    const onDisk = JSON.parse(await fsp.readFile(path.join(dir, 'settings.json'), 'utf8')) as { api: Record<string, unknown> };
    expect(onDisk.api.apiKey).toBe(FAKE_KEY);
    expect('apiKeySource' in onDisk.api).toBe(false);
    if (process.platform !== 'win32') {
      const st = await fsp.stat(path.join(dir, 'settings.json'));
      expect(st.mode & 0o777).toBe(0o600);
    }

    // 新实例从磁盘读取，结果一致
    const fresh = new SettingsStore(dir);
    expect((await fresh.getPublic()).api.models.main).toBe('glm-5.3');
  });
});
