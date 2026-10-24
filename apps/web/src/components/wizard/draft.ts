import type { ModelMap, ProjectCreateInput, ProjectSubject, ResearchMode, TeamMember, TrackId } from '../../types';
import { loadJson, removeKey, saveJson } from '../../lib/storage';

export interface Draft {
  track: TrackId | null;
  group: string;
  name: string;
  teamName: string;
  subject: ProjectSubject;
  team: TeamMember[];
  params: {
    researchMode: ResearchMode;
    pageTarget: number;
    pitchMinutes: 8 | 10 | 15;
    reviewRounds: number;
    targetScore: number;
    models: Partial<ModelMap>;
  };
}

export const DRAFT_KEY = 'bizgold.newProject.draft';
export const DRAFT_STEP_KEY = 'bizgold.newProject.step';

export function emptyDraft(): Draft {
  return {
    track: null,
    group: '',
    name: '',
    teamName: '',
    subject: { type: 'brand', name: '', company: '', isAssignedCase: false, caseTitle: '', brief: '', region: '', industry: '' },
    team: [{ role: '队长 / 总策划', major: '', strength: '' }],
    params: { researchMode: 'synthetic', pageTarget: 60, pitchMinutes: 10, reviewRounds: 3, targetScore: 85, models: {} },
  };
}

export function loadDraft(): Draft {
  const d = loadJson<Partial<Draft> | null>(DRAFT_KEY, null);
  if (!d) return emptyDraft();
  const base = emptyDraft();
  return {
    ...base,
    ...d,
    subject: { ...base.subject, ...(d.subject ?? {}) },
    params: { ...base.params, ...(d.params ?? {}), models: { ...(d.params?.models ?? {}) } },
    team: d.team?.length ? d.team : base.team,
  };
}

export function saveDraft(d: Draft): void {
  saveJson(DRAFT_KEY, d);
}

export function clearDraft(): void {
  removeKey(DRAFT_KEY);
  removeKey(DRAFT_STEP_KEY);
}

export type Errors = Partial<Record<string, string>>;

export function validateStep1(d: Draft): Errors {
  const e: Errors = {};
  if (!d.track) e.track = '请选择赛道';
  if (!d.group) e.group = '请选择组别';
  return e;
}

export function validateStep2(d: Draft): Errors {
  const e: Errors = {};
  if (d.name.trim().length < 2) e.name = '作品名至少 2 个字';
  else if (d.name.trim().length > 60) e.name = '作品名不超过 60 个字';
  if (d.teamName.trim().length < 2) e.teamName = '团队名至少 2 个字';
  if (!d.subject.name.trim()) e['subject.name'] = '请填写品牌 / 项目名';
  if (d.subject.isAssignedCase && !d.subject.caseTitle) e['subject.caseTitle'] = '请选择命题';
  if (d.subject.brief.trim().length < 20) e['subject.brief'] = '项目简介至少 20 个字，便于智能体准确读题';
  if (!d.team.length || d.team.some((m) => !m.role.trim())) e.team = '至少一名成员，且每位成员需填写角色';
  const forbidden = /(大学|学院|学校|同学|老师|教授)/;
  const hit = d.team.find((m) => forbidden.test(`${m.role}${m.major}${m.strength}`));
  if (hit) e.team = '匿名制：团队信息中不得出现院校或人名相关字样';
  if (forbidden.test(d.teamName)) e.teamName = '匿名制：团队名不得包含院校名称';
  return e;
}

export function validateStep3(d: Draft, pageLimit: number): Errors {
  const e: Errors = {};
  if (d.params.pageTarget > pageLimit) e.pageTarget = `目标页数不得超过赛道上限 ${pageLimit} 页`;
  if (!d.params.models.main) e.models = '请选择主代理模型';
  return e;
}

export function toCreateInput(d: Draft): ProjectCreateInput {
  return {
    name: d.name.trim(),
    teamName: d.teamName.trim(),
    track: d.track ?? 'brand',
    group: d.group,
    subject: { ...d.subject, name: d.subject.name.trim(), brief: d.subject.brief.trim() },
    team: d.team.map((m) => ({ role: m.role.trim(), major: m.major.trim(), strength: m.strength.trim() })),
    params: {
      researchMode: d.params.researchMode,
      pageTarget: d.params.pageTarget,
      pitchMinutes: d.params.pitchMinutes,
      reviewRounds: d.params.reviewRounds,
      targetScore: d.params.targetScore,
      language: 'zh-CN',
      models: d.params.models,
    },
  };
}
