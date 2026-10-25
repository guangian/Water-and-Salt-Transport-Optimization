"""五个赛道的内置机器检查配置（契约 §7.8 `tracks/<id>/checks.json` 的兜底默认值）。

`tracks/<track>/checks.json` 可能尚未生成，或只覆盖部分字段；`load_track_config()` 以内置默认为底、
用赛道文件（先找工作区下 `tracks/<id>/checks.json`，再找仓库 `agent/tracks/<id>/checks.json`）逐字段覆盖。

必备要素清单取自 `docs/research/01_赛事情报.md` 第二章各赛道“必备章节”（品牌 10 / 数字营销 9 / 案例 8 / 会奖 10），
创业计划赛道细则未规定章节清单，按 `02_国金级材料标准.md` §1.2 的 13 章 BP 结构；别名用于标题层级匹配（含常见同义写法）。
禁语表与匿名正则取自 `02_国金级材料标准.md` §6.7。

两套字段名：内置默认用扁平字段（`sourceOfficialRatio` / `pitchMinutes` / `deckPages` / `qaMin` / `*Required` 开关），
真实赛道文件用 `sourceOfficialRatioMin` / `sourceFreshRatioMin` / `requiredEvidence{}` / `numberKeys[]` / `pitch{}`；
`load_track_config()` 两套都认，并在返回前由 `normalize_config()` 互相补齐（以赛道文件为准），调用方任选一套读取即可。
"""
from __future__ import annotations

import copy
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from ._util import read_json

# 禁语：纯字符串 = 命中即计；字典 {term, needSource:true} = 同句无来源短语时才命中（“据统计”“填补空白”）
# 广告法绝对化用语只收具体组合，“最 / 第一”单字不入表（“第一年”“最大”会误伤）。
FORBIDDEN_TERMS_DEFAULT: List[Union[str, Dict[str, Any]]] = [
    "颠覆", "独一无二", "蓝海", "众所周知", "前景广阔", "我们学校", "老师说",
    {"term": "据统计", "needSource": True},
    {"term": "填补空白", "needSource": True},
    {"term": "填补国内空白", "needSource": True},
    "国家级", "顶级", "唯一", "最佳", "最高级", "全国第一", "行业第一", "全球第一", "第一品牌", "最先进", "史无前例",
]

# 匿名正则（§6.7）：院校 / 校徽 / 导师 / 学号 / 手机号 / 邮箱；白名单在扫描前用占位符保护
ANONYMITY_REGEX_DEFAULT = r"(大学|学院|University|College|校徽|指导教师|指导老师|导师|学号|1[3-9]\d{9}|[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,})"
# 常见误伤：消费人群“大学生”、出版社、研究机构名
ANONYMITY_WHITELIST_DEFAULT: List[str] = [
    "大学生", "大学生群体", "在校大学生", "大学城", "出版社", "中国社会科学院", "中国科学院", "社科院",
    "国务院", "研究院", "科学院", "工程院", "医院", "法院", "设计学院派",
]

# 来源分级：官方来源的域名 / 类型关键词（统计局、政府、协会、年报招股书、国际组织）
OFFICIAL_DOMAIN_HINTS: List[str] = [
    "gov.cn", "stats.gov", "cnnic", "pbc.gov", "mofcom", "miit", "ndrc", "sse.com.cn", "szse.cn", "cninfo", "hkex",
    "sec.gov", "worldbank", "imf.org", "oecd.org", "un.org", "who.int", "wto.org", "europa.eu", "ceic", "wind",
    "cas.cn", "cass", "chinatax", "customs.gov", "association", "协会", "年报", "招股", "统计年鉴",
]
OFFICIAL_TYPES = {"gov", "annual-report", "standard", "database"}

# 广告法绝对化用语（CLAUDE.md §2.2 / §6.7）。赛道文件里的单字 “最” 与 “第一” 由 `expand_adlaw_terms()` 展开为具体组合，
# 否则 “最低点 / 最近三年 / 第一年 / 第一阶段” 全部误报，硬门槛永远过不了。
ADLAW_TERMS_DEFAULT: List[str] = [
    "最佳", "最好", "最优", "最强", "最美", "最先进", "最高级", "最低价", "最便宜", "最专业", "最权威", "最顶级", "最领先",
    "最受欢迎", "最高端", "最尖端", "最完美", "最安全", "最畅销", "第一品牌", "全网第一", "行业第一", "全国第一", "全球第一",
    "世界第一", "销量第一", "国家级", "世界级", "顶级", "唯一", "独家", "极致", "首选", "史无前例", "100%有效", "绝对",
]
# “第一” 后接这些字视为序数 / 计量而非绝对化宣传
ADLAW_ORDINAL_SUFFIX = "年|阶段|季度|季|步|批|次|期|章|节|轮|天|周|月|部分|层|象限|版|印象|发明人|时间|线|产业|方|类|个|家|款|件|位|条|行|列|页|项|篇|张|场|组|轮次|回|届"
ADLAW_WHITELIST_DEFAULT: List[str] = [
    "最低点", "最高点", "最近", "最终", "最初", "最早", "最晚", "最后", "最小", "最大", "最高值", "最低值", "最多", "最少",
    "最新", "最常", "最主要", "最核心", "最关键", "最大化", "最小化", "最优化", "最快", "最低", "最高", "最长", "最短",
    "最深", "最浅", "最重", "最轻", "最远", "最近似", "最优解", "唯一性", "唯一标识", "首选渠道", "第一发明人",
    "国家级非遗", "国家级非物质文化遗产", "国家级经济技术开发区", "国家级高新区", "国家级新区", "国家级标准", "顶级域名",
    "独家协议", "独家授权", "独家代理", "独家经销", "独家合作", "独家访谈",
]
# 需要 “同句有来源” 才放过的禁语（赛道文件若以纯字符串给出同名词，也按此语义处理）
FORBIDDEN_NEED_SOURCE: List[str] = ["据统计", "填补空白", "填补国内空白", "填补市场空白", "预计将"]

# 弱评价词（CLAUDE.md §2.2：每百字 ≤ 0.1）；`WEAK_TERM_EXCEPTIONS` 为统计口径里应当放过的专业搭配
WEAK_TERMS_DEFAULT: List[str] = ["较好", "较高", "较大", "显著", "良好", "合理", "有效", "明显", "一定程度"]
WEAK_TERM_EXCEPTIONS: List[str] = [
    "有效问卷", "有效样本", "有效率", "有效期", "有效回收", "有效性", "有效数字", "有效覆盖", "显著性", "显著水平", "合理化",
    "合理区间", "有效访谈", "有效应答",
]

# 赛道文件字段名 → 内置字段名（互为别名；赛道文件优先）
FIELD_ALIASES: Dict[str, str] = {
    "sourceOfficialRatioMin": "sourceOfficialRatio",
    "sourceFreshRatioMin": "sourceFreshRatio",
    "sourceDomainMin": "sourceDomainsMin",
    "chartsMin": "chartMin",
    "wordsMin": "wordMin",
    "pagesLimit": "pageLimit",
}
# requiredEvidence 开关 ↔ 内置 *Required 布尔字段
EVIDENCE_FLAGS: Dict[str, str] = {
    "has_tam_sam_som": "tamRequired",
    "has_unit_economics": "unitEconomicsRequired",
    "has_sensitivity": "sensitivityRequired",
    "three_statement_check": "financialModelRequired",
    "budget_recalc": "budgetRequired",
    "kpi_table": "kpiRequired",
    "has_competitor_matrix": "competitorMatrixRequired",
    "has_primary_research": "primaryResearchRequired",
}

_BASE: Dict[str, Any] = {
    "pageLimit": 80,
    "wordMin": 0,
    "chartMin": 15,
    "sourceMin": 15,
    "sourceDomainsMin": 8,
    "sourceOfficialRatio": 0.4,
    "sourceFreshRatio": 0.7,
    "sourceFreshYears": 3,
    "numericDensity": 15,
    "numericDensityChapters": {},
    "pageTargetTolerance": 0.1,
    "pitchMinutes": 10,
    "deckPages": [16, 20],
    "deckPagesAppendixMax": 4,
    "deckChartMin": 8,
    "slideWordsMax": 120,
    "scriptCharsPerMinute": 190,
    "qaMin": 40,
    "qaWarn": 20,
    "qaDimensions": ["项目与创新", "市场与客户", "商业模式与财务", "团队与运营", "风险与合规", "社会价值与成长"],
    "sampleSizeMin": 384,
    "sampleSizeWarn": 200,
    "riskMatrixMin": 6,
    "tamRequired": False,
    "unitEconomicsRequired": False,
    "sensitivityRequired": False,
    "financialModelRequired": False,
    "budgetRequired": False,
    "kpiRequired": False,
    "competitorMatrixRequired": True,
    "primaryResearchRequired": True,
    "numberKeys": [],
    "forbiddenTerms": FORBIDDEN_TERMS_DEFAULT,
    "advertisingLawTerms": ADLAW_TERMS_DEFAULT,
    "advertisingLawWhitelist": ADLAW_WHITELIST_DEFAULT,
    "weakTerms": WEAK_TERMS_DEFAULT,
    "weakTermsPer100": 0.1,
    "anonymityRegex": ANONYMITY_REGEX_DEFAULT,
    "anonymityWhitelist": ANONYMITY_WHITELIST_DEFAULT,
}


def _el(eid: str, name: str, *aliases: str) -> Dict[str, Any]:
    return {"id": eid, "name": name, "aliases": list(aliases)}


TRACK_DEFAULTS: Dict[str, Dict[str, Any]] = {
    "brand": {
        **_BASE,
        "name": "品牌策划赛道", "docName": "品牌策划书",
        "chartMin": 20, "sourceMin": 15, "numericDensity": 15,
        "budgetRequired": True, "kpiRequired": True,
        "requiredElements": [
            _el("toc", "目录"),
            _el("summary", "概要", "执行摘要", "摘要", "方案概要", "内容概要"),
            _el("intro", "引言", "前言", "绪论", "项目背景", "策划背景", "命题解读"),
            _el("status", "品牌现状分析", "现状分析", "品牌诊断", "市场与品牌现状", "环境分析", "现状诊断", "品牌现状"),
            _el("positioning", "品牌定位与核心价值", "品牌定位", "核心价值", "定位与核心价值", "STP", "品牌战略"),
            _el("elements", "品牌元素设计", "品牌元素", "元素设计", "视觉识别", "VI 设计", "VI设计", "品牌识别系统", "创意设计"),
            _el("promotion", "品牌推广与传播", "品牌推广", "传播策略", "推广与传播", "整合营销传播", "传播执行", "推广计划", "传播计划"),
            _el("protection", "品牌资产保护和风险管理", "品牌资产保护", "风险管理", "资产保护与风险", "风险与保护", "风险管理与应急", "品牌保护"),
            _el("budget", "财务预算", "预算", "费用预算", "预算与排期", "预算与效果", "投入预算"),
            _el("appendix", "附录和参考资料", "附录", "参考资料", "参考文献", "附录与参考文献", "附件"),
        ],
    },
    "startup": {
        **_BASE,
        "name": "创业计划赛道", "docName": "商业计划书",
        "wordMin": 3000, "chartMin": 15, "sourceMin": 20, "numericDensity": 20,
        "tamRequired": True, "unitEconomicsRequired": True, "sensitivityRequired": True, "financialModelRequired": True,
        "requiredElements": [
            _el("summary", "执行摘要", "摘要", "概要", "项目概述", "项目概要"),
            _el("background", "项目背景与痛点", "项目背景", "痛点", "背景与痛点", "问题与机会", "行业背景", "痛点分析"),
            _el("product", "产品与技术", "产品与服务", "产品", "技术方案", "解决方案", "产品介绍", "核心技术", "产品与创新"),
            _el("market", "市场分析", "行业分析", "市场研究", "市场规模", "目标市场", "市场机会", "市场与客户"),
            _el("competition", "竞争分析", "竞争格局", "竞品分析", "竞争对手", "竞争优势", "竞争与定位"),
            _el("model", "商业模式", "盈利模式", "商业模式画布", "收入模式", "商业模式与单位经济"),
            _el("marketing", "营销与渠道", "营销策略", "市场推广", "渠道策略", "营销计划", "获客策略", "营销与销售"),
            _el("operations", "运营与实施", "运营计划", "实施计划", "运营模式", "实施路径", "里程碑", "运营与里程碑", "执行计划"),
            _el("team", "团队", "团队介绍", "团队与组织", "核心团队", "组织架构", "团队与股权"),
            _el("finance", "财务与融资", "财务预测", "财务分析", "融资计划", "财务计划", "财务与融资计划", "财务预测与融资"),
            _el("risk", "风险", "风险分析", "风险与对策", "风险控制", "风险管理", "风险评估"),
            _el("social", "社会价值", "社会效益", "产业价值", "社会与产业价值", "可持续发展", "社会责任"),
            _el("appendix", "附录", "附录和参考资料", "参考资料", "参考文献", "附件"),
        ],
    },
    "digital": {
        **_BASE,
        "name": "数字营销赛道", "docName": "数字营销方案",
        "chartMin": 20, "sourceMin": 15, "numericDensity": 15,
        "budgetRequired": True, "kpiRequired": True,
        "requiredElements": [
            _el("toc", "目录"),
            _el("summary", "概要", "执行摘要", "摘要", "方案概要"),
            _el("intro", "引言", "前言", "绪论", "背景", "命题解读"),
            _el("company", "目标企业分析", "企业分析", "品牌分析", "企业现状", "企业与产品分析", "目标企业"),
            _el("market", "市场分析", "市场与用户分析", "用户分析", "行业分析", "市场洞察", "用户洞察"),
            _el("strategy", "数字营销策略制定", "数字营销策略", "营销策略", "策略制定", "传播策略", "内容策略", "策略规划"),
            _el("plan", "执行计划与预算", "执行计划", "预算", "排期与预算", "实施计划", "执行与预算"),
            _el("effect", "数字营销效果", "效果评估", "营销效果", "效果预估", "KPI", "效果与评估", "数据效果"),
            _el("appendix", "附录和参考资料", "附录", "参考资料", "参考文献", "短视频脚本", "附件"),
        ],
    },
    "case": {
        **_BASE,
        "name": "案例分析赛道", "docName": "案例分析报告",
        "pageLimit": 50, "wordMin": 3000, "chartMin": 12, "sourceMin": 12, "numericDensity": 15,
        "competitorMatrixRequired": False, "primaryResearchRequired": False, "financialModelRequired": True,
        "requiredElements": [
            _el("toc", "目录"),
            _el("summary", "概要", "执行摘要", "摘要"),
            _el("intro", "引言", "前言", "绪论"),
            _el("case", "案例简介", "企业简介", "案例企业", "公司简介", "案例背景", "企业概况"),
            _el("environment", "企业内外部环境分析", "环境分析", "内外部环境", "外部环境分析", "内部环境分析", "财务分析", "本企业财务分析", "SWOT"),
            _el("bottleneck", "发展瓶颈分析", "瓶颈分析", "问题分析", "策略分析", "问题与挑战", "发展瓶颈", "关键问题"),
            _el("conclusion", "结论与启示", "结论", "启示", "策略计划", "建议与启示", "对策建议", "实施计划"),
            _el("appendix", "附录与参考资料", "附录", "参考资料", "参考文献", "知识点列表", "附件"),
        ],
    },
    "mice": {
        **_BASE,
        "name": "商务会奖旅游策划竞赛", "docName": "会奖策划书",
        "chartMin": 15, "sourceMin": 10, "numericDensity": 12,
        "budgetRequired": True, "kpiRequired": True, "competitorMatrixRequired": False, "primaryResearchRequired": False,
        "budgetRange": [2000000, 4000000],
        "requiredElements": [
            _el("toc", "目录"),
            _el("client", "客户公司介绍", "客户介绍", "客户公司", "客户分析", "客户企业介绍", "委托方介绍"),
            _el("purpose", "方案目的", "活动目的", "方案目标", "目的与目标", "策划目的"),
            _el("feasibility", "方案可行性分析", "可行性分析", "可行性", "目的地分析", "环境分析"),
            _el("theme", "方案主题", "活动主题", "主题创意", "主题设计", "主题与创意"),
            _el("itinerary", "具体行程设计与规划", "行程设计", "行程安排", "日程安排", "行程规划", "活动安排", "日程"),
            _el("contingency", "突发事件预案", "应急预案", "危机处理", "风险预案", "突发事件", "风险与应急"),
            _el("budget", "费用预算", "预算", "财务预算", "预算明细", "费用明细"),
            _el("features", "方案特色", "方案亮点", "特色与亮点", "创新点", "亮点"),
            _el("effect", "方案预计效果", "预期效果", "预计效果", "效果评估", "效果预估"),
        ],
    },
}

TRACK_IDS = list(TRACK_DEFAULTS)


def default_config(track: str) -> Dict[str, Any]:
    """返回赛道默认配置的深拷贝；未知赛道退化为 brand 并在 `_fallbackFrom` 标记。"""
    cfg = copy.deepcopy(TRACK_DEFAULTS.get(track) or TRACK_DEFAULTS["brand"])
    cfg["track"] = track
    if track not in TRACK_DEFAULTS:
        cfg["_fallbackFrom"] = "brand"
    return cfg


def find_track_file(track: str, workspace: Optional[Path], filename: str = "checks.json") -> Optional[Path]:
    """依次在 工作区/tracks、工作区/agent/tracks、仓库 agent/tracks 下找赛道文件。"""
    cands: List[Path] = []
    if workspace is not None:
        ws = Path(workspace)
        cands += [ws / "tracks" / track / filename, ws / "agent" / "tracks" / track / filename]
    here = Path(__file__).resolve()
    # bizkit/ → toolkit/ → agent/ → agent/tracks/<id>/
    cands.append(here.parents[2] / "tracks" / track / filename)
    for c in cands:
        if c.exists():
            return c
    return None


def _expand_adlaw_term(term: str) -> List[str]:
    """单字 “最” → 具体组合表；其余原样。“第一” 的序数排除在 `adlaw_patterns()` 里用正则处理。"""
    t = str(term).strip()
    if not t:
        return []
    if t == "最":
        return [x for x in ADLAW_TERMS_DEFAULT if x.startswith("最")]
    return [t]


def expand_adlaw_terms(terms: Any) -> List[str]:
    out: List[str] = []
    for t in terms if isinstance(terms, list) else []:
        for x in _expand_adlaw_term(t):
            if x not in out:
                out.append(x)
    return out


def adlaw_patterns(terms: List[str]) -> List[str]:
    """把广告法用语转成正则：“第一” 排除序数后缀，其余按字面匹配。"""
    pats: List[str] = []
    for t in terms:
        if t == "第一":
            pats.append(r"第一(?!" + ADLAW_ORDINAL_SUFFIX + r")")
        else:
            pats.append(re.escape(t))
    return pats


def _as_list(v: Any) -> List[Any]:
    return list(v) if isinstance(v, (list, tuple)) else []


def _union(*lists: Any) -> List[Any]:
    out: List[Any] = []
    for lst in lists:
        for x in _as_list(lst):
            if x and x not in out:
                out.append(x)
    return out


def normalize_config(cfg: Dict[str, Any], override: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """让两套字段名互相补齐；`override` 为赛道文件原文（用于判断哪一套是作者写的，以它为准）。"""
    ov = override or {}
    # 1. 扁平字段别名
    for src, dst in FIELD_ALIASES.items():
        if src in ov and ov[src] is not None:
            cfg[dst] = ov[src]
        if dst in cfg and src not in cfg:
            cfg[src] = cfg[dst]
        elif src in cfg and dst not in cfg:
            cfg[dst] = cfg[src]
    # 2. pitch{}：赛道文件的 pitch 优先，其次扁平字段；讲稿字数区间缺省按 分钟 × 190 × [0.87, 1.03] 推算
    minutes = int(cfg.get("pitchMinutes") or 10)
    pitch: Dict[str, Any] = {
        "minutes": minutes,
        "deckPages": _as_list(cfg.get("deckPages")) or [16, 20],
        "deckPagesAppendixMax": cfg.get("deckPagesAppendixMax", 4),
        "qaMin": cfg.get("qaMin", 40),
        "qaWarn": cfg.get("qaWarn", 20),
        "qaDimensions": _as_list(cfg.get("qaDimensions")),
    }
    file_pitch: Dict[str, Any] = {}
    for src in (cfg.get("pitch"), ov.get("pitch")):
        if isinstance(src, dict):
            file_pitch.update({k: v for k, v in src.items() if v is not None})
    pitch.update(file_pitch)
    if _as_list(file_pitch.get("scriptWords")):
        pitch["scriptWords"] = list(file_pitch["scriptWords"])
    elif _as_list(cfg.get("scriptWords")):
        pitch["scriptWords"] = list(cfg["scriptWords"])
    else:
        m = int(pitch.get("minutes") or minutes)
        pitch["scriptWords"] = [round(m * 165), round(m * 195)]
    cfg["pitch"] = pitch
    cfg["pitchMinutes"] = pitch.get("minutes", minutes)
    cfg["deckPages"] = pitch.get("deckPages")
    cfg["deckPagesAppendixMax"] = pitch.get("deckPagesAppendixMax", 4)
    cfg["scriptWords"] = pitch.get("scriptWords")
    cfg["qaMin"] = pitch.get("qaMin", 40)
    cfg["qaWarn"] = pitch.get("qaWarn", cfg.get("qaWarn", 20))
    cfg["qaDimensions"] = pitch.get("qaDimensions") or []
    # 3. requiredEvidence{} ↔ *Required 开关
    if isinstance(ov.get("requiredEvidence"), dict):
        ev: Dict[str, Any] = dict(ov["requiredEvidence"])
    elif isinstance(cfg.get("requiredEvidence"), dict):
        ev = dict(cfg["requiredEvidence"])
    else:
        ev = {}
        for check_id, flag in EVIDENCE_FLAGS.items():
            if cfg.get(flag):
                ev[check_id] = True
        if cfg.get("financialModelRequired"):
            ev.setdefault("financing_consistency", True)
            ev.setdefault("revenue_shape", True)
        if cfg.get("budgetRange"):
            ev.setdefault("budget_in_range", list(cfg["budgetRange"]))
        ev.setdefault("risk_matrix_min", cfg.get("riskMatrixMin", 6))
    cfg["requiredEvidence"] = ev
    for check_id, flag in EVIDENCE_FLAGS.items():
        cfg[flag] = bool(ev.get(check_id))
    if isinstance(ev.get("budget_in_range"), (list, tuple)) and len(ev["budget_in_range"]) == 2:
        cfg["budgetRange"] = list(ev["budget_in_range"])
    elif cfg.get("budgetRange") and "budget_in_range" not in ev:
        ev["budget_in_range"] = list(cfg["budgetRange"])
    if isinstance(ev.get("risk_matrix_min"), (int, float)):
        cfg["riskMatrixMin"] = int(ev["risk_matrix_min"])
    # 4. 列表类：禁语以赛道文件为准；广告法用语展开单字；白名单取并集（只会减少误报）
    cfg["forbiddenTerms"] = _as_list(cfg.get("forbiddenTerms")) or list(FORBIDDEN_TERMS_DEFAULT)
    cfg["advertisingLawTerms"] = expand_adlaw_terms(cfg.get("advertisingLawTerms")) or list(ADLAW_TERMS_DEFAULT)
    cfg["advertisingLawWhitelist"] = _union(ADLAW_WHITELIST_DEFAULT, cfg.get("advertisingLawWhitelist"))
    cfg["weakTerms"] = _as_list(cfg.get("weakTerms")) or list(WEAK_TERMS_DEFAULT)
    cfg["weakTermsPer100"] = float(cfg.get("weakTermsPer100") or 0.1)
    cfg["anonymityWhitelist"] = _union(ANONYMITY_WHITELIST_DEFAULT, cfg.get("anonymityWhitelist"))
    cfg["anonymityRegex"] = cfg.get("anonymityRegex") or ANONYMITY_REGEX_DEFAULT
    cfg["numberKeys"] = _as_list(cfg.get("numberKeys"))
    if not isinstance(cfg.get("numericDensityChapters"), dict):
        cfg["numericDensityChapters"] = {}
    # 5. 检查器用到的其余字段：示例数据标记覆盖率阈值、封面必填字段（默认取 naming.json，调用方可再覆盖）、按变体的开关表
    try:
        cfg["syntheticMarkerCoverage"] = float(cfg.get("syntheticMarkerCoverage", 1.0))
    except (TypeError, ValueError):
        cfg["syntheticMarkerCoverage"] = 1.0
    cfg["coverFields"] = _as_list(cfg.get("coverFields"))
    for k in ("requiredEvidenceByVariant", "coverFieldsByVariant"):
        if not isinstance(cfg.get(k), dict):
            cfg[k] = {}
    return cfg


def apply_variant(cfg: Dict[str, Any], variant: Optional[str]) -> Dict[str, Any]:
    """按变体（案例赛道 accounting | retail）把 `requiredEvidenceByVariant` / `coverFieldsByVariant` 合并进生效配置。

    变体开关覆盖同名的 `requiredEvidence` 项；`word_min` / `deck_min_pages` 两个特殊键同时写回 `wordMin` / `pitch.deckPages[0]`。
    未知或空变体原样返回。
    """
    v = str(variant or "").strip()
    cfg["variant"] = v or None
    if not v:
        return cfg
    by_ev = cfg.get("requiredEvidenceByVariant") or {}
    if isinstance(by_ev.get(v), dict):
        ev = dict(cfg.get("requiredEvidence") or {})
        ev.update({k: x for k, x in by_ev[v].items() if x is not None})
        cfg["requiredEvidence"] = ev
        if isinstance(ev.get("word_min"), (int, float)):
            cfg["wordMin"] = max(int(cfg.get("wordMin") or 0), int(ev["word_min"]))
        if isinstance(ev.get("deck_min_pages"), (int, float)):
            pitch = dict(cfg.get("pitch") or {})
            pages = list(pitch.get("deckPages") or cfg.get("deckPages") or [16, 20])
            pages[0] = max(int(pages[0]), int(ev["deck_min_pages"]))
            pages[1] = max(int(pages[1]), pages[0])
            pitch["deckPages"] = pages
            cfg["pitch"] = pitch
            cfg["deckPages"] = pages
        for check_id, flag in EVIDENCE_FLAGS.items():
            cfg[flag] = bool(ev.get(check_id))
    by_cf = cfg.get("coverFieldsByVariant") or {}
    if isinstance(by_cf.get(v), list):
        cfg["coverFields"] = list(by_cf[v])
    return cfg


def load_track_config(track: str, workspace: Optional[Union[str, Path]] = None) -> Dict[str, Any]:
    """内置默认 ← tracks/<id>/checks.json 覆盖（字段级；requiredElements 整体替换），再 `normalize_config()` 补齐别名。"""
    cfg = default_config(track)
    path = find_track_file(track, Path(workspace) if workspace else None)
    cfg["_source"] = "builtin"
    override: Dict[str, Any] = {}
    if path is not None:
        try:
            loaded = read_json(path)
        except Exception as e:  # 赛道文件损坏不应让检查器崩溃，退回默认并记录
            cfg["_warning"] = f"读取 {path} 失败：{e}"
            return normalize_config(cfg)
        if isinstance(loaded, dict):
            override = loaded
            for k, v in override.items():
                if v is None:
                    continue
                if k == "forbiddenTerms" and isinstance(v, list):
                    cfg[k] = list(v)
                elif k == "requiredElements" and isinstance(v, list):
                    cfg[k] = [x if isinstance(x, dict) else {"id": str(x), "name": str(x), "aliases": []} for x in v]
                else:
                    cfg[k] = v
            cfg["_source"] = str(path)
    return normalize_config(cfg, override)
