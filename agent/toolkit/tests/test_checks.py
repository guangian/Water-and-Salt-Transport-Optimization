"""bizkit.checks（机器检查器，契约 §7.7）测试。

运行：cd agent/toolkit && python -m pytest -q tests/test_checks.py

策略：
1. examples/ws_demo3 是一个带“故意缺陷”的迷你品牌赛道工作区：一处匿名违规、一处禁语、一处广告法用语、
   一处缺「示例数据」标记的样本量、讲稿里一处与账本不一致的预算；同时含竞品矩阵 / KPI 表 / 定位陈述 / 命题回应表等正向证据。
   逐项断言这些缺陷被命中、正向证据被识别、输出文件与 summary 结构符合契约。
2. 缺文件不崩：空工作区、缺 numbers.json / 06_路演 的工作区全部 skipped 而非异常；坏 JSON 也只是 skipped。
3. 对 examples/ws_demo（brand）与 ws_demo2（startup）整体运行不崩、无 errors、硬门槛结构完整。
4. --chapter 单章模式只跑六项且不写文件；CLI 入口 main() 最后一行是 JSON。
5. 变体：案例赛道 apply_variant 合并 requiredEvidenceByVariant / coverFieldsByVariant。
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import pytest

TOOLKIT = Path(__file__).resolve().parents[1]
if str(TOOLKIT) not in sys.path:
    sys.path.insert(0, str(TOOLKIT))

from bizkit import checks, track_defaults as td  # noqa: E402

EXAMPLES = TOOLKIT / "examples"
PASS, WARN, FAIL, SKIP = checks.PASS, checks.WARN, checks.FAIL, checks.SKIP


def _regen(dst: Path, args: List[str]) -> None:
    """仓库只保存示例工作区的源文件（charts.json / assumptions.json / slides.json），派生产物按需重新生成。"""
    env = dict(os.environ, PYTHONPATH=str(EXAMPLES.parent), PYTHONUTF8="1")
    subprocess.run([sys.executable, "-m", *args], cwd=dst, env=env, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def _copy_ws(name: str, dst: Path) -> Path:
    shutil.copytree(EXAMPLES / name, dst, dirs_exist_ok=True, ignore=shutil.ignore_patterns("__pycache__"))
    if (dst / "04_图表" / "charts.json").exists() and not (dst / "04_图表" / "index.json").exists():
        _regen(dst, ["bizkit.charts", "04_图表/charts.json"])
    if (dst / "03_财务" / "assumptions.json").exists() and not (dst / "03_财务" / "财务模型.xlsx").exists():
        _regen(dst, ["bizkit.finmodel", "03_财务/assumptions.json", "--out", "03_财务"])
    if (dst / "06_路演" / "slides.json").exists() and not (dst / "产出" / "路演.pptx.preview.json").exists() and name != "ws_demo3":
        _regen(dst, ["bizkit.pptx_render", "06_路演/slides.json"])
    return dst


def _by_id(result: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    return {c["id"]: c for c in result["checks"]}


# --------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def demo3(tmp_path_factory: pytest.TempPathFactory) -> Dict[str, Any]:
    """ws_demo3 的干净副本上跑一次第 1 轮全量检查，供多个测试共享。"""
    ws = _copy_ws("ws_demo3", tmp_path_factory.mktemp("ws_demo3"))
    result = checks.run(ws, track="brand", round_=1)
    return {"ws": ws, "result": result, "checks": _by_id(result)}


# --------------------------------------------------------------------------- 1. 迷你工作区：故意缺陷与正向证据
def test_demo3_runs_clean(demo3: Dict[str, Any]) -> None:
    r = demo3["result"]
    assert r["ok"] is True
    assert r["errors"] == []
    assert r["track"] == "brand" and r["round"] == 1
    s = r["summary"]
    assert s["degraded"] is False
    assert s["pass"] + s["warn"] + s["fail"] + s["skipped"] == len(r["checks"])
    assert all(c["status"] in (PASS, WARN, FAIL, SKIP) for c in r["checks"])
    assert len({c["id"] for c in r["checks"]}) == len(r["checks"]), "检查项 id 不得重复"


def test_demo3_writes_contract_files(demo3: Dict[str, Any]) -> None:
    ws: Path = demo3["ws"]
    jpath = ws / "07_评审" / "第1轮" / "机器检查.json"
    mpath = ws / "07_评审" / "第1轮" / "机器检查.md"
    assert jpath.is_file() and mpath.is_file()
    data = json.loads(jpath.read_text(encoding="utf-8"))
    for key in ("checks", "hardGates", "crossDoc", "summary"):
        assert key in data
    for c in data["checks"]:
        assert {"id", "value", "threshold", "status", "detail"} <= set(c)
    md = mpath.read_text(encoding="utf-8")
    assert "## 硬门槛" in md and "## 检查项" in md and "## 跨材料数字一致性" in md
    assert "280" in md, "不一致清单应列出讲稿里的 280 万元"
    assert data["summary"]["degraded"] is False


def test_demo3_anonymity_hit(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]["anonymity_hits"]
    assert c["status"] == FAIL and c["value"] == 1
    loc = c["locations"][0]
    assert loc["file"] == "05_正文/02_品牌现状分析.md" and loc["term"] == "大学"
    assert "华中某大学" in loc["text"]
    # 白名单：正文里的「光谷大学城」不算命中
    assert all("大学城" not in h["text"] or "某大学" in h["text"] for h in c["locations"])


def test_demo3_forbidden_terms(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]["forbidden_terms"]
    assert c["value"] == 1 and c["status"] == WARN
    assert c["locations"][0]["term"] == "蓝海" and "蓝海" in c["locations"][0]["text"]
    # needSource 词「据统计」同句带 [3] 来源标注，不应命中
    assert all(h["term"] != "据统计" for h in c["locations"])


def test_demo3_adlaw_hard_gate(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]["forbidden_terms_adlaw"]
    assert c["status"] == FAIL and c["value"] == 1
    assert c["locations"][0]["term"] == "最佳"
    # 「第一阶段」「第一，」「第 1 位」是序数 / 列举用法，不是广告法用语
    assert all(h["term"] != "第一" for h in c["locations"])
    gates = {g["id"]: g for g in demo3["result"]["hardGates"]}
    assert gates["forbidden_terms_adlaw"]["passed"] is False
    assert gates["anonymity_hits"]["passed"] is False
    assert gates["cover_fields"]["passed"] is True
    assert demo3["result"]["summary"]["hardGatePassed"] is False


def test_demo3_positive_evidence(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]
    assert c["has_competitor_matrix"]["status"] == PASS and c["has_competitor_matrix"]["value"] >= 5
    assert "主要竞品对比" in c["has_competitor_matrix"]["detail"]
    assert c["kpi_table"]["status"] == PASS and c["kpi_table"]["value"] >= 3
    assert c["positioning_statement"]["status"] == PASS
    assert c["big_idea"]["status"] == PASS
    assert c["kol_tiers"]["status"] == PASS
    assert c["risk_matrix"]["status"] == PASS and c["risk_matrix"]["value"] >= 6
    assert c["budget_recalc"]["status"] == PASS and c["budget_recalc"]["value"] == 0
    assert "机动金 8.3%" in c["budget_recalc"]["detail"]
    assert c["has_primary_research"]["status"] == PASS and c["has_primary_research"]["value"] == 412
    assert c["appendix_questionnaire"]["status"] == PASS and c["appendix_questionnaire"]["value"] == 26
    assert c["has_sensitivity"]["status"] == PASS
    assert c["title_pattern"]["status"] == PASS
    assert c["honor_ratio"]["status"] == PASS and c["honor_ratio"]["value"] < 0.05
    # 命题案例（project.subject.isAssignedCase=true）且赛道 brief_response_table="if_assigned" → 检查并命中
    assert c["brief_response_table"]["status"] == PASS and c["brief_response_table"]["value"] == 5


def test_demo3_sources(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]
    assert c["source_count"]["value"] == 6 and c["source_count"]["status"] == FAIL
    assert c["source_official_ratio"]["value"] == 0.5 and c["source_official_ratio"]["status"] == PASS
    assert c["source_freshness"]["value"] == 1.0
    sp = c["source_pairing"]
    assert sp["status"] == WARN and len(sp["locations"]) == 1
    assert "fig_stores" in sp["locations"][0]["text"] or "门店焕新前后对比" in sp["locations"][0]["text"]


def test_demo3_synthetic_marker_missing_once(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]["synthetic_marker_coverage"]
    assert c["status"] == WARN
    assert c["occurrences"] == 9 and len(c["locations"]) == 1
    miss = c["locations"][0]
    assert miss["key"] == "sample_size" and miss["file"] == "05_正文/03_品牌定位与传播策略.md"
    assert "412 份有效问卷" in miss["text"]


def test_demo3_cross_doc_finds_script_mismatch(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]["cross_doc_consistency"]
    cross = demo3["result"]["crossDoc"]
    assert c["status"] == FAIL and c["value"] == 1
    assert cross["consistent"] is False
    assert set(cross["keys"]) == {"budget_total", "sample_size", "kpi_exposure", "roi"}
    (m,) = cross["mismatches"]
    assert m["key"] == "budget_total" and m["expected"] == 300
    assert m["script"] == 280.0 and m["book"] == 300.0 and m["deck"] == 300.0 and m["qa"] == 300.0
    assert m["locations"][0]["file"] == "06_路演/讲稿.md" and "280" in m["locations"][0]["text"]
    gates = {g["id"]: g for g in demo3["result"]["hardGates"]}
    assert gates["cross_doc_consistency"]["passed"] is False
    assert gates["cross_doc_consistency"]["mismatches"][0]["key"] == "budget_total"


def test_demo3_pitch_checks(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]
    st = c["script_timing"]
    assert st["status"] == PASS and 1650 <= st["value"] <= 1950 and st["paragraphs"] >= 8
    assert c["sentence_length"]["status"] == PASS and c["sentence_length"]["max"] <= 40
    qa = c["qa_bank"]
    assert qa["value"] == 12 and qa["status"] == FAIL
    assert qa["missingOwner"] == 0 and qa["missingThreePart"] == 0
    assert all(n == 2 for n in qa["perDimension"].values())
    sm = c["slide_metrics"]
    assert sm["status"] == FAIL and sm["metrics"]["bodyPages"] == 9 and sm["sub"]["pages"] == FAIL
    assert sm["metrics"]["overWordsPages"] == 0
    assert c["deck_font_min"]["status"] == SKIP  # 没有 产出/路演.render.json


def test_demo3_page_and_density(demo3: Dict[str, Any]) -> None:
    c = demo3["checks"]
    pc = c["page_count"]
    assert pc["status"] in (WARN, FAIL) and 0 < pc["value"] < 20 and "估算" in pc["detail"]
    nd = c["numeric_density"]
    assert nd["status"] == PASS and nd["value"] > 15
    assert c["weak_terms"]["value"] == 0.0  # 「有效问卷 / 有效率」属专业搭配不计弱词
    assert c["struct_coverage"]["status"] == FAIL and 0 < c["struct_coverage"]["value"] < 1


# --------------------------------------------------------------------------- 2. 缺文件 → skipped，不崩
def test_empty_workspace_all_skipped(tmp_path: Path) -> None:
    r = checks.run(tmp_path, track="brand", round_=2)
    assert r["ok"] is True and r["errors"] == []
    assert r["summary"]["skipped"] == len(r["checks"]) and r["summary"]["fail"] == 0
    assert r["crossDoc"]["consistent"] is True and r["crossDoc"]["checked"] == 0
    assert all(g["passed"] for g in r["hardGates"]) and all(g.get("skipped") for g in r["hardGates"])
    assert (tmp_path / "07_评审" / "第2轮" / "机器检查.json").is_file()


def test_partial_workspace_and_bad_json(tmp_path: Path) -> None:
    """只有正文、无 numbers.json / 06_路演，且 参考文献.json 是坏 JSON：依赖项 skipped，其余照常。"""
    _copy_ws("ws_demo3", tmp_path)
    shutil.rmtree(tmp_path / "06_路演")
    (tmp_path / "03_财务" / "numbers.json").unlink()
    (tmp_path / "01_调研" / "参考文献.json").write_text("{not json", encoding="utf-8")
    r = checks.run(tmp_path, track="brand", round_=None)
    assert r["ok"] is True and r["errors"] == [] and r["summary"]["degraded"] is False
    c = _by_id(r)
    for cid in ("script_timing", "sentence_length", "qa_bank", "slide_metrics", "source_count", "source_official_ratio",
                "synthetic_marker_coverage", "cross_doc_consistency"):
        assert c[cid]["status"] == SKIP, cid
    assert c["anonymity_hits"]["value"] == 1  # 正文仍被检查
    assert c["budget_recalc"]["status"] == PASS  # 账本缺失时以正文预算表自洽为准
    assert (tmp_path / "07_评审" / "即时" / "机器检查.json").is_file()


def test_unknown_track_falls_back(tmp_path: Path) -> None:
    r = checks.run(tmp_path, track="no_such_track", round_=None, write=False)
    assert r["ok"] is True and r["files"] is None
    assert r["config"].get("warning")


# --------------------------------------------------------------------------- 3. 对已有示例工作区整体运行不崩
@pytest.mark.parametrize("name,track", [("ws_demo", "brand"), ("ws_demo2", "startup")])
def test_examples_do_not_crash(tmp_path: Path, name: str, track: str) -> None:
    _copy_ws(name, tmp_path)
    r = checks.run(tmp_path, track=None, round_=1)  # 赛道从 project.json 推断
    assert r["ok"] is True and r["errors"] == [] and r["track"] == track
    assert r["summary"]["degraded"] is False
    gate_ids = {g["id"] for g in r["hardGates"]}
    assert {"anonymity_hits", "page_limit", "cross_doc_consistency"} <= gate_ids
    c = _by_id(r)
    assert set(c) >= {"struct_coverage", "page_count", "word_count", "anonymity_hits", "numeric_density", "source_count", "forbidden_terms",
                      "forbidden_terms_adlaw", "weak_terms", "synthetic_marker_coverage", "cross_doc_consistency"}
    if name == "ws_demo":
        if ".pdf" in c["page_count"]["detail"]:
            assert c["page_count"]["pages"] == 25  # 产出/*.pdf 真页数（pypdf）；干净克隆没有 PDF 时走估算
        else:
            assert c["page_count"]["pages"] > 0
        assert c["has_competitor_matrix"]["status"] == PASS
        assert c["cross_doc_consistency"]["status"] == PASS and c["cross_doc_consistency"]["checked"] == 2
    else:
        assert c["three_statement_check"]["status"] == PASS  # 财务模型.xlsx「检查」页
        assert c["financing_consistency"]["status"] == PASS
        assert c["cross_doc_consistency"]["status"] == FAIL  # ws_demo2 的 PPT 里本就有两处与账本不一致
        assert {m["key"] for m in r["crossDoc"]["mismatches"]} == {"breakeven_month", "post_money"}


# --------------------------------------------------------------------------- 4. 单章模式与 CLI
def test_chapter_mode(tmp_path: Path) -> None:
    _copy_ws("ws_demo3", tmp_path)
    r = checks.run_chapter(tmp_path, "brand", "05_正文/02_品牌现状分析.md")
    assert r["ok"] is True and r["mode"] == "chapter" and r["files"] is None
    assert [c["id"] for c in r["checks"]] == checks.CHAPTER_CHECK_IDS
    c = _by_id(r)
    assert c["anonymity_hits"]["value"] == 1 and c["forbidden_terms"]["value"] == 1
    assert c["forbidden_terms_adlaw"]["value"] == 0  # 「最佳」在第 3 章，不在本章
    assert c["source_pairing"]["status"] == PASS  # 图经 04_图表/index.json 的 source 配对
    assert c["numeric_density"]["chapters"][0]["threshold"] == 25  # numericDensityChapters 对 02_品牌现状分析 的专项阈值
    assert not (tmp_path / "07_评审").exists()


def test_chapter_mode_missing_file(tmp_path: Path) -> None:
    r = checks.run_chapter(tmp_path, "brand", "05_正文/不存在.md")
    assert r["ok"] is True and all(c["status"] == SKIP for c in r["checks"])


def test_cli_main_prints_json_last_line(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    _copy_ws("ws_demo3", tmp_path)
    code = checks.main(["--track", "brand", "--workspace", str(tmp_path), "--chapter", "05_正文/03_品牌定位与传播策略.md"])
    out = capsys.readouterr().out.strip().splitlines()
    data = json.loads(out[-1])
    assert code == 0 and data["ok"] is True and data["mode"] == "chapter"
    assert _by_id(data)["forbidden_terms_adlaw"]["value"] == 1

    code = checks.main(["--track", "brand", "--workspace", str(tmp_path), "--round", "3"])
    out = capsys.readouterr().out.strip().splitlines()
    data = json.loads(out[-1])
    assert code == 0 and data["round"] == 3
    assert (tmp_path / "07_评审" / "第3轮" / "机器检查.md").is_file()


def test_cli_never_raises(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code = checks.main(["--workspace", str(tmp_path / "nowhere"), "--no-write"])
    data = json.loads(capsys.readouterr().out.strip().splitlines()[-1])
    assert code == 0 and data["ok"] is True  # 空目录 → 全部 skipped，仍是合法结果


# --------------------------------------------------------------------------- 5. 变体与配置
def test_case_track_variant_merge() -> None:
    base = td.load_track_config("case", None)
    retail = td.apply_variant(td.load_track_config("case", None), "retail")
    assert retail["variant"] == "retail"
    assert retail["coverFields"] and retail["coverFields"] != base.get("coverFields")
    assert retail["requiredEvidence"] != base["requiredEvidence"]
    acct = td.apply_variant(td.load_track_config("case", None), "accounting")
    assert acct["requiredEvidence"] != retail["requiredEvidence"]


def test_variant_changes_enabled_checks(tmp_path: Path) -> None:
    (tmp_path / "project.json").write_text(json.dumps({"track": "case", "subject": {"type": "case"}}, ensure_ascii=False), encoding="utf-8")
    r1 = checks.run(tmp_path, track="case", variant="accounting", write=False)
    r2 = checks.run(tmp_path, track="case", variant="retail", write=False)
    assert r1["variant"] == "accounting" and r2["variant"] == "retail"
    assert r1["ok"] and r2["ok"]
    assert {c["id"] for c in r1["checks"]} != {c["id"] for c in r2["checks"]}


# --------------------------------------------------------------------------- 6. 解析器单元
def test_parse_doc_tables_and_figures() -> None:
    text = "# 1 标题\n\n表：示例表\n\n| 维度 | A | B |\n|---|---|---|\n| x | 1 | 2 |\n| y | 3 | 4 |\n\n来源：测试\n\n![图 1](../04_图表/fig_a.png)\n"
    d = checks.parse_doc("05_正文/01_t.md", text)
    assert len(d.tables) == 1 and d.tables[0].title == "示例表" and d.tables[0].nrows == 2 and d.tables[0].heading == "1 标题"
    assert len(d.figures) == 1 and d.figures[0].fig_id == "fig_a"
    assert checks._count_words("中文 abc 123 中") == 5


def test_adlaw_patterns_skip_ordinals() -> None:
    cfg = td.load_track_config("brand", None)
    pats = checks._adlaw_patterns(cfg)
    units = [("f.md", 1, "第一阶段完成后，第一，复盘；第一年收入"), ("f.md", 2, "我们是行业第一的最佳品牌"), ("f.md", 3, "100% 天然，100% 完成率")]
    hits = checks._scan(units, pats, cfg.get("advertisingLawWhitelist") or [])
    terms = sorted(h["term"] for h in hits)
    assert "第一" not in [h["term"] for h in hits if h["line"] == 1]
    assert {"最佳", "行业第一"} <= set(terms) or {"最佳", "第一"} <= set(terms)
    assert sum(1 for h in hits if h["line"] == 3) == 1  # 只有「100% 天然」算宣传搭配


# --------------------------------------------------------------------------- 7. 文风三项：套话 / 示例数据标注（章节级）/ 无源关键数字
def _mini_ws(root: Path, chapter: str, project: Any = None, numbers: Any = None, script: Any = None, qa: Any = None) -> Path:
    """只有一章正文（可选 project.json / numbers.json / 讲稿 / 题库）的最小工作区。"""
    (root / "05_正文").mkdir(parents=True, exist_ok=True)
    (root / "05_正文" / "01_市场分析.md").write_text(chapter, encoding="utf-8")
    if project is not None:
        (root / "project.json").write_text(json.dumps(project, ensure_ascii=False), encoding="utf-8")
    if numbers is not None:
        (root / "03_财务").mkdir(exist_ok=True)
        (root / "03_财务" / "numbers.json").write_text(json.dumps(numbers, ensure_ascii=False), encoding="utf-8")
    if script is not None or qa is not None:
        (root / "06_路演").mkdir(exist_ok=True)
        if script is not None:
            (root / "06_路演" / "讲稿.md").write_text(script, encoding="utf-8")
        if qa is not None:
            (root / "06_路演" / "答辩题库.md").write_text(qa, encoding="utf-8")
    return root


_CLICHE_CHAPTER = (
    "# 2 市场分析\n\n"
    "随着电解质饮料市场的不断发展，功能性零食迎来空前的机遇。赋能渠道是关键抓手，形成生态闭环，引领行业。\n\n"
    "本方案可以说明问题所在；综上所述，强强联合、全方位多层次打造标杆，助力县域经济高质量发展。\n\n"
    "市场规模 120 亿元[1]，不难发现，一站式服务前景广阔。\n"
)


def test_cliche_terms_hits_and_wildcards(tmp_path: Path) -> None:
    r = checks.run(_mini_ws(tmp_path, _CLICHE_CHAPTER, project={"track": "brand"}), track="brand", write=False)
    c = _by_id(r)["cliche_terms"]
    assert c["status"] == FAIL and c["value"] >= 9 and c["hits"] == len(c["locations"]) >= 10
    terms = [h["term"] for h in c["locations"]]
    assert "随着电解质饮料市场的不断发展" in terms  # 「随着…发展」模板句（… 为通配）
    assert "助力县域经济高质量发展" in terms
    assert {"赋能", "抓手", "生态闭环", "引领行业", "综上所述", "强强联合", "打造标杆", "不难发现", "一站式", "前景广阔"} <= set(terms)
    assert "可以说" not in terms  # 「可以说明」在 clicheWhitelist 里
    assert all(h["file"] == "05_正文/01_市场分析.md" and h["line"] in (3, 5, 7) and h["text"] for h in c["locations"])
    assert "高频" in c["detail"] and "/万字" in c["detail"]
    # 不进硬门槛
    assert "cliche_terms" not in {g["id"] for g in r["hardGates"]}


def test_cliche_terms_clean_and_thresholds(demo3: Dict[str, Any], tmp_path: Path) -> None:
    c = demo3["checks"]["cliche_terms"]
    assert c["status"] == PASS and c["value"] == 0.0 and c["hits"] == 0 and c["words"] > 3000  # 正文 + 讲稿一起计字数
    # 每万字 ≤ 3 pass；(3, 9) warn；≥ 9 fail —— 用 1 处套话 + 不同字数验证三档
    filler = "本方案第一年目标销售额 4,860 万元，按客单价 16 元与复购模型测算。" * 60  # ≈ 1 800 字
    one_hit = "# 1 概要\n\n" + filler + "\n\n我们要赋能渠道。\n"
    r = checks.run(_mini_ws(tmp_path / "w", one_hit, project={"track": "brand"}), track="brand", write=False)
    c1 = _by_id(r)["cliche_terms"]
    assert c1["hits"] == 1 and c1["status"] == WARN and 3 < c1["value"] < 9
    r = checks.run(_mini_ws(tmp_path / "p", "# 1 概要\n\n" + filler * 3 + "\n\n我们要赋能渠道。\n", project={"track": "brand"}), track="brand", write=False)
    assert _by_id(r)["cliche_terms"]["status"] == PASS
    # 赛道 clicheWhitelist 可放过专业搭配；clicheTerms 可整体替换
    cfg = td.load_track_config("brand", None)
    assert "赋能" in cfg["clicheTerms"] and "可以说明" in cfg["clicheWhitelist"]
    pats = checks._cliche_patterns({"clicheTerms": ["在…背景下", {"term": r"深度\s*融合", "regex": True}]})
    hits = checks._scan([("f.md", 1, "在消费升级的背景下深度 融合；随着，市场发展")], pats, [])
    assert sorted(h["term"] for h in hits) == ["在消费升级的背景下", "深度 融合"]  # 通配不跨标点（「随着，市场发展」不计）


_RESEARCH_CHAPTER = (
    "# 2 消费者调研\n\n"
    "团队在 6 个商圈发放问卷并完成深度访谈，细节见附录。\n\n"
    "调研显示，受访者对品牌的认知率为 71.2%。\n\n"
    "有效样本 412 份，其中 25—35 岁占 58.3%。\n"
)
_SYN_PROJECT = {"track": "brand", "params": {"researchMode": "synthetic"}}


def test_synthetic_research_markers_missing(tmp_path: Path) -> None:
    script = "[0:00–1:00] 主讲\n\n访谈中一位店长说顾客变老了。有效问卷 412 份（示例数据，须替换）。\n"
    qa = "### Q1【市场与数据】样本够吗\n- 结论：412 份问卷在 95% 置信水平误差 ±4.8%。\n"
    r = checks.run(_mini_ws(tmp_path, _RESEARCH_CHAPTER, project=_SYN_PROJECT, script=script, qa=qa), track="brand", write=False)
    c = _by_id(r)["synthetic_research_markers"]
    assert c["status"] == FAIL and c["value"] == 2 and c["files"] == 3
    assert c["missingFiles"] == ["05_正文/01_市场分析.md", "06_路演/答辩题库.md"]  # 讲稿同文件有标记 → 不缺
    locs = {(h["file"], h["line"]): h for h in c["locations"]}
    assert ("05_正文/01_市场分析.md", 3) not in locs  # 「发放问卷并完成深度访谈」是方法描述，不是结论
    assert locs[("05_正文/01_市场分析.md", 5)]["term"] == "调研显示" and "71.2%" in locs[("05_正文/01_市场分析.md", 5)]["text"]
    assert locs[("05_正文/01_市场分析.md", 7)]["term"] == "样本"
    assert locs[("06_路演/答辩题库.md", 2)]["term"] == "问卷"
    assert "缺标记：05_正文/01_市场分析.md（第 5、7 行）" in c["detail"]
    assert "synthetic_research_markers" not in {g["id"] for g in r["hardGates"]}


def test_synthetic_research_markers_marked_or_not_synthetic(tmp_path: Path, demo3: Dict[str, Any]) -> None:
    marked = _RESEARCH_CHAPTER.replace("调研显示", "【示例数据，须替换为真实调研结果】调研显示")
    r = checks.run(_mini_ws(tmp_path / "ok", marked, project=_SYN_PROJECT), track="brand", write=False)
    c = _by_id(r)["synthetic_research_markers"]
    assert c["status"] == PASS and c["value"] == 0 and c["files"] == 1 and c["locations"] == []
    # provided 模式且账本无 synthetic 条目 → skipped
    r = checks.run(_mini_ws(tmp_path / "prov", _RESEARCH_CHAPTER, project={"track": "brand", "params": {"researchMode": "provided"}}), track="brand", write=False)
    assert _by_id(r)["synthetic_research_markers"]["status"] == SKIP
    # 无 researchMode 但账本有 synthetic:true → 启用；示例数值只在带同类单位时触发（「第 32 类」商标类别不算 32%）
    numbers = {"entries": [{"key": "sample_size", "label": "样本", "value": 412, "unit": "份", "synthetic": True},
                           {"key": "research_scene_sport", "label": "运动场景占比", "value": 32, "unit": "%", "synthetic": True}]}
    ch = "# 6 资产保护\n\n商标按第 30 类核心、第 32 类防御注册，预算 412 元。\n\n运动后场景占 32%。\n"
    r = checks.run(_mini_ws(tmp_path / "led", ch, project={"track": "brand"}, numbers=numbers), track="brand", write=False)
    c = _by_id(r)["synthetic_research_markers"]
    assert c["status"] == FAIL and [(h["line"], h["term"]) for h in c["locations"]] == [(5, "账本 research_scene_sport")]
    # ws_demo3：正文 / 讲稿 / 题库全部带标记，只有 PPT 第 3 页「38.7% 的受访者认为…」没有标记
    d = demo3["checks"]["synthetic_research_markers"]
    assert d["status"] == FAIL and d["missingFiles"] == ["06_路演/slides.json#3"] and d["files"] >= 6
    assert d["locations"][0]["file"] == "06_路演/slides.json" and d["locations"][0]["line"] == 3 and "受访者" in d["locations"][0]["text"]


_NUMBERS_CHAPTER = (
    "# 2 市场分析\n\n"
    "## 2.1 规模\n\n"
    "2025 年中国功能糖果零售规模 120 亿元，较 2023 年增长 18.6%[1]。渗透率 3.2%，见图 2-1 与表 2-3，详见第 3 章第 12 页。\n\n"
    "行业均价 15 元/袋，头部品牌客单价 28.5 元，复购率 35%。\n\n"
    "Z 世代人群规模 2.6 亿人，年均消费 1,200 元。\n来源：团队整理自国家统计局。\n\n"
    "| 指标 | 数值 |\n|---|---|\n| 市场规模 | 120 亿元 |\n| 占比 | 45% |\n\n来源：表内来源行。\n\n"
    "按客单价 20 元与复购模型测算，首年销售额 4,860 万元。预算总额 300 万元，样本量 412 份。\n\n"
    "预计到场 3 万人次，粉丝 100 万以上。\n"
)
_LEDGER = {"entries": [{"key": "budget_total", "label": "预算", "value": 300, "unit": "万元"}, {"key": "revenue_y1", "label": "首年收入", "value": 4860, "unit": "万元"}]}


def test_unsourced_numbers_counts_and_exclusions(tmp_path: Path) -> None:
    r = checks.run(_mini_ws(tmp_path, _NUMBERS_CHAPTER, project={"track": "brand"}, numbers=_LEDGER), track="brand", write=False)
    c = _by_id(r)["unsourced_numbers"]
    # 计入：120 亿元 / 18.6% / 3.2%（[1] 覆盖全段）、15 元 / 28.5 元 / 35%（无源）、2.6 亿人 / 1,200 元（下一行来源行）、
    #       20 元 / 4,860 万元 / 300 万元（“测算”）、3 万人次 / 100 万（无源）；年份、图表号、章节页码、表格单元格、[1] 本身不计
    assert c["total"] == 13 and c["unsourced"] == 5
    assert c["status"] == FAIL and c["value"] == round(5 / 13, 4)
    assert [h["term"] for h in c["locations"]] == ["15 元", "28.5 元", "35%", "3 万人次", "100 万"]
    assert all(h["file"] == "05_正文/01_市场分析.md" and h["line"] in (7, 21) and len(h["text"]) <= 60 for h in c["locations"])
    assert "unsourced_numbers" not in {g["id"] for g in r["hardGates"]}


def test_unsourced_numbers_sources_ledger_and_demo3(tmp_path: Path, demo3: Dict[str, Any]) -> None:
    # 同一段补上文献号 / 账本 → 无源数下降到 warn；账本匹配允许千分位 / 小数位 / 万↔亿 / 20%↔0.2
    fixed = _NUMBERS_CHAPTER.replace("复购率 35%。", "复购率 35%[2]。").replace("粉丝 100 万以上。", "粉丝 100 万以上（示例数据）。")
    r = checks.run(_mini_ws(tmp_path, fixed, project={"track": "brand"}, numbers=_LEDGER), track="brand", write=False)
    c = _by_id(r)["unsourced_numbers"]
    assert c["unsourced"] == 0 and c["status"] == PASS and c["locations"] == []
    ledger = [(4860.0, "万元"), (0.2, "%"), (3547.0, "亿元"), (2.5, "倍")]
    assert checks._in_ledger(4860, "万元", ledger) and checks._in_ledger(0.486, "亿元", ledger) and checks._in_ledger(48600000, "元", ledger)
    assert checks._in_ledger(20, "%", ledger) and checks._in_ledger(3547.2, "亿元", ledger) and checks._in_ledger(2.50, "倍", ledger)
    assert not checks._in_ledger(3600, "亿元", ledger)
    # 读图 / 读表句与紧跟图表之后的段落由图表来源覆盖；一手调研口径（问卷 / 受访者）与「见第 N 章 / 见表」指向也算有源
    read = ("# 2 x\n\n![份额](../04_图表/fig_a.png)\n\n读图：头部占 47%，CR5 合计 93.2%。\n\n"
            "受访者中 58% 愿意多付 2 元。\n\n三阶段预算 126 万元、68 万元，见第 7 章明细。\n\n单场预算 2 万元。\n")
    r = checks.run(_mini_ws(tmp_path / "read", read, project={"track": "brand"}), track="brand", write=False)
    c = _by_id(r)["unsourced_numbers"]
    assert c["total"] == 7 and [h["term"] for h in c["locations"]] == ["2 万元"]  # 47% / 93.2% / 58% / 2 元 / 126 万元 / 68 万元 有源
    # 样例最多 20 条，但总数照记
    many = "# 1 x\n\n" + "\n\n".join(f"第 {i} 类客群客单价 {10 + i}.5 元。" for i in range(30)) + "\n"
    r = checks.run(_mini_ws(tmp_path / "many", many, project={"track": "brand"}), track="brand", write=False)
    c = _by_id(r)["unsourced_numbers"]
    assert c["unsourced"] == 30 and len(c["locations"]) == 20 and c["status"] == FAIL
    # ws_demo3：56 个关键数字里 2 个无源（4%）→ pass；被点名的是「粉丝 100 万以上」「预计到场 3 万人次」两处无测算依据的数字，
    # 带 [n] 的段、读图句、问卷结论句（问卷 / 受访者 口径）、账本数（300 万元 / 5000 万次 / 2.5 倍）都算有源
    d = demo3["checks"]["unsourced_numbers"]
    assert d["status"] == PASS and d["total"] == 56 and d["unsourced"] == 2 and d["value"] < 0.1
    assert [(h["file"], h["line"], h["term"]) for h in d["locations"]] == [("05_正文/03_品牌定位与传播策略.md", 20, "100 万"), ("05_正文/03_品牌定位与传播策略.md", 20, "3 万人次")]
    # 无正文 → skipped；有正文但没有带关键单位的数字 → pass 0
    assert _by_id(checks.run(tmp_path / "empty", track="brand", write=False))["unsourced_numbers"]["status"] == SKIP
    r = checks.run(_mini_ws(tmp_path / "nonum", "# 1 x\n\n本章无数字。\n", project={"track": "brand"}), track="brand", write=False)
    assert _by_id(r)["unsourced_numbers"]["status"] == PASS and _by_id(r)["unsourced_numbers"]["total"] == 0


def test_style_checks_in_chapter_mode(tmp_path: Path) -> None:
    """单章快检也跑三项，且账本（numbers.json）在单章模式下可用。"""
    _mini_ws(tmp_path, _NUMBERS_CHAPTER, project=_SYN_PROJECT, numbers=_LEDGER)
    r = checks.run_chapter(tmp_path, "brand", "05_正文/01_市场分析.md")
    assert {"cliche_terms", "unsourced_numbers", "synthetic_research_markers"} <= set(checks.CHAPTER_CHECK_IDS)
    c = _by_id(r)
    assert c["unsourced_numbers"]["total"] == 13 and c["unsourced_numbers"]["unsourced"] == 5  # 300 万元 / 4,860 万元 经账本放过
    assert c["synthetic_research_markers"]["status"] == FAIL and c["synthetic_research_markers"]["locations"][0]["line"] == 19  # 「样本量 412 份」无标记
    assert c["cliche_terms"]["status"] == PASS and c["cliche_terms"]["hits"] == 0
