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
