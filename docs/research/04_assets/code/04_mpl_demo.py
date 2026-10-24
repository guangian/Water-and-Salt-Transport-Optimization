# -*- coding: utf-8 -*-
"""生成 11 张商业样图到 out/charts/，并打印所用字体。"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import bizcharts as bc  # noqa: E402

OUT = os.path.join(os.path.dirname(__file__), "out", "charts")
os.makedirs(OUT, exist_ok=True)

font = bc.setup_chinese_font()
bc.apply_style()

years = [2021, 2022, 2023, 2024, 2025, 2026]
bc.save_fig(bc.market_size_stacked_bar(
    years, {"硬件": [120, 150, 190, 240, 300, 370], "软件订阅": [30, 55, 90, 140, 210, 300], "服务": [10, 18, 30, 50, 80, 120]},
    title="智慧养老市场规模（亿元）", source="艾瑞咨询、团队测算"), f"{OUT}/01_市场规模堆叠柱.png")

bc.save_fig(bc.share_barh(["A公司", "B公司", "本项目", "C公司", "其他"], [28.5, 21.0, 6.5, 15.0, 29.0],
                          title="细分市场份额（2025）", highlight="本项目", source="公司年报"), f"{OUT}/02_份额条形.png")

bc.save_fig(bc.growth_line(years, {"营业收入": [320, 480, 690, 980, 1350, 1800], "毛利": [96, 158, 248, 372, 540, 756]},
                           title="收入与毛利增长（万元）", unit="万元", yoy_bars=[0, 0.5, 0.44, 0.42, 0.38, 0.33],
                           source="财务模型"), f"{OUT}/03_增长折线.png")

bc.save_fig(bc.waterfall(["营业收入", "原材料", "人工", "营销费用", "研发费用", "管理费用", "税金"],
                         [1800, -620, -310, -220, -180, -120, -90], title="2026E 利润瀑布（万元）", source="财务模型"),
            f"{OUT}/04_瀑布图.png")

bc.save_fig(bc.funnel(["曝光", "访问", "注册", "试用", "付费"], [1_200_000, 360_000, 90_000, 27_000, 8_100],
                      title="用户转化漏斗", source="后台数据"), f"{OUT}/05_漏斗.png")

bc.save_fig(bc.tornado(["客单价", "获客成本", "续费率", "人工成本", "折现率"],
                       low=[1800, 3100, 2200, 2900, 3000], high=[3600, 2300, 3200, 2500, 2400], base=2700,
                       source="敏感性分析"), f"{OUT}/06_龙卷风.png")

bc.save_fig(bc.breakeven(fixed_cost=600, price=0.45, var_cost=0.18, max_units=4000, title="盈亏平衡分析",
                         unit="万元", qty_unit="台"), f"{OUT}/07_盈亏平衡.png")

bc.save_fig(bc.cashflow_combo(["2024", "2025", "2026", "2027", "2028"], [300, 700, 1400, 2300, 3300],
                              [900, 950, 1200, 1700, 2200], title="现金流预测（万元）"), f"{OUT}/08_现金流组合.png")

bc.save_fig(bc.competition_matrix(["A公司", "B公司", "C公司", "本项目", "D公司"], [8, 6, 3, 4, 7], [7, 5, 3, 8, 4],
                                  size=[90, 60, 30, 20, 50], highlight="本项目"), f"{OUT}/09_竞争矩阵.png")

bc.save_fig(bc.gantt([("需求调研", 1, 2, "筹备期"), ("产品研发", 2, 4, "研发期"), ("试点运营", 5, 3, "试运营"),
                      ("市场推广", 7, 4, "推广期"), ("融资与扩张", 9, 4, "推广期")], total_periods=12),
            f"{OUT}/10_甘特图.png")

bc.save_fig(bc.business_model_canvas({
    "重要伙伴": ["社区养老服务中心", "三甲医院老年科", "硬件代工厂"], "关键业务": ["健康监测设备研发", "AI 预警算法"],
    "核心资源": ["专利 6 项", "医疗数据库", "研发团队"], "价值主张": ["24h 无感健康监测", "异常 3 分钟预警", "家属远程可视"],
    "客户关系": ["订阅制服务", "社区驻点"], "渠道通路": ["民政采购", "线上商城", "养老机构"],
    "客户细分": ["独居老人家庭", "养老机构", "社区政府"], "成本结构": ["硬件 BOM 成本", "研发人员薪酬", "获客成本"],
    "收入来源": ["设备销售", "年度订阅费", "数据增值服务"]}), f"{OUT}/11_商业模式画布.png")

print("charts written:", sorted(os.listdir(OUT)))
print("font:", font)
