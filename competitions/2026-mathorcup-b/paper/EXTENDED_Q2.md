# 第二问证据分配与压缩映射

本文以局部饱和基质的起始边界为主线。原始详细记录与结果表均保留于 `Q2_initial_state_20261029_v3/`，以下定位以该目录为根；此文件不是新的计算或独立验收报告。

|主文内容|直接来源|分配原则|
|---|---|---|
|研究对象、阶段几何、速度与外地材料|`Q2_paper.md` §2.1–2.3；`data/BC_geometry_tuples.csv`、`data/cghc_f7b_digitized.csv`、`data/blatten_gsd_summary.csv`|图2-1展示运动和粒度；保留不同来源身份、相对时钟、坡度和渗透率转移。|
|状态、本构、无拉力分层及预算|`Q2_paper.md` §3.1–3.3，`config.json`|保留模型方程、两支连接、初态准入、自由边界与固水预算。|
|四初态和近似诊断|`results/initial_state_boundary.csv`、`initial_state_response.csv`、`initial_state_cumulative_check.csv`；图07/08|主文表2-1及图2-2、2-3，区分根与上邻侧事件的偏差。|
|材料与前缀窗|`results/material_window_sensitivity.csv`；`Q2_paper.md` §4.6|保留临界结构变化与三个窗的累计组合，完整成员和各窗输入仍见来源表。|
|压力形状与外部密土反例|`results/closure_sensitivity.csv`、`iverson2000_directional_comparison.csv`；图05/04|图2-4与正文保留配对闭合和方向反例。|
|条件纯融冰库存|`results/pure_melt_counterfactual.csv`；`Q2_paper.md` §4.4|正文压缩为假设、范围和水源意义；不推断真实局部水量。|
|数值与预算检查|`REVIEW.md`、`verification/`、`results/implementation_checks.csv` 和 `closure_checks.csv`|检查证据保留于原包；正文不把数值自洽写成自然模型验证。|
|外部湿床、鼓试验、2017试验元数据、静床对照|`Q2_paper.md` §2.3–2.5、4.3、5；`data/`及`results/external_process_comparison.json`；静床见 `Q2_delivery_20261029_r2/results/reference_compatibility.csv`|这些材料分别涉及既成来流卷入、流态定义、观测资格或另一操作定义下的对照，未作为本模型预测成绩。|

八幅原图 PDF 与 PNG 均原样复制在 `figures/q2/`：`fig01_new_evidence` → 正文图2-1；`fig07_initial_state_boundary` → 图2-2；`fig08_initial_state_diagnostics` → 图2-3；`fig05_closure_sensitivity` → 图2-4。其余 `fig02_response_control`（参考响应/自由排水）、`fig03_critical_conditions`（厚度与结构）、`fig04_external_process`（外部孔压方向）、`fig06_cumulative_shear`（累计剪切和前缀窗）作为扩展图保留，完整图注见原稿。所有图中的参数与结果沿用原包，未重新绘图。主文用的关键数值和范围可追至上述表；该映射不等同于重新运行或实际自然过程验证。

首轮返修时，正文四幅图的TeX显示窗口从原PDF上边裁去20 pt、下边裁去24 pt；`figures/q2/display/` 中PNG按原像素同范围裁边，用于带图Markdown。原始PDF/PNG保留，裁除图外总标题和底注，坐标、图例、子图标题和数据保留；条件解释由正文图注承载。此显示裁边不是重算或重绘。
