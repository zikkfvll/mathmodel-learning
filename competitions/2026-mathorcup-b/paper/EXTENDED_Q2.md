# 第二问阅读与结果索引

问题二以局部饱和基质的临界边界为主线，依次解释资料、状态方程、初态效应、近似诊断及材料转移。完整论述见[论文](paper.pdf)，细节见[Q2补充分析](Q2_SUPPLEMENTARY_ANALYSES.md)。

| 内容 | 对应材料 | 阅读重点 |
|---|---|---|
| 阶段几何、速度与外地材料 | [data](../reproduction/q2/data/)中的`BC_geometry_tuples.csv`、`cghc_f7b_digitized.csv`、`blatten_gsd_summary.csv` | 区分事件读点、外地材料和设计条件 |
| 状态、本构、分层与预算 | [config.json](../reproduction/q2/config.json)、[sol.py](../reproduction/q2/sol.py) | 状态准入、自由边界、固水预算 |
| 四初态与近似诊断 | [results](../reproduction/q2/results/)中的`initial_state_boundary.csv`、`initial_state_response.csv`、`initial_state_cumulative_check.csv` | 临界根与上邻侧首达偏差属于不同诊断 |
| 材料与时间窗 | [material_window_sensitivity.csv](../reproduction/q2/results/material_window_sensitivity.csv) | 阈值组合依赖初态、材料和规定窗口 |
| 压力剖面与外部反例 | [closure_sensitivity.csv](../reproduction/q2/results/closure_sensitivity.csv)、[iverson2000_directional_comparison.csv](../reproduction/q2/results/iverson2000_directional_comparison.csv) | 配对改变闭合，保留密土方向不符的结果 |
| 条件融冰库存 | [pure_melt_counterfactual.csv](../reproduction/q2/results/pure_melt_counterfactual.csv) | 条件水源比较，不推断真实局部水量 |
| 数值与预算 | [implementation_checks.csv](../reproduction/q2/results/implementation_checks.csv)、[closure_checks.csv](../reproduction/q2/results/closure_checks.csv) | 数值自洽与自然验证分别解释 |

## 八组图

[figures/q2](figures/q2/)保留八组图及其已有格式：

- `fig01_new_evidence`：运动与粒度资料。
- `fig02_response_control`：参考响应与自由排水对照。
- `fig03_critical_conditions`：临界厚度与结构。
- `fig04_external_process`：外部孔压方向对照。
- `fig05_closure_sensitivity`：剖面闭合与条件敏感性。
- `fig06_cumulative_shear`：累计剪切与前缀窗口。
- `fig07_initial_state_boundary`：初态与条件边界。
- `fig08_initial_state_diagnostics`：初态积分与近似诊断。

论文选用其中四组作为主文图。`display/`中的PNG是显示裁边版本，坐标、图例和数据保持原图内容；重新绘制时应从程序和结果表生成，而不是根据截图倒取数据。
