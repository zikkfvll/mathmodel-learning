# 复现指南

本目录提供四问的程序、必要输入和已保存结果。可以只运行需要的一问，或只重绘已有结果，无需重新下载和清理全部原始资料。

## 环境

建议使用 Python 3.12，并为项目建立虚拟环境。各问依赖文件保存了原计算所用版本；按所运行的问题安装即可。以下从本比赛目录开始：

```bash
python3.12 -m venv .venv
source .venv/bin/activate
```

Windows下使用对应的`.venv\Scripts\activate`。Q4的C后端需要C99编译器；中文绘图需要可用中文字体。数据准备中的栅格和矢量工具仅在重新清理原始附件时需要。

下面Q1、Q2、Q3的代码块各自从本比赛目录开始。使用同一终端时，先回到本比赛目录，再执行下一问的`cd`；不要在上一问的子目录中直接粘贴。Q4的两个代码块均从本比赛的`reproduction/`目录运行。

## Q1：环境特征

```bash
cd reproduction/q1
python -m pip install -r requirements.txt
python -B sol.py --mode enhance --d1-csv data/d1_era5/era5_precipitation_hourly_2021_2026.csv --d0-cache results/d0_window_cache.json.gz --baseline-results results/baseline --output-dir new_run
python -B draw.py --results-dir results --output-dir new_figures
```

`new_run`保存重新计算结果；上述绘图命令读取仓库中的`results`。若要绘制新计算结果，应将`--results-dir`指向相应目录。ERA5为格点回溯产品，不能据运行成功把它改称站点实测。

## Q2：局部材料柱

```bash
cd reproduction/q2
python -m pip install -r requirements.txt
python -B sol.py --output-dir new_results
python -B draw.py --results-dir results --output-dir new_figures
```

配置见[config.json](../reproduction/q2/config.json)。四初态响应、临界边界和积分诊断分别保存在`results/initial_state_*.csv`等表中。`data/v2_reused_results/`是明确复用的前期计算输入，保留它才能按当前配置运行。

## Q3：条件首达与概率

```bash
cd reproduction/q3
python -m pip install -r requirements.txt
python -B sol.py --data data --results results_new --workers 2
python -B draw.py --results results --output figures_new
```

该问包括概率积分、代理训练比较、固定模拟成员验证及敏感性分析。随机种子用于复现相同实验，重复同种子不会新增一次独立统计认证。若要看某项已有结果，可直接阅读`results/`，不必先完成整个运行。

## Q4：二维传播与物理核验

以下从`reproduction/`运行：

```bash
python -m pip install -r q4/requirements-lock.txt
python -B q4/code/sol.py decide --out q4_run/new_physics --inputs q4/data/local --backend c --volume 475682.8460010884 --mu .075 --gp-state q4/data/GP_STATE.json
```

`--out`目录必须此前不存在。此命令调用新的30 m物理解，并在新输出中保存结果和运行回执；耗时取决于机器。当前发布整理只核对了这条路径的代码逻辑，没有再次执行该命令，因此复现结论应以实际新输出为准。物理解失败、证据不一致或出域时，接口保留`unknown`，不会以GP估计直接替代核验。

论文附录和历史材料里另有`--cache-entry ... --cache-only`命令。它属于原始归档环境：公开文件移除了本机用户名，改变了`TERRAIN.json`和部分历史回执的字节，旧哈希不再适用于这些文本。请使用上面的新计算命令，省略两个缓存参数。`q4/data/cache/`用于查看原实验的输入、证据和结构，不应当作当前版本可直接通过的缓存认证。

重绘已有Q4结果：

```bash
python -B q4/code/draw.py --root q4 --out q4_guard_figures
python -B q4/followup/code/draw.py --root q4/followup --out q4_followup_figures
python -B q4/base/code/draw.py --root q4/base --out q4_base_figures
```

## 原始数据准备

`reproduction/prepare/`提供清理、下载与转换程序，对应依赖见其中的`requirements.txt`。完整清理需要自行提供官方附件1—8；站点子入口只生成站点及日温结果。取数程序在缺少原响应时会联网获取公开数据。现有四问计算已带必要输入，通常可以跳过这一阶段。

## 运行结果与论文

使用新目录保存输出，便于与仓库结果比较。重新计算后，先核对单位、类别、有效域和关键结果表，再决定是否更新图与论文。当前论文和图件来自已保存的研究结果；本次文档整理未重算模型或重新编译。关于数据身份和使用条件，参见[数据来源与许可](../reproduction/数据来源与许可.txt)。
