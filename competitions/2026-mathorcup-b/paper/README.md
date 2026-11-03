# 论文与图表

**[阅读论文 PDF](paper.pdf)** · **[LaTeX 源码](main.tex)** · **[图件](figures/)** · **[程序与复现](../docs/reproduction.md)**

论文题目为《吉隆灾害链的环境特征、触发条件、预警与仿真》。前30页包括摘要、目录、正文、参考文献及文后声明，后283页是完整程序附录。阅读研究方法与结果时可先看前30页。

## 文件说明

- `main.tex`与`cumcmthesis.cls`：论文入口及类文件。
- `figures/`：主文和补充分析的图件；提供已有的PDF、PNG或SVG格式。
- `tools/appendix_source.tex`与`appendix_code/`：19份程序的源码附录输入。
- `Q1_SUPPORTING_ANALYSES.md`至`Q4_SUPPORTING_ANALYSES.md`及Q2补充文件：正文之外的方法比较、敏感性与诊断。
- `DATA_LICENSES/`：地形资料的许可与归属。

## 编译

在本目录运行，使用XeLaTeX至少编译两遍：

```bash
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

也可使用`tectonic main.tex`。需要相应TeX宏包及中文字体；工程在macOS使用Songti/Heiti，在其他系统回退Fandol。新编译结果名为`main.pdf`，系统字体或宏包版本变化可能影响分页，编译后应核对图表、目录和附录。仓库中的`paper.pdf`保留研究完成时的PDF。

想复用空白结构而不带入本次研究内容，请使用[个人论文模板](../../../templates/paper-latex/README.md)。该模板可独立下载。

## 补充说明

[写作与排版](STYLE_AND_LAYOUT.md)说明章节和格式，[四问衔接](Q4_INTEGRATION_NOTES.md)解释局部模型、概率与传播的关系，[版本记录](VERSION_RECORD.md)给出当前论文文件指纹。计算应使用[复现指南](../docs/reproduction.md)中的命令；源码附录保留的Q4历史缓存命令仅适用于原始归档环境。
