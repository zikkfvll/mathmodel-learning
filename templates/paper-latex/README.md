# 数学建模论文 LaTeX 模板

一个可直接开始写作的中文数学建模论文模板，包含完整章节结构、按问组织的建模段落、图表与公式示例，以及可选的声明和附录。

**[下载独立模板包](https://github.com/zikkfvll/mathmodel-learning/raw/refs/heads/main/downloads/paper-latex-template.zip)** · **[查看七页样张](模板样张.pdf)** · **[详细使用说明](使用说明.md)**

## 快速开始

1. 下载并解压模板包，或复制本目录到新项目。
2. 修改`config.tex`中的题名、队号、题号、关键词和开关。
3. 在`sections/`填入各章内容，按需增减问题数量。
4. 以`main.tex`为入口，使用XeLaTeX编译两遍。在Overleaf中选择XeLaTeX。

```bash
xelatex -interaction=nonstopmode -halt-on-error main.tex
xelatex -interaction=nonstopmode -halt-on-error main.tex
```

## 包含什么

| 文件 | 用途 |
|---|---|
| `main.tex` | 章节顺序与文档入口 |
| `config.tex` | 题名、身份字段、关键词和可选章节开关 |
| `personal-style.tex` | 字体、页边距、行距、标题及图表编号 |
| `sections/` | 摘要、问题重述与分析、假设、符号、各问、评价、文献与附录 |
| `模板样张.pdf` | 完整版式示例 |

模板保留A4、25 mm页边距、宋黑体层级、三线表、连续公式编号及按问编号的图表。默认四问，可在`sections/questions.tex`调整；题名过长时自动换行。

模板来自[本次比赛论文](https://github.com/zikkfvll/mathmodel-learning/tree/main/competitions/2026-mathorcup-b/paper)的结构和样式，不包含原研究数据、结论、图件或个人字体。它是个人模板，正式使用前按当届官方要求调整匿名、封面、目录、页数与声明。第三方类文件的来源见[NOTICE](NOTICE.md)。
