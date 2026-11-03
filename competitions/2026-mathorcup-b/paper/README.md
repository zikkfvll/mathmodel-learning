# 论文工程

[paper.pdf](paper.pdf)是现有最终PDF，实际313页：前30页为摘要、目录、正文、参考文献和文后声明，后283页为完整源码附录。[main.tex](main.tex)是相应源稿，含17个实际引用的图文件，`tools/appendix_source.tex`和`appendix_code/`保留19份程序的附录输入。

本公开副本没有重新编译论文。原工程曾使用XeTeX兼容引擎与本机中文字体；换环境需安装宏包及可用字体，并复核版面。[STYLE_AND_LAYOUT.md](STYLE_AND_LAYOUT.md)保存既有风格说明。实际编译main.tex通常生成main.pdf，与这里保存的paper.pdf文件名不同；不要将新生成PDF冒称与原最终PDF逐字节相同。

[README.source-history.md](README.source-history.md)是原项目的历史写作说明，其中“待审”、旧页数和旧main.pdf文件名不代表当前公开版状态。其他科学补充说明保留原内容和引用；公开整理与指纹变化见上级PUBLICATION_CHANGES.md。

论文PDF及源码附录保留的Q4 `--cache-only` 是原始未脱敏包的历史命令。公开副本脱敏后旧request/receipt哈希不再匹配，不应照抄该命令回放；公开副本的新物理解命令与未实跑说明见上级README及reproduction/软件与运行命令.txt。
