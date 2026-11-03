# 公开副本变更记录

整理日期：2026-11-03。原项目和原始文件保持不变。

## 处理规则

- 文本中的本机绝对路径将 `/Users/<本机用户名>/` 改为 `/Users/USER/`，保留后续目录和相对路径。
- 带明确Codex/ChatGPT会话URL前缀的会话标识替换为公开占位，不改科学运行ID、种子和数据标识。
- 本轮不重算、不重绘、不编译；科学数字、条件、计算逻辑和数组结果保持原有版本。
- 论文PDF仅改保存文件名，字节未改。原论文README保留为README.source-history.md，另新增当前使用说明。
- 不复制签名承诺书、个人工具详情、历史工作日志、编译产物及官方原始大附件。
- 公开运行说明补充：旧Q4 cache-only命令只适用于原始未脱敏包；公开副本给出省略两个cache参数的新物理解命令。本轮只核读代码，未执行该新命令。
- 补入科学支持Markdown实际引用的4张既有图，图像字节未改。

## 修改文件

历史SHA与运行证据绑定原始文件。下列公开副本已改变文本，不能将其标为仍满足原文件逐字节SHA绑定；新的SHA见PUBLICATION_MANIFEST.json。

- `reproduction/q4/base/data/full/TERRAIN.json`：local_username=2
- `reproduction/q4/base/results/inverse/INVERSE_RESULTS.json`：local_username=6
- `reproduction/q4/data/FORMAL_POLICY.json`：local_username=1
- `reproduction/q4/data/MORRIS_QUERIES.json`：local_username=76
- `reproduction/q4/data/cache/new1/CACHE_ENTRY.json`：local_username=2
- `reproduction/q4/data/cache/new1/ORIGINAL_CACHE_ENTRY.json`：local_username=1
- `reproduction/q4/data/cache/new1/evidence/RUN_RECEIPT.json`：local_username=3
- `reproduction/q4/data/cache/new2/CACHE_ENTRY.json`：local_username=2
- `reproduction/q4/data/cache/new2/ORIGINAL_CACHE_ENTRY.json`：local_username=1
- `reproduction/q4/data/cache/new2/evidence/RUN_RECEIPT.json`：local_username=3
- `reproduction/q4/data/cache/old0/CACHE_ENTRY.json`：local_username=4
- `reproduction/q4/data/cache/old0/ORIGINAL_CACHE_ENTRY.json`：local_username=2
- `reproduction/q4/data/cache/old0/evidence/RUN_RECEIPT.json`：local_username=3
- `reproduction/q4/data/cache/old1/CACHE_ENTRY.json`：local_username=4
- `reproduction/q4/data/cache/old1/ORIGINAL_CACHE_ENTRY.json`：local_username=2
- `reproduction/q4/data/cache/old1/evidence/RUN_RECEIPT.json`：local_username=3
- `reproduction/q4/data/cache/old2/CACHE_ENTRY.json`：local_username=4
- `reproduction/q4/data/cache/old2/ORIGINAL_CACHE_ENTRY.json`：local_username=2
- `reproduction/q4/data/cache/old2/evidence/RUN_RECEIPT.json`：local_username=3
- `reproduction/q4/data/cache/old3/CACHE_ENTRY.json`：local_username=4
- `reproduction/q4/data/cache/old3/ORIGINAL_CACHE_ENTRY.json`：local_username=2
- `reproduction/q4/data/cache/old3/evidence/RUN_RECEIPT.json`：local_username=3
- `reproduction/q4/data/local/TERRAIN.json`：local_username=2

- `reproduction/软件与运行命令.txt`：补充公开版Q4运行方式、旧缓存绑定失效说明与本轮未实跑边界。

除上列脱敏和运行说明补充文件外，复制文件与清单中原始来源字节相同。新生成文件仅包括本说明、目录README与PUBLICATION_MANIFEST.json。清单自排除以避免循环哈希。
