# 结果、判定规则与排错

[返回 README](../README.md) · [完整参数](parameters.md) · [输入格式](input-formats.md)

## 先确认运行成功

命令成功返回退出码 `0`。检查 `run.json` 中 `status` 是否为 `success`，再读取最终表格。

- 输入/工具预检查失败：可能尚未创建运行目录，不保证有 `run.json`。
- 运行目录创建后发生错误：写入 `status=failed` 和 `error`，返回非零；中间文件不能作为完整结果使用。
- 搜索成功但没有命中：成功输出空 ID/FASTA 文件及带表头的空表，是合法的零结果。
- 非空输出目录：直接拒绝，须使用新的或空目录；不自动清理旧结果。

## 谁进入证据表

`candidate_evidence.tsv` 包含：

```text
原始/物种 HMMER 已报告命中的所有蛋白
UNION
至少一个 BLAST HSP 通过补充筛选的蛋白
```

HMMER 已报告但未通过筛选的蛋白也可能有一行 `rejected`。只有不合格 BLAST 命中、且没有 HMM 命中的蛋白，不会新增为一行。
因此表格及 `rejected` 都不是完整蛋白组的成员/非成员分类清单。

## accepted、review 和 rejected

先评价搜索证据，再对尚未因显著性不合格被排除的蛋白进行域核查。

### 搜索证据

| 情况 | 搜索阶段结果 |
|---|---|
| 任一 HMM 域通过全序列 E-value、域 independent E-value 和覆盖度三条件 | 合格 HMM 证据，继续域核查 |
| 没有合格 HMM 域，但有合格 BLAST HSP | 合格补充证据，继续域核查 |
| HMM 显著性合格，但覆盖不足，且没有合格 BLAST | 保留 review；继续检查域规则，但域通过也不自动升级 |
| 没有任何 HMM 显著性合格域，也没有合格 BLAST | rejected，原因 `search_threshold_not_met` |

HMM 的全序列与域阈值须在同一域记录上同时通过。BLAST 补充是独立候选入口，不要求该蛋白必须已有 HMM 命中。

### 域核查

完整规则需同时有 `--interpro-tsv` 和至少一个 `--required-domain`：

| 情况 | 最终状态 |
|---|---|
| 没有验证文件 | review |
| 有文件，但没有必需域规则 | review；禁用域规则也未启用 |
| 没有该蛋白的可用验证记录 | review |
| 有记录且含任一禁用域 | rejected |
| 有记录但缺少任一必需域 | rejected |
| 所有必需域存在，无禁用域，且搜索证据完整 | accepted |
| 域规则通过，但 HMM 覆盖不足且无合格 BLAST | review |

先判断禁用域，再判断必需域。覆盖不足的 review 仍可能因违反域规则改为 rejected。
`accepted` 只代表通过所配置的计算规则；功能、催化活性或具体亚家族归属需要其他证据。

### 四个例子

1. HMM 两种 E-value 与覆盖均合格，但没有 InterProScan 文件：`review`。
2. 仅 BLAST 合格，验证记录包含全部必需域且没有禁用域：`accepted`。
3. HMM 显著性合格，覆盖度 0.4，而下限为 0.6；无 BLAST，即使验证通过：`review`。
4. HMM 证据合格，验证中有其他域但缺少要求的家族域：`rejected`。完全没有验证记录则为 `review`，不是相同情况。

## 文件对应关系

| 文件 | 何时生成 | 内容 |
|---|---|---|
| `candidate_evidence.tsv` | 成功运行 | 每个进入证据集合的蛋白一行 |
| `domain_evidence.tsv` | 成功运行 | 全部已报告 HMM 域命中，含未通过阈值的域 |
| `candidates.ids.txt` / `candidates.protein.fasta` | 成功运行 | accepted + review |
| `accepted.ids.txt` / `accepted.protein.fasta` | 成功运行 | accepted；没有时为空文件 |
| `review.ids.txt` / `review.protein.fasta` | 成功运行 | review |
| `rejected.ids.txt` / `rejected.protein.fasta` | 成功运行 | 证据集合中的 rejected |
| 上述各集合的 `.cds.fasta` | 提供 `--cds` | 对应集合中可查到的 CDS，每个 CDS ID 仅一次 |
| `missing_cds.tsv` | 提供 `--cds` | candidates 中未找到 CDS 的蛋白及所用 CDS ID；无缺失时只有表头 |
| `gene_representatives.protein.fasta` | 提供 `--id-map` | candidates 中每个基因座一个代表蛋白 |
| `original.domtblout` | 第一轮完成 | 执行或导入的原始 HMM 域表 |
| `seed_domains.fasta` | 开启 `--refine` | 通过种子筛选的域片段；不足仍保留该文件 |
| `seed_alignment.fasta`、`species.hmm`、`species.domtblout` | 重建实际执行且成功 | 种子比对、物种模型、第二轮结果 |
| `blast.tsv` | 执行或导入 BLAST | 原始 HSP 行，含未通过后续筛选的行 |
| `interpro.tsv` | 导入验证文件 | 验证输入副本 |
| `*.log` | 对应外部命令执行 | 工具 stdout/stderr |
| `run.json` | 运行目录创建后 | 参数、输入 SHA256、命令、工具信息、状态、警告及成功时的计数 |

## 候选证据表字段

| 字段 | 含义 |
|---|---|
| `protein_id` | 输入蛋白 ID |
| `gene_id` | 映射表中的基因座，未提供映射则为空 |
| `cds_id` | 提供 CDS 时的查找 ID；没有 CDS 输入则为空 |
| `source` | 所有已报告 HMM 来源，加上合格 BLAST 来源；`pfam_hmm`、`species_hmm`、`blast` 用 `+` 连接 |
| `best_domain_source` | 表内展示的最佳域来自哪一轮 HMM 搜索 |
| `model_id` | 最佳域的 Pfam accession（去版本号），无 accession 则使用模型名 |
| `sequence_evalue` / `domain_ievalue` | 最佳域所属蛋白整体 E-value / 该域 independent E-value |
| `domain_score` / `hmm_coverage` | 最佳域得分 / 模型覆盖比例 |
| `domain_start` / `domain_end` | 最佳域的蛋白坐标，1-based，两端包含 |
| `blast_query` / `blast_evalue` / `blast_query_coverage` | 最佳合格 HSP 的参考蛋白、E-value、query 覆盖比例 |
| `domains` | 验证文件中该蛋白的 signature 与 InterPro 编号集合，`;` 分隔 |
| `status` / `reason` | 最终状态 / 判定原因，多个原因用 `;` 分隔 |
| `cds_found` | 提供 CDS 时为 true/false，否则为空 |
| `gene_representative` | 提供映射时是否被选为基因座代表，否则为空 |

`source` 表示命中来源，不表示每个来源均通过筛选。应结合 `status`、指标和 `domain_evidence.tsv` 阅读。
仅 BLAST 候选的 HMM 指标为空是正常情况。

展示最佳域时优先在通过全部 HMM 条件的域中选择；若没有，则在显著性合格域中选择；仍没有则从全部域中选择。每层按 independent E-value 最小、得分最高、覆盖度最高排序。
`domain_evidence.tsv` 才是完整域列表；不能只根据单个展示域推断蛋白的完整域架构。

## reason 对照

| 原因 | 解释与处理 |
|---|---|
| `search_threshold_not_met` | 没有合格 HMM 显著性或 BLAST 证据；检查模型、参考序列与阈值 |
| `hmm_coverage_below_threshold` | HMM 显著性合格但域覆盖不足；检查截短注释、片段与家族边界 |
| `domain_validation_not_run` | 没有提供验证文件 |
| `required_domain_rule_not_set` | 提供文件但未指定必需域；补上明确规则 |
| `no_validation_record` | 无可用验证记录；检查是否实际扫描、是否有命中及记录状态 |
| `forbidden_domain:...` | 出现指定禁用域；核查家族定义与真实域组合 |
| `required_domain_missing:...` | 缺少指定必需域；检查编号及是否保留完整扫描记录 |
| `search_and_domain_rules_passed` | 搜索证据与域存在/缺失规则通过 |

## 计数与代表蛋白

`run.json` 的 `counts.accepted/review/rejected` 是证据表中各状态的**蛋白数**。
`candidate_proteins = accepted + review`。
只有提供 `--id-map` 才报告 `candidate_gene_loci`，为 candidates 中不同 gene_id 的数量。

每个候选基因座选择一个代表：accepted 优先于 review，其次最长蛋白，仍同长则取 protein ID 按字典序最前者。
完整候选集不因代表选择而删除转录本；不同基因座即使序列相同也不合并。
多个转录本可映射到同一个 CDS，因此 CDS 数量不一定等于候选蛋白数。

## 常见错误

| 错误或现象 | 检查方向 |
|---|---|
| `executable not found` | 检查 PATH，或传入相应 `--...-exe` 完整路径；MUSCLE 要求 v5 |
| `output directory must be new or empty` | 换输出目录，不重复写入已有运行 |
| `invalid domtblout` / `original HMM model` 不匹配 | 是否误用了 tblout、hmmscan、其他模型或其他蛋白版本 |
| `invalid BLAST TSV` | 是否严格使用规定的自定义 12 列和正确 query/subject 方向 |
| `protein MD5 mismatch` | InterProScan 文件是否来自相同 ID 和相同完整序列；仅长度相同不够 |
| `ID map must cover the input proteome` | 映射须覆盖全部输入蛋白，而非只填写候选 |
| `duplicate FASTA ID` | 先修复数据中的重复 ID，程序不自行覆盖或猜测 |
| `too few qualifying seed proteins` 警告 | 原始搜索仍成功；可检查种子与阈值后决定是否重建 |
| `alignment IDs, lengths or ungapped sequences` 不匹配 | 比对文件的 ID、等长性或去 gap 后序列与种子不一致 |
| 没有 accepted | 首先检查验证文件与必需域规则，再读 review/rejected 的原因 |

## 分享与版本控制

源码、文档、测试和人工示例可纳入版本控制。运行结果由各分析生成，`run.json` 和工具日志含数据路径与环境信息；公开分享结果前应检查其内容。
默认忽略 `results/`、`*.run/`、`demo_result*/`、`verification/`、`run.json` 和 `*.log`；如果使用其他自定义结果目录，也应将该目录加入 `.gitignore`。
