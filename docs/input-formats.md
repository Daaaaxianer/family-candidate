# 输入文件格式

[返回 README](../README.md) · [完整参数](parameters.md) · [结果规则](results.md)

## 蛋白与 CDS FASTA

```fasta
>protein_a.1 description or other metadata
MSTAVK
>protein_a.2
MSTAVKGG
```

程序使用标题的第一个空白分隔字段作为 ID，例如 `protein_a.1`。描述部分不用于匹配，导出完整蛋白时保留原标题。

要求：

- 文件至少有一条非空序列。
- 同一个文件中不允许重复 ID；不同 ID 可以有相同序列。
- ID 精确匹配，`gene1` 不会匹配 `gene10`，`.` 等字符按字面处理。
- 可包含换行的序列；空行忽略。
- 输入蛋白应是可供 HMMER/BLAST 分析的氨基酸序列；输入 CDS 应是对应编码序列。

`--proteins` 是待鉴定全蛋白组，`--reference-proteins` 是已验证参考蛋白，不能互换。
`--cds` 只参与序列导出，不参与蛋白显著性或域规则判定。

## HMM 文件

`--hmm` 接受 HMMER3 profile 文件，开头为 `HMMER3/`，具有 `NAME` 和有效的 `LENG` 信息。
可以使用对应家族的 Pfam HMM 或自行构建且经过核查的模型。多个模型可以放在同一个文件内进行基本搜索；开启 `--refine` 时必须只有一个模型。

仓库附带 `NB-ARC.hmm` 的 accession 为 `PF00931.22`。它是示例模型，不意味着所有含该域的蛋白都属于同一个功能亚类。

## HMMER domtblout

`--domtbl` 接受 **hmmsearch 的域表输出**：

```bash
hmmsearch --cpu 4 --domtblout original.domtblout family.hmm proteome.fasta > original.log
```

导入的是 `original.domtblout`，不是 `original.log`。`hmmscan` 的 query/target 含义不同，其域表不能直接导入；`--tblout` 也不能代替 `--domtblout`。

| 域表列号（从 1 开始） | 内容 | 本工具用途 |
|---:|---|---|
| 1 / 3 | 目标蛋白 ID / 长度 | 核查是否对应输入蛋白 |
| 4 / 5 / 6 | 模型名 / accession / 长度 | 模型信息及覆盖度计算 |
| 7 | 全序列 E-value | 蛋白整体显著性 |
| 13 | 域 independent E-value | 当前域显著性 |
| 14 | 域 bit score | 最佳域排序 |
| 16 / 17 | hmm_from / hmm_to | 模型覆盖度 |
| 18 / 19 | ali_from / ali_to | 蛋白片段坐标，1-based、两端包含 |

忽略空行和 `#` 注释；记录不足 22 列、坐标越界或非有限 E-value 会报错。会核查蛋白 ID、蛋白长度、原始模型名和模型长度。
这些检查不能证明两个同名等长模型完全相同，导入者仍应确认结果来自指定模型和蛋白版本。

导入后只能筛选已经报告的命中。若原始搜索的报告阈值更严格，后续放宽筛选不会恢复被原始工具省略的结果，需要重新搜索。

## BLASTP TSV

**必须使用以下 12 列，不能直接使用默认 `outfmt 6` 输出：**

```text
qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen
```

每行一个 HSP，空白分隔，推荐 TAB；不带表头，允许 `#` 注释及空行。方向固定为参考蛋白 query、待鉴定蛋白 subject。

可按如下方式生成：

```bash
makeblastdb -in proteome.fasta -dbtype prot -out protein_db
blastp \
  -query references.fasta \
  -db protein_db \
  -evalue 1e-5 \
  -num_threads 4 \
  -max_target_seqs 100000 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" \
  -out blast.tsv
```

`100000` 是示例命中数量上限，可换成待鉴定蛋白组记录数；上限太小可能遗漏补充候选。程序自行运行 BLAST 时将它设置为输入蛋白记录数。

导入时核查 subject ID 与长度、query/subject 坐标及 E-value。筛选使用单个 HSP 的 E-value 和 query coverage，不合并多个 HSP，也不按 pident 做额外筛选。
自行搜索时的 `-evalue` 应能包含期望导入的命中。`--blast-tsv` 文件中的低质量行被过滤，但原始文件仍复制到运行目录供检查；仅 BLAST 未通过的蛋白不会因此新增到证据表。

人工格式示例见 [synthetic.blast.tsv](../examples/synthetic.blast.tsv)。

## InterProScan 标准 TSV

`--interpro-tsv` 接受标准的、无表头的 InterProScan TSV，而不是 Excel、JSON、GFF3 或自定义两列表。
程序不运行 InterProScan；可扫描全部输入蛋白，也可扫描候选完整蛋白序列，随后重新运行程序导入结果。

| 列号 | 内容 | 要求 |
|---:|---|---|
| 1 | protein accession | 与输入 FASTA ID 完全相同 |
| 2 | sequence MD5 | 等于输入蛋白序列大写后的 MD5 |
| 3 | sequence length | 等于输入蛋白长度 |
| 4 | analysis | 原数据库/分析方法名称 |
| 5 | signature accession | 必需/禁用域规则可匹配此编号 |
| 6 | signature description | 保留标准字段位置 |
| 7 / 8 | start / stop | 有效的 1-based 蛋白坐标 |
| 9 | score | 保留标准字段位置；本工具不再对它统一设阈值 |
| 10 | status | 仅 `T` 的记录加入域集合 |
| 11 | date | 保留标准字段位置 |
| 12 / 13 | InterPro accession / description | 可选；至少有 13 列且 accession 非 `-` 时加入域集合 |

文件至少需要 11 列。使用 `IPR...` 作为规则时，必须确保输出具有对应的 InterPro 映射字段；没有这些字段时不能匹配 IPR 编号。
含未知蛋白 ID、长度不符或 MD5 不符会报错。不得修改标题或截短序列后沿用之前的扫描。

对每个蛋白，程序将 signature 和 InterPro 编号收集为集合，然后检查全部必需域和任意禁用域。它不统计域拷贝数量，也不核查排列或保守位点。

**应保留每条已扫描蛋白的完整命中记录。** 若只保留目标域行，可能隐藏禁用域；若只保留部分 signature，也可能导致必需域缺失。未扫描、完全无命中或没有可用 `T` 记录，统一保留待复核，程序不能从 TSV 判断具体是哪一种情况。

人工格式示例见 [synthetic.interpro.tsv](../examples/synthetic.interpro.tsv)。

## 蛋白—基因座—CDS 映射

`--id-map` 接受带表头的 TAB 分隔文件，至少包含以下三个字段：

```text
protein_id	gene_id	cds_id
protein_a.1	gene_a	transcript_a.1
protein_a.2	gene_a	transcript_a.2
protein_b.1	gene_b	transcript_b.1
```

上述示例使用实际 TAB 分隔；逗号分隔 CSV 不兼容。可参考 [synthetic.id_map.tsv](../examples/synthetic.id_map.tsv) 的真实文件。

- 必须覆盖 `--proteins` 中**全部蛋白**，不只是预期候选；每个 protein_id 恰好一行，不得出现未知蛋白。
- gene_id 不能为空。同一基因座的多个转录本填相同 gene_id。
- cds_id 字段必须存在；没有对应 CDS 时允许字段留空。
- 没有提供 `--cds` 也能利用映射统计基因座和选择代表蛋白；此时不检查 CDS 是否存在。
- 提供 CDS 文件后，根据 cds_id 精确查找；缺失记录写清单，不改变蛋白候选状态。
- CDS FASTA 保留自身 ID，不为配合蛋白 ID 而改名。

不提供映射时，gene_id 不推断，CDS 默认按相同 protein_id 查找。多个蛋白映射同一 CDS 时，对应输出集合中的 CDS 只导出一次。

## 人工示例文件

仓库的 `examples/synthetic.*` 用于演示各类格式和筛选分支，可直接运行 README 的导入命令。
其中 HMMER、BLAST、InterProScan 行均为人工构造；不能作为真实搜索结果或家族成员注释使用。
