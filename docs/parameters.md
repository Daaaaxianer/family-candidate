# 参数参考

[返回 README](../README.md) · [输入格式](input-formats.md) · [结果规则](results.md)

本页适用于 `family_candidate.py`。参数使用 `--参数名 值`；开关 `--refine` 不带值。查看帮助：

```bash
python3 family_candidate.py --help
```

## 输入与运行设置

| 参数 | 必需/默认 | 作用与约束 |
|---|---|---|
| `--hmm FILE` | 必需 | 原始 HMMER3 家族模型；导入已有命中时仍需提供。普通搜索允许多个模型，`--refine` 仅允许一个模型。 |
| `--proteins FILE` | 必需 | 待鉴定的完整蛋白 FASTA，ID 唯一；不是参考蛋白文件。 |
| `--outdir DIR` | 必需 | 新目录或空目录。相对路径相对于当前工作目录。 |
| `--cds FILE` | 不提供 | 可选 CDS FASTA，仅用于序列导出；不参与候选显著性筛选。 |
| `--id-map FILE` | 不提供 | 带表头的 TAB 分隔映射表，覆盖所有输入蛋白；启用基因座统计与代表蛋白导出。可不提供 CDS 文件。 |
| `--cpu N` | `1` | HMMER `--cpu` 和 BLASTP `-num_threads` 的值，正整数；不传给 MUSCLE。纯导入时不执行搜索。 |
| `--domtbl FILE` | 不提供 | 使用已有原始模型 `hmmsearch --domtblout` 输出，跳过第一轮搜索。不是普通终端日志，也不是 `--tblout`。 |

相对文件路径不会自动相对于程序文件所在目录解释。例如：

```bash
python3 /path/to/family-candidate/family_candidate.py \
  --hmm /path/to/data/family.hmm \
  --proteins /path/to/data/proteome.fasta \
  --outdir /path/to/results/run_01
```

以上路径是占位示例，请换成实际位置。

## 最终 HMM 命中筛选

| 参数 | 默认 | 判定 |
|---|---:|---|
| `--sequence-evalue X` | `1e-10` | 整条蛋白的 E-value ≤ X。 |
| `--domain-evalue X` | `1e-10` | 当前域的 independent E-value ≤ X；不是 conditional E-value。 |
| `--min-hmm-coverage F` | `0.0` | 当前域覆盖 HMM 的比例 ≥ F。范围 0–1，不是 0–100。 |

对**同一个域命中**同时检查：

```text
sequence_evalue <= sequence-evalue
AND domain_ievalue <= domain-evalue
AND (hmm_to - hmm_from + 1) / model_length >= min-hmm-coverage
```

蛋白有至少一个合格域时，具有合格 HMM 搜索证据。不能将一个域的 E-value 与另一个域的覆盖度拼接判定。
原始模型和物种模型都使用这组最终筛选参数，任一来源合格即可提供 HMM 证据。

`0.6` 表示覆盖模型的 60%，不是覆盖蛋白全长的 60%。默认 `0` 表示不设置最终覆盖下限。
显著性合格但覆盖不足的蛋白，在没有合格 BLAST 补充证据时保留为 `review`，即使结构域规则通过也不自动升级。

两个 E-value 参数都要求有限正数；数值越小越严格。搜索结果中的 E-value `0` 可正常解析并通过正阈值。
这些参数是**结果筛选条件**，不同于 HMMER 输出报告阈值。自行导入结果时，须保证原搜索的报告范围能够涵盖所需命中。

## 可选物种 HMM 重建

| 参数 | 默认 | 作用 |
|---|---:|---|
| `--refine` | 关闭 | 从原始 HMM 命中中提取种子域，MUSCLE 比对、hmmbuild 建模，再搜索同一蛋白组；保留原始候选。 |
| `--seed-evalue X` | `1e-20` | 用于重建的种子蛋白整体 E-value 上限。 |
| `--seed-domain-evalue X` | `1e-20` | 用于重建的种子域 independent E-value 上限。 |
| `--seed-min-coverage F` | `0.6` | 种子域覆盖 HMM 的最低比例，范围 0–1。 |
| `--min-seeds N` | `3` | 最小种子数，正整数；单模型重建中，每个蛋白最多贡献一个最佳合格域。 |

未加 `--refine` 时，这四个种子筛选参数不用于候选筛选。仅修改 `--seed-evalue` **不会**同步修改 `--seed-domain-evalue`；二者独立。

最佳种子域依次按 independent E-value 最小、域得分最高、覆盖度最高、坐标靠前选择。种子少于 N 时跳过重建并在 `run.json` 写警告。
启用 `--refine` 后会预先要求 `hmmsearch` 可用；种子足够时才检查 MUSCLE v5 和 `hmmbuild`。启用重建但命令失败时整次运行失败，不静默退回部分结果。

例如，严格种子筛选、相对宽松的最终筛选：

```text
--refine --seed-evalue 1e-20 --seed-domain-evalue 1e-20
--seed-min-coverage 0.6 --min-seeds 3
--sequence-evalue 1e-10 --domain-evalue 1e-10
```

## BLASTP 补充

| 参数 | 默认 | 作用与约束 |
|---|---:|---|
| `--reference-proteins FILE` | 不提供 | 以参考蛋白为 query，新建待鉴定蛋白数据库并执行 BLASTP；需要 BLAST+。 |
| `--blast-tsv FILE` | 不提供 | 导入本工具约定的 12 列 BLASTP 输出，跳过 BLAST 搜索。 |
| `--blast-evalue X` | `1e-5` | 单个 HSP 的 E-value 上限，有限正数。 |
| `--blast-min-query-coverage F` | `0.5` | 单个 HSP 的参考蛋白覆盖比例下限，范围 0–1。 |

`--reference-proteins` 与 `--blast-tsv` 不能一起使用。不提供任何一个参数时，不执行 BLAST；其筛选参数没有实际筛选对象。

```text
query_coverage = (qend - qstart + 1) / qlen
```

例如参考蛋白为 400 aa，单个 HSP 覆盖 200 aa，则比例为 0.5。不是待鉴定蛋白的 subject coverage。
同一目标蛋白任一 HSP 同时通过 E-value 和覆盖规则，即为 BLAST 合格候选。多个 HSP 的覆盖不累加；不额外要求 pident 达到固定阈值。
展示的最佳 BLAST 证据按 E-value 最小、bitscore 最高选择。

## 结构域核查

| 参数 | 默认 | 作用与约束 |
|---|---|---|
| `--interpro-tsv FILE` | 不提供 | 导入完整的 InterProScan 标准 TSV 命中；程序本身不运行 InterProScan。 |
| `--required-domain ID` | 空列表 | 必需 signature accession 或 InterPro accession；可重复指定，所有条目须存在。 |
| `--forbidden-domain ID` | 空列表 | 禁用 signature accession 或 InterPro accession；可重复指定，任一个出现即违反规则。 |

完整验证需要 **`--interpro-tsv` + 至少一个 `--required-domain`**。只有验证文件而没有必需域，或只配置禁用域，均不能形成完整验证规则；候选保留待复核。

```text
--required-domain PF00931
--required-domain PF00931 --required-domain IPR000001
--required-domain PF00931 --forbidden-domain PFXXXXX
```

第二行要求两个编号均被记录。第三行中的 `PFXXXXX` 是占位符，换成实际需要排除的域。
每次参数只接受一个编号；逗号分隔不会被拆成多个域。不得将同一编号同时设为必需和禁用。Pfam 编号的版本后缀会归一化，例如 `PF00931.22` 与 `PF00931` 等价。

规则比较的是 InterProScan 第 5 列 signature accession，以及存在时的第 12 列 InterPro accession。
`IPR...` 规则需要输出含有相应 InterPro 映射字段。序列 ID、长度、MD5 均须对应本次输入。
这套规则不检查域数量、排列、位置先后、保守位点或实验功能。

## 外部工具路径

| 参数 | 默认命令 | 使用条件 |
|---|---|---|
| `--hmmsearch-exe PATH` | `hmmsearch` | 未导入第一轮结果，或开启重建。 |
| `--hmmbuild-exe PATH` | `hmmbuild` | 开启重建且种子足够。 |
| `--muscle-exe PATH` | `muscle` | 开启重建且种子足够；必须是 v5。 |
| `--blastp-exe PATH` | `blastp` | 提供 `--reference-proteins`。 |
| `--makeblastdb-exe PATH` | `makeblastdb` | 提供 `--reference-proteins`。 |

可使用 PATH 中的命令名，也可使用可执行文件的完整路径：

```text
--hmmsearch-exe "/opt/bio tools/hmmer/bin/hmmsearch"
```

参数只接受可执行文件，不要把命令参数附在这个字符串中。

## 参数组合速查

| 组合 | 实际行为 |
|---|---|
| 只有三个必需参数 | 执行原始 HMM 搜索；不运行 BLAST/重建；没有结构域核查，合格候选为 review。 |
| `--domtbl`，未加 `--refine` | 只导入原始 HMM 结果，不运行 hmmsearch。 |
| `--domtbl --refine` | 从导入结果选择种子，足够时执行 MUSCLE/hmmbuild/第二轮 hmmsearch。 |
| `--reference-proteins` | 执行 BLAST，不自动生成结构域验证文件。 |
| `--blast-tsv` | 只导入 BLAST；若未提供 `--domtbl`，第一轮 HMM 搜索仍会执行。 |
| `--interpro-tsv`，没有必需域 | 导入文件，但候选因规则未定义而保持 review。 |
| `--required-domain`，没有验证文件 | 候选因未进行域核查而保持 review。 |
| `--forbidden-domain`，没有必需域 | 未启用完整域规则，不能单独用于确认或按禁用域排除。 |
| `--id-map`，没有 `--cds` | 可以统计基因座及输出代表蛋白，不导出 CDS FASTA。 |
| `--cds`，没有 `--id-map` | 用蛋白 ID 精确查找相同 CDS ID；不统计基因座。 |

## 兼容命令

推荐新分析使用 `family_candidate.py`；旧 Shell 入口保留六位置参数：

```bash
bash batch_pfam_hmmer_to_family_candidate_search.sh \
  NB-ARC.hmm \
  Athaliana_167_TAIR10.protein.fa \
  Athaliana_167_TAIR10.cds.fa \
  1e-20 \
  1e-10 \
  NB-ARC.in.At
```

| 位置 | 含义 | 对应新参数 |
|---:|---|---|
| 1 | 原始 HMM | `--hmm` |
| 2 | 待鉴定蛋白组 | `--proteins` |
| 3 | CDS 文件，旧入口仍必需 | `--cds` |
| 4 | 严格种子阈值 | 同时设置 `--seed-evalue` 与 `--seed-domain-evalue` |
| 5 | 最终显著性阈值 | 同时设置 `--sequence-evalue` 与 `--domain-evalue` |
| 6 | 输出前缀 | 运行目录为 `<前缀>.run` |

旧入口固定启用 `--refine`，种子覆盖下限和最小种子数使用默认值。它不接受新 CLI 的可选模块；需要 BLAST、验证或自定义覆盖度时使用新入口。
旧 `.id.txt`、`.protein.fasta`、`.cds.fasta` 从运行目录的 `candidates.*` 复制，代表候选集。已有旧输出或非空运行目录会被拒绝。

`batch_pfam_hmmer_to_family_candidate_search.perl.sh` 调用同一 Shell 核心；同样使用 Python 3 和 MUSCLE v5。`PYTHON` 环境变量可指定 Python 解释器，默认 `python3`。

辅助命令：

```bash
python3 retrieve.seq.from.all.fasta.py ids.txt proteome.fasta selected.fasta
python3 refine_domain_seq.py original.domtblout proteome.fasta seed_domains.fasta 1e-20 --domain-evalue 1e-20 --min-coverage 0.6
```

序列提取采用完整 ID 匹配，目标 ID 缺失时在写文件前失败；ID 列表去空行、去重复，并保持列表顺序。
域提取命令的第四参数是全序列阈值；未指定 `--domain-evalue` 时，域阈值与它相同。域序列标题为 `protein|model|start-end`。
相应 `.pl` 文件转调 Python 实现。这两个辅助 Python 命令会覆盖指定的输出文件，新分析应使用独立文件名。
