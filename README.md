# family-candidate

基于 HMMER 的蛋白编码基因家族候选鉴定工具，支持同源序列搜索补充、结构域结果核查和成员证据导出。

```text
原始家族 HMM ─────────────────┐
可选物种特异 HMM ─────────────┼→ 合并搜索证据 → 核查结构域 → accepted / review / rejected
可选参考蛋白 BLASTP ──────────┘
```

原始模型的候选始终保留；物种特异 HMM 是补充步骤。没有提供结构域验证文件及明确的必需域规则时，通过搜索阈值的候选进入 `review`。

## 选择运行方式

| 你的数据和目标 | 主要参数 | 外部工具 |
|---|---|---|
| 有 HMM 和全蛋白组，先找候选 | `--hmm --proteins --outdir` | `hmmsearch` |
| 想增加物种特异模型搜索 | 在基本搜索上加 `--refine` | `hmmsearch`；种子足够时还需要 `hmmbuild`、MUSCLE v5 |
| 有参考蛋白，想补查同源成员 | 加 `--reference-proteins` | `makeblastdb`、`blastp` |
| 已有 HMMER/BLAST 结果，只做筛选 | `--domtbl`，可加 `--blast-tsv` | 单纯导入不需要搜索工具 |
| 想输出通过结构域规则的成员 | 加 `--interpro-tsv` 和至少一个 `--required-domain` | 导入已有 InterProScan TSV，不在程序内运行 InterProScan |
| 想统计基因座、对应不同的 CDS ID | 加 `--id-map`；要导出 CDS 再加 `--cds` | 无新增依赖 |

各模块可组合。`--reference-proteins` 与 `--blast-tsv` 二选一；`--domtbl` 与 `--refine` 可以一起使用，但第二轮搜索仍需要 HMMER。

## 环境与输入

- Python 3.9+，Python 程序仅使用标准库。
- 真实 HMM 搜索需要 HMMER 3；物种模型重建使用 MUSCLE **v5**。
- 执行 BLAST 补充需要 NCBI BLAST+。
- Windows 上可将示例中的 `python3` 改为 `python`。多行示例使用 Bash 的 `\`；在 PowerShell 中将命令合并为一行。
- 所有相对数据路径都相对于**运行命令时所在目录**；从其他目录运行时，给程序和数据使用对应的路径。含空格的路径需加引号。

必需输入：HMMER3 格式的家族 HMM、待鉴定物种的蛋白 FASTA，以及新的或空的输出目录。CDS 是可选输入，缺少 CDS 不影响蛋白候选鉴定。

FASTA ID 使用标题中第一个空白符之前的内容，要求唯一。每个转录本可以作为一条蛋白保留；要按基因座统计，须提供映射表。

详细说明：[完整参数](docs/parameters.md) · [输入文件格式](docs/input-formats.md) · [结果与筛选规则](docs/results.md)

## 快速开始：无需外部工具的演示

在仓库根目录执行：

```bash
python3 family_candidate.py \
  --hmm NB-ARC.hmm \
  --proteins examples/synthetic.protein.fasta \
  --cds examples/synthetic.cds.fasta \
  --domtbl examples/synthetic.domtblout \
  --blast-tsv examples/synthetic.blast.tsv \
  --interpro-tsv examples/synthetic.interpro.tsv \
  --required-domain PF00931 \
  --id-map examples/synthetic.id_map.tsv \
  --outdir demo_result
```

`examples/synthetic.*` 是人工构造的程序示例，序列及命中不代表真实生物学结果。预期输出为 2 个 `accepted`、0 个 `review`、1 个 `rejected`、2 个候选基因座。

先查看 `demo_result/candidate_evidence.tsv` 和 `demo_result/run.json`。重复演示时换一个输出目录，例如 `demo_result_2`。

## 场景一：只进行原始 HMM 搜索

安装 HMMER 后，使用仓库附带的拟南芥示例：

```bash
python3 family_candidate.py \
  --hmm NB-ARC.hmm \
  --proteins Athaliana_167_TAIR10.protein.fa \
  --cds Athaliana_167_TAIR10.cds.fa \
  --sequence-evalue 1e-10 \
  --domain-evalue 1e-10 \
  --min-hmm-coverage 0.6 \
  --cpu 4 \
  --outdir results/nb_arc_candidates
```

程序执行原始 HMM 搜索，筛选域命中并导出候选，不重建物种模型，也不执行 BLAST。没有结构域验证时，`accepted` 文件为空是正常行为，候选在 `review` 文件中。

这里的 `0.6` 是示例覆盖阈值，程序默认值为 `0`。覆盖不足但满足显著性阈值的命中保留待复核。NB-ARC 命中支持“含该域的候选”，具体家族或 NLR 亚类还需符合相应定义。

## 场景二：同源搜索补充并核查结构域

准备已验证的参考蛋白 `references.fasta`，以及针对待鉴定蛋白序列生成的完整 InterProScan 标准 TSV `scan.tsv`：

```bash
python3 family_candidate.py \
  --hmm family.hmm \
  --proteins proteome.fasta \
  --reference-proteins references.fasta \
  --blast-evalue 1e-5 \
  --blast-min-query-coverage 0.5 \
  --interpro-tsv scan.tsv \
  --required-domain PF00931 \
  --outdir results/validated_candidates
```

替换 `family.hmm`、蛋白文件和 `PF00931` 为研究家族对应的输入与规则。参考蛋白是 BLAST 的 query，待鉴定全蛋白组是数据库。BLAST 独有候选也进入结构域核查。

程序导入 `scan.tsv`，不自动运行 InterProScan。可扫描完整待鉴定蛋白组，或先搜索候选再扫描候选的完整蛋白序列；必须保留每条已扫描蛋白的全部 signature 命中。

例如，一个定义要求同时具有 `PFXXXXX` 和 `PFYYYYY`，写成：

```text
--required-domain PFXXXXX --required-domain PFYYYYY
```

这是 **AND** 条件；两个域均须存在。`PFXXXXX`、`PFYYYYY` 仅为格式占位符，须换成真实编号。不能写成逗号分隔的单个参数。当前规则不支持“两个域任选其一”。

## 场景三：增加物种特异 HMM

```bash
python3 family_candidate.py \
  --hmm family.hmm \
  --proteins proteome.fasta \
  --refine \
  --seed-evalue 1e-20 \
  --seed-domain-evalue 1e-20 \
  --seed-min-coverage 0.6 \
  --min-seeds 3 \
  --sequence-evalue 1e-10 \
  --domain-evalue 1e-10 \
  --outdir results/refined_candidates
```

种子阈值只控制哪些原始模型命中用于重建 HMM；最终阈值控制原始和重建模型的候选筛选。两组参数相互独立。

每个蛋白选择最佳合格域，MUSCLE 比对后交给 `hmmbuild`。种子少于 `--min-seeds` 时跳过重建并记录警告，继续输出原始模型结果；足够时，重建或第二轮搜索失败会终止运行。

`--refine` 仅接受单模型 HMM。种子数量达到要求也不保证覆盖全部亚群，应检查种子与比对。此场景仍未提供验证文件，符合搜索阈值的结果依然需要复核。

## 场景四：导入已有搜索结果

```bash
python3 family_candidate.py \
  --hmm family.hmm \
  --proteins proteome.fasta \
  --domtbl original.domtblout \
  --blast-tsv blast.tsv \
  --interpro-tsv scan.tsv \
  --required-domain PF00931 \
  --outdir results/imported_candidates
```

`--domtbl` 替代第一轮 `hmmsearch`；`--blast-tsv` 替代 BLAST 搜索。HMM 和蛋白 FASTA 仍然必需，用于核查结果及提取序列。

HMMER 文件必须为 `hmmsearch --domtblout` 输出。BLAST 文件须采用本工具约定的 12 列，默认 BLAST `outfmt 6` 不兼容。导出方法见[输入文件格式](docs/input-formats.md)。

## 输出与常见问题

| 输出 | 用途 |
|---|---|
| `candidate_evidence.tsv` | 每个已报告 HMM 命中蛋白及通过 BLAST 筛选的蛋白，其证据、状态和原因 |
| `candidates.ids.txt / .protein.fasta` | `accepted` 与 `review` 的合并候选集 |
| `accepted.*` | 满足搜索证据与指定结构域规则的成员 |
| `review.*` | 缺少验证或存在覆盖不足等待核查情况 |
| `rejected.*` | 已进入证据表但未通过规则的蛋白 |
| `domain_evidence.tsv` | 全部已报告 HMM 域命中 |
| `run.json` | 参数、输入哈希、工具信息、计数、警告及运行状态 |

添加 `--cds` 后会生成相应 CDS FASTA；添加 `--id-map` 后可统计基因座并导出代表蛋白。[完整结果说明与错误排查](docs/results.md)。

- **没有 accepted？** 检查是否同时提供了 `--interpro-tsv` 和 `--required-domain`，再查看 `reason`。
- **只有蛋白数，没有基因数？** 加入覆盖全部输入蛋白的 `--id-map`；程序不从 ID 后缀猜测基因座。
- **CDS 数量较少？** 查看 `missing_cds.tsv`；多个蛋白映射同一 CDS 时该 CDS 在对应 FASTA 中仅导出一次。
- **想放宽已导入结果的阈值？** 可以调整筛选参数，但原始搜索未报告的命中须重新搜索。
- **输出目录已存在？** 使用新目录；程序不覆盖非空运行目录。

## 测试与旧入口

```bash
python3 -m unittest discover -s tests -v
python3 family_candidate.py --help
```

测试包括解析、筛选、异常处理和模拟外部工具流程；不替代真实工具集成测试或生物学准确率评估。

六位置参数的旧 Shell 命令仍可使用；两个 Shell 文件名均调用同一 Python 核心。参数含义、输出对应关系和 Perl 辅助入口见[兼容命令](docs/parameters.md#兼容命令)。

## 功能范围

支持域存在/缺失规则，尚未实现域数量与排列、催化位点、GA 自动筛选、基因组漏注释补查或伪基因判定。`accepted` 表示通过配置的计算规则，不等同于实验功能验证。

运行目录和 `run.json` 含输入路径及环境信息，属于分析记录；提交源码时无需包含运行结果。默认 `.gitignore` 排除常见结果目录及运行记录。
