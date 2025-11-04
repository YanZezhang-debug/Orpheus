# 同源搜索功能说明

## 概述

Orpheus 现在支持基于同源证据的 ORF 预测，而不是仅仅依赖"最长ORF"的启发式假设。通过与已知蛋白质数据库（如 SwissProt）进行比对，可以：

- **提高预测准确性**：保留有同源证据的 ORF，即使它们不是最长的
- **减少假阳性**：过滤掉长但无功能的假 ORF
- **支持可变剪接**：识别较短但功能重要的转录本亚型
- **基于进化保守性**：优先选择在已知蛋白质中有同源的 ORF

## 为什么需要同源证据？

"最长 ORF" 假设的局限性：

1. **可变剪接**：同一基因的较短剪接体可能是主要功能形式
2. **不完整转录本**：真正的全长转录本可能因5'端缺失而显得"不够长"
3. **假ORF**：UTR区域可能包含偶然出现的长ORF，但无生物学功能
4. **进化保守性被忽略**：最长的ORF可能缺乏同源支持

同源搜索通过与已知蛋白质比对，提供了更可靠的证据。

## 支持的工具

### 1. Diamond（推荐）

- **速度**：比 BLASTP 快 100-10000 倍
- **准确性**：与 BLASTP 相当
- **适用场景**：大规模转录组分析

**安装**：
```bash
# Conda
conda install -c bioconda diamond

# 或从源码安装
wget http://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
sudo mv diamond /usr/local/bin/
```

### 2. NCBI BLASTP

- **标准工具**：经典的序列比对工具
- **高灵敏度**：适合检测远缘同源
- **适用场景**：小规模数据或需要最高灵敏度时

**安装**：
```bash
# Conda
conda install -c bioconda blast

# 或下载预编译版本
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

## 配置方法

### 步骤 1：准备蛋白质数据库

#### 使用 Diamond

如果您已有 `.dmnd` 格式的数据库（如 `swissprot_db.dmnd`），可以直接使用。

如果只有 FASTA 格式，需要先构建数据库：

```bash
# 下载 SwissProt 数据库
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# 构建 Diamond 数据库
diamond makedb --in uniprot_sprot.fasta --db swissprot_db
# 这会生成 swissprot_db.dmnd
```

#### 使用 BLASTP

```bash
# 如果使用 FASTA 格式，构建 BLAST 数据库
makeblastdb -in uniprot_sprot.fasta -dbtype prot -out uniprot_sprot_blast
```

### 步骤 2：修改配置文件

编辑您的配置文件（如 `config/default.yaml` 或自定义配置）：

```yaml
transdecoder:
  # ... 其他配置 ...
  
  # 启用同源搜索
  homology_search:
    enabled: true                    # 设为 true 启用
    
    # 是否只使用完整 ORF (type:complete) 进行同源搜索
    # true: 仅完整 ORF（推荐用于 PCR/克隆，详见下文说明）
    # false: 包括不完整 ORF（适合基因发现/功能注释）
    use_complete_only: true
    
    tool: "diamond"                  # 使用 Diamond（推荐）或 "blastp"
    diamond_executable: "diamond"    # Diamond 命令路径
    blastp_executable: "blastp"      # BLASTP 命令路径（如果使用 blastp）
    
    # 数据库路径（根据您的实际路径修改）
    database: "/mnt/data/01.databases/swissprot_db.dmnd"
    
    evalue: 1e-5                     # E-value 阈值
    max_target_seqs: 1               # 每个查询保留的最大目标序列数
    threads: 8                       # 使用的线程数
    extra_params: ""                 # 其他自定义参数
```

### 步骤 3：运行流程

```bash
python orpheus_cli.py --config your_config.yaml
```

## 工作原理

### 完整流程

1. **TransDecoder.LongOrfs**：识别所有候选 ORF（生成 `longest_orfs.pep`）

2. **同源搜索**（新增）：
   - 使用 Diamond/BLASTP 将候选 ORF 与蛋白质数据库比对
   - 生成 `blastp_results.outfmt6`（表格格式）

3. **TransDecoder.Predict**：
   - 使用 `--retain_blastp_hits` 参数读取同源搜索结果
   - **优先保留有同源证据的 ORF**
   - 同时保留超过长度阈值的长 ORF（作为补充）

### 输出文件

同源搜索会在 TransDecoder 工作目录中生成：

```
<input>.transdecoder_dir/
  ├── longest_orfs.pep              # 候选 ORF 蛋白序列
  ├── blastp_results.outfmt6        # 同源搜索结果（表格格式）
  └── ... 其他 TransDecoder 文件 ...
```

最终预测结果：

```
<input>.transdecoder.pep            # 预测的蛋白质序列（基于同源证据）
<input>.transdecoder.cds            # 预测的 CDS 序列
<input>.transdecoder.gff3           # GFF3 注释文件
<input>.transdecoder.bed            # BED 格式注释
```

## 参数调优

### 完整 ORF 过滤（use_complete_only）

**非常重要**：此选项控制是否只使用完整的 ORF 进行同源搜索。

#### 什么是完整 ORF？

完整 ORF (type:complete) 是指同时满足以下条件的开放阅读框：
- 包含起始密码子（ATG）
- 包含终止密码子（TAA/TAG/TGA）
- 从起始到终止密码子之间无其他终止密码子

#### 配置选项

```yaml
homology_search:
  use_complete_only: true   # 默认值，推荐
```

| 设置 | 说明 | 适用场景 | 优点 | 缺点 |
|-----|------|---------|------|------|
| **true** (默认) | 仅使用完整 ORF | • PCR 引物设计<br>• 基因克隆<br>• 蛋白表达 | • 可扩增完整基因<br>• 提高匹配质量<br>• 减少假阳性 | • 可能遗漏不完整转录本 |
| **false** | 使用所有 ORF | • 基因发现<br>• 功能注释研究<br>• 不完整转录组 | • 覆盖更全面<br>• 适合不完整数据 | • 可能增加假阳性<br>• 包含无法克隆的序列 |

#### 使用建议

**推荐使用 `use_complete_only: true`（默认）**，因为：

1. **PCR 实验需求**：不完整的 ORF 缺少起始/终止密码子，无法设计有效的扩增引物
2. **提高准确性**：完整 ORF 更可能代表真实的蛋白编码序列
3. **蛋白表达**：只有完整 ORF 才能在表达系统中产生完整蛋白
4. **避免混淆**：过滤掉片段序列，减少误导性结果

**何时使用 `use_complete_only: false`**：

- 转录组组装不完整（缺失 5' 端）
- 纯粹用于功能注释，不进行实验验证
- 研究目的是基因发现而非基因克隆

#### 运行时日志

当 `use_complete_only: true` 时，日志会显示：
```
筛选完整的 ORF 用于同源搜索（推荐用于后续 PCR 实验）...
✓ 将使用完整 ORF 文件进行同源搜索: xxx.complete.pep
```

当 `use_complete_only: false` 时，日志会显示：
```
使用所有 ORF（包括不完整的）进行同源搜索
```

### E-value 阈值

- **`1e-5`**（默认）：较严格，高置信度匹配
- **`1e-3`**：较宽松，可能包含更多远缘同源
- **`1e-10`**：非常严格，只保留高度相似的匹配

```yaml
homology_search:
  evalue: 1e-5
```

### 最大目标序列数

- **`1`**（默认）：每个 ORF 只保留最佳匹配
- **`5`**：保留前 5 个匹配（可能增加计算时间）

```yaml
homology_search:
  max_target_seqs: 1
```

### 线程数

根据您的服务器配置调整：

```yaml
homology_search:
  threads: 16  # 使用 16 个线程
```

### 数据库选择

| 数据库 | 特点 | 适用场景 |
|--------|------|----------|
| **SwissProt** | 高质量、人工审核 | 推荐，准确性高 |
| **TrEMBL** | 自动注释、规模大 | 需要更广泛覆盖时 |
| **RefSeq** | NCBI 参考序列 | 特定物种研究 |
| **自定义数据库** | 特定物种蛋白组 | 近缘物种研究 |

## 性能对比

### Diamond vs BLASTP

使用 10,000 个候选 ORF 与 SwissProt (560K 序列) 比对：

| 工具 | 时间 | 内存 | 敏感度 |
|------|------|------|--------|
| Diamond | ~2 分钟 | ~2 GB | 高 |
| BLASTP | ~2 小时 | ~1 GB | 非常高 |

**建议**：对于大多数转录组项目，Diamond 是最佳选择。

## 示例配置

### 配置 1：PCR/克隆项目（推荐）

适用于需要进行 PCR 扩增或基因克隆的项目。

```yaml
transdecoder:
  homology_search:
    enabled: true
    use_complete_only: true  # 只使用完整 ORF（推荐）
    tool: "diamond"
    database: "/mnt/data/01.databases/swissprot_db.dmnd"
    evalue: 1e-5
    max_target_seqs: 1
    threads: 8
```

### 配置 2：功能注释研究

适用于纯粹的功能注释，不需要实验验证。

```yaml
transdecoder:
  homology_search:
    enabled: true
    use_complete_only: false  # 包括不完整 ORF，增加覆盖率
    tool: "diamond"
    database: "/mnt/data/01.databases/swissprot_db.dmnd"
    evalue: 1e-5
    max_target_seqs: 1
    threads: 8
```

### 配置 3：高灵敏度（BLASTP + 宽松阈值）

适用于检测远缘同源或者小数据集。

```yaml
transdecoder:
  homology_search:
    enabled: true
    use_complete_only: true
    tool: "blastp"
    database: "/data/uniprot_sprot.fasta"
    evalue: 1e-3
    max_target_seqs: 5
    threads: 16
```

### 配置 4：仅使用长度标准（不推荐）

不使用同源证据，仅基于长度阈值（传统方法）。

```yaml
transdecoder:
  homology_search:
    enabled: false  # 禁用同源搜索
  
  retain_long_orfs: true
  retain_long_orfs_length: 900
```

## 常见问题

### Q1: 同源搜索失败会怎样？

如果同源搜索失败，Orpheus 会：
1. 记录警告日志
2. 回退到仅使用长度标准（`retain_long_orfs`）
3. 继续完成流程

### Q2: 可以同时使用同源证据和长度标准吗？

可以！TransDecoder 会：
- 保留所有有同源证据的 ORF
- **同时**保留超过长度阈值的长 ORF
- 最终输出是两者的并集

### Q3: 数据库路径在 Windows 下如何配置？

```yaml
# Windows 路径
database: "D:\\databases\\swissprot_db.dmnd"

# 或使用正斜杠（推荐）
database: "D:/databases/swissprot_db.dmnd"
```

### Q4: 如何验证同源搜索是否生效？

查看日志输出：

```
步骤 2/3: 同源搜索 (Diamond/BLASTP)
运行 Diamond BLASTP 同源搜索...
  查询: .../longest_orfs.pep
  数据库: /mnt/data/01.databases/swissprot_db.dmnd
  E-value: 1e-05
✓ Diamond 搜索完成
  找到 8523 个同源匹配

步骤 3/3: 预测编码区域 (TransDecoder.Predict)
使用同源搜索证据: .../blastp_results.outfmt6
```

### Q5: 为什么推荐使用 `use_complete_only: true`？

完整 ORF 对于实验验证至关重要：

1. **PCR 扩增**：需要起始密码子（ATG）和终止密码子设计引物
2. **蛋白表达**：不完整的 ORF 无法在表达系统中产生功能蛋白
3. **避免浪费**：筛选出的基因都可以直接用于实验，不会在实验阶段才发现序列不完整
4. **提高质量**：完整 ORF 更可能是真实的蛋白编码基因

**如果日志显示"筛选完整的 ORF"**，说明此功能已启用。

### Q6: 我的转录组不完整，是否应该设置 `use_complete_only: false`？

**建议**：
- 如果**主要目的是 PCR/克隆**：保持 `true`，宁可少一些结果，但确保都是可用的
- 如果**主要目的是功能注释**：可以设置为 `false`，获得更全面的注释结果
- 如果**不确定**：保持默认值 `true`

即使转录组不完整，完整 ORF 仍然是最有价值的部分。

### Q7: 同源搜索需要多少时间？

取决于：
- 候选 ORF 数量（通常 10K-100K）
- 数据库大小（SwissProt ~560K 序列）
- 工具选择（Diamond vs BLASTP）
- 线程数

**典型时间**（50K 候选 ORF + SwissProt）：
- Diamond (8 threads): 5-10 分钟
- BLASTP (8 threads): 2-5 小时

## 参考文献

1. TransDecoder: https://github.com/TransDecoder/TransDecoder
2. Diamond: Buchfink B, Xie C, Huson DH. (2015) Nature Methods 12, 59–60
3. BLAST+: Camacho C, et al. (2009) BMC Bioinformatics 10:421

## 相关文档

- [TransDecoder 文档](TRANSDECODER.md)
- [流程说明](PIPELINE_FLOW.md)
- [配置文件说明](../config/default.yaml)

