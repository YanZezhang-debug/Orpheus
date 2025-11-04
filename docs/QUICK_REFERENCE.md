# Orpheus 快速参考

## 命令行参数速查表

| 参数 | 长格式 | 说明 | 示例 |
|------|--------|------|------|
| `-h` | `--help` | 显示帮助信息 | `orpheus -h` |
| `-v` | `--version` | 显示版本信息 | `orpheus -v` |
| `-i FILE` | `--input FILE` | **必需** 输入 FASTA 文件 | `-i trinity.fasta` |
| `-c FILE` | `--config FILE` | 配置文件路径 | `-c my_config.yaml` |
| `-o FILE` | `--output FILE` | 输出文件路径 | `-o result.fasta` |
| `-t N` | `--threads N` | 线程数（覆盖配置文件） | `-t 32` |
| | `--start-from STEP` | 从指定步骤开始 | `--start-from cdhit` |

## 快速使用示例

### 1️⃣ 查看帮助

```bash
orpheus -h
```

### 2️⃣ 基本使用

```bash
# 使用默认配置从头运行
orpheus -i trinity_assembly.fasta

# 使用自定义配置
orpheus -i trinity_assembly.fasta -c my_config.yaml
```

### 3️⃣ 指定线程数

```bash
# 使用 32 个线程
orpheus -i trinity_assembly.fasta -t 32

# 使用 16 个线程，仅运行 TransDecoder
orpheus -i input.fasta --start-from transdecoder -t 16
```

### 4️⃣ 步骤控制（支持自动文件查找）

```bash
# 从头开始（默认）
orpheus -i trinity_assembly.fasta

# 从 CD-HIT 开始
orpheus -i trinity_assembly.fasta --start-from cdhit

# 跳过 CD-HIT，从 TransDecoder 开始
# 方式1：手动指定输入文件
orpheus -i cdhit_result.fasta --start-from transdecoder

# 方式2：自动查找工作目录中的 cdhit_result.fasta（推荐）
orpheus --start-from transdecoder
# ✓ 自动找到 ./orpheus_output/cdhit_result.fasta

# 直接运行评分步骤（自动查找所有需要的文件）
orpheus --start-from busco_scoring
# ✓ 自动找到 GFF3 和 CD-HIT 结果
```

**自动文件查找**：
- 跳过步骤时，会自动在工作目录中查找中间文件
- 无需手动指定输入文件（除非使用自定义文件）
- 详见：[步骤跳过与文件自动查找](PIPELINE_FLOW.md#步骤跳过与文件自动查找)

### 5️⃣ 完整示例

```bash
# 指定所有参数
orpheus \
  -i trinity_assembly.fasta \
  -c my_config.yaml \
  -t 24 \
  -o final_result.fasta \
  --start-from cdhit
```

## 步骤选项说明

| 步骤 | 说明 | 输入要求 |
|------|------|----------|
| `beginning` | 从头开始执行所有步骤（默认） | Trinity 组装结果 |
| `busco_before` | BUSCO 质量评估（处理前） | FASTA 文件 |
| `cdhit` | 从 CD-HIT 去冗余开始 | Trinity 组装结果 |
| `transdecoder` | TransDecoder ORF 预测 | FASTA 文件（可以是 CD-HIT 结果） |
| `busco_after` | BUSCO 质量评估（处理后） | FASTA 文件 |
| `busco_scoring` | **新增** BUSCO ORF 评分和转录本筛选 | GFF3 + CD-HIT 结果（自动查找） |

## 配置文件快速设置

### 最小配置（仅必需项）

```yaml
cdhit:
  identity: 0.95

transdecoder:
  min_protein_length: 100
```

### 启用同源搜索（推荐用于 PCR/克隆）

```yaml
transdecoder:
  homology_search:
    enabled: true
    use_complete_only: true  # 仅使用完整 ORF（推荐）
    tool: "diamond"
    database: "/path/to/swissprot_db.dmnd"
    threads: 16
    evalue: 1e-5
```

**关键参数**：
- `use_complete_only: true` - 只使用带起始/终止密码子的完整 ORF
  - ✅ 推荐用于 PCR 引物设计和基因克隆
  - ✅ 确保所有结果都可以直接用于实验
- `use_complete_only: false` - 包括不完整 ORF
  - 适合纯功能注释研究
  - 适合不完整的转录组数据

### 配置评分系统（新增）

```yaml
scoring:
  # 评分权重（总和应为 1.0）
  completeness_weight: 0.4  # ORF 完整性
  length_weight: 0.3        # ORF 长度
  homology_weight: 0.3      # 同源性
  
  # 最小分数阈值
  min_score_threshold: 0.5
  
  # 每个基因保留的最大转录本数（0 = 保留所有）
  max_transcripts_per_gene: 1
  
  # ORF 长度归一化参数（bp）
  max_orf_length: 3000
```

### 设置线程数

```yaml
cdhit:
  threads: 16

transdecoder:
  homology_search:
    threads: 32
```

**注意**：命令行 `-t` 参数会覆盖配置文件中的线程设置。

## 输出文件说明

### CD-HIT 输出

- `cdhit_output.fasta` - 去冗余后的序列
- `cdhit_output.fasta.clstr` - 聚类信息

### TransDecoder 输出

- `*.transdecoder.pep` - 预测的蛋白质序列
- `*.transdecoder.cds` - 预测的编码序列
- `*.transdecoder.gff3` - GFF3 格式的 ORF 注释
- `*.transdecoder.bed` - BED 格式的 ORF 位置

### 如果启用同源搜索

- `blastp.outfmt6` 或 `diamond.outfmt6` - 比对结果

### BUSCO ORF 评分输出（新增）

- `scored_transcripts.tsv` - 转录本评分表（TSV 格式）
  - 列：transcript_id, gene_id, orf_type, orf_length, homology, score, quality
  - 按基因分组，每组保留评分最高的转录本

## 常见问题速查

### Q1: 如何查看所有参数？

```bash
python orpheus_cli.py -h
```

### Q2: 如何指定线程数？

使用 `-t` 参数：
```bash
python orpheus_cli.py -i input.fasta -t 32
```

### Q3: 如何跳过某些步骤？

使用 `--start-from` 参数：
```bash
# 跳过 CD-HIT，仅运行 TransDecoder
python orpheus_cli.py -i input.fasta --start-from transdecoder
```

### Q4: 配置文件在哪里？

- 默认：`config/default.yaml`
- 示例：`example_config.yaml`
- 自定义：使用 `-c` 参数指定

### Q5: 如何启用同源搜索？

在配置文件中设置：
```yaml
transdecoder:
  homology_search:
    enabled: true
    use_complete_only: true  # 推荐：仅完整 ORF
    tool: "diamond"
    database: "/path/to/database.dmnd"
```

**重要**：设置 `use_complete_only: true` 确保只使用完整 ORF（有起始和终止密码子），这对于后续 PCR 和克隆实验至关重要。

### Q6: Diamond 和 BLASTP 如何选择？

- **推荐使用 Diamond**：速度快（100-10000x），适合大数据集
- **使用 BLASTP**：如果需要更高灵敏度或已有 BLAST 数据库

配置文件中设置：
```yaml
tool: "diamond"  # 或 "blastp"
```

## 性能调优建议

### 小数据集 (<1GB)

```bash
python orpheus_cli.py -i small.fasta -t 8
```

### 中等数据集 (1-10GB)

```bash
python orpheus_cli.py -i medium.fasta -t 24
```

### 大数据集 (>10GB)

```bash
python orpheus_cli.py -i large.fasta -t 48
```

### 服务器环境推荐

```bash
# 使用 75% 可用核心数
# 例如 64 核服务器
python orpheus_cli.py -i trinity.fasta -t 48
```

## 更多信息

- 详细文档：[README.md](../README.md)
- 流程说明：[PIPELINE_FLOW.md](PIPELINE_FLOW.md)
- 同源搜索：[HOMOLOGY_SEARCH.md](HOMOLOGY_SEARCH.md)
- 配置示例：[example_config.yaml](../example_config.yaml)

