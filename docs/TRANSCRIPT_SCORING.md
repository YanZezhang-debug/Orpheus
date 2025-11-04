# 转录本评分和筛选功能

## 概述

Orpheus v0.2.0+ 提供了综合评分系统，用于评估转录本质量并筛选出最可信的转录本。该系统整合了多种证据：**BUSCO 保守基因匹配**、**ORF 完整性**、**同源性**和 **ORF 长度**，通过加权评分为每个转录本计算质量得分。

## 评分依据

### 1. BUSCO 基因匹配（默认权重：40%）⭐ 最强证据

**BUSCO**（Benchmarking Universal Single-Copy Orthologs）检测保守的单拷贝直系同源基因。如果一个转录本匹配到 BUSCO 基因，说明它是重要的核心功能基因。

- **匹配 BUSCO 基因**（Complete/Duplicated 状态）
  - `busco_score = 1.0`
  - ✅✅✅ **最强证据**：高度保守的核心基因
  - 这些基因在进化过程中高度保守，通常以单拷贝形式存在
  
- **未匹配 BUSCO 基因**
  - `busco_score = 0.0`
  - ℹ️ 不意味着质量差，可能是物种特异性基因或非核心基因

**注意**：
- BUSCO 完整性评估的是**整个转录组**的质量（有多少核心基因被组装出来）
- BUSCO 基因匹配评估的是**单个转录本**是否为核心保守基因
- 只有 `Complete` 和 `Duplicated` 状态的基因被计入（`Fragmented` 基因不计入）

### 2. ORF 完整性（默认权重：30%）

完整的 ORF（同时具有起始密码子和终止密码子）更可能是真实的编码区域：

- **complete**（完整）：起始密码子 + CDS + 终止密码子
  - `completeness_score = 1.0`
  - ✅ 推荐用于 PCR 实验
  
- **5prime_partial**（5' 端不完整）：缺少起始密码子
  - `completeness_score = 0.6`
  - ⚠️ 可能是组装不完整或真实的 N 端截短蛋白
  
- **3prime_partial**（3' 端不完整）：缺少终止密码子
  - `completeness_score = 0.6`
  - ⚠️ 可能是组装不完整
  
- **internal**（内部片段）：两端都不完整
  - `completeness_score = 0.3`
  - ❌ 不推荐使用

### 3. 同源性（默认权重：20%）

有同源蛋白支持的 ORF 更可能是真实的：

- **有同源证据**：在蛋白质数据库中找到显著匹配
  - `homology_score = 1.0`
  - ✅ 强证据支持，功能可推断
  
- **无同源证据**：未找到显著匹配
  - `homology_score = 0.0`
  - ⚠️ 可能是新基因、物种特异性基因或假基因

### 4. ORF 长度（默认权重：10%）

较长的 ORF 通常更可靠，因为随机形成的长 ORF 概率极低：

- 长度归一化公式：`length_score = orf_length / max_orf_length`
- `max_orf_length` 为数据集中最长的 ORF
- 示例（假设最长 ORF 为 3000 bp）：
  - 1500 bp → score = 0.5
  - 3000 bp → score = 1.0

## 评分公式

### 默认权重（有 BUSCO 数据时）

```
总分 = busco_score × 0.4 + completeness_score × 0.3 + homology_score × 0.2 + length_score × 0.1
```

### 备用权重（无 BUSCO 数据时）

```
总分 = completeness_score × 0.5 + homology_score × 0.3 + length_score × 0.2
```

**分数范围**：0.0 - 1.0

## 质量分级

根据总分自动分级：

| 质量等级 | 分数范围 | 说明 |
|---------|---------|------|
| **high** | ≥ 0.7 | 高质量转录本，推荐使用 |
| **medium** | 0.5 - 0.7 | 中等质量，可能需要验证 |
| **low** | < 0.5 | 低质量，不推荐使用 |

## 使用方法

### 方法一：运行完整流程（推荐）

```bash
# 运行包含评分步骤的完整流程
orpheus -i trinity_assembly.fasta -t 16
```

这将依次执行：
1. CD-HIT 去冗余
2. TransDecoder ORF 预测（含同源搜索）
3. BUSCO ORF 评分和转录本筛选

### 方法二：单独运行评分步骤

如果您已经运行了前面的步骤，可以只运行评分：

```bash
# 从评分步骤开始
orpheus -i trinity_assembly.fasta --start-from busco_scoring
```

**注意**：此时需要确保以下文件已存在：
- CD-HIT 输出文件（默认：`orpheus_output/cdhit_result.fasta`）
- TransDecoder GFF3 文件（默认：`orpheus_output/transdecoder_results/*.transdecoder.gff3`）

## 配置选项

在配置文件中自定义评分参数：

```yaml
scoring:
  # 评分权重（总和应为 1.0）
  weights:
    busco: 0.4          # BUSCO 基因匹配权重
    completeness: 0.3   # ORF 完整性权重
    homology: 0.2       # 同源性权重
    length: 0.1         # ORF 长度权重
  
  # 筛选参数
  threshold: 0.5        # 最小分数阈值
  top_n: null           # 保留得分最高的 N 个（null = 不限制）
  
  # 同源搜索结果文件（可选）
  homology_file: null   # DIAMOND/BLASTP 结果文件路径
```

### 权重调整建议

根据您的研究目标调整权重：

#### 场景 1：核心基因研究（推荐）
```yaml
scoring:
  weights:
    busco: 0.5          # 最重视进化保守的核心基因
    completeness: 0.3
    homology: 0.15
    length: 0.05
```

#### 场景 2：PCR 克隆优先
```yaml
scoring:
  weights:
    busco: 0.3
    completeness: 0.5   # 更重视完整性（需要完整 CDS）
    homology: 0.15
    length: 0.05
```

#### 场景 3：功能注释优先
```yaml
scoring:
  weights:
    busco: 0.3
    completeness: 0.2
    homology: 0.4       # 更重视同源性（便于功能推断）
    length: 0.1
```

#### 场景 4：新基因发现（无 BUSCO 约束）
```yaml
scoring:
  weights:
    busco: 0.1
    completeness: 0.5   # 更重视完整性
    homology: 0.1       # 降低同源性要求（允许新基因）
    length: 0.3         # 长度作为主要可靠性指标
```

## 输出文件

### 评分结果文件（scored_transcripts.tsv）

TSV 格式，包含以下列：

| 列名 | 说明 | 示例值 |
|-----|------|--------|
| transcript_id | 转录本 ID | TRINITY_DN100_c0_g1_i1 |
| gene_id | 基因 ID | TRINITY_DN100_c0_g1 |
| orf_type | ORF 类型 | complete, 5prime_partial, etc. |
| orf_length | ORF 长度（bp） | 1200 |
| homology | 是否有同源证据 | yes / no |
| busco_gene | BUSCO 基因名（如有） | 12345at33090 或空 |
| busco_status | BUSCO 基因状态 | Complete / Duplicated 或空 |
| score | 总评分 | 0.85 |
| quality | 质量等级 | high / medium / low |

### 示例输出

```tsv
transcript_id	gene_id	orf_type	orf_length	homology	busco_gene	busco_status	score	quality
TRINITY_DN100_c0_g1_i1	TRINITY_DN100_c0_g1	complete	1200	yes	12345at33090	Complete	1.00	high
TRINITY_DN100_c0_g1_i2	TRINITY_DN100_c0_g1	5prime_partial	800	no			0.38	low
TRINITY_DN200_c0_g1_i1	TRINITY_DN200_c0_g1	complete	2400	yes	67890at33090	Duplicated	1.00	high
TRINITY_DN300_c0_g1_i1	TRINITY_DN300_c0_g1	complete	900	no			0.52	medium
```

**说明**：
- 匹配到 BUSCO 基因的转录本会显示 `busco_gene` 和 `busco_status`
- 未匹配的转录本这两列为空
- BUSCO 基因通常是核心保守基因，具有很高的可信度

## 转录本筛选

### 基因级别筛选

对于每个基因（gene_id），可以选择保留：

1. **最高分转录本**（默认）：`max_transcripts_per_gene: 1`
   - 每个基因只保留评分最高的 1 个转录本
   - 适合：需要非冗余转录本集、PCR 验证
   
2. **前 N 个转录本**：`max_transcripts_per_gene: N`
   - 每个基因保留评分前 N 高的转录本
   - 适合：研究可变剪接
   
3. **所有转录本**：`max_transcripts_per_gene: 0`
   - 保留所有转录本及其评分
   - 适合：手动筛选、进一步分析

### 质量阈值筛选

低于 `min_score_threshold` 的转录本会被标记为 `low` 质量：

```yaml
min_score_threshold: 0.5  # 默认值
```

建议设置：
- **严格筛选**：0.7（只保留高质量）
- **平衡筛选**：0.5（默认）
- **宽松筛选**：0.3（保留更多候选）

## 使用场景示例

### 场景 1：选择用于 PCR 克隆的转录本

**目标**：找到完整、可靠的转录本用于克隆

**配置**：
```yaml
scoring:
  completeness_weight: 0.6
  length_weight: 0.2
  homology_weight: 0.2
  min_score_threshold: 0.7
  max_transcripts_per_gene: 1
```

**筛选标准**：
```bash
# 从结果中选择 complete + high quality
awk '$3=="complete" && $7=="high"' scored_transcripts.tsv
```

### 场景 2：研究可变剪接

**目标**：保留每个基因的多个转录本

**配置**：
```yaml
scoring:
  max_transcripts_per_gene: 3  # 保留前 3 个
```

### 场景 3：功能注释

**目标**：优先考虑有同源性的转录本

**配置**：
```yaml
scoring:
  completeness_weight: 0.3
  length_weight: 0.2
  homology_weight: 0.5
  min_score_threshold: 0.5
```

**筛选标准**：
```bash
# 优先使用有同源证据的转录本
awk '$5=="yes"' scored_transcripts.tsv | sort -k6 -rn
```

## 评分统计

运行评分后，日志会显示统计信息：

```
转录本评分统计:
  总转录本数: 50000
  complete ORF: 35000 (70.0%)
  5prime_partial: 8000 (16.0%)
  3prime_partial: 5000 (10.0%)
  internal: 2000 (4.0%)
  
  有同源证据: 30000 (60.0%)
  
  质量分布:
    high (≥0.7): 28000 (56.0%)
    medium (0.5-0.7): 15000 (30.0%)
    low (<0.5): 7000 (14.0%)
  
  基因级别筛选:
    原始基因数: 20000
    保留转录本数: 20000
    平均每基因转录本数: 1.0
```

## 注意事项

1. **同源搜索是可选的**，但强烈推荐启用
   - 未启用时，`homology_score` 对所有转录本都是 0
   - 评分主要依赖完整性和长度

2. **权重总和应为 1.0**
   - 系统会自动归一化，但建议手动确保总和为 1.0

3. **ORF 长度是核苷酸数，不是氨基酸数**
   - 1200 bp ≈ 400 aa

4. **评分不是绝对标准**
   - 低分转录本可能是真实的（如短肽、新基因）
   - 高分转录本也需要实验验证

## 相关文档

- [快速参考手册](QUICK_REFERENCE.md)
- [同源搜索文档](HOMOLOGY_SEARCH.md)
- [主文档](../README.md)

