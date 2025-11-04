# 新功能：转录本综合评分和筛选系统

## 功能概述

Orpheus v0.2.0 新增了转录本综合评分系统，实现了您提出的需求：

> "我们的步骤应该可以从BUSCO开始的，输入gff文件与CD-HIT去冗余文件，输出每一个开放阅读框的完整程度，这一条要和有没有预测出同源性那一步综合投票，决出最可能的转录本"

## 主要特性

### 1. 综合评分机制

整合三种证据类型：
- **ORF 完整性**（权重 40%）：complete ORF 优先
- **ORF 长度**（权重 30%）：长 ORF 更可靠
- **同源性**（权重 30%）：有同源蛋白支持优先

### 2. 智能筛选

- 按基因分组
- 每个基因保留评分最高的转录本（可配置）
- 自动质量分级（high/medium/low）

### 3. 灵活配置

```yaml
scoring:
  completeness_weight: 0.4
  length_weight: 0.3
  homology_weight: 0.3
  min_score_threshold: 0.5
  max_transcripts_per_gene: 1  # 每个基因保留的转录本数
  max_orf_length: 3000
```

## 技术实现

### 新增模块

1. **orpheus/tools/scorer.py**
   - `TranscriptScorer` 类
   - 解析 GFF3 文件
   - 解析同源搜索结果
   - 计算综合评分
   - 筛选和排序转录本

2. **Pipeline 集成**
   - 新增 `busco_scoring` 步骤
   - 自动查找 GFF3 和 CD-HIT 文件
   - 支持从此步骤开始执行

### 工作流程

```
TransDecoder GFF3 → 解析 ORF 信息
                          ↓
同源搜索结果      → 解析同源性信息
                          ↓
CD-HIT 文件       → 识别基因分组
                          ↓
                     综合评分
                          ↓
                   筛选和排序
                          ↓
                  输出 TSV 报告
```

## 使用示例

### 完整流程

```bash
# 运行所有步骤，包括评分
orpheus -i trinity.fasta -t 16
```

### 单独运行评分

```bash
# 如果已有 GFF3 和 CD-HIT 文件
orpheus -i trinity.fasta --start-from busco_scoring
```

### 自定义权重

创建配置文件 `my_config.yaml`：

```yaml
scoring:
  completeness_weight: 0.5  # 更重视完整性
  length_weight: 0.3
  homology_weight: 0.2
  max_transcripts_per_gene: 2  # 每个基因保留前2个
```

运行：
```bash
orpheus -i trinity.fasta -c my_config.yaml
```

## 输出文件

### scored_transcripts.tsv

```tsv
transcript_id	gene_id	orf_type	orf_length	homology	score	quality
TRINITY_DN100_c0_g1_i1	TRINITY_DN100_c0_g1	complete	1200	yes	1.00	high
TRINITY_DN200_c0_g1_i1	TRINITY_DN200_c0_g1	complete	2400	yes	1.00	high
TRINITY_DN300_c0_g1_i1	TRINITY_DN300_c0_g1	complete	900	no	0.70	high
TRINITY_DN400_c0_g1_i2	TRINITY_DN400_c0_g1	5prime_partial	800	no	0.38	low
```

### 列说明

- **transcript_id**: 转录本 ID
- **gene_id**: 基因 ID（从 Trinity ID 提取）
- **orf_type**: ORF 类型（complete/5prime_partial/3prime_partial/internal）
- **orf_length**: ORF 长度（bp）
- **homology**: 是否有同源证据（yes/no）
- **score**: 综合评分（0.0-1.0）
- **quality**: 质量等级（high≥0.7, medium 0.5-0.7, low<0.5）

## 评分规则

### ORF 完整性评分

```python
if orf_type == "complete":
    completeness_score = 1.0
else:
    completeness_score = 0.0
```

### ORF 长度评分

```python
length_score = min(orf_length / max_orf_length, 1.0)
# 默认 max_orf_length = 3000 bp
```

### 同源性评分

```python
if has_homology:
    homology_score = 1.0
else:
    homology_score = 0.0
```

### 总分计算

```python
total_score = (
    completeness_score × completeness_weight +
    length_score × length_weight +
    homology_score × homology_weight
)
```

## 筛选策略

### 基因级别筛选

1. 按基因分组（gene_id）
2. 每组内按评分降序排序
3. 保留前 N 个（默认 N=1）

### 质量阈值筛选

所有转录本都会被评分，但低于阈值的会被标记为 `low` 质量：

```python
if score >= 0.7:
    quality = "high"
elif score >= min_score_threshold:  # 默认 0.5
    quality = "medium"
else:
    quality = "low"
```

## 应用场景

### 场景 1：PCR 克隆

**需求**：选择完整、可靠的转录本用于克隆

**配置**：
```yaml
scoring:
  completeness_weight: 0.6  # 高完整性要求
  length_weight: 0.2
  homology_weight: 0.2
  min_score_threshold: 0.7  # 高标准
```

**筛选**：
```bash
# 只选择 complete + high quality
awk '$3=="complete" && $7=="high"' scored_transcripts.tsv
```

### 场景 2：功能注释

**需求**：优先考虑有同源性的转录本

**配置**：
```yaml
scoring:
  completeness_weight: 0.3
  length_weight: 0.2
  homology_weight: 0.5  # 高同源性要求
```

**筛选**：
```bash
# 有同源证据的转录本
awk '$5=="yes"' scored_transcripts.tsv
```

### 场景 3：可变剪接研究

**需求**：保留每个基因的多个转录本

**配置**：
```yaml
scoring:
  max_transcripts_per_gene: 3  # 保留前3个
```

## 技术细节

### GFF3 解析

提取信息：
- 转录本 ID
- ORF 类型（从 `type:` 字段）
- ORF 长度（计算 CDS 区域总长度）

### 同源搜索结果解析

支持格式：
- Diamond outfmt6
- BLASTP outfmt6

提取：有显著匹配的转录本列表

### 基因 ID 提取

从 Trinity ID 提取基因级别标识：
```python
# TRINITY_DN100_c0_g1_i1 → TRINITY_DN100_c0_g1
gene_id = "_".join(transcript_id.split("_")[:-1])
```

## 性能特点

- **快速**：纯 Python 实现，文件解析高效
- **内存友好**：逐行解析，不一次性加载全部数据
- **健壮**：完善的错误处理和日志记录

## 未来改进

可能的扩展方向：

1. **表达量权重**：整合 Salmon/Kallisto 的定量信息
2. **覆盖度评分**：考虑序列覆盖度均匀性
3. **长读长支持**：整合 PacBio/Nanopore 全长证据
4. **可变权重**：根据不同基因类型使用不同权重
5. **交互式可视化**：生成 HTML 报告

## 相关文档

- [转录本评分详细文档](docs/TRANSCRIPT_SCORING.md)
- [快速参考手册](docs/QUICK_REFERENCE.md)
- [主文档](README.md)

## 更新日志

### v0.2.0-dev
- ✨ 新增：转录本综合评分系统
- ✨ 新增：基于评分的智能筛选
- ✨ 新增：`busco_scoring` 步骤
- ✨ 新增：质量分级（high/medium/low）
- 📝 新增：评分系统详细文档

