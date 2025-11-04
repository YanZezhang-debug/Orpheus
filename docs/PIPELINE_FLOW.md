# Orpheus 流程架构

## 流程图

```
输入: Trinity 组装结果
         |
         v
┌────────────────────┐
│   CD-HIT 去冗余     │ <-- 可以从这里开始 (--start-from cdhit)
│  (去除相似序列)     │
└────────────────────┘
         |
         v
   去冗余的转录本
         |
         v
┌────────────────────┐
│ TransDecoder ORF   │ <-- 可以从这里开始 (--start-from transdecoder)
│      预测          │
│ (识别编码区域)      │
└────────────────────┘
         |
         v
    所有候选 ORF
  (longest_orfs.pep)
         |
         v
┌────────────────────┐
│  完整 ORF 筛选     │ (可选，如果 use_complete_only=true)
│ (type:complete)    │
└────────────────────┘
         |
         v
  完整 ORF 子集
(longest_orfs.complete.pep)
         |
         v
┌────────────────────┐
│   同源搜索 (可选)   │
│ Diamond/BLASTP     │
└────────────────────┘
         |
         v
    同源证据文件
         |
         v
┌────────────────────┐
│ TransDecoder.Predict│
│   (最终预测)       │
└────────────────────┘
         |
         v
  最终 CDS/蛋白预测
   (.pep/.cds/.gff3)
         |
         v
┌────────────────────┐
│  BUSCO 质量评估    │ <-- 可以从这里开始 (--start-from busco)
│  (核心基因完整性)   │
└────────────────────┘
         |
         v
  BUSCO 评估结果
 (完整/重复/片段/缺失)
         |
         v
┌────────────────────┐
│   转录本评分       │ <-- 可以从这里开始 (--start-from scoring)
│ (综合质量评估)     │
└────────────────────┘
         |
         v
  评分结果文件
(scored_transcripts.tsv)
         |
         v
   最终筛选结果
```

## 步骤详细说明

### 步骤 1: CD-HIT 去冗余 (`cdhit`)

**目的**：去除高度相似的冗余转录本

**输入**：
- Trinity 组装的原始转录本 FASTA 文件

**输出**：
- 去冗余后的转录本 FASTA 文件
- 聚类信息文件 (.clstr)

**关键参数**：
- `identity`: 序列相似度阈值 (默认 0.95)
- `coverage`: 覆盖度阈值 (默认 0.9)
- `threads`: 并行线程数

**统计信息**：
- 输入序列数量
- 输出序列数量
- 去冗余率

### 步骤 2: TransDecoder ORF 预测 (`transdecoder`)

**目的**：预测开放阅读框（ORF）并评估完整度

**输入**：
- CD-HIT 输出的去冗余转录本（或任何转录本文件）

**输出**：
- `.transdecoder.pep`: 预测的蛋白质序列
- `.transdecoder.cds`: 预测的 CDS 序列
- `.transdecoder.gff3`: GFF3 注释文件（含完整度信息）

**ORF 完整度分类**：
- `complete`: 含起始密码子和终止密码子
- `5prime_partial`: 缺少起始密码子
- `3prime_partial`: 缺少终止密码子
- `internal`: 两端都不完整

**关键参数**：
- `min_protein_length`: 最小蛋白长度（氨基酸数）
- `genetic_code`: 遗传密码表
- `single_best_only`: 是否只保留最佳 ORF

**同源搜索配置**（提高预测准确性）：
- `homology_search.enabled`: 是否启用同源搜索（强烈推荐）
- `homology_search.use_complete_only`: 是否仅用完整 ORF 进行同源搜索（推荐，对 PCR 实验重要）
- `homology_search.tool`: 搜索工具（diamond 或 blastp）
- `homology_search.database`: 蛋白质数据库路径（如 Swiss-Prot）
- `homology_search.evalue`: E-value 阈值

**工作流程**：
1. **ORF 识别**: `TransDecoder.LongOrfs` 识别所有候选 ORF
   - 扫描所有六个阅读框
   - 根据 `min_protein_length` 筛选候选 ORF
   - 输出 `longest_orfs.pep` (所有候选 ORF)
   
2. **完整 ORF 筛选**（如果启用 `use_complete_only`）：
   - 自动从 `longest_orfs.pep` 中筛选完整 ORF
   - 筛选标准：序列头包含 `type:complete`
   - 输出 `longest_orfs.complete.pep`（仅完整 ORF）
   - 提供统计信息：总 ORF 数、完整 ORF 数和比例
   - **重要性**：完整 ORF 对 PCR 实验设计至关重要，确保能扩增完整基因
   
3. **同源搜索**（可选但推荐）：
   - 如果 `use_complete_only=True`：仅使用完整 ORF 进行同源比对
   - 使用 Diamond/BLASTP 在数据库中搜索同源蛋白
   - 生成同源证据文件 (`blastp.outfmt6` 或 `diamond.outfmt6`)
   
4. **最终预测**: `TransDecoder.Predict` 结合多种证据
   - 整合长度标准
   - 整合同源搜索证据
   - 输出最终的 `.pep`、`.cds` 和 `.gff3` 文件

**统计信息**：
- 预测的候选 ORF 总数
- 完整 ORF 数量和比例
- 不完整 ORF 的分类统计（5' partial, 3' partial, internal）
- 同源搜索匹配数量（如果启用）
- 最终保留的 CDS 数量

### 步骤 3: BUSCO 质量评估 (`busco`)

**目的**：评估转录组的核心基因完整性

**输入**：
- TransDecoder 预测的蛋白质序列 (`.transdecoder.pep`)

**输出**：
- BUSCO 评估报告 (`short_summary.*.txt`)
- 完整基因列表
- 重复基因列表
- 片段基因列表
- 缺失基因列表

**关键参数**：
- `enabled`: 是否启用 BUSCO 评估（默认 false）
- `lineage`: BUSCO 数据集路径或名称（如 `poales_odb12`）
- `mode`: 运行模式（`transcriptome`、`genome` 或 `proteins`）
- `threads`: 线程数

**BUSCO 基因状态分类**：
- **Complete**: 完整的单拷贝基因（最佳质量）
- **Duplicated**: 完整但有多个拷贝（可能是基因家族或组装冗余）
- **Fragmented**: 部分匹配（组装不完整）
- **Missing**: 未找到（可能真的缺失或组装质量问题）

**统计信息**：
- 总 BUSCO 基因数
- 完整基因数和百分比
- 重复基因数和百分比
- 片段基因数和百分比
- 缺失基因数和百分比

**配置示例**：
```yaml
busco:
  enabled: true
  lineage: "/mnt/data/01.databases/poales_odb12"
  mode: "transcriptome"
  threads: 8
```

### 步骤 4: 转录本评分 (`scoring`)

**目的**：整合多种证据，为每个转录本计算综合质量评分

**输入**：
- TransDecoder GFF3 文件（ORF 完整性信息）
- 同源搜索结果（如果启用）
- BUSCO 结果（如果启用）

**输出**：
- `scored_transcripts.tsv`: 包含所有转录本的详细评分
- `high_quality_transcripts.fasta`: 高质量转录本序列（可选）

**评分依据**：
1. **BUSCO 基因匹配** (权重 0.4)
   - 匹配到 BUSCO 核心基因的转录本得最高分
   - 这些是最可靠的转录本
   
2. **ORF 完整性** (权重 0.3)
   - `complete`: 1.0
   - `5prime_partial`: 0.5
   - `3prime_partial`: 0.5
   - `internal`: 0.0
   
3. **同源性** (权重 0.2)
   - 有同源匹配: 1.0
   - 无同源匹配: 0.0
   
4. **ORF 长度** (权重 0.1)
   - 归一化长度（基于 `max_orf_length` 参数）

**关键参数**：
- `weights`: 各项评分的权重配置
- `min_score_threshold`: 最小质量阈值（默认 0.5）
- `max_transcripts_per_gene`: 每个基因保留的最大转录本数

**质量等级**：
- **high**: score ≥ 0.7（推荐用于下游分析）
- **medium**: 0.5 ≤ score < 0.7（可以保留）
- **low**: score < 0.5（建议过滤）

**输出文件列**：
- `transcript_id`: 转录本 ID
- `gene_id`: 基因 ID
- `orf_type`: ORF 完整性类型
- `orf_length`: ORF 长度（bp）
- `homology`: 是否有同源证据
- `busco_gene`: BUSCO 基因名（如有）
- `busco_status`: BUSCO 基因状态（如有）
- `score`: 综合评分
- `quality`: 质量等级

**配置示例**：
```yaml
scoring:
  weights:
    busco: 0.4
    completeness: 0.3
    homology: 0.2
    length: 0.1
  min_score_threshold: 0.5
  max_transcripts_per_gene: 1
```

详细说明请参考：[TRANSCRIPT_SCORING.md](TRANSCRIPT_SCORING.md)

### 未来步骤 (计划中)

#### 覆盖度分析 (`coverage`)

**目的**：评估短读覆盖度的均匀性

**关键输出**：
- 覆盖度分布
- 均匀性指标

## 步骤控制

### 从头开始运行所有步骤

```bash
orpheus -i trinity.fasta
# 等价于
orpheus -i trinity.fasta --start-from beginning
```

执行步骤：cdhit → transdecoder → busco → scoring

### 从 CD-HIT 开始

```bash
orpheus -i trinity.fasta --start-from cdhit
```

执行步骤：cdhit → transdecoder → busco → scoring

### 仅运行 TransDecoder

```bash
orpheus -i cdhit_result.fasta --start-from transdecoder
```

执行步骤：transdecoder → busco → scoring

### 从 BUSCO 开始

```bash
orpheus -i trinity.fasta --start-from busco
```

执行步骤：busco → scoring

（需要已经运行过 TransDecoder，生成了 `.transdecoder.pep` 文件）

### 仅运行评分

```bash
orpheus -i trinity.fasta --start-from scoring
```

执行步骤：scoring

（需要已经有 TransDecoder 输出，可选 BUSCO 结果）

### 使用场景

| 场景 | 命令 | 说明 |
|------|------|------|
| 首次运行 | `--start-from beginning` | 执行所有步骤 |
| CD-HIT 完成，想重跑 TransDecoder | `--start-from transdecoder` | 跳过 CD-HIT |
| 测试不同 TransDecoder 参数 | `--start-from transdecoder` | 避免重复运行 CD-HIT |
| 调整 BUSCO 参数重新评估 | `--start-from busco` | 跳过前面的步骤 |
| 调整评分权重重新评分 | `--start-from scoring` | 只重新评分，不重跑分析 |
| 某步骤失败后继续 | `--start-from <失败的步骤>` | 从失败处继续 |

## 扩展新步骤

Orpheus 的架构支持轻松添加新步骤：

1. **定义步骤**：在 `PIPELINE_STEPS` 列表中添加
2. **实现方法**：添加 `run_<步骤名>()` 方法
3. **步骤调度**：在 `run()` 方法中添加调用逻辑
4. **更新文档**：更新帮助信息和文档

示例：

```python
# orpheus/pipeline.py
PIPELINE_STEPS = [
    'cdhit',
    'transdecoder',
    'busco',
    'scoring',
    'coverage',  # 新步骤
]

def run_coverage(self, input_file=None, output_file=None):
    """运行覆盖度分析"""
    # 实现逻辑
    pass
```

CLI 自动支持：

```bash
orpheus -i input.fasta --start-from coverage
```

## 中间文件管理

每个步骤的输出默认保存在工作目录（`orpheus_output/`）：

```
orpheus_output/
├── cdhit_result.fasta                      # CD-HIT 输出
├── cdhit_result.fasta.clstr                # CD-HIT 聚类信息
├── transdecoder_results/                   # TransDecoder 输出目录
│   ├── input.fasta.transdecoder_dir/       # TransDecoder 工作目录
│   │   ├── longest_orfs.pep                # 所有候选 ORF
│   │   ├── longest_orfs.complete.pep       # 完整 ORF（如果启用）
│   │   ├── longest_orfs.cds                # 候选 CDS
│   │   ├── blastp.outfmt6                  # BLASTP 同源结果
│   │   └── diamond.outfmt6                 # Diamond 同源结果
│   ├── input.fasta.transdecoder.pep        # 最终预测的蛋白质序列
│   ├── input.fasta.transdecoder.cds        # 最终预测的 CDS 序列
│   └── input.fasta.transdecoder.gff3       # 最终 GFF3 注释文件
├── busco_results/                          # BUSCO 输出目录（如果启用）
│   ├── short_summary.*.txt                 # BUSCO 摘要报告
│   ├── full_table.tsv                      # 完整结果表
│   ├── missing_busco_list.tsv              # 缺失基因列表
│   └── run_*/                              # BUSCO 详细结果
│       ├── busco_sequences/
│       │   ├── single_copy_busco_sequences/
│       │   └── multi_copy_busco_sequences/
│       └── hmmer_output/
├── scored_transcripts.tsv                  # 转录本评分结果
├── high_quality_transcripts.fasta          # 高质量转录本（可选）
└── orpheus.log                             # 日志文件
```

### 文件说明

**CD-HIT 输出**：
- `cdhit_result.fasta`: 去冗余后的转录本序列
- `cdhit_result.fasta.clstr`: 聚类信息文件

**TransDecoder 中间文件**：
- `longest_orfs.pep`: TransDecoder.LongOrfs 识别的所有候选 ORF
- `longest_orfs.complete.pep`: 筛选出的完整 ORF（仅当 `use_complete_only=True` 时生成）
- `blastp.outfmt6` / `diamond.outfmt6`: 同源搜索结果（仅当启用同源搜索时生成）

**TransDecoder 最终输出**：
- `.transdecoder.pep`: 最终预测的蛋白质序列
- `.transdecoder.cds`: 最终预测的 CDS 序列  
- `.transdecoder.gff3`: 包含完整度信息的注释文件

**BUSCO 输出**：
- `short_summary.*.txt`: BUSCO 评估摘要（包含完整性百分比）
- `full_table.tsv`: 每个 BUSCO 基因的详细匹配信息
- `busco_sequences/`: 匹配到的 BUSCO 基因序列

**评分输出**：
- `scored_transcripts.tsv`: 包含所有转录本的综合评分（包含 BUSCO、完整性、同源性、长度评分）
- `high_quality_transcripts.fasta`: 高质量转录本序列（可选）

## 步骤跳过与文件自动查找

### 文件查找逻辑

当使用 `--start-from` 跳过某些步骤时，Orpheus 会**自动查找**工作目录中的中间文件，无需手动指定。

#### TransDecoder 步骤的输入文件查找

当执行 `--start-from transdecoder` 时，输入文件按以下优先级查找：

1. **命令行指定的输入文件** (`-i` 参数)
   ```bash
   orpheus -i my_sequences.fasta --start-from transdecoder
   ```

2. **工作目录中的 CD-HIT 结果**（自动查找）
   ```
   ./orpheus_output/cdhit_result.fasta
   ```
   如果找到此文件，会自动使用它作为输入

3. **配置文件中的原始输入**
   ```yaml
   io:
     trinity_assembly: /path/to/trinity.fasta
   ```

**示例 1：使用已有的 CD-HIT 结果**
```bash
# 第一次运行：完整流程
orpheus -i trinity.fasta

# 第二次运行：跳过 CD-HIT，使用之前的结果
orpheus --start-from transdecoder
# ✓ 自动找到 ./orpheus_output/cdhit_result.fasta
```

**示例 2：指定其他输入文件**
```bash
# 使用自定义的 FASTA 文件
orpheus -i my_custom.fasta --start-from transdecoder
```

#### BUSCO Scoring 步骤的文件查找

当执行 `--start-from busco_scoring` 时，需要以下文件：

**1. GFF3 文件**（TransDecoder 输出，必需）

查找优先级：
- 优先使用：上一步 TransDecoder 的输出 (`self.transdecoder_output_dir/*.transdecoder.gff3`)
- 自动查找：`./orpheus_output/transdecoder_results/*.transdecoder.gff3`

**2. CD-HIT FASTA 文件**（必需）

查找优先级：
- 优先使用：上一步 CD-HIT 的输出 (`self.cdhit_output_file`)
- 自动查找：`./orpheus_output/cdhit_result.fasta`

**3. TransDecoder 目录**（自动识别）

当找到 GFF3 文件时，会自动查找对应的 TransDecoder 工作目录：
```
./orpheus_output/transdecoder_results/cdhit_result.fasta.transdecoder_dir/
```

此目录包含：
- `longest_orfs.pep`: 所有候选 ORF
- `longest_orfs.complete.pep`: 完整 ORF（如果启用了 `use_complete_only`）
- 同源搜索结果（`blastp.outfmt6` 或 `diamond.outfmt6`）

**示例 1：直接运行 BUSCO Scoring**
```bash
# 确保工作目录包含必要文件：
# ./orpheus_output/cdhit_result.fasta
# ./orpheus_output/transdecoder_results/*.transdecoder.gff3
# ./orpheus_output/transdecoder_results/*.transdecoder_dir/

orpheus --start-from busco_scoring
# ✓ 自动找到 GFF3 文件
# ✓ 自动找到 CD-HIT 文件
# ✓ 自动找到 TransDecoder 目录
```

**示例 2：使用不同工作目录的结果**
```bash
# 假设之前的结果在 /path/to/results/
orpheus --work-dir /path/to/results --start-from busco_scoring
# ✓ 在 /path/to/results/orpheus_output/ 中查找所有必需文件
```

### 错误处理

如果找不到必需的文件，会显示清晰的错误信息：

```
错误: 输入文件不存在: None
请确保以下之一存在:
  1. 工作目录中的 CD-HIT 结果: ./orpheus_output/cdhit_result.fasta
  2. 配置文件中指定的输入文件
  3. 通过命令行 -i 参数指定的输入文件
```

### 最佳实践

1. **保持工作目录结构**：
   - 不要移动或删除中间文件（`cdhit_result.fasta` 等）
   - 使用默认的输出目录结构

2. **使用相同的工作目录**：
   ```yaml
   io:
     work_dir: ./orpheus_output  # 保持一致
   ```

3. **重跑特定步骤**：
   ```bash
   # 重跑 TransDecoder，使用之前的 CD-HIT 结果
   orpheus --start-from transdecoder
   
   # 重跑评分，使用之前的所有结果
   orpheus --start-from busco_scoring
   ```

## 完整 ORF 筛选详解

### 为什么需要完整 ORF？

在使用同源搜索辅助 ORF 预测时，**完整 ORF**（含起始和终止密码子）比不完整 ORF 更可靠：

1. **生物学意义**：
   - 完整 ORF 代表完整的蛋白质编码序列
   - 包含完整的功能域，同源比对更准确
   
2. **实验应用**：
   - **PCR 引物设计**：需要知道完整的起始和终止位置
   - **基因克隆**：完整 ORF 可以直接用于表达
   - **功能研究**：不完整序列可能缺失关键功能域

3. **预测质量**：
   - 完整 ORF 通常代表高质量的转录本组装
   - 减少假阳性和低质量匹配
   - 提高同源搜索的特异性

### 筛选流程

当配置 `homology_search.use_complete_only: true` 时：

1. TransDecoder.LongOrfs 生成所有候选 ORF → `longest_orfs.pep`
2. 自动筛选完整 ORF（`type:complete`）→ `longest_orfs.complete.pep`
3. 仅使用完整 ORF 进行同源搜索
4. 将同源证据反馈给 TransDecoder.Predict

### 完整度分类

TransDecoder 将 ORF 分为四类：

| 类型 | 描述 | 起始密码子 | 终止密码子 | 适合同源搜索 |
|------|------|-----------|-----------|------------|
| `complete` | 完整 ORF | ✓ | ✓ | **推荐** |
| `5prime_partial` | 5' 端不完整 | ✗ | ✓ | 不推荐 |
| `3prime_partial` | 3' 端不完整 | ✓ | ✗ | 不推荐 |
| `internal` | 两端都不完整 | ✗ | ✗ | 不推荐 |

### 配置示例

```yaml
transdecoder:
  min_protein_length: 100
  homology_search:
    enabled: true
    use_complete_only: true  # 仅使用完整 ORF
    tool: diamond
    database: /path/to/swissprot
    evalue: 1e-5
```

### 统计输出示例

```
✓ 完整 ORF 筛选完成:
  总 ORF 数: 15234
  完整 ORF: 8956 (58.78%)
  不完整 ORF (已排除): 6278
```

**解读**：
- 高比例（>50%）的完整 ORF 表明转录本组装质量较好
- 低比例（<30%）可能提示需要优化组装参数或数据质量

