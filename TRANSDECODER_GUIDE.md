# TransDecoder 使用指南

## 什么是 TransDecoder？

TransDecoder 是一个用于识别转录本序列中候选编码区域（CDS）的工具。它可以：
1. 识别长开放阅读框（ORF）
2. 预测最可能的编码区域
3. 评估 ORF 的完整度

## Orpheus 中的 TransDecoder 流程

### 自动运行流程

当您运行 Orpheus 时，会自动执行以下步骤：

```
Trinity 组装结果 
    ↓
CD-HIT 去冗余
    ↓
TransDecoder.LongOrfs (识别长 ORF)
    ↓
TransDecoder.Predict (预测编码区域)
    ↓
输出：蛋白序列、CDS 序列、GFF3 注释
```

### 输出文件说明

运行后在 `orpheus_output/transdecoder_results/` 目录下会生成：

1. **`*.transdecoder.pep`** - 预测的蛋白序列（FASTA 格式）
   - 用于后续的同源比对（DIAMOND/BLAST）

2. **`*.transdecoder.cds`** - 预测的 CDS 序列（FASTA 格式）
   - 核苷酸序列，用于 PCR 引物设计

3. **`*.transdecoder.gff3`** - GFF3 格式注释文件 ⭐ 重要！
   - 包含每个 ORF 的位置信息
   - **包含 ORF 完整度信息**

4. **`*.transdecoder.bed`** - BED 格式注释文件
   - 可在基因组浏览器中查看

## ORF 完整度类型

在 GFF3 文件中，查看第 9 列的 `type:` 字段：

### 1. `type:complete` - 完整 ORF ✅
- 有起始密码子（ATG）
- 有终止密码子（TAA/TAG/TGA）
- **ORF_score = 1.0**
- **推荐用于 PCR 克隆**

### 2. `type:5prime_partial` - 5'端部分 ORF ⚠️
- 缺少起始密码子
- 可能是组装不完整或真实的 N端截断
- **ORF_score = 0.0**（除非有强补证据）
- **需要 5' RACE 验证**

### 3. `type:3prime_partial` - 3'端部分 ORF ⚠️
- 缺少终止密码子
- 可能是组装不完整
- **ORF_score = 0.0**（除非有强补证据）
- **需要 3' RACE 验证**

### 4. `type:internal` - 内部片段 ❌
- 两端都缺失
- **ORF_score = 0.0**
- **不推荐用于 PCR**

## 如何查看 ORF 完整度

### 方法 1: 直接查看 GFF3 文件

```bash
# 查看完整 ORF 数量
grep "type:complete" *.transdecoder.gff3 | wc -l

# 查看不完整 ORF 数量
grep "type:5prime_partial" *.transdecoder.gff3 | wc -l
grep "type:3prime_partial" *.transdecoder.gff3 | wc -l
grep "type:internal" *.transdecoder.gff3 | wc -l
```

### 方法 2: 查看 Orpheus 日志

Orpheus 会自动统计并在日志中显示：

```
ORF 完整度统计:
  总 ORF 数: 10000
  完整 ORF: 7500 (75.00%)
  不完整 ORF: 2500 (25.00%)
    - 5'端缺失: 1200
    - 3'端缺失: 1000
    - 内部片段: 300
```

## ORF 评分规则（Orpheus 评估系统）

### 完整 ORF
```
if type:complete
    → ORF_score = 1.0
    → 状态: ✅ Recommended for PCR
```

### 不完整 ORF（无强证据）
```
if type != complete AND no_strong_evidence
    → ORF_score = 0.0
    → 状态: ❌ Not recommended
    → 建议: 不要用于 PCR
```

### 不完整 ORF（有强补证据）
```
if type != complete AND has_strong_evidence
    → ORF_score = 0.5
    → 状态: ⚠️ Manual review required
    → 建议: 只做短片段 PCR (100-250 bp) 或 RACE 验证
```

**强补证据包括**（满足任一项）：
1. **长读长全长支持**
   - ≥1 条长读（PacBio/ONT）跨越大部分转录本
   - 支持当前的拼接结构

2. **强同源命中**
   - DIAMOND/BLAST 比对结果：
     - `align_len / transcript_len ≥ 0.7`
     - `e-value ≤ 1e-20`
     - 命中序列有明确的蛋白域

3. **高且均匀的短读覆盖**
   - `percent_bases_covered ≥ 0.95`
   - `mean_depth ≥ 20`（可配置）

*注：强补证据的检测将在 Orpheus 后续版本中实现*

## 配置 TransDecoder 参数

在 `config/default.yaml` 中调整：

```yaml
transdecoder:
  # 最小蛋白长度（氨基酸数）
  min_protein_length: 100  # 默认 100，可改为 50 以获得更多短 ORF
  
  # 遗传密码表
  genetic_code: "universal"  # 真核生物通常用 universal
                             # 线粒体可用 Mitochondrial-Canonical
  
  # 是否只保留单一最佳 ORF
  single_best_only: false    # 改为 true 可减少冗余，但可能遗漏可变剪切
```

## 实际应用示例

### 场景 1: 选择用于 PCR 克隆的转录本

```bash
# 1. 运行 Orpheus
python3 orpheus_cli.py -i trinity_assembly.fasta

# 2. 提取完整 ORF 的转录本 ID
grep "type:complete" orpheus_output/transdecoder_results/*.gff3 | \
    cut -f1 | sort -u > complete_orf_transcripts.txt

# 3. 根据这个列表筛选转录本用于 PCR
```

### 场景 2: 评估组装质量

完整 ORF 的比例可以反映组装质量：
- **>80%**: 优秀的组装质量
- **60-80%**: 良好的组装质量
- **40-60%**: 中等组装质量，可能需要优化
- **<40%**: 较差的组装质量，建议重新组装

### 场景 3: 识别需要 RACE 验证的转录本

```bash
# 提取 5'端不完整的转录本（需要 5' RACE）
grep "type:5prime_partial" *.gff3 | cut -f1 | sort -u > need_5prime_race.txt

# 提取 3'端不完整的转录本（需要 3' RACE）
grep "type:3prime_partial" *.gff3 | cut -f1 | sort -u > need_3prime_race.txt
```

## 下一步分析

1. **同源比对**（即将在 Orpheus 中实现）
   - 使用 `*.transdecoder.pep` 进行 DIAMOND/BLAST 比对
   - 验证预测的 ORF 是否有同源蛋白支持

2. **表达量分析**（即将在 Orpheus 中实现）
   - 使用 Salmon/Kallisto 计算转录本表达量
   - 高表达转录本更可能是真实基因

3. **功能注释**
   - 使用 InterProScan 或 eggNOG-mapper
   - 识别蛋白域和功能

## 常见问题

### Q: 为什么我的完整 ORF 比例很低？
A: 可能的原因：
1. 测序深度不足
2. RNA 样品质量差
3. Trinity 参数需要优化
4. 物种进化距离较远，基因较新颖

### Q: single_best_only 应该设为 true 还是 false？
A: 
- **false（推荐）**: 保留所有可能的 ORF，适合研究可变剪切
- **true**: 只保留最佳 ORF，减少冗余，适合快速筛选

### Q: 不完整的 ORF 是否就是错误的？
A: 不一定！可能原因：
1. 真实的组装不完整（需要 RACE）
2. 真实的可变剪切（N端或C端截断）
3. 测序或组装错误

需要结合其他证据（同源比对、表达量、长读长）综合判断。

## 参考资料

- TransDecoder 官方文档: https://github.com/TransDecoder/TransDecoder/wiki
- Trinity 转录组组装指南: https://github.com/trinityrnaseq/trinityrnaseq/wiki

