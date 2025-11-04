# Orpheus 使用示例

## 1. 基本用法 - 完整流程

从头开始运行完整流程（CD-HIT + TransDecoder）：

```bash
python orpheus_cli.py -i trinity_assembly.fasta
```

这将：
1. 使用 CD-HIT 去除冗余转录本
2. 使用 TransDecoder 预测 ORF 并评估完整度

## 2. 使用自定义配置

```bash
python orpheus_cli.py -i trinity_assembly.fasta -c my_config.yaml
```

## 3. 跳过步骤执行

### 场景 1：仅运行 TransDecoder（跳过 CD-HIT）

当你已经有去冗余的序列，或者想测试不同的 TransDecoder 参数：

```bash
python orpheus_cli.py -i cdhit_result.fasta --start-from transdecoder
```

### 场景 2：从 CD-HIT 开始

如果有之前步骤的结果，想从 CD-HIT 重新开始：

```bash
python orpheus_cli.py -i trinity_assembly.fasta --start-from cdhit
```

### 场景 3：步骤失败后重新运行

假设 TransDecoder 步骤失败了，修复问题后从该步骤继续：

```bash
# 修复 TransDecoder 配置或安装问题后
python orpheus_cli.py -i cdhit_output.fasta --start-from transdecoder
```

## 4. 参数组合使用

自定义配置 + 跳过步骤：

```bash
python orpheus_cli.py -i cdhit_result.fasta \
    -c custom_transdecoder.yaml \
    --start-from transdecoder
```

## 5. 工作流示例

### 典型工作流

```bash
# 步骤 1: 运行完整流程
python orpheus_cli.py -i trinity.fasta

# 步骤 2: 如果需要调整 TransDecoder 参数，编辑配置文件后重新运行
# 编辑 config/default.yaml 中的 transdecoder 部分
python orpheus_cli.py -i orpheus_output/cdhit_result.fasta --start-from transdecoder
```

### 调试工作流

```bash
# 1. 先只测试 CD-HIT（修改配置，将后续步骤注释）
python orpheus_cli.py -i trinity.fasta

# 2. CD-HIT 成功后，测试 TransDecoder
python orpheus_cli.py -i orpheus_output/cdhit_result.fasta --start-from transdecoder
```

## 6. 可用步骤列表

当前版本支持的步骤：

| 步骤名称        | 说明                    | 依赖工具        |
|----------------|------------------------|----------------|
| `cdhit`        | CD-HIT 去冗余           | CD-HIT-EST     |
| `transdecoder` | TransDecoder ORF 预测   | TransDecoder   |

**未来计划的步骤**（示例）：
- `diamond` - DIAMOND 同源比对
- `blast` - BLAST 同源比对
- `evaluation` - 综合评估和打分

## 7. 常见问题

### Q: 如何知道当前版本有哪些步骤？

```bash
python orpheus_cli.py --help
```

查看 `--start-from` 参数的可选值。

### Q: --start-from 和配置文件中的步骤控制有什么区别？

- `--start-from`: 控制从哪个步骤**开始**执行到最后
- 配置文件: 控制每个步骤的具体参数

### Q: 可以只运行某一个步骤（不运行后续步骤）吗？

当前版本不支持。`--start-from transdecoder` 会运行 TransDecoder 及之后的所有步骤。

如果未来添加了更多步骤（如 diamond），但你只想运行到 transdecoder，可以临时修改配置文件禁用后续步骤。

## 8. 扩展性说明

Orpheus 的步骤控制系统设计为可扩展的。开发者可以轻松添加新步骤：

1. 在 `orpheus/pipeline.py` 的 `PIPELINE_STEPS` 列表中添加步骤名称
2. 在 `_get_step_name()` 中添加友好名称映射
3. 在 `run()` 方法中添加步骤的执行逻辑
4. 更新 CLI 的 `--start-from` choices

示例（未来添加 DIAMOND 步骤）：

```python
# orpheus/pipeline.py
PIPELINE_STEPS = [
    'cdhit',
    'transdecoder',
    'diamond',  # 新步骤
]
```

命令行自动支持：

```bash
python orpheus_cli.py -i input.fasta --start-from diamond
```


