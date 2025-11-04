#!/usr/bin/env python3
"""
Orpheus (奥菲斯) 命令行入口
无参转录组组装转录本可信性评估工具
"""

import argparse
import sys
from pathlib import Path

from orpheus.pipeline import OrpheusPipeline
from orpheus import __version__


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='Orpheus - De novo Transcriptome Assembly Reliability Assessment Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run with default configuration from the beginning
  orpheus -i trinity_assembly.fasta
  
  # Use custom configuration file
  orpheus -i trinity_assembly.fasta -c my_config.yaml
  
  # Specify thread count (overrides configuration file defaults)
  orpheus -i trinity_assembly.fasta -t 32
  
  # Start from CD-HIT step (skip previous steps)
  orpheus -i trinity_assembly.fasta --start-from cdhit
  
  # Run only TransDecoder ORF prediction with 16 threads
  orpheus -i cdhit_result.fasta --start-from transdecoder -t 16
  
  # Run BUSCO scoring step (auto-finds required files in work directory)
  orpheus --start-from busco_scoring
  
  # Full example: specify config, threads, and output
  orpheus -i trinity_assembly.fasta -c my_config.yaml -t 24 -o result.fasta

For more information, visit: https://github.com/yourusername/orpheus
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=False,
        default=None,
        metavar='FILE',
        help='Input transcriptome assembly file in FASTA format (e.g., Trinity.fasta)\n'
             'Optional when using --start-from with auto-file finding'
    )
    
    parser.add_argument(
        '-c', '--config',
        default=None,
        metavar='FILE',
        help='Path to configuration file (default: config/default.yaml)'
    )
    
    parser.add_argument(
        '-o', '--output',
        default=None,
        metavar='FILE',
        help='Output file path (default: read from configuration file)'
    )
    
    parser.add_argument(
        '--start-from',
        choices=['beginning', 'busco_before', 'cdhit', 'transdecoder', 'busco_after', 'busco_scoring'],
        default='beginning',
        metavar='STEP',
        help='Starting step for pipeline execution (default: beginning)\n'
             '  beginning     - Execute all steps from the start\n'
             '  busco_before  - BUSCO quality assessment (before processing)\n'
             '  cdhit         - Start from CD-HIT redundancy removal\n'
             '  transdecoder  - Run TransDecoder ORF prediction\n'
             '  busco_after   - BUSCO quality assessment (after processing)\n'
             '  busco_scoring - BUSCO-based ORF scoring and transcript filtering'
    )
    
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=None,
        metavar='N',
        help='Number of threads to use (overrides configuration file settings)\n'
             'Applies to both CD-HIT clustering and homology searches (Diamond/BLASTP)'
    )
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f'Orpheus {__version__}'
    )
    
    args = parser.parse_args()
    
    # 检查输入文件是否存在（仅当提供了输入文件时）
    if args.input is not None and not Path(args.input).exists():
        print(f"错误: 输入文件不存在: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # 验证：除了 busco_scoring 步骤外，其他所有步骤都需要提供输入文件
    # busco_scoring 可以自动查找所需文件（GFF 和 CD-HIT 输出）
    if args.input is None and args.start_from != 'busco_scoring':
        print(f"错误: 步骤 '{args.start_from}' 需要提供输入文件 (-i 参数)", file=sys.stderr)
        print(f"提示: 只有 'busco_scoring' 步骤可以不提供输入文件（自动查找所需文件）", file=sys.stderr)
        sys.exit(1)
    
    try:
        # 初始化流程
        pipeline = OrpheusPipeline(config_path=args.config, threads=args.threads)
        
        # 运行流程
        success = pipeline.run(input_file=args.input, start_from=args.start_from)
        
        # 返回退出码
        sys.exit(0 if success else 1)
        
    except Exception as e:
        print(f"错误: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()

