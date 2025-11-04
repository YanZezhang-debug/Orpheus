"""
Orpheus 主流程管道
"""

import os
from pathlib import Path
from typing import Optional
import logging

from .config import Config
from .utils import setup_logger, ensure_dir, validate_fasta
from .tools import CDHitRunner, TransDecoderRunner, BUSCORunner, TranscriptScorer


class OrpheusPipeline:
    """Orpheus 流程管道"""
    
    # 定义流程步骤（按执行顺序）
    PIPELINE_STEPS = [
        'busco_before',    # 处理前的 BUSCO 评估（可选）
        'cdhit',
        'transdecoder',
        'busco_after',     # 处理后的 BUSCO 评估（可选）
        'busco_scoring',   # BUSCO ORF 评分和转录本筛选
    ]
    
    def __init__(self, config_path: Optional[str] = None, threads: Optional[int] = None):
        """
        初始化流程管道
        
        Args:
            config_path: 配置文件路径
            threads: 线程数（如果指定，将覆盖配置文件中的线程设置）
        """
        # 加载配置
        self.config = Config(config_path)
        
        # 保存线程参数
        self.threads_override = threads
        
        # 设置日志
        log_level = self.config.get('logging.level', 'INFO')
        log_file = self.config.get('logging.log_file', 'orpheus.log')
        console = self.config.get('logging.console', True)
        
        self.logger = setup_logger('Orpheus', log_level, log_file, console)
        self.logger.info("=" * 60)
        self.logger.info("Orpheus (奥菲斯) - 转录本可信性评估工具")
        self.logger.info("版本: 0.2.0-dev (支持 TransDecoder)")
        self.logger.info("=" * 60)
        
        # 如果指定了线程数，显示覆盖信息
        if self.threads_override is not None:
            self.logger.info(f"✓ 使用命令行指定的线程数: {self.threads_override}")
            self.logger.info("  （将覆盖配置文件中的线程设置）")
        
        # 初始化工作目录
        self.work_dir = self.config.get('io.work_dir', './orpheus_output')
        ensure_dir(self.work_dir)
        self.logger.info(f"工作目录: {os.path.abspath(self.work_dir)}")
        
        # 用于存储中间文件路径
        self.cdhit_output_file = None
        self.transdecoder_output_dir = None
        self.busco_before_dir = None
        self.busco_after_dir = None
        self.scoring_output_file = None
    
    def run_cdhit(self, input_file: Optional[str] = None, 
                  output_file: Optional[str] = None) -> bool:
        """
        运行CD-HIT去冗余
        
        Args:
            input_file: 输入FASTA文件，如果为None则从配置读取
            output_file: 输出FASTA文件，如果为None则从配置读取
        
        Returns:
            True如果成功，False否则
        """
        # 获取输入输出文件
        if input_file is None:
            input_file = self.config.get('io.trinity_assembly')
        
        if not input_file:
            self.logger.error("未指定输入文件！请在配置文件中设置 io.trinity_assembly 或作为参数传入")
            return False
        
        if output_file is None:
            output_file = self.config.get('io.cdhit_output', 'cdhit_result.fasta')
        
        # 如果输出文件是相对路径，放到工作目录下
        if not os.path.isabs(output_file):
            output_file = os.path.join(self.work_dir, output_file)
        
        # 验证输入文件
        if not validate_fasta(input_file):
            self.logger.error(f"输入文件不是有效的FASTA文件: {input_file}")
            return False
        
        # 创建CD-HIT运行器
        cdhit_config = self.config.get('cdhit', {})
        
        # 如果命令行指定了线程数，覆盖配置文件中的值
        if self.threads_override is not None:
            cdhit_config = cdhit_config.copy()  # 避免修改原配置
            cdhit_config['threads'] = self.threads_override
            self.logger.info(f"  CD-HIT 线程数: {self.threads_override} (命令行覆盖)")
        
        runner = CDHitRunner(cdhit_config, self.logger)
        
        # 检查CD-HIT安装
        if not runner.check_installation():
            self.logger.error("CD-HIT-EST 未安装或不可用！")
            self.logger.error("请安装 CD-HIT: https://github.com/weizhongli/cdhit")
            return False
        
        # 运行CD-HIT
        success = runner.run(input_file, output_file)
        
        if success:
            self.logger.info("CD-HIT 去冗余步骤完成！")
            # 保存输出文件路径供后续步骤使用
            self.cdhit_output_file = output_file
        else:
            self.logger.error("CD-HIT 去冗余步骤失败！")
        
        return success
    
    def run_transdecoder(self, input_file: Optional[str] = None,
                        output_dir: Optional[str] = None) -> bool:
        """
        运行 TransDecoder ORF 预测
        
        Args:
            input_file: 输入转录本 FASTA 文件，如果为 None 则自动查找
            output_dir: 输出目录，如果为 None 则从配置读取
        
        Returns:
            True 如果成功，False 否则
        """
        # 获取输入文件（优先级顺序）
        if input_file is None:
            # 1. 首先尝试使用上一步 CD-HIT 的输出
            if self.cdhit_output_file is not None:
                input_file = self.cdhit_output_file
                self.logger.info(f"使用 CD-HIT 输出文件: {input_file}")
            else:
                # 2. 尝试在工作目录中查找 cdhit_result.fasta
                cdhit_result = os.path.join(self.work_dir, 'cdhit_result.fasta')
                if os.path.exists(cdhit_result):
                    input_file = cdhit_result
                    self.logger.info(f"找到已存在的 CD-HIT 结果文件: {input_file}")
                else:
                    # 3. 最后使用配置文件中的原始输入文件
                    input_file = self.config.get('io.trinity_assembly')
                    if input_file:
                        self.logger.info(f"使用配置文件中的输入文件: {input_file}")
        
        if not input_file or not os.path.exists(input_file):
            self.logger.error(f"输入文件不存在: {input_file}")
            self.logger.error("请确保以下之一存在:")
            self.logger.error(f"  1. 工作目录中的 CD-HIT 结果: {os.path.join(self.work_dir, 'cdhit_result.fasta')}")
            self.logger.error("  2. 配置文件中指定的输入文件")
            self.logger.error("  3. 通过命令行 -i 参数指定的输入文件")
            return False
        
        # 获取输出目录
        if output_dir is None:
            output_dir = self.config.get('io.transdecoder_output', 'transdecoder_results')
        
        # 如果输出目录是相对路径，放到工作目录下
        if not os.path.isabs(output_dir):
            output_dir = os.path.join(self.work_dir, output_dir)
        
        # 创建 TransDecoder 运行器
        transdecoder_config = self.config.get('transdecoder', {})
        
        # 如果命令行指定了线程数，覆盖同源搜索的线程数
        if self.threads_override is not None:
            transdecoder_config = transdecoder_config.copy()  # 避免修改原配置
            # 深拷贝 homology_search 配置
            if 'homology_search' in transdecoder_config:
                transdecoder_config['homology_search'] = transdecoder_config['homology_search'].copy()
                transdecoder_config['homology_search']['threads'] = self.threads_override
                self.logger.info(f"  同源搜索线程数: {self.threads_override} (命令行覆盖)")
        
        runner = TransDecoderRunner(transdecoder_config, self.logger)
        
        # 检查 TransDecoder 安装
        if not runner.check_installation():
            self.logger.error("TransDecoder 未安装或不可用！")
            self.logger.error("请安装 TransDecoder:")
            self.logger.error("  conda install -c bioconda transdecoder")
            self.logger.error("  或访问: https://github.com/TransDecoder/TransDecoder")
            return False
        
        # 运行 TransDecoder
        success = runner.run(input_file, output_dir)
        
        if success:
            self.logger.info("TransDecoder ORF 预测步骤完成！")
            self.transdecoder_output_dir = output_dir
        else:
            self.logger.error("TransDecoder ORF 预测步骤失败！")
        
        return success
    
    def run_busco(self, input_file: str, stage: str = 'before') -> bool:
        """
        运行 BUSCO 质量评估
        
        Args:
            input_file: 输入序列文件
            stage: 阶段标识 ('before' 或 'after')
        
        Returns:
            True 如果成功或跳过，False 如果失败
        """
        # 创建 BUSCO 运行器
        busco_config = self.config.get('busco', {})
        
        # 如果未启用，直接返回成功
        if not busco_config.get('enabled', False):
            self.logger.info(f"BUSCO 评估（{stage}）未启用，跳过")
            return True
        
        # 如果命令行指定了线程数，覆盖配置文件中的值
        if self.threads_override is not None:
            busco_config = busco_config.copy()  # 避免修改原配置
            busco_config['threads'] = self.threads_override
            self.logger.info(f"  BUSCO 线程数: {self.threads_override} (命令行覆盖)")
        
        runner = BUSCORunner(busco_config, self.logger)
        
        # 检查 BUSCO 安装
        if not runner.check_installation():
            self.logger.warning("BUSCO 未安装或不可用，跳过质量评估")
            self.logger.info("提示: 安装 BUSCO 可以评估转录组完整性")
            self.logger.info("  conda install -c bioconda busco")
            return True  # 不作为失败，仅跳过
        
        # 设置输出目录
        output_dir = os.path.join(self.work_dir, f'busco_{stage}')
        
        # 运行 BUSCO
        run_name = f"orpheus_{stage}"
        success = runner.run(input_file, output_dir, run_name)
        
        if success:
            if stage == 'before':
                self.busco_before_dir = output_dir
            else:
                self.busco_after_dir = output_dir
            
            # 如果两个阶段都完成了，比较结果
            if self.busco_before_dir and self.busco_after_dir and stage == 'after':
                before_summary = runner._find_summary_file(self.busco_before_dir, "orpheus_before")
                after_summary = runner._find_summary_file(self.busco_after_dir, "orpheus_after")
                if before_summary and after_summary:
                    runner.compare_results(before_summary, after_summary)
        
        return success
    
    def run_busco_scoring(self, gff_file: Optional[str] = None, 
                          cdhit_file: Optional[str] = None,
                          transdecoder_dir: Optional[str] = None,
                          output_file: Optional[str] = None) -> bool:
        """
        运行 BUSCO ORF 评分和转录本筛选
        
        Args:
            gff_file: TransDecoder 的 GFF3 文件，如果为 None 则自动查找
            cdhit_file: CD-HIT 的输出文件，如果为 None 则自动查找
            transdecoder_dir: TransDecoder 目录，如果为 None 则自动查找
            output_file: 输出文件路径，如果为 None 则从配置读取
        
        Returns:
            True 如果成功，False 否则
        """
        # 获取 GFF 文件（优先级顺序）
        if gff_file is None:
            # 1. 首先尝试使用 self.transdecoder_output_dir
            if self.transdecoder_output_dir:
                td_dir = Path(self.transdecoder_output_dir)
                gff_files = list(td_dir.glob("*.transdecoder.gff3"))
                if gff_files:
                    gff_file = str(gff_files[0])
                    self.logger.info(f"使用 TransDecoder 输出的 GFF3 文件: {gff_file}")
            
            # 2. 如果没有，尝试在工作目录的 transdecoder_results 中查找
            if gff_file is None:
                td_dir = Path(self.work_dir) / 'transdecoder_results'
                if td_dir.exists():
                    # 首先尝试查找最终的 .transdecoder.gff3 文件
                    gff_files = list(td_dir.glob("*.transdecoder.gff3"))
                    if gff_files:
                        gff_file = str(gff_files[0])
                        self.transdecoder_output_dir = str(td_dir)
                        self.logger.info(f"✓ 自动找到 GFF3 文件: {gff_file}")
                        
                        # 同时设置 transdecoder_dir（如果未指定）
                        if transdecoder_dir is None:
                            # 查找 transdecoder_dir 目录
                            td_dir_path = td_dir / (gff_files[0].stem.replace('.transdecoder', '') + '.transdecoder_dir')
                            if td_dir_path.exists():
                                transdecoder_dir = str(td_dir_path)
                                self.logger.info(f"✓ 自动找到 TransDecoder 目录: {transdecoder_dir}")
                    else:
                        # 如果没有找到最终文件，尝试查找 .transdecoder_dir 中的 longest_orfs.gff3
                        transdecoder_dirs = list(td_dir.glob("*.transdecoder_dir"))
                        for td_dir_candidate in transdecoder_dirs:
                            longest_orfs_gff = td_dir_candidate / "longest_orfs.gff3"
                            if longest_orfs_gff.exists():
                                gff_file = str(longest_orfs_gff)
                                self.transdecoder_output_dir = str(td_dir)
                                self.logger.warning(f"⚠ 找到中间 GFF3 文件: {gff_file}")
                                self.logger.warning("  建议完整运行 TransDecoder.Predict 步骤以生成最终的 .transdecoder.gff3 文件")
                                
                                # 设置 transdecoder_dir
                                if transdecoder_dir is None:
                                    transdecoder_dir = str(td_dir_candidate)
                                    self.logger.info(f"✓ 自动找到 TransDecoder 目录: {transdecoder_dir}")
                                break
            
            if gff_file is None:
                self.logger.error("未找到 TransDecoder GFF3 文件")
                self.logger.error("请确保以下之一存在:")
                self.logger.error(f"  1. 工作目录中的 TransDecoder 结果: {Path(self.work_dir) / 'transdecoder_results'}")
                self.logger.error("  2. 已运行 TransDecoder 步骤")
                self.logger.error("  3. 手动指定 GFF3 文件路径")
                return False
        
        if not os.path.exists(gff_file):
            self.logger.error(f"GFF3 文件不存在: {gff_file}")
            return False
        
        # 获取 CD-HIT 文件（优先级顺序）
        if cdhit_file is None:
            # 1. 首先尝试使用上一步 CD-HIT 的输出
            if self.cdhit_output_file is not None:
                cdhit_file = self.cdhit_output_file
                self.logger.info(f"使用 CD-HIT 输出文件: {cdhit_file}")
            else:
                # 2. 尝试在工作目录中查找 cdhit_result.fasta
                cdhit_result = os.path.join(self.work_dir, 'cdhit_result.fasta')
                if os.path.exists(cdhit_result):
                    cdhit_file = cdhit_result
                    self.cdhit_output_file = cdhit_file
                    self.logger.info(f"✓ 自动找到 CD-HIT 结果文件: {cdhit_file}")
        
        if not cdhit_file or not os.path.exists(cdhit_file):
            self.logger.error(f"CD-HIT 文件不存在: {cdhit_file or '(未指定)'}")
            self.logger.error("请确保以下之一存在:")
            self.logger.error(f"  1. 工作目录中的 CD-HIT 结果: {os.path.join(self.work_dir, 'cdhit_result.fasta')}")
            self.logger.error("  2. 已运行 CD-HIT 步骤")
            self.logger.error("  3. 手动指定 CD-HIT 文件路径")
            return False
        
        # 获取输出文件
        if output_file is None:
            output_file = self.config.get('io.scoring_output', 'scored_transcripts.tsv')
        
        # 如果输出文件是相对路径，放到工作目录下
        if not os.path.isabs(output_file):
            output_file = os.path.join(self.work_dir, output_file)
        
        # 创建评分器
        scorer = TranscriptScorer(self.logger)
        
        # 获取 BUSCO 输出目录（如果有）
        busco_dir = None
        # 首先检查正确的属性名
        if hasattr(self, 'busco_after_dir') and self.busco_after_dir:
            busco_dir = self.busco_after_dir
            self.logger.info(f"使用 BUSCO after 结果目录: {busco_dir}")
        elif hasattr(self, 'busco_before_dir') and self.busco_before_dir:
            busco_dir = self.busco_before_dir
            self.logger.info(f"使用 BUSCO before 结果目录: {busco_dir}")
        else:
            # 尝试在工作目录中自动查找 BUSCO 结果
            busco_after = os.path.join(self.work_dir, 'busco_after')
            busco_before = os.path.join(self.work_dir, 'busco_before')
            
            if os.path.exists(busco_after):
                busco_dir = busco_after
                self.busco_after_dir = busco_dir
                self.logger.info(f"✓ 自动找到 BUSCO after 结果目录: {busco_dir}")
            elif os.path.exists(busco_before):
                busco_dir = busco_before
                self.busco_before_dir = busco_dir
                self.logger.info(f"✓ 自动找到 BUSCO before 结果目录: {busco_dir}")
        
        # 获取评分配置
        scoring_config = self.config.get('scoring', {})
        threshold = scoring_config.get('threshold', 0.5)
        top_n = scoring_config.get('top_n', None)
        weights = scoring_config.get('weights', None)
        
        # 获取同源搜索结果文件（如果有）
        homology_file = scoring_config.get('homology_file', None)
        
        # 如果配置文件中未指定同源搜索文件，尝试自动查找
        if homology_file is None:
            # 尝试从 transdecoder_dir 或 gff_file 路径推断
            search_dir = None
            
            if transdecoder_dir is not None:
                search_dir = Path(transdecoder_dir)
            elif gff_file is not None:
                # 从 GFF3 文件路径推断 transdecoder_dir
                # 例如: /path/to/input.fasta.transdecoder.gff3 -> /path/to/input.fasta.transdecoder_dir
                gff_path = Path(gff_file)
                if '.transdecoder.gff3' in gff_path.name:
                    base_name = gff_path.name.replace('.transdecoder.gff3', '')
                    search_dir = gff_path.parent / f"{base_name}.transdecoder_dir"
            
            if search_dir is not None:
                blastp_file = search_dir / "blastp_results.outfmt6"
                if blastp_file.exists():
                    homology_file = str(blastp_file)
                    self.logger.info(f"✓ 自动找到同源搜索结果文件: {homology_file}")
                else:
                    self.logger.debug(f"未找到同源搜索结果文件: {blastp_file}")
        
        # 运行评分
        success = scorer.score_and_filter(
            gff3_file=gff_file,
            cdhit_file=cdhit_file,
            output_dir=os.path.dirname(output_file),
            homology_file=homology_file,
            busco_dir=busco_dir,
            threshold=threshold,
            top_n=top_n,
            weights=weights
        )
        
        if success:
            self.logger.info("BUSCO ORF 评分和转录本筛选完成！")
            self.scoring_output_file = output_file
        else:
            self.logger.error("BUSCO ORF 评分和转录本筛选失败！")
        
        return success
    
    def get_steps_to_run(self, start_from: str = 'beginning') -> list:
        """
        根据 start_from 参数获取需要执行的步骤列表
        
        Args:
            start_from: 从哪一步开始，可以是 'beginning' 或步骤名称
        
        Returns:
            需要执行的步骤列表
        """
        if start_from == 'beginning':
            return self.PIPELINE_STEPS.copy()
        
        if start_from not in self.PIPELINE_STEPS:
            raise ValueError(f"未知的步骤: {start_from}. 可用步骤: {', '.join(self.PIPELINE_STEPS)}")
        
        # 从指定步骤开始到结束的所有步骤
        start_index = self.PIPELINE_STEPS.index(start_from)
        return self.PIPELINE_STEPS[start_index:]
    
    def run(self, input_file: Optional[str] = None, start_from: str = 'beginning'):
        """
        运行完整流程
        
        Args:
            input_file: Trinity组装结果文件（或中间文件，取决于start_from）
            start_from: 从哪一步开始执行 ('beginning' 或步骤名称)
        """
        self.logger.info("开始运行 Orpheus 流程...")
        
        # 获取需要执行的步骤
        steps_to_run = self.get_steps_to_run(start_from)
        total_steps = len(steps_to_run)
        
        if start_from != 'beginning':
            self.logger.info(f"从步骤 '{start_from}' 开始执行")
            skipped_steps = [s for s in self.PIPELINE_STEPS if s not in steps_to_run]
            if skipped_steps:
                self.logger.info(f"跳过步骤: {', '.join(skipped_steps)}")
        
        # 步骤执行结果
        step_results = {}
        
        # 当前输入文件（会在步骤间传递）
        current_input = input_file
        
        # 执行每个步骤
        for idx, step in enumerate(steps_to_run, 1):
            self.logger.info(f"\n步骤 {idx}/{total_steps}: {self._get_step_name(step)}")
            
            # 根据步骤名称调用相应的方法
            if step == 'busco_before':
                # 处理前的 BUSCO 评估（使用原始输入文件）
                success = self.run_busco(current_input, stage='before')
            elif step == 'cdhit':
                success = self.run_cdhit(current_input)
                # CD-HIT 完成后，更新当前输入为 CD-HIT 输出
                if success:
                    current_input = self.cdhit_output_file
            elif step == 'transdecoder':
                # 如果是直接从 transdecoder 开始，使用命令行传入的 input_file
                # 否则使用上一步（CD-HIT）的输出
                success = self.run_transdecoder(current_input)
            elif step == 'busco_after':
                # 处理后的 BUSCO 评估（使用当前输入文件，通常是 CD-HIT 的输出）
                success = self.run_busco(current_input, stage='after')
            elif step == 'busco_scoring':
                # BUSCO ORF 评分和转录本筛选
                # 这一步需要 GFF 文件和 CD-HIT 文件
                success = self.run_busco_scoring()
            else:
                # 未来添加的步骤可以在这里处理
                self.logger.warning(f"步骤 '{step}' 尚未实现，跳过")
                success = True
            
            step_results[step] = success
            
            # 如果步骤失败，停止执行
            if not success:
                self.logger.error("\n" + "=" * 60)
                self.logger.error(f"Orpheus 流程失败！（{self._get_step_name(step)} 步骤失败）")
                self.logger.error("=" * 60)
                self._print_step_results(step_results)
                return False
        
        # 所有步骤完成
        all_success = all(step_results.values())
        
        if all_success:
            self.logger.info("\n" + "=" * 60)
            self.logger.info("Orpheus 流程完成！")
            self.logger.info("=" * 60)
            self._print_summary()
        else:
            self.logger.error("\n" + "=" * 60)
            self.logger.error("Orpheus 流程部分失败！")
            self.logger.error("=" * 60)
            self._print_step_results(step_results)
        
        return all_success
    
    def _get_step_name(self, step: str) -> str:
        """
        获取步骤的友好名称
        
        Args:
            step: 步骤标识符
        
        Returns:
            步骤的友好名称
        """
        step_names = {
            'busco_before': 'BUSCO 质量评估（处理前）',
            'cdhit': 'CD-HIT 去冗余',
            'transdecoder': 'TransDecoder ORF 预测',
            'busco_after': 'BUSCO 质量评估（处理后）',
            'busco_scoring': 'BUSCO ORF 评分和转录本筛选',
        }
        return step_names.get(step, step)
    
    def _print_step_results(self, step_results: dict):
        """
        打印步骤执行结果
        
        Args:
            step_results: 步骤执行结果字典
        """
        self.logger.info("\n步骤执行结果:")
        for step, success in step_results.items():
            status = "✓ 成功" if success else "✗ 失败"
            self.logger.info(f"  {self._get_step_name(step)}: {status}")
    
    def _print_summary(self):
        """打印结果摘要"""
        self.logger.info("\n结果摘要:")
        
        if self.cdhit_output_file:
            self.logger.info(f"  CD-HIT 结果: {self.cdhit_output_file}")
        
        if self.transdecoder_output_dir:
            self.logger.info(f"  TransDecoder 结果目录: {self.transdecoder_output_dir}")
            
            # 列出主要输出文件
            td_dir = Path(self.transdecoder_output_dir)
            pep_file = None
            cds_file = None
            gff3_file = None
            
            for f in td_dir.glob("*.transdecoder.pep"):
                pep_file = f
                break
            for f in td_dir.glob("*.transdecoder.cds"):
                cds_file = f
                break
            for f in td_dir.glob("*.transdecoder.gff3"):
                gff3_file = f
                break
            
            if pep_file:
                self.logger.info(f"    - 预测蛋白序列: {pep_file.name}")
            if cds_file:
                self.logger.info(f"    - 预测 CDS 序列: {cds_file.name}")
            if gff3_file:
                self.logger.info(f"    - GFF3 注释: {gff3_file.name}")
        
        if self.scoring_output_file:
            self.logger.info(f"  评分和筛选结果: {self.scoring_output_file}")
        
        self.logger.info("\n建议:")
        if self.scoring_output_file:
            self.logger.info("  1. 查看评分结果文件了解每个转录本的质量得分")
            self.logger.info("  2. 使用筛选后的高质量转录本进行下游分析")
        else:
            self.logger.info("  1. 查看 GFF3 文件了解 ORF 完整度信息")
            self.logger.info("  2. 使用预测的蛋白序列进行同源比对（DIAMOND/BLAST）")
            self.logger.info("  3. 运行 busco_scoring 步骤进行转录本质量评估")

