"""
TransDecoder 调用模块
用于预测转录本的开放阅读框（ORF）并评估完整度
支持基于同源搜索（Diamond/BLASTP）的证据辅助预测
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple
import logging
import re


class TransDecoderRunner:
    """TransDecoder 运行器"""
    
    def __init__(self, config: Dict[str, Any], logger: Optional[logging.Logger] = None):
        """
        初始化 TransDecoder 运行器
        
        Args:
            config: TransDecoder 配置字典
            logger: 日志记录器
        """
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        self.longorfs_exec = config.get('longorfs_executable', 'TransDecoder.LongOrfs')
        self.predict_exec = config.get('predict_executable', 'TransDecoder.Predict')
    
    def check_installation(self) -> bool:
        """
        检查 TransDecoder 是否已安装
        
        Returns:
            True 如果已安装，False 否则
        """
        try:
            # 方法1: 使用 shutil.which 检查命令是否存在
            longorfs_path = shutil.which(self.longorfs_exec)
            predict_path = shutil.which(self.predict_exec)
            
            if longorfs_path and predict_path:
                self.logger.info(f"TransDecoder 已找到")
                self.logger.info(f"  - LongOrfs: {longorfs_path}")
                self.logger.info(f"  - Predict: {predict_path}")
                return True
            
            # 方法2: 尝试运行命令检查
            # 检查 TransDecoder.LongOrfs
            result1 = subprocess.run(
                [self.longorfs_exec, '--help'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # 检查 TransDecoder.Predict
            result2 = subprocess.run(
                [self.predict_exec, '--help'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # 检查是否有输出且不是"找不到命令"错误
            def check_command(result, exec_name):
                has_output = bool(result.stdout or result.stderr)
                error_output = result.stderr.lower() if result.stderr else ""
                is_not_found = any(phrase in error_output for phrase in [
                    'command not found',
                    'no such file',
                    'cannot find',
                    '不是内部或外部命令'
                ])
                
                # returncode 为 0 或 255 通常表示命令存在
                if result.returncode in [0, 255]:
                    return True
                # 即使返回码不是0/255，只要有输出且不是"找不到"错误，也认为已安装
                elif has_output and not is_not_found:
                    return True
                return False
            
            installed = check_command(result1, self.longorfs_exec) and \
                       check_command(result2, self.predict_exec)
            
            if installed:
                self.logger.info(f"TransDecoder 已找到")
                self.logger.info(f"  - LongOrfs: {self.longorfs_exec}")
                self.logger.info(f"  - Predict: {self.predict_exec}")
            else:
                self.logger.error(f"TransDecoder 未正确安装")
                if not check_command(result1, self.longorfs_exec):
                    self.logger.error(f"  - LongOrfs 未找到: {self.longorfs_exec}")
                if not check_command(result2, self.predict_exec):
                    self.logger.error(f"  - Predict 未找到: {self.predict_exec}")
            
            return installed
            
        except FileNotFoundError:
            self.logger.error(f"TransDecoder 未找到")
            self.logger.error("请确保 TransDecoder 已安装并在 PATH 中")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"检查 TransDecoder 时超时")
            return False
        except Exception as e:
            self.logger.error(f"检查 TransDecoder 安装时出错: {e}")
            return False
    
    def build_longorfs_command(self, input_file: str, output_dir: str) -> list:
        """
        构建 TransDecoder.LongOrfs 命令
        
        Args:
            input_file: 输入转录本 FASTA 文件
            output_dir: 输出目录
        
        Returns:
            命令列表
        """
        cmd = [
            self.longorfs_exec,
            '-t', input_file
        ]
        
        # 最小蛋白长度
        min_protein_len = self.config.get('min_protein_length', 100)
        if min_protein_len:
            cmd.extend(['-m', str(min_protein_len)])
        
        # 遗传密码表
        genetic_code = self.config.get('genetic_code', 'universal')
        if genetic_code and genetic_code != 'universal':
            cmd.extend(['-G', genetic_code])
        
        # 输出目录（通过工作目录设置）
        # TransDecoder 会在当前目录创建 .transdecoder_dir
        
        return cmd
    
    def build_predict_command(self, input_file: str, blastp_output: Optional[str] = None) -> list:
        """
        构建 TransDecoder.Predict 命令
        
        Args:
            input_file: 输入转录本 FASTA 文件
            blastp_output: BLASTP/Diamond 搜索结果文件路径（可选）
        
        Returns:
            命令列表
        """
        cmd = [
            self.predict_exec,
            '-t', input_file
        ]
        
        # 如果提供了同源搜索结果，使用它作为证据
        if blastp_output and os.path.exists(blastp_output):
            cmd.extend(['--retain_blastp_hits', blastp_output])
            self.logger.info(f"使用同源搜索证据: {blastp_output}")
        
        # 保留长 ORF（核苷酸长度，默认 900）
        # retain_long_orfs_length 配置值单位为核苷酸（nt）
        if self.config.get('retain_long_orfs', True):
            retain_nt_len = self.config.get('retain_long_orfs_length', 900)
            cmd.extend(['--retain_long_orfs', str(retain_nt_len)])
        
        # 单一最佳 ORF
        if self.config.get('single_best_only', False):
            cmd.append('--single_best_only')
        
        return cmd
    
    def run(self, input_file: str, output_dir: str) -> bool:
        """
        运行 TransDecoder 完整流程
        
        Args:
            input_file: 输入转录本 FASTA 文件
            output_dir: 输出目录
        
        Returns:
            True 如果成功，False 否则
        """
        # 检查输入文件
        if not os.path.exists(input_file):
            self.logger.error(f"输入文件不存在: {input_file}")
            return False
        
        # 确保输出目录存在
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # 获取输入文件的绝对路径
        input_file_abs = os.path.abspath(input_file)
        
        self.logger.info("=" * 60)
        self.logger.info("开始运行 TransDecoder")
        self.logger.info(f"输入文件: {input_file_abs}")
        self.logger.info(f"输出目录: {os.path.abspath(output_dir)}")
        self.logger.info("=" * 60)
        
        # 清理旧的 transdecoder_dir（避免残留文件导致错误）
        transdecoder_dir = Path(output_dir) / f"{Path(input_file_abs).name}.transdecoder_dir"
        if transdecoder_dir.exists():
            self.logger.info(f"检测到旧的 TransDecoder 工作目录，正在清理: {transdecoder_dir}")
            try:
                shutil.rmtree(transdecoder_dir)
                self.logger.info("✓ 旧工作目录已清理")
            except Exception as e:
                self.logger.warning(f"清理旧工作目录失败: {e}")
                self.logger.warning("这可能导致 TransDecoder 运行出错")
        
        # 步骤 1: TransDecoder.LongOrfs
        self.logger.info("\n步骤 1/3: 识别长 ORF (TransDecoder.LongOrfs)")
        if not self._run_longorfs(input_file_abs, output_dir):
            return False
        
        # 步骤 2: 同源搜索（可选）
        blastp_output = None
        homology_config = self.config.get('homology_search', {})
        if homology_config.get('enabled', False):
            self.logger.info("\n步骤 2/3: 同源搜索 (Diamond/BLASTP)")
            
            # 找到 longest_orfs.pep 文件
            transdecoder_dir = Path(output_dir) / f"{Path(input_file_abs).name}.transdecoder_dir"
            pep_file = transdecoder_dir / "longest_orfs.pep"
            
            if not pep_file.exists():
                self.logger.error(f"未找到 longest_orfs.pep 文件: {pep_file}")
                return False
            
            # 设置同源搜索输出文件
            blastp_output = str(transdecoder_dir / "blastp_results.outfmt6")
            
            # 获取是否只使用完整 ORF 的配置
            use_complete_only = homology_config.get('use_complete_only', True)
            
            if not self.run_homology_search(str(pep_file), blastp_output, use_complete_only):
                self.logger.warning("同源搜索失败，将继续使用长度标准")
                blastp_output = None
        else:
            self.logger.info("\n步骤 2/3: 同源搜索 (跳过)")
        
        # 步骤 3: TransDecoder.Predict
        self.logger.info("\n步骤 3/3: 预测编码区域 (TransDecoder.Predict)")
        if not self._run_predict(input_file_abs, output_dir, blastp_output):
            return False
        
        self.logger.info("\n" + "=" * 60)
        self.logger.info("TransDecoder 成功完成！")
        self.logger.info("=" * 60)
        
        # 分析结果
        self._analyze_results(input_file_abs, output_dir)
        
        return True
    
    def filter_complete_orfs(self, input_pep: str, output_pep: str) -> bool:
        """
        从 longest_orfs.pep 中筛选出完整的 ORF（带有起始密码子和终止密码子）
        
        TransDecoder 的 ORF 头部格式示例：
        >TRINITY_DN1000_c0_g1_i1.p1 TRINITY_DN1000_c0_g1_i1 type:complete len:300 (+)
        >TRINITY_DN1001_c0_g1_i1.p1 TRINITY_DN1001_c0_g1_i1 type:5prime_partial len:250 (+)
        
        Args:
            input_pep: 输入的 longest_orfs.pep 文件
            output_pep: 输出的完整 ORF 文件
        
        Returns:
            True 如果成功，False 否则
        """
        try:
            complete_count = 0
            total_count = 0
            
            with open(input_pep, 'r') as infile, open(output_pep, 'w') as outfile:
                current_header = None
                current_seq = []
                
                for line in infile:
                    if line.startswith('>'):
                        # 处理上一个序列
                        if current_header is not None:
                            total_count += 1
                            # 检查是否为完整 ORF
                            if 'type:complete' in current_header:
                                outfile.write(current_header + '\n')
                                outfile.write(''.join(current_seq))
                                complete_count += 1
                        
                        # 开始新序列
                        current_header = line.strip()
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                # 处理最后一个序列
                if current_header is not None:
                    total_count += 1
                    if 'type:complete' in current_header:
                        outfile.write(current_header + '\n')
                        outfile.write(''.join(current_seq))
                        complete_count += 1
            
            self.logger.info(f"✓ 完整 ORF 筛选完成:")
            self.logger.info(f"  总 ORF 数: {total_count}")
            self.logger.info(f"  完整 ORF: {complete_count} ({complete_count/total_count*100:.2f}%)" if total_count > 0 else "  完整 ORF: 0")
            self.logger.info(f"  不完整 ORF (已排除): {total_count - complete_count}")
            
            if complete_count == 0:
                self.logger.warning("⚠ 未找到任何完整的 ORF！")
                self.logger.warning("  这可能表明转录本质量较差或组装不完整")
                return False
            
            return True
            
        except Exception as e:
            self.logger.error(f"筛选完整 ORF 时出错: {e}")
            return False
    
    def run_homology_search(self, pep_file: str, output_file: str, use_complete_only: bool = True) -> bool:
        """
        运行同源搜索（Diamond 或 BLASTP）
        
        Args:
            pep_file: 输入蛋白序列文件（longest_orfs.pep）
            output_file: 输出文件路径
            use_complete_only: 是否只使用完整的 ORF 进行同源搜索（推荐）
        
        Returns:
            True 如果成功，False 否则
        """
        homology_config = self.config.get('homology_search', {})
        
        # 检查是否启用
        if not homology_config.get('enabled', False):
            self.logger.info("同源搜索未启用，跳过")
            return True
        
        # 检查数据库
        database = homology_config.get('database', '')
        if not database:
            self.logger.warning("未配置同源搜索数据库，跳过")
            return True
        
        if not os.path.exists(database):
            self.logger.error(f"数据库文件不存在: {database}")
            return False
        
        # 检查输入文件
        if not os.path.exists(pep_file):
            self.logger.error(f"蛋白序列文件不存在: {pep_file}")
            return False
        
        # 如果配置了只使用完整 ORF，先筛选
        query_file = pep_file
        if use_complete_only:
            self.logger.info("筛选完整的 ORF 用于同源搜索（推荐用于后续 PCR 实验）...")
            complete_pep = pep_file.replace('.pep', '.complete.pep')
            
            if self.filter_complete_orfs(pep_file, complete_pep):
                query_file = complete_pep
                self.logger.info(f"✓ 将使用完整 ORF 文件进行同源搜索: {complete_pep}")
            else:
                self.logger.warning("⚠ 完整 ORF 筛选失败，将使用所有 ORF")
                query_file = pep_file
        else:
            self.logger.info("使用所有 ORF（包括不完整的）进行同源搜索")
        
        # 选择工具
        tool = homology_config.get('tool', 'diamond').lower()
        
        if tool == 'diamond':
            return self._run_diamond_blastp(query_file, output_file, homology_config)
        elif tool == 'blastp':
            return self._run_blastp(query_file, output_file, homology_config)
        else:
            self.logger.error(f"不支持的同源搜索工具: {tool}")
            return False
    
    def _run_diamond_blastp(self, query_file: str, output_file: str, config: Dict[str, Any]) -> bool:
        """
        运行 Diamond BLASTP
        
        Args:
            query_file: 查询序列文件
            output_file: 输出文件
            config: 同源搜索配置
        
        Returns:
            True 如果成功，False 否则
        """
        diamond_exec = config.get('diamond_executable', 'diamond')
        database = config['database']
        evalue = config.get('evalue', 1e-5)
        max_target_seqs = config.get('max_target_seqs', 1)
        threads = config.get('threads', 8)
        extra_params = config.get('extra_params', '')
        
        self.logger.info("运行 Diamond BLASTP 同源搜索...")
        self.logger.info(f"  查询: {query_file}")
        self.logger.info(f"  数据库: {database}")
        self.logger.info(f"  E-value: {evalue}")
        
        cmd = [
            diamond_exec, 'blastp',
            '--query', query_file,
            '--db', database,
            '--out', output_file,
            '--outfmt', '6',  # 表格格式
            '--evalue', str(evalue),
            '--max-target-seqs', str(max_target_seqs),
            '--threads', str(threads)
        ]
        
        if extra_params:
            cmd.extend(extra_params.split())
        
        self.logger.info(f"命令: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout:
                self.logger.debug(f"Diamond 输出: {result.stdout[:500]}")
            
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                self.logger.info(f"✓ Diamond 搜索完成")
                self.logger.info(f"  输出文件: {output_file} ({file_size} bytes)")
                
                # 统计匹配数
                with open(output_file, 'r') as f:
                    hit_count = sum(1 for _ in f)
                self.logger.info(f"  找到 {hit_count} 个同源匹配")
                return True
            else:
                self.logger.error("Diamond 未生成输出文件")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Diamond 运行失败: {e}")
            if e.stdout:
                self.logger.error(f"标准输出: {e.stdout}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"Diamond 运行出错: {e}")
            return False
    
    def _run_blastp(self, query_file: str, output_file: str, config: Dict[str, Any]) -> bool:
        """
        运行 NCBI BLASTP
        
        Args:
            query_file: 查询序列文件
            output_file: 输出文件
            config: 同源搜索配置
        
        Returns:
            True 如果成功，False 否则
        """
        blastp_exec = config.get('blastp_executable', 'blastp')
        database = config['database']
        evalue = config.get('evalue', 1e-5)
        max_target_seqs = config.get('max_target_seqs', 1)
        threads = config.get('threads', 8)
        extra_params = config.get('extra_params', '')
        
        self.logger.info("运行 BLASTP 同源搜索...")
        self.logger.info(f"  查询: {query_file}")
        self.logger.info(f"  数据库: {database}")
        self.logger.info(f"  E-value: {evalue}")
        
        cmd = [
            blastp_exec,
            '-query', query_file,
            '-db', database,
            '-out', output_file,
            '-outfmt', '6',  # 表格格式
            '-evalue', str(evalue),
            '-max_target_seqs', str(max_target_seqs),
            '-num_threads', str(threads)
        ]
        
        if extra_params:
            cmd.extend(extra_params.split())
        
        self.logger.info(f"命令: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout:
                self.logger.debug(f"BLASTP 输出: {result.stdout[:500]}")
            
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                self.logger.info(f"✓ BLASTP 搜索完成")
                self.logger.info(f"  输出文件: {output_file} ({file_size} bytes)")
                
                # 统计匹配数
                with open(output_file, 'r') as f:
                    hit_count = sum(1 for _ in f)
                self.logger.info(f"  找到 {hit_count} 个同源匹配")
                return True
            else:
                self.logger.error("BLASTP 未生成输出文件")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"BLASTP 运行失败: {e}")
            if e.stdout:
                self.logger.error(f"标准输出: {e.stdout}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"BLASTP 运行出错: {e}")
            return False
    
    def _run_longorfs(self, input_file: str, output_dir: str) -> bool:
        """运行 TransDecoder.LongOrfs"""
        cmd = self.build_longorfs_command(input_file, output_dir)
        
        self.logger.info(f"命令: {' '.join(cmd)}")
        
        try:
            # TransDecoder 需要在输出目录运行
            result = subprocess.run(
                cmd,
                cwd=output_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout:
                self.logger.info("TransDecoder.LongOrfs 输出:")
                for line in result.stdout.splitlines():
                    if line.strip():
                        self.logger.info(f"  {line}")
            
            if result.stderr:
                self.logger.debug("TransDecoder.LongOrfs 调试信息:")
                for line in result.stderr.splitlines():
                    if line.strip():
                        self.logger.debug(f"  {line}")
            
            # 检查输出目录是否生成
            transdecoder_dir = Path(output_dir) / f"{Path(input_file).name}.transdecoder_dir"
            if transdecoder_dir.exists():
                self.logger.info(f"✓ 长 ORF 识别完成")
                self.logger.info(f"  输出目录: {transdecoder_dir}")
                return True
            else:
                self.logger.error("TransDecoder.LongOrfs 未生成预期输出目录")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"TransDecoder.LongOrfs 运行失败: {e}")
            if e.stdout:
                self.logger.error(f"标准输出: {e.stdout}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"运行 TransDecoder.LongOrfs 时发生错误: {e}")
            return False
    
    def _run_predict(self, input_file: str, output_dir: str, blastp_output: Optional[str] = None) -> bool:
        """
        运行 TransDecoder.Predict
        
        Args:
            input_file: 输入转录本文件
            output_dir: 输出目录
            blastp_output: 同源搜索结果文件（可选，绝对路径或相对于当前目录的路径）
        
        Returns:
            True 如果成功，False 否则
        """
        # 详细日志和路径处理
        if blastp_output:
            # 检查文件是否存在
            if os.path.exists(blastp_output):
                file_size = os.path.getsize(blastp_output)
                self.logger.info(f"使用同源搜索证据: {blastp_output}")
                
                # 计算相对于 output_dir 的相对路径（因为 TransDecoder 会在 output_dir 中运行）
                try:
                    blastp_abs = Path(blastp_output).resolve()
                    output_abs = Path(output_dir).resolve()
                    blastp_relative = os.path.relpath(blastp_abs, output_abs)
                    self.logger.debug(f"  绝对路径: {blastp_abs}")
                    self.logger.debug(f"  相对于运行目录的路径: {blastp_relative}")
                    blastp_output = blastp_relative
                except ValueError:
                    # Windows 上可能因为不同盘符导致无法计算相对路径，使用绝对路径
                    self.logger.debug(f"  无法计算相对路径，使用绝对路径")
                    blastp_output = str(Path(blastp_output).resolve())
            else:
                self.logger.warning(f"同源搜索结果文件不存在: {blastp_output}")
                self.logger.warning("将不使用同源搜索证据")
                blastp_output = None
        else:
            self.logger.info("未提供同源搜索结果，仅使用长度标准")
        
        cmd = self.build_predict_command(input_file, blastp_output)
        
        self.logger.info(f"命令: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                cwd=output_dir,
                capture_output=True,
                text=True,
                check=True
            )
            
            if result.stdout:
                self.logger.info("TransDecoder.Predict 输出:")
                for line in result.stdout.splitlines():
                    if line.strip():
                        self.logger.info(f"  {line}")
            
            if result.stderr:
                self.logger.debug("TransDecoder.Predict 调试信息:")
                for line in result.stderr.splitlines():
                    if line.strip():
                        self.logger.debug(f"  {line}")
            
            # 检查输出文件
            input_basename = Path(input_file).name
            pep_file = Path(output_dir) / f"{input_basename}.transdecoder.pep"
            cds_file = Path(output_dir) / f"{input_basename}.transdecoder.cds"
            gff3_file = Path(output_dir) / f"{input_basename}.transdecoder.gff3"
            bed_file = Path(output_dir) / f"{input_basename}.transdecoder.bed"
            
            if pep_file.exists() and cds_file.exists():
                self.logger.info(f"✓ ORF 预测完成")
                self.logger.info(f"  蛋白序列: {pep_file}")
                self.logger.info(f"  CDS 序列: {cds_file}")
                if gff3_file.exists():
                    self.logger.info(f"  GFF3 注释: {gff3_file}")
                if bed_file.exists():
                    self.logger.info(f"  BED 文件: {bed_file}")
                return True
            else:
                self.logger.error("TransDecoder.Predict 未生成预期输出文件")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"TransDecoder.Predict 运行失败: {e}")
            if e.stdout:
                self.logger.error(f"标准输出: {e.stdout}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"运行 TransDecoder.Predict 时发生错误: {e}")
            return False
    
    def _analyze_results(self, input_file: str, output_dir: str):
        """分析 TransDecoder 结果"""
        try:
            input_basename = Path(input_file).name
            gff3_file = Path(output_dir) / f"{input_basename}.transdecoder.gff3"
            
            if not gff3_file.exists():
                self.logger.warning("GFF3 文件不存在，跳过结果分析")
                return
            
            # 统计 ORF 信息
            complete_orfs = 0
            incomplete_orfs = 0
            orfs_5prime_partial = 0
            orfs_3prime_partial = 0
            orfs_internal = 0
            
            with open(gff3_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    feature_type = fields[2]
                    if feature_type != 'CDS':
                        continue
                    
                    attributes = fields[8]
                    
                    # 检查是否为完整 ORF
                    if 'type:complete' in attributes:
                        complete_orfs += 1
                    elif 'type:5prime_partial' in attributes:
                        orfs_5prime_partial += 1
                        incomplete_orfs += 1
                    elif 'type:3prime_partial' in attributes:
                        orfs_3prime_partial += 1
                        incomplete_orfs += 1
                    elif 'type:internal' in attributes:
                        orfs_internal += 1
                        incomplete_orfs += 1
                    else:
                        incomplete_orfs += 1
            
            total_orfs = complete_orfs + incomplete_orfs
            
            self.logger.info("\n" + "=" * 60)
            self.logger.info("ORF 完整度统计:")
            self.logger.info(f"  总 ORF 数: {total_orfs}")
            self.logger.info(f"  完整 ORF: {complete_orfs} ({complete_orfs/total_orfs*100:.2f}%)" if total_orfs > 0 else "  完整 ORF: 0")
            self.logger.info(f"  不完整 ORF: {incomplete_orfs} ({incomplete_orfs/total_orfs*100:.2f}%)" if total_orfs > 0 else "  不完整 ORF: 0")
            if orfs_5prime_partial > 0:
                self.logger.info(f"    - 5'端缺失: {orfs_5prime_partial}")
            if orfs_3prime_partial > 0:
                self.logger.info(f"    - 3'端缺失: {orfs_3prime_partial}")
            if orfs_internal > 0:
                self.logger.info(f"    - 内部片段: {orfs_internal}")
            self.logger.info("=" * 60)
            
        except Exception as e:
            self.logger.warning(f"分析结果时出错: {e}")
    
    def parse_orf_completeness(self, gff3_file: str) -> Dict[str, Dict[str, Any]]:
        """
        解析 GFF3 文件，获取每个转录本的 ORF 完整度信息（仅第一个 ORF，向后兼容）
        
        Args:
            gff3_file: TransDecoder 生成的 GFF3 文件路径
        
        Returns:
            字典：{transcript_id: {'is_complete': bool, 'type': str, 'length': int, ...}}
        """
        orf_info = {}
        
        if not os.path.exists(gff3_file):
            self.logger.warning(f"GFF3 文件不存在: {gff3_file}")
            return orf_info
        
        try:
            with open(gff3_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # 只处理 CDS 特征
                    if fields[2] != 'CDS':
                        continue
                    
                    transcript_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    attributes = fields[8]
                    
                    # 解析 ORF 类型
                    orf_type = 'unknown'
                    is_complete = False
                    
                    if 'type:complete' in attributes:
                        orf_type = 'complete'
                        is_complete = True
                    elif 'type:5prime_partial' in attributes:
                        orf_type = '5prime_partial'
                    elif 'type:3prime_partial' in attributes:
                        orf_type = '3prime_partial'
                    elif 'type:internal' in attributes:
                        orf_type = 'internal'
                    
                    # 存储信息（仅第一个 ORF）
                    if transcript_id not in orf_info:
                        orf_info[transcript_id] = {
                            'is_complete': is_complete,
                            'type': orf_type,
                            'start': start,
                            'end': end,
                            'length': end - start + 1
                        }
            
            self.logger.info(f"解析了 {len(orf_info)} 个转录本的 ORF 信息")
            
        except Exception as e:
            self.logger.error(f"解析 GFF3 文件时出错: {e}")
        
        return orf_info
    
    def parse_all_orfs(self, gff3_file: str, complete_only: bool = True) -> Dict[str, List[Dict[str, Any]]]:
        """
        解析 GFF3 文件，获取每个转录本的所有 ORF 信息（支持多 ORF）
        
        Args:
            gff3_file: TransDecoder 生成的 GFF3 文件路径
            complete_only: 是否只返回完整的 ORF（type:complete），默认 True
        
        Returns:
            字典：{transcript_id: [orf1, orf2, ...]}
            每个 ORF 包含：{'type': str, 'start': int, 'end': int, 'length': int, 
                          'strand': str, 'attributes': str, 'orf_id': str}
        """
        all_orfs = {}
        
        if not os.path.exists(gff3_file):
            self.logger.warning(f"GFF3 文件不存在: {gff3_file}")
            return all_orfs
        
        try:
            with open(gff3_file, 'r') as f:
                for line in f:
                    if line.startswith('#') or not line.strip():
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    
                    # 只处理 CDS 特征
                    if fields[2] != 'CDS':
                        continue
                    
                    transcript_id = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    attributes = fields[8]
                    
                    # 解析 ORF 类型
                    orf_type = 'unknown'
                    if 'type:complete' in attributes:
                        orf_type = 'complete'
                    elif 'type:5prime_partial' in attributes:
                        orf_type = '5prime_partial'
                    elif 'type:3prime_partial' in attributes:
                        orf_type = '3prime_partial'
                    elif 'type:internal' in attributes:
                        orf_type = 'internal'
                    
                    # 如果只要完整 ORF，跳过非完整的
                    if complete_only and orf_type != 'complete':
                        continue
                    
                    # 提取 ORF ID（如果有）
                    orf_id = transcript_id
                    for attr in attributes.split(';'):
                        if attr.startswith('ID='):
                            orf_id = attr.split('=')[1]
                            break
                    
                    # 创建 ORF 信息字典
                    orf_info = {
                        'type': orf_type,
                        'start': start,
                        'end': end,
                        'length': end - start + 1,
                        'strand': strand,
                        'attributes': attributes,
                        'orf_id': orf_id
                    }
                    
                    # 添加到对应转录本的 ORF 列表
                    if transcript_id not in all_orfs:
                        all_orfs[transcript_id] = []
                    all_orfs[transcript_id].append(orf_info)
            
            # 统计信息
            total_orfs = sum(len(orfs) for orfs in all_orfs.values())
            multi_orf_count = sum(1 for orfs in all_orfs.values() if len(orfs) > 1)
            
            self.logger.info(f"解析了 {len(all_orfs)} 个转录本的 {total_orfs} 个 ORF")
            if complete_only:
                self.logger.info(f"  - 仅包含完整 ORF (type:complete)")
            self.logger.info(f"  - 单 ORF 转录本: {len(all_orfs) - multi_orf_count}")
            self.logger.info(f"  - 多 ORF 转录本: {multi_orf_count}")
            
        except Exception as e:
            self.logger.error(f"解析 GFF3 文件时出错: {e}")
        
        return all_orfs

