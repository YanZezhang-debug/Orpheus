"""
BUSCO 质量评估模块
用于评估转录组组装的完整性和质量
"""

import os
import re
import subprocess
import logging
from pathlib import Path
from typing import Dict, Any, Optional, Tuple


class BUSCORunner:
    """
    BUSCO 运行器
    
    BUSCO (Benchmarking Universal Single-Copy Orthologs) 用于评估基因组/转录组的完整性。
    它通过检查单拷贝直系同源基因的存在情况来评估组装质量。
    """
    
    def __init__(self, config: Dict[str, Any], logger: Optional[logging.Logger] = None):
        """
        初始化 BUSCO 运行器
        
        Args:
            config: BUSCO 配置字典
            logger: 日志记录器
        """
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        
        # BUSCO 配置
        self.executable = config.get('executable', 'busco')
        self.lineage = config.get('lineage', '')
        self.mode = config.get('mode', 'transcriptome')
        self.threads = config.get('threads', 8)
        self.extra_params = config.get('extra_params', '')
        
    def check_installation(self) -> bool:
        """
        检查 BUSCO 是否已安装
        
        Returns:
            True 如果 BUSCO 可用，False 否则
        """
        try:
            result = subprocess.run(
                [self.executable, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            
            if result.returncode == 0:
                # 提取版本号
                version_match = re.search(r'(\d+\.\d+\.\d+)', result.stdout)
                version = version_match.group(1) if version_match else "未知版本"
                self.logger.info(f"✓ 检测到 BUSCO {version}")
                return True
            else:
                self.logger.warning(f"BUSCO 可能未正确安装")
                return False
                
        except FileNotFoundError:
            self.logger.error(f"未找到 BUSCO: {self.executable}")
            self.logger.error("请确保 BUSCO 已安装并在 PATH 中，或在配置文件中指定正确路径")
            return False
        except Exception as e:
            self.logger.error(f"检查 BUSCO 安装时出错: {e}")
            return False
    
    def list_lineages(self) -> bool:
        """
        列出可用的 BUSCO lineage 数据集
        
        Returns:
            True 如果成功，False 否则
        """
        self.logger.info("查询可用的 BUSCO lineage 数据集...")
        
        try:
            result = subprocess.run(
                [self.executable, '--list-datasets'],
                capture_output=True,
                text=True,
                timeout=30
            )
            
            if result.returncode == 0:
                self.logger.info("可用的 BUSCO lineage 数据集:")
                self.logger.info(result.stdout)
                return True
            else:
                self.logger.error("无法列出 BUSCO lineage 数据集")
                if result.stderr:
                    self.logger.error(result.stderr)
                return False
                
        except Exception as e:
            self.logger.error(f"列出 BUSCO lineage 时出错: {e}")
            return False
    
    def run(self, input_file: str, output_dir: str, run_name: Optional[str] = None) -> bool:
        """
        运行 BUSCO 质量评估
        
        Args:
            input_file: 输入序列文件（FASTA 格式）
            output_dir: 输出目录
            run_name: 运行名称（用于标识结果，默认自动生成）
        
        Returns:
            True 如果成功，False 否则
        """
        if not self.config.get('enabled', False):
            self.logger.info("BUSCO 评估未启用，跳过")
            return True
        
        # 检查 lineage 配置
        if not self.lineage:
            self.logger.warning("未配置 BUSCO lineage，跳过质量评估")
            self.logger.info("提示: 请在配置文件中设置 busco.lineage")
            self.logger.info("      例如: eukaryota_odb10, metazoa_odb10, viridiplantae_odb10 等")
            return True
        
        # 检查输入文件
        if not os.path.exists(input_file):
            self.logger.error(f"输入文件不存在: {input_file}")
            return False
        
        # 创建输出目录
        os.makedirs(output_dir, exist_ok=True)
        
        # 生成运行名称
        if not run_name:
            run_name = Path(input_file).stem
        
        self.logger.info("\n" + "=" * 60)
        self.logger.info("开始 BUSCO 质量评估")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件: {input_file}")
        self.logger.info(f"输出目录: {output_dir}")
        self.logger.info(f"Lineage: {self.lineage}")
        self.logger.info(f"模式: {self.mode}")
        self.logger.info(f"线程数: {self.threads}")
        
        # 构建命令
        cmd = [
            self.executable,
            '-i', input_file,
            '-o', run_name,
            '-l', self.lineage,
            '-m', self.mode,
            '--cpu', str(self.threads),
            '--out_path', output_dir,
            '--force'  # 覆盖已有结果
        ]
        
        # 添加额外参数
        if self.extra_params:
            cmd.extend(self.extra_params.split())
        
        self.logger.info(f"命令: {' '.join(cmd)}")
        
        try:
            # 运行 BUSCO
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=3600  # 1小时超时
            )
            
            # 记录输出
            if result.stdout:
                self.logger.debug(f"BUSCO 标准输出:\n{result.stdout}")
            
            if result.returncode == 0:
                self.logger.info("✓ BUSCO 评估完成")
                
                # 解析并显示结果
                summary_file = self._find_summary_file(output_dir, run_name)
                if summary_file:
                    self._display_summary(summary_file)
                
                return True
            else:
                self.logger.error("BUSCO 运行失败")
                if result.stderr:
                    self.logger.error(f"错误信息: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            self.logger.error("BUSCO 运行超时（超过1小时）")
            return False
        except Exception as e:
            self.logger.error(f"运行 BUSCO 时出错: {e}")
            return False
    
    def _find_summary_file(self, output_dir: str, run_name: str) -> Optional[str]:
        """
        查找 BUSCO 结果摘要文件
        
        Args:
            output_dir: 输出目录
            run_name: 运行名称
        
        Returns:
            摘要文件路径，如果未找到则返回 None
        """
        # BUSCO v5+ 的摘要文件路径
        summary_patterns = [
            f"{run_name}/short_summary.specific.*.{run_name}.txt",
            f"{run_name}/short_summary.txt",
            f"short_summary.specific.*.{run_name}.txt"
        ]
        
        for pattern in summary_patterns:
            import glob
            matches = glob.glob(os.path.join(output_dir, pattern))
            if matches:
                return matches[0]
        
        return None
    
    def _display_summary(self, summary_file: str) -> None:
        """
        解析并显示 BUSCO 结果摘要
        
        Args:
            summary_file: 摘要文件路径
        """
        try:
            with open(summary_file, 'r') as f:
                content = f.read()
            
            self.logger.info("\n" + "=" * 60)
            self.logger.info("BUSCO 评估结果")
            self.logger.info("=" * 60)
            
            # 提取关键统计信息
            stats = self._parse_summary(content)
            
            if stats:
                total = stats.get('total', 0)
                complete = stats.get('complete', 0)
                complete_single = stats.get('complete_single', 0)
                complete_duplicated = stats.get('complete_duplicated', 0)
                fragmented = stats.get('fragmented', 0)
                missing = stats.get('missing', 0)
                
                self.logger.info(f"完整 BUSCOs (C): {complete}/{total} ({complete/total*100:.1f}%)")
                if complete > 0:
                    self.logger.info(f"  - 单拷贝 (S): {complete_single}")
                    self.logger.info(f"  - 多拷贝 (D): {complete_duplicated}")
                self.logger.info(f"片段化 (F): {fragmented}/{total} ({fragmented/total*100:.1f}%)")
                self.logger.info(f"缺失 (M): {missing}/{total} ({missing/total*100:.1f}%)")
                
                # 质量评价
                completeness = complete / total * 100 if total > 0 else 0
                self.logger.info("")
                if completeness >= 90:
                    self.logger.info("✓ 组装质量: 优秀")
                elif completeness >= 80:
                    self.logger.info("✓ 组装质量: 良好")
                elif completeness >= 70:
                    self.logger.info("✓ 组装质量: 中等")
                else:
                    self.logger.info("⚠ 组装质量: 需要改进")
                
                self.logger.info(f"\n详细报告: {summary_file}")
            else:
                # 如果无法解析，直接显示文件内容
                self.logger.info(content)
                
        except Exception as e:
            self.logger.warning(f"无法读取 BUSCO 摘要文件: {e}")
    
    def _parse_summary(self, content: str) -> Optional[Dict[str, int]]:
        """
        解析 BUSCO 摘要内容
        
        Args:
            content: 摘要文件内容
        
        Returns:
            统计字典，如果解析失败返回 None
        """
        try:
            stats = {}
            
            # 匹配 BUSCO 结果行
            # 例如: C:90.5%[S:85.2%,D:5.3%],F:5.2%,M:4.3%,n:1440
            pattern = r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
            match = re.search(pattern, content)
            
            if match:
                c_pct, s_pct, d_pct, f_pct, m_pct, total = match.groups()
                total = int(total)
                
                stats['total'] = total
                stats['complete'] = int(float(c_pct) * total / 100)
                stats['complete_single'] = int(float(s_pct) * total / 100)
                stats['complete_duplicated'] = int(float(d_pct) * total / 100)
                stats['fragmented'] = int(float(f_pct) * total / 100)
                stats['missing'] = int(float(m_pct) * total / 100)
                
                return stats
            
            # 尝试匹配其他格式
            # 例如: 1234 Complete BUSCOs (C)
            complete_match = re.search(r'(\d+)\s+Complete BUSCOs', content)
            single_match = re.search(r'(\d+)\s+Complete and single-copy', content)
            dup_match = re.search(r'(\d+)\s+Complete and duplicated', content)
            frag_match = re.search(r'(\d+)\s+Fragmented', content)
            miss_match = re.search(r'(\d+)\s+Missing', content)
            total_match = re.search(r'(\d+)\s+Total BUSCO groups', content)
            
            if all([complete_match, total_match]):
                stats['complete'] = int(complete_match.group(1))
                stats['complete_single'] = int(single_match.group(1)) if single_match else 0
                stats['complete_duplicated'] = int(dup_match.group(1)) if dup_match else 0
                stats['fragmented'] = int(frag_match.group(1)) if frag_match else 0
                stats['missing'] = int(miss_match.group(1)) if miss_match else 0
                stats['total'] = int(total_match.group(1))
                
                return stats
            
            return None
            
        except Exception as e:
            self.logger.debug(f"解析 BUSCO 摘要时出错: {e}")
            return None
    
    def compare_results(self, before_file: str, after_file: str) -> Tuple[Optional[Dict], Optional[Dict]]:
        """
        比较处理前后的 BUSCO 结果
        
        Args:
            before_file: 处理前的摘要文件
            after_file: 处理后的摘要文件
        
        Returns:
            (before_stats, after_stats) 元组
        """
        before_stats = None
        after_stats = None
        
        try:
            if os.path.exists(before_file):
                with open(before_file, 'r') as f:
                    before_stats = self._parse_summary(f.read())
            
            if os.path.exists(after_file):
                with open(after_file, 'r') as f:
                    after_stats = self._parse_summary(f.read())
            
            if before_stats and after_stats:
                self.logger.info("\n" + "=" * 60)
                self.logger.info("BUSCO 结果对比")
                self.logger.info("=" * 60)
                
                total = after_stats.get('total', 1)
                
                for key in ['complete', 'fragmented', 'missing']:
                    before_val = before_stats.get(key, 0)
                    after_val = after_stats.get(key, 0)
                    diff = after_val - before_val
                    
                    label = {'complete': '完整', 'fragmented': '片段化', 'missing': '缺失'}[key]
                    
                    before_pct = before_val / total * 100
                    after_pct = after_val / total * 100
                    diff_pct = diff / total * 100
                    
                    sign = '+' if diff > 0 else ''
                    self.logger.info(
                        f"{label}: {before_val} ({before_pct:.1f}%) → "
                        f"{after_val} ({after_pct:.1f}%) [{sign}{diff}, {sign}{diff_pct:.1f}%]"
                    )
                
        except Exception as e:
            self.logger.error(f"比较 BUSCO 结果时出错: {e}")
        
        return before_stats, after_stats

