"""
CD-HIT-EST 调用模块
用于合并相似的转录本序列
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import Optional, Dict, Any
import logging


class CDHitRunner:
    """CD-HIT-EST 运行器"""
    
    def __init__(self, config: Dict[str, Any], logger: Optional[logging.Logger] = None):
        """
        初始化CD-HIT运行器
        
        Args:
            config: CD-HIT配置字典
            logger: 日志记录器
        """
        self.config = config
        self.logger = logger or logging.getLogger(__name__)
        self.executable = config.get('executable', 'cd-hit-est')
    
    def check_installation(self) -> bool:
        """
        检查CD-HIT-EST是否已安装
        
        Returns:
            True如果已安装，False否则
        """
        try:
            # 方法1: 使用 which/shutil.which 检查命令是否存在
            executable_path = shutil.which(self.executable)
            if executable_path:
                self.logger.info(f"CD-HIT-EST 已找到: {executable_path}")
                return True
            
            # 方法2: 尝试运行命令检查
            # 注意：很多工具在显示帮助信息时返回非零退出码（如255），这是正常的
            result = subprocess.run(
                [self.executable, '-h'],
                capture_output=True,
                text=True,
                timeout=5
            )
            
            # 检查是否有输出（说明命令存在并能运行）
            has_output = bool(result.stdout or result.stderr)
            
            # 检查是否为常见的"命令不存在"错误
            error_output = result.stderr.lower() if result.stderr else ""
            is_not_found = any(phrase in error_output for phrase in [
                'command not found',
                'no such file',
                'cannot find',
                '不是内部或外部命令'
            ])
            
            # 判断：有输出且不是"找不到命令"错误，或者返回码为0/255（常见帮助信息返回码）
            if has_output and not is_not_found:
                # returncode 为 0 或 255 通常表示命令存在（255是很多工具帮助信息的标准返回码）
                if result.returncode in [0, 255]:
                    self.logger.info(f"CD-HIT-EST 已找到: {self.executable}")
                    return True
                # 即使返回码不是0/255，只要有输出且不是"找不到"错误，也认为已安装
                elif has_output:
                    self.logger.info(f"CD-HIT-EST 已找到: {self.executable} (返回码: {result.returncode})")
                    return True
            
            # 如果都失败了
            self.logger.error(f"CD-HIT-EST 未找到: {self.executable}")
            if result.stderr:
                self.logger.debug(f"错误信息: {result.stderr[:200]}")
            return False
            
        except FileNotFoundError:
            self.logger.error(f"CD-HIT-EST 未找到: {self.executable}")
            self.logger.error("请确保 CD-HIT 已安装并在 PATH 中")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"检查 CD-HIT-EST 时超时")
            return False
        except Exception as e:
            self.logger.error(f"检查CD-HIT-EST安装时出错: {e}")
            return False
    
    def build_command(self, input_file: str, output_file: str) -> list:
        """
        构建CD-HIT-EST命令
        
        Args:
            input_file: 输入FASTA文件
            output_file: 输出FASTA文件
        
        Returns:
            命令列表
        """
        cmd = [
            self.executable,
            '-i', input_file,
            '-o', output_file,
            '-c', str(self.config.get('identity', 0.95)),
            '-n', str(self.config.get('word_size', 10)),
            '-M', str(self.config.get('memory', 0)),
            '-T', str(self.config.get('threads', 8)),
            '-aS', str(self.config.get('coverage', 0.9)),
            '-g', str(self.config.get('coverage_mode', 0)),
            '-d', str(self.config.get('description', 0))
        ]
        
        # 添加额外参数
        extra_params = self.config.get('extra_params', '')
        if extra_params:
            cmd.extend(extra_params.split())
        
        return cmd
    
    def run(self, input_file: str, output_file: str) -> bool:
        """
        运行CD-HIT-EST
        
        Args:
            input_file: 输入FASTA文件路径
            output_file: 输出FASTA文件路径
        
        Returns:
            True如果成功，False否则
        """
        # 检查输入文件
        if not os.path.exists(input_file):
            self.logger.error(f"输入文件不存在: {input_file}")
            return False
        
        # 确保输出目录存在
        output_dir = Path(output_file).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建命令
        cmd = self.build_command(input_file, output_file)
        
        self.logger.info("=" * 60)
        self.logger.info("开始运行 CD-HIT-EST")
        self.logger.info(f"输入文件: {input_file}")
        self.logger.info(f"输出文件: {output_file}")
        self.logger.info(f"命令: {' '.join(cmd)}")
        self.logger.info("=" * 60)
        
        try:
            # 运行CD-HIT-EST
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            # 输出标准输出和标准错误
            if result.stdout:
                self.logger.info("CD-HIT-EST 输出:")
                for line in result.stdout.splitlines():
                    self.logger.info(f"  {line}")
            
            if result.stderr:
                self.logger.warning("CD-HIT-EST 错误输出:")
                for line in result.stderr.splitlines():
                    self.logger.warning(f"  {line}")
            
            # 检查输出文件是否生成
            if os.path.exists(output_file):
                self.logger.info(f"CD-HIT-EST 成功完成！")
                self.logger.info(f"结果文件: {output_file}")
                self.logger.info(f"聚类信息: {output_file}.clstr")
                
                # 统计结果
                self._report_stats(input_file, output_file)
                return True
            else:
                self.logger.error("CD-HIT-EST 运行失败：输出文件未生成")
                return False
                
        except subprocess.CalledProcessError as e:
            self.logger.error(f"CD-HIT-EST 运行失败: {e}")
            if e.stdout:
                self.logger.error(f"标准输出: {e.stdout}")
            if e.stderr:
                self.logger.error(f"错误输出: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"运行CD-HIT-EST时发生错误: {e}")
            return False
    
    def _report_stats(self, input_file: str, output_file: str):
        """
        报告统计信息
        
        Args:
            input_file: 输入文件
            output_file: 输出文件
        """
        try:
            input_count = self._count_sequences(input_file)
            output_count = self._count_sequences(output_file)
            
            self.logger.info("=" * 60)
            self.logger.info("统计结果:")
            self.logger.info(f"  输入序列数: {input_count}")
            self.logger.info(f"  输出序列数: {output_count}")
            self.logger.info(f"  去除冗余数: {input_count - output_count}")
            self.logger.info(f"  保留比例: {output_count/input_count*100:.2f}%")
            self.logger.info("=" * 60)
        except Exception as e:
            self.logger.warning(f"统计序列数时出错: {e}")
    
    @staticmethod
    def _count_sequences(fasta_file: str) -> int:
        """
        统计FASTA文件中的序列数量
        
        Args:
            fasta_file: FASTA文件路径
        
        Returns:
            序列数量
        """
        count = 0
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
        return count

