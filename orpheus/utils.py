"""
工具函数模块
"""

import os
import logging
import subprocess
from pathlib import Path
from typing import Optional


def setup_logger(name: str, level: str = "INFO", 
                log_file: Optional[str] = None, 
                console: bool = True) -> logging.Logger:
    """
    设置日志记录器
    
    Args:
        name: 日志记录器名称
        level: 日志级别
        log_file: 日志文件路径
        console: 是否输出到控制台
    
    Returns:
        配置好的日志记录器
    """
    logger = logging.getLogger(name)
    logger.setLevel(getattr(logging, level.upper()))
    
    # 清除已有的处理器
    logger.handlers.clear()
    
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    if console:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    if log_file:
        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger


def check_executable(executable: str) -> bool:
    """
    检查可执行文件是否存在
    
    Args:
        executable: 可执行文件名或路径
    
    Returns:
        True如果存在，False否则
    """
    try:
        result = subprocess.run(
            [executable, '--version'],
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        # 尝试使用which/where命令
        try:
            cmd = 'where' if os.name == 'nt' else 'which'
            result = subprocess.run(
                [cmd, executable],
                capture_output=True,
                timeout=5
            )
            return result.returncode == 0
        except:
            return False


def ensure_dir(directory: str) -> Path:
    """
    确保目录存在，如果不存在则创建
    
    Args:
        directory: 目录路径
    
    Returns:
        Path对象
    """
    dir_path = Path(directory)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def validate_fasta(fasta_file: str) -> bool:
    """
    验证FASTA文件格式
    
    Args:
        fasta_file: FASTA文件路径
    
    Returns:
        True如果格式正确，False否则
    """
    if not os.path.exists(fasta_file):
        return False
    
    try:
        with open(fasta_file, 'r') as f:
            first_line = f.readline().strip()
            return first_line.startswith('>')
    except:
        return False

