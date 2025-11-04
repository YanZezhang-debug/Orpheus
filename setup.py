#!/usr/bin/env python3
"""
Orpheus 安装脚本
"""

from setuptools import setup, find_packages
from pathlib import Path

# 读取 README
readme_file = Path(__file__).parent / "README.md"
long_description = ""
if readme_file.exists():
    with open(readme_file, "r", encoding="utf-8") as f:
        long_description = f.read()

# 读取版本号
version = {}
with open("orpheus/__init__.py", "r", encoding="utf-8") as f:
    for line in f:
        if line.startswith("__version__"):
            exec(line, version)
            break

setup(
    name="orpheus-bio",
    version=version.get("__version__", "0.1.0"),
    author="Your Name",
    author_email="your.email@example.com",
    description="无参转录组组装转录本可信性评估工具",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/orpheus",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=[
        "PyYAML>=5.1",
        "biopython>=1.78",
    ],
    entry_points={
        "console_scripts": [
            "orpheus=orpheus_cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["config/*.yaml"],
    },
    zip_safe=False,
    keywords="bioinformatics transcriptome assembly orf prediction quality-assessment",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/orpheus/issues",
        "Source": "https://github.com/yourusername/orpheus",
    },
)

