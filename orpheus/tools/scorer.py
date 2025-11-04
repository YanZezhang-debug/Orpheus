"""
转录本可信性评分模块
综合 BUSCO 完整性评估和同源性证据，对转录本进行评分筛选
"""

import os
import re
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict


class TranscriptScorer:
    """
    转录本评分器
    
    综合多种证据对转录本进行评分：
    1. BUSCO 基因匹配（是否匹配到保守的单拷贝直系同源基因）- 最强权重
    2. 同源性证据（来自 DIAMOND/BLASTP 结果）- 中等权重
    3. ORF 完整性（起始密码子 + 终止密码子）- 中等权重
    4. ORF 长度 - 弱权重（基础过滤）
    
    评分策略：
    - BUSCO 基因：如果转录本匹配到 BUSCO 保守基因，说明它是重要的核心基因
    - 同源性：有已知蛋白匹配，功能可推断
    - ORF 完整性：有起始和终止密码子，可以翻译完整蛋白
    - ORF 长度：足够长度降低随机 ORF 的可能性
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """
        初始化评分器
        
        Args:
            logger: 日志记录器
        """
        self.logger = logger or logging.getLogger(__name__)
    
    def parse_busco_results(self, busco_dir: str) -> Tuple[Set[str], Dict[str, Dict[str, str]]]:
        """
        解析 BUSCO 结果，提取匹配到 BUSCO 保守基因的转录本及其详细信息
        
        Args:
            busco_dir: BUSCO 输出目录
            
        Returns:
            元组：(匹配到 BUSCO 基因的转录本 ID 集合, BUSCO 详细信息字典)
            BUSCO 详细信息字典格式: {transcript_id: {'busco_id': str, 'status': str, 'score': str, 'length': str}}
        """
        self.logger.info(f"解析 BUSCO 结果: {busco_dir}")
        
        busco_genes = set()
        busco_details = {}
        
        # 查找 full_table.tsv 文件
        # BUSCO v5+ 路径可能是: 
        #   - busco_output/run_name/full_table.tsv (旧版)
        #   - busco_output/run_name/run_lineage/full_table.tsv (新版)
        full_table_patterns = [
            os.path.join(busco_dir, "**/full_table.tsv"),  # 递归查找所有子目录
            os.path.join(busco_dir, "*/full_table.tsv"),    # 一层子目录
            os.path.join(busco_dir, "full_table.tsv"),       # 直接在根目录
        ]
        
        full_table_file = None
        import glob
        for pattern in full_table_patterns:
            matches = glob.glob(pattern, recursive=True)  # 启用递归匹配
            if matches:
                full_table_file = matches[0]
                break
        
        if not full_table_file:
            self.logger.warning(f"未找到 BUSCO full_table.tsv 文件: {busco_dir}")
            self.logger.warning("将不使用 BUSCO 基因匹配信息进行评分")
            return busco_genes, busco_details
        
        self.logger.info(f"找到 BUSCO 详细表: {full_table_file}")
        
        # 解析 full_table.tsv
        # 格式：
        # # Busco id	Status	Sequence	Score	Length
        # EOG09360001	Complete	TRINITY_DN1000_c0_g1_i1.p1	1234.5	567
        # EOG09360002	Fragmented	TRINITY_DN1001_c0_g1_i1.p1	456.7	234
        # EOG09360003	Missing	-	-	-
        
        try:
            with open(full_table_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    
                    # 跳过注释和空行
                    if not line or line.startswith('#'):
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) < 3:
                        continue
                    
                    busco_id = fields[0]
                    status = fields[1]
                    sequence = fields[2]
                    score = fields[3] if len(fields) > 3 else '-'
                    length = fields[4] if len(fields) > 4 else '-'
                    
                    # 保留 Complete, Duplicated 和 Fragmented 状态的 BUSCO 基因
                    # Complete/Duplicated: 完整匹配，最高分
                    # Fragmented: 部分匹配，中等分（仍然是有价值的同源性证据）
                    if status in ['Complete', 'Duplicated', 'Fragmented'] and sequence != '-':
                        # 从 TRINITY_DN1000_c0_g1_i1.p1 提取 TRINITY_DN1000_c0_g1_i1
                        transcript_id = sequence
                        if transcript_id.endswith('.p1') or transcript_id.endswith('.p2'):
                            transcript_id = transcript_id.rsplit('.', 1)[0]
                        
                        busco_genes.add(transcript_id)
                        
                        # 存储详细信息
                        busco_details[transcript_id] = {
                            'busco_id': busco_id,
                            'status': status,
                            'score': score,
                            'length': length
                        }
            
            self.logger.info(f"找到 {len(busco_genes)} 个匹配 BUSCO 保守基因的转录本")
            
            # 统计各状态的数量
            status_counts = {}
            for detail in busco_details.values():
                status = detail.get('status', 'Unknown')
                status_counts[status] = status_counts.get(status, 0) + 1
            
            if status_counts:
                self.logger.info("  BUSCO 匹配状态分布:")
                for status, count in sorted(status_counts.items()):
                    self.logger.info(f"    {status}: {count}")
            
            # 显示部分示例（前5个）
            if busco_genes:
                examples = list(busco_genes)[:5]
                self.logger.debug(f"示例: {', '.join(examples)}")
            
        except Exception as e:
            self.logger.error(f"解析 BUSCO 结果时出错: {e}")
            return set(), {}
        
        return busco_genes, busco_details
        
    def parse_gff3(self, gff3_file: str) -> Dict[str, Dict]:
        """
        解析 TransDecoder 生成的 GFF3 文件，提取 ORF 信息
        
        Args:
            gff3_file: GFF3 文件路径
            
        Returns:
            字典，键为转录本 ID，值为 ORF 信息
            {
                'transcript_id': {
                    'gene': 'gene_id',
                    'type': 'complete' | 'internal' | '5prime_partial' | '3prime_partial',
                    'length': int,  # 氨基酸长度
                    'score': float,
                    'coords': (start, end),
                }
            }
        """
        self.logger.info(f"解析 GFF3 文件: {gff3_file}")
        
        if not os.path.exists(gff3_file):
            self.logger.error(f"GFF3 文件不存在: {gff3_file}")
            return {}
        
        orf_info = {}
        mrna_info = {}  # 存储 mRNA 的类型信息，key 是 mRNA ID
        
        # 第一遍：解析 mRNA 行，提取类型信息
        with open(gff3_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                
                seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
                
                # 只处理 mRNA 特征
                if feature_type != 'mRNA':
                    continue
                
                # 解析属性
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                mrna_id = attr_dict.get('ID', '')
                if not mrna_id:
                    continue
                
                # 从 Name 字段提取 ORF 类型
                # Name 格式: ORF%20...%20type%3Acomplete%20len%3A407%20...
                name = attr_dict.get('Name', '')
                orf_type = 'unknown'
                if 'type%3Acomplete' in name or 'type:complete' in name:
                    orf_type = 'complete'
                elif 'type%3Ainternal' in name or 'type:internal' in name:
                    orf_type = 'internal'
                elif 'type%3A5prime_partial' in name or 'type:5prime_partial' in name:
                    orf_type = '5prime_partial'
                elif 'type%3A3prime_partial' in name or 'type:3prime_partial' in name:
                    orf_type = '3prime_partial'
                
                mrna_info[mrna_id] = {
                    'type': orf_type,
                    'transcript_id': seqid
                }
        
        # 第二遍：解析 CDS 行，提取坐标信息，并关联 mRNA 类型
        with open(gff3_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                if not line or line.startswith('#'):
                    continue
                
                fields = line.split('\t')
                if len(fields) < 9:
                    continue
                
                seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
                
                # 只处理 CDS 特征
                if feature_type != 'CDS':
                    continue
                
                # 解析属性
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                # Parent 指向 mRNA ID
                mrna_id = attr_dict.get('Parent', '')
                transcript_id = seqid  # 转录本 ID 是序列 ID
                
                if not mrna_id:
                    continue
                
                # 从 mRNA 信息中获取 ORF 类型
                orf_type = 'unknown'
                if mrna_id in mrna_info:
                    orf_type = mrna_info[mrna_id]['type']
                
                # 计算 ORF 长度（氨基酸数）
                orf_length = (int(end) - int(start) + 1) // 3
                
                # 保存 ORF 信息（如果同一个转录本有多个 ORF，保留最长的）
                if transcript_id not in orf_info or orf_length > orf_info[transcript_id]['length']:
                    orf_info[transcript_id] = {
                        'gene': mrna_id,
                        'type': orf_type,
                        'length': orf_length,
                        'score': float(score) if score != '.' else 0.0,
                        'coords': (int(start), int(end)),
                        'strand': strand,
                    }
        
        self.logger.info(f"✓ 解析完成，找到 {len(orf_info)} 个转录本的 ORF 信息")
        
        # 统计 ORF 类型分布
        type_counts = defaultdict(int)
        for info in orf_info.values():
            type_counts[info['type']] += 1
        
        self.logger.info("  ORF 类型分布:")
        for orf_type, count in sorted(type_counts.items()):
            self.logger.info(f"    {orf_type}: {count}")
        
        return orf_info
    
    def parse_homology_results(self, homology_file: str) -> Tuple[Set[str], Dict[str, Dict[str, str]]]:
        """
        解析同源搜索结果（DIAMOND/BLASTP outfmt6 格式）
        
        Args:
            homology_file: 同源搜索结果文件
            
        Returns:
            元组：(有同源证据的基因 ID 集合, 同源详细信息字典)
            同源详细信息字典格式: {gene_id: {'subject': str, 'identity': str, 'evalue': str, 'bitscore': str}}
        """
        self.logger.info(f"解析同源搜索结果: {homology_file}")
        
        if not os.path.exists(homology_file):
            self.logger.warning(f"同源搜索结果文件不存在: {homology_file}")
            return set(), {}
        
        homology_genes = set()
        homology_details = {}
        
        with open(homology_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # outfmt6 格式: query, subject, identity, length, ...
                fields = line.split('\t')
                if len(fields) < 2:
                    continue
                
                query_id = fields[0]
                subject_id = fields[1]
                identity = fields[2] if len(fields) > 2 else '-'
                evalue = fields[10] if len(fields) > 10 else '-'
                bitscore = fields[11] if len(fields) > 11 else '-'
                
                # TransDecoder 输出的 query ID 格式: GENE.p1, GENE.p2 等
                # 提取基因 ID
                gene_id = re.sub(r'\.p\d+$', '', query_id)
                homology_genes.add(gene_id)
                
                # 只保留最佳匹配（第一个出现的，因为通常已按 bitscore 排序）
                if gene_id not in homology_details:
                    homology_details[gene_id] = {
                        'subject': subject_id,
                        'identity': identity,
                        'evalue': evalue,
                        'bitscore': bitscore
                    }
        
        self.logger.info(f"✓ 找到 {len(homology_genes)} 个基因有同源证据")
        
        return homology_genes, homology_details
    
    def calculate_scores(self, 
                        orf_info: Dict[str, Dict],
                        homology_genes: Set[str],
                        busco_genes: Optional[Set[str]] = None,
                        busco_details: Optional[Dict[str, Dict[str, str]]] = None,
                        weights: Optional[Dict[str, float]] = None) -> Dict[str, float]:
        """
        计算转录本综合得分
        
        评分标准：
        1. BUSCO 基因匹配：
           - Complete/Duplicated: 1.0（完整匹配保守基因）
           - Fragmented: 0.5（部分匹配，仍是有价值的同源性证据）
           - 无匹配: 0.0
        2. ORF 完整性：complete (1.0) > 5'/3' partial (0.6) > internal (0.3) > unknown (0.0)
        3. 同源性证据：有同源匹配 (1.0) vs 无 (0.0)
        4. ORF 长度：标准化到 [0, 1]
        
        Args:
            orf_info: ORF 信息字典（来自 parse_gff3）
            homology_genes: 有同源证据的基因 ID 集合（来自 parse_homology_results）
            busco_genes: 匹配 BUSCO 保守基因的转录本 ID 集合（来自 parse_busco_results）
            busco_details: BUSCO 详细信息字典（包含状态信息）
            weights: 权重字典 {'busco': float, 'completeness': float, 'homology': float, 'length': float}
            
        Returns:
            字典，键为转录本 ID，值为综合得分
        """
        self.logger.info("\n" + "=" * 60)
        self.logger.info("计算转录本综合得分")
        self.logger.info("=" * 60)
        
        # 默认权重
        if weights is None:
            weights = {
                'busco': 0.4,          # BUSCO 基因匹配权重（最强）
                'completeness': 0.3,   # ORF 完整性权重
                'homology': 0.2,       # 同源性权重
                'length': 0.1,         # 长度权重（最弱）
            }
        
        # 如果没有 BUSCO 数据，重新分配权重
        if not busco_genes:
            self.logger.warning("未提供 BUSCO 基因匹配数据，重新分配权重")
            weights = {
                'busco': 0.0,
                'completeness': 0.5,
                'homology': 0.3,
                'length': 0.2,
            }
        
        self.logger.info("评分权重:")
        for criterion, weight in weights.items():
            self.logger.info(f"  {criterion}: {weight}")
        
        # 初始化 busco_genes 为空集合（如果没有提供）
        if busco_genes is None:
            busco_genes = set()
        if busco_details is None:
            busco_details = {}
        
        # BUSCO 状态得分
        busco_status_scores = {
            'Complete': 1.0,      # 完整匹配保守基因
            'Duplicated': 1.0,    # 完整匹配但有多个拷贝
            'Fragmented': 0.5,    # 部分匹配，仍是有价值的同源性证据
        }
        
        # ORF 完整性得分
        completeness_scores = {
            'complete': 1.0,
            '5prime_partial': 0.6,
            '3prime_partial': 0.6,
            'internal': 0.3,
            'unknown': 0.0,
        }
        
        # 找出最长的 ORF 用于归一化
        max_length = max((info['length'] for info in orf_info.values()), default=1)
        
        scores = {}
        
        for transcript_id, info in orf_info.items():
            # 1. BUSCO 基因匹配得分（根据状态给不同的分数）
            busco_score = 0.0
            if transcript_id in busco_genes and transcript_id in busco_details:
                status = busco_details[transcript_id].get('status', '')
                busco_score = busco_status_scores.get(status, 0.0)
            
            # 2. ORF 完整性得分
            comp_score = completeness_scores.get(info['type'], 0.0)
            
            # 3. 同源性得分
            gene_id = info['gene']
            homology_score = 1.0 if gene_id in homology_genes else 0.0
            
            # 4. 长度得分（归一化）
            length_score = info['length'] / max_length
            
            # 综合得分
            total_score = (
                weights['busco'] * busco_score +
                weights['completeness'] * comp_score +
                weights['homology'] * homology_score +
                weights['length'] * length_score
            )
            
            scores[transcript_id] = total_score
        
        self.logger.info(f"\n✓ 计算完成，共 {len(scores)} 个转录本")
        
        # 统计得分分布
        score_ranges = {
            'excellent (≥0.8)': sum(1 for s in scores.values() if s >= 0.8),
            'good (0.6-0.8)': sum(1 for s in scores.values() if 0.6 <= s < 0.8),
            'fair (0.4-0.6)': sum(1 for s in scores.values() if 0.4 <= s < 0.6),
            'poor (<0.4)': sum(1 for s in scores.values() if s < 0.4),
        }
        
        self.logger.info("\n得分分布:")
        for range_name, count in score_ranges.items():
            percentage = count / len(scores) * 100 if scores else 0
            self.logger.info(f"  {range_name}: {count} ({percentage:.1f}%)")
        
        return scores
    
    def filter_transcripts(self,
                          scores: Dict[str, float],
                          threshold: float = 0.5,
                          top_n: Optional[int] = None) -> List[str]:
        """
        根据得分筛选转录本
        
        Args:
            scores: 转录本得分字典
            threshold: 最低得分阈值（0-1）
            top_n: 保留得分最高的 N 个转录本（None = 不限制）
            
        Returns:
            筛选后的转录本 ID 列表（按得分降序排列）
        """
        self.logger.info("\n" + "=" * 60)
        self.logger.info("筛选高可信度转录本")
        self.logger.info("=" * 60)
        self.logger.info(f"阈值: {threshold}")
        if top_n:
            self.logger.info(f"最多保留: {top_n} 个转录本")
        
        # 按得分降序排序
        sorted_transcripts = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        
        # 应用阈值筛选
        filtered = [(tid, score) for tid, score in sorted_transcripts if score >= threshold]
        
        # 应用数量限制
        if top_n and len(filtered) > top_n:
            filtered = filtered[:top_n]
        
        selected_ids = [tid for tid, score in filtered]
        
        self.logger.info(f"\n✓ 筛选完成")
        self.logger.info(f"  原始转录本数: {len(scores)}")
        self.logger.info(f"  筛选后数量: {len(selected_ids)}")
        
        if len(filtered) > 0:
            avg_score = sum(score for _, score in filtered) / len(filtered)
            min_score = filtered[-1][1]
            max_score = filtered[0][1]
            self.logger.info(f"  得分范围: {min_score:.3f} - {max_score:.3f}")
            self.logger.info(f"  平均得分: {avg_score:.3f}")
        
        return selected_ids
    
    def export_filtered_fasta(self,
                             input_fasta: str,
                             output_fasta: str,
                             selected_ids: List[str]) -> bool:
        """
        导出筛选后的转录本序列
        
        Args:
            input_fasta: 输入 FASTA 文件（CD-HIT 去冗余结果）
            output_fasta: 输出 FASTA 文件
            selected_ids: 要保留的转录本 ID 列表
            
        Returns:
            True 如果成功，False 否则
        """
        self.logger.info(f"\n导出筛选后的序列到: {output_fasta}")
        
        if not os.path.exists(input_fasta):
            self.logger.error(f"输入文件不存在: {input_fasta}")
            return False
        
        selected_set = set(selected_ids)
        exported_count = 0
        
        try:
            with open(input_fasta, 'r') as fin, open(output_fasta, 'w') as fout:
                write_seq = False
                
                for line in fin:
                    if line.startswith('>'):
                        # 提取序列 ID
                        seq_id = line[1:].split()[0]
                        
                        # 检查是否在筛选列表中
                        if seq_id in selected_set:
                            write_seq = True
                            exported_count += 1
                            fout.write(line)
                        else:
                            write_seq = False
                    elif write_seq:
                        fout.write(line)
            
            self.logger.info(f"✓ 成功导出 {exported_count} 个转录本序列")
            return True
            
        except Exception as e:
            self.logger.error(f"导出序列时出错: {e}")
            return False
    
    def score_and_filter(self, 
                        gff3_file: str,
                        cdhit_file: str,
                        output_dir: str,
                        homology_file: Optional[str] = None,
                        busco_dir: Optional[str] = None,
                        threshold: float = 0.5,
                        top_n: Optional[int] = None,
                        weights: Optional[Dict[str, float]] = None) -> bool:
        """
        一站式评分和筛选接口（便捷方法）
        
        Args:
            gff3_file: TransDecoder 生成的 GFF3 文件
            cdhit_file: CD-HIT 去冗余结果文件
            output_dir: 输出目录
            homology_file: 同源搜索结果文件（可选）
            busco_dir: BUSCO 输出目录（可选）
            threshold: 评分阈值
            top_n: 保留的最大数量
            weights: 评分权重
            
        Returns:
            True 如果成功，False 否则
        """
        # 1. 解析 GFF3 文件
        orf_info = self.parse_gff3(gff3_file)
        if not orf_info:
            self.logger.error("无法从 GFF3 文件中提取 ORF 信息")
            return False
        
        # 2. 解析 BUSCO 结果
        busco_genes = set()
        busco_details = {}
        if busco_dir and os.path.exists(busco_dir):
            busco_genes, busco_details = self.parse_busco_results(busco_dir)
        else:
            self.logger.warning("未提供 BUSCO 结果，将不使用 BUSCO 基因匹配信息")
        
        # 3. 解析同源搜索结果
        homology_genes = set()
        homology_details = {}
        if homology_file and os.path.exists(homology_file):
            homology_genes, homology_details = self.parse_homology_results(homology_file)
        else:
            self.logger.warning("未提供同源搜索结果")
        
        # 4. 计算综合得分
        scores = self.calculate_scores(orf_info, homology_genes, busco_genes, busco_details, weights)
        
        # 5. 筛选转录本
        selected_ids = self.filter_transcripts(scores, threshold, top_n)
        
        if not selected_ids:
            self.logger.warning("没有转录本通过筛选阈值")
            return False
        
        # 6. 导出结果
        os.makedirs(output_dir, exist_ok=True)
        
        # 导出筛选后的序列
        output_fasta = os.path.join(output_dir, "high_confidence_transcripts.fasta")
        if not self.export_filtered_fasta(cdhit_file, output_fasta, selected_ids):
            return False
        
        # 导出评分表
        score_table = os.path.join(output_dir, "transcript_scores.tsv")
        self.export_score_table(scores, orf_info, homology_genes, busco_genes, score_table)
        
        # 生成整合报告
        integrated_report = os.path.join(output_dir, "integrated_report.tsv")
        self.generate_integrated_report(scores, orf_info, homology_details, busco_details, integrated_report)
        
        self.logger.info("\n" + "=" * 60)
        self.logger.info("✓ 评分和筛选完成！")
        self.logger.info("=" * 60)
        self.logger.info(f"高可信度转录本: {output_fasta}")
        self.logger.info(f"详细评分表: {score_table}")
        self.logger.info(f"整合报告: {integrated_report}")
        
        return True
    
    def export_score_table(self,
                          scores: Dict[str, float],
                          orf_info: Dict[str, Dict],
                          homology_genes: Set[str],
                          busco_genes: Optional[Set[str]],
                          output_file: str) -> bool:
        """
        导出详细的评分表
        
        Args:
            scores: 转录本得分字典
            orf_info: ORF 信息字典
            homology_genes: 有同源证据的基因 ID 集合
            busco_genes: 匹配 BUSCO 保守基因的转录本 ID 集合
            output_file: 输出文件路径
            
        Returns:
            True 如果成功，False 否则
        """
        self.logger.info(f"\n导出评分详情表到: {output_file}")
        
        # 初始化 busco_genes 为空集合（如果没有提供）
        if busco_genes is None:
            busco_genes = set()
        
        try:
            # 按得分降序排序
            sorted_items = sorted(scores.items(), key=lambda x: x[1], reverse=True)
            
            with open(output_file, 'w') as f:
                # 写入表头
                f.write("Transcript_ID\tScore\tBUSCO_Gene\tORF_Type\tORF_Length\tHomology_Evidence\tStrand\n")
                
                # 写入数据
                for transcript_id, score in sorted_items:
                    info = orf_info.get(transcript_id, {})
                    gene_id = info.get('gene', '')
                    orf_type = info.get('type', 'unknown')
                    orf_length = info.get('length', 0)
                    strand = info.get('strand', '.')
                    is_busco = 'Yes' if transcript_id in busco_genes else 'No'
                    has_homology = 'Yes' if gene_id in homology_genes else 'No'
                    
                    f.write(f"{transcript_id}\t{score:.4f}\t{is_busco}\t{orf_type}\t{orf_length}\t{has_homology}\t{strand}\n")
            
            self.logger.info(f"✓ 成功导出 {len(sorted_items)} 行评分数据")
            return True
            
        except Exception as e:
            self.logger.error(f"导出评分表时出错: {e}")
            return False
    
    def generate_integrated_report(self,
                                   scores: Dict[str, float],
                                   orf_info: Dict[str, Dict],
                                   homology_details: Dict[str, Dict[str, str]],
                                   busco_details: Dict[str, Dict[str, str]],
                                   output_file: str) -> bool:
        """
        生成整合报告，包含所有基因的详细信息
        
        Args:
            scores: 转录本得分字典
            orf_info: ORF 信息字典（从 parse_gff3 返回）
            homology_details: 同源注释详细信息字典
            busco_details: BUSCO 详细信息字典
            output_file: 输出文件路径
            
        Returns:
            True 如果成功，False 否则
        """
        self.logger.info(f"\n生成整合报告: {output_file}")
        
        try:
            # 按基因组织数据
            gene_data = defaultdict(list)
            
            for transcript_id, score in scores.items():
                info = orf_info.get(transcript_id, {})
                gene_id = info.get('gene', transcript_id)
                
                gene_data[gene_id].append({
                    'transcript_id': transcript_id,
                    'score': score,
                    'orf_info': info
                })
            
            # 按得分降序排序
            sorted_genes = sorted(gene_data.items(), 
                                 key=lambda x: max(t['score'] for t in x[1]), 
                                 reverse=True)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                # 写入表头
                header = [
                    "Gene_ID",
                    "Transcript_Count", 
                    "Best_Score",
                    "BUSCO_Status",
                    "BUSCO_Gene_ID",
                    "ORF1_Type",
                    "ORF1_Length",
                    "ORF2_Type",
                    "ORF2_Length",
                    "ORF3_Type",
                    "ORF3_Length",
                    "Homology_Subject",
                    "Homology_Identity",
                    "Homology_Evalue"
                ]
                f.write('\t'.join(header) + '\n')
                
                # 写入数据
                for gene_id, transcripts in sorted_genes:
                    # 按得分排序转录本
                    transcripts.sort(key=lambda x: x['score'], reverse=True)
                    
                    # 获取最佳转录本
                    best_transcript = transcripts[0]
                    best_transcript_id = best_transcript['transcript_id']
                    best_score = best_transcript['score']
                    
                    # BUSCO 信息
                    busco_info = busco_details.get(best_transcript_id, {})
                    busco_status = busco_info.get('status', '-')
                    busco_gene_id = busco_info.get('busco_id', '-')
                    
                    # 同源注释信息
                    homology_info = homology_details.get(gene_id, {})
                    homology_subject = homology_info.get('subject', '-')
                    homology_identity = homology_info.get('identity', '-')
                    homology_evalue = homology_info.get('evalue', '-')
                    
                    # ORF 信息（最多显示前3个）
                    orf_types = []
                    orf_lengths = []
                    for t in transcripts[:3]:
                        orf_types.append(t['orf_info'].get('type', '-'))
                        orf_lengths.append(str(t['orf_info'].get('length', '-')))
                    
                    # 补齐到3列
                    while len(orf_types) < 3:
                        orf_types.append('-')
                        orf_lengths.append('-')
                    
                    # 写入行
                    row = [
                        gene_id,
                        str(len(transcripts)),
                        f"{best_score:.4f}",
                        busco_status,
                        busco_gene_id,
                        orf_types[0],
                        orf_lengths[0],
                        orf_types[1],
                        orf_lengths[1],
                        orf_types[2],
                        orf_lengths[2],
                        homology_subject,
                        homology_identity,
                        homology_evalue
                    ]
                    f.write('\t'.join(row) + '\n')
            
            self.logger.info(f"✓ 成功生成整合报告，包含 {len(sorted_genes)} 个基因")
            return True
            
        except Exception as e:
            self.logger.error(f"生成整合报告时出错: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False


def run_scoring(cdhit_file: str,
                gff3_file: str,
                homology_file: Optional[str],
                busco_dir: Optional[str],
                output_dir: str,
                threshold: float = 0.5,
                top_n: Optional[int] = None,
                weights: Optional[Dict[str, float]] = None,
                logger: Optional[logging.Logger] = None) -> Tuple[bool, Optional[str]]:
    """
    运行转录本评分和筛选流程
    
    Args:
        cdhit_file: CD-HIT 去冗余结果文件
        gff3_file: TransDecoder 生成的 GFF3 文件
        homology_file: 同源搜索结果文件（可选）
        busco_dir: BUSCO 输出目录（可选）
        output_dir: 输出目录
        threshold: 评分阈值
        top_n: 保留的最大数量
        weights: 评分权重
        logger: 日志记录器
        
    Returns:
        (成功标志, 输出文件路径)
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info("\n" + "=" * 60)
    logger.info("转录本可信性评分与筛选")
    logger.info("=" * 60)
    
    # 创建评分器
    scorer = TranscriptScorer(logger)
    
    # 1. 解析 GFF3 文件
    orf_info = scorer.parse_gff3(gff3_file)
    if not orf_info:
        logger.error("无法从 GFF3 文件中提取 ORF 信息")
        return False, None
    
    # 2. 解析 BUSCO 结果
    busco_genes = set()
    busco_details = {}
    if busco_dir and os.path.exists(busco_dir):
        busco_genes, busco_details = scorer.parse_busco_results(busco_dir)
    else:
        logger.warning("未提供 BUSCO 结果，将不使用 BUSCO 基因匹配信息")
    
    # 3. 解析同源搜索结果
    homology_genes = set()
    homology_details = {}
    if homology_file and os.path.exists(homology_file):
        homology_genes, homology_details = scorer.parse_homology_results(homology_file)
    else:
        logger.warning("未提供同源搜索结果，将仅基于 BUSCO、ORF 完整性和长度进行评分")
    
    # 4. 计算综合得分
    scores = scorer.calculate_scores(orf_info, homology_genes, busco_genes, busco_details, weights)
    
    # 4. 筛选转录本
    selected_ids = scorer.filter_transcripts(scores, threshold, top_n)
    
    if not selected_ids:
        logger.warning("没有转录本通过筛选阈值")
        return False, None
    
    # 5. 导出结果
    os.makedirs(output_dir, exist_ok=True)
    
    # 导出筛选后的序列
    output_fasta = os.path.join(output_dir, "high_confidence_transcripts.fasta")
    if not scorer.export_filtered_fasta(cdhit_file, output_fasta, selected_ids):
        return False, None
    
    # 导出评分表
    score_table = os.path.join(output_dir, "transcript_scores.tsv")
    scorer.export_score_table(scores, orf_info, homology_genes, busco_genes, score_table)
    
    # 生成整合报告
    integrated_report = os.path.join(output_dir, "integrated_report.tsv")
    scorer.generate_integrated_report(scores, orf_info, homology_details, busco_details, integrated_report)
    
    logger.info("\n" + "=" * 60)
    logger.info("✓ 评分和筛选完成！")
    logger.info("=" * 60)
    logger.info(f"高可信度转录本: {output_fasta}")
    logger.info(f"详细评分表: {score_table}")
    logger.info(f"整合报告: {integrated_report}")
    
    return True, output_fasta

