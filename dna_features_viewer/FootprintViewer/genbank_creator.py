"""
GenBank文件创建模块

从基因组FASTA和GFF3文件中提取指定区域并创建GenBank文件
"""

import tempfile
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation


class GenBankCreator:
    """GenBank文件创建器"""
    
    def __init__(self):
        # 目标特征类型映射
        self.target_feature_types = {
            'five_prime_UTR': '5UTR',
            'five_prime_utr': '5UTR', 
            '5UTR': '5UTR',
            'CDS': 'CDS',
            'cds': 'CDS',
            'three_prime_UTR': '3UTR',
            'three_prime_utr': '3UTR',
            '3UTR': '3UTR'
        }
    
    def create_from_region(self, genome_fasta, gff3_file, chrom, start, end, output_path=None):
        """
        从基因组FASTA和GFF3文件中提取指定区域并创建GenBank文件
        只读取转录本中的5UTR、CDS、3UTR记录，label形式为转录本-XX
        
        参数:
        - genome_fasta: 基因组FASTA文件路径
        - gff3_file: GFF3注释文件路径  
        - chrom: 染色体名称 (如 'Chr1')
        - start: 起始位置 (1-based)
        - end: 结束位置 (1-based)
        - output_path: 输出GenBank文件路径，如果为None则创建临时文件
        
        返回:
        - GenBank文件路径
        """
        
        # 读取基因组序列
        genome_dict = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
        
        if chrom not in genome_dict:
            raise ValueError(f"染色体 {chrom} 在基因组文件中未找到")
        
        # 提取目标区域序列 (转换为0-based索引)
        region_seq = genome_dict[chrom].seq[start-1:end]
        
        # 创建SeqRecord
        record = SeqRecord(
            region_seq,
            id=f"{chrom}_{start}_{end}",
            description=f"Region {chrom}:{start}-{end}",
            annotations={
                "molecule_type": "DNA",
                "topology": "linear",
                "data_file_division": "PLN",
                "date": "01-JAN-2024",
                "accessions": [f"{chrom}_{start}_{end}"],
                "organism": "Arabidopsis thaliana"
            }
        )
        
        # 解析GFF3文件获取特征
        features = self._parse_gff3_features(gff3_file, chrom, start, end)
        record.features = features
        
        # 保存为GenBank格式
        if output_path is None:
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.gb', delete=False)
            output_path = temp_file.name
            temp_file.close()
        
        SeqIO.write(record, output_path, "genbank")
        return output_path
    
    def _parse_gff3_features(self, gff3_file, chrom, start, end):
        """解析GFF3文件中的特征"""
        features = []
        transcript_map = {}  # 用于映射转录本ID到转录本名称
        
        print(f"正在解析GFF3文件中 {chrom}:{start}-{end} 区域的特征...")
        
        # 第一遍：收集转录本信息
        with open(gff3_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or line.strip() == '':
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                seqname, source, feature_type, start_pos, end_pos, score, strand, phase, attributes = parts
                
                # 只处理目标染色体的mRNA/transcript
                if seqname != chrom or feature_type not in ['mRNA', 'transcript']:
                    continue
                
                # 转换坐标为整数
                try:
                    feat_start = int(start_pos)  # GFF3是1-based
                    feat_end = int(end_pos)
                except ValueError:
                    continue
                
                # 检查是否与目标区域重叠
                if feat_start > end or feat_end < start:
                    continue
                    
                # 解析attributes
                attr_dict = self._parse_gff_attributes(attributes)
                
                # 收集转录本信息
                transcript_id = attr_dict.get('ID', '')
                transcript_name = attr_dict.get('Name', attr_dict.get('ID', ''))
                if transcript_id and transcript_name:
                    transcript_map[transcript_id] = transcript_name
        
        print(f"收集到 {len(transcript_map)} 个转录本映射")
        
        # 第二遍：处理目标特征类型
        feature_count = 0
        target_feature_count = 0
        
        with open(gff3_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or line.strip() == '':
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                seqname, source, feature_type, start_pos, end_pos, score, strand, phase, attributes = parts
                
                # 只处理目标染色体
                if seqname != chrom:
                    continue
                    
                feature_count += 1
                
                # 只处理目标特征类型
                if feature_type not in self.target_feature_types:
                    continue
                    
                target_feature_count += 1
                
                # 转换坐标为整数
                try:
                    feat_start = int(start_pos)  # GFF3是1-based
                    feat_end = int(end_pos)
                except ValueError:
                    continue
                
                # 检查是否与目标区域重叠
                if feat_start > end or feat_end < start:
                    continue
                    
                # 解析attributes
                attr_dict = self._parse_gff_attributes(attributes)
                
                # 计算在目标区域内的相对坐标，处理截断情况
                relative_start = max(0, feat_start - start)
                relative_end = min(end - start + 1, feat_end - start + 1)
                
                # 确保有效的坐标范围
                if relative_end > relative_start:
                    # 创建SeqFeature
                    feature_obj = self._create_seq_feature(
                        feature_type, relative_start, relative_end, strand,
                        attr_dict, transcript_map
                    )
                    features.append(feature_obj)
        
        print(f"处理了 {feature_count} 个特征，找到 {target_feature_count} 个目标类型特征")
        print(f"最终添加到GenBank的特征数: {len(features)}")
        
        return features
    
    def _parse_gff_attributes(self, attributes):
        """解析GFF3属性字符串"""
        attr_dict = {}
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key] = value
        return attr_dict
    
    def _create_seq_feature(self, feature_type, relative_start, relative_end, 
                           strand, attr_dict, transcript_map):
        """创建SeqFeature对象"""
        # 确定转录本名称
        transcript_name = "Unknown"
        parent_attr = attr_dict.get('Parent', '')
        
        if parent_attr:
            # Parent可能包含多个ID，用逗号分隔
            parent_ids = parent_attr.split(',')
            for parent_id in parent_ids:
                if parent_id in transcript_map:
                    transcript_name = transcript_map[parent_id]
                    break
            else:
                # 如果在映射中没找到，使用第一个Parent ID
                transcript_name = parent_ids[0].split('.')[0] if '.' in parent_ids[0] else parent_ids[0]
        
        # 创建标签：转录本-特征类型
        feature_type_label = self.target_feature_types[feature_type]
        label = f"{transcript_name}-{feature_type_label}"
        
        # 处理链方向
        feature_strand = 1 if strand == '+' else -1 if strand == '-' else None
        
        # 创建SeqFeature
        qualifiers = {
            'label': [label]
        }
        
        # 添加其他有用的attributes
        if 'ID' in attr_dict:
            qualifiers['ID'] = [attr_dict['ID']]
        if 'Name' in attr_dict:
            qualifiers['Name'] = [attr_dict['Name']]
        if 'Parent' in attr_dict:
            qualifiers['Parent'] = [attr_dict['Parent']]
        
        return SeqFeature(
            FeatureLocation(relative_start, relative_end, feature_strand),
            type=feature_type,
            qualifiers=qualifiers
        )
