"""
FootPrint数据处理模块

处理footprint分数数据，包括数据加载和热图矩阵转换
"""

import pandas as pd
import numpy as np


class FootprintDataProcessor:
    """FootPrint数据处理器"""
    
    def __init__(self, default_radius_range=(2, 100)):
        """
        初始化数据处理器
        
        参数:
        - default_radius_range: 默认的footprint半径范围 (min, max)
        """
        self.default_radius_range = default_radius_range
    
    def load_footprint_scores(self, fp_score_file, chrom, start, end):
        """
        从footprint分数文件中加载指定区域的数据
        
        参数:
        - fp_score_file: footprint分数parquet文件路径
        - chrom: 染色体名称
        - start: 起始位置 (1-based)  
        - end: 结束位置 (1-based)
        
        返回:
        - DataFrame: 包含该区域的footprint分数数据
        """
        
        # 读取parquet文件
        try:
            df = pd.read_parquet(fp_score_file, engine="fastparquet")
        except ImportError:
            # 如果fastparquet不可用，尝试使用pyarrow
            try:
                df = pd.read_parquet(fp_score_file, engine="pyarrow")
            except ImportError:
                # 如果两个都不可用，使用默认引擎
                df = pd.read_parquet(fp_score_file)
        
        # 检查数据格式
        # 新格式：chrom和pos为列，r2到r100为footprint分数列
        if 'chrom' in df.columns and 'pos' in df.columns and 'r2' in df.columns:
            # 新格式：宽格式数据
            # 筛选指定区域的数据
            region_data = df[
                (df['chrom'] == chrom) & 
                (df['pos'] >= start) & 
                (df['pos'] <= end)
            ].copy()
        else:
            # 旧格式：长格式数据（chrom, pos, radius, score）
            # 筛选指定区域的数据
            region_data = df[
                (df['chrom'] == chrom) & 
                (df['pos'] >= start) & 
                (df['pos'] <= end)
            ].copy()
        
        return region_data
    
    def create_heatmap_data(self, fp_data, start, end, radius_range=None):
        """
        将footprint数据转换为热图矩阵格式
        
        参数:
        - fp_data: footprint数据DataFrame
        - start: 区域起始位置
        - end: 区域结束位置  
        - radius_range: footprint半径范围 (min, max)，None时使用默认值
        
        返回:
        - heatmap_data: 热图数据矩阵 (radius数量 x 位置数量)
        - positions: 位置数组
        - radii: 半径数组
        """
        
        if radius_range is None:
            radius_range = self.default_radius_range
            
        region_length = end - start + 1
        positions = np.arange(start, end + 1)
        radii = np.arange(radius_range[0], radius_range[1] + 1)
        
        # 初始化热图数据矩阵
        heatmap_data = np.zeros((len(radii), region_length))
        
        # 检查数据格式并填充数据
        if 'r2' in fp_data.columns:
            # 新格式：宽格式数据，每个半径是一个列（r2, r3, ..., r100）
            for _, row in fp_data.iterrows():
                pos_idx = int(row['pos']) - start  # 转换为相对位置索引
                
                if 0 <= pos_idx < region_length:
                    # 遍历所有半径列
                    for radius in radii:
                        col_name = f'r{radius}'
                        if col_name in fp_data.columns:
                            radius_idx = radius - radius_range[0]  # 转换为半径索引
                            if 0 <= radius_idx < len(radii):
                                score = row[col_name]
                                # 处理可能的NaN值
                                if pd.notna(score):
                                    heatmap_data[radius_idx, pos_idx] = score
        else:
            # 旧格式：长格式数据（chrom, pos, radius, score）
            for _, row in fp_data.iterrows():
                pos_idx = int(row['pos']) - start  # 转换为相对位置索引
                radius_idx = int(row['radius']) - radius_range[0]  # 转换为半径索引
                
                if 0 <= pos_idx < region_length and 0 <= radius_idx < len(radii):
                    heatmap_data[radius_idx, pos_idx] = row['score']
        
        return heatmap_data, positions, radii
    
    def get_data_statistics(self, heatmap_data):
        """
        获取热图数据的统计信息
        
        参数:
        - heatmap_data: 热图数据矩阵
        
        返回:
        - dict: 包含统计信息的字典
        """
        non_zero_data = heatmap_data[heatmap_data > 0]
        
        stats = {
            'total_points': heatmap_data.size,
            'non_zero_points': len(non_zero_data),
            'zero_ratio': (heatmap_data.size - len(non_zero_data)) / heatmap_data.size,
            'min_value': float(np.min(heatmap_data)),
            'max_value': float(np.max(heatmap_data)),
            'mean_value': float(np.mean(non_zero_data)) if len(non_zero_data) > 0 else 0.0,
            'std_value': float(np.std(non_zero_data)) if len(non_zero_data) > 0 else 0.0
        }
        
        return stats
