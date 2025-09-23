"""
FootprintViewer - DNA FootPrint 数据可视化工具

这个模块提供了用于可视化DNA footprint数据的工具,包括:
- 从基因组FASTA和GFF3文件创建GenBank记录
- 加载和处理footprint分数数据 
- 绘制基因注释图和footprint热图
- 支持多组织数据比较

主要类和函数:
- FootprintDataProcessor: footprint数据处理
- GenBankCreator: GenBank文件创建
- FootprintVisualizer: 可视化工具
- plot_region_with_footprints: 单组织分析
- plot_multi_tissue_comparison: 多组织比较
"""

from .data_processor import FootprintDataProcessor
from .genbank_creator import GenBankCreator
from .visualizer import FootprintVisualizer
from .api import plot_region_with_footprints, plot_multi_tissue_comparison

__version__ = "1.0.0"

__all__ = [
    'FootprintDataProcessor',
    'GenBankCreator', 
    'FootprintVisualizer',
    'plot_region_with_footprints',
    'plot_multi_tissue_comparison'
]
