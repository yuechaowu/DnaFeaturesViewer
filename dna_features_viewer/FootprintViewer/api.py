"""
FootPrint可视化主要API接口

提供高级的FootPrint数据可视化函数，整合所有模块功能
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

from .genbank_creator import GenBankCreator
from .data_processor import FootprintDataProcessor
from .visualizer import FootprintVisualizer, TypeOnlyTranslator


def plot_region_with_footprints(genome_fasta, gff3_file, fp_score_file, chrom, start, end, 
                               figsize=None, highlight_regions=None, output_file=None, genbank_file=None):
    """
    绘制指定区域的基因注释和footprint分数热图
    
    参数:
    - genome_fasta: 基因组FASTA文件路径
    - gff3_file: GFF3注释文件路径
    - fp_score_file: footprint分数parquet文件路径
    - chrom: 染色体名称
    - start: 起始位置 (1-based)
    - end: 结束位置 (1-based)
    - figsize: 图片大小 (如果为None则自动调整)
    - highlight_regions: 高亮区域列表 [(start1, end1), (start2, end2), ...]
    - output_file: 输出文件路径
    - genbank_file: GenBank文件保存路径 (如果指定则不会删除该文件)
    
    返回:
    - fig: matplotlib图形对象
    """
    
    # 初始化组件
    gb_creator = GenBankCreator()
    data_processor = FootprintDataProcessor()
    visualizer = FootprintVisualizer()
    
    # 创建GenBank文件
    print(f"正在为区域 {chrom}:{start}-{end} 创建GenBank文件...")
    gb_file = gb_creator.create_from_region(genome_fasta, gff3_file, chrom, start, end, genbank_file)
    is_temp_file = genbank_file is None
    
    # 加载footprint数据
    print("正在加载footprint分数数据...")
    fp_data = data_processor.load_footprint_scores(fp_score_file, chrom, start, end)
    
    if fp_data.empty:
        print(f"警告: 在区域 {chrom}:{start}-{end} 中未找到footprint数据")
        heatmap_data = np.zeros((99, end - start + 1))  # 创建空的热图数据
        radii = np.arange(2, 101)
        max_score = 1.0  # 默认最大值
    else:
        print(f"找到 {len(fp_data)} 个footprint数据点")
        heatmap_data, positions, radii = data_processor.create_heatmap_data(fp_data, start, end)
        raw_max = np.max(heatmap_data) if np.max(heatmap_data) > 0 else 1.0
        max_score = min(raw_max, 5.0)  # 限制最大值不超过5
        print(f"数据最大值: {raw_max:.3f}, colorbar最大值: {max_score:.3f}")
    
    # 读取GenBank记录
    record = SeqIO.read(gb_file, "genbank")
    
    # 自动调整画布大小（2个轨道：基因注释 + footprint热图）
    if figsize is None:
        figsize = visualizer._auto_figsize(2)
    
    # 创建图形
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=figsize, 
        gridspec_kw={"height_ratios": [3, 2]}
    )
    
    # 绘制基因注释图 - 使用智能布局
    print("正在绘制基因注释（使用智能布局）...")
    
    # 如果记录中有特征，使用智能布局
    if record.features:
        transcript_rows, transcript_ranges, transcripts = visualizer.create_transcript_layout_visualization(record)
        
        if len(transcript_rows) > 1:
            # 多行布局 - 重新创建子图
            plt.close(fig)
            
            # 重新计算子图数量：转录本行数 + 1个footprint热图
            total_rows = len(transcript_rows) + 1
            height_ratios = [1.2] * len(transcript_rows) + [2]  # 转录本行较矮，热图行较高
            
            if figsize is None:
                figsize = visualizer._auto_figsize(total_rows)
            
            fig, axes = plt.subplots(
                total_rows, 1, figsize=figsize,
                gridspec_kw={"height_ratios": height_ratios}
            )
            
            if total_rows == 1:
                axes = [axes]
            
            # 绘制每一行的转录本，隐藏所有转录本轨道的x轴标签
            region_len = end - start + 1
            for row_idx, transcript_row in enumerate(transcript_rows):
                visualizer.draw_transcript_row(axes[row_idx], transcript_row, transcript_ranges, transcripts, record, start, region_len=region_len)
                # 隐藏所有转录本轨道的x轴标签（只有最后的footprint热图轨道显示x轴标签）
                axes[row_idx].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            
            # footprint热图使用最后一个轴
            ax2 = axes[-1]
            
            # 不再设置标题
        else:
            # 单行布局 - 使用原有方式，但需要调整坐标系
            translator = TypeOnlyTranslator()
            graphic_record = translator.translate_record(record)
            graphic_record.first_index = start  # 设置起始索引为实际基因组位置
            graphic_record.plot(ax=ax1, with_ruler=False, draw_line=False, strand_in_label_threshold=4)
            # 为单行布局添加虚线连接
            visualizer._add_intron_connections_simple(ax1, record, start)
            ax1.set_xlim(start, end)
            # 隐藏基因注释轨道的x轴标签
            ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    else:
        # 没有特征，使用原有方式
        translator = TypeOnlyTranslator()
        graphic_record = translator.translate_record(record)
        graphic_record.first_index = start  # 设置起始索引为实际基因组位置
        graphic_record.plot(ax=ax1, with_ruler=False, draw_line=False, strand_in_label_threshold=4)
        # 为无特征情况添加虚线连接（如果有合适的特征）
        visualizer._add_intron_connections_simple(ax1, record, start)
        ax1.set_xlim(start, end)
        # 隐藏基因注释轨道的x轴标签
        ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
    # 绘制footprint分数热图
    print("正在绘制footprint分数热图...")
    region_len = end - start + 1
    
    im = visualizer.plot_heatmap(ax2, heatmap_data, radii, region_len, max_score, start)
    ax2.set_xlabel(f"{chrom} position (bp)")
    
    # 添加高亮区域
    axes_to_highlight = axes if 'axes' in locals() and len(axes) > 2 else [ax1, ax2]
    visualizer.add_highlight_regions(axes_to_highlight, highlight_regions, start, region_len)
    
    # 添加颜色条（右上角横置）
    cbar_ax = fig.add_axes([0.85, 0.9, 0.1, 0.03])  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    cbar.set_label("FootPrint Score", fontsize=10)
    cbar.ax.tick_params(labelsize=8)
    
    plt.tight_layout()
    
    # 保存图片
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"图片已保存到: {output_file}")
    
    # 清理临时文件（仅当未指定genbank_file时）
    if is_temp_file and os.path.exists(gb_file):
        os.unlink(gb_file)
        print("临时GenBank文件已清理")
    elif not is_temp_file:
        print(f"GenBank文件已保存到: {gb_file}")
    
    return fig


def plot_multi_tissue_comparison(genome_fasta, gff3_file, fp_files_dict, chrom, start, end,
                                figsize=None, highlight_regions=None, output_file=None, genbank_file=None):
    """
    比较多个组织的footprint分数
    
    参数:
    - genome_fasta: 基因组FASTA文件路径
    - gff3_file: GFF3注释文件路径
    - fp_files_dict: footprint分数文件字典 {"tissue_name": "file_path", ...}
    - chrom: 染色体名称
    - start: 起始位置 (1-based)
    - end: 结束位置 (1-based)
    - figsize: 图片大小 (如果为None则自动调整)
    - highlight_regions: 高亮区域列表
    - output_file: 输出文件路径
    - genbank_file: GenBank文件保存路径 (如果指定则不会删除该文件)
    
    返回:
    - fig: matplotlib图形对象
    """
    
    # 初始化组件
    gb_creator = GenBankCreator()
    data_processor = FootprintDataProcessor()
    visualizer = FootprintVisualizer()
    
    n_tissues = len(fp_files_dict)
    
    # 创建GenBank文件
    print(f"正在为区域 {chrom}:{start}-{end} 创建GenBank文件...")
    gb_file = gb_creator.create_from_region(genome_fasta, gff3_file, chrom, start, end, genbank_file)
    is_temp_file = genbank_file is None
    record = SeqIO.read(gb_file, "genbank")
    region_len = end - start + 1
    
    # 智能绘制基因注释图（与单组织保持一致风格）
    transcript_rows = []
    base_idx = 1  # 热图起始轴索引（默认基因注释占1行）
    
    if record.features:
        transcript_rows, transcript_ranges, transcripts = visualizer.create_transcript_layout_visualization(record)
        if len(transcript_rows) > 1:
            # 多行转录本布局：总行数 = 转录本行 + n个组织热图
            total_rows = len(transcript_rows) + n_tissues
            height_ratios = [1.2] * len(transcript_rows) + [2] * n_tissues
            if figsize is None:
                figsize = visualizer._auto_figsize(total_rows)
            fig, axes = plt.subplots(
                total_rows, 1, figsize=figsize,
                gridspec_kw={"height_ratios": height_ratios}
            )
            if total_rows == 1:
                axes = [axes]
            # 绘制每一行的转录本，隐藏x轴标签
            for row_idx, transcript_row in enumerate(transcript_rows):
                visualizer.draw_transcript_row(axes[row_idx], transcript_row, transcript_ranges, transcripts, record, start, region_len=region_len)
                axes[row_idx].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            base_idx = len(transcript_rows)
        else:
            # 单行转录本：1个基因注释轨道 + n个组织热图
            total_rows = 1 + n_tissues
            height_ratios = [3] + [2] * n_tissues
            if figsize is None:
                figsize = visualizer._auto_figsize(total_rows)
            fig, axes = plt.subplots(
                total_rows, 1, figsize=figsize,
                gridspec_kw={"height_ratios": height_ratios}
            )
            if total_rows == 1:
                axes = [axes]
            translator = TypeOnlyTranslator()
            graphic_record = translator.translate_record(record)
            graphic_record.first_index = start
            graphic_record.plot(ax=axes[0], with_ruler=False, draw_line=False, strand_in_label_threshold=4)
            visualizer._add_intron_connections_simple(axes[0], record, start)
            axes[0].set_xlim(start, end)
            axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            base_idx = 1
    else:
        # 无特征：退回简单单行基因注释
        total_rows = 1 + n_tissues
        height_ratios = [3] + [2] * n_tissues
        if figsize is None:
            figsize = visualizer._auto_figsize(total_rows)
        fig, axes = plt.subplots(
            total_rows, 1, figsize=figsize,
            gridspec_kw={"height_ratios": height_ratios}
        )
        if total_rows == 1:
            axes = [axes]
        translator = TypeOnlyTranslator()
        graphic_record = translator.translate_record(record)
        graphic_record.first_index = start
        graphic_record.plot(ax=axes[0], with_ruler=False, draw_line=False, strand_in_label_threshold=4)
        visualizer._add_intron_connections_simple(axes[0], record, start)
        axes[0].set_xlim(start, end)
        axes[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        base_idx = 1
    
    # 为每个组织绘制footprint热图
    tissue_names = list(fp_files_dict.keys())
    all_max_scores = []  # 收集所有组织的最大值
    tissue_data = []     # 存储所有组织的数据
    
    # 第一遍：加载所有数据并计算全局最大值
    for i, (tissue, fp_file) in enumerate(fp_files_dict.items()):
        print(f"正在处理 {tissue} 组织的footprint数据...")
        
        # 加载footprint数据
        fp_data = data_processor.load_footprint_scores(fp_file, chrom, start, end)
        
        if fp_data.empty:
            print(f"警告: {tissue} 组织在区域 {chrom}:{start}-{end} 中未找到footprint数据")
            heatmap_data = np.zeros((99, region_len))
            radii = np.arange(2, 101)
            max_score = 0
        else:
            print(f"{tissue}: 找到 {len(fp_data)} 个footprint数据点")
            heatmap_data, positions, radii = data_processor.create_heatmap_data(fp_data, start, end)
            raw_max = np.max(heatmap_data) if np.max(heatmap_data) > 0 else 0
            max_score = min(raw_max, 5.0)  # 限制每个组织的最大值不超过5
            print(f"{tissue} 数据最大值: {raw_max:.3f}")
        
        tissue_data.append((tissue, heatmap_data, radii, max_score))
        all_max_scores.append(max_score)
    
    # 计算全局最大值，限制不超过5
    raw_global_max = max(all_max_scores) if all_max_scores and max(all_max_scores) > 0 else 1.0
    global_max = min(raw_global_max, 5.0)
    print(f"所有组织最大值: {raw_global_max:.3f}, colorbar最大值: {global_max:.3f}")
    
    # 第二遍：绘制热图
    im = None
    for i, (tissue, heatmap_data, radii, _) in enumerate(tissue_data):
        ax = axes[i + base_idx]
        im = visualizer.plot_heatmap(ax, heatmap_data, radii, region_len, global_max, start, title=f"{tissue}")
        
        # 与单组织风格一致：y轴仅显示数值含义，将组织名放在左侧标题位置
        ax.set_ylabel("FootPrint Size (bp)")
        ax.set_title(tissue, loc='left', fontsize=10)
        
        # 隐藏除了最后一个子图外的x轴刻度和标签
        if i < n_tissues - 1:
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        else:
            # 只在最后一个子图显示x轴标签
            ax.set_xlabel(f"{chrom} position (bp)")
    
    # 添加高亮区域
    visualizer.add_highlight_regions(axes, highlight_regions, start, region_len)
    
    # 添加颜色条（右上角横置）
    if im is not None:
        # 与单组织保持一致的颜色条风格与位置
        cbar_ax = fig.add_axes([0.85, 0.9, 0.1, 0.03])
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar.set_label("FootPrint Score", fontsize=10)
        cbar.ax.tick_params(labelsize=8)
    
    plt.tight_layout()
    
    # 保存图片（避免在PDF后端中因bbox_inches='tight'导致的超大栅格区域问题）
    if output_file:
        plt.savefig(output_file, dpi=300)
        print(f"比较图已保存到: {output_file}")
    
    # 清理临时文件（仅当未指定genbank_file时）
    if is_temp_file and os.path.exists(gb_file):
        os.unlink(gb_file)
        print("临时GenBank文件已清理")
    elif not is_temp_file:
        print(f"GenBank文件已保存到: {gb_file}")
    
    return fig


def create_example_usage():
    """
    创建使用示例
    
    返回:
    - dict: 包含示例参数的字典
    """
    
    example_params = {
        "genome_fasta": "/mnt/data/wyc/project/FP2LLM/data/genome/arabidopsis.fa",
        "gff3_file": "/mnt/data/wyc/project/FP2LLM/data/genome/arabidopsis.gff3",
        "fp_files": {
            "inflorescence": "/mnt/data/wyc/project/FP2LLM/data/FP_score_outputs/inflorescence_exp8_FTScore.parquet",
            "leaf": "/mnt/data/wyc/project/FP2LLM/data/FP_score_outputs/leaf_exp3_FTScore.parquet", 
            "root": "/mnt/data/wyc/project/FP2LLM/data/FP_score_outputs/root_exp3_FTScore.parquet"
        },
        "example_region": {
            "chrom": "Chr1",
            "start": 3631,
            "end": 5899
        }
    }
    
    usage_doc = """
    # FootPrint 分析工具使用示例
    
    ## 单组织分析
    ```python
    from dna_features_viewer.FootprintViewer import plot_region_with_footprints
    
    fig = plot_region_with_footprints(
        genome_fasta="./data/genome/arabidopsis.fa",
        gff3_file="./data/genome/arabidopsis.gff3",
        fp_score_file="./data/FP_score_outputs/leaf_exp3_FTScore.parquet",
        chrom="Chr1",
        start=3631,
        end=5899,
        highlight_regions=[(4000, 4200)],
        output_file="single_tissue_analysis.pdf"
    )
    ```
    
    ## 多组织比较
    ```python
    from dna_features_viewer.FootprintViewer import plot_multi_tissue_comparison
    
    fp_files = {
        "leaf": "./data/FP_score_outputs/leaf_exp3_FTScore.parquet",
        "root": "./data/FP_score_outputs/root_exp3_FTScore.parquet",
        "inflorescence": "./data/FP_score_outputs/inflorescence_exp8_FTScore.parquet"
    }
    
    fig = plot_multi_tissue_comparison(
        genome_fasta="./data/genome/arabidopsis.fa",
        gff3_file="./data/genome/arabidopsis.gff3",
        fp_files_dict=fp_files,
        chrom="Chr1",
        start=3631,
        end=5899,
        highlight_regions=[(4000, 4200), (4800, 5000)],
        output_file="multi_tissue_comparison.pdf"
    )
    ```
    """
    
    return example_params, usage_doc
