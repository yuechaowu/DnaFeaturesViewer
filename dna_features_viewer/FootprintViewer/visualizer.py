"""
FootPrint可视化模块

包含用于绘制基因注释和footprint热图的可视化工具
"""

import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord


class TypeOnlyTranslator(BiopythonTranslator):
    """自定义的转录本标签翻译器"""
    
    def compute_feature_label(self, feature):
        """
        计算特征标签，使用现有的label如果存在，否则生成默认标签
        """
        # 优先使用现有的label
        if 'label' in feature.qualifiers and feature.qualifiers['label']:
            return feature.qualifiers['label'][0]
        
        t = (feature.type or "").lower()
        
        # 基本类型映射
        type_map = {
            "five_prime_utr": "5UTR",
            "five_prime_untranslated_region": "5UTR",
            "5utr": "5UTR",
            "three_prime_utr": "3UTR", 
            "three_prime_untranslated_region": "3UTR",
            "3utr": "3UTR",
            "cds": "CDS",
            "exon": "exon",
            "mrna": "mRNA",
            "gene": "gene"
        }
        
        base_label = type_map.get(t, feature.type or "feature")
        
        # 尝试添加转录本信息
        qualifiers = feature.qualifiers
        transcript_info = ""
        
        # 优先级：Name > ID > Parent
        if 'Name' in qualifiers and qualifiers['Name']:
            transcript_info = f":{qualifiers['Name'][0]}"
        elif 'ID' in qualifiers and qualifiers['ID']:
            transcript_info = f":{qualifiers['ID'][0]}"
        elif 'Parent' in qualifiers and qualifiers['Parent']:
            parent_id = qualifiers['Parent'][0]
            # 简化Parent ID显示
            if '.' in parent_id:
                parent_id = parent_id.split('.')[-1]
            transcript_info = f":{parent_id}"
        
        return base_label + transcript_info
    
    def compute_feature_color(self, feature):
        """
        为不同类型的特征设置不同颜色
        """
        t = (feature.type or "").lower()
        
        color_map = {
            "gene": "#ff9999",           # 浅红色
            "mrna": "#99ff99",           # 浅绿色
            "cds": "#9999ff",            # 浅蓝色
            "exon": "#ffff99",           # 浅黄色
            "five_prime_utr": "#ff99ff", # 浅紫色
            "5utr": "#ff99ff",
            "three_prime_utr": "#99ffff", # 浅青色
            "3utr": "#99ffff",
            "utr": "#f0f0f0"             # 浅灰色
        }
        
        return color_map.get(t, "#cccccc")  # 默认浅灰色


class FootprintVisualizer:
    """FootPrint数据可视化器"""
    
    def __init__(self):
        """初始化可视化器"""
        # 设置matplotlib参数
        import matplotlib as mpl
        mpl.rcParams['svg.fonttype'] = 'none'
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        mpl.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
        mpl.rcParams['axes.unicode_minus'] = False
        
        # 转录本组件的颜色
        self.component_colors = {
            '5UTR': "#FFB6C1",   # 浅粉色
            'CDS': "#87CEEB",    # 天蓝色
            '3UTR': "#98FB98",   # 浅绿色
            'EXON': "#F0E68C",   # 卡其色
            'INTRON': "#D3D3D3", # 浅灰色
            'MOTIF': "#DDA0DD",  # 浅紫色
            'source': "#F5F5F5"  # 浅灰色
        }
    
    def _auto_figsize(self, num_tracks):
        """
        根据轨道数量自动调整画布大小
        
        参数:
        - num_tracks: 轨道数量（包括基因注释轨道和footprint热图轨道）
        
        返回:
        - (width, height): 画布大小元组
        """
        if num_tracks <= 2:  # 单个组织：1个基因注释轨道 + 1个热图轨道
            return (12, 4)
        elif num_tracks <= 4:  # 2-3个组织
            return (12, 6)
        else:  # 更多组织
            base_h = 1.5  # 基础高度（每个轨道）
            pad = 1.0     # 额外边距
            h = base_h * num_tracks + pad
            return (14, h)
    
    def create_transcript_layout_visualization(self, record):
        """
        根据label中'-'前的转录本名称分组，
        将不重叠的转录本放在同一行，重叠的转录本放在不同行
        """
        
        # 1. 解析所有特征，按转录本/特征分组
        transcripts = {}
        
        for feature in record.features:
            if 'label' in feature.qualifiers:
                label = feature.qualifiers['label'][0]
                if '-' in label:
                    transcript_name = label.split('-')[0]  # 取'-'前的部分作为转录本名称
                    component_type = label.split('-')[1]   # 取'-'后的部分作为组件类型
                    
                    if transcript_name not in transcripts:
                        transcripts[transcript_name] = []
                    
                    transcripts[transcript_name].append({
                        'feature': feature,
                        'component_type': component_type,
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'strand': feature.location.strand,
                        'label': label
                    })
                else:
                    # 将没有'-'的特征作为单独的"转录本"处理
                    transcript_name = label  # 直接使用label作为名称
                    transcripts[transcript_name] = [{
                        'feature': feature,
                        'component_type': feature.type,  # 使用特征类型作为组件类型
                        'start': int(feature.location.start),
                        'end': int(feature.location.end),
                        'strand': feature.location.strand,
                        'label': label
                    }]
        
        print("识别到的转录本:")
        for transcript_name, components in transcripts.items():
            print(f"  {transcript_name}: {len(components)}个组件")
            for comp in components:
                print(f"    - {comp['component_type']}: {comp['start']}-{comp['end']}")
        
        # 2. 计算转录本的整体范围
        transcript_ranges = {}
        for transcript_name, components in transcripts.items():
            starts = [comp['start'] for comp in components]
            ends = [comp['end'] for comp in components]
            transcript_ranges[transcript_name] = {
                'start': min(starts),
                'end': max(ends),
                'components': components
            }
        
        # 3. 智能布局算法：将不重叠的转录本放在同一行
        transcript_rows = self._arrange_transcript_rows(transcript_ranges)
        
        print(f"\n布局结果（共{len(transcript_rows)}行）:")
        for i, row in enumerate(transcript_rows):
            print(f"  第{i+1}行: {', '.join(row)}")
        
        return transcript_rows, transcript_ranges, transcripts
    
    def _arrange_transcript_rows(self, transcript_ranges):
        """智能布局算法：将不重叠的转录本放在同一行"""
        def transcripts_overlap(t1_range, t2_range):
            """检查两个转录本是否重叠"""
            return not (t1_range['end'] < t2_range['start'] or t2_range['end'] < t1_range['start'])
        
        # 贪心算法分配行
        transcript_rows = []
        transcript_names = list(transcript_ranges.keys())
        
        for transcript_name in transcript_names:
            transcript_range = transcript_ranges[transcript_name]
            
            # 尝试放入现有行
            placed = False
            for row_idx, row in enumerate(transcript_rows):
                # 检查与该行中所有转录本是否重叠
                overlaps = False
                for existing_transcript in row:
                    if transcripts_overlap(transcript_range, transcript_ranges[existing_transcript]):
                        overlaps = True
                        break
                
                if not overlaps:
                    row.append(transcript_name)
                    placed = True
                    break
            
            # 如果无法放入现有行，创建新行
            if not placed:
                transcript_rows.append([transcript_name])
        
        return transcript_rows
    
    def draw_transcript_row(self, ax, transcript_row, transcript_ranges, transcripts, record):
        """
        绘制一行转录本
        """
        # 为每个转录本分配独特的基础颜色
        transcript_names = list(transcript_ranges.keys())
        transcript_base_colors = plt.cm.Set3(np.linspace(0, 1, len(transcript_names)))
        transcript_color_map = dict(zip(transcript_names, transcript_base_colors))
        
        graphic_features = []
        
        for transcript_name in transcript_row:
            transcript_range = transcript_ranges[transcript_name]
            components = transcript_range['components']
            
            # 判断是否为多组件转录本
            if len(components) > 1:
                # 多组件转录本：先绘制整体框架，然后在框架内绘制组件
                base_color = transcript_color_map[transcript_name]
                
                # 绘制转录本整体框架（有填充的背景）
                background_color = (*base_color[:3], 0.3)  # 添加透明度
                graphic_features.append(GraphicFeature(
                    start=transcript_range['start'],
                    end=transcript_range['end'],
                    strand=components[0]['strand'],
                    color=background_color,  # 有填充
                    linewidth=2,
                    linecolor=base_color,
                    label=f"{transcript_name}",
                    thickness=20
                ))
                
                # 在框架内绘制各个组件（较细，位置稍微偏移以显示在框架内）
                for comp in components:
                    comp_color = self.component_colors.get(comp['component_type'], "#CCCCCC")
                    
                    graphic_features.append(GraphicFeature(
                        start=comp['start'],
                        end=comp['end'],
                        strand=comp['strand'],
                        color=comp_color,
                        label=f"{comp['component_type']}",
                        thickness=12,  # 较细的组件
                        linewidth=1
                    ))
            else:
                # 单组件特征：直接显示
                comp = components[0]
                comp_color = self.component_colors.get(comp['component_type'], "#CCCCCC")
                
                graphic_features.append(GraphicFeature(
                    start=comp['start'],
                    end=comp['end'],
                    strand=comp['strand'],
                    color=comp_color,
                    label=f"{transcript_name}",
                    thickness=15
                ))
        
        # 创建并绘制图形记录
        graphic_record = GraphicRecord(
            sequence_length=len(record.seq),
            features=graphic_features
        )
        
        graphic_record.plot(ax=ax, with_ruler=False)
        
        # 设置行标签
        row_labels = ', '.join(transcript_row)
        ax.set_ylabel(f"{row_labels}", fontsize=9)
        ax.set_xlim(0, len(record.seq))
        
        # 减小y轴标签的边距
        ax.tick_params(axis='y', labelsize=8)
        ax.yaxis.set_label_coords(-0.05, 0.5)
    
    def plot_heatmap(self, ax, heatmap_data, radii, seq_len, max_score, title=""):
        """
        绘制footprint分数热图
        
        参数:
        - ax: matplotlib轴对象
        - heatmap_data: 热图数据矩阵
        - radii: 半径数组
        - seq_len: 序列长度
        - max_score: 最大分数值
        - title: 图表标题
        """
        im = ax.imshow(
            heatmap_data,
            aspect="auto",
            origin="lower",
            extent=[0, seq_len, radii[0]-0.5, radii[-1]+0.5],
            interpolation="nearest",
            cmap=plt.cm.Blues,
            vmin=0, vmax=max_score
        )
        
        ax.set_xlim(0, seq_len)
        ax.set_ylabel("FootPrint Size (bp)")
        if title:
            ax.set_title(title)
        
        return im
    
    def add_highlight_regions(self, axes, highlight_regions, start, seq_len):
        """
        在图表中添加高亮区域
        
        参数:
        - axes: matplotlib轴对象列表
        - highlight_regions: 高亮区域列表 [(start1, end1), (start2, end2), ...]
        - start: 区域起始位置
        - seq_len: 序列长度
        """
        if highlight_regions:
            for region_start, region_end in highlight_regions:
                # 转换为相对坐标
                rel_start = max(0, region_start - start)
                rel_end = min(seq_len, region_end - start)
                
                # 在所有轴上添加高亮
                for ax in axes:
                    ax.axvspan(rel_start, rel_end, color="red", alpha=0.3, linewidth=0.5)
