# FootprintViewer 模块

FootprintViewer 是 dna_features_viewer 的一个扩展模块，专门用于可视化 DNA FootPrint 数据。

## 功能特性

- **基因注释可视化**: 从基因组FASTA和GFF3文件自动创建基因注释图
- **FootPrint热图**: 可视化不同大小的FootPrint分数数据
- **智能布局**: 自动排列重叠的转录本，避免视觉冲突
- **多组织比较**: 支持同时比较多个组织的FootPrint数据
- **高亮区域**: 可以标记感兴趣的特定区域
- **高质量输出**: 支持PDF等矢量格式输出

## 模块结构

```
FootprintViewer/
├── __init__.py           # 模块入口
├── genbank_creator.py    # GenBank文件创建
├── data_processor.py     # FootPrint数据处理
├── visualizer.py         # 可视化工具
├── api.py               # 主要API接口
└── README.md            # 说明文档
```

## 快速开始

### 基本导入

```python
from dna_features_viewer.FootprintViewer import (
    plot_region_with_footprints,
    plot_multi_tissue_comparison
)
```

### 单组织分析

```python
# 绘制单个组织的FootPrint分析图
fig = plot_region_with_footprints(
    genome_fasta="path/to/genome.fa",
    gff3_file="path/to/annotations.gff3", 
    fp_score_file="path/to/tissue_scores.parquet",
    chrom="Chr1",
    start=1000,
    end=5000,
    highlight_regions=[(2000, 2500)],  # 可选：高亮区域
    output_file="single_tissue_analysis.pdf"
)
```

### 多组织比较

```python
# 比较多个组织的FootPrint数据
fp_files = {
    "leaf": "path/to/leaf_scores.parquet",
    "root": "path/to/root_scores.parquet", 
    "stem": "path/to/stem_scores.parquet"
}

fig = plot_multi_tissue_comparison(
    genome_fasta="path/to/genome.fa",
    gff3_file="path/to/annotations.gff3",
    fp_files_dict=fp_files,
    chrom="Chr1", 
    start=1000,
    end=5000,
    highlight_regions=[(2000, 2500), (3000, 3500)],
    output_file="multi_tissue_comparison.pdf"
)
```

## 高级使用

### 使用独立组件

```python
from dna_features_viewer.FootprintViewer import (
    GenBankCreator,
    FootprintDataProcessor,
    FootprintVisualizer
)

# 创建GenBank文件
gb_creator = GenBankCreator()
gb_file = gb_creator.create_from_region(
    genome_fasta, gff3_file, "Chr1", 1000, 5000
)

# 处理FootPrint数据
processor = FootprintDataProcessor()
fp_data = processor.load_footprint_scores(fp_file, "Chr1", 1000, 5000)
heatmap_data, positions, radii = processor.create_heatmap_data(fp_data, 1000, 5000)

# 自定义可视化
visualizer = FootprintVisualizer()
# ... 自定义绘图逻辑
```

## 参数说明

### 主要函数参数

- `genome_fasta`: 基因组FASTA文件路径
- `gff3_file`: GFF3注释文件路径  
- `fp_score_file` / `fp_files_dict`: FootPrint分数文件（Parquet格式）
- `chrom`: 染色体名称（如 "Chr1"）
- `start`, `end`: 区域起止位置（1-based坐标）
- `figsize`: 图形大小，为None时自动调整
- `highlight_regions`: 高亮区域列表 `[(start1, end1), (start2, end2), ...]`
- `output_file`: 输出文件路径（支持PDF、PNG、SVG等格式）
- `genbank_file`: GenBank文件保存路径（可选，用于调试）

## 数据格式要求

### FootPrint分数文件格式

Parquet文件应包含以下列：
- `chrom`: 染色体名称
- `pos`: 基因组位置
- `radius`: FootPrint大小/半径
- `score`: FootPrint分数值

### GFF3文件格式

支持标准的GFF3格式，主要处理以下特征类型：
- `five_prime_UTR` / `5UTR`
- `CDS`
- `three_prime_UTR` / `3UTR`
- `mRNA` / `transcript`

## 输出示例

生成的图表包含：
1. **基因注释轨道**: 显示基因结构（UTR、CDS等）
2. **FootPrint热图轨道**: 以热图形式显示不同大小FootPrint的分数
3. **颜色条**: 显示FootPrint分数的数值范围
4. **高亮区域**: 用红色半透明区域标记感兴趣的位置

## 依赖库

- `pandas`: 数据处理
- `numpy`: 数值计算
- `matplotlib`: 绘图
- `biopython`: 生物信息学数据处理
- `dna_features_viewer`: 基因注释可视化（父包）

## 注意事项

1. 确保所有输入文件路径正确且文件存在
2. FootPrint分数文件必须是Parquet格式
3. 坐标系统使用1-based（与GFF3一致）
4. 图表会自动调整大小，但可以通过`figsize`参数自定义
5. 输出的PDF文件适合用于科研出版

## 错误处理

模块包含完善的错误处理机制：
- 文件不存在时会抛出清晰的错误信息
- 空数据时会生成空白热图并给出警告
- 解析错误时会显示详细的错误堆栈

## 版本历史

- v1.0.0: 初始版本，基于FP_heatmap.ipynb重构
