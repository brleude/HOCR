import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import ticker
import math
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import load_scores
from statistical_index import select_top_proteins
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


# 绘制各方法的Jackknife 方法关键蛋白质累加数量曲线
def cumulative_performance_visualization(all_rankings, key_proteins, max_range):
    # 定义排名范围
    ranking_ranges = list(range(0, max_range + 1, 1))

    # 累积计算函数
    def compute_cumulative_tp(sorted_proteins, key_proteins):
        """支持动态间隔的累积计算"""
        cumulative = []
        count = 0
        current_index = 0  # 当前处理的蛋白质索引

        # 预先生成有效排序列表（处理不足max_range的情况）
        valid_sorted = sorted_proteins[:max_range] if len(sorted_proteins) >= max_range else sorted_proteins

        for target_rank in ranking_ranges:
            # 处理超出实际数据范围的情况
            while current_index < target_rank:
                if current_index < len(valid_sorted):
                    protein = valid_sorted[current_index]
                    if protein in key_proteins:
                        count += 1
                current_index += 1
            cumulative.append(count)
        return cumulative

    # 生成足够长的排序列表
    method_rankings = {
        method: select_top_proteins(scores, max_range)
        for method, scores in all_rankings.items()
    }

    # 计算累积结果
    cumulative_results = {}
    for method, ranking in method_rankings.items():
        cumulative_results[method] = compute_cumulative_tp(ranking, key_proteins)

    # 优化可视化设置
    plt.figure(figsize=(11.69, 7), dpi=300)  # 修改图片大小和分辨率

    # 颜色方案
    color_palette = {
        'DC': '#84ba42',
        'BC': '#4485c7',
        'CC': '#7abbdb',
        'EC': '#dbb428',
        'PageRank': '#d4562e',
        'CDR': '#a51c36',
        'CSR': '#682487'
    }

    # 绘制曲线
    for method in method_rankings.keys():
        plt.plot(
            ranking_ranges,
            cumulative_results[method],
            markersize=4,
            linewidth=1.4,  # 修改线条宽度
            color=color_palette[method],
            label=method
        )

    # 优化坐标轴
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))  # 主刻度每50单位
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))  # 次要刻度每10单位
    ax.yaxis.set_major_locator(ticker.MaxNLocator(10))  # Y轴自动刻度

    # 添加辅助网格
    # ax.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.4)  # 修改网格线宽度
    # ax.grid(which='minor', linestyle=':', alpha=0.3, linewidth=1.4)  # 修改网格线宽度

    # 调整标签
    plt.xlabel("Top Ranges", fontsize=20, labelpad=8)  # 修改x轴标签字体大小
    plt.ylabel("Cumulative Key Proteins", fontsize=20, labelpad=8)  #修改y轴标签字体大小
    # plt.title("Key Protein Prediction Performance", 
    #          fontsize=14, pad=15, fontweight='bold')

    # 设置刻度字体大小
    ax.tick_params(axis='x', labelsize=15)  # 修改x轴刻度字体大小
    ax.tick_params(axis='y', labelsize=15)  # 修改y轴刻度字体大小

    # 去掉右方和上方的框
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # 紧凑图例
    legend = plt.legend(
        loc='lower right',
        ncol=1,  # 一列布局
        fontsize=15,  # 修改图例字体大小
        framealpha=0.9
    )
    legend.get_frame().set_linewidth(1.4)  # 修改图例边框线宽度

    plt.tight_layout()
    plt.savefig('酿酒酵母PPI网络Jackknife关键蛋白质累加数量曲线.png', dpi=300)  # 保存图片
    plt.savefig('酿酒酵母PPI网络Jackknife关键蛋白质累加数量曲线.pdf', dpi=300)  # 保存图片
    plt.show()


# 准备数据
def read_key_and_no_key_proteins(file_path):
    """
    读取关键蛋白和非关键蛋白数据
    :param file_path: xls 文件路径
    :return: 关键蛋白列表，非关键蛋白列表
    """
    # 读取关键蛋白
    key_proteins_df = pd.read_excel(file_path, sheet_name=0, header=None)
    key_proteins = key_proteins_df[0].tolist()  # 第一列为关键蛋白
    
    # 读取非关键蛋白
    non_key_proteins_df = pd.read_excel(file_path, sheet_name=1, header=None)
    non_key_proteins = non_key_proteins_df[0].tolist()  # 第一列为非关键蛋白
    
    return key_proteins, non_key_proteins
# 关键蛋白质和非关键蛋白质文件路径
current_file = Path(__file__).resolve()
parent3 = current_file.parent.parent.parent
key_nonkey_file_path = parent3/'data/sc/Key and Non-key Proteins.xls'
# 读取关键蛋白质和非关键蛋白质
key_proteins, non_key_proteins = read_key_and_no_key_proteins(key_nonkey_file_path)
# 读取网络数据
def load_ppi_edges(ppi_file):
    """加载PPI边列表"""
    ppi_edges = []
    with open(ppi_file, 'r') as f:
        for line in f:
            u, v = line.strip().split()
            ppi_edges.append((u, v))
    return ppi_edges
# PPI网络文件路径
ppi_file = parent3/'data/sc/PPI Network.txt'
# 加载PPI网络
ppi_edges = load_ppi_edges(ppi_file)
G = nx.Graph()
G.add_edges_from(ppi_edges)

# 读取各方法得分
dc_scores = load_scores('results/dc_scores_sc.txt')
cc_scores = load_scores('results/cc_scores_sc.txt')
bc_scores = load_scores('results/bc_scores_sc.txt')
ec_scores = load_scores('results/ec_scores_sc.txt')
pr_scores = load_scores('results/pr_scores_sc.txt')
cdr_scores = load_scores('results/cdr_scores_sc.txt')
csr_scores = load_scores('results/csr_scores_sc.txt')

# 整合得分
all_rankings = {
    'DC': dc_scores,
    'BC': bc_scores,
    'CC': cc_scores,
    'EC': ec_scores,
    'PageRank': pr_scores,
    'CDR': cdr_scores,
    'CSR': csr_scores,
}
max_range = math.floor(len(G.nodes)/5) # 进行比较的排名范围上界
# 绘制Jackknife 方法关键蛋白质累加数量曲线
cumulative_performance_visualization(all_rankings, key_proteins, max_range)