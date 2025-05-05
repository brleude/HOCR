import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List
from matplotlib import gridspec
from matplotlib.lines import Line2D
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import load_scores
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 绘制各方法的关键蛋白质排名分布
def plot_multiple_ranking_distribution(
    rankings: Dict[str, Dict[str, float]], 
    key_proteins: List[str], 
    bin_size: int = 100):
    """
    比较多种方法的排名分布
    
    参数:
        rankings: 包含所有方法得分的字典，格式为 {'方法名': {蛋白质: 得分}}
        key_proteins: 关键蛋白列表
        bin_size: 排名区间大小
    """
    # 设置图表布局
    plt.figure(figsize=(11.69, 7), dpi=300)  # 修改图片大小和分辨率
    fig = plt.gcf()
    gs = gridspec.GridSpec(2, 4, figure=fig)

    # 方法名称和对应的颜色及位置
    methods = {
        'DC': ('#84ba42', 0, 0),
        'BC': ('#4485c7', 0, 1),
        'CC': ('#7abbdb', 0, 2),
        'EC': ('#dbb428', 0, 3),
        'PageRank': ('#d4562e', 1, 0),
        'CDR': ('#a51c36', 1, 1),
        'CSR': ('#682487', 1, 2)
    }

    def process_rankings(scores):
        sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        ranks = [i + 1 for i, (p, _) in enumerate(sorted_scores) if p in key_proteins]
        if not ranks:
            return [], np.array([]), np.array([])
        bins = range(0, max(ranks) + bin_size, bin_size)
        hist, edges = np.histogram(ranks, bins=bins)
        x = (edges[:-1] + edges[1:]) / 2
        window_size = 3
        if len(hist) >= window_size:
            y = np.convolve(hist, np.ones(window_size) / window_size, mode='valid')
            x = x[window_size - 1:]
        else:
            y = hist
            x = x[:len(y)]
        return ranks, x, y

    # 预处理数据并计算全局y轴最大值
    processed_data = {}
    max_y = 0
    for method_name in methods:
        scores = rankings.get(method_name, {})
        ranks, x, y = process_rankings(scores)
        processed_data[method_name] = (ranks, x, y)
        if len(y) > 0:
            current_max = np.max(y)
            max_y = current_max if current_max > max_y else max_y

    # 统一设置纵坐标范围
    y_max = max_y * 1.1 if max_y > 0 else 1

    # 存储平均排名
    avg_rankings = {}

    # 用于收集平均排名线条
    avg_line = None
    # 用于收集各方法曲线线条
    method_lines = []

    # 绘制子图
    for method_name, (color, row, col) in methods.items():
        ranks, x, y = processed_data[method_name]
        avg_rank = np.mean(ranks) if len(ranks) > 0 else 0
        avg_rankings[method_name] = avg_rank

        ax = fig.add_subplot(gs[row, col])
        line = None
        if len(x) > 0 and len(y) > 0:
            line, = ax.plot(x, y, color=color, linewidth=1.4)  # 修改线条宽度
            ax.fill_between(x, y, alpha=0.2, color=color)
            method_lines.append(line)
        line = ax.axvline(x=avg_rank, color='red', linestyle='--')
        # ax.set_title(f'{method_name}', fontsize=20)

        # 根据位置设置x轴和y轴
        if row == 0 and col in [0, 1, 2]:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        elif row == 0 and col == 3:
            ax.set_xlabel('Rank', fontsize=20)  # 修改x轴标签字体大小
        elif row == 1:
            ax.set_xlabel('Rank', fontsize=20)  # 修改x轴标签字体大小
        if row == 0 and col in [1, 2, 3]:
            ax.set_yticklabels([])
            ax.set_ylabel('')
        elif row == 1 and col in [1, 2]:
            ax.set_yticklabels([])
            ax.set_ylabel('')
        else:
            ax.set_ylabel('# of Key Proteins', fontsize=20)  # 修改y轴标签字体大小

        ax.set_ylim(0, y_max)  # 统一纵坐标
        # ax.grid(True, alpha=0.3)
        ax.tick_params(axis='x', labelsize=15)  # 修改x轴刻度字体大小
        ax.tick_params(axis='y', labelsize=15)  # 修改y轴刻度字体大小

        # 设置边框线宽
        for spine in ax.spines.values():
            spine.set_linewidth(1.4)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False) 

        # 只收集一次平均排名线条
        if avg_line is None:
            avg_line = line

    # 生成图例标签
    legend_lines = [avg_line]
    legend_labels = ['Average Ranking\nof Key Proteins']
    for method_name, avg_rank in avg_rankings.items():
        legend_labels.append(f'{method_name}: {avg_rank:.2f}')

    # 创建代理艺术家
    proxy_lines = [Line2D([0], [0], linestyle="", marker="") for _ in range(len(legend_labels) - 1)]
    all_lines = legend_lines + method_lines + proxy_lines

    # 在右下角（1, 3）处绘制统一的图例
    ax = fig.add_subplot(gs[1, 3])
    ax.axis('off')  # 隐藏坐标轴
    ax.legend(all_lines, legend_labels, loc='lower right', fontsize=13)  # 修改图例字体大小

    plt.tight_layout()
    plt.savefig('酿酒酵母PPI网络关键蛋白质排名分布.png', dpi=300)
    plt.savefig('酿酒酵母PPI网络关键蛋白质排名分布.pdf', dpi=300)
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
# 读取各方法得分
dc_scores = load_scores('results/dc_scores_sc.txt')
cc_scores = load_scores('results/cc_scores_sc.txt')
bc_scores = load_scores('results/bc_scores_sc.txt')
ec_scores = load_scores('results/ec_scores_sc.txt')
pr_scores = load_scores('results/pr_scores_sc.txt')
cdr_scores = load_scores('results/cdr_scores_sc.txt')
csr_scores = load_scores('results/csr_scores_sc.txt')
# 整合所有方法的得分
all_rankings = {
    'DC': dc_scores,
    'BC': bc_scores,
    'CC': cc_scores,
    'EC': ec_scores,
    'PageRank': pr_scores,
    'CDR': cdr_scores,
    'CSR': csr_scores,
}

# 绘制关键蛋白质排名分布图
plot_multiple_ranking_distribution(all_rankings, key_proteins, bin_size=100)