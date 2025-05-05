import numpy as np
import json
import csv
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import ticker
from typing import Dict, List
from matplotlib import gridspec
from matplotlib.lines import Line2D
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent.parent)
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
    plt.figure(figsize=(11.69, 21), dpi=300)  # 修改图片大小和分辨率
    fig = plt.gcf()
    gs = gridspec.GridSpec(4, 2, figure=fig)

    # 方法名称和对应的颜色及位置
    methods = {
        'DC': ('#84ba42', 0, 0),
        'BC': ('#4485c7', 0, 1),
        'CC': ('#7abbdb', 1, 0),
        'EC': ('#dbb428', 1, 1),
        'PageRank': ('#d4562e', 2, 0),
        'CDR_5': ('#a51c36', 2, 1),
        'CSR_3': ('#682487', 3, 0)
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
            line, = ax.plot(x, y, color=color, linewidth=2.4)  # 修改线条宽度
            ax.fill_between(x, y, alpha=0.2, color=color)
            method_lines.append(line)
        line = ax.axvline(x=avg_rank, color='red', linestyle='--')
        # ax.set_title(f'{method_name}', fontsize=20)  # 标题设置

        # 根据位置设置x轴和y轴
        if col == 0 and row in [0, 1, 2]:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        elif col == 0 and row == 3:
            ax.set_xlabel('Rank', fontsize=40)
        elif col == 1:
            ax.set_yticklabels([])
            ax.set_ylabel('')
        if col == 1 and row in [0, 1]:
            ax.set_xticklabels([])
            ax.set_xlabel('')
        elif col == 1 and row in [2]:
            ax.set_xlabel('Rank', fontsize=40)
        else:
            ax.set_ylabel('# of Key\n Proteins', fontsize=40)  # 修改y轴标签字体大小

        ax.set_ylim(0, y_max)  # 统一纵坐标
        # ax.grid(True, alpha=0.3)
        ax.tick_params(axis='x', labelsize=36)  # 修改x轴刻度字体大小
        ax.tick_params(axis='y', labelsize=36)  # 修改y轴刻度字体大小

        ax = plt.gca()
        ax.xaxis.set_major_locator(ticker.MultipleLocator(300))  # 主刻度每50单位
        # ax.xaxis.set_minor_locator(ticker.MultipleLocator(200))  # 次要刻度每10单位
        ax.yaxis.set_major_locator(ticker.MaxNLocator(3))  # Y轴自动刻度

        # 设置边框线宽
        for spine in ax.spines.values():
            spine.set_linewidth(1.4)
        ax.spines['right'].set_visible(False)  # 去掉右方的框
        ax.spines['top'].set_visible(False)  # 去掉上方的框

        # 只收集一次平均排名线条
        if avg_line is None:
            avg_line = line

    # 生成图例标签
    legend_lines = [avg_line]
    legend_labels = ['(Average Ranking\nof Key Proteins)']
    for method_name, avg_rank in avg_rankings.items():
        legend_labels.append(f'{method_name}: {avg_rank:.2f}')

    # 创建代理艺术家
    proxy_lines = [Line2D([0], [0], linestyle="", marker="") for _ in range(len(legend_labels) - 1)]
    all_lines = legend_lines + method_lines + proxy_lines

    # 在右下角（1, 3）处绘制统一的图例
    ax = fig.add_subplot(gs[3, 1])
    ax.axis('off')  # 隐藏坐标轴
    ax.legend(all_lines, legend_labels, fontsize=20, loc='lower right', title_fontsize=15)

    plt.tight_layout()
    plt.savefig('大肠杆菌PPI网络关键蛋白质排名分布.png', dpi=300)  # 保存图片
    plt.savefig('大肠杆菌PPI网络关键蛋白质排名分布.pdf', dpi=300)  # 保存图片
    plt.show()


# 提取原始PPI网络
current_file = Path(__file__).resolve()
parent4 = current_file.parent.parent.parent.parent
ec_ppi_path = parent4/'data/ec/Ec PPI Network.csv'
ppi_edges_ec = []
with open(ec_ppi_path, 'r', encoding='utf-8') as f:
    reader = csv.reader(f)
    next(reader)  # 跳过标题行
    for row in reader:
        ppi_edges_ec.append( (row[0], row[1]) )

# 提取关键蛋白质
def extract_gene_names_from_file(filename):
    """
    从TSV文件中提取蛋白质名称
    :param filename: TSV文件路径
    :return: 包含所有蛋白质名称的列表，包含拆分后的双名称
    """
    gene_names = []
    
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            # 读取标题行
            headers = f.readline().strip().split('\t')
            try:
                attr_index = headers.index('Attributes')
            except ValueError:
                return []  # 无Attributes列直接返回空列表

            # 逐行处理数据
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) <= attr_index:
                    continue

                # 解析JSON
                try:
                    attr = json.loads(parts[attr_index])
                    if 'GeneName' in attr:
                        names = attr['GeneName'].split('/')
                        gene_names.extend([name.strip() for name in names if name.strip()])
                except json.JSONDecodeError:
                    continue  # 跳过无效JSON

    except FileNotFoundError:
        print(f"错误：文件 {filename} 未找到")
        return []
    
    return list(set(gene_names))
ec_key_genes = extract_gene_names_from_file(parent4/"data/ec/Ec Key Proteins.tsv")

# 对网络中蛋白质进行分类
def get_key_protein_sets(ppi_edges, key_genes):
    """
    根据PPI网络和关键基因列表生成关键蛋白和非关键蛋白集合
    
    参数：
    ppi_edges -- 来自get_ppi_edges的边列表，格式[(name1,name2),...]
    key_genes -- 来自extract_gene_names_from_file的基因名称列表
    
    返回：
    (key_proteins, non_key_proteins) 两个集合
    """
    # 提取PPI网络中所有节点（自动去重）
    ppi_nodes = set()
    for u, v in ppi_edges:
        ppi_nodes.add(u)
        ppi_nodes.add(v)
    
    # 转换关键基因为集合
    key_gene_set = set(key_genes)
    
    # 关键蛋白：在PPI网络中的关键基因
    key_proteins = key_gene_set & ppi_nodes
    
    # 非关键蛋白：网络中存在但不在关键集合中
    non_key_proteins = ppi_nodes - key_proteins
    
    return list(key_proteins), list(non_key_proteins)
key_proteins_ec, non_key_proteins_ec = get_key_protein_sets(ppi_edges_ec, ec_key_genes)

# 构建大肠杆菌PPI网络
G_ec = nx.Graph()
G_ec.add_edges_from(ppi_edges_ec)

# 加载待比较得分
dc_scores_ec = load_scores('results/dc_scores_ec.txt')
bc_scores_ec = load_scores('results/bc_scores_ec.txt')
cc_scores_ec = load_scores('results/cc_scores_ec.txt')
ec_scores_ec = load_scores('results/ec_scores_ec.txt')
pr_scores_ec = load_scores('results/pr_scores_ec.txt')
cdr_scores_ec = load_scores('results/cdr_scores_ec.txt')
csr_scores_ec = load_scores('results/csr_scores_ec.txt')

# 整合所有方法的得分
all_rankings_ec = {
    'DC': dc_scores_ec,
    'BC': bc_scores_ec,
    'CC': cc_scores_ec,
    'EC': ec_scores_ec,
    'PageRank': pr_scores_ec,
    'CDR_5': cdr_scores_ec,
    'CSR_3': csr_scores_ec,
}

# 绘制关键蛋白质排名分布图
plot_multiple_ranking_distribution(all_rankings_ec, key_proteins_ec, bin_size=50)