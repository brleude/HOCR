import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 计算网络中关键节点和非关键节点的各大小模体参与次数
def compute_motif_participation(G, key_nodes, non_key_nodes=None, motifs=[3,4,5]):
    """
    计算网络中关键节点和非关键节点的模体参与情况
    
    参数:
        G: 网络图
        key_nodes: 关键节点列表
        non_key_nodes: 非关键节点列表(可选，如果为None则自动计算)
        motifs: 要计算的模体尺寸列表
        
    返回:
        (key_counts, non_key_counts): 两个字典，分别包含关键节点和非关键节点的各大小模体参与次数计数
    """
    all_motifs = sorted(set(motifs))
    max_k = max(all_motifs) if all_motifs else 0
    if max_k < 3:
        return {}, {}

    # 处理节点分类
    all_nodes = set(G.nodes())
    key_set = set(key_nodes) & all_nodes
    
    if non_key_nodes is None:
        # 如果未指定非关键节点，则自动计算
        non_key_set = all_nodes - key_set
    else:
        # 确保非关键节点不包含任何关键节点
        non_key_set = (set(non_key_nodes) & all_nodes) - key_set
    
    # 初始化计数结构
    counters = {
        node: {k:0 for k in all_motifs} 
        for node in key_set | non_key_set  # 只为关注的节点初始化计数器
    }

    # 分层处理clique
    prev_cliques = []
    for k in range(3, max_k+1):
        current_cliques = []
        
        if k == 3:
            # 优化3-clique生成
            for u in G:
                if u not in counters:  # 跳过不关注的节点
                    continue
                    
                neighbors = sorted([v for v in G[u] if v > u and v in counters])
                for i, v in enumerate(neighbors):
                    common = set(neighbors[i+1:]) & set(G[v])
                    for w in common:
                        if w > v and w in counters:
                            clique = (u, v, w)
                            current_cliques.append(clique)
                            if 3 in all_motifs:
                                for node in clique:
                                    counters[node][3] += 1
        else:
            # 流式处理大尺寸clique
            for clique in prev_cliques:
                last_node = clique[-1]
                common_neighbors = set(G[clique[0]])
                for node in clique[1:]:
                    common_neighbors &= set(G[node])
                
                for neighbor in common_neighbors:
                    if neighbor > last_node and neighbor in counters:
                        new_clique = (*clique, neighbor)
                        current_cliques.append(new_clique)
                        if k in all_motifs:
                            for node in new_clique:
                                counters[node][k] += 1
        
        # 仅保留必要数据
        prev_cliques = current_cliques if k+1 <= max_k else []
    
    # 分类统计结果
    key_counts = {
        node: {k: counters[node][k] for k in all_motifs} 
        for node in key_set
    }
    
    non_key_counts = {
        node: {k: counters[node][k] for k in all_motifs} 
        for node in non_key_set
    }
    
    return key_counts, non_key_counts

########################################## 准备数据 ##########################################
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
parent4 = current_file.parent.parent.parent.parent
key_nonkey_file_path = parent4/'data/sc/Key and Non-key Proteins.xls'
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
ppi_file = parent4/'data/sc/PPI Network.txt'
# 加载PPI网络
ppi_edges = load_ppi_edges(ppi_file)
G = nx.Graph()
G.add_edges_from(ppi_edges)

# 绘图数据
motif_sizes = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
key_counts, non_key_counts = compute_motif_participation(G, key_proteins, non_key_proteins, 
                                                         motifs=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12])

# 转换数据结构
def prepare_plot_data(count_dict, motif_sizes):
    plot_data = {k: [] for k in motif_sizes}
    for node_data in count_dict.values():
        for k in motif_sizes:
            plot_data[k].append(node_data[k])
    return plot_data

key_plot_data = prepare_plot_data(key_counts, motif_sizes)
non_key_plot_data = prepare_plot_data(non_key_counts, motif_sizes)

########################################## 绘图 ##########################################
# 创建图形
plt.figure(figsize=(11.69, 4), dpi=300)  #图片尺寸和分辨率

# 设置颜色
key_color = '#FF6B6B'
non_key_color = '#4D96FF'

# 创建自定义小提琴图绘制函数
def plot_side_violin(pos, data, color, side='left'):
    """
    side: 'left' or 'right'
    """
    # 创建完整的小提琴图
    violin = plt.violinplot(data, positions=[pos], widths=0.8, showmeans=False, showmedians=False)

    # 修改小提琴图形状为半边形
    for pc in violin['bodies']:
        path = pc.get_paths()[0]
        verts = path.vertices

        if side == 'left':
            # 只保留x <= pos的部分
            verts[:, 0] = np.clip(verts[:, 0], -np.inf, pos)
        else:
            # 只保留x >= pos的部分
            verts[:, 0] = np.clip(verts[:, 0], pos, np.inf)

        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.8)
        pc.set_linewidth(1.4000)  # 设置线条宽度

    # 添加均值线
    mean = np.mean(data)

    if side == 'left':
        plt.hlines(mean, pos - 0.4, pos, colors='black', linestyles='--', lw=1.4000, zorder=10)  # 设置线条宽度
    else:
        plt.hlines(mean, pos, pos + 0.4, colors='black', linestyles='--', lw=1.4000, zorder=10)  # 设置线条宽度

# 绘制组合小提琴图
for i, k in enumerate(motif_sizes):
    pos = i + 1
    # 左侧：关键节点
    plot_side_violin(pos, key_plot_data[k], key_color, 'left')
    # 右侧：非关键节点
    plot_side_violin(pos, non_key_plot_data[k], non_key_color, 'right')

# 设置x轴标签和标题
plt.xlabel('Clique Size', fontsize=18)  # 设置x轴标签字体大小
plt.ylabel('Participation Count', fontsize=18)  # 设置y轴标签字体大小
plt.xticks(range(1, len(motif_sizes) + 1), [f'{k}' for k in motif_sizes], fontsize=16)  # 设置x轴刻度字体大小
# 设置y轴范围
plt.ylim(0, 200)
plt.yticks(np.arange(0, 201, 50), fontsize=16)  # 设置y轴刻度字体大小

# 添加自定义图例
legend_elements = [
    Line2D([0], [0], color=key_color, lw=8, label='Key Proteins'),
    Line2D([0], [0], color=non_key_color, lw=8, label='Non-key Proteins'),
    Line2D([0], [0], color='black', linestyle='--', lw=1.4000, label='Mean')  # 设置线条宽度
]
plt.legend(handles=legend_elements, loc='upper right', fontsize=15)  # 设置图例字体大小

# 去掉图片上方和右边的框
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 调整布局
plt.tight_layout()

# 设置分辨率为300 dpi并保存图片
plt.savefig('酿酒酵母PPI网络关键节点和非关键节点模体参与计数分布.png', dpi=300)
plt.savefig('酿酒酵母PPI网络关键节点和非关键节点模体参与计数分布.pdf', dpi=300)

# 显示图形
plt.show()

