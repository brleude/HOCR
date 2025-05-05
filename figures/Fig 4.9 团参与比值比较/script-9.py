import numpy as np
import pandas as pd
import json
import csv
import matplotlib.pyplot as plt
from pathlib import Path
import networkx as nx
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
# 统计平均参与次数
def get_average_counts(count_dict):
    avg = {k: np.mean([v[k] for v in count_dict.values()]) for k in count_dict[next(iter(count_dict))]}
    return avg


###################################### 准备数据 ######################################
######  酿酒酵母相关数据  ######
# 读取关键，非关键蛋白质
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
# 计算团参与情况
key_counts, non_key_counts = compute_motif_participation(G, key_proteins, non_key_proteins, 
                                                         motifs=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
key_avg = get_average_counts(key_counts)
non_key_avg = get_average_counts(non_key_counts)

######  大肠杆菌相关数据  ######
# 提取原始PPI网络
parent3 = current_file.parent.parent.parent
ec_ppi_path = parent3/'data/ec/Ec PPI Network.csv'
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
ec_key_genes = extract_gene_names_from_file(parent3/"data/ec/Ec Key Proteins.tsv")
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
# 计算团参与情况
key_counts_ec, non_key_counts_ec = compute_motif_participation(G_ec, key_proteins_ec, non_key_proteins_ec, 
                                                         motifs=[3, 4, 5, 6, 7, 8, 9, 10])
key_avg_ec = get_average_counts(key_counts_ec)
non_key_avg_ec = get_average_counts(non_key_counts_ec)
# 提取clique大小和对应的平均参与次数
clique_sizes_net1 = np.array(list(key_avg_ec.keys()))
critical_net1 = np.array([key_avg_ec[k] for k in clique_sizes_net1])
noncritical_net1 = np.array([non_key_avg_ec[k] for k in clique_sizes_net1])
ratios_net1 = critical_net1 / noncritical_net1

clique_sizes_net2 = np.array(list(key_avg.keys()))
critical_net2 = np.array([key_avg[k] for k in clique_sizes_net2])
noncritical_net2 = np.array([non_key_avg[k] for k in clique_sizes_net2])
noncritical_net2[noncritical_net2 == 0] = 0.05  # 处理除零情况
ratios_net2 = critical_net2 / noncritical_net2


###################################### 绘制图形 ######################################
plt.figure(figsize=(11.69, 6), dpi=300)  #图片尺寸和分辨率
ax = plt.gca()

# 使用所有可能的clique大小（3 - 12）
all_clique_sizes = np.arange(3, 13)
x = np.arange(len(all_clique_sizes))
width = 0.35

# 创建完整数据数组（用NaN填充缺失值）
full_ratios_net1 = np.full(len(all_clique_sizes), np.nan)
full_ratios_net2 = np.full(len(all_clique_sizes), np.nan)

# 填充网络1数据
for i, size in enumerate(clique_sizes_net1):
    idx = np.where(all_clique_sizes == size)[0][0]
    full_ratios_net1[idx] = ratios_net1[i]

# 填充网络2数据
for i, size in enumerate(clique_sizes_net2):
    idx = np.where(all_clique_sizes == size)[0][0]
    full_ratios_net2[idx] = ratios_net2[i]

# 绘图顺序，E.coli 在右侧
rects1 = ax.bar(x + width / 2, full_ratios_net1, width,
                label='E.coli', color='#4D96FF', edgecolor='black', alpha=0.8, linewidth=1.4, zorder=2)
rects2 = ax.bar(x - width / 2, full_ratios_net2, width,
                label='S.cerevisiae', color='#FF6B6B', edgecolor='black', alpha=0.8, linewidth=1.4, zorder=2)

# 设置坐标轴
ax.set_ylabel('Average Clique Participation Ratio\n(Key Proteins / Non-key Proteins)', fontsize=20) # y轴标签字体大小
ax.set_xlabel('Clique Size', fontsize=20) # x轴标签字体大小
ax.set_xticks(x)
ax.set_xticklabels(all_clique_sizes)
ax.tick_params(axis='both', which='major', labelsize=18) # 坐标轴刻度字体大小
# ax.set_title('Comparison of Clique Participation', fontsize=12, pad=15)
# 添加网格
ax.grid(True, axis='y', ls="--", alpha=0.7, zorder=1)
# 添加图例
ax.legend(frameon=True, loc='upper left', fontsize=18)  # 图例字体大小
# 添加数据标签
def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        if not np.isnan(height):  # 只标注有效数据
            ax.annotate(f'{height:.1f}',
                        xy=(rect.get_x() + rect.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=15)
autolabel(rects1)
autolabel(rects2)
# 调整y轴范围，避免顶部空白
ax.set_ylim(0, max(np.nanmax(full_ratios_net1), np.nanmax(full_ratios_net2)) * 1.1)
# 去掉上方和右侧的边框
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# 调整布局
plt.tight_layout()
# 保存图像
plt.savefig('团参与比值比较.png', dpi=300)
plt.savefig('团参与比值比较.pdf', dpi=300)
plt.show()