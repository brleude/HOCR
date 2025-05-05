import pandas as pd
import networkx as nx
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
# 计算节点参与各大小的团的情况
def calculate_clique_stats(counts):
    clique_sizes = set()
    for protein in counts:
        clique_sizes.update(counts[protein].keys())
    clique_sizes = sorted(clique_sizes)

    total_nodes = len(counts)
    node_counts = {size: 0 for size in clique_sizes}
    node_proportions = {size: 0 for size in clique_sizes}

    for protein in counts:
        for size in clique_sizes:
            if counts[protein].get(size, 0) > 0:
                node_counts[size] += 1

    for size in clique_sizes:
        node_proportions[size] = node_counts[size] / total_nodes

    return node_counts, node_proportions
# 将结果保存为表格
def save_to_table(key_proportions, non_key_proportions, all_proportions, file_path):
    """
    将关键节点、非关键节点和所有节点的团比例信息保存为表格
    :param key_proportions: 关键节点的团比例字典
    :param non_key_proportions: 非关键节点的团比例字典
    :param all_proportions: 所有节点的团比例字典
    :param file_path: 保存文件的路径
    """
    # 创建 DataFrame
    data = {
        'Node Type': ['Key Node Proportions', 'Non-key Node Proportions', 'All Nodes Proportions'],
        **key_proportions,
        **non_key_proportions,
        **all_proportions
    }
    df = pd.DataFrame({
        'Node Type': ['Key Node Proportions', 'Non-key Node Proportions', 'All Nodes Proportions']
    })
    clique_sizes = sorted(set(key_proportions.keys()) | set(non_key_proportions.keys()) | set(all_proportions.keys()))
    for size in clique_sizes:
        df[size] = [key_proportions.get(size), non_key_proportions.get(size), all_proportions.get(size)]

    # 保存为 CSV 文件
    df.to_csv(file_path, index=False)
    print('数据已保存')

    return df


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
# 加载网络边信息
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
# 构建PPI网络
ppi_edges = load_ppi_edges(ppi_file)
G = nx.Graph()
G.add_edges_from(ppi_edges)


# 计算相关信息
key_counts, non_key_counts = compute_motif_participation(G, key_proteins, non_key_proteins, 
                                                         motifs=[3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
# 合并 key_counts 和 non_key_counts
all_counts = {**key_counts, **non_key_counts}
# 计算 key_counts 的统计信息
key_node_counts, key_node_proportions = calculate_clique_stats(key_counts)
# 计算 non_key_counts 的统计信息
non_key_node_counts, non_key_node_proportions = calculate_clique_stats(non_key_counts)
# 计算所有节点的统计信息
all_node_counts, all_node_proportions = calculate_clique_stats(all_counts)

# 保存为表格
file_path = '酿酒酵母PPI网络中参与各大小的团的节点所占百分比.csv'
df = save_to_table(key_node_proportions, non_key_node_proportions, all_node_proportions, file_path)
print(df)