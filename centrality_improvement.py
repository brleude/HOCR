import networkx as nx
from collections import defaultdict
from itertools import combinations
import random
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

def find_motifs(G, max_clique=None):
    """查找网络中的所有团

    参数:
        G: 网络图对象
        max_clique: 最大团大小，如果为None则自动计算网络中的最大团

    返回:
        包含所有motif的列表，每个motif是一个排序后的元组
    """
    all_motifs = []
    if max_clique is None:
        max_clique = max(len(c) for c in nx.find_cliques(G))
    for size in range(3, max_clique + 1):
        cliques = []
        for clique in nx.find_cliques(G):
            if len(clique) >= size:
                for sub_clique in combinations(clique, size):
                    if all(G.has_edge(u, v) for u, v in combinations(sub_clique, 2)):
                        cliques.append(tuple(sorted(sub_clique)))
        all_motifs.extend(list(set(cliques)))
    return all_motifs


def compute_combined_motif_network(edges, motifs):
    """构建考虑权重的高阶加权网络

    参数:
        edges: 原始网络的边列表
        motifs: 所有motif列表

    返回:
        高阶加权后的网络G_prime
    """
    G = nx.Graph()
    G.add_edges_from(edges)
    edge_weights = defaultdict(int)
    for u, v in G.edges():
        edge_weights[(u, v)] = 1

    clique_sizes = set(len(motif) for motif in motifs)
    for size in clique_sizes:
        size_motifs = [motif for motif in motifs if len(motif) == size]
        for clique in size_motifs:
            for u, v in combinations(clique, 2):
                edge_weights[(u, v)] += size

    G_prime = nx.Graph()
    G_prime.add_nodes_from(G.nodes())
    for (u, v), weight in edge_weights.items():
        G_prime.add_edge(u, v, weight=weight)

    return G_prime


def count_motif(node, motifs):
    """计算节点参与的motif数量

    参数:
        node: 目标节点
        motifs: 所有motif列表

    返回:
        该节点参与的motif总数
    """
    return sum(1 for motif in motifs if node in motif)


def improved_centrality(G, base_scores='dc', max_clique=None, params=None):
    """计算改进的中心性指标(CDR和CSR)

    参数:
        G: 网络图对象
        base_scores: 基础中心性类型，可选'dc'(度中心性，默认)、'bc'(中介中心性)、
                    'cc'(接近中心性)、'ec'(特征向量中心性)、'pr'(PageRank)
                    或直接提供预计算的中心性分数字典
        max_clique: 最大团大小，None则自动计算网络中的最大团
        params: 包含theta和lambda参数的字典，None则随机生成

    返回:
        两个字典: CDR分数字典, CSR分数字典
    """
    if base_scores == 'bc':
        base_scores = nx.betweenness_centrality(G)
    elif base_scores == 'cc':
        base_scores = nx.closeness_centrality(G)
    elif base_scores == 'dc':
        base_scores = nx.degree_centrality(G)
    elif base_scores == 'ec':
        base_scores = nx.eigenvector_centrality(G)
    elif base_scores == 'pr':
        base_scores = nx.pagerank(G)
    else:
        base_scores = base_scores

    if params is None:
        # 设置随机种子
        random.seed(42)
        # 自动生成参数
        params = {}
        params['theta'] = random.uniform(0, 1)
        if max_clique is None:
            max_clique = max(len(c) for c in nx.find_cliques(G))
        for size in range(3, max_clique + 1):
            params[f'lambda_{size}'] = random.uniform(0, 1)
            
    motifs = find_motifs(G, max_clique)
    G_prime = compute_combined_motif_network(G.edges, motifs)

    motif_degree = dict(G_prime.degree())
    motif_strength = dict(nx.degree(G_prime, weight='weight'))

    node_motif_counts = {}
    clique_sizes = set(len(motif) for motif in motifs)
    for size in clique_sizes:
        size_motifs = [motif for motif in motifs if len(motif) == size]
        for node in G.nodes():
            count = count_motif(node, size_motifs)
            if node not in node_motif_counts:
                node_motif_counts[node] = {}
            node_motif_counts[node][size] = count

    total_degree = sum(motif_degree.values()) + 1e-10
    norm_degree = {n: d / total_degree for n, d in motif_degree.items()}

    total_strength = sum(motif_strength.values()) + 1e-10
    norm_strength = {n: s / total_strength for n, s in motif_strength.items()}

    def _calculate_scores(base_scores, norm_motif_vals):
        adjusted = {}
        for n in G.nodes():
            base_effect = (norm_motif_vals[n]) ** params.get('theta', 1.0)
            high_order_correction = 1
            for size in node_motif_counts[n]:
                lambda_n = params.get(f'lambda_{size}', 0)
                total_size = sum(node_motif_counts[m][size] for m in G.nodes()) + 1e-10
                high_order_correction += lambda_n * (node_motif_counts[n][size] / total_size)
            adjusted[n] = base_scores[n] * base_effect * high_order_correction

        total = sum(adjusted.values())
        return {n: score / total for n, score in adjusted.items()}

    cdr = _calculate_scores(base_scores, norm_degree)
    csr = _calculate_scores(base_scores, norm_strength)

    return cdr, csr

# # 使用示例
# if __name__ == "__main__":
#     # 创建一个示例图
#     G = nx.Graph()
#     # 添加边
#     G.add_edges_from([(1, 2), (2, 3), (3, 1), (1, 4), (4, 5), (3,4),(4,2)])
#     base_scores = 'dc'
#     max_clique = None
#     params = {'theta': 0.9, 'lambda_3': 0.38, 'lambda_4': 0.74}
#     cdr, csr = improved_centrality(G, base_scores, max_clique, params)
#     print("CDR:", cdr)
#     print("CSR:", csr)