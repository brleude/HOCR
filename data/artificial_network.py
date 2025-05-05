import numpy as np
import networkx as nx
from tqdm import tqdm
import pickle
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 生成无标度网络
def generate_network(N=1000, avg_degree=10, initial_nodes=5):
    """生成无标度无向网络"""
    # 生成对数正态分布
    q = np.random.lognormal(mean=-1, sigma=0.8, size=N)
    q = (q - q.min()) / (q.max() - q.min()) * 0.9 + 0.1  # 映射到[0.1,1.0]
    Q=q
    
    # 创建初始完全连接的种子网络
    G = nx.complete_graph(initial_nodes)
    degrees = np.zeros(N, dtype=int)  # 预分配所有节点的度数数组
    degrees[:initial_nodes] = initial_nodes - 1  # 初始节点的度数
    
    # 预分配边列表（减少动态添加的开销）
    edge_list = list(G.edges())
    
    # 添加剩余节点和边（向量化优化）
    for i in tqdm(range(initial_nodes, N), desc="Add nodes"):
        # 计算现有节点的连接概率
        existing_nodes = np.arange(i)
        weights = Q[existing_nodes] * (degrees[existing_nodes] + 1)
        prob = weights / weights.sum()
        
        # 一次性选择m个连接目标（无放回抽样）
        m = 2  # 每个新节点的初始连接数
        targets = np.random.choice(existing_nodes, size=m, p=prob, replace=False)
        
        # 更新边列表和度数
        for j in targets:
            edge_list.append((i, j))
            degrees[i] += 1
            degrees[j] += 1
    
    # 批量添加所有边
    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edge_list)
    
    # 继续添加剩余边（优化概率计算）
    current_edges = G.number_of_edges()
    total_edges = int(N * avg_degree / 2)
    
    # 预计算动态归一化因子
    max_Q = Q.max()
    max_degree = degrees.max()
    max_p = (max_Q * (max_degree + 1)) ** 2
    
    with tqdm(total=total_edges - current_edges, desc="Add other edges") as pbar:
        while current_edges < total_edges:
            # 随机选择两个不同节点（避免重复检查）
            i, j = np.random.choice(N, size=2, replace=False)
            if not G.has_edge(i, j):
                # 快速计算连接概率
                p_ij = (Q[i] * (degrees[i] + 1)) * (Q[j] * (degrees[j] + 1))
                if np.random.rand() < p_ij / max_p:
                    G.add_edge(i, j)
                    degrees[i] += 1
                    degrees[j] += 1
                    current_edges += 1
                    pbar.update(1)
                    
                    # 动态更新max_p（每100条边更新一次）
                    if current_edges % 100 == 0:
                        max_degree = degrees.max()
                        max_p = (max_Q * (max_degree + 1)) ** 2
    
    return G, Q

# 生成网络并获取内在质量分数
G_test, Q = generate_network(N=1000, avg_degree=10)

# 保存图对象 G_test 和内在质量 Q
with open('artificial_network.pkl', 'wb') as f:
    pickle.dump({'G': G_test, 'Q': Q}, f)

# 保存图对象 G_test 为边列表文件
nx.write_edgelist(G_test, 'graph.edgelist')

# 获取内在质量前top_n的关键节点并保存
top_n = 100 # 选择100个关键节点
key_nodes_100 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
with open('key_nodes_100.txt', 'w') as file:
    file.write(','.join(str(item) for item in key_nodes_100))

top_n = 200 # 选择200个关键节点
key_nodes_200 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
with open('key_nodes_200.txt', 'w') as file:
    file.write(','.join(str(item) for item in key_nodes_200))
    
top_n = 300 # 选择300个关键节点
key_nodes_300 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
with open('key_nodes_300.txt', 'w') as file:
    file.write(','.join(str(item) for item in key_nodes_300))

top_n = 400 # 选择400个关键节点
key_nodes_400 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
with open('key_nodes_400.txt', 'w') as file:
    file.write(','.join(str(item) for item in key_nodes_400))

top_n = 500 # 选择500个关键节点
key_nodes_500 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
with open('key_nodes_500.txt', 'w') as file:
    file.write(','.join(str(item) for item in key_nodes_500))