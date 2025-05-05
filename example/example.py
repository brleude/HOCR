import sys
import networkx as nx
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from centrality_improvement import improved_centrality
from bayesian_optimization import optimize_method
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


# 构建原始网络
ppi_test = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4),(3,5),(4,5),(4,6)]
G_test = nx.Graph()
G_test.add_edges_from(ppi_test)

# 计算CDR,CSR
# 如果有关键节点数据
base_scores = 'cc'
max_clique = None
key_nodes = [3, 4]
non_key_nodes = [1, 2, 5, 6]
best_params, _ = optimize_method(G_test, base_scores=base_scores, key_nodes=key_nodes, non_key_nodes=non_key_nodes, max_clique=max_clique, verbose=False)
cdr1, csr1 = improved_centrality(G_test, base_scores=base_scores, max_clique=max_clique, params=best_params)
print('如果有关键节点数据：')
print("CDR:", cdr1)
print("CSR:", csr1)

# 如果没有关键节点数据
params = None
cdr2, csr2 = improved_centrality(G_test, base_scores=base_scores, max_clique=max_clique, params=params)
print('如果没有关键节点数据：')
print('CDR:', cdr2)
print('CSR:', csr2)