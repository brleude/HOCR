import pandas as pd
import networkx as nx
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import save_scores
from centrality_improvement import improved_centrality
from bayesian_optimization import optimize_method
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


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
# print(key_proteins,non_key_proteins)

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

# 计算并保存5中经典中心性得分
pr_scores = nx.pagerank(G)
dc_scores = nx.degree_centrality(G)
bc_scores = nx.betweenness_centrality(G)
ec_scores = nx.eigenvector_centrality(G)
cc_scores = nx.closeness_centrality(G)
save_scores(pr_scores, 'pr_scores_sc.txt')
save_scores(dc_scores, 'dc_scores_sc.txt')
save_scores(bc_scores, 'bc_scores_sc.txt')
save_scores(ec_scores, 'ec_scores_sc.txt')
save_scores(cc_scores, 'cc_scores_sc.txt')

# 优化CDR,CSR的参数
best_params_cdr, best_score_cdr = optimize_method(G, 'dc', key_proteins, non_key_proteins, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr)
print("Best score:", best_score_cdr)
best_params_csr, best_score_csr = optimize_method(G, 'cc', key_proteins, non_key_proteins, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr)
print("Best score:", best_score_csr)

# 计算CDR,CSR得分(此处采用论文中的最优参数组合)
params_cdr = {'theta': 0.06, 'lambda_3': 0.38, 'lambda_4': 0.78, 'lambda_5': 0.55, 'lambda_6': 0.39, 
              'lambda_7': 0.89, 'lambda_8': 0.75, 'lambda_9': 0.83, 'lambda_10': 0.48, 'lambda_11': 0.83, 
              'lambda_12': 0.51}
params_csr = {'theta': 1, 'lambda_3': 0.65, 'lambda_4': 0.14, 'lambda_5': 0.95, 'lambda_6': 0.15, 
              'lambda_7': 0.73, 'lambda_8': 0.77, 'lambda_9': 0.27, 'lambda_10': 0.39, 'lambda_11': 0.71, 
              'lambda_12': 0.95}
cdr_scores, _ = improved_centrality(G, 'dc', max_clique=None, params=params_cdr)
_, csr_scores = improved_centrality(G, 'cc', max_clique=None, params=params_csr)

# 保存CDR,CSR得分
save_scores(cdr_scores, 'cdr_scores_sc.txt')
save_scores(csr_scores, 'csr_scores_sc.txt')