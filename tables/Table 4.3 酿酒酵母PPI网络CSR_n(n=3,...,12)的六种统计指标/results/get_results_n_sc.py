import pandas as pd
import networkx as nx
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import save_scores
from centrality_improvement import improved_centrality
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


params_csr_3 = {'theta': 1, 'lambda_3': 0.43}
_, csr_scores_3 = improved_centrality(G, 'cc', max_clique=3, params=params_csr_3)

params_csr_4 = {'theta': 1, 'lambda_3': 0.68, 'lambda_4': 0.85}
_, csr_scores_4 = improved_centrality(G, 'cc', max_clique=4, params=params_csr_4)

params_csr_5 = {'theta': 1, 'lambda_3': 0.81, 'lambda_4': 0.59, 'lambda_5':0.76}
_, csr_scores_5 = improved_centrality(G, 'cc', max_clique=5, params=params_csr_5)

params_csr_6 = {'theta': 1, 'lambda_3': 0.37, 'lambda_4': 0.18, 'lambda_5':0.35, 'lambda_6': 0.76}
_, csr_scores_6 = improved_centrality(G, 'cc', max_clique=6, params=params_csr_6)

params_csr_7 = {'theta': 1, 'lambda_3': 0.37, 'lambda_4': 0.76, 'lambda_5':0.49, 'lambda_6': 0.12, 
                'lambda_7': 0.00}
_, csr_scores_7 = improved_centrality(G, 'cc', max_clique=7, params=params_csr_7)

params_csr_8 = {'theta': 1, 'lambda_3': 0.71, 'lambda_4': 0.24, 'lambda_5':0.47, 'lambda_6': 0.74, 
                'lambda_7': 0.48, 'lambda_8': 0.93}
_, csr_scores_8 = improved_centrality(G, 'cc', max_clique=8, params=params_csr_8)

params_csr_9 = {'theta': 1, 'lambda_3': 0.45, 'lambda_4': 0.57, 'lambda_5':0.86, 'lambda_6': 0.65, 
                'lambda_7': 0.56, 'lambda_8': 0.33, 'lambda_9': 0.16}
_, csr_scores_9 = improved_centrality(G, 'cc', max_clique=9, params=params_csr_9)

params_csr_10 = {'theta': 1, 'lambda_3': 0.85, 'lambda_4': 0.01, 'lambda_5':0.53, 'lambda_6': 0.26, 
                'lambda_7': 0.53, 'lambda_8': 0.35, 'lambda_9': 0.28, 'lambda_10': 0.55}
_, csr_scores_10 = improved_centrality(G, 'cc', max_clique=10, params=params_csr_10)

params_csr_11 = {'theta': 1, 'lambda_3': 0.49, 'lambda_4': 0.32, 'lambda_5':0.66, 'lambda_6': 0.57, 
                'lambda_7': 0.47, 'lambda_8': 0.47, 'lambda_9': 0.60, 'lambda_10': 0.55, 'lambda_11': 0.33}
_, csr_scores_11 = improved_centrality(G, 'cc', max_clique=11, params=params_csr_11)

params_csr_12 = {'theta': 1, 'lambda_3': 0.65, 'lambda_4': 0.14, 'lambda_5':0.95, 'lambda_6': 0.15, 
                'lambda_7': 0.73, 'lambda_8': 0.77, 'lambda_9': 0.27, 'lambda_10': 0.39, 'lambda_11': 0.71, 'lambda_12': 0.95}
_, csr_scores_12 = improved_centrality(G, 'cc', max_clique=12, params=params_csr_12)


save_scores(csr_scores_3, 'csr_scores_3_sc.txt')
save_scores(csr_scores_4, 'csr_scores_4_sc.txt')
save_scores(csr_scores_5, 'csr_scores_5_sc.txt')
save_scores(csr_scores_6, 'csr_scores_6_sc.txt')
save_scores(csr_scores_7, 'csr_scores_7_sc.txt')
save_scores(csr_scores_8, 'csr_scores_8_sc.txt')
save_scores(csr_scores_9, 'csr_scores_9_sc.txt')
save_scores(csr_scores_10, 'csr_scores_10_sc.txt')
save_scores(csr_scores_11, 'csr_scores_11_sc.txt')
save_scores(csr_scores_12, 'csr_scores_12_sc.txt')