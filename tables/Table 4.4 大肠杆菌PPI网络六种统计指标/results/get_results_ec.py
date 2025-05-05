import pandas as pd
import networkx as nx
import json
import csv
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

# 计算传统中心性
dc_scores_ec = nx.degree_centrality(G_ec)
bc_scores_ec = nx.betweenness_centrality(G_ec)
cc_scores_ec = nx.closeness_centrality(G_ec)
ec_scores_ec = nx.eigenvector_centrality(G_ec)
pr_scores_ec = nx.pagerank(G_ec)

save_scores(dc_scores_ec, 'dc_scores_ec.txt')
save_scores(bc_scores_ec, 'bc_scores_ec.txt')
save_scores(cc_scores_ec, 'cc_scores_ec.txt')
save_scores(ec_scores_ec, 'ec_scores_ec.txt')
save_scores(pr_scores_ec, 'pr_scores_ec.txt')

# 优化CDR_5,CSR_3的参数
best_params_cdr, best_score_cdr = optimize_method(G_ec, 'cc', key_proteins_ec, non_key_proteins_ec, max_clique=5,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr)
print("Best score:", best_score_cdr)
best_params_csr, best_score_csr = optimize_method(G_ec, 'cc', key_proteins_ec, non_key_proteins_ec, max_clique=3,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr)
print("Best score:", best_score_csr)

# 计算CDR_5,CSR_3得分(此处采用论文中的最优参数组合)
params_cdr = {'theta': 0.15, 'lambda_3': 0.03, 'lambda_4': 0.06, 'lambda_5': 0.06}
params_csr = {'theta': 0.03, 'lambda_3': 0.03}
cdr_scores_ec, _ = improved_centrality(G_ec, 'cc', max_clique=5, params=params_cdr)
_, csr_scores_ec = improved_centrality(G_ec, 'cc', max_clique=3, params=params_csr)

# 保存CDR_5,CSR_3得分
save_scores(cdr_scores_ec, 'cdr_scores_ec.txt')
save_scores(csr_scores_ec, 'csr_scores_ec.txt')