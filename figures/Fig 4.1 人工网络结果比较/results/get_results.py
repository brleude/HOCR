import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from centrality_improvement import improved_centrality
from bayesian_optimization import evaluate_rank, optimize_method
import pickle
import numpy as np
import networkx as nx
from pathlib import Path
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 获取当前脚本的绝对路径
current_file = Path(__file__).resolve()
# 上三级目录
parent4 = current_file.parent.parent.parent.parent
data_path = parent4/'data/artificial_network.pkl'
# 加载图对象 G_test 和内在质量 Q
with open(data_path, 'rb') as f:
    data = pickle.load(f)
    G_test = data['G']
    Q = data['Q']

ppi_edges_test = list(G_test.edges())

# 获取内在质量前top_n的关键节点
top_n = 100 # 选择100个关键节点
key_nodes_100 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
non_key_nodes_100 = np.setdiff1d(range(G_test.number_of_nodes()), key_nodes_100)

top_n = 200 # 选择200个关键节点
key_nodes_200 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
non_key_nodes_200 = np.setdiff1d(range(G_test.number_of_nodes()), key_nodes_200)

top_n = 300 # 选择300个关键节点
key_nodes_300 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
non_key_nodes_300 = np.setdiff1d(range(G_test.number_of_nodes()), key_nodes_300)

top_n = 400 # 选择400个关键节点
key_nodes_400 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
non_key_nodes_400 = np.setdiff1d(range(G_test.number_of_nodes()), key_nodes_400)

top_n = 500 # 选择500个关键节点
key_nodes_500 = np.argsort(Q)[-top_n:][::-1]  # 按Q值降序排列的前top_n个节点
non_key_nodes_500 = np.setdiff1d(range(G_test.number_of_nodes()), key_nodes_500)

# 计算人工网络传统中心性指标
pr = nx.pagerank(G_test)
bc = nx.betweenness_centrality(G_test)
cc = nx.closeness_centrality(G_test)
dc = nx.degree_centrality(G_test)
ec = nx.eigenvector_centrality(G_test)

# 优化5种参数情况下CDR,CSR的参数组合
best_params_cdr_100, best_score_cdr_100 = optimize_method(G_test, 'ec', key_nodes_100, non_key_nodes_100, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr_100)
print("Best score:", best_score_cdr_100)
best_params_csr_100, best_score_csr_100 = optimize_method(G_test, 'ec', key_nodes_100, non_key_nodes_100, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr_100)
print("Best score:", best_score_csr_100)

best_params_cdr_200, best_score_cdr_200 = optimize_method(G_test, 'ec', key_nodes_200, non_key_nodes_200, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr_200)
print("Best score:", best_score_cdr_200)
best_params_csr_200, best_score_csr_200 = optimize_method(G_test, 'dc', key_nodes_200, non_key_nodes_200, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr_200)
print("Best score:", best_score_csr_200)

best_params_cdr_300, best_score_cdr_300 = optimize_method(G_test, 'cc', key_nodes_300, non_key_nodes_300, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr_300)
print("Best score:", best_score_cdr_300)
best_params_csr_300, best_score_csr_300 = optimize_method(G_test, 'dc', key_nodes_300, non_key_nodes_300, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr_300)
print("Best score:", best_score_csr_300)

best_params_cdr_400, best_score_cdr_400 = optimize_method(G_test, 'cc', key_nodes_400, non_key_nodes_400, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr_400)
print("Best score:", best_score_cdr_400)
best_params_csr_400, best_score_csr_400 = optimize_method(G_test, 'dc', key_nodes_400, non_key_nodes_400, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr_400)
print("Best score:", best_score_csr_400)

best_params_cdr_500, best_score_cdr_500 = optimize_method(G_test, 'cc', key_nodes_500, non_key_nodes_500, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='cdr')
print("Best parameters:", best_params_cdr_500)
print("Best score:", best_score_cdr_500)
best_params_csr_500, best_score_csr_500 = optimize_method(G_test, 'dc', key_nodes_500, non_key_nodes_500, max_clique=None,
                                          n_trials=50, metric='ap', threshold=None, rank_type='csr')
print("Best parameters:", best_params_csr_500)
print("Best score:", best_score_csr_500)


# 计算选取不同个数(100,200,300,400,500)关键节点下人工网络的最优CDR和CSR(最优参数组合经过调参选取,这里使用Fig4.1结果的对应参数)
params_cdr_100 = {'theta': 0.27, 'lambda_3': 0.53, 'lambda_4': 0.97, 'lambda_5': 0.98, 'lambda_6': 0.96, 'lambda_7': 0.68, 'lambda_8': 0.27}
params_csr_100 = {'theta': 0.57, 'lambda_3': 0.45, 'lambda_4': 0.42, 'lambda_5': 0.83, 'lambda_6': 0.09, 'lambda_7': 0.03, 'lambda_8': 0.14}
cdr_100, _ = improved_centrality(G_test, 'ec', max_clique=None, params=params_cdr_100)
_, csr_100 = improved_centrality(G_test, 'ec', max_clique=None, params=params_csr_100)

params_cdr_200 = {'theta': 0.59, 'lambda_3': 0.71, 'lambda_4': 0.93, 'lambda_5': 0.26, 'lambda_6': 0.75, 'lambda_7': 0.08, 'lambda_8': 0.04}
params_csr_200 = {'theta': 0.44, 'lambda_3': 0.38, 'lambda_4': 0.10, 'lambda_5': 0.10, 'lambda_6': 0.02, 'lambda_7': 0.04, 'lambda_8': 0.15}
cdr_200, _ = improved_centrality(G_test, 'ec', max_clique=None, params=params_cdr_200)
_, csr_200 = improved_centrality(G_test, 'dc', max_clique=None, params=params_csr_200)

params_cdr_300 = {'theta': 0.50, 'lambda_3': 0.80, 'lambda_4': 0.27, 'lambda_5': 0.05, 'lambda_6': 0.42, 'lambda_7': 0.03, 'lambda_8': 0.64}
params_csr_300 = {'theta': 0.11, 'lambda_3': 0.53, 'lambda_4': 0.45, 'lambda_5': 0.20, 'lambda_6': 0.50, 'lambda_7': 0.32, 'lambda_8': 0.90}
cdr_300, _ = improved_centrality(G_test, 'cc', max_clique=None, params=params_cdr_300)
_, csr_300 = improved_centrality(G_test, 'dc', max_clique=None, params=params_csr_300)

params_cdr_400 = {'theta': 0.48, 'lambda_3': 0.33, 'lambda_4': 0.74, 'lambda_5': 0.35, 'lambda_6': 0.40, 'lambda_7': 0.55, 'lambda_8': 0.75}
params_csr_400 = {'theta': 0.22, 'lambda_3': 0.56, 'lambda_4': 0.92, 'lambda_5': 0.86, 'lambda_6': 0.99, 'lambda_7': 0.50, 'lambda_8': 0.76}
cdr_400, _ = improved_centrality(G_test, 'cc', max_clique=None, params=params_cdr_400)
_, csr_400 = improved_centrality(G_test, 'dc', max_clique=None, params=params_csr_400)

params_cdr_500 = {'theta': 0.53, 'lambda_3': 0.17, 'lambda_4': 0.94, 'lambda_5': 0.28, 'lambda_6': 0.51, 'lambda_7': 0.12, 'lambda_8': 0.93}
params_csr_500 = {'theta': 0.23, 'lambda_3': 0.36, 'lambda_4': 0.80, 'lambda_5': 0.15, 'lambda_6': 0.41, 'lambda_7': 0.49, 'lambda_8': 0.78}
cdr_500, _ = improved_centrality(G_test, 'cc', max_clique=None, params=params_cdr_500)
_, csr_500 = improved_centrality(G_test, 'dc', max_clique=None, params=params_csr_500)

# 评估多种中心性指标下人工网络不同版本的AUC、平均排名和平均精度
def evaluate_clique_performance(true_cores, dc_scores, cc_scores, bc_scores, ec_scores, pr_scores, cdr_scores, csr_scores):
    """评估多种中心性指标下的AUC、平均排名和平均精度"""
    all_proteins = list(dc_scores.keys())
    key_nodes = true_cores
    non_key_nodes = [node for node in all_proteins if node not in key_nodes]

    score_dicts = {
        'DC': dc_scores,
        'CC': cc_scores,
        'BC': bc_scores,
        'EC': ec_scores,
        'PR': pr_scores,
        'CDR': cdr_scores,
        'CSR': csr_scores
    }

    results = {}
    for score_name, score_dict in score_dicts.items():
        # 计算AUC
        auc_value = evaluate_rank(score_dict, key_nodes, non_key_nodes, all_proteins, metric='auc')

        # 计算平均排名
        def get_avg_rank(scores_dict):
            sorted_nodes = sorted(scores_dict, key=scores_dict.get, reverse=True)
            return np.mean([sorted_nodes.index(n) + 1 for n in key_nodes if n in scores_dict])

        avg_rank = get_avg_rank(score_dict)

        # 计算平均精度
        ap_value = evaluate_rank(score_dict, key_nodes, non_key_nodes, all_proteins, metric='ap')

        results[f'{score_name} AUC'] = auc_value
        results[f'{score_name} AvgRank'] = avg_rank
        results[f'{score_name} AP'] = ap_value

    return results
       
results_100 = evaluate_clique_performance(key_nodes_100, dc, cc, bc, ec, pr, cdr_100, csr_100)
print(results_100)

results_200 = evaluate_clique_performance(key_nodes_200, dc, cc, bc, ec, pr, cdr_200, csr_200)
print(results_200)

results_300 = evaluate_clique_performance(key_nodes_300, dc, cc, bc, ec, pr, cdr_300, csr_300)
print(results_300)

results_400 = evaluate_clique_performance(key_nodes_400, dc, cc, bc, ec, pr, cdr_400, csr_400)
print(results_400)

results_500 = evaluate_clique_performance(key_nodes_500, dc, cc, bc, ec, pr, cdr_500, csr_500)
print(results_500)


with open('results_100.pkl', 'wb') as f:
    pickle.dump({'results': results_100}, f)
with open('results_200.pkl', 'wb') as f:
    pickle.dump({'results': results_200}, f)
with open('results_300.pkl', 'wb') as f:
    pickle.dump({'results': results_300}, f)
with open('results_400.pkl', 'wb') as f:
    pickle.dump({'results': results_400}, f)
with open('results_500.pkl', 'wb') as f:
    pickle.dump({'results': results_500}, f)