import optuna
import networkx as nx
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score
from centrality_improvement import improved_centrality
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

def evaluate_rank(rank_scores, key_nodes, non_key_nodes, all_nodes,
                  metric='ap', threshold=None):
    """综合评估排名得分

    参数:
        rank_scores: 排名得分字典 {nodes: score}
        key_nodes: 关键节点列表
        non_key_nodes: 非关键节点列表
        all_nodes: 网络中所有节点列表
        metric: 评估指标，可选 'auc', 'ap'(默认), 'f1'
        threshold: 仅当metric='f1'时有效，分类阈值，None表示自动选择最佳阈值

    返回:
        对应指标的评估分数
    """
    # 获取共有的节点
    common_key = set(key_nodes) & set(all_nodes)
    common_non_key = set(non_key_nodes) & set(all_nodes)

    # 准备标签和得分
    y_true = []
    y_score = []

    for protein in common_key:
        y_true.append(1)
        y_score.append(rank_scores.get(protein, 0))

    for protein in common_non_key:
        y_true.append(0)
        y_score.append(rank_scores.get(protein, 0))

    # 检查是否有两类样本
    if len(set(y_true)) < 2:
        return 0.5 if metric == 'auc' else 0.0

    # 根据metric参数选择评估指标
    if metric == 'auc':
        return roc_auc_score(y_true, y_score)
    elif metric == 'ap':
        return average_precision_score(y_true, y_score)
    elif metric == 'f1':
        if threshold is None:
            # 自动选择最佳阈值
            thresholds = sorted(set(y_score), reverse=True)
            best_f1 = 0
            for t in thresholds:
                y_pred = [1 if score >= t else 0 for score in y_score]
                current_f1 = f1_score(y_true, y_pred)
                if current_f1 > best_f1:
                    best_f1 = current_f1
            return best_f1
        else:
            y_pred = [1 if score >= threshold else 0 for score in y_score]
            return f1_score(y_true, y_pred)
    else:
        raise ValueError(f"Invalid metric '{metric}'. Choose from 'auc', 'ap', 'f1'")

def optimize_method(G, base_scores, key_nodes, non_key_nodes=None, max_clique=None,
                   n_trials=50, metric='ap', threshold=None, rank_type='csr',
                   verbose=True):
    """使用贝叶斯优化方法

    参数:
        G: 网络对象
        base_scores: 基础中心性得分字典 {node: score}
        key_nodes: 关键节点列表
        non_key_nodes: 非关键节点列表，若为None则认为其是G中的节点去掉key_nodes的节点
        max_clique: 最大团大小，若为None则考虑网络中最大团的大小
        n_trials: 迭代次数
        metric: 评估指标，可选 'auc', 'ap'(默认), 'f1'
        threshold: 仅当metric='f1'时有效，分类阈值，None表示自动选择最佳阈值
        rank_type: 排名类型，可选 'csr'(默认) 或 'cdr'
        verbose: 是否输出训练过程，默认为True

    返回:
        最佳参数和最佳得分
    """
    all_nodes = list(G.nodes())
    if non_key_nodes is None:
        non_key_nodes = list(set(all_nodes) - set(key_nodes))

    if max_clique is None:
        max_clique = max(len(c) for c in nx.find_cliques(G))

    def objective(trial):
        # 动态生成参数
        params = {}
        params['theta'] = trial.suggest_float('theta', 0, 1.0)
        for size in range(3, max_clique + 1):
            params[f'lambda_{size}'] = trial.suggest_float(f'lambda_{size}', 0, 1.0)

        cdr, csr = improved_centrality(G, base_scores, max_clique, params)
        rank_scores = csr if rank_type == 'csr' else cdr
        return evaluate_rank(rank_scores, key_nodes, non_key_nodes, all_nodes, metric, threshold)

    # 保存当前日志级别
    original_log_level = optuna.logging.get_verbosity()
    
    if not verbose:
        # 完全禁用所有输出
        optuna.logging.set_verbosity(optuna.logging.ERROR)
        optuna.logging.disable_default_handler()
        optuna.logging.disable_propagation()

    study = optuna.create_study(direction='maximize')
    study.optimize(objective, n_trials=n_trials, show_progress_bar=verbose)

    if not verbose:
        # 恢复日志设置
        optuna.logging.set_verbosity(original_log_level)
        optuna.logging.enable_default_handler()
        optuna.logging.enable_propagation()

    return study.best_params, study.best_value