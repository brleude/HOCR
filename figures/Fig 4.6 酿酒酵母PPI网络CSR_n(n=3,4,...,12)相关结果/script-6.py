import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import ticker
from sklearn.metrics import roc_curve, auc, average_precision_score, precision_recall_curve
from matplotlib import gridspec
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import load_scores
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 酿酒酵母PPI网络上CSR_n(n=3,4,...,12)相关结果
def evaluate_parameter_combinations_4in1(score_dicts, key_proteins, non_key_proteins, G, 
                                       top_n_range=None, step=1):
    """
    四合一评估函数 (2x2布局)
    
    参数:
        score_dicts: [{'name':str, 'scores':dict}] 方法得分列表
        key_proteins: 关键蛋白列表
        non_key_proteins: 非关键蛋白列表
        G: PPI网络图
        top_n_range: 累积曲线最大范围（蛋白数）
        step: 累积曲线步长
    
    返回:
        各方法评估指标字典
    """
    # ====== 数据预处理 ======
    network_nodes = set(G.nodes())
    valid_key = [p for p in key_proteins if p in network_nodes]
    valid_non_key = [p for p in non_key_proteins if p in network_nodes]
    
    y_true = np.array([1]*len(valid_key) + [0]*len(valid_non_key))
    all_proteins = np.array(valid_key + valid_non_key)

    # ====== 可视化配置 ======
    plt.style.use('default')
    plt.rcParams.update({
        'figure.facecolor': 'white',
        # 'axes.grid': True,
        # 'grid.color': '#e6e6e6',
        # 'grid.linewidth': 0.8,
        'axes.edgecolor': '#333333',
        'axes.linewidth': 1.4,
        'font.size': 16,  # 默认字体大小
    })
    
    color_palette = ['#dbb428', '#7abbdb', '#bcbd22', '#84ba42', '#17becf', 
                     '#7f7f7f', '#4485c7', '#d4562e', '#a51c36', '#682487']
    
    # 创建2x2的子图布局，固定figsize
    fig = plt.figure(figsize=(11.69, 8.27), dpi=300)
    gs = gridspec.GridSpec(2, 2, figure=fig)
    
    # 存储所有结果
    all_results = {}
    
    # ====== 1. ROC曲线 ======
    ax1 = fig.add_subplot(gs[0, 0])
    roc_results = {}
    
    for idx, score_dict in enumerate(score_dicts):
        name = score_dict['name']
        scores = score_dict['scores']
        y_score = np.array([scores.get(p, 0) for p in all_proteins])
        
        fpr, tpr, _ = roc_curve(y_true, y_score)
        metric = auc(fpr, tpr)
        ax1.plot(fpr, tpr, color=color_palette[idx], linewidth=1.4,  # 调整线条宽度
               label=f'{name} (AUC={metric:.4f})')
        roc_results[name] = metric
    
    ax1.set_xlabel('False Positive Rate\n(a)', fontsize=20)  # 单独设置x轴标签字体大小
    ax1.set_ylabel('True Positive Rate', fontsize=20)  # 单独设置y轴标签字体大小
    ax1.tick_params(axis='x', labelsize=15)  # 调整x轴刻度字体大小
    ax1.tick_params(axis='y', labelsize=15)  # 修改y轴刻度字体大小
    ax1.plot([0, 1], [0, 1], 'k--', linewidth=1.4)  # 调整线条宽度
    ax1.legend(loc='lower right', frameon=True, fontsize=8)  # 调整图例字体大小
    # ax1.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.0)  # 调整网格线宽度
    # ax1.grid(which='minor', linestyle=':', alpha=0.3, linewidth=0.8)  # 调整网格线宽度
    ax1.spines['right'].set_visible(False)  # 去掉右侧的框
    ax1.spines['top'].set_visible(False)  # 去掉上方的框
    
    all_results['roc'] = roc_results
    
    # ====== 2. PR曲线 ======
    ax2 = fig.add_subplot(gs[0, 1])
    pr_results = {}
    
    for idx, score_dict in enumerate(score_dicts):
        name = score_dict['name']
        scores = score_dict['scores']
        y_score = np.array([scores.get(p, 0) for p in all_proteins])
        
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        ap = average_precision_score(y_true, y_score)
        ax2.plot(recall, precision, color=color_palette[idx], linewidth=1.4,  # 调整线条宽度
               label=f'{name} (AP={ap:.4f})')
        pr_results[name] = ap
    
    baseline = len(valid_key)/len(all_proteins)
    ax2.set_xlabel('Recall\n(b)', fontsize=20)  # 单独设置x轴标签字体大小
    ax2.set_ylabel('Precision', fontsize=20)  # 单独设置y轴标签字体大小
    ax2.tick_params(axis='x', labelsize=15)  # 调整x轴刻度字体大小
    ax2.tick_params(axis='y', labelsize=15)  # 修改y轴刻度字体大小
    ax2.axhline(y=baseline, color='k', linestyle='--', linewidth=1.4)  # 调整线条宽度
    ax2.legend(loc='upper right', frameon=True, fontsize=8)  # 调整图例字体大小
    # ax2.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.0)  # 调整网格线宽度
    # ax2.grid(which='minor', linestyle=':', alpha=0.3, linewidth=0.8)  # 调整网格线宽度
    ax2.spines['right'].set_visible(False)  # 去掉右侧的框
    ax2.spines['top'].set_visible(False)  # 去掉上方的框
    
    all_results['pr'] = pr_results
    
    # ====== 3. 排名比较柱状图 ======
    ax3 = fig.add_subplot(gs[1, 0])
    rank_results = {}
    
    # 计算各方法平均排名
    method_names = []
    avg_ranks = []
    colors = []
    
    for idx, score_dict in enumerate(score_dicts):
        name = score_dict['name']
        scores = score_dict['scores']
        
        # 统一排名计算方式
        sorted_proteins = sorted(scores.items(), key=lambda x: x[1], reverse=True)
        ranks = [i+1 for i, (p, _) in enumerate(sorted_proteins) if p in valid_key]
        avg_rank = np.mean(ranks) if len(ranks) > 0 else np.nan
        
        method_names.append(name)
        avg_ranks.append(avg_rank)
        colors.append(color_palette[idx % len(color_palette)])
        rank_results[name] = avg_rank
    
    # 绘制柱状图
    bars = ax3.bar(method_names, avg_ranks, color=colors, alpha=0.7, zorder=2,
                   width=0.6)  # 调整柱状图宽度
    
    # 自动调整Y轴范围
    min_rank = np.nanmin(avg_ranks)
    max_rank = np.nanmax(avg_ranks)
    rank_range = max_rank - min_rank
    ax3.set_ylim(
        max(0, min_rank - 0.1 * rank_range),
        min(max_rank + 0.2 * rank_range, len(network_nodes)))
    
    # 添加数据标签
    for bar in bars:
        height = bar.get_height()
        if not np.isnan(height):
            ax3.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}',
                    ha='center', va='bottom',
                    fontsize=7)  # 调整标签字体大小
    
    ax3.set_xlabel('Methods\n(c)', fontsize=20)  # 调整标签字体大小
    ax3.set_ylabel('Average Rank', fontsize=20)  # 调整标签字体大小
    ax3.tick_params(axis='x', rotation=45, labelsize=15)  # 调整x轴刻度字体大小
    ax3.tick_params(axis='y', labelsize=15)  # 修改y轴刻度字体大小
    # ax3.grid(axis='y', linestyle='--', alpha=0.6, linewidth=1.0, zorder=1)  # 调整网格线宽度
    ax3.spines['right'].set_visible(False)  # 去掉右侧的框
    ax3.spines['top'].set_visible(False)  # 去掉上方的框
    
    all_results['rank'] = rank_results
    
    # ====== 4. 累积关键蛋白曲线 ======
    ax4 = fig.add_subplot(gs[1, 1])
    cumulative_results = {}
    
    # 计算实际要显示的蛋白数量
    total_nodes = len(network_nodes)
    if top_n_range is None:
        top_n = total_nodes
    elif top_n_range <= 1:  # 百分比模式
        top_n = int(total_nodes * top_n_range)
    else:  # 绝对数模式
        top_n = min(int(top_n_range), total_nodes)
    
    step = max(1, int(step))  # 确保步长至少为1
    
    for idx, score_dict in enumerate(score_dicts):
        name = score_dict['name']
        scores = score_dict['scores']
        
        # 按得分降序排序
        sorted_proteins = [p for _, p in sorted(
            zip([scores.get(p, 0) for p in all_proteins], all_proteins),
            reverse=True
        )]
        
        # 动态计算累积关键蛋白
        x_vals = list(range(0, top_n + 1, step))
        y_vals = []
        count = 0
        
        for i, protein in enumerate(sorted_proteins[:top_n]):
            if protein in valid_key:
                count += 1
            if i % step == 0:
                y_vals.append(count)
        
        # 补全最后数据点
        if top_n % step != 0:
            y_vals.append(count)
        
        cumulative_results[name] = {
            'max_found': count,
            'total_key': len(valid_key),
            'ratio': count / len(valid_key)
        }
        
        ax4.plot(x_vals[:len(y_vals)], y_vals,
                color=color_palette[idx], 
                linewidth=1.4,  # 调整线条宽度
                label=f'{name}')
    
    # 图表美化
    ax4.set_xlabel('Top Ranges\n(d)', fontsize=20)  # 调整标签字体大小
    ax4.set_ylabel('Cumulative Key Proteins', fontsize=20)  # 调整标签字体大小
    ax4.tick_params(axis='x', labelsize=15)  # 调整x轴刻度字体大小
    ax4.tick_params(axis='y', labelsize=15)  # 修改y轴刻度字体大小
    ax4.xaxis.set_major_locator(ticker.MultipleLocator(200))
    ax4.xaxis.set_minor_locator(ticker.MultipleLocator(100))
    ax4.yaxis.set_major_locator(ticker.MaxNLocator(7))
    ax4.legend(loc='lower right', frameon=True, fontsize=8)  # 调整图例字体大小，确保不遮挡图像
    # ax4.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.0)  # 调整网格线宽度
    # ax4.grid(which='minor', linestyle=':', alpha=0.3, linewidth=0.8)  # 调整网格线宽度
    ax4.spines['right'].set_visible(False)  # 去掉右侧的框
    ax4.spines['top'].set_visible(False)  # 去掉上方的框
    
    all_results['cumulative'] = cumulative_results
    
    # 调整布局
    plt.tight_layout()
    plt.savefig('酿酒酵母PPI网络CSR_n(n=3,4,...,12)相关结果.png', dpi=300)  # 保存图片
    plt.savefig('酿酒酵母PPI网络CSR_n(n=3,4,...,12)相关结果.pdf', dpi=300)  # 保存图片
    plt.show()
    
    return all_results

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

csr_scores_3 = load_scores('results/csr_scores_3_sc.txt')
csr_scores_4 = load_scores('results/csr_scores_4_sc.txt')
csr_scores_5 = load_scores('results/csr_scores_5_sc.txt')
csr_scores_6 = load_scores('results/csr_scores_6_sc.txt')
csr_scores_7 = load_scores('results/csr_scores_7_sc.txt')
csr_scores_8 = load_scores('results/csr_scores_8_sc.txt')
csr_scores_9 = load_scores('results/csr_scores_9_sc.txt')
csr_scores_10 = load_scores('results/csr_scores_10_sc.txt')
csr_scores_11 = load_scores('results/csr_scores_11_sc.txt')
csr_scores_12 = load_scores('results/csr_scores_12_sc.txt')

# 数据准备
parameter_combinations = [
    {'name': 'CSR_3', 'scores': csr_scores_3},
    {'name': 'CSR_4', 'scores': csr_scores_4},
    {'name': 'CSR_5', 'scores': csr_scores_5},
    {'name': 'CSR_6', 'scores': csr_scores_6},
    {'name': 'CSR_7', 'scores': csr_scores_7},
    {'name': 'CSR_8', 'scores': csr_scores_8},
    {'name': 'CSR_9', 'scores': csr_scores_9},
    {'name': 'CSR_10', 'scores': csr_scores_10},
    {'name': 'CSR_11', 'scores': csr_scores_11},
    {'name': 'CSR_12', 'scores': csr_scores_12}
]

# 调用四合一函数
results = evaluate_parameter_combinations_4in1(
    parameter_combinations, 
    key_proteins,
    non_key_proteins,
    G,
    top_n_range=0.2,
    step=1
)