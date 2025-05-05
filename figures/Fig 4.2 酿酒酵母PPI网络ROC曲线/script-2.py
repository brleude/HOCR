import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from matplotlib import ticker
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


# 绘制各方法的ROC曲线
def plot_roc_curves(key_proteins, non_key_proteins, dicts_list, labels):
    """
    ROC曲线对比

    :param key_proteins: 关键蛋白列表（正样本，ground truth）
    :param non_key_proteins: 非关键蛋白列表（负样本，必须与正样本互斥）
    :param dicts_list: 多个评分字典列表，每个字典形如{protein: score}
    :param labels: 方法标签列表，与dicts_list顺序一致

    :return: None
    """
    # ====== 图形初始化 ======
    plt.figure(figsize=(11.69, 7), dpi=300)
    ax = plt.gca()

    # ====== 颜色方案 ======
    color_palette = {
            'DC': '#84ba42',
            'BC': '#4485c7',
            'CC': '#7abbdb',
            'EC': '#dbb428',
            'PageRank': '#d4562e',
            'CDR': '#a51c36',
            'CSR': '#682487'
    }

    # ====== 数据验证 ======
    # 检查正负样本重叠
    overlap = set(key_proteins) & set(non_key_proteins)
    if overlap:
        raise ValueError(f"正负样本存在重叠蛋白: {overlap}")

    # ====== 绘制曲线 ======
    for i, (score_dict, label) in enumerate(zip(dicts_list, labels)):
        # --- 数据过滤 ---
        valid_scores = []
        valid_labels = []

        for protein, score in score_dict.items():
            if protein in key_proteins:
                valid_scores.append(score)
                valid_labels.append(1)
            elif protein in non_key_proteins:
                valid_scores.append(score)
                valid_labels.append(0)
            # 忽略不在两个列表中的蛋白

        # --- 数据检查 ---
        if len(valid_labels) == 0:
            raise ValueError(f"'{label}'方法无有效样本")

        y_true = np.array(valid_labels)
        if len(np.unique(y_true)) < 2:
            class_ratio = np.bincount(y_true)
            raise ValueError(
                f"'{label}'方法样本不均衡: 正样本={class_ratio[1]}, 负样本={class_ratio[0]}"
            )

        # --- 计算指标 ---
        fpr, tpr, _ = roc_curve(y_true, valid_scores)
        roc_auc = auc(fpr, tpr)

        # --- 可视化配置 ---
        line_color = color_palette.get(label, f'C{i}')
        line_style = '-'

        plt.plot(
            fpr, tpr,
            color=line_color,
            linestyle=line_style,
            linewidth=1.4,  
            label=f'{label} (AUC={roc_auc:.4f})'
        )

    # ====== 参考线 ======
    plt.plot([0, 1], [0, 1],
             color='#666666',
             linewidth=1.4,  
             linestyle='-.',
             alpha=0.6)

    # ====== 学术图表优化 ======
    # 坐标轴刻度
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))

    # 字体与标签
    plt.xlabel('False Positive Rate', fontsize=20)  # 修改x轴标签字体大小
    plt.ylabel('True Positive Rate', fontsize=20)  # 修改y轴标签字体大小
    # plt.title('Comparative ROC Analysis',
    #          fontsize=15,
    #          pad=18,
    #          fontweight='bold')

    # 设置刻度字体大小
    ax.tick_params(axis='x', labelsize=18)  # 修改x轴刻度字体大小
    ax.tick_params(axis='y', labelsize=18)  # 修改y轴刻度字体大小

    # 网格与边框
    # ax.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.4)
    # ax.grid(which='minor', linestyle=':', alpha=0.3, linewidth=1.4)
    # for spine in ax.spines.values():
    #     spine.set_linewidth(1.4)
    #     spine.set_color('#333333')

    # 去掉右方和上方的框
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # 图例设置
    legend = plt.legend(
        loc='lower right',
        frameon=True,
        edgecolor='black',
        fancybox=False,
        fontsize=18,  # 设置图例字体大小
        borderpad=0.8
    )
    legend.get_frame().set_alpha(0.9)

    # 边框优化
    for spine in ax.spines.values():
        spine.set_linewidth(1.4)

    plt.tight_layout()

    # 设置分辨率为300 dpi并保存图片
    plt.savefig('酿酒酵母PPI网络ROC曲线.png', dpi=300)
    plt.savefig('酿酒酵母PPI网络ROC曲线.pdf', dpi=300)

    plt.show()


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
# 读取各方法得分
dc_scores = load_scores('results/dc_scores_sc.txt')
cc_scores = load_scores('results/cc_scores_sc.txt')
bc_scores = load_scores('results/bc_scores_sc.txt')
ec_scores = load_scores('results/ec_scores_sc.txt')
pr_scores = load_scores('results/pr_scores_sc.txt')
cdr_scores = load_scores('results/cdr_scores_sc.txt')
csr_scores = load_scores('results/csr_scores_sc.txt')

dicts_list = [dc_scores, bc_scores, cc_scores, ec_scores, pr_scores, cdr_scores, csr_scores]  # 包含每个算法的中心性值字典
labels = ["DC", "BC", "CC", "EC", "PageRank", "CDR", "CSR"]  # 每条曲线的标签

# 绘制 ROC 曲线
plot_roc_curves(key_proteins, non_key_proteins, dicts_list, labels)
