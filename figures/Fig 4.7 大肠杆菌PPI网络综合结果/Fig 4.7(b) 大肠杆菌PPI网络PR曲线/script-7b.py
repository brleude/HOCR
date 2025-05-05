import numpy as np
import json
import csv
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import ticker
from sklearn.metrics import precision_recall_curve, average_precision_score
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import load_scores
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 绘制各方法的PR曲线
def plot_precision_recall_curves(key_proteins, non_key_proteins, dicts_list, labels):
    """
    PR曲线对比

    :param key_proteins: 已验证关键蛋白列表（正样本）
    :param non_key_proteins: 已验证非关键蛋白列表（负样本）
    :param dicts_list: 方法得分字典列表 [{蛋白:得分}, ...]
    :param labels: 方法名称列表

    :return: None
    """
    # ====== 初始化检查 ======
    # 检查正负样本重叠
    overlap = set(key_proteins) & set(non_key_proteins)
    if overlap:
        raise ValueError(f"正负样本存在{len(overlap)}个重叠蛋白，首5个：{list(overlap)[:5]}")

    # ====== 可视化初始化 ======
    plt.figure(figsize=(11.69, 7), dpi=300)  # 图片大小和分辨率设置
    ax = plt.gca()

    # ====== 颜色编码 ======
    color_palette = {
        'DC': '#84ba42',
        'BC': '#4485c7',
        'CC': '#7abbdb',
        'EC': '#dbb428',
        'PageRank': '#d4562e',
        'CDR_5': '#a51c36',
        'CSR_3': '#682487'
    }  # 迁移颜色设置

    # ====== 遍历各方法 ======
    for idx, (score_dict, method_label) in enumerate(zip(dicts_list, labels)):
        # --- 数据过滤 ---
        valid_scores, valid_labels = [], []
        for protein, score in score_dict.items():
            if protein in key_proteins:
                valid_scores.append(score)
                valid_labels.append(1)
            elif protein in non_key_proteins:
                valid_scores.append(score)
                valid_labels.append(0)

        # --- 数据验证 ---
        if not valid_scores:
            raise ValueError(f"'{method_label}'无有效数据（未匹配任何正负样本）")

        y_true = np.array(valid_labels)
        if len(np.unique(y_true)) < 2:
            class_counts = np.bincount(y_true)
            raise ValueError(f"'{method_label}'样本不均衡：正样本={class_counts[1]}, 负样本={class_counts[0]}")

        # --- 指标计算 ---
        precision, recall, _ = precision_recall_curve(y_true, valid_scores)
        ap_score = average_precision_score(y_true, valid_scores)

        # --- 可视化配置 ---
        line_config = {
            'color': color_palette.get(method_label, f'C{idx}'),
            'linewidth': 2.4  # 线条宽度设置
        }

        # --- 绘制曲线 ---
        plt.plot(
            recall,
            precision,
            **line_config,
            label=f'{method_label} (AP={ap_score:.4f})'
        )

    # ====== 图表优化 ======
    # 坐标轴
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))

    # 标签与标题
    plt.xlabel('Recall', fontsize=40)  # x轴标签字体大小设置
    plt.ylabel('Precision', fontsize=40)  # y轴标签字体大小设置
    # plt.title('Precision-Recall Analysis',
    #          fontsize=15,
    #          pad=18,
    #          fontweight='bold',
    #          )

    # 设置刻度字体大小
    ax.tick_params(axis='x', labelsize=36)  # x轴刻度字体大小设置
    ax.tick_params(axis='y', labelsize=36)  # y轴刻度字体大小设置

    # # 网格与边框
    # ax.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.4)  # 网格线宽度设置
    # ax.grid(which='minor', linestyle=':', alpha=0.3, linewidth=1.4)  # 网格线宽度设置
    # for spine in ax.spines.values():
    #     spine.set_linewidth(1.4)  # 边框线宽度设置
    #     spine.set_color('#333333')

    # 去掉右方和上方的框
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # 图例
    legend = plt.legend(
        loc='upper right',
        frameon=True,
        edgecolor='black',
        fancybox=False,
        fontsize=19.5,  # 图例字体大小设置
        borderpad=0.1
    )
    legend.get_frame().set_facecolor('#FFFFFF')
    legend.get_frame().set_alpha(0.9) 

    plt.tight_layout()
    plt.savefig('大肠杆菌PPI网络PR曲线.png', dpi=300)  # 保存图片
    plt.savefig('大肠杆菌PPI网络PR曲线.pdf', dpi=300)  # 保存图片
    plt.show()


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

# 加载待比较得分
dc_scores_ec = load_scores('results/dc_scores_ec.txt')
bc_scores_ec = load_scores('results/bc_scores_ec.txt')
cc_scores_ec = load_scores('results/cc_scores_ec.txt')
ec_scores_ec = load_scores('results/ec_scores_ec.txt')
pr_scores_ec = load_scores('results/pr_scores_ec.txt')
cdr_scores_ec = load_scores('results/cdr_scores_ec.txt')
csr_scores_ec = load_scores('results/csr_scores_ec.txt')

# 准备数据
dicts_list_ec = [dc_scores_ec, bc_scores_ec, cc_scores_ec, ec_scores_ec, pr_scores_ec, cdr_scores_ec, csr_scores_ec]  # 包含每个算法的中心性值字典
labels = ["DC", "BC", "CC", "EC", "PageRank", "CDR_5", "CSR_3"]  # 每条曲线的标签

# 绘制 PR 曲线
plot_precision_recall_curves(key_proteins_ec, non_key_proteins_ec, dicts_list_ec, labels)