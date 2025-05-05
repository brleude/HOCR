import math
import json
import csv
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import ticker
import sys
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import load_scores
from statistical_index import select_top_proteins
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


# 绘制各方法的Jackknife 方法关键蛋白质累加数量曲线
def cumulative_performance_visualization(all_rankings, key_proteins, max_range):
    # 定义排名范围
    ranking_ranges = list(range(0, max_range + 1, 1))

    # 累积计算函数
    def compute_cumulative_tp(sorted_proteins, key_proteins):
        """支持动态间隔的累积计算"""
        cumulative = []
        count = 0
        current_index = 0  # 当前处理的蛋白质索引

        # 预先生成有效排序列表（处理不足max_range的情况）
        valid_sorted = sorted_proteins[:max_range] if len(sorted_proteins) >= max_range else sorted_proteins

        for target_rank in ranking_ranges:
            # 处理超出实际数据范围的情况
            while current_index < target_rank:
                if current_index < len(valid_sorted):
                    protein = valid_sorted[current_index]
                    if protein in key_proteins:
                        count += 1
                current_index += 1
            cumulative.append(count)
        return cumulative

    # 生成足够长的排序列表
    method_rankings = {
        method: select_top_proteins(scores, max_range)
        for method, scores in all_rankings.items()
    }

    # 计算累积结果
    cumulative_results = {}
    for method, ranking in method_rankings.items():
        cumulative_results[method] = compute_cumulative_tp(ranking, key_proteins)

    # 优化可视化设置
    plt.figure(figsize=(11.69, 7), dpi=300)  # 修改图片大小和分辨率

    # 颜色方案
    color_palette = {
        'DC': '#84ba42',
        'BC': '#4485c7',
        'CC': '#7abbdb',
        'EC': '#dbb428',
        'PageRank': '#d4562e',
        'CDR_5': '#a51c36',
        'CSR_3': '#682487'
    }

    # 绘制曲线
    for method in method_rankings.keys():
        plt.plot(
            ranking_ranges,
            cumulative_results[method],
            markersize=4,
            linewidth=2.4,  # 修改线条宽度
            color=color_palette[method],
            label=method
        )

    # 优化坐标轴
    ax = plt.gca()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100))  # 主刻度每50单位
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(50))  # 次要刻度每10单位
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))  # Y轴自动刻度

    # 添加辅助网格
    # ax.grid(which='major', linestyle='--', alpha=0.6, linewidth=1.4)  # 修改网格线宽度
    # ax.grid(which='minor', linestyle=':', alpha=0.3, linewidth=1.4)  # 修改网格线宽度

    # 调整标签
    plt.xlabel("Top Ranges", fontsize=40, labelpad=8)  # 修改x轴标签字体大小
    plt.ylabel("Cumulative Key\n Proteins", fontsize=40, labelpad=8)  #修改y轴标签字体大小
    # plt.title("Key Protein Prediction Performance", 
    #          fontsize=14, pad=15, fontweight='bold')

    # 设置刻度字体大小
    ax.tick_params(axis='x', labelsize=36)  # 修改x轴刻度字体大小
    ax.tick_params(axis='y', labelsize=36)  # 修改y轴刻度字体大小

    # 去掉右方和上方的框
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # 紧凑图例
    legend = plt.legend(
        loc='lower right',
        ncol=1,  # 一列布局
        fontsize=19.5  # 修改图例字体大小
    )
    legend.get_frame().set_linewidth(1.4)  # 修改图例边框线宽度

    plt.tight_layout()
    plt.savefig('大肠杆菌PPI网络Jackknife关键蛋白质累加数量曲线.png', dpi=300)  # 保存图片
    plt.savefig('大肠杆菌PPI网络Jackknife关键蛋白质累加数量曲线.pdf', dpi=300)  # 保存图片
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

# 整合得分
all_rankings_ec = {
    'DC': dc_scores_ec,
    'BC': bc_scores_ec,
    'CC': cc_scores_ec,
    'EC': ec_scores_ec,
    'PageRank': pr_scores_ec,
    'CDR_5': cdr_scores_ec,
    'CSR_3': csr_scores_ec,
}
max_range_ec = math.floor(len(G_ec.nodes)/2) # 进行比较的排名范围上界
# 绘制Jackknife 方法关键蛋白质累加数量曲线
cumulative_performance_visualization(all_rankings_ec, key_proteins_ec, max_range_ec)