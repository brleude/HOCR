import pandas as pd
import sys
import csv
import json
from pathlib import Path
python_code_path = str(Path(__file__).parent.parent.parent)
if python_code_path not in sys.path:
    sys.path.append(python_code_path)
from save_load_scores import load_scores
from statistical_index import calculate_all_metrics
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 将数据保存为表格
def data_to_table(data, file_path, file_format='csv'):
    """
    将字典数据转换为 pandas DataFrame 表格并保存到文件
    :param data: 包含方法和指标数据的字典
    :param file_path: 保存文件的路径
    :param file_format: 文件格式，默认为 'csv'，也支持 'excel'
    :return: 转换后的 pandas DataFrame
    """
    df = pd.DataFrame.from_dict(data, orient='index')
    if file_format == 'csv':
        df.to_csv(file_path, index_label='Method')
        print("数据已保存")
    elif file_format == 'excel':
        df.to_excel(file_path, index_label='Method')
        print("数据已保存")
    else:
        raise ValueError("不支持的文件格式，仅支持 'csv' 和 'excel'。")
    return df

# 提取原始PPI网络
current_file = Path(__file__).resolve()
parent3 = current_file.parent.parent.parent
ec_ppi_path = parent3/'data/ec/Ec PPI Network.csv'
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
ec_key_genes = extract_gene_names_from_file(parent3/"data/ec/Ec Key Proteins.tsv")

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

# 加载待比较得分
dc_scores_ec = load_scores('results/dc_scores_ec.txt')
bc_scores_ec = load_scores('results/bc_scores_ec.txt')
cc_scores_ec = load_scores('results/cc_scores_ec.txt')
ec_scores_ec = load_scores('results/ec_scores_ec.txt')
pr_scores_ec = load_scores('results/pr_scores_ec.txt')
cdr_scores_ec = load_scores('results/cdr_scores_ec.txt')
csr_scores_ec = load_scores('results/csr_scores_ec.txt')

# 统计指标计算
top_m = 200
all_rankings_ec = {
    'DC': dc_scores_ec,
    'BC': bc_scores_ec,
    'CC': cc_scores_ec,
    'EC': ec_scores_ec,
    'PageRank': pr_scores_ec,
    'CDR_5': cdr_scores_ec,
    'CSR_3': csr_scores_ec,
}

all_results_ec = calculate_all_metrics(all_rankings_ec, key_proteins_ec, non_key_proteins_ec, top_m)

# 将数据转换为表格并保存为 CSV 文件
csv_file_path = '大肠杆菌PPI网络六种统计指标.csv'
table = data_to_table(all_results_ec, csv_file_path, file_format='csv')
print(table)