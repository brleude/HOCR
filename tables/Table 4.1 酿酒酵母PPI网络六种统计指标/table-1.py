import pandas as pd
import sys
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

# 加载得分数据
dc_scores = load_scores('results/dc_scores_sc.txt')
cc_scores = load_scores('results/cc_scores_sc.txt')
bc_scores = load_scores('results/bc_scores_sc.txt')
ec_scores = load_scores('results/ec_scores_sc.txt')
pr_scores = load_scores('results/pr_scores_sc.txt')
cdr_scores = load_scores('results/cdr_scores_sc.txt')
csr_scores = load_scores('results/csr_scores_sc.txt')

# 统计指标计算
top_n = 500
all_rankings = {
    'DC': dc_scores,
    'BC': bc_scores,
    'CC': cc_scores,
    'EC': ec_scores,
    'PageRank': pr_scores,
    'CDR': cdr_scores,
    'CSR': csr_scores,
}
all_results = calculate_all_metrics(all_rankings, key_proteins, non_key_proteins, top_n=top_n)

# 将数据转换为表格并保存为 CSV 文件
csv_file_path = '酿酒酵母PPI网络六种统计指标.csv'
table = data_to_table(all_results, csv_file_path, file_format='csv')
print(table)