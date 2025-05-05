import argparse
import os
import ast
import networkx as nx
from centrality_improvement import improved_centrality
from bayesian_optimization import optimize_method

def parse_arguments():
    """
    解析命令行参数
    """
    parser = argparse.ArgumentParser(description="网络分析工具")
    # 定义必需的参数
    required_args = parser.add_argument_group('必需的参数')
    required_args.add_argument('-G', '--network', type=str, help="网络的 edgelist 文件路径", required=True)
    def base_scores_type(value):
        valid_strings = ['bc', 'cc', 'dc', 'ec', 'pr']
        if value in valid_strings:
            return value
        try:
            return ast.literal_eval(value)
        except (SyntaxError, ValueError):
            raise argparse.ArgumentTypeError(f"base_scores 必须是 'bc', 'cc', 'dc', 'ec', 'pr' 之一或有效的字典字符串")
    required_args.add_argument('-b', '--base_scores', type=base_scores_type,
                               help="基础中心性得分，可选值为 'bc', 'cc', 'dc', 'ec', 'pr' 或得分字典，如 {'a':0.12, 'b':0.43}",
                               required=True)
    required_args.add_argument('-o', '--output', type=str, help="输出文件夹路径", required=True)

    # 定义可选参数
    optional_args = parser.add_argument_group('可选的参数')
    def read_nodes_from_txt(file_path):
        try:
            with open(file_path, 'r') as f:
                content = f.read().strip()
                if content:
                    return [int(item) for item in content.split(',')]
                return []
        except FileNotFoundError:
            raise argparse.ArgumentTypeError(f"文件 {file_path} 未找到")

    def node_list_type(value):
        try:
            return [int(item) for item in value.split(',')]
        except ValueError:
            return read_nodes_from_txt(value)
    optional_args.add_argument('-k', '--key_nodes', type=node_list_type,
                               help="关键节点列表，用逗号分隔（例如：1,2,3）或 txt 文件路径", default=None)
    optional_args.add_argument('-n', '--non_key_nodes', type=node_list_type,
                               help="非关键节点列表，用逗号分隔（例如：4,5,6）或 txt 文件路径", default=None)
    required_args.add_argument('-m', '--max_clique', type=int, help="所考虑的最大团的大小，若为None则自动识别网络中最大的团", required=False)
    required_args.add_argument('-r', '--rank_type', type=str, help="选择对CDR或CSR进行调参(输入cdr或csr)，若为None则默认为csr", required=False)
    optional_args.add_argument('-v', '--verbose', type=lambda x: (str(x).lower() == 'true'), help="是否输出调参日志 (True 或 False)", 
                               default=False)

    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = parse_arguments()
    output_folder = args.output
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_file_path = os.path.join(output_folder, 'output.txt')
    G = nx.read_edgelist(args.network)
    if args.key_nodes:
        best_params, _ = optimize_method(G, base_scores=args.base_scores, key_nodes=args.key_nodes, 
                                         non_key_nodes=args.non_key_nodes, max_clique=args.max_clique,
                                         rank_type=args.rank_type, verbose=args.verbose)
        cdr1, csr1 = improved_centrality(G, base_scores=args.base_scores, max_clique=args.max_clique, 
                                         params=best_params)
        # 将结果保存到txt文件
        with open(output_file_path, 'w') as f:
            f.write(f"cdr: {cdr1}\n")
            f.write(f"csr: {csr1}\n")
        
    else:
        cdr2, csr2 = improved_centrality(G, base_scores=args.base_scores, max_clique=args.max_clique,
                                         params=None)
        # 将结果保存到txt文件
        with open(output_file_path, 'w') as f:
            f.write(f"cdr: {cdr2}\n")
            f.write(f"csr: {csr2}\n")