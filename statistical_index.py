import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 选择排名前 top_n 个蛋白质
def select_top_proteins(dicts, top_n=20):
    """
    根据中心性值排序并选择前 top_n 个蛋白质
    :param dicts: 字典，键为蛋白质名称，值为中心性值
    :param top_n: 选择前 top_n 个蛋白质
    :return: 前 top_n 个蛋白质列表
    """
    # 按分值降序排序
    sorted_proteins = sorted(dicts.keys(), key=lambda x: dicts[x], reverse=True)
    
    # 选择前 top_n 个蛋白质
    top_proteins = sorted_proteins[:top_n]
    
    return top_proteins
# 计算混淆矩阵
def calculate_confusion_matrix(top_proteins, key_proteins, non_key_proteins):
    """
    计算混淆矩阵
    :param top_proteins: 前 top_n 个蛋白质列表
    :param key_proteins: 关键蛋白列表
    :param non_key_proteins: 非关键蛋白列表
    :return: TP, FP, TN, FN
    """
    TP = len([protein for protein in top_proteins if protein in key_proteins])  # 真阳性
    FP = len([protein for protein in top_proteins if protein in non_key_proteins])  # 假阳性
    TN = len([protein for protein in non_key_proteins if protein not in top_proteins])  # 真阴性
    FN = len([protein for protein in key_proteins if protein not in top_proteins])  # 假阴性
    
    return TP, FP, TN, FN

# 计算六个统计指标（SN, SP, PPV, NPV, F-score, ACC）
def calculate_metrics(TP, FP, TN, FN):
    """
    计算敏感性、特异性、阳性预测值、阴性预测值、F 值、准确率
    :param TP: 真阳性
    :param FP: 假阳性
    :param TN: 真阴性
    :param FN: 假阴性
    :return: SN, SP, PPV, NPV, F, ACC
    """
    SN = TP / (TP + FN) if (TP + FN) != 0 else 0  # 敏感性
    SP = TN / (TN + FP) if (TN + FP) != 0 else 0  # 特异性
    PPV = TP / (TP + FP) if (TP + FP) != 0 else 0  # 阳性预测值
    NPV = TN / (TN + FN) if (TN + FN) != 0 else 0  # 阴性预测值
    F = 2 * (PPV * SN) / (PPV + SN) if (PPV + SN) != 0 else 0  # F 值
    ACC = (TP + TN) / (TP + FP + TN + FN) if (TP + FP + TN + FN) != 0 else 0  # 准确率
    
    return SN, SP, PPV, NPV, F, ACC

# 整合函数，计算各个方法的六个评价指标得分
def calculate_all_metrics(all_rankings, key_proteins, non_key_proteins, top_n=500):
    """
    整合函数，计算各个方法的六个评价指标得分
    :param all_rankings: 包含各个方法得分字典的字典，键为方法名，值为得分字典
    :param key_proteins: 关键蛋白列表
    :param non_key_proteins: 非关键蛋白列表
    :param top_n: 选择前 top_n 个蛋白质
    :return: 包含各个方法六个评价指标得分的字典
    """
    result = {}
    for method, scores in all_rankings.items():
        top_proteins = select_top_proteins(scores, top_n)
        TP, FP, TN, FN = calculate_confusion_matrix(top_proteins, key_proteins, non_key_proteins)
        SN, SP, PPV, NPV, F, ACC = calculate_metrics(TP, FP, TN, FN)
        result[method] = {
            'SN': SN,
            'SP': SP,
            'PPV': PPV,
            'NPV': NPV,
            'F': F,
            'ACC': ACC
        }
    return result