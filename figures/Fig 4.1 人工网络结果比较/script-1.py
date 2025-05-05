import numpy as np
import matplotlib.pyplot as plt
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 数据准备
node_counts = [100, 200, 300, 400, 500]
categories = [f"{n}" for n in node_counts]

# 人工网络选择不同数量关键节点时，CDR，CSR与其他五种传统方法的AUC、AP、关键节点平均排名得分数据
# 每种指标为一个字典，键为方法名称，值为对应的分数列表（分别对应不同的关键节点数量）
auc_metrics = {
    'DC': [0.8330, 0.7701, 0.7544, 0.7331, 0.7098],
    'CC': [0.8124, 0.7591, 0.7214, 0.6973, 0.6838],
    'BC': [0.7952, 0.7480, 0.7224, 0.7157, 0.6919],
    'EC': [0.8362, 0.7761, 0.7393, 0.7128, 0.6974],
    'PageRank': [0.8274, 0.7664, 0.7505, 0.7313, 0.7071],
    'CDR': [0.8405, 0.7813, 0.7589, 0.7359, 0.7162], 
    'CSR': [0.8402, 0.7798, 0.7587, 0.7368, 0.7161]  
}
ap = {
    'DC': [0.4564, 0.5371, 0.6155, 0.6771, 0.7291],
    'CC': [0.4556, 0.5304, 0.5925, 0.6541,0.7130],
    'BC': [0.4202, 0.5148, 0.5860, 0.6679, 0.7195],
    'EC': [0.4781, 0.5544, 0.6206, 0.6776, 0.7317],
    'PageRank': [0.4552, 0.5429, 0.6214, 0.6865, 0.7387],
    'CDR': [0.4814, 0.5608, 0.6329, 0.6955, 0.7465], 
    'CSR': [0.4865, 0.5552, 0.6314, 0.6948, 0.7452]  
}
avgrank_metrics = {
    'DC': [206.33, 288.94, 327.69, 363.89, 397.34],
    'CC': [219.43, 293.30, 345.56, 382.16, 408.63],
    'BC': [234.86, 302.10, 344.86, 371.08, 404.55],
    'EC': [197.90, 279.65, 333.02, 372.85, 401.82],
    'PageRank': [205.82, 287.36, 325.11, 361.74, 396.97],
    'CDR': [194.04, 275.46, 319.26, 358.95, 392.40],  
    'CSR': [194.35, 276.99, 320.19, 358.71, 392.31] 
}

# 定义颜色
colors = {
    'DC': '#84ba42',
    'BC': '#4485c7',
    'CC': '#7abbdb',
    'EC': '#dbb428',
    'PageRank': '#d4562e',
    'CDR': '#a51c36',
    'CSR': '#682487'
}

# 创建子图
fig, axes = plt.subplots(3, 1, figsize=(18, 13))  
plt.subplots_adjust(hspace=0.1, wspace=0.2)  

bar_width = 0.1
x = np.arange(len(categories))

# 绘制AUC对比图 (子图a)
for i, (metric, values) in enumerate(auc_metrics.items()):
    offset = (i - len(auc_metrics) / 2) * bar_width
    axes[0].bar(x + offset, values, width=bar_width,
                color=colors[metric], edgecolor='black', label=metric, zorder=2)

# 设置坐标刻度和标签
axes[0].set_xticks(x)
axes[0].set_xticklabels(categories, fontsize=20)
axes[0].set_xlabel('# of Selected Key Nodes\n(a)', fontsize=28)
axes[0].set_ylabel('AUC', fontsize=28)

min_auc = min([min(vals) for vals in auc_metrics.values()])
max_auc = max([max(vals) for vals in auc_metrics.values()])
margin_auc = (max_auc - min_auc) * 0.1
axes[0].set_ylim(min_auc - margin_auc, max_auc + margin_auc)

axes[0].legend(loc='upper center', ncol=7, fontsize=20,
               frameon=True, framealpha=0.9, edgecolor='black', bbox_to_anchor=(0.5, 1.30))
axes[0].grid(axis='y', linestyle='--', alpha=0.7, zorder=1)
# axes[0].text(0.5, -0.24, '(a)', transform=axes[0].transAxes, ha='center', fontsize=28)
axes[0].tick_params(axis='y', labelsize=20) 

# 绘制平均精度AP对比图 (子图b)
for i, (metric, values) in enumerate(ap.items()):
    offset = (i - len(ap) / 2) * bar_width
    axes[1].bar(x + offset, values, width=bar_width,
                color=colors[metric], edgecolor='black', label=metric,zorder=2)

# 设置坐标刻度和标签
axes[1].set_xticks(x)
axes[1].set_xticklabels(categories, fontsize=20)
axes[1].set_xlabel('# of Selected Key Nodes\n(b)', fontsize=28)
axes[1].set_ylabel('AP', fontsize=28)

min_ap = min([min(vals) for vals in ap.values()])
max_ap = max([max(vals) for vals in ap.values()])
margin_ap = (max_ap - min_ap) * 0.1
axes[1].set_ylim(min_ap - margin_ap, max_ap + margin_ap)

# axes[1].legend(loc='upper center', ncol=7, fontsize=20,
#                frameon=True, framealpha=0.9, edgecolor='black')
axes[1].grid(axis='y', linestyle='--', alpha=0.7, zorder=1)
# axes[1].text(0.5, -0.24, '(b)', transform=axes[1].transAxes, ha='center', fontsize=28)
axes[1].tick_params(axis='y', labelsize=20) 

# 绘制平均排名对比图 (子图c)
for i, (metric, values) in enumerate(avgrank_metrics.items()):
    offset = (i - len(avgrank_metrics) / 2) * bar_width
    axes[2].bar(x + offset, values, width=bar_width,
                color=colors[metric], edgecolor='black', label=metric, zorder=2)

# 设置坐标刻度和标签
axes[2].set_xticks(x)
axes[2].set_xticklabels(categories, fontsize=20)
axes[2].set_xlabel('# of Selected Key Nodes\n(c)', fontsize=28)
axes[2].set_ylabel('Average Rank', fontsize=28)

all_avg_rank_values = [val for sublist in avgrank_metrics.values() for val in sublist]
min_avg_rank = min(all_avg_rank_values)
max_avg_rank = max(all_avg_rank_values)
margin = (max_avg_rank - min_avg_rank) * 0.1
axes[2].set_ylim(min_avg_rank - margin, max_avg_rank + margin)

# axes[2].legend(loc='upper center', ncol=7, fontsize=20,
#                frameon=True, framealpha=0.9, edgecolor='black')
axes[2].grid(axis='y', linestyle='--', alpha=0.7, zorder=1)
# axes[2].text(0.5, -0.24, '(c)', transform=axes[2].transAxes, ha='center', fontsize=28)
axes[2].tick_params(axis='y', labelsize=20) 

# 只显示有刻度的x轴和y轴
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

plt.tight_layout()
# 设置分辨率为300 dpi并保存图片
plt.savefig('人工网络结果.png', dpi=300)
plt.savefig('人工网络结果.pdf', dpi=300)
plt.show()