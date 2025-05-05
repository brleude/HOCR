import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 读取图片
p1 = mpimg.imread('Fig 4.8(a) 酿酒酵母PPI网络关键节点和非关键节点模体参与计数分布/酿酒酵母PPI网络关键节点和非关键节点模体参与计数分布.png')
p2 = mpimg.imread('Fig 4.8(b) 大肠杆菌PPI网络关键节点和非关键节点模体参与计数分布/大肠杆菌PPI网络关键节点和非关键节点模体参与计数分布.png')

# 创建一个 2x1 的子图布局，使用 gridspec 来控制布局
fig = plt.figure(figsize=(11.69, 4 * 2), dpi=300)  # 调整总大小，宽度为 11.69，高度为两张 4 高的图
gs = plt.GridSpec(2, 1, width_ratios=[1], height_ratios=[4, 4])

# 绘制 p1 在上方
ax1 = fig.add_subplot(gs[0, 0])
ax1.imshow(p1)
ax1.axis('off')
# 在 p1 下方居中添加编号 (a)
ax1.text(0.5, -0.05, '(a)', transform=ax1.transAxes, ha='center', va='center', fontsize=18)

# 绘制 p2 在下方
ax2 = fig.add_subplot(gs[1, 0])
ax2.imshow(p2)
ax2.axis('off')
# 在 p2 下方居中添加编号 (b)
ax2.text(0.5, -0.05, '(b)', transform=ax2.transAxes, ha='center', va='center', fontsize=18)

# 调整子图之间的间距
plt.tight_layout()
# 保存图片
plt.savefig('模体参与计数分布整合.png', dpi=300) 
plt.savefig('模体参与计数分布整合.pdf', dpi=300)  
# 显示图片
plt.show()