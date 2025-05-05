import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)


# 读取图片
p1 = mpimg.imread('Fig 4.7(a) 大肠杆菌PPI网络ROC曲线/大肠杆菌PPI网络ROC曲线.png')
p2 = mpimg.imread('Fig 4.7(b) 大肠杆菌PPI网络PR曲线/大肠杆菌PPI网络PR曲线.png')
p3 = mpimg.imread('Fig 4.7(c) 大肠杆菌PPI网络Jackknife关键蛋白质累加数量曲线/大肠杆菌PPI网络Jackknife关键蛋白质累加数量曲线.png')
p4 = mpimg.imread('Fig 4.7(d) 大肠杆菌PPI网络关键蛋白质排名分布/大肠杆菌PPI网络关键蛋白质排名分布.png')

# 创建一个 3x2 的子图布局，使用 gridspec 来控制布局
fig = plt.figure(figsize=(11.69 * 2, 7 * 3), dpi=300)  # 调整总大小，宽度为两张 11.69 宽的图，高度为三张 7 高的图
gs = plt.GridSpec(3, 2, width_ratios=[1, 1], height_ratios=[7, 7, 7])

# 绘制 p1、p2、p3 在左列
ax1 = fig.add_subplot(gs[0, 0])
ax1.imshow(p1)
ax1.axis('off')
# 在 p1 下方居中添加编号 (a)
ax1.text(0.5, -0.05, '(a)', transform=ax1.transAxes, ha='center', va='center', fontsize=40)

ax2 = fig.add_subplot(gs[1, 0])
ax2.imshow(p2)
ax2.axis('off')
# 在 p2 下方居中添加编号 (b)
ax2.text(0.5, -0.05, '(b)', transform=ax2.transAxes, ha='center', va='center', fontsize=40)

ax3 = fig.add_subplot(gs[2, 0])
ax3.imshow(p3)
ax3.axis('off')
# 在 p3 下方居中添加编号 (c)
ax3.text(0.5, -0.05, '(c)', transform=ax3.transAxes, ha='center', va='center', fontsize=40)

# 绘制 p4 在右列，跨越三行
ax4 = fig.add_subplot(gs[:, 1])
ax4.imshow(p4)
ax4.axis('off')
# 在 p4 下方居中添加编号 (d)
ax4.text(0.5, -0.02, '(d)', transform=ax4.transAxes, ha='center', va='center', fontsize=40)

# 调整子图之间的间距
plt.tight_layout()

# 保存图片
plt.savefig('大肠杆菌PPI网络综合结果.png', dpi=300)
plt.savefig('大肠杆菌PPI网络综合结果.pdf', dpi=300)

# 显示图片
plt.show()