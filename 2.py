
import matplotlib.pyplot as plt
from dragon_positon import *

points_x,points_y,point_p_x,point_p_y ,crash_head_theta= dragon_position(410, 415,0.55 / (2 * np.pi),True,0.01)


# 绘制标记的点
if len(point_p_x) > 0:
    plt.scatter(point_p_x[-1], point_p_y[-1], color='red', label='Point Position')
plt.scatter(points_x, points_y, color='green', label='Meter Position')
plt.plot(points_x, points_y, color='yellow', label='Meter Position')
# 添加标题和标签
plt.title('Archimedean Spiral with Meter Positions Marked')
plt.xlabel('X Axis')
plt.ylabel('Y Axis')
# 添加网格
plt.grid()
plt.axis('equal')
plt.legend()
plt.show()







