
from dragon_positon import *
currect_head_theta =1000
crash_head_theta = 0
#遍历螺线间隔使得与圆相交的点刚好碰撞
new_b= 0.4504 / (2 * np.pi)
while(crash_head_theta <= currect_head_theta):
    new_b-=0.0002/(2 * np.pi)
    currect_head_theta,currect_head_x,currect_head_y=find_intersection(new_b)
    print(currect_head_x,currect_head_y)
    points_x,points_y,point_p_x,point_p_y,crash_head_theta = dragon_position(200, 260,new_b,True,0.01)
    print(crash_head_theta/(2 * np.pi),currect_head_theta/(2 * np.pi))
    print()

print(new_b*2*np.pi)

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







