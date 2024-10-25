import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.optimize as spo
from matplotlib.animation import FuncAnimation

# 设置常数
b = 0.55 / (2 * np.pi)
theta_start = 64 * np.pi  # 起点为16圈
theta_head = 32 * np.pi
theta_end = 0  # 终点为源点
bench_lengths = [2.86] + [1.65] * 222  # 每节板凳的长度
points_x = []
points_y = []

# 做螺线图
theta = np.linspace(theta_end, theta_start, 20000)  # 生成点模拟曲线
r = b * theta  # 阿基米德螺线方程

# 定义阿基米德螺线的导数
def integrand(theta):
    r_value = b * theta  # 当前的 r
    dr_dtheta = b  # dr/dtheta
    return np.sqrt(dr_dtheta**2 + r_value**2)

# 计算弧长函数
def arc_length(theta):
    return spi.quad(integrand, theta, theta_head)[0]  # 从0到theta的弧长

# 动画初始化函数
fig, ax = plt.subplots()

def init():
    ax.clear()
    ax.plot(r * np.cos(theta), r * np.sin(theta), label='Archimedean Spiral')
    ax.set_title('Archimedean Spiral with Meter Positions Marked')
    ax.set_xlabel('X Axis')
    ax.set_ylabel('Y Axis')
    ax.grid()
    ax.axis('equal')
    ax.legend()

# 更新函数：每次调用更新一帧
def update(item):
    points_x.clear()
    points_y.clear()
    target_length = item

    # 寻找使得 arc_length(theta) = target_length 的 θ
    theta_solution = spo.brentq(lambda x: arc_length(x) - target_length, theta_end, theta_head)

    # 计算对应的 r 和 (x, y) 坐标
    r_solution = b * theta_solution
    x_solution = r_solution * np.cos(theta_solution)
    y_solution = r_solution * np.sin(theta_solution)
    points_x.append(x_solution)
    points_y.append(y_solution)

    for length in bench_lengths:
        # 将圆的方程转化为一个函数
        def circle_equation(theta_val):
            r_ed = b * theta_val
            x_eq = r_ed * np.cos(theta_val)
            y_eq = r_ed * np.sin(theta_val)
            return np.abs((x_eq - x_solution)**2 + (y_eq - y_solution)**2 - length**2)

        # 寻找解，使用 fsolve
        theta_guess = theta_solution + 0.3 * np.pi  # 使用当前解作为初始猜测
        theta_solutions = spo.fsolve(circle_equation, theta_guess)

        for sol in theta_solutions:
            current_r = b * sol
            current_theta = sol

            if current_theta > theta_solution and (current_theta - theta_solution) < 0.5 * np.pi:
                r_solution = current_r
                x_solution = r_solution * np.cos(current_theta)
                y_solution = r_solution * np.sin(current_theta)
                points_x.append(x_solution)
                points_y.append(y_solution)

                theta_solution = current_theta
                break

    # 每帧更新图中的标记点和连线
    ax.plot(r * np.cos(theta), r * np.sin(theta), label='Archimedean Spiral')
    ax.scatter(points_x, points_y, color='green', label='Meter Position')
    ax.plot(points_x, points_y, color='yellow', label='Meter Position')

# 创建动画对象
ani = FuncAnimation(fig, update, frames=np.arange(100, 300), init_func=init, interval=500, repeat=False)

# 显示动画
plt.show()

# 保存为 GIF 动图文件（可选）
# ani.save('spiral_animation.gif', writer='imagemagick')
