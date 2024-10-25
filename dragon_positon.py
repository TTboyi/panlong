import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as spi
import scipy.optimize as spo
import math
import sympy as sp  # 确保导入 sympy

# 设置常数

theta_start = 48 * np.pi  # 起点为16圈
theta_head = 32 * np.pi
theta_end = 0  # 终点为源点
bench_lengths = [2.86] + [1.65] * 222  # 每节板凳的长度
points_x = []
points_y = []
point_p_x = []
point_p_y = []
head_theta=0
m = 0.275**2 + 0.15**2
n = 3.135**2 + 0.15**2
# 做螺线图
theta = np.linspace(theta_end, theta_start, 20000)  # 生成点模拟曲线
# r = b * theta  # 阿基米德螺线方程

def change_b(new_b =  0.55 / (2 * np.pi)):
   global b
   b= new_b


# 定义阿基米德螺线的导数
def integrand(theta, b):
    r_value = b * theta  # 当前的 r
    dr_dtheta = b  # dr/dtheta
    return np.sqrt(dr_dtheta**2 + r_value**2)

# 计算弧长函数，接受 b 作为参数
def arc_length(theta, b):
    return spi.quad(integrand, theta, theta_head, args=(b,))[0]  # 从0到theta的弧长



import numpy as np
from scipy.optimize import fsolve


def find_intersection(b):
    def equations(vars):
        theta = vars[0]  # 只需要 theta 来求解交点
        x = b * theta * np.cos(theta)  # 螺线的 x 坐标
        y = b * theta * np.sin(theta)  # 螺线的 y 坐标

        # 圆的方程
        circle_eq = x ** 2 + y ** 2 - 4.5 ** 2

        return [circle_eq]

    # 初始猜测，假设 theta 在 0 到 2π 之间
    initial_guess = [1.0]  # 可以根据需要调整
    theta_solution = fsolve(equations, initial_guess)

    # 计算交点的坐标
    intersection_x = b * theta_solution[0] * np.cos(theta_solution[0])
    intersection_y = b * theta_solution[0] * np.sin(theta_solution[0])

    return theta_solution[0],intersection_x, intersection_y





def calculate_distance_to_line(points_x, points_y, point_p_x, point_p_y):
    """计算点(px, py)到由点(x1, y1)和(x2, y2)构成的直线的距离"""
    min_distance = float('inf')  # 使用无穷大初始化最小距离

    for i in range(len(points_x) - 1):
        # 直线的系数 A, B, C
        A = points_y[i + 1] - points_y[i]
        B = points_x[i] - points_x[i + 1]
        C = (points_x[i + 1] * points_y[i]) - (points_x[i] * points_y[i + 1])

        # 计算垂足的 x 坐标
        if A != 0:  # 确保不除以零
            foot_x = (B * (B * point_p_x - A * point_p_y) - A * C) / (A ** 2 + B ** 2)
            foot_y = (A * (A * point_p_y - B * point_p_x) - B * C) / (A ** 2 + B ** 2)
        else:
            # 如果 A 为 0，表示水平线，垂足的 y 值等于直线的 y 值
            foot_x = point_p_x
            foot_y = points_y[i]  # 选择线段的 y 值

        # 判断垂足是否在线段之间
        if (min(points_x[i], points_x[i + 1]) <= foot_x <= max(points_x[i], points_x[i + 1]) and
                min(points_y[i], points_y[i + 1]) <= foot_y <= max(points_y[i], points_y[i + 1])):
            # 计算到直线的距离
            distance = abs(A * point_p_x + B * point_p_y + C) / np.sqrt(A ** 2 + B ** 2)

            # 更新最小距离
            if distance < min_distance:
                min_distance = distance

    return min_distance if min_distance != float('inf') else None  # 返回最小距离或 None








# 计算初始龙身位置
def dragon_position(start_time = 0,end_time = 1,new_b = 0.55 / (2 * np.pi), if_p_point = False, interval= 1.0):
    b= new_b
    crash_head_theta= 1000

    for item in np.arange(start_time, end_time, interval):
        plt.clf()
        points_x.clear()
        points_y.clear()
        target_length = item
        # 寻找使得 arc_length(theta) = target_length 的 θ
        theta_solution = spo.brentq(lambda x: arc_length(x,b) - target_length, theta_end, theta_head)

        # 计算对应的 r 和 (x, y) 坐标
        r_solution = b * theta_solution
        x_solution = r_solution * np.cos(theta_solution)
        y_solution = r_solution * np.sin(theta_solution)

        #计算出龙头位置
        points_x.append(x_solution)
        points_y.append(y_solution)
        head_theta = theta_solution
        for length in bench_lengths:
            # 将圆的方程转化为一个函数
            def circle_equation(theta_val):
                r_ed = b * theta_val
                x_eq = r_ed * np.cos(theta_val)
                y_eq = r_ed * np.sin(theta_val)
                # print((x_eq - x_solution)**2 + (y_eq - y_solution)**2 - length**2)
                return np.abs((x_eq - x_solution)**2 + (y_eq - y_solution)**2 - length**2)

            # 寻找解，使用 fsolve
            theta_guess = theta_solution+0.3*np.pi  # 使用当前解作为初始猜测
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

        if if_p_point:
            # 使用 fsolve 找到一个解
            def equations(vars):
                x_val, y_val = vars
                return [
                    (x_val - points_x[0]) ** 2 + (y_val - points_y[0]) ** 2 - m,
                    (x_val - points_x[1]) ** 2 + (y_val - points_y[1]) ** 2 - n
                ]

            initial_guess = [points_x[0]*1.5, points_y[0]*1.5]  # 初始猜测
            solution = spo.fsolve(equations, initial_guess)
            p_x = solution[0]
            p_y = solution[1]
            point_p_x.append(p_x)
            point_p_y.append(p_y)

            min_dintance = calculate_distance_to_line(points_x, points_y, p_x , p_y)
            if min_dintance <= 0.15:
                crash_head_theta = head_theta
                print(item)
                print(b*2*np.pi)

                break


    plt.plot(b * theta * np.cos(theta), b * theta * np.sin(theta), label='Archimedean Spiral')
    return points_x, points_y, point_p_x, point_p_y ,crash_head_theta




