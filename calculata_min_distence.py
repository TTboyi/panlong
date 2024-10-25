

def calculate_distance_to_line(points_x, points_y,  point_p_x, point_p_y):
    """计算点(px, py)到由点(x1, y1)和(x2, y2)构成的直线的距离"""
    min_distance = 1000
    for i in range(len(points_x)):
        A = points_y[i + 1] - points_y[i]
        B = points_x[i] - points_x[i + 1]
        C = points_x[i + 1] * points_y[i] - points_x[i] * points_y[i + 1]

        distance = abs(A * point_p_x + B * point_p_y + C) / np.sqrt(A ** 2 + B ** 2)

        if distance < min_distance:
            min_distance = distance

    return min_distance

