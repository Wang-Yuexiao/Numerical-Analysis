import numpy as np

# 参数 (从程序输出中获取)
a1, a2, a3 = 0.981888, 0.002541, -0.375174
a4, a5, a6 = 0.001250, 0.982163, 1.157715

# 从文件读取数据
data = np.loadtxt("fitdata1.dat")
x, y, x_actual, y_actual = data[:, 0], data[:, 1], data[:, 2], data[:, 3]

# 计算预测值
x_predicted = a1 * x + a2 * y + a3
y_predicted = a4 * x + a5 * y + a6

# 计算误差
mse_x = np.mean((x_predicted - x_actual) ** 2)
mse_y = np.mean((y_predicted - y_actual) ** 2)

print(f"MSE for x': {mse_x}")
print(f"MSE for y': {mse_y}")