import numpy as np


# 谱线矫正
# 参数data: 图像数据(numpy标准格式, 二维数组)
# 参数x0: 曲线矫正对称轴
# 参数C: 曲线矫正二次系数
# 输出：矫正后的图像数据
def Curve_correction(data, x0, C):
    H, W = data.shape
    # print(H,W)
    for x in range(W):
        stdx = np.arange(0, H, 1)
        stdx = stdx / (C * (x - x0) * (x - x0) + 1)
        stdy = data[:, x]
        now = 1
        for y in range(H):
            while now < H - 1 and stdx[now] < y:
                now += 1
            # data[y][x] = max(stdy[now],stdy[now-1])
            if y > stdx[now]:
                data[y][x] = stdx[now]
            else:
                if stdx[now] - stdx[now - 1] < 2:
                    data[y][x] = stdy[now - 1]
                else:
                    data[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                            y - stdx[now - 1])
    return data


# 平场计算
# 参数data: 图像数据(numpy标准格式, 二维数组)
# 输出: 平场数据
def Pingchang(data):
    H, W = data.shape
    mean = np.sum(data, axis=1)
    for i in range(H):
        data[i, :] = data[i, :] / mean[i]
    return data


# 红蓝移补偿参考线拟合
# 参数x,y: numpy数组
# 输出: 直线参数k,b
def linefit(x, y):  # np.polyfit(X, Y, 1)
    N = len(x)
    sx, sy, sxx, syy, sxy = 0, 0, 0, 0, 0
    for i in range(N):
        sx += x[i]
        sy += y[i]
        sxx += x[i] * x[i]
        syy += y[i] * y[i]
        sxy += x[i] * y[i]
    k = (sy * sx / N - sxy) / (sx * sx / N - sxx)
    b = (sy - k * sx) / N
    return k, b


# 红蓝移矫正
def RB_repair(data):
    H, W = data.shape
    k, b = linefit(sum(data, axis=1), data)
    for i in range(H):
        data[i, :] = data[i, :] / (k * i + b)
    return data
