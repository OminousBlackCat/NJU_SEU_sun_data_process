import numpy as np
import re
import os
import time
from astropy.io import fits
import suntools
import matplotlib.pyplot as plt
from PIL import Image
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from numpy import polyfit, poly1d
import sys


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
def GetFlat(data):
    H, W = data.shape
    mean = np.sum(data, axis=1)/W
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
    k = (float)(sy * sx / N - sxy) / (sx * sx / N - sxx)
    b = (float)(sy - k * sx) / N
    return k, b


# 标准太阳光谱获取
# 参数filepath：储存标准光谱的文件路径
# 输出：文件包含数据
def get_Sunstd(filepath):
    print("标准光谱数据地址：" + filepath)
    data = []
    with open(filepath) as f:
        line = f.readline()
        if len(re.findall(r"\d+\.?\d*", line)) > 0:
            data.append(float(re.findall(r"\d+\.?\d*", line)[1]))
        while line:
            line = f.readline()
            if len(re.findall(r"\d+\.?\d*", line)) > 0:
                data.append(float(re.findall(r"\d+\.?\d*", line)[1]))
    return np.array(data)


# 红蓝移矫正
def RB_repair(data, sun_std):
    H, W = data.shape
    cov = np.polyfit( np.arange(0, H, 1),np. sum(data, axis=1) / W / sun_std[0:H] , 1)
    k, b = cov[0],cov[1]
    #print(k,b)
    for i in range(H):
        data[i, :] = data[i, :] / (k * i + b)
    return data

# 矩阵减法
def subtract(data_A,data_B):
    H, W = data_A.shape
    for i in range(H):
        for j in range(W):
            if data_A[i][j]>data_B[i][j]:
                data_A[i][j]=data_A[i][j]-data_B[i][j]
            else:
                data_A[i][j]=data_B[i][j]-data_A[i][j]
    return data_A

#平滑操作
def Smooth(data):
    H, W = data.shape
    SmoothData = data
    win = 3
    for x in range(H-win):
        for y in range(W-win):
            SmoothData[x][y] = np.median(data[x:x+win,y:y+win].reshape(-1))
    return SmoothData

if __name__ == "__main__":

    filepath_result = "testResult\\"
    filepath_test = "testData\\"
    filepath_bash="bass2000.txt"
    base = np.array(get_Sunstd(filepath_bash),dtype=float)
    filelist = os.listdir(filepath_test)
    #print(filelist)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file),dtype=float)

    image_file = get_pkg_data_filename(filepath_test + 'for_flat.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)

    data = Curve_correction(flat_data, 2321.26, 1.92909e-011)
    data=GetFlat(data)
    # plt.figure()
    # plt.imshow(data, cmap="gray")
    # plt.show()
    print("Ping is over")
    time_start = time.time()
    image_file = get_pkg_data_filename(filepath_test + filelist[2022])
    test_data = np.array(fits.getdata(image_file), dtype=float)
    test_data = Curve_correction(test_data-dark_data, 2321.26, 1.92909e-011)
    test_data = RB_repair(test_data,base)
    test_data = test_data/data
    time_end1 = time.time()
    test_data = Smooth(test_data)
    time_end = time.time()
    print(time_end - time_start)
    print(time_end - time_end1)
    plt.figure()
    plt.imshow(test_data, cmap="gray")
    plt.show()

