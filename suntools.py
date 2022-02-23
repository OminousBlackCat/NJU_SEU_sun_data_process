import numpy as np
import re
import os
import time
import matplotlib
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from PIL import Image
import scipy.signal as signal

# 谱线矫正
# 参数data: 图像数据(numpy标准格式, 二维数组)
# 参数x0: 曲线矫正对称轴
# 参数C: 曲线矫正二次系数
# 输出：矫正后的图像数据
def curve_correction(imgData, x0, C):
    H, W = imgData.shape
    # print(H,W)
    for x in range(W):
        stdx = np.arange(0, H, 1)
        stdx = stdx / (C * (x - x0) * (x - x0) + 1)
        stdy = imgData[:, x]
        now = 1
        for y in range(H):
            while now < H - 1 and stdx[now] < y:
                now += 1
            # data[y][x] = max(stdy[now],stdy[now-1])
            if y > stdx[now]:
                imgData[y][x] = stdx[now]
            else:
                if stdx[now] - stdx[now - 1] < 2:
                    imgData[y][x] = stdy[now - 1]
                else:
                    imgData[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                            y - stdx[now - 1])
    return imgData


# 平场计算
# 参数flatData: 平场图像(numpy标准格式, 二维数组)
# 输出: 平场数据
def getFlat(flatData):
    H, W = flatData.shape
    mean = np.sum(flatData, axis=1) / W
    for i in range(H):
        flatData[i, :] = flatData[i, :] / mean[i]
    return flatData


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
    dataX = []
    dataY = []
    with open(filepath) as f:
        line = f.readline()
        while line:
            line = re.findall(r"\d+\.?\d*", line)
            if len(line) > 0:
                dataX.append(float(line[0]))
                dataY.append(float(line[1]))
            line = f.readline()
    dataX = np.array(dataX)
    dataX = (dataX - 6559.5804) / 0.024202301
    now = 1
    ansY = []
    for i in range(400):
        while(dataX[now]<i):
            now += 1
        ansY.append(dataY[now-1] + (dataY[now] - dataY[now-1]) / (dataX[now] - dataX[now-1]) * (i - dataX[now-1]))
    ansY = np.array(ansY)
    return ansY/np.max(ansY)


# 红蓝移矫正
def RB_repair(imgData, sun_std):
    H, W = imgData.shape
    sun_image = np.sum(imgData, axis=1) / W
    sun_image /= np.max(sun_image)
    cov = np.polyfit(np.arange(0, H, 1), sun_image / sun_std[0:H], 1)
    k, b = cov[0], cov[1]
    # print(k,b)
    for i in range(H):
        imgData[i, :] = imgData[i, :] / (k * i + b)
    return imgData


# 矩阵减法
def subtract(data_A, data_B):
    H, W = data_A.shape
    for i in range(H):
        for j in range(W):
            if data_A[i][j] > data_B[i][j]:
                data_A[i][j] = data_A[i][j] - data_B[i][j]
            else:
                data_A[i][j] = data_B[i][j] - data_A[i][j]
    return data_A


# 图像均值平滑操作
# 参数
def MeanSmooth(imgData, winSize=4):
    H, W = imgData.shape
    offset = int(winSize / 2)
    OffsetData = np.zeros((H + offset * 2, W + offset * 2))
    SmoothData = np.zeros((H + offset * 2, W + offset * 2))
    OffsetData[offset:offset+H,offset:offset+W] = imgData[0:H,0:W]
    for i in range(winSize):
        for j in range(winSize):
            SmoothData[offset:offset+H,offset:offset+W] += OffsetData[i : H + i,j:W + j]
    return SmoothData[offset: offset + H, offset: offset + W]/winSize/winSize


# 图像中值平滑操作
# 参数
def MedSmooth(imgData, winSize=4):
    imgData = signal.medfilt(imgData, kernel_size=3)
    # H, W = imgData.shape
    # offset = int(winSize / 2)
    # SmoothData = np.zeros((H + offset * 2, W + offset * 2))
    # for i in range(H-winSize):
    #     for j in range(W-winSize):
    #         SmoothData[offset+i][offset+j]=np.median(imgData[i:i+winSize,j:j+winSize])
    return imgData


if __name__ == "__main__":
    matplotlib.rcParams['font.sans-serif'] = ['KaiTi']
    filepath_result = "data/"
    filepath_test = "testData/"
    filepath_bash = "bass2000.txt"
    base = get_Sunstd(filepath_bash)
    # print(base)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file), dtype=float)

    image_file = get_pkg_data_filename(filepath_test + 'for_flat.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)

    filelist = os.listdir(filepath_test)

    # flat_data = np.zeros(dark_data.shape)
    # for i in range (200):
    #     image_file = get_pkg_data_filename(filepath_test + filelist[2712+i])
    #     flat_data += np.array(fits.getdata(image_file), dtype=float)
    # data = curve_correction(flat_data/400 - dark_data, 2321.26, 1.92909e-011)
    # data = smooth(data)
    data = curve_correction(flat_data , 2321.26, 1.92909e-011)
    data = getFlat(data)
    # plt.figure()
    # plt.imshow(data, cmap="gray")
    # plt.show()
    # plt.figure()
    # plt.imshow(data, cmap="gray")
    # plt.show()
    print("Ping is over")
    plt.figure()
    time_start = time.time()
    image_file = get_pkg_data_filename(filepath_test + filelist[1700])
    test_data = np.array(fits.getdata(image_file), dtype=float)
    plt.subplot(5, 1, 1)
    plt.imshow(test_data - dark_data, cmap="gray")
    plt.title('去暗场')
    test_data = curve_correction(test_data - dark_data, 2321.26, 1.92909e-011)
    plt.subplot(5, 1, 2)
    plt.imshow(test_data, cmap="gray")
    plt.title('谱线矫正')
    test_data = test_data / data
    plt.subplot(5, 1, 3)
    plt.imshow(test_data, cmap="gray")
    plt.title('去平场')
    test_data = RB_repair(test_data, base)
    plt.subplot(5, 1, 4)
    plt.imshow(test_data, cmap="gray")
    plt.title('红蓝翼矫正')
    time_end1 = time.time()
    test_data = MedSmooth(test_data)
    plt.subplot(5, 1, 5)
    plt.imshow(test_data, cmap="gray")
    plt.title('滤波')
    time_end = time.time()
    print(time_end - time_start)
    print(time_end - time_end1)
    plt.figure()
    plt.imshow(test_data, cmap="gray")
    plt.show()
