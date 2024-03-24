"""
本py文件内部包含了所有在对图像预处理时使用的处理函数与工具函数
所有与算法相关的代码都在本文件内

@author: seu_lcl
@editor: seu_wxy
"""
import copy
# import jplephem
import datetime
import os.path
import random
import re
import time
from math import *

import sys
import astropy
import cv2
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
from PIL import Image
from astropy import coordinates
from astropy.io import fits
from astropy.time import Time
from astropy.utils.data import get_pkg_data_filename
from numba import jit
from scipy import interpolate

import config
import sim


# 调试用
np.set_printoptions(threshold=sys.maxsize)

cv2.setNumThreads(1)

# 定义参数
bin_count = config.bin_count
if bin_count == 1:
    height_Ha = config.height_Ha  # ha窗口的长度
    height_Fe = config.height_Fe  # fe窗口的长度
    HA = config.HA_start  # 红蓝移HA参数
    FE = config.FE_start  # 红蓝移FE参数
    K = config.wavelength_resolution_bin_1  # 红蓝移K参数
    # K = K * bin_count
    x0 = config.curve_cor_x0_bin_1
    C = config.curve_cor_C_bin_1
else:
    height_Ha = config.height_Ha  # ha窗口的长度
    height_Fe = config.height_Fe  # fe窗口的长度
    HA = config.HA_start  # 红蓝移HA参数
    FE = config.FE_start  # 红蓝移FE参数
    K = config.wavelength_resolution_bin_2  # 红蓝移K参数
    # K = K * bin_count
    x0 = config.curve_cor_x0_bin_2
    C = config.curve_cor_C_bin_2

Ha_Lower = int(config.Ha_lower / bin_count)
Ha_Upper = int(config.Ha_Upper / bin_count)
Fe_Lower = int(config.Fe_lower / bin_count)
Fe_Upper = int(config.Fe_Upper / bin_count)

HA_lineCore = config.HA_lineCore
FE_lineCore = config.FE_lineCore

txt_size = config.date_font_size  # 字体大小
txt_thick = config.date_font_thick  # 字体粗细
Interpolation_parameter = config.interpolation_parameter  # 插值算法次数
# 以第一个大于目标坐标为起点
Interpolation_front = int((Interpolation_parameter + 1) / 2)  # 插值起点
Interpolation_back = int(Interpolation_parameter / 2) + 1  # 插值终点


# 修改平场
def FlatNormalization(flatData):
    H, W = flatData.shape
    flatData[:, 0:int(20 / bin_count)] = np.ones(
        (H, int(20 / bin_count)))  # np.clip(flatData[:,0:int(20/bin_count)],1,200000)
    flatData[:, -int(68 / bin_count):] = np.ones(
        (H, int(68 / bin_count)))  # np.clip(flatData[:, -int(68 / bin_count):],1,200000)
    return flatData * 1.4


# 多次插值
@jit(nopython=True)
def Interpolation(X_data: np.array, Y_data: np.array, x: float) -> float:
    N = len(X_data)
    if N == 1:
        return Y_data[0]
    ans = np.ones(N)
    output = 0
    for i in range(N):
        for j in range(N):
            if i != j:
                ans[i] = ans[i] * (x - X_data[j]) / (X_data[i] - X_data[j])
        output = output + ans[i] * Y_data[i]
    return output


# 谱线矫正
# 参数data: 图像数据(numpy标准格式, 二维数组)
# 参数x0: 曲线矫正对称轴
# 参数C: 曲线矫正二次系数
# 参数bin: 模式参数
# 输出：矫正后的图像数据
@jit(nopython=True)
def curve_correction(imgData, x0, C):
    # 获取图片高度和宽度
    H, W = imgData.shape
    x0 = x0 - 1
    # 定义两个窗口的高度 H窗口height_Ha行 Fe窗口height_Fe行
    bad_Ha = int(height_Ha / bin_count) - int(24 / bin_count)
    bad_Fe = int(height_Fe / bin_count) - int(24 / bin_count)
    ansData = np.zeros((bad_Fe + bad_Ha, W))
    # 进行矫正操作
    # 分两个窗口分别操作
    for x in range(W):
        # 对于H窗口进行操作
        # 计算原坐标经过变换之后对应的坐标
        stdx = np.arange(0, int(height_Ha / bin_count), 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * K + HA) / (C * (x - x0) * (x - x0) + 1) - HA) / K
        # 获取原数据值
        stdy = imgData[:, x]
        # 确定插值的坐标
        now = 1
        for y in range(bad_Ha):
            # 移动到第一个大于该坐标的地方
            while now < int(height_Ha / bin_count) - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now]:
                ansData[y][x] = stdy[now]
            else:
                # 计算插值
                ansData[y][x] = Interpolation(
                    stdx[max(0, now - Interpolation_front):min(now + Interpolation_back, int(height_Ha / bin_count))],
                    stdy[max(0, now - Interpolation_front):min(now + Interpolation_back, int(height_Ha / bin_count))],
                    y)
                # 计算插值
                # if now > 1:
                #     imgData[y][x] = QuadraticInterpolation(stdx[now - 2], stdy[now - 2] ,stdx[now - 1], stdy[now - 1] ,
                #                                            stdx[now], stdy[now] , y)
                # else:
                #     imgData[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                #         y - stdx[now - 1])

        # 对于Fe窗口进行操作
        # 计算原坐标经过变换之后对应的坐标
        stdx = np.arange(0, int(height_Fe / bin_count), 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * K + FE) / (C * (x - x0) * (x - x0) + 1) - FE) / K
        # 获取原数据值
        stdy = imgData[int(height_Ha / bin_count):int(height_Fe / bin_count) + int(height_Ha / bin_count), x]
        # 确定插值的坐标
        now = 1
        for y in range(bad_Fe):
            # 移动到第一个大于该坐标的地方
            while now < int(height_Fe / bin_count) - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now]:
                ansData[y + bad_Ha][x] = stdy[now]
            else:
                # 计算插值
                ansData[y + bad_Ha][x] = Interpolation(
                    stdx[max(0, now - Interpolation_front):min(now + Interpolation_back, int(height_Ha / bin_count))],
                    stdy[max(0, now - Interpolation_front):min(now + Interpolation_back, int(height_Ha / bin_count))],
                    y)
                # 计算插值
                # if now > 1:
                #     imgData[y][x] = QuadraticInterpolation(stdx[now - 2], stdy[now - 2], stdx[now - 1], stdy[now - 1],
                #                                            stdx[now], stdy[now], y)
                # else:
                #     imgData[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                #             y - stdx[now - 1])

    # print(bad_Ha,bad_Fe)
    # 删除坏行 并输出两窗口最后的行数

    return ansData, bad_Ha, bad_Fe


# 平场计算
# 参数flatData: 平场图像(numpy标准格式, 二维数组)
# 输出: 平场数据
def getFlat(flatData):
    # 获取图片高度和宽度
    H, W = flatData.shape
    # 计算平均轮廓
    mean = np.sum(flatData, axis=1) / W
    # 计算平场
    for i in range(H):
        flatData[i, :] = flatData[i, :] / mean[i]
    return flatData


# 标准太阳光谱获取
# 参数filepath：储存标准光谱的文件路径
# 输出：文件包含数据
def get_Sunstd(filepath):
    dataX = []
    dataY = []
    # 获取太阳标准数据
    with open(filepath) as f:
        line = f.readline()
        while line:
            line = re.findall(r"\d+\.?\d*", line)
            if len(line) > 0:
                dataX.append(float(line[0]))
                dataY.append(float(line[1]))
            line = f.readline()

    # 数据类型转换
    dataX = np.array(dataX)
    now = 1
    ansY = []
    stdx = np.zeros(int(height_Fe / bin_count + height_Ha / bin_count))
    # 根据窗口转换成需要的数据
    # 坐标变换
    stdx[0:int(height_Ha / bin_count)] = np.arange(0, int(height_Ha / bin_count), 1) * K + HA
    stdx[int(height_Ha / bin_count):] = np.arange(0, int(height_Fe / bin_count), 1) * K + FE
    # 插值
    for i in range(int(height_Fe / bin_count + height_Ha / bin_count)):
        # 找到插值所需的两侧点
        while dataX[now] < stdx[i] and now < len(dataX) - 1:
            now += 1
            # 计算插值
        ansY.append(Interpolation(
            dataX[max(0, now - Interpolation_front):min(now + Interpolation_back, len(dataX))],
            dataY[max(0, now - Interpolation_front):min(now + Interpolation_back, len(dataX))],
            stdx[i]))
        # 进行插值操作
        # if now > 1:
        #     ansY.append(QuadraticInterpolation(dataX[now - 2], dataY[now - 2], dataX[now - 1], dataY[now - 1],
        #                                            dataX[now], dataY[now], stdx[i]))
        # else:
        #     ansY.append(
        #         dataY[now - 1] + (dataY[now] - dataY[now - 1]) / (dataX[now] - dataX[now - 1]) * (stdx[i] - dataX[now - 1]))
    # 类型转换
    ansY = np.array(ansY)
    # 归一化输出
    return ansY / np.max(ansY)


# 吸收系数获取
# 参数filepathFE/HA：储存标准光谱的文件路径
# 输出：文件包含数据
def get_Absorstd(filepathHA, filepathFE, HofHa, HofFe):
    ansY = np.zeros(height_Ha + height_Fe)
    # 获取太阳标准数据
    i = 0
    with open(filepathHA) as f:
        line = f.readline()
        while line:
            line = re.findall(r"\d+\.?\d*", line)
            if len(line) > 0:
                ansY[i] = float(line[0])
                i += 1
            line = f.readline()
    f.close()
    i = height_Ha

    with open(filepathFE) as f:
        line = f.readline()
        while line:
            line = re.findall(r"\d+\.?\d*", line)
            if len(line) > 0:
                ansY[i] = float(line[0])
                i += 1
            line = f.readline()
    f.close()

    if bin_count == 2:
        for i in range(int((height_Ha + height_Fe) / 2)):
            ansY[i] = (ansY[i * 2] + ansY[i * 2 + 1]) / 2
        ansY[HofHa:HofHa + HofFe] = ansY[int(height_Ha / 2):int(height_Ha / 2) + HofFe]
    else:
        ansY[HofHa:HofHa + HofFe] = ansY[height_Ha:height_Ha + HofFe]
    ansY = ansY[0:HofHa + HofFe]
    # 归一化输出
    return ansY / np.max(ansY)


# 红蓝移矫正
# 参数bin: 模式参数
def RB_getdata(imgData, sun_std, HofHa, HofFe):
    temp_std = np.array(sun_std)
    ans = np.zeros(HofHa + HofFe)
    # 获取图片尺寸
    H, W = imgData.shape
    # 获取波长强度和
    sun_image = np.sum(imgData, axis=1) / W
    # 归一化
    # print(np.max(sun_image))
    sun_image /= np.max(sun_image)
    # # 提取所需要的对应数据 for i in range(HofHa): ans[i] = (temp_std[i * bin_count] + temp_std[i * bin_count + bin_count -
    # 1]) / 2 for i in range(HofFe): ans[HofHa + i] = (temp_std[i * bin_count + height_Ha] + temp_std[i * bin_count +
    # height_Ha + bin_count - 1]) / 2
    stdx = np.zeros(H)
    ans[0:HofHa] = temp_std[0:HofHa]
    ans[HofHa:] = temp_std[int(height_Ha / bin_count):int(height_Ha / bin_count) + HofFe]
    # 坐标转化为波长
    stdx[0:HofHa] = np.arange(0, HofHa, 1) * K + HA
    stdx[HofHa:] = np.arange(0, HofFe, 1) * K + FE
    # 拟合一次函数
    cov = np.polyfit(stdx, sun_image / ans[0:H], 1)
    k, b = cov[0], cov[1]
    # 红蓝翼矫正
    RB_absorb = np.zeros(H)
    for i in range(H):
        RB_absorb[i] = k * stdx[i] + b
    return RB_absorb


def RB_repair(imgData, RB_absorb):
    H, W = imgData.shape
    for i in range(H):
        imgData[i, :] = imgData[i, :] / RB_absorb[i]
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
    OffsetData[offset:offset + H, offset:offset + W] = imgData[0:H, 0:W]
    for i in range(winSize):
        for j in range(winSize):
            SmoothData[offset:offset + H, offset:offset + W] += OffsetData[i: H + i, j:W + j]
    return SmoothData[offset: offset + H, offset: offset + W] / winSize / winSize


def DivFlat(imgData, flatData):
    np.seterr(divide='ignore', invalid='ignore')
    return imgData / flatData


# 图像中值平滑操作
# 参数
def MedSmooth(imgData, HofHa, HofFe, winSize=3):
    H, W = imgData.shape
    img = np.zeros([HofHa + HofFe + 4, W])
    img[1:1 + HofHa] = imgData[0:HofHa]
    img[0] = imgData[0]
    img[1 + HofHa] = imgData[HofHa]

    img[3 + HofHa:3 + HofHa + HofFe] = imgData[HofHa:]
    img[2 + HofHa] = imgData[HofHa]
    img[3 + HofHa + HofFe] = imgData[HofHa + HofFe - 1]

    if bin_count == 1:
        img = signal.medfilt(img, kernel_size=winSize)
    if bin_count == 2:
        img = signal.medfilt(img, kernel_size=winSize - 2)

    imgData[:HofHa] = img[1:1 + HofHa]
    imgData[HofHa:] = img[HofHa + 3: 3 + HofHa + HofFe]
    return imgData


def amplify(Data):
    H, W = Data.shape
    amplify_size = 10
    Data1 = np.zeros((H, (W - 1) * amplify_size + 1))
    for j in range(W - 1):
        for k in range(amplify_size):
            Data1[:, j * amplify_size + k] = Data[:, j] * (1 - k / amplify_size) + Data[:, j + 1] * (k / amplify_size)
    Data1[:, (W - 1) * amplify_size] = Data[:, W - 1]
    return Data1


# 计算偏差
def getFlatOffset(flatData, imgData):
    flatTempData = np.array(flatData)
    HofHa = int(height_Ha / bin_count)
    HofFe = int(height_Fe / bin_count)
    # 获取图片尺寸
    H, W = imgData.shape
    # 计算中心点

    cx = int(H / 2)
    cy = int(W / 2)
    # 获取序列
    # img = imgData[cx :cx + int(H / 4), cy - int(W / 8):cy + int(W / 8)]
    # flat = flatData[cx :cx + int(H / 4), cy - int(W / 8):cy + int(W / 8)]
    img = imgData[HofHa - int(H / 8):HofHa + int(H / 8), cy - int(W / 8):cy + int(W / 8) + 1]
    flat = flatTempData[HofHa - int(H / 8):HofHa + int(H / 8), cy - int(W / 8):cy + int(W / 8) + 1]
    img = amplify(img)
    flat = amplify(flat)
    # FFT变化
    imgFFT = np.fft.fft2(img)
    flatFFT = np.fft.fft2(flat)

    # plt.figure()
    # plt.plot(img[50, :].reshape(-1))
    # plt.plot(flat[50, :].reshape(-1))
    # plt.show()
    # print(np.unravel_index(np.argmax(np.abs(img[10, :650])), img[10, :650].shape))
    # print(np.unravel_index(np.argmax(np.abs(flat[10, :650])), flat[10, :650].shape))
    # IFFT
    FR = imgFFT * np.conj(flatFFT)
    R = np.fft.ifft2(FR)
    R = np.fft.fftshift(R)
    h, w = img.shape
    # plt.figure()
    # plt.imshow(np.abs(R), cmap="gray", aspect='auto')
    # plt.show()

    # 获取最大值坐标
    pos = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    # print(pos)
    # 计算偏移量
    mx = int((pos[1] - int((w + 1) / 2)) / 10)
    my = (pos[1] - int((w + 1) / 2)) - mx * 10
    log("偏移量：" + str(mx) + str(my))

    # 偏移操作
    if mx * 10 + my < 0:
        my *= -1
        flatTempData[:, 0:W + mx - 1] = flatTempData[:, -mx:W - 1] * (1 - my / 10) + flatTempData[:, -mx + 1:W] * (
                my / 10)
    else:
        flatTempData[:, mx + 1:W - 1] = flatTempData[:, 0:W - mx - 2] * my / 10 + flatTempData[:, 1:W - mx - 1] * (
                1 - my / 10)
    return flatTempData


# 对图像进行横向上的平移
def moveImg(imgdata, offset):
    H, W = imgdata.shape
    if offset < 0:
        imgdata[int(height_Ha / bin_count) + 1:, 0:W + int(offset / bin_count)] = imgdata[
                                                                                  int(height_Ha / bin_count) + 1:,
                                                                                  -int(offset / bin_count):W]
    else:
        imgdata[int(height_Ha / bin_count) + 1:, int(offset / bin_count):W] = imgdata[int(height_Ha / bin_count) + 1:,
                                                                              0:W - int(offset / bin_count)]
    return imgdata


# 从色谱文件中获得输出png的色谱
def get_color_map(fname):
    colors = []
    with open(fname) as f:
        line = f.readline()
        while line:
            line = re.findall(r"\d+\.?\d*", line)
            if len(line) > 0:
                # 对颜色值进行归一化
                colors.append([float(line[0]) / 255, float(line[1]) / 255, float(line[2]) / 255])
            line = f.readline()

    clrmap = matplotlib.colors.LinearSegmentedColormap.from_list("mycmap", colors)
    return clrmap


# 将bin=1的图转换为bin=2（测试时使用）
def change(img):
    if bin_count == 1:
        return img
    H, W = img.shape
    ans = np.zeros([int(H / 2), int(W / 2)])
    for i in range(int(H / 2)):
        for j in range(int(W / 2)):
            ans[i][j] = (img[i * 2][j * 2] + img[i * 2 + 1][j * 2] + img[i * 2][j * 2 + 1] + img[i * 2 + 1][
                j * 2 + 1]) / 4
    return ans


# 获取图片所在的bin值
def getBin(imgData):
    H, W = imgData.shape
    log(H, W, height_Ha + height_Fe)
    if H >= height_Ha + height_Fe:
        return 1
    return 2


# 完整的工作流程
def entireWork(filename, darkDate, flatData, abortion):
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    # imgData = change(imgData)
    imgData = moveImg(imgData, -2)
    T1 = time.perf_counter()
    imgData, HofHa, HofFe = curve_correction(imgData - darkDate, x0, C)
    T2 = time.perf_counter()
    print('程序运行时间:%s毫秒' % ((T2 - T1) * 1000))
    print("去暗场")
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    print("去平场")
    plt.figure()
    plt.imshow(flatData, cmap="gray", aspect='auto')
    plt.show()
    plt.figure()
    plt.imshow(DivFlat(imgData, flatData), cmap="gray", aspect='auto')
    plt.show()
    # print(HofHa, HofFe)
    print(flatData.min())
    imgData = DivFlat(imgData, flatData)
    # plt.figure()
    # plt.imshow(imgData, cmap="gray", aspect='auto')
    # plt.show()
    # print(HofHa)
    plt.figure()
    plt.plot(imgData[HofHa:, 1100].reshape(-1))
    imgDataRB = RB_repair(imgData, abortion)
    plt.plot(imgDataRB[HofHa:, 1100].reshape(-1))
    imgData = MedSmooth(np.array(imgDataRB), HofHa, HofFe, 5)

    plt.plot(imgData[HofHa:, 1100].reshape(-1))
    plt.show()
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    return imgDataRB, imgData


def circle(x1, y1, x2, y2, x3, y3):
    """
    :return:  x0 and y0 is center of a circle, r is radius of a circle
    """
    a = x1 - x2
    b = y1 - y2
    c = x1 - x3
    d = y1 - y3
    a1 = ((x1 * x1 - x2 * x2) + (y1 * y1 - y2 * y2)) / 2.0
    a2 = ((x1 * x1 - x3 * x3) + (y1 * y1 - y3 * y3)) / 2.0
    theta = b * c - a * d
    if abs(theta) < 1e-7:
        return 0, 0, 0
    x0 = (b * a2 - d * a1) / theta
    y0 = (c * a1 - a * a2) / theta
    r = np.sqrt(pow((x1 - x0), 2) + pow((y1 - y0), 2))
    return x0, y0, r


# 通过灰度图拟合图中的圆
def getCircle(image, idd=0):
    # 二值化
    image_max = np.max(image)
    image = np.clip((image - image_max * 0.05) * 10, 0, 1)
    # plt.figure()
    # plt.imshow(image)
    # plt.show()

    # 通过卷积，使用Sobel算子提取边界
    conv1 = np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])
    conv2 = np.array([[1, 1, 1], [1, -8, 1], [1, 1, 1]])
    conv3 = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]]) / 16
    image = signal.convolve2d(image, conv3, "valid")
    # gradient_Y = signal.convolve2d(image, conv1, "valid")
    gradient_X = signal.convolve2d(image, conv2, "valid")
    gradient = np.abs(gradient_X)  # + np.abs(gradient_Y)
    gradient = np.clip(gradient, 0, np.max(gradient) * 0.3)
    gradient_max = np.max(gradient)
    gradient = np.clip((gradient - gradient * 0.05) * 10, 0, 1)

    # plt.figure()
    # plt.imshow(gradient)
    # plt.show()
    # H, W = gradient.shape
    # gradient_reverse = np.array(gradient)
    # for i in range(W):
    #     gradient_reverse[:, i] = gradient[:, W - 1 - i]
    #
    # # plt.figure()
    # # plt.imshow(gradient_reverse)
    # # plt.show()
    #
    # # FFT变化
    # imgFFT = np.fft.fft2(gradient)
    # flatFFT = np.fft.fft2(gradient_reverse)
    #
    # # IFFT
    # FR = imgFFT * np.conj(flatFFT)
    # R = np.fft.ifft2(FR)
    # R = np.fft.fftshift(R)
    # # 获取最大值坐标
    # pos = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    # # print(pos)
    # # print((pos[1] + W / 2) / 2)  # 460 340
    # Offset = int(pos[1] - W / 2)
    # if Offset > 0:
    #     gradient_reverse[:, Offset:] = gradient_reverse[:, 0:-Offset]
    # else:
    #     Offset *= -1
    #     gradient_reverse[:, 0:-Offset] = gradient_reverse[:, Offset:]
    # gradient *= gradient_reverse

    # plt.figure()
    # plt.imshow(gradient)
    # plt.show()

    gradient = signal.medfilt(gradient, kernel_size=7)

    H, W = gradient.shape

    points = []
    for i in range(H):
        for j in range(W):
            if gradient[i][j] > 0.98:
                points.append([i, j])
    L = len(points)

    gradient *= 255  # 变换为0-255的灰度值
    im = Image.fromarray(gradient)
    im = im.convert('L')
    im = np.array(im)
    circles = cv2.HoughCircles(im, cv2.HOUGH_GRADIENT, 1, 500, param1=100, param2=10,
                               minRadius=int(800 * 2 / bin_count),
                               maxRadius=int(1000 * 2 / bin_count))
    img = np.array(im)
    if circles is None:
        raise ValueError
    if len(circles) == 0:
        raise ValueError
    # id = -1
    # goal = -1
    # now = 0

    # print(circles)
    # for circle in circles[0]:
    #     x = circle[0]
    #     y = circle[1]
    #     r = circle[2]
    #     goali = 0
    #     for i in range(L):
    #         # goali += abs((points[i][0]-x)**2 + (points[i][1]-y)**2 - r**2)
    #         if abs(math.sqrt((points[i][0]-x)**2 + (points[i][1]-y)**2) - r) < 1 and im[points[i][0]][points[i][1]]>254:
    #             goali += 1
    #     if goali > goal or id == -1:
    #         id = now
    #         goal = goali
    #     now += 1
    #     cv2.circle(img, (int(circles[0][now-1][0]), int(circles[0][now-1][1])), int(circles[0][now-1][2]), (255), 10)
    # print(now-1,goali)
    # print(id)
    id = 0
    cv2.circle(gradient, (int(circles[0][id][0]), int(circles[0][id][1])), int(circles[0][id][2]), (100), 3)
    # plt.figure()
    # plt.imshow(img)
    # plt.show()
    # plt.figure()
    # plt.imshow(im)
    # plt.show()
    #
    # gradient = signal.convolve2d(gradient, conv1, "valid")
    # gradient = np.clip(gradient - 8, 0, 1)
    # print(times)
    # print(x,y,r*0.52)
    return circles[0][id][0] + 2, circles[0][id][1] + 2, circles[0][id][2]


# 辅助计算软件的运算
def rotation_matrix3(xyz, theta):
    # theta in radian
    if xyz == 'x':
        R = np.array([[1, 0, 0], [0, cos(theta), sin(theta)], [0, -sin(theta), cos(theta)]])
    if xyz == 'y':
        R = np.array([[cos(theta), 0, -sin(theta)], [0, 1, 0], [sin(theta), 0, cos(theta)]])
    if xyz == 'z':
        R = np.array([[cos(theta), sin(theta), 0], [-sin(theta), cos(theta), 0], [0, 0, 1]])
    return R


# 计算需求数据，输入来自于头文件
def getB0P0(q0, q1, q2, q3, strtime):
    # astropy.utils.data.import_file_to_cache(config.de_file_url,
    #                                         config.de_file_name,
    #                                         remove_original=False, pkgname='astropy', replace=True)
    astropy.coordinates.solar_system_ephemeris.set(config.de_file_url)
    t = Time(strtime)

    sun_6 = astropy.coordinates.get_body_barycentric_posvel('sun', t)
    earth_6 = astropy.coordinates.get_body_barycentric_posvel('earth', t)

    sun_pos = sun_6[0]
    earth_pos = earth_6[0]

    # position is in ICRF (aberration not corrected)
    earth2sun_rx = sun_pos.x.value - earth_pos.x.value
    earth2sun_ry = sun_pos.y.value - earth_pos.y.value
    earth2sun_rz = sun_pos.z.value - earth_pos.z.value
    earth2sun_pos = np.array([[earth2sun_rx], [earth2sun_ry], [earth2sun_rz]])

    normalize_factor = sqrt(earth2sun_pos[0] ** 2 + earth2sun_pos[1] ** 2 + earth2sun_pos[2] ** 2)
    earth2sun_pos_normalize = earth2sun_pos / normalize_factor  # the satellite position to the sun position (point from satellite to sun)
    g1 = 2 * (q1 * q3 - q0 * q2)
    g2 = 2 * (q2 * q3 + q0 * q1)
    g3 = q0 ** 2 - q1 ** 2 - q2 ** 2 + q3 ** 2
    g4 = 2 * (q1 * q2 + q0 * q3)
    g5 = q0 ** 2 + q1 ** 2 - q2 ** 2 - q3 ** 2

    rad2deg = 180 / pi

    # theta, psi, gamma in degrees
    theta = atan(g2 / g3) * rad2deg
    if g3 < 0:
        theta = theta + 180

    psi = atan(g4 / g5) * rad2deg
    if g5 < 0:
        psi = psi + 180

    cosgamma = g2 / sin(theta / rad2deg)
    if abs(sin(theta / rad2deg)) < 0.5:
        cosgamma = g3 / cos(theta / rad2deg)
    singamma = -g1
    gamma = atan(singamma / cosgamma) * rad2deg
    if cosgamma < 0:
        gamma = gamma + 180

    # print('ephemeris:\n', earth2sun_pos_normalize)

    deg2rad = pi / 180
    Chase_arr_b = np.array([[0], [1], [0]])  # Chase y-axis pointing in Chase coordinates
    R1 = np.dot(rotation_matrix3('x', -theta * deg2rad), Chase_arr_b)
    R2 = np.dot(rotation_matrix3('y', -gamma * deg2rad), R1)
    R3 = np.dot(rotation_matrix3('z', -psi * deg2rad), R2)
    Chase_arr_E = R3  # Chase pointing in Equatorial coordinates
    # print('quaternion:\n', Chase_arr_E)

    arr_mul = np.multiply(earth2sun_pos_normalize, Chase_arr_E)
    # print('bias:\n', acos(arr_mul[0] + arr_mul[1] + arr_mul[2]) * 180 / pi * 60 * 60, 'arcsec')

    # the direction of the solar rotation axis in J2000
    alpha0 = 286.13
    delta0 = 63.87
    R_a0d0 = np.array([[cos(alpha0 * deg2rad) * cos(delta0 * deg2rad)], [sin(alpha0 * deg2rad) * cos(delta0 * deg2rad)],
                       [sin(delta0 * deg2rad)]])

    vn = np.multiply(earth2sun_pos_normalize, R_a0d0)
    sigvn = float(vn[0]) + float(vn[1]) + float(vn[2])
    B0 = 90 - acos(sigvn) / deg2rad
    B0 = -B0

    R_E0 = R_a0d0
    R1 = np.dot(rotation_matrix3('z', psi / rad2deg), R_E0)
    R2 = np.dot(rotation_matrix3('y', gamma / rad2deg), R1)
    R3 = np.dot(rotation_matrix3('x', theta / rad2deg), R2)
    R_b0 = R3

    # solar rotation axis projection to the sky plane (xz plane)
    R_b0xz = np.array([[float(R_b0[0])], [float(R_b0[2])]])

    # picture top vector
    # R_top = np.array([[0], [1]])

    # the angle between R_top and R_b0xz, different with p0 by a constant
    p0 = atan(R_b0xz[0] / R_b0xz[1]) * 180 / pi
    if R_b0xz[1] < 0:
        p0 = 180 + p0

    # instrument rotation clockwise
    INST_ROT = -p0
    return B0, INST_ROT


# 下采样 根据不同的bin_count 下采样不同的scale
# 对于bin = 1 下采样1/4(16个像素->1个像素)
# 对于bin = 2 下采样1/2(4个像素->1个像素)
# 顺便将输出的数组上下翻转 以适配习惯
def down_sample(data: np.array):
    scale_factor = int(4 / bin_count)
    return_array = np.zeros((int(data.shape[0] / scale_factor), int(data.shape[1] / scale_factor)))
    for i in range(return_array.shape[0]):
        for j in range(return_array.shape[1]):
            scaled_i = i * scale_factor
            scaled_j = j * scale_factor
            if scaled_i < data.shape[0] and scaled_j < data.shape[1]:
                return_array[return_array.shape[0] - i - 1][j] = np.mean(
                    data[scaled_i: scaled_i + scale_factor, scaled_j: scaled_j + scale_factor])
            else:
                return_array[return_array.shape[0] - i - 1][j] = 0
    return return_array

# 灰度图添加文字
# I_array为目标图 time_txt为目标文字
# 字体和粗细可通过config调整
def add_time(Input_array, time_txt, max_value):
    I_array = np.array(Input_array)
    H, W = I_array.shape
    h1 = int(H * 0.55)
    w1 = int(W * 0.97)
    cv2.putText(I_array, time_txt, (h1, w1), cv2.FONT_HERSHEY_PLAIN, txt_size, (max_value, max_value, max_value),
                txt_thick)
    return I_array


def cal_center_mean(Input_np):
    data = Input_np.mean(axis=1).mean(axis=1)
    # print(data.shape)
    height_ha = int(height_Ha / bin_count) - int(24 / bin_count)
    height_fe = int(height_Fe / bin_count) - int(24 / bin_count)
    log("HA 与FE 窗口大小为:", height_ha, height_fe)
    Ha_B = (data[Ha_Lower] + data[Ha_Upper]) / 2 * 0.9
    Fe_B = (data[height_ha + Fe_Lower] + data[height_ha + Fe_Upper]) / 2 * 0.9
    log("HA与FE的窗口平均值: ", Ha_B, Fe_B)
    Ha_U = Ha_L = 0
    Fe_U = Fe_L = 0

    for i in range(Ha_Lower, Ha_Upper + 1):
        if data[i] < Ha_B:
            Ha_U = Ha_U + (data[i] - Ha_B) * i
            Ha_L = Ha_L + (data[i] - Ha_B)
    for i in range(Fe_Lower, Fe_Upper + 1):
        if data[height_ha + i] < Fe_B:
            Fe_U = Fe_U + (data[height_ha + i] - Fe_B) * i
            Fe_L = Fe_L + (data[height_ha + i] - Fe_B)

    Ha_mean = Ha_U / Ha_L
    Fe_mean = Fe_U / Fe_L + (5430 - 5096) / bin_count
    log("HA_U 与 Ha_L:", Ha_U, Ha_L)
    log("Fe_U 与 Fe_L:", Fe_U, Fe_L)

    mean_K = (FE_lineCore - HA_lineCore) / (Fe_mean - Ha_mean)
    mean_b = FE_lineCore - Fe_mean * mean_K
    log("线心像元为: ", Fe_mean, Ha_mean)
    log("结果为:", mean_K, mean_b)
    return mean_K, mean_b, mean_b + mean_K * (5430 - 5096) / bin_count


# log所用函数
def log(*args):
    print(datetime.datetime.now().strftime("[%Y-%m-%d-%H:%M:%S]"), end='')
    print(*args)


# 整合南大2023年6月份提出的需求，新增太阳图像的去抖动和旋转矫正
# 原始代码提供者:绕世豪（南大天文系）
# 整合人：褚有骋

# 四元数旋转函数
def quaternion_rot(q, p):
    """
    @param q:
    @param p:
    @return
    """
    s_q = q[0]
    lambda_q = sqrt(1 - s_q ** 2)
    if lambda_q != 0:
        v_q = [q[i] / lambda_q for i in range(1, 4)]
    else:
        v_q = [0, 0, 0]
    return (2 * lambda_q ** 2 * np.dot(np.dot(v_q, p), v_q) + (s_q ** 2 - lambda_q ** 2) * p + \
            2 * lambda_q * s_q * np.cross(v_q, p))


# 旋转矫正函数
def rotate_fits(width, center_x, center_y, se00xx_imwing, rsun, data, angle, isHa, time_series_data_array=None):
    """
    对三维矩阵内的每个光谱维度的全日面像以日心为中心旋转对应角度
    @param width
    @param center_x
    @param center_y
    @param se00xx_imwing
    @param rsun
    @param data
    @param angle
    @param isHa
    @param time_series_data_array
    """
    se00xx_center = sim.circle_center(se00xx_imwing)
    centerx, centery = se00xx_center[0], se00xx_center[1]
    se00xx_rotate1 = np.zeros([center_x, center_y])
    se00xx_rotate2 = np.zeros([center_x, center_y])

    for j in range(center_x):
        for k in range(center_y):
            original_coor = [k - centerx, 0, j - centery]
            local_rot = np.dot(sim.rotation_matrix3('y', angle / 180 * pi), \
                               original_coor)
            se00xx_rotate1[j, k] = local_rot[0] + centerx
            se00xx_rotate2[j, k] = local_rot[2] + centery

    se00xx_rotate1 = se00xx_rotate1.astype('float32')
    se00xx_rotate2 = se00xx_rotate2.astype('float32')

    for j in range(width):
        data[j, :, :] = cv2.remap(data[j, :, :], se00xx_rotate1, se00xx_rotate2, \
                                  borderMode=cv2.BORDER_CONSTANT, \
                                  interpolation=cv2.INTER_LINEAR)

    # 对时间矩阵做相同的操作
    if time_series_data_array is not None:
        time_series_data_array[:, :] = cv2.remap(time_series_data_array, se00xx_rotate1, se00xx_rotate2, \
                                                 borderMode=cv2.BORDER_CONSTANT, \
                                                 borderValue=-1, \
                                                 interpolation=cv2.INTER_LINEAR)

    se00xx_mask = np.ones([center_x, center_y])
    for j in range(center_x):
        for k in range(center_y):
            if (j - centery) ** 2 + (k - centerx) ** 2 >= (1.1 * rsun) ** 2:
                se00xx_mask[j, k] = 0

    if isHa:
        rotated_data = copy.deepcopy(np.multiply(data[110 * 2 // bin_count, :, :], se00xx_mask))
    else:
        rotated_data = copy.deepcopy(np.multiply(data[10 * 2 // bin_count, :, :], se00xx_mask))

    se00xx_center = sim.circle_center(rotated_data)
    se00xx_centerx, se00xx_centery = se00xx_center[0], se00xx_center[1]
    return se00xx_centerx, se00xx_centery


def getEphemerisPos(strtime):
    """

    @param strtime
    @return
    """
    t = Time(strtime)  # 修改时间格式，以便输入星表（ephemeris）查询太阳和球的位置（J2000坐标系）
    astropy.coordinates.solar_system_ephemeris.set(config.de_file_url)
    sun_6 = astropy.coordinates.get_body_barycentric_posvel('sun', t)
    earth_6 = astropy.coordinates.get_body_barycentric_posvel('earth', t)
    return sun_6[0], earth_6[0]


def head_distortion_correction(spec_win, axis_width, biasx, biasz, sequence_data_array, time_series_data_array=None):
    """
    对头部0000序列进行畸变矫正函数, 直接修改传入的sequence_data_array np数组
    当传入的时间矩阵不为None时, 对时间矩阵做相同的畸变矫正

    @param spec_win
    @param axis_width
    @param biasx
    @param biasz
    @param sequence_data_array 传入的完整全日面太阳光谱像(引用) 函数将直接修改nparray内部数据
    @param time_series_data_array
    @return
    """
    a1, a2, a3 = axis_width
    if len(biasx) < 2313 * 2 // bin_count:
        a2 = len(biasx)
    se0000_x, se0000_z = np.zeros((a2, a3)), np.zeros((a2, a3))
    lse0000_za, lse0000_zb, lse0000_zc = [], [], []
    for j in range(a2):
        for k in range(a3):
            se0000_x[j, k] = k - biasx[j]  # 拍摄时，若镜头向右偏移，重映射时，采用相应图像向左平移
            if k == 0:
                lse0000_za.append(j)  # 完整的z方向2312帧坐标
                if j + biasz[j] not in lse0000_zb:  # 新图像j + l_se0000_biasz[j]位置由j位置替换
                    lse0000_zb.append(j + biasz[j])  # 因为不均匀，故不涵盖所有2312个z方向坐标，该位置作为插值的自变量
                    lse0000_zc.append(j)  # 该位置作为插值的因变量

    se0000_maxtestz = max(lse0000_zb)
    se0000_mintestz = min(lse0000_zb)
    se0000_tck = interpolate.splrep(lse0000_zb, lse0000_zc)
    lse0000_z = interpolate.splev(lse0000_za, se0000_tck, der=0)  # 插值
    lse0000_z = list(lse0000_z)

    se0000_shiftmask = np.ones([a2, a3])
    for j in range(a2):
        for k in range(a3):
            se0000_z[j, k] = lse0000_z[j]  # 插值之后，每个点都对应了
            if j > se0000_maxtestz or j < se0000_mintestz:
                se0000_shiftmask[j, k] = 0  # 外插的点不可信，需要擦除

    se0000_x = se0000_x.astype('float32')  # 改变数据类型，否则remap时候会报错（下同）
    se0000_z = se0000_z.astype('float32')

    for wl in range(a1):
        se0000_im_wl = sequence_data_array[wl, :, :]
        se0000_imout = cv2.remap(se0000_im_wl, se0000_x, se0000_z, borderMode=cv2.BORDER_CONSTANT, \
                                 interpolation=cv2.INTER_LINEAR)  # 非刚性位移改正畸变，边界值默认为常数0以免假信号造成误解
        se0000_imout = np.multiply(se0000_imout, se0000_shiftmask)
        sequence_data_array[wl, :, :] = se0000_imout

    se0000_imwing = None
    if spec_win == 'HA':
        se0000_imwing = sequence_data_array[110 * 2 // bin_count, :, :]
    if spec_win == 'FE':
        se0000_imwing = sequence_data_array[10 * 2 // bin_count, :, :]

    # 如果传进来了时间矩阵, 那就对时间矩阵做相同的变换
    if time_series_data_array is not None:
        time_series_data_array[:, :] = cv2.remap(time_series_data_array, se0000_x, se0000_z,
                                                 borderMode=cv2.BORDER_CONSTANT, \
                                                 borderValue=-1, \
                                                 interpolation=cv2.INTER_LINEAR)
        time_series_data_array[:, :] = np.multiply(time_series_data_array, se0000_shiftmask)

    return se0000_imwing


def non_head_distortion_correction(spec_win, se00xx_data, se0000_center, se00xx_RSUN, axis_width, se00xx_hacore,
                                   hacore0, width, height, time_series_data_array=None):
    """
    对非头部序列的序列进行畸变矫正, 需要额外输入头部序列所获日心
    当时间矩阵不为None时, 对时间矩阵做相同的畸变矫正

    @param spec_win:
    @param se00xx_data:
    @param se0000_center:
    @param se00xx_RSUN:
    @param axis_width:
    @param se00xx_hacore:
    @param hacore0:
    @param width:
    @param height:
    @param time_series_data_array:
    @return
    """
    a1, a2, a3 = axis_width
    se0000_centerX0, se0000_centerY0 = se0000_center

    if spec_win == 'HA':  # 读取光球图像用于定日心坐标
        se00xx_wing = se00xx_data[110 * 2 // bin_count, :, :]
    if spec_win == 'FE':
        se00xx_wing = se00xx_data[10 * 2 // bin_count, :, :]
    se00xx_center = sim.circle_center(se00xx_wing)
    se00xx_centerx, se00xx_centery = se00xx_center[0], se00xx_center[1]
    dx, dy = -(se00xx_centerx - se0000_centerX0), -(se00xx_centery - se0000_centerY0)
    se00xx_hacore2 = sim.imshift(se00xx_hacore, [int(dy), int(dx)])  # 刚性对齐

    se00xx_flow = cv2.calcOpticalFlowFarneback(hacore0, se00xx_hacore2, flow=None, pyr_scale=0.5, \
                                               levels=3, winsize=config.winsize * 2 // bin_count, iterations=5, poly_n=5, \
                                               poly_sigma=1.2, flags=0)  # 计算畸变

    se00xx_x1, se00xx_y1 = np.meshgrid(np.arange(width), np.arange(height))

    se00xx_x2 = copy.deepcopy(se00xx_x1)
    se00xx_y2 = copy.deepcopy(se00xx_y1)
    se00xx_flow2 = copy.deepcopy(se00xx_flow)

    for j in range(a2):  # 只选取日面范围的区域（2.1日面半径的正方形）计算每步扫描的指向偏差，或许可以自适应调节区域大小
        if abs(j - se00xx_centery) <= 1.05 * se00xx_RSUN:
            xrange = [se00xx_centerx - sqrt((1.05 * se00xx_RSUN) ** 2 - (se00xx_centery - j) ** 2), \
                      se00xx_centerx + sqrt((1.05 * se00xx_RSUN) ** 2 - (se00xx_centery - j) ** 2)]
            flowx_mean = se00xx_flow[j, int(xrange[0]):int(xrange[1]), 0].mean()
            flowy_mean = se00xx_flow[j, int(xrange[0]):int(xrange[1]), 1].mean()
            for k in range(a3):
                se00xx_flow2[j, k, 0] = flowx_mean  # 同一步的图像需要刚性平移改正畸变
                se00xx_flow2[j, k, 1] = flowy_mean

    se00xx_x2 = se00xx_x2 + se00xx_flow2[:, :, 0]
    se00xx_y2 = se00xx_y2 + se00xx_flow2[:, :, 1]
    se00xx_x2 = se00xx_x2.astype('float32')
    se00xx_y2 = se00xx_y2.astype('float32')

    se00xx_imwing = None

    for j in range(a1):
        if spec_win == 'HA':
            se00xx_im3 = se00xx_data[j, :, :]
            if j == 110 * 2 // bin_count:
                se00xx_imwing = se00xx_data[j, :, :]
        if spec_win == 'FE':
            se00xx_im3 = se00xx_data[j, :, :]
            if j == 10 * 2 // bin_count:
                se00xx_imwing = se00xx_data[j, :, :]
        se00xx_imout = cv2.remap(se00xx_im3, se00xx_x2, se00xx_y2, borderMode=cv2.BORDER_CONSTANT, \
                                 interpolation=cv2.INTER_LINEAR)  # 非刚性位移改正畸变
        se00xx_data[j, :, :] = se00xx_imout

    if time_series_data_array is not None:
        time_series_data_array[:, :] = cv2.remap(time_series_data_array, se00xx_x2, se00xx_y2, borderMode=cv2.BORDER_CONSTANT, \
                                           interpolation=cv2.INTER_LINEAR)  # 对时间矩阵也应用相同的非刚性位移改正畸变

    return se00xx_imwing


def cropPNG(data: np.array, centerX: int, centerY: int):
    """
    直接修改data数组做剪裁，取日心附近1040长度的数据
    @param data: 二维数组, Ha线心全日面矩阵
    @param centerX: X轴日心
    @param centerY: Y轴日心
    @return 返回剪裁后的数组
    """
    return_data = data[centerY - 1040:centerY + 1040, centerX - 1040:centerX + 1040]
    return return_data


def judgeBinMode(fits_list: list[str], root_path: str):
    """
    根据传入的fits_list文件绝对路径地址列表, 判断程序运行时的bin模式
    将随机读取10张fits文件, 若超过半数的fits文件分辨率均 < 2500,  则设置bin模式为2
    若超过半数的分辨率 > 2500, 则设置bin模式为1
    否则出错
    """
    big_than_2500_cnt = 0
    sml_than_2500_cnt = 0
    for i in range(10):
        current_fit_name = random.choice(fits_list)
        img_data = fits.open(os.path.join(root_path, current_fit_name))[0].data
        if img_data.shape[1] > 2500:
            big_than_2500_cnt += 1
        else:
            sml_than_2500_cnt += 1
    if big_than_2500_cnt > 5:
        return 1
    if sml_than_2500_cnt > 5:
        return 2
    return -1