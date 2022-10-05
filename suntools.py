"""
本py文件内部包含了所有在对图像预处理时使用的处理函数与工具函数
所有与算法相关的代码都在本文件内

@author: seu_lcl
@editor: seu_wxy
"""

import numpy as np
import re
import os
import time
import matplotlib
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
from PIL import Image
import math
import scipy.signal as signal
import config
import astropy
# import jplephem
import datetime
import random
from math import *
from astropy.io import fits
from astropy.time import Time
from astropy import coordinates
import multiprocessing as mp
import cv2

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


# 多次插值
def Interpolation(X_data, Y_data, x):
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
def curve_correction(imgData, x0, C):
    # 获取图片高度和宽度
    H, W = imgData.shape
    x0 = x0 - 1
    # 定义两个窗口的高度 H窗口height_Ha行 Fe窗口height_Fe行
    bad_Ha = height_Ha / bin_count
    bad_Fe = height_Fe / bin_count

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
        for y in range(int(height_Ha / bin_count)):
            # 移动到第一个大于该坐标的地方
            while now < int(height_Ha / bin_count) - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now] and x > 400 / bin_count and x < W - 400 / bin_count:
                imgData[y][x] = stdy[now]
                if y < bad_Ha:
                    bad_Ha = y
            else:
                # 计算插值
                imgData[y][x] = Interpolation(
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
        for y in range(int(height_Fe / bin_count)):
            # 移动到第一个大于该坐标的地方
            while now < int(height_Fe / bin_count) - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now] and x > 400 / bin_count and x < W - 400 / bin_count:
                imgData[y + int(height_Ha / bin_count)][x] = stdy[now]
                if y < bad_Fe:
                    bad_Fe = y
            else:
                # 计算插值
                imgData[y + int(height_Ha / bin_count)][x] = Interpolation(
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
    # if bad_Ha > int(height_Ha / bin_count) - int(24 / bin_count):
    bad_Ha = int(height_Ha / bin_count) - int(24 / bin_count)
    # if bad_Fe > int(height_Fe / bin_count) - int(24 / bin_count):
    bad_Fe = int(height_Fe / bin_count) - int(24 / bin_count)
    # print(bad_Ha,bad_Fe,int(height_Ha / bin_count) - int(24 / bin_count),int(height_Fe / bin_count) - int(24 / bin_count))
    # 删除坏行 并输出两窗口最后的行数
    imgData[bad_Ha:bad_Ha + int(height_Fe / bin_count)] = imgData[
                                                          int(height_Ha / bin_count):int(height_Fe / bin_count) + int(
                                                              height_Ha / bin_count)]
    return imgData[0:bad_Ha + bad_Fe], bad_Ha, bad_Fe


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


# 将bin=2的图转换为bin=1（测试时使用）
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
    imgData, HofHa, HofFe = curve_correction(imgData - darkDate, x0, C)
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


def test():
    matplotlib.rcParams['font.sans-serif'] = ['KaiTi']
    filepath_result = "testResult/"
    filepath_test = "testData/"
    filepath_bash = "bass2000.txt"

    # print(base)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file), dtype=float)
    dark_data = change(dark_data)
    image_file = get_pkg_data_filename(filepath_test + 'for_flat_binning2.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)
    # RSM20211222T215254-0010-2313-基准.fts    RSM20211222T215555-0013-2367-测试.fts
    # RSM20220120T062536-0017-1081.fts
    H, W = flat_data.shape
    print(H, W)
    print(bin_count)
    flat_data, b, d = curve_correction(flat_data - dark_data, x0, C)
    flat_data = getFlat(flat_data)
    filename = filepath_test + 'RSM20220808T005114-0000-1837.fts'
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    imgData = moveImg(imgData, -1)
    imgData, HofHa, HofFe = curve_correction(imgData - dark_data, x0, C)
    print("pic1")
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flat_data)
    base = get_Sunstd(filepath_bash)
    filepathHA = "HA_absorption.txt"
    filepathFE = "FE_absorption.txt"
    abortion = get_Absorstd(filepathHA, filepathFE, HofHa, HofFe)
    plt.figure()
    plt.plot(base)
    plt.show()
    plt.figure()
    plt.imshow(flat_data, cmap="gray", aspect='auto')
    plt.show()
    # filelist = os.listdir(filepath_test)
    image_file, imgData = entireWork(filepath_test + 'RSM20220808T005114-0000-1837.fts', dark_data, flat_data, abortion)
    #
    print("OK")
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()

    # grey = fits.PrimaryHDU(image_file)


# 灰度图添加文字
# I_array为目标图 time_txt为目标文字
# 字体和粗细可通过config调整
def add_time(Input_array, time_txt):
    I_array = np.array(Input_array)
    H, W = I_array.shape
    h1 = int(H * 0.55)
    w1 = int(W * 0.95)
    cv2.putText(I_array, time_txt, (h1, w1), cv2.FONT_HERSHEY_PLAIN, txt_size, (I_array.max()), txt_thick)
    return I_array


def cal_center_mean(Input_np):
    data = Input_np.mean(axis=1).mean(axis=1)
    # print(data.shape)
    height_ha = int(height_Ha / bin_count) - int(24 / bin_count)
    height_fe = int(height_Fe / bin_count) - int(24 / bin_count)
    Ha_B = (data[Ha_Lower] + data[Ha_Upper]) / 2
    print(height_ha,Fe_Upper)
    Fe_B = (data[height_ha + Fe_Lower] + data[height_ha + Fe_Upper]) / 2
    Ha_U = Ha_L = 0
    Fe_U = Fe_L = 0

    for i in range(Ha_Lower, Ha_Upper + 1):
        Ha_U = Ha_U + (data[i] - Ha_B) * i
        Ha_L = Ha_L + (data[i] - Ha_B)
    for i in range(Fe_Lower, Fe_Upper + 1):
        Fe_U = Fe_U + (data[height_ha + i] - Fe_B) * i
        Fe_L = Fe_L + (data[height_ha + i] - Fe_B)

    Ha_mean = Ha_U / Ha_L
    Fe_mean = Fe_U / Fe_L + (5430 - 5096) / bin_count

    mean_K = (FE_lineCore - HA_lineCore) / (Fe_mean - Ha_mean)
    mean_b = FE_lineCore - Fe_mean * mean_K
    print(Fe_mean,Ha_mean)
    return mean_K, mean_b, mean_b + mean_K * (5430 - 5096) / bin_count


# log所用函数
def log(*args):
    print(datetime.datetime.now().strftime("[%Y-%m-%d-%H:%M:%S]"), end='')
    print(args)


if __name__ == "__main__":
    A = config.FE_start - config.HA_start
    B = A / K
    C = (5430 - 5096) / bin_count
    print(A, B, C)
    A = height_Ha
    B = 5430 - 5096
    print(A, B)
    Q =  np.array([1,2,3,4,5,6,7])
    D = 3
    print(Q[3:int(D)+3])
    print(cal_center_mean(np.zeros((164,3,3))))
    # cal_center_mean(np.array([[[1,2],[3,4]],[[5,6],[7,8]]]))
    # test()
    # I = Image.open("123.png")

    # I_array = np.array(I.convert('L'))
    # H, W = I_array.shape
    # text1 = str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    # I_array = add_time(I_array, text1)
    # plt.figure()
    # plt.imshow(I_array, cmap="gray")
    # plt.show()
    # testPath = "circle/circle/"
    # type = "check"
    # if type == "test":
    #     Filelist = os.listdir(testPath)
    #     if True:
    #         id = 105
    #         Filelist = os.listdir(testPath)
    #         I = Image.open(testPath + Filelist[1 + id])
    #         I_array = np.array(I.convert('L'))
    #
    #         # image_file = get_pkg_data_filename(testPath + 'sum8.fts')
    #         # I_array = np.array(fits.getdata(image_file), dtype=float)
    #         # # print(np.shape(I_array))
    #         rx, ry, r = getCircle(I_array)
    #         print(id, rx, ry, r)
    #         H, W = I_array.shape
    #         for i in range(H):
    #             for j in range(W):
    #                 if abs((i - rx) * (i - rx) + (j - ry) * (j - ry) - r * r) < 10000:
    #                     I_array[i][j] = 240
    #         # for point in points:
    #         #     for i in range(20):
    #         #         for j in range(20):
    #         #             I_array[point[0]+i-8][point[1]+j-8] = 240
    #         # print(rx, ry, r * 0.52)
    #         plt.figure()
    #         plt.imshow(I_array)
    #         plt.show()
    # H,W = I_array.shape
    # print(rx,ry,r)
    # for i in range(H):
    #     for j in range(W):
    #         if abs((i-rx)*(i-rx) + (j-ry)*(j-ry) -r*r) <10000:
    #             I_array[i][j]=240
    # plt.figure()
    # plt.imshow(I_array)
    # plt.show()
    # if_first_print = True
    # for i in range(100):
    #     remaining_count = mp.Value('i', int(i))
    #     if if_first_print:
    #         print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value), end='')
    #         if_first_print = False
    #     else:
    #         print('\b' * (5 + len(str(remaining_count)) + 1 + len(str(file_count.value))) + '当前进度:' + str(
    #             remaining_count.value) + '/' + str(file_count.value), end='')
    #     time.sleep(0.5)
    # test()
    # if type == "check":
    #     Filelist = os.listdir(testPath)
    #     print(Filelist)
    #     L = len(Filelist) - 1
    #     err = []
    #     id = 0
    #     while id < L:
    #         I = Image.open(testPath + Filelist[id + 1])
    #         I_array = np.array(I.convert('L'))
    #
    #         # image_file = get_pkg_data_filename(testPath + 'sum8.fts')
    #         # I_array = np.array(fits.getdata(image_file), dtype=float)
    #         # # print(np.shape(I_array))
    #         ry, rx, r = getCircle(I_array, id)
    #         print(id, rx, ry, r)
    #         H, W = I_array.shape
    #         for i in range(H):
    #             for j in range(W):
    #                 if abs((i - rx) * (i - rx) + (j - ry) * (j - ry) - r * r) < 10000 / bin_count:
    #                     I_array[i][j] = 240
    #         plt.imsave("Result/result/" + str(id) + ".jpg", I_array)
    #         # for point in points:
    #         #     for i in range(20):
    #         #         for j in range(20):
    #         #             I_array[point[0]+i-8][point[1]+j-8] = 240
    #         # print(rx, ry, r * 0.52)=
    #         id += 1
    # test()
