import numpy as np
import re
import os
import time
import matplotlib
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
from PIL import Image
import math
import scipy.signal as signal
import config
import astropy
import jplephem
import datetime
import random
import numpy as np
from math import *
from astropy.io import fits
from astropy.time import Time
from astropy import coordinates
import multiprocessing as mp

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
            if y > stdx[now]:
                imgData[y][x] = stdx[now]
                if y < bad_Ha:
                    bad_Ha = y
            else:
                # 计算插值
                imgData[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                        y - stdx[now - 1])

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
            if y > stdx[now]:
                imgData[y + int(height_Ha / bin_count)][x] = stdx[now]
                if y < bad_Fe:
                    bad_Fe = y
            else:
                # 计算插值
                imgData[y + int(height_Ha / bin_count)][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (
                        stdx[now] - stdx[now - 1]) * (
                                                                     y - stdx[now - 1])
    if bad_Ha < int(height_Ha / bin_count) - int(29 / bin_count):
        bad_Ha = int(height_Ha / bin_count) - int(29 / bin_count)
    if bad_Fe < int(height_Fe / bin_count) - int(29 / bin_count):
        bad_Fe = int(height_Fe / bin_count) - int(29 / bin_count)
    # print(bad_Ha,bad_Fe)
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
        # 进行插值操作
        ansY.append(
            dataY[now - 1] + (dataY[now] - dataY[now - 1]) / (dataX[now] - dataX[now - 1]) * (stdx[i] - dataX[now - 1]))
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
    return imgData / flatData


# 图像中值平滑操作
# 参数
def MedSmooth(imgData, winSize=3):
    zero_range = 100
    if bin_count == 1:
        imgData = signal.medfilt(imgData, kernel_size=winSize)
    if bin_count == 2:
        imgData = signal.medfilt(imgData, kernel_size=winSize - 2)
    imgData[:, imgData.shape[1] - int(zero_range / bin_count): imgData.shape[1] - 1] = 0
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
    print("偏移量：", end='')
    print(mx, my)

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


def getBin(imgData):
    H, W = imgData.shape
    print(H, W, height_Ha + height_Fe)
    if H >= height_Ha + height_Fe:
        return 1
    return 2


def entireWork(filename, darkDate, flatData, abortion):
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    # imgData = change(imgData)
    imgData = moveImg(imgData, -2)
    imgData, HofHa, HofFe = curve_correction(imgData - darkDate, x0, C)
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flatData)
    # plt.figure()
    # plt.imshow(imgData, cmap="gray", aspect='auto')
    # plt.show()
    plt.figure()
    plt.plot(imgData[:, 2200].reshape(-1))
    imgDataRB = RB_repair(imgData, abortion)
    imgDataRB = MedSmooth(imgDataRB, 3)
    plt.plot(imgDataRB[:, 2200].reshape(-1))
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
def getCircle(image):
    # 二值化
    image_max = np.max(image)
    image = np.clip((image - image_max * 0.05) * 10, 0, 1)
    # plt.figure()
    # plt.imshow(image)
    # plt.show()

    # 通过卷积，使用Sobel算子提取边界
    # conv1 = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    conv2 = np.array([[1, 1, 1], [1, -8, 1], [1, 1, 1]])
    conv3 = np.array([[1, 2, 1], [2, 4, 2], [1, 2, 1]]) / 16
    image = signal.convolve2d(image, conv3, "valid")
    # gradient_Y = signal.convolve2d(image, conv1, "valid")
    gradient_X = signal.convolve2d(image, conv2, "valid")
    gradient = np.abs(gradient_X)  # + np.abs(gradient_Y)
    gradient = np.clip(gradient, 0, np.max(gradient) * 0.6)
    gradient_max = np.max(gradient)
    gradient = np.clip((gradient - gradient * 0.05) * 10, 0, 1)

    # plt.figure()
    # plt.imshow(gradient)
    # plt.show()

    gradient_reverse = np.array(gradient)
    H, W = gradient.shape
    for i in range(W):
        gradient_reverse[:, i] = gradient[:, W - 1 - i]

    # plt.figure()
    # plt.imshow(gradient_reverse)
    # plt.show()

    # FFT变化
    imgFFT = np.fft.fft2(gradient)
    flatFFT = np.fft.fft2(gradient_reverse)

    # IFFT
    FR = imgFFT * np.conj(flatFFT)
    R = np.fft.ifft2(FR)
    R = np.fft.fftshift(R)
    # 获取最大值坐标
    pos = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    # print(pos)
    # print((pos[1] + W / 2) / 2)  # 460 340
    Offset = int(pos[1] - W / 2)
    if Offset > 0:
        gradient_reverse[:, Offset:] = gradient_reverse[:, 0:-Offset]
    else:
        Offset *= -1
        gradient_reverse[:, 0:-Offset] = gradient_reverse[:, Offset:]
    gradient *= gradient_reverse
    # plt.figure()
    # plt.imshow(gradient)
    # plt.show()
    points = []
    for i in range(H):
        for j in range(W):
            if gradient[i][j] > 0.8:
                points.append([i, j])
    L = len(points)
    if L < 200:
        print("圆检测失败")
        return -1,-1,-1
    flag = True
    times = 0
    while flag:
        p1, p2, p3 = random.sample(range(1, int(L / 2)), 3)
        # print(p1,p2,p3)
        y, x, r = circle(points[p1][1], points[p1][0], points[p2][1], points[p2][0], points[p3][1], points[p3][0])
        s = 0
        for i in range(L):
            if abs((points[i][0] - x) ** 2 + (points[i][1] - y) ** 2 - r ** 2) < 10000:
                s += 1
        times += 1
        if times > 100 * 2.5:
            return -1, -1, -1
        if s > L * 0.3 or (times > 100 and s > L * (0.3 - times / 1000)):
            flag = False
    # print(times)
    # print(x,y,r*0.52)
    return x + 4, y + 4, r - 20 / bin_count


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


def test():
    matplotlib.rcParams['font.sans-serif'] = ['KaiTi']
    filepath_result = "testResult/"
    filepath_test = "testData/"
    filepath_bash = "bass2000.txt"

    # print(base)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file), dtype=float)
    image_file = get_pkg_data_filename(filepath_test + 'for_flat.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)
    # RSM20211222T215254-0010-2313-基准.fts    RSM20211222T215555-0013-2367-测试.fts
    # RSM20220120T062536-0017-1081.fts
    H, W = flat_data.shape
    print(H, W)
    filelist = os.listdir(filepath_test)
    image_file = get_pkg_data_filename(filepath_test + 'RSM20211222T060119-0008-1353.fts')
    img_data = np.array(fits.getdata(image_file), dtype=float)
    # img_data = change(img_data)
    # bin = getBin(img_data)
    print(bin)
    img_data = moveImg(img_data, -2)
    # flat_data = change(flat_data)
    dark_data = change(dark_data)
    flat_data, b, d = curve_correction(flat_data - dark_data, x0, C)
    img_data, HofHa, HofFe = curve_correction(img_data - dark_data, x0, C)
    flat_data = getFlatOffset(flat_data, img_data)
    # print(flat_data)
    flat_data = getFlat(flat_data)

    filename = filepath_test + 'RSM20211222T060119-0008-1353.fts'
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    # imgData = change(imgData)
    imgData = moveImg(imgData, -2)
    imgData, HofHa, HofFe = curve_correction(imgData - dark_data, x0, C)
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
    image_file, imgData = entireWork(filepath_test + 'RSM20211222T060119-0008-1353.fts', dark_data, flat_data, abortion)
    #
    print("OK")
    plt.figure()
    plt.imshow(image_file, cmap="gray", aspect='auto')
    plt.show()

    # grey = fits.PrimaryHDU(image_file)


# 全局进度控制
file_count = mp.Value('i', 20)
if_first_print = mp.Value('b', True)
remaining_count = mp.Value('i', 0)

if __name__ == "__main__":
    testPath = "circle/circle/"
    type = "check"
    if type == "test":
        Filelist = os.listdir(testPath)
        I = Image.open(testPath + Filelist[178])
        I_array = np.array(I.convert('L'))
        # image_file = get_pkg_data_filename(testPath + 'sum8.fts')
        # I_array = np.array(fits.getdata(image_file), dtype=float)
        # # print(np.shape(I_array))
        rx, ry, r = getCircle(I_array)
        H, W = I_array.shape
        for i in range(H):
            for j in range(W):
                if abs((i - rx) * (i - rx) + (j - ry) * (j - ry) - r * r) < 10000:
                    I_array[i][j] = 240
        print(rx, ry, r * 0.52)
        plt.figure()
        plt.imshow(I_array)
        plt.show()
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
    else:
        Filelist = os.listdir(testPath)
        print(Filelist)
        L = len(Filelist) - 1
        err = []
        i = 0
        while i < L:
            I = Image.open(testPath + Filelist[i + 1])
            I_array = np.array(I.convert('L'))
            print("测试图像" + str(i) + ":", end="")
            rx, ry, r = getCircle(I_array)
            print(rx, ry, r * 0.52 * 2)
            if r * 0.52 > 985 or r * 0.52 < 955:
                err.append([i, Filelist[i + 1], r * 0.52])
            i += 1
        print(err)
