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
import config

# 定义参数
height_Ha = config.height_Ha  # ha窗口的长度
height_Fe = config.height_Fe  # he窗口的长度
HA = config.HA   # 红蓝移HA参数
FE = config.FE   # 红蓝移FE参数
K = config.K     # 红蓝移K参数


# 谱线矫正
# 参数data: 图像数据(numpy标准格式, 二维数组)
# 参数x0: 曲线矫正对称轴
# 参数C: 曲线矫正二次系数
# 输出：矫正后的图像数据
def curve_correction(imgData, x0, C):
    # 获取图片高度和宽度
    H, W = imgData.shape
    x0 = x0 - 1
    # 定义两个窗口的高度 H窗口height_Ha行 Fe窗口height_Fe行
    bad_Ha = height_Ha
    bad_Fe = height_Fe

    # 进行矫正操作
    # 分两个窗口分别操作
    for x in range(W):
        # 对于H窗口进行操作
        # 计算原坐标经过变换之后对应的坐标
        stdx = np.arange(0, height_Ha, 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * K + HA) / (C * (x - x0) * (x - x0) + 1) - HA) / K
        # 获取原数据值
        stdy = imgData[:, x]
        # 确定插值的坐标
        now = 1
        for y in range(height_Ha):
            # 移动到第一个大于该坐标的地方
            while now < height_Ha - 1 and stdx[now] < y:
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
        stdx = np.arange(0, height_Fe, 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * K + FE) / (C * (x - x0) * (x - x0) + 1) - FE) / K
        # 获取原数据值
        stdy = imgData[height_Ha:height_Fe + height_Ha, x]
        # 确定插值的坐标
        now = 1
        for y in range(height_Fe):
            # 移动到第一个大于该坐标的地方
            while now < height_Fe - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now]:
                imgData[y + height_Ha][x] = stdx[now]
                if y < bad_Fe:
                    bad_Fe = y
            else:
                # 计算插值
                imgData[y + height_Ha][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (
                            stdx[now] - stdx[now - 1]) * (
                                                    y - stdx[now - 1])
    if bad_Ha < height_Ha - 29:
        bad_Ha = height_Ha - 29
    if bad_Fe < height_Fe - 29:
        bad_Fe = height_Fe - 29
    # 删除坏行 并输出两窗口最后的行数
    imgData[bad_Ha:bad_Ha + height_Fe] = imgData[height_Ha:height_Fe + height_Ha]
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
    stdx = np.zeros(height_Fe + height_Ha)
    # 根据窗口转换成需要的数据
    # 坐标变换
    stdx[0:height_Ha] = np.arange(0, height_Ha, 1) * K + HA
    stdx[height_Ha:] = np.arange(0, height_Fe, 1) * K + FE
    # 插值
    for i in range(height_Fe + height_Ha):
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


# 红蓝移矫正
def RB_getdata(imgData, sun_std, HofHa, HofFe):
    temp_std = np.array(sun_std)
    # 获取图片尺寸
    H, W = imgData.shape
    # 获取波长强度和
    sun_image = np.sum(imgData, axis=1) / W
    # 归一化
    # print(np.max(sun_image))
    sun_image /= np.max(sun_image)
    # 提取所需要的对应数据
    # temp_std[HofHa:HofHa + HofFe] = temp_std[height_Ha:height_Ha + HofFe]
    stdx = np.zeros(H)
    # 坐标转化为波长
    stdx[0:HofHa] = np.arange(0, HofHa, 1) * K + HA
    stdx[HofHa:] = np.arange(0, HofFe, 1) * K + FE
    # 拟合一次函数
    cov = np.polyfit(stdx, sun_image / temp_std[0:H], 1)
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
    return imgData / np.clip(flatData, 5e-5, 200)


# 图像中值平滑操作
# 参数
def MedSmooth(imgData, winSize=4):
    imgData = signal.medfilt(imgData, kernel_size=int(winSize / 2) * 2 + 1)
    return imgData


# 计算偏差
def getFlatOffset(flatData, imgData):
    # 获取图片尺寸
    H, W = imgData.shape
    # 计算中心点

    cx = int(H / 2)
    cy = int(W / 2)
    # 获取序列
    img = imgData[cx - int(H / 8):cx + int(H / 8), cy - int(W / 8):cy + int(W / 8)]
    flat = flatData[cx - int(H / 8):cx + int(H / 8), cy - int(W / 8):cy + int(W / 8)]

    # FFT变化
    imgFFT = np.fft.fft2(img)
    flatFFT = np.fft.fft2(flat)

    # IFFT
    FR = imgFFT * np.conj(flatFFT)
    R = np.fft.ifft2(FR)
    R = np.fft.fftshift(R)
    # plt.figure()
    # plt.imshow(np.abs(R), cmap="gray", aspect='auto')
    # plt.show()

    # 获取最大值坐标
    pos = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    # 计算偏移量
    mx = pos[1] - int(W / 8)
    print(mx)

    # 偏移操作
    if mx < 0:
        flatData[:, 0:W + mx] = flatData[:, -mx:W]
    else:
        flatData[:, mx:W] = flatData[:, 0:W - mx]

    return flatData


# 对图像进行横向上的平移
def moveImg(imgdata, offset):
    H, W = imgdata.shape
    if offset < 0:
        imgdata[height_Ha:, 0:W + offset] = imgdata[height_Ha:, -offset:W]
    else:
        imgdata[height_Ha:, offset:W] = imgdata[height_Ha:, 0:W - offset]
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


def entireWork(filename, darkDate, flatData, abortion):
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    imgData = moveImg(imgData, -2)
    imgData, HofHa, HofFe = curve_correction(imgData - darkDate, 2321.26, 1.92909e-011)
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flatData)
    plt.figure()
    plt.plot(imgData[:, 2200].reshape(-1))
    imgDataRB = RB_repair(imgData, abortion)
    imgDataRB = MedSmooth(imgDataRB, 3)
    plt.plot(imgDataRB[:, 2200].reshape(-1))
    plt.show()
    return imgDataRB, imgData


if __name__ == "__main__":
    matplotlib.rcParams['font.sans-serif'] = ['KaiTi']
    filepath_result = "testResult/"
    filepath_test = "testData/"
    filepath_bash = "bass2000.txt"

    base = get_Sunstd(filepath_bash)
    # print(base)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file), dtype=float)

    image_file = get_pkg_data_filename(filepath_test + 'for_flat.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)
    H, W = flat_data.shape
    # print(H,W)
    filelist = os.listdir(filepath_test)
    image_file = get_pkg_data_filename(filepath_test + filelist[2314])
    img_data = np.array(fits.getdata(image_file), dtype=float)
    flat_data = getFlatOffset(flat_data, img_data)
    flat_data, b, d = curve_correction(flat_data - dark_data, 2321.26, 1.92909e-011)
    # print(flat_data)
    flat_data = getFlat(flat_data)

    filename = filepath_test + filelist[2314]
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    imgData = moveImg(imgData, -2)
    imgData, HofHa, HofFe = curve_correction(imgData - dark_data, 2321.26, 1.92909e-011)
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flat_data)
    abortion = RB_getdata(imgData, base, HofHa, HofFe)

    # filelist = os.listdir(filepath_test)
    image_file, imgData = entireWork(filepath_test + filelist[631], dark_data, flat_data, abortion)
    # flat_data = np.zeros(dark_data.shape)
    # for i in range (200):
    #     image_file = get_pkg_data_filename(filepath_test + filelist[2712+i])
    #     flat_data += np.array(fits.getdata(image_file), dtype=float)
    # data = curve_correction(flat_data/400 - dark_data, 2321.26, 1.92909e-011)
    # data = smooth(data)
    # plt.figure()
    # plt.imshow(flat_data, cmap="gray", aspect='auto')
    # plt.show()
    # data, H, F = curve_correction(flat_data, 2321.26, 1.92909e-011)
    # plt.figure()
    # plt.imshow(data, cmap="gray",aspect='auto')
    # plt.show()
    # data = getFlat(data)
    # plt.figure()
    # plt.imshow(data, cmap="gray",aspect='auto')
    # plt.show()
    # print("Ping is over")
    # time_start = time.time()
    # image_file = entireWork(filepath_test + filelist[3000], dark_data, data, base)

    # time_end = time.time()
    # print(time_end - time_start)

    # plt.figure()
    # plt.plot(image_file[:,2200].reshape(-1))
    # plt.show()

    plt.figure()
    plt.imshow(image_file, cmap="gray", aspect='auto')
    plt.show()

    grey = fits.PrimaryHDU(image_file)
    greyHDU = fits.HDUList([grey])
    greyHDU.writeto(filepath_result + 'result.fits')
