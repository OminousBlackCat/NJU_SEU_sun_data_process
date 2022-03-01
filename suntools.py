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
    # 获取图片高度和宽度
    H, W = imgData.shape
    x0 = x0 - 1
    # 定义两个窗口的高度 H窗口260行 Fe窗口116行
    bad_H = 260
    bad_Fe = 116

    # 进行矫正操作
    # 分两个窗口分别操作
    for x in range(W):
        # 对于H窗口进行操作
        # 计算原坐标经过变换之后对应的坐标
        stdx = np.arange(0, 260, 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * 0.024202301 + 6562.82) / (C * (x - x0) * (x - x0) + 1) - 6562.82) / 0.024202301
        # 获取原数据值
        stdy = imgData[:, x]
        # 确定插值的坐标
        now = 1
        for y in range(260):
            # 移动到第一个大于该坐标的地方
            while now < 260 - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now]:
                imgData[y][x] = stdx[now]
                if y < bad_H:
                    bad_H = y
            else:
                # 计算插值
                imgData[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                        y - stdx[now - 1])

        # 对于Fe窗口进行操作
        # 计算原坐标经过变换之后对应的坐标
        stdx = np.arange(0, 116, 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * 0.024202301 + 6569.22) / (C * (x - x0) * (x - x0) + 1) - 6569.22) / 0.024202301
        # 获取原数据值
        stdy = imgData[260:376, x]
        # 确定插值的坐标
        now = 1
        for y in range(116):
            # 移动到第一个大于该坐标的地方
            while now < 116 - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now]:
                imgData[y + 260][x] = stdx[now]
                if y < bad_Fe:
                    bad_Fe = y
            else:
                # 计算插值
                imgData[y + 260][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                        y - stdx[now - 1])
    # 删除坏行 并输出两窗口最后的行数
    imgData[bad_H:bad_H + 116] = imgData[260:376]
    return imgData[0:bad_H + bad_Fe], bad_H, bad_Fe


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
    stdx = np.zeros(376)
    # 根据窗口转换成需要的数据
    # 坐标变换
    stdx[0:260] = np.arange(0, 260, 1) * 0.024202301 + 6562.82
    stdx[260:] = np.arange(0, 116, 1) * 0.024202301 + 6569.22
    # 插值
    for i in range(376):
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
def RB_repair(imgData, sun_std, HofH, HofFe):
    temp_std = np.array(sun_std)
    # 获取图片尺寸
    H, W = imgData.shape
    # 获取波长强度和
    sun_image = np.sum(imgData, axis=1) / W
    # 归一化
    #print(np.max(sun_image))
    sun_image /= np.max(sun_image)
    # 提取所需要的对应数据
    temp_std[HofH:HofH + HofFe] = temp_std[260:260 + HofFe]
    stdx = np.zeros(H)
    # 坐标转化为波长
    stdx[0:HofH] = np.arange(0, HofH, 1) * 0.024202301 + 6562.82
    stdx[HofH:] = np.arange(0, HofFe, 1) * 0.024202301 + 6569.22
    # 拟合一次函数
    cov = np.polyfit(stdx, sun_image / temp_std[0:H], 1)
    k, b = cov[0], cov[1]
    # 红蓝翼矫正
    for i in range(H):
        imgData[i, :] = imgData[i, :] / (k * stdx[i] + b)
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
    return imgData / np.clip(flatData,5e-5,200)


# 图像中值平滑操作
# 参数
def MedSmooth(imgData, winSize=4):
    imgData = signal.medfilt(imgData, kernel_size=int(winSize/2)*2+1)
    return imgData


def entireWork(filename, darkDate, flatData, sun_std):
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    imgData, HofH, HofFe = curve_correction(imgData - darkDate, 2321.26, 1.92909e-011)
    imgData = DivFlat(imgData, flatData)
    imgData = RB_repair(imgData, sun_std, HofH, HofFe)
    imgData = MedSmooth(imgData)
    return imgData

# 计算偏差
def getFlatOffset(flatData,imgData):
    # 获取图片尺寸
    H, W = imgData.shape
    # 计算中心点

    cx = int(H/2)
    cy = int(W/2)
    # 获取序列
    img = imgData[50:150,cy - 1000:cy + 1000]
    flat = flatData[50:150,cy - 1000:cy + 1000]
    Img = np.zeros((100,20000))
    Flat = np.zeros((100, 20000))
    # 图片拉伸
    for i in range (2000):
        for j in range(10):
            Img[:, i * 10 + j] = img[:, i]*(1 - j / 10)+img[:, i]* j / 10
            Flat[:, i * 10 + j] = flat[:, i] * (1 - j / 10) + flat[:, i] * j / 10

    # FFT变化
    imgFFT = np.fft.fft2(Img)
    flatFFT = np.fft.fft2(Flat)

    # IFFT
    FR = imgFFT*np.conj(flatFFT)
    R = np.fft.ifft2(FR)
    R = np.fft.fftshift(R)
    # plt.figure()
    # plt.imshow(np.abs(R), cmap="gray", aspect='auto')
    # plt.show()

    # 获取最大值坐标
    pos = np.unravel_index(np.argmax(np.abs(R)),R.shape)

    # 计算偏移量
    mx = int(pos[1] / 10) - 1000
    my = pos[1] - int(pos[1] / 10) * 10
    print(mx,my)

    # 偏移操作
    front = flatData
    next = flatData
    if mx < 0 :
        front[0:260,0:W + mx] = flatData[0:260,-mx:W]
        next[260:,-mx + 1:W] = flatData[260:,0:W + mx - 1]
        flatData[0:260] = front[0:260] + (next[0:260] - front[0:260])*my/10
    else:
        front[0:260, mx:W] = flatData[0:260, 0:W-mx]
        next[0:260, mx+1:W] = flatData[0:260, 0:W-mx-1]
        flatData[0:260] = front[0:260] + (next[0:260] - front[0:260])*my/10

    # 获取序列
    img = imgData[250:280,cy - 1000:cy + 1000]
    flat = flatData[250:280,cy - 1000:cy + 1000]
    Img = np.zeros((30, 20000))
    Flat = np.zeros((30, 20000))
    for i in range(2000):
        for j in range(10):
            Img[:, i * 10 + j] = img[:, i] * (1 - j / 10) + img[:, i] * j / 10
            Flat[:, i * 10 + j] = flat[:, i] * (1 - j / 10) + flat[:, i] * j / 10

    # FFT变化
    imgFFT = np.fft.fft2(Img)
    flatFFT = np.fft.fft2(Flat)

    # iFFT变化
    FR = imgFFT*np.conj(flatFFT)
    R = np.fft.ifft2(FR)
    R = np.fft.fftshift(R)
    # plt.figure()
    # plt.imshow(np.abs(R), cmap="gray", aspect='auto')
    # plt.show()
    pos = np.unravel_index(np.argmax(np.abs(R)),R.shape)
    print(pos)

    # 计算偏移量
    mx = int(pos[1] / 10) - 1000
    my = pos[1] - int(pos[1] / 10) * 10
    print(mx,my)
    front = flatData
    next = flatData

    # 偏移操作
    if mx < 0 :
        front[260:,0:W + mx] = flatData[260:,-mx:W]
        next[260:,-mx + 1:W] = flatData[260:,0:W + mx - 1]
        flatData[260:] = front[260:] + (next[260:] - front[260:])*my/10
    else:
        front[260:, mx:W] = flatData[260:, 0:W-mx]
        next[260:, mx+1:W] = flatData[260:, 0:W-mx-1]
        flatData[260:] = front[260:] + (next[260:] - front[260:])*my/10
    return flatData



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
    filelist = os.listdir(filepath_test)
    image_file = get_pkg_data_filename(filepath_test + filelist[2312])
    img_data = np.array(fits.getdata(image_file), dtype=float)
    flat_data = getFlatOffset(flat_data,img_data)
    flat_data,b,d = curve_correction(flat_data - dark_data, 2321.26, 1.92909e-011)
    #print(flat_data)
    flat_data = getFlat(flat_data)
    # filelist = os.listdir(filepath_test)
    image_file = entireWork(filepath_test + filelist[3000], dark_data, flat_data, base)
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
    plt.figure()
    plt.imshow(image_file, cmap="gray", aspect='auto')
    plt.show()

    grey = fits.PrimaryHDU(image_file)
    greyHDU = fits.HDUList([grey])
    greyHDU.writeto(filepath_result+'result.fits')
