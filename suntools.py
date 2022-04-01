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
HA = config.HA  # 红蓝移HA参数
FE = config.FE  # 红蓝移FE参数
K = config.K  # 红蓝移K参数
bin = config.bin_count


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
    bad_Ha = height_Ha / bin
    bad_Fe = height_Fe / bin

    # 进行矫正操作
    # 分两个窗口分别操作
    for x in range(W):
        # 对于H窗口进行操作
        # 计算原坐标经过变换之后对应的坐标
        stdx = np.arange(0, int(height_Ha / bin), 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * K + HA) / (C * (x - x0) * (x - x0) + 1) - HA) / K
        # 获取原数据值
        stdy = imgData[:, x]
        # 确定插值的坐标
        now = 1
        for y in range(int(height_Ha / bin)):
            # 移动到第一个大于该坐标的地方
            while now < int(height_Ha / bin) - 1 and stdx[now] < y:
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
        stdx = np.arange(0, int(height_Fe / bin), 1)
        # 先转成波长 根据波长公式进行变换后反变化回坐标
        stdx = ((stdx * K + FE) / (C * (x - x0) * (x - x0) + 1) - FE) / K
        # 获取原数据值
        stdy = imgData[int(height_Ha / bin):int(height_Fe / bin) + int(height_Ha / bin), x]
        # 确定插值的坐标
        now = 1
        for y in range(int(height_Fe / bin)):
            # 移动到第一个大于该坐标的地方
            while now < int(height_Fe / bin) - 1 and stdx[now] < y:
                now += 1
            # 若越界则标记为坏点
            if y > stdx[now]:
                imgData[y + int(height_Ha / bin)][x] = stdx[now]
                if y < bad_Fe:
                    bad_Fe = y
            else:
                # 计算插值
                imgData[y + int(height_Ha / bin)][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (
                        stdx[now] - stdx[now - 1]) * (
                                                               y - stdx[now - 1])
    if bad_Ha < int(height_Ha / bin) - int(29 / bin):
        bad_Ha = int(height_Ha / bin) - int(29 / bin)
    if bad_Fe < int(height_Fe / bin) - int(29 / bin):
        bad_Fe = int(height_Fe / bin) - int(29 / bin)
    # 删除坏行 并输出两窗口最后的行数
    imgData[bad_Ha:bad_Ha + int(height_Fe / bin)] = imgData[
                                                    int(height_Ha / bin):int(height_Fe / bin) + int(height_Ha / bin)]
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
    # 提取所需要的对应数据
    for i in range(HofHa):
        ans[i] = (temp_std[i * bin] + temp_std[i * bin + bin - 1]) / 2
    for i in range(HofFe):
        ans[HofHa + i] = (temp_std[i * bin + height_Ha] + temp_std[i * bin + height_Ha + bin - 1]) / 2
    stdx = np.zeros(H)
    # 坐标转化为波长
    stdx[0:HofHa] = np.arange(0, HofHa, 1) * K * bin + HA
    stdx[HofHa:] = np.arange(0, HofFe, 1) * K * bin + FE
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
    return imgData / np.clip(flatData, 0.01, 200)


# 图像中值平滑操作
# 参数
def MedSmooth(imgData, winSize=4):
    imgData = signal.medfilt(imgData, kernel_size=int(winSize / 2) * 2 + 1)
    return imgData

def amplify(Data):
    H, W = Data.shape
    amplify_size = 10
    Data1 = np.zeros((H,(W-1)*amplify_size+1))
    for j in range(W-1):
        for k in range(amplify_size):
            Data1[:,j*amplify_size+k] = Data[:,j]*(1 - k/amplify_size) + Data[:,j + 1]*(k/amplify_size)
    Data1[:, (W-1)*amplify_size] = Data[:,W - 1]
    return Data1


# 计算偏差
def getFlatOffset(flatData, imgData):
    HofHa = int(height_Ha / bin)
    HofFe = int(height_Fe / bin)
    # 获取图片尺寸
    H, W = imgData.shape
    # 计算中心点

    cx = int(H / 2)
    cy = int(W / 2)
    # 获取序列
    # img = imgData[cx :cx + int(H / 4), cy - int(W / 8):cy + int(W / 8)]
    # flat = flatData[cx :cx + int(H / 4), cy - int(W / 8):cy + int(W / 8)]
    img = imgData[HofHa - int(H / 8):HofHa + int(H / 8), cy - int(W / 8):cy + int(W / 8) + 1]
    flat = flatData[HofHa - int(H / 8):HofHa + int(H / 8), cy - int(W / 8):cy + int(W / 8) + 1]
    img = amplify(img)
    flat = amplify(flat)
    # FFT变化
    imgFFT = np.fft.fft2(img)
    flatFFT = np.fft.fft2(flat)

    # plt.figure()
    # plt.plot(img[50, :].reshape(-1))
    # plt.plot(flat[50, :].reshape(-1))
    # plt.show()
    print(np.unravel_index(np.argmax(np.abs(img[10, :650])), img[10, :650].shape))
    print(np.unravel_index(np.argmax(np.abs(flat[10, :650])), flat[10, :650].shape))
    # IFFT
    FR = imgFFT * np.conj(flatFFT)
    R = np.fft.ifft2(FR)
    R = np.fft.fftshift(R)
    h,w = img.shape
    # plt.figure()
    # plt.imshow(np.abs(R), cmap="gray", aspect='auto')
    # plt.show()

    # 获取最大值坐标
    pos = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    print(pos)
    # 计算偏移量
    mx = int((pos[1] - int((w+1) / 2))/10)
    my = (pos[1] - int((w+1) / 2)) - mx * 10
    print(mx,my)

    # 偏移操作
    if mx * 10 + my < 0:
        my *= -1
        flatData[:, 0:W + mx - 1] = flatData[:, -mx:W - 1]*(1 - my/10) + flatData[:, -mx + 1:W]*(my/10)
    else:
        flatData[:, mx:W-1] = flatData[:, 0:W - mx-1]*(1 - my/10) + flatData[:, 1:W - mx]*my/10
    return flatData


# 对图像进行横向上的平移
def moveImg(imgdata, offset):
    H, W = imgdata.shape
    if offset < 0:
        imgdata[int(height_Ha / bin) + 1:, 0:W + int(offset / bin)] = imgdata[int(height_Ha / bin) + 1:, -int(offset / bin):W]
    else:
        imgdata[int(height_Ha / bin) + 1:, int(offset / bin):W] = imgdata[int(height_Ha / bin) + 1:, 0:W - int(offset / bin)]
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
    if bin == 1:
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
    imgData = change(imgData)
    bin = getBin(imgData)
    imgData = moveImg(imgData, -2 )
    imgData, HofHa, HofFe = curve_correction(imgData - darkDate, 2321.26, 1.92909e-011)
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flatData)
    # plt.figure()
    # plt.imshow(imgData, cmap="gray", aspect='auto')
    # plt.show()
    # plt.figure()
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

    # print(base)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file), dtype=float)

    image_file = get_pkg_data_filename(filepath_test + 'for_flat.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)
#RSM20211222T215254-0010-2313-基准.fts    RSM20211222T215555-0013-2367-测试.fts
    H, W = flat_data.shape
    filelist = os.listdir(filepath_test)
    image_file = get_pkg_data_filename(filepath_test + 'RSM20211222T215555-0013-2367-测试.fts')
    img_data = np.array(fits.getdata(image_file), dtype=float)
    img_data = change(img_data)
    # bin = getBin(img_data)
    print(bin)
    img_data = moveImg(img_data, -2)
    flat_data = change(flat_data)
    dark_data = change(dark_data)
    flat_data, b, d = curve_correction(flat_data - dark_data, 2321.26, 1.92909e-011)
    img_data, HofHa, HofFe = curve_correction(img_data - dark_data, 2321.26, 1.92909e-011)
    flat_data = getFlatOffset(flat_data, img_data)
    # print(flat_data)
    flat_data = getFlat(flat_data)

    filename = filepath_test + 'RSM20211222T215555-0013-2367-测试.fts'
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    imgData = change(imgData)
    imgData = moveImg(imgData, -2)
    imgData, HofHa, HofFe = curve_correction(imgData - dark_data, 2321.26, 1.92909e-011)
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flat_data)
    base = get_Sunstd(filepath_bash)
    abortion = RB_getdata(imgData, base, HofHa, HofFe)

    plt.figure()
    plt.imshow(flat_data, cmap="gray", aspect='auto')
    plt.show()

    # filelist = os.listdir(filepath_test)
    image_file, imgData = entireWork(filepath_test + 'RSM20211222T215555-0013-2367-测试.fts', dark_data, flat_data, abortion)
    #
    plt.figure()
    plt.imshow(image_file, cmap="gray", aspect='auto')
    plt.show()

    # grey = fits.PrimaryHDU(image_file)
