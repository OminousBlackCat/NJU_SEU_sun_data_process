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
    # 260 116分开做
    bad_H=260
    bad_Fe=116
    for x in range(W):
        stdx = np.arange(0, 260, 1)
        stdx = ((stdx * 0.024202301 + 6562.82)  / (C * (x - x0) * (x - x0) + 1) - 6562.82) / 0.024202301
        stdy = imgData[:, x]
        now = 1
        for y in range(260):
            while now < 260 - 1 and stdx[now] < y:
                now += 1
            # data[y][x] = max(stdy[now],stdy[now-1])
            if y > stdx[now]:
                imgData[y][x] = stdx[now]
                if y<bad_H:
                    bad_H = y
            else:
                imgData[y][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                            y - stdx[now - 1])
        stdx = np.arange(0, 116, 1)
        stdx = ((stdx * 0.024202301 + 6569.22)  / (C * (x - x0) * (x - x0) + 1) - 6569.22) / 0.024202301
        stdy = imgData[260:376, x]
        now = 1
        for y in range(116):
            while now < 116 - 1 and stdx[now] < y:
                now += 1
            # data[y][x] = max(stdy[now],stdy[now-1])
            if y > stdx[now]:
                imgData[y+260][x] = stdx[now]
                if y<bad_Fe:
                    bad_Fe = y
            else:
                imgData[y+260][x] = stdy[now - 1] + (stdy[now] - stdy[now - 1]) / (stdx[now] - stdx[now - 1]) * (
                            y - stdx[now - 1])
    imgData[bad_H:bad_H+116] = imgData[260:376]
    # print(bad_H,bad_Fe)# 230 80
    return imgData[0:bad_H+bad_Fe]


# 平场计算
# 参数flatData: 平场图像(numpy标准格式, 二维数组)
# 输出: 平场数据
def getFlat(flatData):
    H, W = flatData.shape
    mean = np.sum(flatData, axis=1) / W
    for i in range(H):
        flatData[i, :] = flatData[i, :] / mean[i]
    return flatData


# 标准太阳光谱获取
# 参数filepath：储存标准光谱的文件路径
# 输出：文件包含数据
def get_Sunstd(filepath):
    # print("标准光谱数据地址：" + filepath)
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
    now = 1
    ansY = []
    stdx = np.zeros(376)
    stdx[0:260] = np.arange(0, 260, 1) * 0.024202301 + 6562.82
    stdx[260:] = np.arange(0, 116, 1) * 0.024202301 + 6569.22
    for i in range(376):
        while dataX[now]<stdx[i] and now<len(dataX)-1 :
            now += 1
        ansY.append(dataY[now-1] + (dataY[now] - dataY[now-1]) / (dataX[now] - dataX[now-1]) * (stdx[i] - dataX[now-1]))
    ansY = np.array(ansY)
    return ansY/np.max(ansY)


# 红蓝移矫正
def RB_repair(imgData, sun_std):
    H, W = imgData.shape
    sun_image = np.sum(imgData, axis=1) / W
    sun_image /= np.max(sun_image)
    stdx = np.zeros(318)
    sun_std[231:231+87] = sun_std[260:260+87]
    stdx[0:231] = np.arange(0, 231, 1) * 0.024202301 + 6562.82
    stdx[231:] = np.arange(0, 87, 1) * 0.024202301 + 6569.22
    # print(H)
    cov = np.polyfit(stdx, sun_image / sun_std[0:H], 1)
    k, b = cov[0], cov[1]
    # print(k,b)
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
    OffsetData[offset:offset+H,offset:offset+W] = imgData[0:H,0:W]
    for i in range(winSize):
        for j in range(winSize):
            SmoothData[offset:offset+H,offset:offset+W] += OffsetData[i : H + i,j:W + j]
    return SmoothData[offset: offset + H, offset: offset + W]/winSize/winSize


def DivFlat(imgData,flatData):
    H, W = imgData.shape
    # flatList = np.sum(flatData,axis=0)
    # maxFlatindex = np.argmax(flatList)
    # imgList = np.sum(imgData,axis=0)
    # maxImgindex = np.argmax(imgList)
    # offset = maxImgindex - maxFlatindex
    # 230 80
    # print(offset)
    imgHa = imgData[0:230]
    flatHa = flatData[0:230]
    flatList = np.sum(flatHa, axis=0)
    maxFlatindex = np.argmax(flatList)
    imgList = np.sum(imgHa, axis=0)
    maxImgindex = np.argmax(imgList)
    offset = maxImgindex - maxFlatindex
    if offset < 0:
        offset *= -1
        flatHa[0:230 - offset, 0:W - offset] = flatHa[offset:230,offset:W]
    else:
        flatHa[offset:230,offset:W] = flatHa[0:230 - offset, 0:W - offset]
    imgFe = imgData[230:]
    flatFe = flatData[230:]
    flatList = np.sum(flatFe, axis=0)
    maxFlatindex = np.argmax(flatList)
    imgList = np.sum(imgFe, axis=0)
    maxImgindex = np.argmax(imgList)
    offset = maxImgindex - maxFlatindex
    if offset < 0:
        offset *= -1
        flatFe[0:80 - offset, 0:W - offset] = flatFe[offset:80, offset:W]
    else:
        flatFe[offset:80, offset:W] = flatFe[0:80 - offset, 0:W - offset]
    flatData[0:230] = flatHa
    flatData[230:] = flatFe
    return imgData/flatData


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

    # flat_data = np.zeros(dark_data.shape)
    # for i in range (200):
    #     image_file = get_pkg_data_filename(filepath_test + filelist[2712+i])
    #     flat_data += np.array(fits.getdata(image_file), dtype=float)
    # data = curve_correction(flat_data/400 - dark_data, 2321.26, 1.92909e-011)
    # data = smooth(data)
    # plt.figure()
    # plt.imshow(flat_data, cmap="gray", aspect='auto')
    # plt.show()
    data = curve_correction(flat_data , 2321.26, 1.92909e-011)
    # plt.figure()
    # plt.imshow(data, cmap="gray",aspect='auto')
    # plt.show()
    data = getFlat(data)
    # plt.figure()
    # plt.imshow(data, cmap="gray",aspect='auto')
    # plt.show()
    print("Ping is over")
    plt.figure()
    time_start = time.time()
    image_file = get_pkg_data_filename(filepath_test + filelist[3000])
    test_data = np.array(fits.getdata(image_file), dtype=float)
    plt.subplot(5, 1, 1)
    plt.imshow(test_data - dark_data, cmap="gray",aspect='auto')
    #plt.title('去暗场')
    test_data = curve_correction(test_data - dark_data, 2321.26, 1.92909e-011)
    plt.subplot(5, 1, 2)
    plt.imshow(test_data, cmap="gray",aspect='auto')
    #plt.title('谱线矫正')
    H, W = data.shape
    # data[2:H,2:W] = data[0:H-2,0:W-2]
    test_data = DivFlat(test_data,data)
    plt.subplot(5, 1, 3)
    plt.imshow(test_data, cmap="gray",aspect='auto')
    #plt.title('去平场')
    test_data = RB_repair(test_data, base)
    plt.subplot(5, 1, 4)
    test_data = np.array(test_data, dtype=np.int16)
    plt.imshow(test_data, cmap="gray",aspect='auto')
    #plt.title('红蓝翼矫正')
    time_end1 = time.time()
    test_data = MedSmooth(test_data)
    plt.subplot(5, 1, 5)
    plt.imshow(test_data, cmap="gray",aspect='auto')
    #plt.title('滤波')
    time_end = time.time()
    print(time_end - time_start)
    print(time_end - time_end1)
    plt.figure()
    plt.imshow(test_data, cmap="gray",aspect='auto')
    plt.show()

    grey = fits.PrimaryHDU(test_data)
    greyHDU = fits.HDUList([grey])
    greyHDU.writeto(filepath_result+'result.fits')
