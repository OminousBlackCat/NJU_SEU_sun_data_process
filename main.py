import numpy as np
import os
import time
from astropy.io import fits
import suntools
import matplotlib.pyplot as plt
from PIL import Image
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

mode = 'a'
#'a'谱线修复测试  'b'时间操作（只包含修复和合成） 'c'平场测试
filepath_result="testResult/"
filepath_test="testData/"

if mode=='a':
    image_file = get_pkg_data_filename('testData/RSM20211206T021759-0006-2088.fts')
    fits.info(image_file)
    # 读取fits图像数据（2D numpy 数组）
    image_data = fits.getdata(image_file)
    image_data = suntools.Curve_correction(image_data,2225,0.06/(2225-770)/(2225-770))

    plt.figure()
    plt.imshow(image_data, cmap="gray")
    plt.colorbar()
    plt.show()
    plt.imsave(filepath_result+'out_repair.jpg', image_data)

#12631.2s 未作红蓝翼矫正
if mode=='b':
    filePath = 'testData'
    N=len(os.listdir(filePath))
    filelist=os.listdir(filePath)
    now=0
    data=np.zeros((N,4608))
    time_start=time.time()

    for filename in os.listdir(filePath):
        image_file = get_pkg_data_filename(filePath+"/"+filename)
        image_data = fits.getdata(image_file)
        data[now,:]=suntools.Curve_correction(image_data, 2225, 0.06 / (2225 - 770) / (2225 - 770))[220,:]
        now+=1
        if(now>=N):
            break;
        print(str(now/N*100)+'%')
        time_end = time.time()
        print('time cost',time_end-time_start,'s')
        print('预计时间',(time_end-time_start)/now*N,'s')
    plt.figure()
    plt.imshow(data)
    plt.show()
    plt.imsave(filepath_result+'sun.jpg', data)


if mode=='c':
    for i in range(400):
        filename=filelist[2200+i]
        image_file = get_pkg_data_filename(filePath + "/" + filename)
        image_data = fits.getdata(image_file)
        print(filename)
        data=data+image_data
    N=400
    data=data/N
    mean=np.sum(data,axis=1)
    for i in range(376):
        data[i,:]=data[i,:]/mean[i]
    plt.figure()
    plt.imshow(data)
    plt.show()
    plt.imsave(filepath_result+'sun.jpg', data)

if mode=='d':
    I = Image.open(filepath_result+'sun.jpg')
    I_array = np.array(I)
    H,W,D = I_array.shape
    data=np.zeros((H,W))
    # plt.imshow(I_array, cmap="gray")
    # plt.show()

    for i in range(H):
        for j in range(W):
            data[i][j]=I_array[i][j][0]*0.3+I_array[i][j][1]*0.59+I_array[i][j][2]*0.11
            # I_array[i][j][0]=I_array[i][j][0]*0.3+I_array[i][j][1]*0.59+I_array[i][j][2]*0.11
            # I_array[i][j][1]=I_array[i][j][2]=0
    #print(I_array)
    plt.figure()
    plt.imshow(data, cmap="gray")
    plt.colorbar()
    plt.show()