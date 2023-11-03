"""
此py文件用来封装生成视频的函数
主要用途: 读取当前文件夹内的所有PNG文件, 将其分为两组(FE窗口图像与HA窗口图片)
同时生成两个视频, 每张图片占据视频的10帧
视频帧率固定为 10 frame/s, 格式固定为avi
视频像素大小固定为config内的序列长度与图像宽度

@author: seu_wxy
"""

import cv2
import os
import numpy as np
import glob

from PIL import Image

import config
import suntools
import traceback

from math import *
from astropy.io import fits
import matplotlib.pyplot as plt



framePerSec = config.frame_pre_sec
pic_count = config.write_to_video_count
pic_bias = config.write_to_video_bias
fileDir = config.save_dir_path
saveDir = config.video_dir_path
ha_list = []
fe_list = []


def getPNGList():
    arr = os.listdir(fileDir)
    for filename in arr:
        if filename.split('.')[-1] == 'png':
            if filename.split('_')[-1].split('.')[0] == 'HA':
                ha_list.append(filename)
            if filename.split('_')[-1].split('.')[0] == 'FE':
                fe_list.append(filename)
    ha_list.sort()
    fe_list.sort()
# 数组补全函数，将裁剪后的数组不足1040长度的元素用0填充
def completeArray(list_args):
    my_len = 2080
    result = []

    for my_list in list_args:
        if len(my_list) < my_len:
            for i in range(my_len - len(my_list)):
                my_list.append(0)
        result.append(my_list)
    return result

def createVideo():
    # 检查输出文件夹是否存在 不存在则创建
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    getPNGList()
    ha_videoOut = cv2.VideoWriter(saveDir + 'RSM' + fileDate.strftime('%Y-%m-%d') + '_HA.avi',
                                  cv2.VideoWriter_fourcc(*'XVID'), framePerSec, (2080, 2080), True)
    # fe_videoOut = cv2.VideoWriter(saveDir + 'RSM' + fileDate.strftime('%Y-%m-%d') + '_FE.avi',
    #                               cv2.VideoWriter_fourcc(*'XVID'), framePerSec, (frameShape[1], frameShape[0]), True)
    for cnt in range(len(ha_list)):
        ha_img = cv2.imread(fileDir + ha_list[cnt])
        im = Image.open(fileDir + ha_list[cnt])
        metadata = im.info
        centerx = metadata['CenterX']
        centery = metadata['CenterY']
        ha_img = ha_img[int(centery)-1040:int(centery)+1040,int(centerx)-1040:int(centerx)+1040]
        # completeArray(ha_img)
        bias_tmp = cnt % (pic_bias + 1)
        if bias_tmp == 0:
            for i in range(pic_count):
                ha_videoOut.write(ha_img)  # 写入对应帧(1s)
    # for fe in fe_list:
    #     fe_img = cv2.imread(fileDir + fe)
    #     for i in range(framePerSec):
    #         fe_videoOut.write(fe_img)
    ha_videoOut.release()
    # fe_videoOut.release()

def createVideoMp4(img_array,filename):
    size = (0,0)
    for indexf in range(len(img_array)):  # 这个循环是为了读取所有要用的图片文件
        img1 = cv2.imread(fileDir + img_array[indexf], 1)
        if img1 is None:
            print(img_array[indexf] + " is error!")
            continue
        size = (img1.shape[1], img1.shape[0])
        img_array.append(img1)

    videowrite = cv2.VideoWriter(filename, -1, 6, size)

    for i in range(len(img_array)):  # 把读取的图片文件写进去
        videowrite.write(img_array[i])

    videowrite.release()

def createVideoNJU():
    getPNGList()
    createVideoMp4(ha_list,saveDir + 'RSM' + fileDate.strftime('%Y-%m-%d') + '_HA.mp4')
    createVideoMp4(fe_list, saveDir + 'RSM' + fileDate.strftime('%Y-%m-%d') + '_FE.mp4')

def monographNJU(hadata_path0,fedata_path0):
    ang_res = 0.5218 * 2
    tick_pixel = [1040 - 1000 / ang_res, 1040 - 500 / ang_res, 1040, 1040 + 500 / ang_res, 1040 + 1000 / ang_res]
    tick_arcsec = [-1000, -500, 0, 500, 1000]

    for rsm_path in hadata_path0:
        hdu = fits.open(rsm_path)
        hacore = hdu[1].data[68, :, :]
        cx = int(hdu[1].header['CRPIX1'])
        cy = int(hdu[1].header['CRPIX2'])
        datetime = rsm_path.split('/')[-1][3:7] + '-' + rsm_path.split('/')[-1][7:9] + '-' + rsm_path.split('/')[-1][
                                                                                             9:11] + \
                   ' UT ' + rsm_path.split('/')[-1][12:14] + ':' + rsm_path.split('/')[-1][14:16] + ':' + \
                   rsm_path.split('/')[-1][16:18]
        plt.figure(figsize=(20, 20))
        plt.imshow(hacore[cy - 1040:cy + 1040, cx - 1040:cx + 1040], origin='lower', vmin=0, vmax=4 * hacore.mean(),
                   cmap='afmhot')
        plt.text(1950, 50, datetime, horizontalalignment='right', verticalalignment='bottom', color='white', size=28)
        plt.xticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.yticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.xlabel('Solar X (arcsec)', fontsize=30)
        plt.ylabel('Solar Y (arcsec)', fontsize=30)
        plt.savefig('/data/home/MikeRao/data/Time_array/ha' + rsm_path.split('/')[-1][3:18] + 'full_sun.png', dpi=100)
        plt.close()
        print('/data/home/MikeRao/data/Time_array/ha' + rsm_path.split('/')[-1][3:18] + 'full_sun.png')

    ang_res = 0.5218 * 2
    tick_pixel = [1040 - 1000 / ang_res, 1040 - 500 / ang_res, 1040, 1040 + 500 / ang_res, 1040 + 1000 / ang_res]
    tick_arcsec = [-1000, -500, 0, 500, 1000]

    for rsm_path in fedata_path0:
        hdu = fits.open(rsm_path)
        hacore = hdu[1].data[10, :, :]
        cx = int(hdu[1].header['CRPIX1'])
        cy = int(hdu[1].header['CRPIX2'])
        datetime = rsm_path.split('/')[-1][3:7] + '-' + rsm_path.split('/')[-1][7:9] + '-' + rsm_path.split('/')[-1][
                                                                                             9:11] + \
                   ' UT ' + rsm_path.split('/')[-1][12:14] + ':' + rsm_path.split('/')[-1][14:16] + ':' + \
                   rsm_path.split('/')[-1][16:18]
        plt.figure(figsize=(20, 20))
        plt.imshow(hacore[cy - 1040:cy + 1040, cx - 1040:cx + 1040], origin='lower', vmin=0, vmax=4 * hacore.mean(),
                   cmap='afmhot')
        plt.text(1950, 50, datetime, horizontalalignment='right', verticalalignment='bottom', color='white', size=28)
        plt.xticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.yticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.xlabel('Solar X (arcsec)', fontsize=30)
        plt.ylabel('Solar Y (arcsec)', fontsize=30)
        plt.savefig('/data/home/MikeRao/data/Time_array/fe' + rsm_path.split('/')[-1][3:18] + 'full_sun.png', dpi=100)
        plt.close()
        print('/data/home/MikeRao/data/Time_array/fe' + rsm_path.split('/')[-1][3:18] + 'full_sun.png')



def main():
    createVideo()


if __name__ == '__main__':
    main()
