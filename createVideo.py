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
import config
import datetime

framePerSec = 10
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


def createVideo(fileDate: datetime.datetime):
    # 检查输出文件夹是否存在 不存在则创建
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    getPNGList()
    frameShape = cv2.imread(fileDir + ha_list[0]).shape
    ha_videoOut = cv2.VideoWriter(saveDir + 'RSM' + fileDate.strftime('%Y-%m-%d') + '_HA.avi',
                                  cv2.VideoWriter_fourcc(*'XVID'), framePerSec, (frameShape[1], frameShape[0]), True)
    # fe_videoOut = cv2.VideoWriter(saveDir + 'RSM' + fileDate.strftime('%Y-%m-%d') + '_FE.avi',
    #                               cv2.VideoWriter_fourcc(*'XVID'), framePerSec, (frameShape[1], frameShape[0]), True)
    for ha in ha_list:
        ha_img = cv2.imread(fileDir + ha)
        for i in range(framePerSec):
            ha_videoOut.write(ha_img)  # 写入对应帧(1s)
    # for fe in fe_list:
    #     fe_img = cv2.imread(fileDir + fe)
    #     for i in range(framePerSec):
    #         fe_videoOut.write(fe_img)
    ha_videoOut.release()
    # fe_videoOut.release()


def main():
    createVideo(datetime.datetime.now())


if __name__ == '__main__':
    main()
