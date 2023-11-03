"""
此py文件用来封装生成png预览图与预览视频的函数
主要用途: 读取当前文件夹内的所有fits文件, 将其分为两组(FE窗口图像与HA窗口图片)
随后生成对应的png图像，并且生成两个MP4预览视频

@origin: nju_
@author: seu_wxy
"""
from datetime import datetime
import cv2
import os
from astropy.io import fits
import matplotlib.pyplot as plt


def monographNJU(target_path: str, image_dpi = 100):
    """
    读取传入的目标路径, 列出目标文件夹内所有的ha fe文件
    再根据ha/fe fits文件列表在target_dir的文件夹内生成相应的png文件
    生成的png文件与fits文件的文件名相同, 仅格式不同
    @param target_path 传入的目标路径
    @param image_dpi 图像dpi, 默认100
    """
    # 读取文件路径, 将fit分为ha和fe两组存入list
    # list内 存的是完整路径
    ha_fits_list = []
    fe_fits_list = []
    file_arr = os.listdir(target_path)
    for filename in file_arr:
        file_format = filename.split('.')[-1]
        if file_format == 'fits' or file_format == 'fts' or file_format == 'FITS' or file_format == 'FTS':
            if filename.split('_')[-1].split('.')[0] == 'HA' or filename.split('_')[-1].split('.')[0] == 'ha':
                ha_fits_list.append(os.path.join(target_path, filename))
            if filename.split('_')[-1].split('.')[0] == 'FE' or filename.split('_')[-1].split('.')[0] == 'fe':
                fe_fits_list.append(os.path.join(target_path, filename))
    # 一些参数
    ang_res = 0.5218 * 2
    tick_pixel = [1040 - 1000 / ang_res, 1040 - 500 / ang_res, 1040, 1040 + 500 / ang_res, 1040 + 1000 / ang_res]
    tick_arcsec = [-1000, -500, 0, 500, 1000]

    for rsm_path in ha_fits_list:
        hdu = fits.open(rsm_path)
        hacore = hdu[1].data[68, :, :]
        cx = int(hdu[1].header['CRPIX1'])
        cy = int(hdu[1].header['CRPIX2'])
        current_filename = rsm_path.split('/')[-1]
        current_file_datetime = current_filename[3:7] + '-' + current_filename[7:9] + '-' + current_filename[9:11] + \
                   ' UT ' + current_filename[12:14] + ':' + current_filename[14:16] + ':' + current_filename[16:18]
        plt.figure(figsize=(20, 20))
        # plt.imshow(hacore[cy - 1040:cy + 1040, cx - 1040:cx + 1040], origin='lower', vmin=0, vmax=4 * hacore.mean(),
        #            cmap='afmhot')
        plt.text(1950, 50, current_file_datetime, horizontalalignment='right', verticalalignment='bottom', color='white', size=28)
        plt.xticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.yticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.xlabel('Solar X (arcsec)', fontsize=30)
        plt.ylabel('Solar Y (arcsec)', fontsize=30)
        plt.savefig(target_path + current_filename.split('.')[0]  + ".png", dpi=image_dpi)
        plt.close()
        # print('/data/home/MikeRao/data/Time_array/ha' + rsm_path.split('/')[-1][3:18] + 'full_sun.png')

    ang_res = 0.5218 * 2
    tick_pixel = [1040 - 1000 / ang_res, 1040 - 500 / ang_res, 1040, 1040 + 500 / ang_res, 1040 + 1000 / ang_res]
    tick_arcsec = [-1000, -500, 0, 500, 1000]

    for rsm_path in fe_fits_list:
        hdu = fits.open(rsm_path)
        hacore = hdu[1].data[10, :, :]
        cx = int(hdu[1].header['CRPIX1'])
        cy = int(hdu[1].header['CRPIX2'])
        current_filename = rsm_path.split('/')[-1]
        current_file_datetime = current_filename[3:7] + '-' + current_filename[7:9] + '-' + current_filename[9:11] + \
                                ' UT ' + current_filename[12:14] + ':' + current_filename[14:16] + ':' + current_filename[16:18]
        plt.figure(figsize=(20, 20))
        # plt.imshow(hacore[cy - 1040:cy + 1040, cx - 1040:cx + 1040], origin='lower', vmin=0, vmax=4 * hacore.mean(),
        #            cmap='afmhot')
        plt.text(1950, 50, current_file_datetime, horizontalalignment='right', verticalalignment='bottom', color='white', size=28)
        plt.xticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.yticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.xlabel('Solar X (arcsec)', fontsize=30)
        plt.ylabel('Solar Y (arcsec)', fontsize=30)
        plt.savefig(target_path + current_filename.split('.')[0] + ".png", dpi=image_dpi)
        plt.close()
        # print('/data/home/MikeRao/data/Time_array/fe' + rsm_path.split('/')[-1][3:18] + 'full_sun.png')


def createVideoMp4(target_dir: str, img_array: list, filename: str):
    size = (0,0)
    for indexf in range(len(img_array)):  # 这个循环是为了读取所有要用的图片文件
        img1 = cv2.imread(target_dir + img_array[indexf], 1)
        if img1 is None:
            print(img_array[indexf] + " is error!")
            continue
        size = (img1.shape[1], img1.shape[0])
        img_array.append(img1)

    videowrite = cv2.VideoWriter(filename, -1, 6, size)

    for i in range(len(img_array)):  # 把读取的图片文件写进去
        videowrite.write(img_array[i])

    videowrite.release()

def createVideoNJU(target_dir: str, save_dir: str, video_date: datetime):
    arr = os.listdir(target_dir)
    ha_list = []
    fe_list = []
    for filename in arr:
        if filename.split('.')[-1] == 'png':
            if filename.split('_')[-1].split('.')[0] == 'HA':
                ha_list.append(filename)
            if filename.split('_')[-1].split('.')[0] == 'FE':
                fe_list.append(filename)
    ha_list.sort()
    fe_list.sort()
    createVideoMp4(target_dir, ha_list, save_dir + 'RSM' + video_date.strftime('%Y-%m-%d') + '_HA.mp4')
    createVideoMp4(target_dir, fe_list, save_dir + 'RSM' + video_date.strftime('%Y-%m-%d') + '_FE.mp4')


if __name__ == '__main__':
    pass
