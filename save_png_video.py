"""
此py文件用来封装生成png预览图与预览视频的函数
主要用途: 读取当前文件夹内的所有fits文件, 将其分为两组(FE窗口图像与HA窗口图片)
随后生成对应的png图像，并且生成两个AVI预览视频

@origin: nju_rsh
@editor: seu_wxy
"""
import os
import subprocess
from datetime import datetime

import cv2
import matplotlib.pyplot as plt
import requests
from astropy.io import fits

import video_config


def monographNJU(target_path: str, color_map, image_dpi=100):
    """
    读取传入的目标路径, 列出目标文件夹内所有的ha fe文件
    再根据ha/fe fits文件列表在target_dir的文件夹内生成相应的png文件
    生成的png文件与fits文件的文件名相同, 仅格式不同
    @param target_path 传入的目标路径
    @param color_map 传入的色谱类型
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

    for rsm_path in ha_fits_list:
        # 从路径内获取文件名
        current_filename = rsm_path.split('/')[-1]
        current_file_datetime = current_filename[3:7] + '-' + current_filename[7:9] + '-' + current_filename[9:11] + \
                                ' UT ' + current_filename[12:14] + ':' + current_filename[
                                                                         14:16] + ':' + current_filename[16:18]

        # 查看png是否存在
        if os.path.exists(target_path + current_filename.split('.')[0] + ".png"):
            # 如果存在 直接跳过
            continue

        hdu = fits.open(rsm_path)
        # 读取bin mode
        bin_mode = int(hdu[1].header['BIN'])

        hacore = hdu[1].data[68 * 2 // bin_mode, :, :]

        # 一些参数 (与bin相关)
        ang_res = 0.5218 * bin_mode
        tick_pixel_window_size = 1040 * 2 // bin_mode
        tick_pixel = [tick_pixel_window_size - 1000 / ang_res, tick_pixel_window_size - 500 / ang_res,
                      tick_pixel_window_size, tick_pixel_window_size + 500 / ang_res,
                      tick_pixel_window_size + 1000 / ang_res]
        tick_arcsec = [-1000, -500, 0, 500, 1000]
        # 从中间剪裁的窗口大小
        window_scale = 1040 * 2 // bin_mode

        cx = int(hdu[1].header['CRPIX1'])
        cy = int(hdu[1].header['CRPIX2'])

        plt.figure(figsize=(20, 20))
        plt.imshow(hacore[cy - window_scale:cy + window_scale, cx - window_scale:cx + window_scale], origin='lower',
                   vmin=0, vmax=3 * hacore.mean(),
                   cmap=color_map)
        plt.text(1950, 50, current_file_datetime, horizontalalignment='right', verticalalignment='bottom',
                 color='white', size=28)
        plt.xticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.yticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.xlabel('Solar X (arcsec)', fontsize=30)
        plt.ylabel('Solar Y (arcsec)', fontsize=30)
        plt.savefig(target_path + current_filename.split('.')[0] + ".png", dpi=image_dpi, pad_inches=0.0,
                    bbox_inches='tight')
        plt.close()
        # print('/data/home/MikeRao/data/Time_array/ha' + rsm_path.split('/')[-1][3:18] + 'full_sun.png')

    for rsm_path in fe_fits_list:
        # 从路径内获取文件名
        current_filename = rsm_path.split('/')[-1]
        current_file_datetime = current_filename[3:7] + '-' + current_filename[7:9] + '-' + current_filename[9:11] + \
                                ' UT ' + current_filename[12:14] + ':' + current_filename[
                                                                         14:16] + ':' + current_filename[16:18]

        # 查看png是否存在
        if os.path.exists(target_path + current_filename.split('.')[0] + ".png"):
            # 如果存在 直接跳过
            continue

        hdu = fits.open(rsm_path)
        # 读取bin mode
        bin_mode = int(hdu[1].header['BIN'])
        hacore = hdu[1].data[10 * 2 // bin_mode, :, :]

        ang_res = 0.5218 * bin_mode
        tick_pixel_window_size = 1040 * 2 // bin_mode
        tick_pixel = [tick_pixel_window_size - 1000 / ang_res, tick_pixel_window_size - 500 / ang_res,
                      tick_pixel_window_size, tick_pixel_window_size + 500 / ang_res,
                      tick_pixel_window_size + 1000 / ang_res]
        tick_arcsec = [-1000, -500, 0, 500, 1000]

        cx = int(hdu[1].header['CRPIX1'])
        cy = int(hdu[1].header['CRPIX2'])

        plt.figure(figsize=(20, 20))
        plt.imshow(hacore[cy - window_scale:cy + window_scale, cx - window_scale:cx + window_scale], origin='lower',
                   vmin=0, vmax=3 * hacore.mean(),
                   cmap=color_map)
        plt.text(1950, 50, current_file_datetime, horizontalalignment='right', verticalalignment='bottom',
                 color='white', size=28)
        plt.xticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.yticks(ticks=tick_pixel, labels=tick_arcsec, fontsize=22)
        plt.xlabel('Solar X (arcsec)', fontsize=30)
        plt.ylabel('Solar Y (arcsec)', fontsize=30)
        plt.savefig(target_path + current_filename.split('.')[0] + ".png", dpi=image_dpi, pad_inches=0.0,
                    bbox_inches='tight')
        plt.close()
        # print('/data/home/MikeRao/data/Time_array/fe' + rsm_path.split('/')[-1][3:18] + 'full_sun.png')


def createVideoMp4(target_dir: str, img_array: list, filename: str):
    """
    传入png文件列表, 生成mpeg-4编码的avi格式视频
    @param target_dir 生成视频的存放文件夹
    @param img_array 文件列表, 内容为绝对路径
    @param filename 生成视频的名称
    """
    size = video_config.video_frame_size
    imgs = []
    for indexf in range(len(img_array)):  # 这个循环是为了读取所有要用的图片文件
        if indexf % video_config.target_png_interval == 0:
            img1 = cv2.imread(target_dir + img_array[indexf], 1)
            if img1 is None:
                print(img_array[indexf] + " is error!")
                continue
            resized = cv2.resize(img1, size, interpolation=cv2.INTER_AREA)
            imgs.append(resized)

    videowrite = cv2.VideoWriter(filename, cv2.VideoWriter_fourcc(*'XVID'), video_config.video_frame_per_sec, size)

    for i in range(len(imgs)):  # 把读取的图片文件写进去
        videowrite.write(imgs[i])

    videowrite.release()


def createVideoNJU(target_dir: str, save_dir: str):
    """
    读取png列表, 创建视频, 主要调用createVideoMp4函数执行生成视频
    @param target_dir 目标png所在的文件路径
    @param save_dir 视频存储路径
    """
    if not os.path.exists(save_dir):
        # 如果文件夹不存在，则创建文件夹
        try:
            os.makedirs(save_dir)
            print("文件夹创建成功")
        except OSError as e:
            print(f"文件夹创建失败: {e.filename} - {e.strerror}")
    else:
        print("文件夹已经存在")

    arr = os.listdir(target_dir)
    ha_list = []
    fe_list = []
    for i in range(len(arr)):
        filename = arr[i]
        if filename.split('.')[-1] == 'png':
            if filename.split('_')[-1].split('.')[0] == 'HA':
                ha_list.append(filename)
            if filename.split('_')[-1].split('.')[0] == 'FE':
                fe_list.append(filename)
    ha_list.sort()
    fe_list.sort()
    if target_dir.split('/')[-1] == "":
        current_date_str = target_dir.split('/')[-4] + target_dir.split('/')[-3].zfill(2) + target_dir.split('/')[
            -2].zfill(2)
    else:
        current_date_str = target_dir.split('/')[-3] + target_dir.split('/')[-2].zfill(2) + target_dir.split('/')[
            -1].zfill(2)
    HA_filename = 'RSM' + current_date_str + '_HA.avi'
    FE_filename = 'RSM' + current_date_str + '_FE.avi'
    createVideoMp4(target_dir, ha_list, save_dir + HA_filename)
    createVideoMp4(target_dir, fe_list, save_dir + FE_filename)
    return HA_filename, FE_filename, current_date_str;


# def test_monoGraphNJU():
#     suntools.log("重新输出png中")
#     png_target_dir = config.save_dir_path
#     color_map = suntools.get_color_map(config.color_camp_name)
#     monographNJU(png_target_dir, color_map)

def scpVideoToRemote(remote_ssh_url: str, source_file: str):
    """
    通过scp传送单个目标视频文件至远端服务器
    @param remote_ssh_url: 远端服务器sshURL 具有格式 username@ipaddr:PATH/filename
    @param source_file: 目标文件路径
    """
    current_filename = source_file.split("/")[-1]

    subprocess.run(["scp", source_file, remote_ssh_url])


def format_date(input_date):
    # 将输入的时间字符串解析为datetime对象
    date_obj = datetime.strptime(input_date, "%Y%m%d")
    # 格式化为所需的格式
    formatted_date = date_obj.strftime("/%Y-%m-%d/")
    return formatted_date


def convert_video_to_h264(input_file, output_file):
    # Use ffmpeg to convert video to H.264 format
    temp_output_file = output_file.replace('.avi', '.mp4');
    cmd = ['ffmpeg', '-y', '-i', input_file, '-c:v', 'libx264', '-crf', '23', '-c:a', 'copy', temp_output_file]
    subprocess.run(cmd)
    # Replace original file with temp file
    print(input_file)
    os.remove(input_file);

def upload_file(file_path):
    try:
        upload_url = 'http://114.212.184.2:8089/console-service/v1/nd_video/dailyVideo/upload'
        # 打开文件
        with open(file_path, 'rb') as file:
            # 构建请求数据
            files = {'file': file}
            # 发送POST请求
            response = requests.post(upload_url, files=files)
            # 检查响应状态码
            if response.status_code == 200:
                print("文件上传成功")
            else:
                print(f"文件上传失败: {response.status_code}")
    except Exception as e:
        print(f"文件上传失败: {e}")

if __name__ == '__main__':
    # test_monoGraphNJU()
    png_target_dir = video_config.path_to_target_png
    video_save_dir = video_config.path_to_video_save

    print(
        f"视频生成参数如下: \n帧大小:{video_config.video_frame_size.__str__()} \nFPS: {video_config.video_frame_per_sec.__str__()}  \n生成间隔: {video_config.target_png_interval.__str__()}")
    print(f"存放PNG的目标文件夹为:{png_target_dir}")
    Ha, Fe, Date = createVideoNJU(png_target_dir, video_save_dir)

    remote_host = video_config.remote_host
    destPath = video_config.remote_path + format_date(Date)
    convert_video_to_h264(video_save_dir + Ha, video_save_dir + Ha)
    convert_video_to_h264(video_save_dir + Fe, video_save_dir + Fe)
    user = video_config.user;
    # command = f"ssh {user}@{remote_host} 'mkdir -p {destPath}'"
    # # 执行远程命令
    # subprocess.run(command, shell=True)
    #
    # command = f"ssh {user}@{remote_host} 'mkdir -p {destPath}'"
    # # 执行远程命令
    # subprocess.run(command, shell=True)
    # scpVideoToRemote(
    #     "root@114.212.184.2:/home/aiocloud/oradiag_aiocloud/app/bg/ui-console/ui/ui-console/src/main/resources/static/video/dailyVideo" + format_date(
    #         Date) + Ha.replace('.avi', '.mp4'), video_save_dir + Ha.replace('.avi', '.mp4'))
    # scpVideoToRemote(
    #     "root@114.212.184.2:/home/aiocloud/oradiag_aiocloud/app/bg/ui-console/ui/ui-console/src/main/resources/static/video/dailyVideo" + format_date(
    #         Date) + Fe.replace('.avi', '.mp4'), video_save_dir + Fe.replace('.avi', '.mp4'))
    upload_file( video_save_dir + Ha.replace('.avi', '.mp4') )
    upload_file( video_save_dir + Fe.replace('.avi', '.mp4') )
    print(f"生成的视频已保存在: {video_save_dir}内")
    print(f"生成的HA预览视频文件名为: {Ha.replace('.avi', '.mp4')}")
    print(f"生成的FE预览视频文件名为: {Fe.replace('.avi', '.mp4')}")
