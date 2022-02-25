import multiprocessing as mp
import os
import suntools
import time
import matplotlib.pyplot as plt
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import scipy.signal as signal
import config


# 读入配置文件 引入参数
read_dir = config.data_dir_path
out_dir = config.save_dir_path
multiprocess_count = 1
if config.multiprocess_count != 'default':
    multiprocess_count = config.multiprocess_count
else:
    multiprocess_count = mp.cpu_count() - 4

# 读取暗场文件
# 此处的数据均未做共享处理，因为共享数据量并不是很大，在LINUX环境下使用multiprocess并fork()将直接复制此些全局变量
# ***在Windows环境下并不适用!!!!***
temp_img = get_pkg_data_filename(config.dark_fits_name)
dark_img = np.array(fits.getdata(temp_img), dtype=float)
# 平场需要进行谱线弯曲矫正并计算才可以得到真正的平场
temp_img = get_pkg_data_filename(config.flat_fits_name)
flat_img = np.array(fits.getdata(temp_img), dtype=float)
flat_img, temp1, temp2 = suntools.curve_correction(flat_img, config.curve_cor_x0, config.curve_cor_C)
flat_img = suntools.getFlat(flat_img)
# 读取标准太阳光谱数据
sun_std = suntools.get_Sunstd(config.sun_std_name)


# 读取数据文件夹所有文件
def read_fits_directory():
    arr = os.listdir(read_dir)
    return arr


# 全局进度控制
file_count = mp.Value('i', len(read_fits_directory()))
remaining_count = mp.Value('i', 0)


# 定义target task
# 传入一个文件名，读取此文件名对应的fits文件并对其做曲线矫正
def target_task(filename):
    filePath = read_dir + "/" + filename
    file_data = get_pkg_data_filename(filePath)
    image_data = np.array(fits.getdata(file_data), dtype=float)
    # 去暗场
    image_data = image_data - dark_img
    # 谱线弯曲矫正
    image_data, HofH, HofFe = suntools.curve_correction(image_data, config.curve_cor_x0, config.curve_cor_C)
    # 去平场
    image_data = suntools.DivFlat(image_data, flat_img, HofH, HofFe)
    # 红蓝移矫正
    image_data = suntools.RB_repair(image_data, sun_std, HofH, HofFe)
    # 转为整型
    image_data = np.array(image_data, dtype=np.int16)
    # 滤波
    image_data = signal.medfilt(image_data, kernel_size=config.filter_kernel_size)
    # 存储fits
    primaryHDU = fits.PrimaryHDU(image_data)
    greyHDU = fits.HDUList([primaryHDU])
    greyHDU.writeto(out_dir + filename)
    # 进度输出
    remaining_count.value += 1
    print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value))


def main():
    # 测试消耗时间 时间起点
    time_start = time.time()
    # 获得文件夹列表 读取相关参数
    data_file_lst = read_fits_directory()
    # 并行处理
    pool = mp.Pool(processes=multiprocess_count)
    print('多核并行数:' + str(multiprocess_count))
    pool.map(target_task, data_file_lst)
    time_end = time.time()
    print('并行进度已完成，所花费时间为：', (time_end - time_start) / 60, 'min(分钟)')
    # 汇总处理结果
    sum_file_path = config.sum_dir_path
    now = 0
    N = len(os.listdir(out_dir))
    data = np.zeros((N, 4608))
    for filename in os.listdir(out_dir):
        image_file = get_pkg_data_filename(out_dir + "/" + filename)
        image_data = fits.getdata(image_file)
        count = int(filename[-8:-4]) - 1
        data[count, :] = image_data[0, :]
        now += 1
        if now >= N:
            break
    plt.imsave(sum_file_path + 'sum.jpg', data)


if __name__ == "__main__":
    main()
