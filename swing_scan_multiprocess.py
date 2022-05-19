import multiprocessing as mp
import datetime
import os
import header
import suntools
import time
import numpy as np
from astropy.io import fits
import urllib.error as uEr
import config
import matplotlib.pyplot as plt
import ctypes as c
import sys

# 读入配置文件 引入参数
GLOBAL_BINNING = config.bin_count  # binning 数值
READ_DIR = config.data_dir_path  # 读文件的文件夹
OUT_DIR = config.save_dir_path  # 输出文件夹
SUM_DIR = config.sum_dir_path  # 汇总结果文件夹
DARK_FITS_FILE = config.dark_fits_name  # 暗场文件路径
HA_ABSORPTION_FILE = config.HA_absorption_path  # HA吸收系数文件路径
FE_ABSORPTION_FILE = config.FE_absorption_path  # FE吸收系数文件路径
COLOR_CAMP_FILE = config.color_camp_name  # 色彩盘文件路径
HEADER_FILE = config.header_file  # 读取的标准头部文件路径
HA_START = config.HA_start  # HA窗口起始波长
FE_START = config.FE_start  # FE窗口起始波长
HA_LINE_CORE = config.HA_lineCore  # HA线心波长
FE_LINE_CORE = config.FE_lineCore  # FE线心波长
FILTER_KERNEL_SIZE = config.filter_kernel_size  # 滤波窗口
FLAT_FITS_FILE = ''  # 平场文件路径(与bin相关)
SUN_ROW_COUNT = 0  # 太阳序列数(与bin相关)
STANDARD_FILE_INDEX = 0  # 标准文件在序列中的位置(与bin相关)
CURVE_X0 = 0  # 谱线矫正参数x0(与bin相关)
CURVE_C = 0  # 谱线矫正参数c(与bin相关)
WAVE_RESOLUTION = 0  # 波长分辨率(与bin相关)
SUM_ROW_INDEX_HA = 0  # 合并HA日像所选的行数(与bin相关)
SUM_ROW_INDEX_FE = 0  # 合并FE日像所选的行数(与bin相关)
SCAN_TIME_OFFSET = config.scan_time_offset  # 时间偏差
SIT_STARE_MODE = config.sit_stare_mode  # sit stare模式
PIXEL_RESOLUTION = config.pixel_resolution  # 像素分辨率
PIXEL_ZERO_COUNT = config.pixel_to_zero_count  # 置零区间
if GLOBAL_BINNING == 1:
    FLAT_FITS_FILE = config.flat_fits_name_bin_1
    SUN_ROW_COUNT = config.sun_row_count_bin_1
    STANDARD_FILE_INDEX = config.standard_offset_index_bin_1
    CURVE_X0 = config.curve_cor_x0_bin_1
    CURVE_C = config.curve_cor_C_bin_1
    WAVE_RESOLUTION = config.wavelength_resolution_bin_1
    SUM_ROW_INDEX_HA = config.sum_row_index_HA_bin_1
    SUM_ROW_INDEX_FE = config.sum_row_index_FE_bin_1
if GLOBAL_BINNING == 2:
    FLAT_FITS_FILE = config.flat_fits_name_bin_2
    SUN_ROW_COUNT = config.sun_row_count_bin_2
    STANDARD_FILE_INDEX = config.standard_offset_index_bin_2
    CURVE_X0 = config.curve_cor_x0_bin_2
    CURVE_C = config.curve_cor_C_bin_2
    WAVE_RESOLUTION = config.wavelength_resolution_bin_2
    SUM_ROW_INDEX_HA = config.sum_row_index_HA_bin_2
    SUM_ROW_INDEX_FE = config.sum_row_index_FE_bin_2

multiprocess_count = 1
if config.multiprocess_count != 'default':
    multiprocess_count = config.multiprocess_count
else:
    multiprocess_count = mp.cpu_count() - 4
print('多核并行数:' + str(multiprocess_count))


# 读取数据文件夹所有文件
def read_fits_directory():
    arr = []
    arr = os.listdir(READ_DIR)
    if len(arr) == 0:
        raise OSError
    return arr


# 此处的数据均未做共享处理，因为共享数据量并不是很大，在LINUX环境下使用multiprocess并fork()将直接复制此些全局变量
# 预读输入目录
try:
    data_file_lst = read_fits_directory()
except OSError:
    print('没有获得原始数据文件，请检查config中的读入数据目录')
    sys.exit("程序终止")
print('文件总数为: ' + str(len(data_file_lst)))
print('共包含:' + str(len(data_file_lst) / SUN_ROW_COUNT) + '个序列')

# 将读入的文件按照序列分成不同的组
# 此处坑比较大
# 首先需要按照文件前18位的时间来进行排序
# 相同时间的按照后面的序列号和帧号来排序
# 并且将其分为不同的组 不能按照序列号进行索引 可以直接将3-27位拉出来字符串排序
# 一个标准文件名 如下:
# RSM 2021   12     22     T 060105   -   0008-     0001       .fts
# 012 3456   78     90     1 234567   8   90123     4567       8901
#     [year] [mon]  [day] [T hhMMSS]      [index]   [frame]
global_multiprocess_list = []  # 存放序列dict的全局数组
# 对list内的文件名排序
data_file_lst.sort(key=lambda x: x.split('-')[0] + x.split('-')[1] + str(int(x.split('-')[2].split('.'))).zfill(6))
global_wave_line_strength_list = []
# 读取每个文件某一行的像素强度并记录在list内
for filename in data_file_lst:
    temp_img = fits.open(READ_DIR + filename)
    temp_data = np.array(temp_img[0].data, dtype=float)
    temp_mean = np.mean(temp_data[SUM_ROW_INDEX_HA, :])
    global_wave_line_strength_list.append(temp_mean)
for gwl in global_wave_line_strength_list:
    print(gwl)

