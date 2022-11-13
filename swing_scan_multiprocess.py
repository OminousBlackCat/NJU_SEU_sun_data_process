"""
本py文件对应安装后的run_scan_solar文件 使用source指令可以执行该主程序
本程序对应摆扫扫描模式下对图像数据的预处理操作
与multiprocess.py的主要区别体现在对读入文件列表的处理 增添了对每一个扫描序列的判断
注意在使用前修改config.py内的参数!

@author: seu_wxy
"""


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
import createVideo
import sys

time_start = time.time()

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
CENTER_MEAN_THRESHOLD = config.center_mean_threshold  # (摆扫)计算序列中心的阈值
CENTER_MEAN_INDEX = config.center_mean_index  # (摆扫)使用的线心位置
CENTER_MEAN_COUNT = config.center_mean_count  # (摆扫)使用的线心数量
REVERSAL_MODE = config.reversal_mode  # (摆扫)翻转模式
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

# 检查输出文件夹是否存在 不存在则创建
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

multiprocess_count = 1
if config.multiprocess_count != 'default':
    multiprocess_count = config.multiprocess_count
else:
    multiprocess_count = mp.cpu_count() - 4
suntools.log('多核并行数:' + str(multiprocess_count))


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
    suntools.log('没有获得原始数据文件，请检查config中的读入数据目录')
    sys.exit("程序终止")
suntools.log("程序目标文件夹为：" + READ_DIR)
suntools.log('文件总数为: ' + str(len(data_file_lst)))
suntools.log('当前运行处在 摆扫序列处理模式')

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
real_data_file_lst = []  # 存放没有非法文件名的文件列表
# 此处对文件合法性进行判断
suntools.log('判断文件名合法性中...')
for filename in data_file_lst:
    try:
        if filename.split('-')[0][0:3] != 'RSM':
            raise ValueError
        if len(filename.split('-')[1]) > 4:
            raise ValueError
        sort_test = filename.split('-')[0] + filename.split('-')[1] + str(int(filename.split('-')[2].split('.')[0])).zfill(8)
    except BaseException as e:
        suntools.log('<文件:' + filename + '>非法, 已剔除(并未删除硬盘上的文件)')
    real_data_file_lst.append(filename)

data_file_lst = real_data_file_lst  # 将剔除后的列表赋过去
# 对list内的文件名排序
suntools.log('对文件进行排序中...')
# 文件名排序关键字: yyyymmddThhMMSS0000-00000001(帧数补零8位) 可以直接按照时间与帧数排为正序

data_file_lst.sort(key=lambda x: x.split('-')[0] + x.split('-')[1] + str(int(x.split('-')[2].split('.')[0])).zfill(8))
global_wave_line_strength_list = []
# 读取每个文件某一行的像素强度并记录在list内
suntools.log('读取图像像素中...')
have_read_count = 0
if_read_first_print = True
for filename in data_file_lst:
    # if if_read_first_print:
    #     print('当前进度:' + str(have_read_count) + '/' + str(len(data_file_lst)), end='')
    #     sys.stdout.flush()
    #     if_read_first_print = False
    # else:
    #     print('\b' * (9 + len(str(have_read_count)) + 1 + len(str(len(data_file_lst)))), end='')
    #     print('当前进度:' + str(have_read_count) + '/' + str(len(data_file_lst)), end='')
    #     sys.stdout.flush()
    temp_img = fits.open(READ_DIR + filename)
    temp_data = np.array(temp_img[0].data, dtype=float)
    temp_mean = np.mean(temp_data[CENTER_MEAN_INDEX: CENTER_MEAN_INDEX + CENTER_MEAN_COUNT, :])
    global_wave_line_strength_list.append(temp_mean)
    have_read_count += 1
last_wave_line_strength = 0
significant_point_list = []
symmetry_axis_list = []
# 以150为分界线寻找对称轴 记录这些关键点
# 标记0为上升点 标记1为下降点
for i in range(len(global_wave_line_strength_list)):
    if global_wave_line_strength_list[i] >= CENTER_MEAN_THRESHOLD >= last_wave_line_strength:
        significant_point_list.append([i, 0])
    if global_wave_line_strength_list[i] <= CENTER_MEAN_THRESHOLD <= last_wave_line_strength:
        significant_point_list.append([i, 1])
    last_wave_line_strength = global_wave_line_strength_list[i]
for point in significant_point_list:
    if point[1] == 0:
        if significant_point_list.index(point) + 1 >= len(significant_point_list):
            break
        temp_frame_index = int(data_file_lst[point[0]].split('-')[-1].split('.')[0])
        next_frame_index = int(data_file_lst[significant_point_list[significant_point_list.index(point) + 1][0]]
                               .split('-')[-1].split('.')[0])
        # 如果下一个关键点的key为1 且 下一个关键点对应的帧数序号比当前关键点的序号高 则可以获得一个对称轴
        if significant_point_list[significant_point_list.index(point) + 1][1] == 1 and\
                next_frame_index > temp_frame_index:
            symmetry_axis_list.append(int((point[0] +
                                           significant_point_list[significant_point_list.index(point) + 1][0]) / 2))
suntools.log('此文件夹共找到' + str(len(symmetry_axis_list)) + '个序列')
last_symmetry_axis_frame_index = 0
current_scan_index = 0
current_track_index = 0
for axis in symmetry_axis_list:
    suntools.log('*************************************************************************************')
    suntools.log('文件名:' + data_file_lst[axis] + '/ 平均强度为:' + str(global_wave_line_strength_list[axis]))
    # 当前对称轴的帧数
    current_axis_frame_index = int(data_file_lst[axis].split('-')[-1].split('.')[0])
    # 如果当前的帧数小于上次的帧数 说明是新的轨道
    if current_axis_frame_index < last_symmetry_axis_frame_index:
        current_track_index += 1  # 轨道数加一
        current_scan_index = 0  # 将扫描序列序号置0
    temp_start_file_index = axis - int(SUN_ROW_COUNT / 2)
    temp_last_file_index = axis + int(SUN_ROW_COUNT / 2)
    if temp_start_file_index < 0:
        temp_start_file_index = 0
    if temp_last_file_index > len(data_file_lst) - 1:
        temp_last_file_index = len(data_file_lst) - 1
    global_multiprocess_list.append({
        'key_index': symmetry_axis_list.index(axis),  # 唯一序号
        'track_index': current_track_index,  # 轨道序号
        'scan_index': str(current_scan_index).zfill(4),  # 扫描序列序号
        'file_list': data_file_lst[temp_start_file_index: temp_last_file_index + 1],  # 文件名列表
        'file_count': temp_last_file_index - temp_start_file_index + 1,  # 包含的文件数
        'standard_filename': data_file_lst[axis],  # 标准日心文件名
        'first_filename': data_file_lst[temp_start_file_index],  # 序列开始文件
        'last_filename': data_file_lst[temp_last_file_index],  # 序列结束文件
        'flat_data': None,  # 此序列的校正后平场数据
        'abortion_data': None,  # 此序列的校正后红蓝移数据
        'header': fits.header.Header(),  # 此序列的头部, 构造了一个新的header
        'start_time': datetime.datetime.now()
    })
    suntools.log('此序列处于第:' + str(current_track_index) + '轨')
    suntools.log('对应序列序号为:' + str(current_scan_index).zfill(4))
    suntools.log('序列文件总数为:' + str(temp_last_file_index - temp_start_file_index + 1))
    suntools.log('起始文件名:' + data_file_lst[temp_start_file_index])
    suntools.log('结束文件名:' + data_file_lst[temp_last_file_index])
    suntools.log('*************************************************************************************')
    current_scan_index += 1
    last_symmetry_axis_frame_index = current_axis_frame_index

# 读取头部参数文件
# 为每个序列都创建头
global_header_list = header.read_header_from_txt(HEADER_FILE)
# 将读入的list构造成为header
for h in global_header_list:
    for temp_dict in global_multiprocess_list:
        temp_dict['header'].set(h['key'], value=h['value'], comment=h['comment'])

# 读取暗场文件
temp_img = None
dark_img = None
try:
    suntools.log("正在读取原始暗场文件")
    temp_img = fits.open(DARK_FITS_FILE)
except uEr.URLError:
    suntools.log("Error: 暗场文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    suntools.log("Error: 暗场文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")
if temp_img is not None:
    dark_img = np.array(temp_img[0].data, dtype=float)
    dark_img = suntools.change(dark_img)
temp_img.close()

# 平场需要以日心图片作为基准进行平移矫正 再进行谱线弯曲矫正
flat_img = None
standard_HA_width, standard_FE_width = None, None
try:
    suntools.log("正在读取原始平场文件")
    temp_img = fits.open(FLAT_FITS_FILE)
except uEr.URLError:
    suntools.log("Error: 原始平场文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    suntools.log("Error: 原始平场文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")
if temp_img is not None:
    flat_img = np.array(temp_img[0].data, dtype=float)
    flat_img, standard_HA_width, standard_FE_width = suntools.curve_correction(flat_img - dark_img, CURVE_X0,
                                                                               CURVE_C)
    # 对平场进行归一化
    flat_img = suntools.FlatNormalization(flat_img)
temp_img.close()

# 读取经过日心的图片 作为基准
# 读取标准太阳光谱数据
# sun_std = suntools.get_Sunstd(config.sun_std_name)
global_absorption = suntools.get_Absorstd(
    HA_ABSORPTION_FILE, FE_ABSORPTION_FILE, standard_HA_width, standard_FE_width)
sample_from_standard = None
try:
    for temp_dict in global_multiprocess_list:
        # 对每个序列进行校正
        suntools.log('校正扫描第' + str(temp_dict['track_index']) + '轨, 序列' + temp_dict['scan_index'] + '中...使用标准校正文件为:' + temp_dict['standard_filename'])
        suntools.log('此序列首文件为:' + temp_dict['first_filename'])
        suntools.log('此序列末文件为:' + temp_dict['last_filename'])
        suntools.log("校正平场中...")
        standard_name = None
        if temp_dict['standard_filename'] == '' and SIT_STARE_MODE:
            suntools.log('此序列不完整且不包含标准序列文件, 将使用最靠近中心的文件作为矫正基准')
            if int(temp_dict['first_filename'][24: 28]) < 100:
                standard_name = temp_dict['last_filename']
            if int(temp_dict['last_filename'][24: 28]) > 4000 / GLOBAL_BINNING:
                standard_name = temp_dict['last_filename']
        else:
            standard_name = temp_dict['standard_filename']
        temp_img = fits.open(READ_DIR + standard_name)
        standard_header = temp_img[0].header
        for item in header.copy_header_items:
            temp_dict['header'].set(item['key'], standard_header[item['key']])

        standard_img = np.array(temp_img[0].data, dtype=float)
        standard_img = suntools.moveImg(standard_img, -2)
        standard_img, standard_HA_width, standard_FE_width = suntools.curve_correction(standard_img - dark_img,
                                                                                       CURVE_X0,
                                                                                       CURVE_C)
        sample_from_standard = standard_img
        # 先平移矫正 减去暗场 再谱线弯曲矫正
        flatTemp = suntools.getFlatOffset(flat_img, standard_img)
        flatTemp = suntools.getFlat(flatTemp)
        suntools.log("获得标准太阳光谱数据中...")
        # 以标准文件作为基准 计算红蓝移吸收系数
        # 需要先对标注文件进行一系列操作 去暗场 去平场 再进行红蓝移修正
        standard_img = suntools.DivFlat(standard_img, flatTemp)
        # 获得标准吸收系数
        # abortion = suntools.RB_getdata(standard_img, sun_std, standard_HA_width, standard_FE_width)
        temp_dict['flat_data'] = flatTemp
        # temp_dict['abortion_data'] = abortion
        temp_img.close()
        suntools.log("序列:" + temp_dict['scan_index'] + "矫正完成")
        suntools.log('计算B0, INST_ROT中....')
        temp_B0, temp_INST_ROT = suntools.getB0P0(standard_header['Q0'], standard_header['Q1'], standard_header['Q2'],
                                                  standard_header['Q3'], standard_header['STR_TIME'])

        time_offset = datetime.timedelta(seconds=SCAN_TIME_OFFSET)
        first_name = temp_dict['first_filename']
        last_name = temp_dict['last_filename']
        start_temp_time = datetime.datetime(year=int(first_name[3: 7]), month=int(first_name[7: 9]),
                                            day=int(first_name[9: 11]), hour=int(first_name[12: 14]),
                                            minute=int(first_name[14:16]), second=int(first_name[16: 18]))
        start_temp_time = start_temp_time + time_offset
        temp_dict['start_time'] = start_temp_time
        end_temp_time = datetime.datetime(year=int(last_name[3: 7]), month=int(last_name[7: 9]),
                                          day=int(last_name[9: 11]), hour=int(last_name[12: 14]),
                                          minute=int(last_name[14:16]), second=int(last_name[16: 18]))
        end_temp_time = end_temp_time + time_offset
        # 对头进行赋值
        temp_dict['header'].set('BIN', GLOBAL_BINNING)
        temp_dict['header'].set('B0', temp_B0)
        temp_dict['header'].set('INST_ROT', temp_INST_ROT)
        temp_dict['header'].set('DATE_OBS', start_temp_time.strftime('%Y-%m-%dT%H:%M:%S'))
        temp_dict['header'].set('STR_TIME', start_temp_time.strftime('%Y-%m-%dT%H:%M:%S'))
        temp_dict['header'].set('END_TIME', end_temp_time.strftime('%Y-%m-%dT%H:%M:%S'))
        temp_dict['header'].set('FRM_NUM', '1~' + str(temp_dict['file_count']))

except uEr.URLError as error:
    suntools.log("Error: 标准日心校准文件未找到, 请检查config文件或存放目录")
    suntools.log(error)
    sys.exit("程序终止")
except OSError as error:
    suntools.log("Error: 标准日心校准文件读取发生错误, 请检查文件读取权限")
    suntools.log(error)
    sys.exit("程序终止")

# 读取输出色谱
color_map = suntools.get_color_map(COLOR_CAMP_FILE)


# 全局进度控制
file_count = mp.Value('i', len(read_fits_directory()))
remaining_count = mp.Value('i', 1)
if_first_print = mp.Value('b', True)

# 全局共享内存
# 需要创建一个四维数组
# 对于最后写入的文件 NAXIS顺序为: 狭缝宽度 文件序号(扫描序号)  波长深度
# 因此np array的shape应为(波长深度, 文件序号, 狭缝宽度)
# 每个单独文件对应的是 xy平面的一个二维数组
GLOBAL_ARRAY_X_COUNT = sample_from_standard.shape[0]
GLOBAL_ARRAY_Y_COUNT = SUN_ROW_COUNT
GLOBAL_ARRAY_Z_COUNT = sample_from_standard.shape[1]
suntools.log('SHAPE:' + str(GLOBAL_ARRAY_X_COUNT) + ',' + str(GLOBAL_ARRAY_Y_COUNT) + ',' + str(GLOBAL_ARRAY_Z_COUNT))
# 创建共享内存 大小为 x*y*z*sizeof(int16)
GLOBAL_SHARED_MEM = None
try:
    GLOBAL_SHARED_MEM = mp.Array(c.c_int16, GLOBAL_ARRAY_X_COUNT * GLOBAL_ARRAY_Y_COUNT * GLOBAL_ARRAY_Z_COUNT)
except BaseException as e:
    suntools.log("内存不足, 无法创建共享数组")
    sys.exit("程序结束")
GLOBAL_DICT_INDEX = 0  # 全局dict序号控制 为了在pool.map之后让每个子进程知道自己的dict序号


def INCREASE_DICT_INDEX():
    global GLOBAL_DICT_INDEX
    GLOBAL_DICT_INDEX += 1


# 定义target task
# 传入一个文件名，读取此文件名对应的fits文件并对其做曲线矫正
def target_task(filename):
    try:
        # 一个标准文件名 如下:
        # RSM 2021   12     22T060105   -   0008-     0001       .fts
        # 012 3456   78     901234567   8   90123     4567       8901
        #     [year] [mon]  [day_seq]       [index]   [position]
        # 计算此文件在序列中的相对位置 给后续放入全局数组做准备
        fileRelativePosition = int(filename.split('-')[2].split('.')[0]) - int(
            global_multiprocess_list[GLOBAL_DICT_INDEX]['first_filename'].split('-')[2].split('.')[0])
        filePath = READ_DIR + filename
        file_data = fits.open(filePath)
        image_data = np.array(file_data[0].data, dtype=float)
        # 对fe窗口进行平移
        image_data = suntools.moveImg(image_data, -2)
        # 去暗场
        image_data = image_data - dark_img
        # 谱线弯曲矫正
        image_data, HofH, HofFe = suntools.curve_correction(image_data, CURVE_X0, CURVE_C)
        # 搜索list
        currentFlat = global_multiprocess_list[GLOBAL_DICT_INDEX]['flat_data']
        currentScanIndex = int(global_multiprocess_list[GLOBAL_DICT_INDEX]['scan_index'])
        if currentFlat is None:
            suntools.log("文件：" + filename + "未找到平场数据, 请检查文件夹")
            return
        # 去平场
        image_data = suntools.DivFlat(image_data, currentFlat)
        # 红蓝移矫正
        image_data = suntools.RB_repair(image_data, global_absorption)
        # 滤波
        image_data = suntools.MedSmooth(image_data, HofH, HofFe, winSize=FILTER_KERNEL_SIZE)
        # 转为整型, 并将每行的最后部分置零
        image_data = np.array(image_data, dtype=np.int16)
        image_data[:, image_data.shape[1] - PIXEL_ZERO_COUNT:] = 0
        global_shared_array = np.frombuffer(GLOBAL_SHARED_MEM.get_obj(), dtype=np.int16)
        global_shared_array = global_shared_array.reshape(GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT,
                                                          GLOBAL_ARRAY_Z_COUNT)
        if REVERSAL_MODE == 'odd' and currentScanIndex % 2 == 1:
            global_shared_array[:, SUN_ROW_COUNT - 1 - fileRelativePosition, :] = image_data
        elif REVERSAL_MODE == 'even' and currentScanIndex % 2 == 0:
            global_shared_array[:, SUN_ROW_COUNT - 1 - fileRelativePosition, :] = image_data
        else:
            global_shared_array[:, fileRelativePosition, :] = image_data
        # 进度输出
        remaining_count.value += 1
        file_data.close()
        # if if_first_print.value:
        #     print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value), end='')
        #     sys.stdout.flush()
        #     if_first_print.value = False
        # else:
        #     print('\b' * (9 + len(str(remaining_count.value)) + 1 + len(str(file_count.value))), end='')
        #     print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value), end='')
        #     sys.stdout.flush()
    except BaseException as e:
        suntools.log(e)
        suntools.log('文件:' + filename + '处理失败, 请检查此文件')


def main():
    # 测试消耗时间 时间起点
    # 获得文件夹列表 读取相关参数
    # 并行处理
    suntools.log('开启多核并行处理...')
    for temp_dict in global_multiprocess_list:
        suntools.log('正在处理第' + str(temp_dict['track_index']) + '轨,  扫描序列:' + temp_dict['scan_index'] + '...')
        file_count.value = temp_dict['file_count']
        remaining_count.value = 0
        pool = mp.Pool(processes=multiprocess_count)
        pool.map(target_task, temp_dict['file_list'])
        pool.close()
        pool.join()
        suntools.log('扫描序列' + temp_dict['scan_index'] + '预处理完成...')
        suntools.log('生成完整日像中...')
        try:
            sum_data_HA = np.zeros((SUN_ROW_COUNT, sample_from_standard.shape[1]))
            sum_data_FE = np.zeros((SUN_ROW_COUNT, sample_from_standard.shape[1]))
            global_shared_array = np.frombuffer(GLOBAL_SHARED_MEM.get_obj(), dtype=np.int16)
            global_shared_array = global_shared_array.reshape(GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT,
                                                              GLOBAL_ARRAY_Z_COUNT)
            # 将小于0的值全部赋为0
            global_shared_array[global_shared_array < 0] = 0
            suntools.log("SHAPE为：" + str(global_shared_array.shape))
            # 输出太阳像
            for i in range(global_shared_array.shape[1]):
                sum_data_HA[i] = global_shared_array[SUM_ROW_INDEX_HA, i, :].reshape(sample_from_standard.shape[1])
                sum_data_FE[i] = global_shared_array[standard_HA_width + SUM_ROW_INDEX_FE, i, :].reshape(
                    sample_from_standard.shape[1])
            suntools.log('计算CCD太阳像半径中...')
            R_y, R_x, radius = suntools.getCircle(sum_data_FE)
            OBS_Radius = radius * PIXEL_RESOLUTION * GLOBAL_BINNING
            suntools.log('波长定标中...')
            wavelength_calibrate_input = global_shared_array[:, int(R_y) - 50: int(R_y) + 49, int(R_x) - 50: int(R_x) + 49]
            cdel_t3, crval_l3_ha, crval_l3_fe = suntools.cal_center_mean(wavelength_calibrate_input)
            temp_dict['header'].set('CRPIX1', R_x)
            temp_dict['header'].set('CRPIX2', R_y)
            temp_dict['header'].set('R_SUN', radius)
            temp_dict['header'].set('RSUN_OBS', OBS_Radius)
            temp_dict['header'].set('CDELT1', PIXEL_RESOLUTION * GLOBAL_BINNING)
            temp_dict['header'].set('CDELT2', PIXEL_RESOLUTION * GLOBAL_BINNING)
            temp_dict['header'].set('CDELT3', cdel_t3)
            temp_dict['header'].set('CRVAL3', crval_l3_ha)
            # 下采样 1/4
            suntools.log('下采样中...')
            sum_data_HA_save = suntools.down_sample(sum_data_HA)
            sum_data_FE_save = suntools.down_sample(sum_data_FE)
            sum_data_HA_save = suntools.add_time(sum_data_HA_save, temp_dict['start_time'].strftime('%Y-%m-%d '
                                                                                                    '%H:%M:%S UT'))
            sum_data_FE_save = suntools.add_time(sum_data_FE_save, temp_dict['start_time'].strftime('%Y-%m-%d '
                                                                                                    '%H:%M:%S UT'))
            if config.save_img_form == 'default':
                # 使用读取的色谱进行输出 imsave函数将自动对data进行归一化
                suntools.log('输出序号为' + temp_dict['scan_index'] + '的png...')
                sum_mean_ha = np.mean(sum_data_HA)
                sum_mean_fe = np.mean(sum_data_FE)
                plt.imsave(SUM_DIR + 'RSM' + temp_dict['start_time'].strftime('%Y%m%dT%H%M%S')
                           + '_' + temp_dict['scan_index'] + '_HA' + ".png",
                           sum_data_HA_save, cmap=color_map, vmin=0, vmax=3 * sum_mean_ha)
                plt.imsave(SUM_DIR + 'RSM' + temp_dict['start_time'].strftime('%Y%m%dT%H%M%S')
                           + '_' + temp_dict['scan_index'] + '_FE' + ".png",
                           sum_data_FE_save, cmap=color_map, vmin=0, vmax=3 * sum_mean_fe)
            if config.save_img_form == 'fts':
                # 不对data进行任何操作 直接输出为fts文件
                suntools.log('输出序号为' + temp_dict['scan_index'] + '的fits...')
                primaryHDU = fits.PrimaryHDU(sum_data_HA)
                greyHDU = fits.HDUList([primaryHDU])
                greyHDU.writeto(SUM_DIR + 'SUM' + temp_dict['start_time'].strftime('%Y%m%dT%H%M%S')
                                + '_' + temp_dict['scan_index'] + '_HA' + '.fts', overwrite=True)
                greyHDU.close()
                primaryHDU = fits.PrimaryHDU(sum_data_FE)
                greyHDU = fits.HDUList([primaryHDU])
                greyHDU.writeto(SUM_DIR + 'SUM' + temp_dict['start_time'].strftime('%Y%m%dT%H%M%S')
                                + '_' + temp_dict['scan_index'] + '_FE' + '.fts', overwrite=True)
                greyHDU.close()
            suntools.log('生成HA文件中...')
            temp_dict['header'].set('SPECLINE', 'HA')
            temp_dict['header'].set('WAVE_LEN', HA_LINE_CORE)
            temp_dict['header'].set('PRODATE', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
            primaryHDU = fits.CompImageHDU(global_shared_array[0: standard_HA_width, :, :]
                                           .reshape((standard_HA_width, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT))
                                           , header=temp_dict['header'], compression_type='RICE_1')
            primaryHDU.header.set('NAXIS', comment='Number of data axes')
            primaryHDU.header.set('NAXIS1', comment='Length of data axis 1 (slit dimension)')
            primaryHDU.header.set('NAXIS2', comment='Length of data axis 2 (scanning steps)')
            primaryHDU.header.set('NAXIS3', comment='Length of data axis 3 (wavelength dimension)')
            primaryHDU.header.add_comment('Spectral curvature corrected')
            primaryHDU.header.add_comment('Dark subtracted')
            primaryHDU.header.add_comment('Flat-field corrected')
            primaryHDU.header.add_comment('Processed by RSM_prep')
            suntools.log(repr(primaryHDU.header))
            primaryHDU.writeto(OUT_DIR + 'RSM' + temp_dict['start_time'].strftime('%Y%m%dT%H%M%S') + '_' +
                               temp_dict['scan_index'] + '_HA.fits', overwrite=True)
            suntools.log('生成FE文件中...')
            # 修改header内的SPECLINE与WAVELNTH
            temp_dict['header'].set('SPECLINE', 'FEI')
            temp_dict['header'].set('WAVE_LEN', FE_LINE_CORE)
            temp_dict['header'].set('PRODATE', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
            temp_dict['header'].set('CRVAL3', crval_l3_fe)
            primaryHDU = fits.CompImageHDU(global_shared_array[standard_HA_width:, :, :]
                                           .reshape((standard_FE_width, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT))
                                           , header=temp_dict['header'], compression_type='RICE_1')
            primaryHDU.header.set('NAXIS', comment='Number of data axes')
            primaryHDU.header.set('NAXIS1', comment='Length of data axis 1 (slit dimension)')
            primaryHDU.header.set('NAXIS2', comment='Length of data axis 2 (scanning steps)')
            primaryHDU.header.set('NAXIS3', comment='Length of data axis 3 (wavelength dimension)')
            primaryHDU.header.add_comment('Spectral curvature corrected')
            primaryHDU.header.add_comment('Dark subtracted')
            primaryHDU.header.add_comment('Flat-field corrected')
            primaryHDU.header.add_comment('Processed by RSM_prep')
            primaryHDU.writeto(OUT_DIR + 'RSM' + temp_dict['start_time'].strftime('%Y%m%dT%H%M%S') + '_' +
                               temp_dict['scan_index'] + '_FE.fits', overwrite=True)
        except BaseException as uniformException:
            suntools.log(uniformException)
            suntools.log("当前序列输出错误, 已跳过")
        if_first_print.value = True
        INCREASE_DICT_INDEX()

    time_end = time.time()
    suntools.log('并行进度已完成，所花费时间为：', (time_end - time_start) / 60, 'min(分钟)')
    suntools.log('生成视频中...')
    createVideo.createVideo(global_multiprocess_list[0]['start_time'])
    suntools.log('程序结束！')


if __name__ == "__main__":
    main()
