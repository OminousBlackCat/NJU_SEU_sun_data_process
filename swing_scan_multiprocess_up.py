"""
本py文件对应安装后的run_scan_solar文件 使用source指令可以执行该主程序
本程序对应摆扫扫描模式下对图像数据的预处理操作
与multiprocess.py的主要区别体现在对读入文件列表的处理 增添了对每一个扫描序列的判断
注意在使用前修改config.py内的参数!

@author: seu_wxy
"""

import multiprocessing as mp
import datetime
from math import pi, sqrt, acos, atan, ceil, floor
from sympy.external.tests.test_scipy import scipy

import header
import sim
import suntools
import numpy as np
from astropy.io import fits
import urllib.error as uEr
import config
import matplotlib.pyplot as plt
import save_png_video
import traceback
import sys
import os
import time

# 调试用
np.set_printoptions(threshold=sys.maxsize)

# （调试用）控制台输出记录到文件
# class Logger(object):
#     def __init__(self, file_name="Default.log", stream=sys.stdout):
#         self.terminal = stream
#         self.log = open(file_name, "a")
#
#     def write(self, message):
#         self.terminal.write(message)
#         self.log.write(message)
#
#     def flush(self):
#         pass
#
#
#
#     # 自定义目录存放日志文件
# log_path = './Logs/'
# if not os.path.exists(log_path):
#     os.makedirs(log_path)
#     # 日志文件名按照程序运行时间设置
# log_file_name = log_path + 'log-pipeline' + time.strftime("%Y%m%d-%H%M%S", time.localtime()) + '.log'
# # 记录正常的 print 信息
# sys.stdout = Logger(log_file_name)
# # 记录 traceback 异常信息
# sys.stderr = Logger(log_file_name)


time_start = time.time()

# 读入配置文件 引入参数
GLOBAL_BINNING = config.bin_count  # binning 数值
READ_DIR = config.data_dir_path  # 读文件的文件夹
OUT_DIR = config.save_dir_path  # 输出文件夹
SUM_DIR = config.sum_dir_path  # 汇总结果文件夹
DARK_FITS_FILE = ''  # 暗场文件路径
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
PIXEL_ZERO_RIGHT_COUNT = config.pixel_to_zero_right_count  # 右侧置零区间
PIXEL_ZERO_LEFT_COUNT = config.pixel_to_zero_left_count  # 左侧置零区间
CENTER_MEAN_THRESHOLD = config.center_mean_threshold  # (摆扫)计算序列中心的阈值
CENTER_MEAN_INDEX = config.center_mean_index  # (摆扫)使用的线心位置
CENTER_MEAN_COUNT = config.center_mean_count  # (摆扫)使用的线心数量
REVERSAL_MODE = config.reversal_mode  # (摆扫)翻转模式


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
    has_fts = False
    for filename in arr:
        if filename.split('.')[-1] == 'fts' or filename.split('.')[-1] == 'fits':
            has_fts = True
            # break;
    if has_fts == False:
        raise OSError
    return arr


# 此处的数据均未做共享处理，因为共享数据量并不是很大，在LINUX环境下使用multiprocess并fork()将直接复制此些全局变量
# 预读输入目录
try:
    data_file_lst = read_fits_directory()
except OSError:
    suntools.log(traceback.print_exc())
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
        sort_test = filename.split('-')[0] + filename.split('-')[1] + str(
            int(filename.split('-')[2].split('.')[0])).zfill(8)
        real_data_file_lst.append(filename)
    except BaseException as e:
        suntools.log('<文件:' + filename + '>非法, 已剔除(并未删除硬盘上的文件)')

# 对bin模式进行判断
suntools.log("判断bin模式中...")
tmp_bin_mode = suntools.judgeBinMode(data_file_lst, READ_DIR)
if tmp_bin_mode != -1:
    suntools.log(f"判断bin模式成功, 值为: {tmp_bin_mode}")
    GLOBAL_BINNING = tmp_bin_mode
    config.bin_count = tmp_bin_mode
    import importlib
    importlib.reload(suntools)
else:
    suntools.log("判断bin模式失败...将采用config内读取的bin模式作为默认值尝试继续运行程序")

# 再对依赖bin的参数进行初始化
suntools.log("初始化参数中...")
if GLOBAL_BINNING == 1:
    DARK_FITS_FILE = config.dark_fits_name
    FLAT_FITS_FILE = config.flat_fits_name_bin_1
    SUN_ROW_COUNT = config.sun_row_count_bin_1
    STANDARD_FILE_INDEX = config.standard_offset_index_bin_1
    CURVE_X0 = config.curve_cor_x0_bin_1
    CURVE_C = config.curve_cor_C_bin_1
    WAVE_RESOLUTION = config.wavelength_resolution_bin_1
    SUM_ROW_INDEX_HA = config.sum_row_index_HA_bin_1
    SUM_ROW_INDEX_FE = config.sum_row_index_FE_bin_1
if GLOBAL_BINNING == 2:
    DARK_FITS_FILE = config.dark_fits_name_bin_2
    FLAT_FITS_FILE = config.flat_fits_name_bin_2
    SUN_ROW_COUNT = config.sun_row_count_bin_2
    STANDARD_FILE_INDEX = config.standard_offset_index_bin_2
    CURVE_X0 = config.curve_cor_x0_bin_2
    CURVE_C = config.curve_cor_C_bin_2
    WAVE_RESOLUTION = config.wavelength_resolution_bin_2
    SUM_ROW_INDEX_HA = config.sum_row_index_HA_bin_2
    SUM_ROW_INDEX_FE = config.sum_row_index_FE_bin_2


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
plt.plot(np.array(global_wave_line_strength_list[0: 15000]))
ax = plt.subplot()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.xlabel("frame")
plt.ylabel("average strength")
plt.savefig(os.path.join(OUT_DIR, "fig.png"))
suntools.log("输出plot")
for i in range(len(global_wave_line_strength_list)):
    if last_wave_line_strength != 0 and abs(global_wave_line_strength_list[i] - last_wave_line_strength) > 50:
        suntools.log("出现异常波动值, 已跳过")
        continue
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
        if significant_point_list[significant_point_list.index(point) + 1][1] == 1 and \
                next_frame_index > temp_frame_index:
            symmetry_axis_list.append(int((point[0] +
                                           significant_point_list[significant_point_list.index(point) + 1][0]) / 2))
suntools.log('此文件夹共找到' + str(len(symmetry_axis_list)) + '个序列')
last_symmetry_axis_frame_index = 0
current_scan_index = 0
current_track_index = 0
head_series_of_track = dict()
for axis in symmetry_axis_list:
    l_bias = []
    l_sunpos = []
    l_slitpos = []
    count_miss = []
    l_badqua = []
    count = 0
    time0 = time.time()
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
    # 判断当前序列是否被文件夹截断，如果截断了，则继续处理下一个，不管其他的
    if temp_start_file_index < 0:
        temp_start_file_index = 0
        continue
    if temp_last_file_index > len(data_file_lst) - 1:
        temp_last_file_index = len(data_file_lst) - 1
        continue
    step_start = int(data_file_lst[temp_start_file_index].split('-')[-1].split('.')[0])
    step_start = int(data_file_lst[temp_start_file_index].split('-')[-1].split('.')[0])
    step_mid = int(step_start) + 2312 / config.bin_count  # 对非binning模式的数据，需要修改
    step_end = int(step_start) + floor(4625 / config.bin_count)  # 对非binning模式的数据，需要修改
    # if step_start > axis:
    #     continue
    if abs(step_mid-current_axis_frame_index) > 50:
        continue
    if current_scan_index == 0:
        head_series_of_track[current_track_index] = data_file_lst[temp_start_file_index]
    if str(current_scan_index).zfill(4) == '0000':
        for i in range(temp_start_file_index, temp_last_file_index + 1):
            step = int(data_file_lst[i].split('-')[-1].split('.')[0])  # 读取帧序号，以便判断有无缺帧
            # if step >= step_start and step <= step_end and '-0000-' in data_file_lst[i]:  #
            # 第一个条件判断数据是否为这一个序列的，第二个条件判断是否为0级数据的fits文件,第二个条件根据文件名格式可能需要改动，或者可以删除，整合进pipeline时可以检查一下）
            if count == 0:
                dif_step = step - count  # 后续用于判断是否缺帧（如果有更好的方法可以修改）
            # read fits
            hdu = fits.open(READ_DIR + data_file_lst[i])  # 开始读取0级数据
            q0 = hdu[0].header['Q0']  # 读取四元数（下同）
            q1 = hdu[0].header['Q1']
            q2 = hdu[0].header['Q2']
            q3 = hdu[0].header['Q3']
            strtime = hdu[0].header['STR_TIME']  # 读取该帧的时间
            satx_J2000 = hdu[0].header['SAT_POS1']  # 读取卫星相对地球位置（下同）
            saty_J2000 = hdu[0].header['SAT_POS2']
            satz_J2000 = hdu[0].header['SAT_POS3']

            if abs(q0 ** 2 + q1 ** 2 + q2 ** 2 + q3 ** 2 - 1) >= 0.000001:  # 判断四元数是否错误，若有误，执行以下步骤
                l_q = [q0, q1, q2, q3]  # 记录错误的四元数
                l_badqua.append(count)  # 记录错误的四元数所在的位置
                l_qd = [abs(q0 - old_q0), abs(q1 - old_q1), abs(q2 - old_q2),
                        abs(q3 - old_q3)]  # 计算这组四元数与上一组（没有错误的四元数）的相差（理论上相差不大，相差最大的那个数即为出问题的数）
                qbad_i = l_qd.index(max(l_qd))  # 找出出现问题的四元数
                q_squaresum = 0
                for qi in range(4):
                    if qi != qbad_i:
                        q_squaresum += l_q[qi] ** 2
                l_q[qbad_i] = sqrt(1 - q_squaresum)  # 把出问题的四元数暂时替换，使得四元数的平方和为1
                q0, q1, q2, q3 = l_q[0], l_q[1], l_q[2], l_q[3]
            else:  # 四元数没问题时，记录下该组四元数（默认第一帧的四元数不存在问题，后续若出问题，再调整）
                old_q0 = q0
                old_q1 = q1
                old_q2 = q2
                old_q3 = q3

            Chase_arr_Ea = suntools.quaternion_rot(np.array([q0, q1, q2, q3]),
                                                   np.array(
                                                       [0, 1,
                                                        0]))  # 用四元数实现旋转，原始的方向是[0, 1, 0]，旋转后即为卫星在坐标系（J2000坐标系）中的真实指向
            Chase_arr_E = np.array([[Chase_arr_Ea[0]], [Chase_arr_Ea[1]], [Chase_arr_Ea[2]]])

            # calculate real Sun position
            sun_pos, earth_pos = suntools.getEphemerisPos(strtime)

            R_sat_eq = np.array([[satx_J2000], [saty_J2000], [satz_J2000]])  # 卫星相对于地球的位置（J2000坐标系）
            satx_eq = float(R_sat_eq[0])
            saty_eq = float(R_sat_eq[1])
            satz_eq = float(R_sat_eq[2])

            # position is in ICRF (aberration not corrected)
            sat2sun_rx = sun_pos.x.value - earth_pos.x.value - satx_eq
            sat2sun_ry = sun_pos.y.value - earth_pos.y.value - saty_eq
            sat2sun_rz = sun_pos.z.value - earth_pos.z.value - satz_eq
            sat2sun_pos = np.array([[sat2sun_rx], [sat2sun_ry], [sat2sun_rz]])  # 太阳相对于卫星的位置（J2000坐标系）

            normalize_factor = sqrt(
                sat2sun_pos[0] ** 2 + sat2sun_pos[1] ** 2 + sat2sun_pos[2] ** 2)  # 归一化因子，也是太阳与卫星的距离
            sat2sun_pos_normalize = sat2sun_pos / normalize_factor  # 由卫星指向太阳的单位向量

            arr_mul = np.multiply(sat2sun_pos_normalize, Chase_arr_E)  # 对应位置相乘

            try:
                bias = acos(arr_mul[0] + arr_mul[1] + arr_mul[2]) * 180 / pi * 60 * 60  # 换算为角秒的，太阳实际位置与卫星指向的偏差角
            except:
                pass

            l_bias.append(bias)
            l_sunpos.append(sat2sun_pos)  # 记录太阳相对卫星的位置
            l_slitpos.append(Chase_arr_E)  # 记录卫星指向（J2000坐标系）
            count += 1

            if step - count + 1 != dif_step:  # 判断是否缺帧
                count_miss.append(count)  # 若缺帧，则记录缺帧的位置
                count += 1

            if step == step_mid:  # 由于后续都是以中间帧位置作为零点，此处记录中间帧的一些参数
                # cf: center frame
                cf_q0, cf_q1, cf_q2, cf_q3 = q0, q1, q2, q3
        slitposa_first = l_slitpos[0]
        slitpos_first = np.array([slitposa_first[0][0], slitposa_first[1][0], slitposa_first[2][0]])
        R3_first = suntools.quaternion_rot(np.array([cf_q0, -cf_q1, -cf_q2, -cf_q3]), slitpos_first)
        slitposa_last = l_slitpos[-1]
        slitpos_last = np.array([slitposa_last[0][0], slitposa_last[1][0], slitposa_last[2][0]])
        R3_last = suntools.quaternion_rot(np.array([cf_q0, -cf_q1, -cf_q2, -cf_q3]), slitpos_last)
        if R3_first[2] > R3_last[2]:
            reverse_scan = 1
        else:
            reverse_scan = 0

        l_biasx = []  # x轴方向上（slit方向），卫星实际指向与理想指向的偏差
        l_biasz = []  # z轴方向上（scan方向），卫星实际指向与理想指向的偏差
        for step in range(len(l_slitpos)):

            ave_res_ang = 0.5218 * GLOBAL_BINNING * (1156 * 2 // GLOBAL_BINNING - step)  # 卫星理想指向（中心帧坐标系）
            slitposa = l_slitpos[step]  # 卫星实际指向（J2000坐标系）
            slitpos = np.array([slitposa[0][0], slitposa[1][0], slitposa[2][0]])

            R3 = suntools.quaternion_rot(np.array([cf_q0, -cf_q1, -cf_q2, -cf_q3]), slitpos)  # 将卫星指向在中心帧坐标系上表示
            slitpos_b = R3

            biasx = atan(slitpos_b[0] / slitpos_b[1]) * 180 / pi * 60 * 60  # x轴方向上（slit方向），卫星实际指向与理想指向的偏差

            if reverse_scan == 1:  # 实际z指向与理想z指向相减，即为实际与理想的偏差，扫描获取的图像，绘图时，采用了理想指向，因此表现为畸变
                biasz = atan(slitpos_b[2] / sqrt(slitpos_b[0] ** 2 + slitpos_b[1] ** 2)) * 180 / pi * 60 * 60 \
                        - ave_res_ang  # in arcsec
            else:
                biasz = atan(slitpos_b[2] / sqrt(slitpos_b[0] ** 2 + slitpos_b[1] ** 2)) * 180 / pi * 60 * 60 \
                        + ave_res_ang  # in arcsec

            l_biasx.append(biasx)  # 记录下指向的偏差（下同）
            l_biasz.append(biasz)

        count = 0
        # for cm in count_miss:  # 指向偏差列表的长度不一定是2313，有缺帧时，前一帧和后一帧的平均值填充进去
        #     l_biasx.insert(cm + count, (l_biasx[cm + count - 1] + l_biasx[cm + count + 1]) / 2)
        #     l_biasz.insert(cm + count, (l_biasz[cm + count - 1] + l_biasz[cm + count + 1]) / 2)
        #     count += 1
        for bq in l_badqua:  # 处理四元数中存在错误的点
            countless = 0
            countmore = 0
            while bq - countless in l_badqua:
                countless += 1
            while bq + countmore in l_badqua:
                countmore += 1
            l_biasx[bq] = (l_biasx[bq - countless] + l_biasx[bq + countmore]) / 2  # 用前后四元数正确的时刻算出的指向的平均值填充
            l_biasz[bq] = (l_biasz[bq - countless] + l_biasz[bq + countmore]) / 2
        l_biasx_smooth = scipy.signal.savgol_filter(l_biasx, 299, 3)  # 因为四元数的时间分辨率为11帧或14帧，需要将指向平滑一下（下同）
        l_biasz_smooth = scipy.signal.savgol_filter(l_biasz, 299, 3)
        plt.figure(figsize=(20, 20))
        label = "Bias of Track " + str(current_track_index) + " Series " + str(
            current_scan_index).zfill(4) + " \n Time of series 0000 of this track: " + \
                head_series_of_track[current_track_index].split('T')[0][-8:] + 'T' + \
                head_series_of_track[current_track_index].split('T')[1][
                :6]
        plt.suptitle(label)
        plt.subplot(321)
        plt.plot(l_bias)  # 太阳实际位置和卫星指向的角度差
        plt.xlabel('Step')
        plt.ylabel('Slit Pointing (arcsec)')
        plt.grid()

        plt.subplot(323)
        plt.plot(l_biasx, color='red')  # 卫星的实时指向和卫星中间帧指向的差别
        plt.plot(l_biasx_smooth, color='blue')
        plt.xlabel('Step')
        plt.ylabel('Slit-x Pointing (arcsec)')
        plt.grid()

        plt.subplot(325)
        plt.plot(l_biasz, color='red')  # 卫星的实时指向和卫星中间帧指向的差别
        plt.plot(l_biasz_smooth, color='blue')
        plt.xlabel('Step')
        plt.ylabel('Slit-z Pointing (arcsec)')
        plt.grid()

        plt.subplot(324)
        plt.plot(l_biasx_smooth, color='blue')  # 卫星的实时指向和卫星中间帧指向的差别
        plt.xlabel('Step')
        plt.ylabel('Slit-x Pointing (arcsec)')
        plt.grid()

        plt.subplot(326)
        plt.plot(l_biasz_smooth, color='blue')  # 卫星的实时指向和卫星中间帧指向的差别
        plt.xlabel('Step')
        plt.ylabel('Slit-z Pointing (arcsec)')
        plt.grid()

        # if not os.path.exists(config.bias_dir_path):  # 检查目录是否存在
        #     os.makedirs(config.bias_dir_path)  # 如果不存在则创建目录
        # plt.savefig(config.bias_dir_path + 'T' + str(current_track_index).zfill(4) + 'S' +
        #             str(current_scan_index).zfill(4) + data_file_lst[temp_start_file_index].split('T')[0][-8:] + 'T' +
        #             data_file_lst[temp_start_file_index].split('T')[1][
        #             :6] + '.png')  # 需要将四元数指向示意图暂存在某一路径，可供后续检查数据质量
        if reverse_scan == 1:
            for i in range(floor(len(l_biasx) / 2)):
                temp = l_biasx_smooth[-i - 1]
                l_biasx_smooth[-i - 1] = l_biasx_smooth[i]
                l_biasx_smooth[i] = temp
                temp = l_biasz_smooth[-i - 1]
                l_biasz_smooth[-i - 1] = l_biasz_smooth[i]
                l_biasz_smooth[i] = temp

    if str(current_scan_index).zfill(4) != '0000':
        for i in range(temp_start_file_index, temp_last_file_index + 1):
            if i != temp_start_file_index and i != temp_last_file_index and i != axis:
                continue
            # read fits
            hdu = fits.open(READ_DIR + data_file_lst[i])  # 开始读取0级数据
            q0 = hdu[0].header['Q0']  # 读取四元数（下同）
            q1 = hdu[0].header['Q1']
            q2 = hdu[0].header['Q2']
            q3 = hdu[0].header['Q3']
            strtime = hdu[0].header['STR_TIME']  # 读取该帧的时间
            satx_J2000 = hdu[0].header['SAT_POS1']  # 读取卫星相对地球位置（下同）
            saty_J2000 = hdu[0].header['SAT_POS2']
            satz_J2000 = hdu[0].header['SAT_POS3']

            if abs(q0 ** 2 + q1 ** 2 + q2 ** 2 + q3 ** 2 - 1) >= 0.000001:  # 判断四元数是否错误，若有误，执行以下步骤
                l_q = [q0, q1, q2, q3]  # 记录错误的四元数
                l_badqua.append(count)  # 记录错误的四元数所在的位置
                l_qd = [abs(q0 - old_q0), abs(q1 - old_q1), abs(q2 - old_q2),
                        abs(q3 - old_q3)]  # 计算这组四元数与上一组（没有错误的四元数）的相差（理论上相差不大，相差最大的那个数即为出问题的数）
                qbad_i = l_qd.index(max(l_qd))  # 找出出现问题的四元数
                q_squaresum = 0
                for qi in range(4):
                    if qi != qbad_i:
                        q_squaresum += l_q[qi] ** 2
                l_q[qbad_i] = sqrt(1 - q_squaresum)  # 把出问题的四元数暂时替换，使得四元数的平方和为1
                q0, q1, q2, q3 = l_q[0], l_q[1], l_q[2], l_q[3]
            else:  # 四元数没问题时，记录下该组四元数（默认第一帧的四元数不存在问题，后续若出问题，再调整）
                old_q0 = q0
                old_q1 = q1
                old_q2 = q2
                old_q3 = q3

            Chase_arr_Ea = suntools.quaternion_rot(np.array([q0, q1, q2, q3]),
                                                   np.array(
                                                       [0, 1,
                                                        0]))  # 用四元数实现旋转，原始的方向是[0, 1, 0]，旋转后即为卫星在坐标系（J2000坐标系）中的真实指向
            Chase_arr_E = np.array([[Chase_arr_Ea[0]], [Chase_arr_Ea[1]], [Chase_arr_Ea[2]]])

            # calculate real Sun position
            sun_pos, earth_pos = suntools.getEphemerisPos(strtime)

            R_sat_eq = np.array([[satx_J2000], [saty_J2000], [satz_J2000]])  # 卫星相对于地球的位置（J2000坐标系）
            satx_eq = float(R_sat_eq[0])
            saty_eq = float(R_sat_eq[1])
            satz_eq = float(R_sat_eq[2])

            # position is in ICRF (aberration not corrected)
            sat2sun_rx = sun_pos.x.value - earth_pos.x.value - satx_eq
            sat2sun_ry = sun_pos.y.value - earth_pos.y.value - saty_eq
            sat2sun_rz = sun_pos.z.value - earth_pos.z.value - satz_eq
            sat2sun_pos = np.array([[sat2sun_rx], [sat2sun_ry], [sat2sun_rz]])  # 太阳相对于卫星的位置（J2000坐标系）

            normalize_factor = sqrt(
                sat2sun_pos[0] ** 2 + sat2sun_pos[1] ** 2 + sat2sun_pos[2] ** 2)  # 归一化因子，也是太阳与卫星的距离
            sat2sun_pos_normalize = sat2sun_pos / normalize_factor  # 由卫星指向太阳的单位向量

            arr_mul = np.multiply(sat2sun_pos_normalize, Chase_arr_E)  # 对应位置相乘

            try:
                bias = acos(arr_mul[0] + arr_mul[1] + arr_mul[2]) * 180 / pi * 60 * 60  # 换算为角秒的，太阳实际位置与卫星指向的偏差角
            except:
                pass

            l_bias.append(bias)
            l_sunpos.append(sat2sun_pos)  # 记录太阳相对卫星的位置
            l_slitpos.append(Chase_arr_E)  # 记录卫星指向（J2000坐标系）
            count += 1

            if i == axis:  # 由于后续都是以中间帧位置作为零点，此处记录中间帧的一些参数
                # cf: center frame
                cf_q0, cf_q1, cf_q2, cf_q3 = q0, q1, q2, q3
        slitposa_first = l_slitpos[0]
        slitpos_first = np.array([slitposa_first[0][0], slitposa_first[1][0], slitposa_first[2][0]])
        R3_first = suntools.quaternion_rot(np.array([cf_q0, -cf_q1, -cf_q2, -cf_q3]), slitpos_first)
        slitposa_last = l_slitpos[-1]
        slitpos_last = np.array([slitposa_last[0][0], slitposa_last[1][0], slitposa_last[2][0]])
        R3_last = suntools.quaternion_rot(np.array([cf_q0, -cf_q1, -cf_q2, -cf_q3]), slitpos_last)
        if R3_first[2] > R3_last[2]:
            reverse_scan = 1
        else:
            reverse_scan = 0

    if str(current_scan_index).zfill(4) == '0000':
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
            'start_time': datetime.datetime.now(),
            'bias_x': l_biasx_smooth,
            'bias_z': l_biasz_smooth,
            'reverse_scan': reverse_scan
        })

    if str(current_scan_index).zfill(4) != '0000':
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
            'start_time': datetime.datetime.now(),
            'reverse_scan': reverse_scan
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
    suntools.log(traceback.print_exc())
    suntools.log("Error: 暗场文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    suntools.log(traceback.print_exc())
    suntools.log("Error: 暗场文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")
if temp_img is not None:
    dark_img = np.array(temp_img[0].data, dtype=float)
    # dark_img = suntools.change(dark_img)
temp_img.close()

# 平场需要以日心图片作为基准进行平移矫正 再进行谱线弯曲矫正
flat_img = None
standard_HA_width, standard_FE_width = None, None
try:
    suntools.log("正在读取原始平场文件")
    temp_img = fits.open(FLAT_FITS_FILE)
except uEr.URLError:
    suntools.log(traceback.print_exc())
    suntools.log("Error: 原始平场文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    suntools.log(traceback.print_exc())
    suntools.log("Error: 原始平场文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")
if temp_img is not None:
    flat_img = np.array(temp_img[0].data, dtype=float)
    flat_img, standard_HA_width, standard_FE_width = suntools.curve_correction(flat_img - dark_img, CURVE_X0,
                                                                               CURVE_C)
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
        suntools.log(
            '校正扫描第' + str(temp_dict['track_index']) + '轨, 序列' + temp_dict[
                'scan_index'] + '中...使用标准校正文件为:' + temp_dict[
                'standard_filename'])
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
        # 对平场进行归一化
        flatTemp = suntools.FlatNormalization(flatTemp)
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
    suntools.log(traceback.print_exc())
    suntools.log("Error: 标准日心校准文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError as error:
    suntools.log(traceback.print_exc())
    suntools.log("Error: 标准日心校准文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")

# 读取输出色谱
color_map = suntools.get_color_map(COLOR_CAMP_FILE)

# 需要创建一个四维数组
# 对于最后写入的文件 NAXIS顺序为: 狭缝宽度 文件序号(扫描序号)  波长深度
# 因此np array的shape应为(波长深度, 文件序号, 狭缝宽度)
# 每个单独文件对应的是 xy平面的一个二维数组
GLOBAL_ARRAY_X_COUNT = sample_from_standard.shape[0]
GLOBAL_ARRAY_Y_COUNT = SUN_ROW_COUNT
GLOBAL_ARRAY_Z_COUNT = sample_from_standard.shape[1]
suntools.log('SHAPE:' + str(GLOBAL_ARRAY_X_COUNT) + ',' + str(GLOBAL_ARRAY_Y_COUNT) + ',' + str(GLOBAL_ARRAY_Z_COUNT))


# 多进程并行，对于每一序列的0000和0001要分开来处理
def multiprocess_task(parameter_dic: dict):
    """
    按照序列并行, 每次都传进一个dict作为参数列表, 保存了处理此序列所需的全部信息
    每个进程都会开辟一块内存空间用来存储数据
    处理结束之后便会输出文件
    """
    # 单个序列生成的全日面太阳像, shape = [ 波长深度 , 文件序号 , 狭缝宽度 ]
    sequence_data_array = None
    # 单个序列生成的时间矩阵, shape = [ 文件序号, 狭缝宽度 ]
    time_series_data_array = None
    suntools.log(
        '正在处理第' + str(parameter_dic['track_index']) + '轨,  扫描序列:' + parameter_dic['scan_index'] + '...')
    try:
        sequence_data_array = np.zeros((GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT),
                                       dtype=np.int16)
        time_series_data_array = np.zeros((GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT),
                                          dtype=np.int16)
    except BaseException as uniformException:
        suntools.log(traceback.print_exc())
        suntools.log("内存不足, 无法创建内存空间")
        suntools.log(str(uniformException))
        sys.exit("程序终止")

    # 一个标准文件名 如下:
    # RSM 2021   12     22T060105   -   0008-     0001       .fts
    # 012 3456   78     901234567   8   90123     4567       8901
    #     [year] [mon]  [day_seq]       [index]   [position]
    try:
        # 年无需设置, 其余均需要设置以防止序列时间跨日
        start_datetime = datetime.datetime(
            year=int(parameter_dic['first_filename'][3: 7]),
            month=int(parameter_dic['first_filename'][7: 9]),
            day=int(parameter_dic['first_filename'][9: 11]),
            hour=int(parameter_dic['first_filename'][12: 14]),
            minute=int(parameter_dic['first_filename'][14: 16]),
            second=int(parameter_dic['first_filename'][16: 18]),
            microsecond=0
        )
    except BaseException as uniformException:
        suntools.log(traceback.print_exc())
        suntools.log("格式化文件名出错")
        sys.exit("程序终止")

    for sequence_filename in parameter_dic['file_list']:
        # 对序列内每个狭缝文件进行预处理, 组织为一个全日面太阳像与时间矩阵
        try:
            # 求当前狭缝文件的日期
            current_datetime = datetime.datetime(
                year=int(sequence_filename[3: 7]),
                month=int(sequence_filename[7: 9]),
                day=int(sequence_filename[9: 11]),
                hour=int(sequence_filename[12: 14]),
                minute=int(sequence_filename[14: 16]),
                second=int(sequence_filename[16: 18]),
                microsecond=0
            )
            # 求时间差值, 获取相对值
            relative_time_value = (current_datetime - start_datetime).total_seconds()
            # 计算此文件在序列中的相对位置 给后续放入全局数组做准备
            fileRelativePosition = int(sequence_filename.split('-')[2].split('.')[0]) - int(
                parameter_dic['first_filename'].split('-')[2].split('.')[0])
            if fileRelativePosition < 0:
                fileRelativePosition = 0
            if fileRelativePosition >= SUN_ROW_COUNT:
                fileRelativePosition = SUN_ROW_COUNT - 1
            filePath = READ_DIR + sequence_filename
            file_data = fits.open(filePath)
            image_data = np.array(file_data[0].data, dtype=float)
            # 对fe窗口进行平移
            image_data = suntools.moveImg(image_data, -2)
            # 去暗场
            image_data = image_data - dark_img
            # 谱线弯曲矫正
            image_data, HofH, HofFe = suntools.curve_correction(image_data, CURVE_X0, CURVE_C)
            # 搜索list
            currentFlat = parameter_dic['flat_data']
            currentScanIndex = int(parameter_dic['scan_index'])
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
            image_data[:, image_data.shape[1] - PIXEL_ZERO_RIGHT_COUNT:] = 0
            image_data[:, 0: PIXEL_ZERO_LEFT_COUNT] = 0
            # 原来的上下偏转，用南大提供的新方法替换
            if parameter_dic['reverse_scan'] == 1:
                reverse_index = SUN_ROW_COUNT - 1 if SUN_ROW_COUNT - 1 - fileRelativePosition >= SUN_ROW_COUNT else SUN_ROW_COUNT - 1 - fileRelativePosition
                sequence_data_array[:, reverse_index, :] = image_data
                time_series_data_array[reverse_index, :] = int(relative_time_value)
            # elif REVERSAL_MODE == 'even' and currentScanIndex % 2 == 0:
            #     sequence_data_array[:, SUN_ROW_COUNT - 1 - fileRelativePosition, :] = image_data
            else:
                sequence_data_array[:, fileRelativePosition, :] = image_data
                time_series_data_array[fileRelativePosition, :] = int(relative_time_value)

            file_data.close()
        except BaseException:
            suntools.log(traceback.print_exc())
            suntools.log(sequence_filename, parameter_dic['first_filename'], fileRelativePosition)

            suntools.log('文件:' + filename + '处理失败, 请检查此文件')

    suntools.log(
        '第' + str(parameter_dic['track_index']) + '轨, 扫描序列' + parameter_dic['scan_index'] + '预处理完成...')
    suntools.log('生成完整日像中...')
    try:
        # sum_data_HA = np.zeros((SUN_ROW_COUNT, sample_from_standard.shape[1]))
        # sum_data_FE = np.zeros((SUN_ROW_COUNT, sample_from_standard.shape[1]))
        p0 = parameter_dic['header']['INST_ROT']
        strtime = parameter_dic['header']['STR_TIME']
        se0000_hacore = np.array(sequence_data_array[SUM_ROW_INDEX_HA, :, :])
        h, w = se0000_hacore.shape  # 读取图片高度和宽度
        se0000_rx = parameter_dic['header']['SAT_POS1']
        se0000_ry = parameter_dic['header']['SAT_POS2']
        se0000_rz = parameter_dic['header']['SAT_POS3']
        satpos = [se0000_rx, se0000_ry, se0000_rz]
        radius_ref = sim.theory_rsun(strtime, satpos, GLOBAL_BINNING)
        parameter_dic['header'].set('RSUN_REF', radius_ref)
        suntools.log('计算获得太阳日心坐标理论半径为:', radius_ref)
        x_width = sequence_data_array.shape[1]
        z_width = sequence_data_array.shape[2]
        axis_width_ha = [standard_HA_width, x_width, z_width]
        axis_width_fe = [standard_FE_width, x_width, z_width]
        if parameter_dic['is_head_of_track']:
            biasx = parameter_dic['bias_x'] / (0.5218 * GLOBAL_BINNING)  # 读取卫星指向偏差（此时是角秒单位）
            biasz = parameter_dic['bias_z'] / (0.5218 * GLOBAL_BINNING)  # 读取卫星指向偏差（此时是角秒单位）
            se00xx_imwing_ha = suntools.head_distortion_correction('HA', axis_width_ha, biasx, biasz,
                                                                   sequence_data_array[0:standard_HA_width, :, :],
                                                                   time_series_data_array=time_series_data_array)
            se00xx_imwing_fe = suntools.head_distortion_correction('FE', axis_width_fe, biasx, biasz,
                                                                   sequence_data_array[standard_HA_width:, :, :])
            se0000_hacore0 = np.array(sequence_data_array[SUM_ROW_INDEX_HA, :, :])
            se0000_hawing0 = np.array(sequence_data_array[SUM_ROW_INDEX_FE, :, :])
            se0000_center = sim.circle_center(se0000_hawing0)
            if parameter_dic['scan_index'] == '0000':
                parameter_dic['global_track_se0000_center'][parameter_dic['track_index']] = se0000_center
                parameter_dic['global_track_se0000_hacore'][parameter_dic['track_index']] = se0000_hacore0

        else:
            hacore0 = parameter_dic['global_track_se0000_hacore'][parameter_dic['track_index']]
            # 注意使用deep copy
            se00xx_hacore = np.array(sequence_data_array[SUM_ROW_INDEX_HA, :, :])
            se00xx_RSUN = sim.theory_rsun(strtime, satpos, GLOBAL_BINNING)
            se0000_center = parameter_dic['global_track_se0000_center'][parameter_dic['track_index']]
            se00xx_imwing_ha = suntools.non_head_distortion_correction('HA',
                                                                       sequence_data_array[0:standard_HA_width, :, :],
                                                                       se0000_center, se00xx_RSUN, axis_width_ha,
                                                                       se00xx_hacore, hacore0, w, h,
                                                                       time_series_data_array=time_series_data_array)
            se00xx_imwing_fe = suntools.non_head_distortion_correction('FE',
                                                                       sequence_data_array[standard_HA_width:, :, :],
                                                                       se0000_center, se00xx_RSUN, axis_width_fe,
                                                                       se00xx_hacore, hacore0, w, h)
        suntools.log('开始旋转图像...')
        # 对Ha图像进行旋转
        centerx_ha, centery_ha = suntools.rotate_fits(standard_HA_width,
                                                      x_width, z_width,
                                                      se00xx_imwing_ha,
                                                      radius_ref,
                                                      sequence_data_array[
                                                      0:standard_HA_width,
                                                      :, :], p0, True, time_series_data_array=time_series_data_array)
        # 对Fe图像进行旋转
        centerx_fe, centery_fe = suntools.rotate_fits(
            sequence_data_array.shape[0] - standard_HA_width, x_width, z_width, se00xx_imwing_fe,
            radius_ref, sequence_data_array[standard_HA_width:, :, :], p0, False)
        parameter_dic['header'].set('CRPIX1', (centerx_ha + centerx_fe) / 2)
        parameter_dic['header'].set('CRPIX2', (centery_ha + centery_fe) / 2)
        parameter_dic['header'].set('INST_ROT', 0)

        suntools.log('旋转完成！')

        # # 输出太阳像
        # for seq_index in range(SUN_ROW_COUNT):
        #     sum_data_HA = np.array(sequence_data_array[SUM_ROW_INDEX_HA, :, :])
        #     sum_data_FE = np.array(sequence_data_array[standard_HA_width + SUM_ROW_INDEX_FE, :, :])
        # suntools.log('计算CCD太阳像半径中...')
        R_y, R_x, radius = suntools.getCircle(sequence_data_array[standard_HA_width + SUM_ROW_INDEX_FE, :, :])
        # OBS_Radius = radius * PIXEL_RESOLUTION * GLOBAL_BINNING
        suntools.log('波长定标中...')
        wavelength_calibrate_input = np.array(sequence_data_array[:, int(R_y) - 50: int(R_y) + 49,
                                     int(R_x) - 50: int(R_x) + 49])
        cdel_t3, crval_l3_ha, crval_l3_fe = suntools.cal_center_mean(wavelength_calibrate_input)
        parameter_dic['header'].set('R_SUN', radius)
        parameter_dic['header'].set('RSUN_OBS', radius * 0.5218 * GLOBAL_BINNING)
        parameter_dic['header'].set('CDELT1', PIXEL_RESOLUTION * GLOBAL_BINNING)
        parameter_dic['header'].set('CDELT2', PIXEL_RESOLUTION * GLOBAL_BINNING)
        parameter_dic['header'].set('CDELT3', cdel_t3)
        parameter_dic['header'].set('CRVAL3', crval_l3_ha)
        # 下采样 1/4
        # suntools.log('下采样中...')
        # sum_data_HA_save = suntools.down_sample(sum_data_HA)
        # sum_data_FE_save = suntools.down_sample(sum_data_FE)
        # 旋转INST_ROT度（已改用南大提供的方法，这里目前弃用）
        # sum_data_HA_save = ndimage.rotate(sum_data_HA_save, -parameter_dic['header']['INST_ROT'], reshape=False)
        # sum_data_FE_save = ndimage.rotate(sum_data_FE_save, -parameter_dic['header']['INST_ROT'], reshape=False)
        # sum_mean_ha = np.mean(sum_data_HA)
        # sum_mean_fe = np.mean(sum_data_FE)
        # sum_data_HA = suntools.cropPNG(sum_data_HA, int(R_x), int(R_y))
        # sum_data_FE = suntools.cropPNG(sum_data_FE, int(R_x), int(R_y))
        # sum_data_HA_save = suntools.add_time(sum_data_HA, parameter_dic['start_time'].
        #                                      strftime('%Y-%m-%d ''%H:%M:%S UT'), 3 * sum_mean_ha)
        # sum_data_FE_save = suntools.add_time(sum_data_FE, parameter_dic['start_time'].
        #                                      strftime('%Y-%m-%d ''%H:%M:%S UT'), 3 * sum_mean_fe)

        # 将小于0的值全部赋为0
        sequence_data_array[sequence_data_array < 0] = 0
        suntools.log("SHAPE为：" + str(sequence_data_array.shape))
        # if config.save_img_form == 'default':
        #     # 使用读取的色谱进行输出 imsave函数将自动对data进行归一化
        #     suntools.log('输出序号为' + parameter_dic['scan_index'] + '的png...')
        #     METADATA = {'CenterX': str(round(centerx_ha)), 'CenterY': str(round(centery_ha))}
        #     plt.imsave(SUM_DIR + 'RSM' + parameter_dic['start_time'].strftime('%Y%m%dT%H%M%S')
        #                + '_' + parameter_dic['scan_index'] + '_HA' + ".png",
        #                sum_data_HA_save, cmap=color_map, vmin=0, vmax=3 * sum_mean_ha, metadata=METADATA)
        #     plt.imsave(SUM_DIR + 'RSM' + parameter_dic['start_time'].strftime('%Y%m%dT%H%M%S')
        #                + '_' + parameter_dic['scan_index'] + '_FE' + ".png",
        #                sum_data_FE_save, cmap=color_map, vmin=0, vmax=3 * sum_mean_fe, metadata=METADATA)
        # if config.save_img_form == 'fts':
        #     # 不对data进行任何操作 直接输出为fts文件
        #     suntools.log('输出序号为' + parameter_dic['scan_index'] + '的fits...')
        #     primaryHDU = fits.PrimaryHDU(sum_data_HA)
        #     greyHDU = fits.HDUList([primaryHDU])
        #     greyHDU.writeto(SUM_DIR + 'SUM' + parameter_dic['start_time'].strftime('%Y%m%dT%H%M%S')
        #                     + '_' + parameter_dic['scan_index'] + '_HA' + '.fts', overwrite=True)
        #     greyHDU.close()
        #     primaryHDU = fits.PrimaryHDU(sum_data_FE)
        #     greyHDU = fits.HDUList([primaryHDU])
        #     greyHDU.writeto(SUM_DIR + 'SUM' + parameter_dic['start_time'].strftime('%Y%m%dT%H%M%S')
        #                     + '_' + parameter_dic['scan_index'] + '_FE' + '.fts', overwrite=True)
        #     greyHDU.close()
        suntools.log('生成HA文件中...')
        parameter_dic['header'].set('SPECLINE', 'HA')
        parameter_dic['header'].set('WAVE_LEN', HA_LINE_CORE)
        parameter_dic['header'].set('PRODATE', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))

        primaryHDU = fits.CompImageHDU(sequence_data_array[0: standard_HA_width, :, :],
                                       header=parameter_dic['header'], compression_type='RICE_1')
        time_data_HDU = fits.PrimaryHDU(time_series_data_array)
        primaryHDU.header.set('NAXIS', comment='Number of data axes')
        primaryHDU.header.set('NAXIS1', comment='Length of data axis 1 (slit dimension)')
        primaryHDU.header.set('NAXIS2', comment='Length of data axis 2 (scanning steps)')
        primaryHDU.header.set('NAXIS3', comment='Length of data axis 3 (wavelength dimension)')
        primaryHDU.header.add_comment('Spectral curvature corrected')
        primaryHDU.header.add_comment('Dark subtracted')
        primaryHDU.header.add_comment('Flat-field corrected')
        primaryHDU.header.add_comment('Processed by RSM_prep')
        ha_hdu_list = fits.HDUList([time_data_HDU, primaryHDU])
        ha_hdu_list.writeto(OUT_DIR + 'RSM' + parameter_dic['start_time'].strftime('%Y%m%dT%H%M%S') + '_' +
                            parameter_dic['scan_index'] + '_HA.fits', overwrite=True)
        suntools.log('生成FE文件中...')
        # 修改header内的SPECLINE与WAVELNTH
        parameter_dic['header'].set('SPECLINE', 'FEI')
        parameter_dic['header'].set('WAVE_LEN', FE_LINE_CORE)
        parameter_dic['header'].set('PRODATE', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
        parameter_dic['header'].set('CRVAL3', crval_l3_fe)
        primaryHDU = fits.CompImageHDU(sequence_data_array[standard_HA_width:, :, :],
                                       header=parameter_dic['header'], compression_type='RICE_1')
        primaryHDU.header.set('NAXIS', comment='Number of data axes')
        primaryHDU.header.set('NAXIS1', comment='Length of data axis 1 (slit dimension)')
        primaryHDU.header.set('NAXIS2', comment='Length of data axis 2 (scanning steps)')
        primaryHDU.header.set('NAXIS3', comment='Length of data axis 3 (wavelength dimension)')
        primaryHDU.header.add_comment('Spectral curvature corrected')
        primaryHDU.header.add_comment('Dark subtracted')
        primaryHDU.header.add_comment('Flat-field corrected')
        primaryHDU.header.add_comment('Processed by RSM_prep')
        fe_hdu_list = fits.HDUList([time_data_HDU, primaryHDU])
        fe_hdu_list.writeto(OUT_DIR + 'RSM' + parameter_dic['start_time'].strftime('%Y%m%dT%H%M%S') + '_' +
                            parameter_dic['scan_index'] + '_FE.fits', overwrite=True)
    except BaseException as uniformException:
        traceback.print_exc()
        suntools.log("当前序列输出错误, 已跳过")


def join(pool):
    while True:  # 主进程轮询子进程状态，代替pool.join()
        finish_flag = True  # Pool is finish(True) or not
        for sub_p in pool._pool:
            if sub_p.is_alive():
                finish_flag = False
                break
        if finish_flag:
            error_flag = False  # error(True) or not
            for sub_p in pool._pool:
                if sub_p.exitcode != 0:
                    print('Sub-processes exit error')
                    pool.terminate()
                    flag = True
                    break
            break

def main():
    """
    主函数, 使用pool函数对全局dict进行并行处理
    并在pool.join()之后生成视频
    """
    suntools.log('开启多核并行处理...')
    headseries_of_track = []
    non_headseries_of_track = []
    global_track_se0000_center = mp.Manager().dict()
    global_track_se0000_hacore = mp.Manager().dict()
    # mp.set_start_method("spawn")

    for i in range(len(global_multiprocess_list)):
        if global_multiprocess_list[i]['scan_index'] == '0000':
            global_multiprocess_list[i]['is_head_of_track'] = True
            global_multiprocess_list[i]['global_track_se0000_center'] = global_track_se0000_center
            global_multiprocess_list[i]['global_track_se0000_hacore'] = global_track_se0000_hacore
            headseries_of_track.append((global_multiprocess_list[i]))
        else:
            global_multiprocess_list[i]['is_head_of_track'] = False
            global_multiprocess_list[i]['global_track_se0000_center'] = global_track_se0000_center
            global_multiprocess_list[i]['global_track_se0000_hacore'] = global_track_se0000_hacore
            non_headseries_of_track.append((global_multiprocess_list[i]))

    pool = mp.Pool(processes=multiprocess_count)
    # pool = ProcessPoolExecutor(max_workers=multiprocess_count)
    pool.map(multiprocess_task, headseries_of_track)
    pool.close()
    join(pool)
    # pool.shutdown(wait=True)
    # pool1 = ProcessPoolExecutor(max_workers=multiprocess_count)
    pool1 = mp.Pool(processes=multiprocess_count)
    pool1.map(multiprocess_task, non_headseries_of_track)
    pool1.close()
    # pool.join()
    join(pool1)
    # pool1.shutdown(wait=True)
    time_end = time.time()
    suntools.log('并行进度已完成，所花费时间为：', (time_end - time_start) / 60, 'min(分钟)')
    suntools.log('生成预览图像中...')
    suntools.log('视频请使用脚本生成...')
    save_png_video.monographNJU(SUM_DIR, color_map, image_dpi=config.png_dpi_value)
    # save_png_video.createVideoNJU(SUM_DIR, config.video_dir_path)
    suntools.log('程序结束！')

if __name__ == "__main__":
    main()
