import multiprocessing as mp
import os
import header
import sys
import suntools
import time
import numpy as np
from astropy.io import fits
import scipy.signal as signal
import urllib.error as uEr
import config
import matplotlib.pyplot as plt
import ctypes as c
import sys

# 读入配置文件 引入参数
read_dir = config.data_dir_path
out_dir = config.save_dir_path
multiprocess_count = 1
if config.multiprocess_count != 'default':
    multiprocess_count = config.multiprocess_count
else:
    multiprocess_count = mp.cpu_count() - 4
print('多核并行数:' + str(multiprocess_count))


# 读取数据文件夹所有文件
def read_fits_directory():
    arr = []
    arr = os.listdir(read_dir)
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


# 将读入的文件按照序列分成不同的组
global_multiprocess_list = []
for i in range(len(data_file_lst)):
    filename = data_file_lst[i]
    # 选取图像文件名的最后四个字符作为index
    temp_index = int(filename[-13: -9])
    # 在list中寻找对应的dict
    ifFind = False
    for j in range(len(global_multiprocess_list)):
        if global_multiprocess_list[j]['scan_index'] == temp_index:
            global_multiprocess_list[j]['file_list'].append(filename)
            global_multiprocess_list[j]['file_count'] += 1
            ifFind = True
            if int(filename[24:28]) == config.standard_offset_index:
                global_multiprocess_list[j]['standard_filename'] = filename
            break
    if not ifFind:
        global_multiprocess_list.append({
            'scan_index': temp_index,       # 扫描序号
            'file_list': [],                # 文件名列表
            'file_count': 0,                # 包含的文件数
            'standard_filename': '',        # 标准日心文件名
            'flat_data': None,              # 此序列的校正后平场数据
            'abortion_data': None           # 此序列的校正后红蓝移数据
        })
        global_multiprocess_list[len(global_multiprocess_list) - 1]['file_list'].append(filename)
        global_multiprocess_list[len(global_multiprocess_list) - 1]['file_count'] += 1
        if int(filename[24:28]) == config.standard_offset_index:
            global_multiprocess_list[len(global_multiprocess_list) - 1]['standard_filename'] = filename

# 剔除不完整序列
for i in range(len(global_multiprocess_list)):
    if global_multiprocess_list[i]['file_count'] < config.sun_row_count - 300 or global_multiprocess_list[i]['standard_filename'] is None:
        print('文件夹中包含不完整序列, 序列序号为:' + str(global_multiprocess_list[i]['scan_index']))
        print('本次数据处理将不此序列进行处理.....')
        global_multiprocess_list.remove(global_multiprocess_list[i])

# 读取暗场文件
temp_img = None
dark_img = None
try:
    print("正在读取原始暗场文件")
    temp_img = fits.open(config.dark_fits_name)
except uEr.URLError:
    print("Error: 暗场文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    print("Error: 暗场文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")
if temp_img is not None:
    dark_img = np.array(temp_img[0].data, dtype=float)
    dark_img = suntools.change(dark_img)
temp_img.close()

# 平场需要以日心图片作为基准进行平移矫正 再进行谱线弯曲矫正
flat_img = None
try:
    print("正在读取原始平场文件")
    temp_img = fits.open(config.flat_fits_name)
except uEr.URLError:
    print("Error: 原始平场文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    print("Error: 原始平场文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")
if temp_img is not None:
    flat_img = np.array(temp_img[0].data, dtype=float)
    flat_img, temp1, temp2 = suntools.curve_correction(flat_img - dark_img, config.curve_cor_x0, config.curve_cor_C)
temp_img.close()

# 读取经过日心的图片 作为基准
# 读取标准太阳光谱数据
sun_std = suntools.get_Sunstd(config.sun_std_name)
sample_from_standard = None
try:
    for temp_dict in global_multiprocess_list:
        # 对每个序列进行校正
        print('校正扫描序列' + str(temp_dict['scan_index']) + '中...使用标准校正文件为:' + temp_dict['standard_filename'])
        print("校正平场中...")
        standard_name = temp_dict['standard_filename']
        temp_img = fits.open(read_dir + '/' + standard_name)
        standard_img = np.array(temp_img[0].data, dtype=float)
        standard_img = suntools.moveImg(standard_img, -2)
        standard_img, temp1, temp2 = suntools.curve_correction(standard_img - dark_img, config.curve_cor_x0,
                                                               config.curve_cor_C)
        sample_from_standard = standard_img
        # 先平移矫正 减去暗场 再谱线弯曲矫正
        flatTemp = suntools.getFlatOffset(flat_img, standard_img)
        flatTemp = suntools.getFlat(flatTemp)
        print("获得标准太阳光谱数据中...")
        # 以标准文件作为基准 计算红蓝移吸收系数
        # 需要先对标注文件进行一系列操作 去暗场 去平场 再进行红蓝移修正
        standard_img = suntools.DivFlat(standard_img, flatTemp)
        # 获得标准吸收系数
        abortion = suntools.RB_getdata(standard_img, sun_std, temp1, temp2)
        temp_dict['flat_data'] = flatTemp
        temp_dict['abortion_data'] = abortion
        temp_img.close()
        print("序列:" + str(int(standard_name[19:23])) + "矫正完成")
except uEr.URLError:
    print("Error: 标准日心校准文件未找到, 请检查config文件或存放目录")
    sys.exit("程序终止")
except OSError:
    print("Error: 标准日心校准文件读取发生错误, 请检查文件读取权限")
    sys.exit("程序终止")

# 读取输出色谱
color_map = suntools.get_color_map(config.color_camp_name)

# 检查输出文件夹是否存在 不存在则创建
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# 读取头部参数文件
HA_HEADER = fits.header.Header()
FE_HEADER = fits.header.Header()
HA_list = header.read_header_from_txt(config.HA_header_file)
FE_list = header.read_header_from_txt(config.FE_header_file)
# 将读入的list构造成为header
for h in HA_list:
    HA_HEADER.set(h['key'], value=0, comment=h['comment'])
for h in FE_list:
    FE_HEADER.set(h['key'], value=0, comment=h['comment'])

# 全局进度控制
file_count = mp.Value('i', len(read_fits_directory()))
remaining_count = mp.Value('i', 0)
if_first_print = mp.Value('b', True)

# 全局共享内存
# 需要创建一个四维数组
# 三维 xyz 分别为文件序号(扫描序号) 狭缝宽度 与 波长深度
# 每个单独文件对应的是 yz平面的一个二维数组
GLOBAL_ARRAY_X_COUNT = config.sun_row_count
GLOBAL_ARRAY_Y_COUNT = sample_from_standard.shape[0]
GLOBAL_ARRAY_Z_COUNT = sample_from_standard.shape[1]
# 创建共享内存 大小为 x*y*z*sizeof(int16)
GLOBAL_SHARED_MEM = mp.Array(c.c_int16, GLOBAL_ARRAY_X_COUNT * GLOBAL_ARRAY_Y_COUNT * GLOBAL_ARRAY_Z_COUNT)


# 定义target task
# 传入一个文件名，读取此文件名对应的fits文件并对其做曲线矫正
def target_task(filename):
    # 一个标准文件名 如下:
    # RSM 2021   12     22T060105   -   0008-     0001       .fts
    # 012 3456   78     901234567   8   90123     4567       8901
    #     [year] [mon]  [day_seq]       [index]   [position]
    file_index = filename[19:23]
    file_position = filename[24:28]
    filePath = read_dir + "/" + filename
    file_data = fits.open(filePath)
    image_data = np.array(file_data[0].data, dtype=float)
    # 对fe窗口进行平移
    image_data = suntools.moveImg(image_data, -2)
    # 去暗场
    image_data = image_data - dark_img
    # 谱线弯曲矫正
    image_data, HofH, HofFe = suntools.curve_correction(image_data, config.curve_cor_x0, config.curve_cor_C)
    # 搜索list
    currentFlat = None
    currentAbortion = None
    for dataTemp in global_multiprocess_list:
        if dataTemp['scan_index'] == int(file_index):
            currentFlat = dataTemp['flat_data']
            currentAbortion = dataTemp['abortion_data']
            break
    if currentAbortion is None or currentFlat is None:
        print("文件：" + filename + "未找到平场数据与吸收系数, 请检查文件夹")
        return
    # 去平场
    image_data = suntools.DivFlat(image_data, currentFlat)
    # 红蓝移矫正
    image_data = suntools.RB_repair(image_data, currentAbortion)
    # 滤波
    image_data = signal.medfilt(image_data, kernel_size=config.filter_kernel_size)
    # 转为整型
    image_data = np.array(image_data, dtype=np.int16)
    global_shared_array = np.frombuffer(GLOBAL_SHARED_MEM.get_obj(), dtype=np.int16)
    global_shared_array = global_shared_array.reshape(GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT)
    global_shared_array[int(file_position) - 1] = image_data
    # 进度输出
    remaining_count.value += 1
    file_data.close()
    if if_first_print.value:
        print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value), end='')
        sys.stdout.flush()
        if_first_print.value = False
    else:
        print('\b' * (9 + len(str(remaining_count.value)) + 1 + len(str(file_count.value))), end='')
        print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value), end='')
        sys.stdout.flush()


def main():
    # 测试消耗时间 时间起点
    time_start = time.time()
    # 获得文件夹列表 读取相关参数
    # 并行处理
    print('开启多核并行处理...')
    for temp_dict in global_multiprocess_list:
        print('正在处理扫描序列:' + str(temp_dict['scan_index']) + '...')
        file_count.value = temp_dict['file_count']
        remaining_count.value = 0
        pool = mp.Pool(processes=multiprocess_count)
        pool.map(target_task, temp_dict['file_list'])
        pool.close()
        pool.join()
        print('\n扫描序列' + str(temp_dict['scan_index']) + '预处理完成...')
        print('生成完整日像中...')
        sum_data = np.zeros((config.sun_row_count, sample_from_standard.shape[0]))
        global_shared_array = np.frombuffer(GLOBAL_SHARED_MEM.get_obj(), dtype=np.int16)
        global_shared_array = global_shared_array.reshape(GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT,
                                                          GLOBAL_ARRAY_Z_COUNT)
        print("SHAPE为：" + str(global_shared_array.shape))
        for i in range(global_shared_array.shape[0]):
            sum_data[i] = global_shared_array[i, :, config.sum_row_index].reshape(sample_from_standard.shape[0])
        sum_data[sum_data < 0] = 0
        if config.save_img_form == 'default':
            # 使用读取的色谱进行输出 imsave函数将自动对data进行归一化
            print('输出序号为' + str(temp_dict['scan_index']) + '的png...')
            plt.imsave(config.sum_dir_path + 'sum' + str(temp_dict['scan_index']) + ".png", sum_data, cmap=color_map)
        if config.save_img_form == 'fts':
            # 不对data进行任何操作 直接输出为fts文件
            print('输出序号为' + str(temp_dict['scan_index']) + '的fits...')
            primaryHDU = fits.PrimaryHDU(sum_data)
            greyHDU = fits.HDUList([primaryHDU])
            greyHDU.writeto(config.sum_dir_path + 'sum' + str(temp_dict['scan_index']) + '.fts')
            greyHDU.close()
        print('生成HA文件中...')
        file_year = temp_dict['standard_filename'][3:7]
        file_mon = temp_dict['standard_filename'][7:9]
        file_day_seq = temp_dict['standard_filename'][9:18]
        primaryHDU = fits.PrimaryHDU(global_shared_array[:, :, 0: config.height_Ha]
                                     .reshape((GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT, config.height_Ha)))
        greyHDU = fits.HDUList([primaryHDU])
        greyHDU.writeto(config.save_dir_path + 'RSM' + file_year + '-' + file_mon + '-' + file_day_seq + '_' + str(
            temp_dict['scan_index']) + '_HA.fits')
        greyHDU.close()
        print('生成FE文件中...')
        primaryHDU = fits.PrimaryHDU(global_shared_array[:, :, config.height_Ha:]
                                     .reshape((GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT, config.height_Fe)))
        greyHDU = fits.HDUList([primaryHDU])
        greyHDU.writeto(config.save_dir_path + 'RSM' + file_year + '-' + file_mon + '-' + file_day_seq + '_' + str(
            temp_dict['scan_index']) + '_FE.fits')
        greyHDU.close()

    time_end = time.time()
    print('\n并行进度已完成，所花费时间为：', (time_end - time_start) / 60, 'min(分钟)')
    print('程序结束！')


if __name__ == "__main__":
    main()
