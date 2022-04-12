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
            if int(filename[24:28]) == 1:
                global_multiprocess_list[j]['first_filename'] = filename
            if int(filename[24:28]) > int(global_multiprocess_list[j]['last_filename'][24:28]):
                global_multiprocess_list[j]['last_filename'] = filename
            break
    if not ifFind:
        global_multiprocess_list.append({
            'scan_index': temp_index,  # 扫描序号
            'file_list': [],  # 文件名列表
            'file_count': 0,  # 包含的文件数
            'standard_filename': '',  # 标准日心文件名
            'first_filename': '',  # 序列开始文件
            'last_filename': '',  # 序列结束文件
            'flat_data': None,  # 此序列的校正后平场数据
            'abortion_data': None,  # 此序列的校正后红蓝移数据
            'header': fits.header.Header()  # 此序列的头部, 构造了一个新的header
        })
        global_multiprocess_list[len(global_multiprocess_list) - 1]['file_list'].append(filename)
        global_multiprocess_list[len(global_multiprocess_list) - 1]['file_count'] += 1
        if int(filename[24:28]) == config.standard_offset_index:
            global_multiprocess_list[len(global_multiprocess_list) - 1]['standard_filename'] = filename
        if int(filename[24:28]) == 1:
            global_multiprocess_list[len(global_multiprocess_list) - 1]['first_filename'] = filename
        global_multiprocess_list[len(global_multiprocess_list) - 1]['last_filename'] = filename

# 剔除不完整序列
for i in range(len(global_multiprocess_list)):
    if global_multiprocess_list[i]['file_count'] < config.sun_row_count - 5000 or \
            global_multiprocess_list[i]['standard_filename'] is None:
        print('文件夹中包含不完整序列, 序列序号为:' + str(global_multiprocess_list[i]['scan_index']))
        print('本次数据处理将不此序列进行处理.....')
        global_multiprocess_list.remove(global_multiprocess_list[i])

# 读取头部参数文件
# 为每个序列都创建头
global_header_list = header.read_header_from_txt(config.HA_header_file)
# 将读入的list构造成为header
for h in global_header_list:
    for temp_dict in global_multiprocess_list:
        temp_dict['header'].set(h['key'], value=h['value'], comment=h['comment'])
# 将静态的值赋入header
for item in header.static_header_items:
    for temp_dict in global_multiprocess_list:
        temp_dict['header'].set(item['key'], item['value'])

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
    flat_img, standard_HA_width, standard_FE_width = suntools.curve_correction(flat_img - dark_img, config.curve_cor_x0,
                                                                               config.curve_cor_C)
temp_img.close()

# 读取经过日心的图片 作为基准
# 读取标准太阳光谱数据
sun_std = suntools.get_Sunstd(config.sun_std_name)
sample_from_standard = None
HA_LINECORE = None
FEI_LINECORE = None
try:
    for temp_dict in global_multiprocess_list:
        # 对每个序列进行校正
        print('校正扫描序列' + str(temp_dict['scan_index']) + '中...使用标准校正文件为:' + temp_dict['standard_filename'])
        print("校正平场中...")
        standard_name = temp_dict['standard_filename']
        temp_img = fits.open(read_dir + '/' + standard_name)
        standard_header = temp_img[0].header
        # 复制一些值去头部
        HA_LINECORE = float(standard_header['WAVE_LEN'].split('and')[0].strip())
        FEI_LINECORE = float(standard_header['WAVE_LEN'].split('and')[1].strip())
        for item in header.copy_header_items:
            temp_dict['header'].set(item['key'], standard_header[item['key']])
        temp_dict['header'].set('BIN', config.bin_count)
        temp_dict['header'].set('DATE_OBS', standard_header['STR_TIME'])
        standard_img = np.array(temp_img[0].data, dtype=float)
        standard_img = suntools.moveImg(standard_img, -2)
        standard_img, standard_HA_width, standard_FE_width = suntools.curve_correction(standard_img - dark_img,
                                                                                       config.curve_cor_x0,
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
        abortion = suntools.RB_getdata(standard_img, sun_std, standard_HA_width, standard_FE_width)
        temp_dict['flat_data'] = flatTemp
        temp_dict['abortion_data'] = abortion
        temp_img.close()
        print("序列:" + str(int(standard_name[19:23])) + "矫正完成")
        print('计算B0, INST_ROT中....')
        # temp_B0, temp_INST_ROT = suntools.getB0P0(standard_header['Q0'], standard_header['Q1'], standard_header['Q2'],
        #                                           standard_header['Q3'], standard_header['STR_TIME'])
        # temp_dict['header'].set('B0', temp_B0)
        # temp_dict['header'].set('INST_ROT', temp_INST_ROT)
        first_name = temp_dict['first_filename']
        last_name = temp_dict['last_filename']
        temp_dict['header'].set('STR_TIME', first_name[3:7] + '-' + first_name[7:9] + '-' + first_name[9:11] + 'T'
                                + first_name[12:14] + ':' + first_name[14:16] + ':' + first_name[16:18])
        temp_dict['header'].set('END_TIME', last_name[3:7] + '-' + last_name[7:9] + '-' + last_name[9:11] + 'T'
                                + last_name[12:14] + ':' + last_name[14:16] + ':' + last_name[16:18])
        temp_dict['header'].set('FRM_NUM', '1~' + str(temp_dict['file_count']))

except uEr.URLError as error:
    print("Error: 标准日心校准文件未找到, 请检查config文件或存放目录")
    print(error)
    sys.exit("程序终止")
except OSError as error:
    print("Error: 标准日心校准文件读取发生错误, 请检查文件读取权限")
    print(error)
    sys.exit("程序终止")

# 读取输出色谱
color_map = suntools.get_color_map(config.color_camp_name)

# 检查输出文件夹是否存在 不存在则创建
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# 全局进度控制
file_count = mp.Value('i', len(read_fits_directory()))
remaining_count = mp.Value('i', 1)
if_first_print = mp.Value('b', True)

# 全局共享内存
# 需要创建一个四维数组
# 对于最后写入的文件 NAXIS顺序为: 狭缝宽度 文件序号(扫描序号)  波长深度
# 因此np array的shape应为(波长深度, 文件序号, 狭缝宽度)
# 每个单独文件对应的是 xy平面的一个二维数组
# TODO: 是否需要将shape调转顺序？
GLOBAL_ARRAY_X_COUNT = sample_from_standard.shape[0]
GLOBAL_ARRAY_Y_COUNT = config.sun_row_count
GLOBAL_ARRAY_Z_COUNT = sample_from_standard.shape[1]
print('SHAPE:' + str(GLOBAL_ARRAY_X_COUNT) + ',' + str(GLOBAL_ARRAY_Y_COUNT) + ',' + str(GLOBAL_ARRAY_Z_COUNT))
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
    image_data = suntools.MedSmooth(image_data)
    # 转为整型, 并将每行的最后部分置零
    image_data = np.array(image_data, dtype=np.int16)
    global_shared_array = np.frombuffer(GLOBAL_SHARED_MEM.get_obj(), dtype=np.int16)
    global_shared_array = global_shared_array.reshape(GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT)
    global_shared_array[:, int(file_position) - 1, :] = image_data
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
        sum_data = np.zeros((config.sun_row_count, sample_from_standard.shape[1]))
        global_shared_array = np.frombuffer(GLOBAL_SHARED_MEM.get_obj(), dtype=np.int16)
        global_shared_array = global_shared_array.reshape(GLOBAL_ARRAY_X_COUNT, GLOBAL_ARRAY_Y_COUNT,
                                                          GLOBAL_ARRAY_Z_COUNT)
        # 将小于0的值全部赋为0
        global_shared_array[global_shared_array < 0] = 0
        print("SHAPE为：" + str(global_shared_array.shape))
        # 输出太阳像
        for i in range(global_shared_array.shape[1]):
            sum_data[i] = global_shared_array[config.sum_row_index, i, :].reshape(sample_from_standard.shape[1])
        print('计算CCD太阳像半径中...')
        R_x, R_y, radius = suntools.getCircle(sum_data)
        OBS_Radius = radius * temp_dict['header']['CDELT1']
        temp_dict['header'].set('CRPIX1', R_x)
        temp_dict['header'].set('CRPIX2', R_y)
        temp_dict['header'].set('R_SUN', radius)
        temp_dict['header'].set('RSUN_OBS', OBS_Radius)
        if config.save_img_form == 'default':
            # 使用读取的色谱进行输出 imsave函数将自动对data进行归一化
            print('输出序号为' + str(temp_dict['scan_index']) + '的png...')
            plt.imsave(config.sum_dir_path + 'sum' + str(temp_dict['scan_index']) + ".png", sum_data, cmap=color_map)
        if config.save_img_form == 'fts':
            # 不对data进行任何操作 直接输出为fts文件
            print('输出序号为' + str(temp_dict['scan_index']) + '的fits...')
            primaryHDU = fits.PrimaryHDU(sum_data)
            greyHDU = fits.HDUList([primaryHDU])
            greyHDU.writeto(config.sum_dir_path + 'sum' + str(temp_dict['scan_index']) + '.fts', overwrite=True)
            greyHDU.close()
        print('生成HA文件中...')
        temp_dict['header'].set('SPECLINE', 'HA')
        temp_dict['header'].set('LINECORE', HA_LINECORE)
        temp_dict['header'].set('PRODATE', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
        file_year = temp_dict['standard_filename'][3:7]
        file_mon = temp_dict['standard_filename'][7:9]
        file_day_seq = temp_dict['standard_filename'][9:18]
        print(repr(temp_dict['header']))
        primaryHDU = fits.CompImageHDU(global_shared_array[0: standard_HA_width, :, :]
                                       .reshape((standard_HA_width, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT))
                                       , header=temp_dict['header'], compression_type='GZIP_1')
        primaryHDU.header.set('NAXIS', comment='Number of data axes')
        primaryHDU.header.set('NAXIS1', comment='Length of data axis 1 (slit dimension)')
        primaryHDU.header.set('NAXIS2', comment='Length of data axis 2 (scanning steps)')
        primaryHDU.header.set('NAXIS3', comment='Length of data axis 3 (wavelength dimension)')
        primaryHDU.header.set('BZERO', comment='Data is Unsigned Integer')
        primaryHDU.header.set('BSCALE', comment='default scaling factor')
        primaryHDU.header.add_comment('Spectral curvature corrected')
        primaryHDU.header.add_comment('Dark subtracted')
        primaryHDU.header.add_comment('Flat-field corrected')
        primaryHDU.header.add_comment('Processed by RSM_prep')
        primaryHDU.writeto(config.save_dir_path + 'RSM' + file_year + '-' + file_mon + '-' + file_day_seq + '_' + str(
            temp_dict['scan_index']).zfill(4) + '_HA.fits', overwrite=True)
        print('生成FE文件中...')
        # 修改header内的SPECLINE与LINECORE
        temp_dict['header'].set('SPECLINE', 'FEI')
        temp_dict['header'].set('LINECORE', FEI_LINECORE)
        temp_dict['header'].set('PRODATE', datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
        primaryHDU = fits.CompImageHDU(global_shared_array[standard_HA_width:, :, :]
                                       .reshape((standard_FE_width, GLOBAL_ARRAY_Y_COUNT, GLOBAL_ARRAY_Z_COUNT))
                                       , header=temp_dict['header'], compression_type='GZIP_1')
        primaryHDU.header.set('NAXIS', comment='Number of data axes')
        primaryHDU.header.set('NAXIS1', comment='Length of data axis 1 (slit dimension)')
        primaryHDU.header.set('NAXIS2', comment='Length of data axis 2 (scanning steps)')
        primaryHDU.header.set('NAXIS3', comment='Length of data axis 3 (wavelength dimension)')
        primaryHDU.header.set('BZERO', comment='Data is Unsigned Integer')
        primaryHDU.header.set('BSCALE', comment='default scaling factor')
        primaryHDU.header.add_comment('Spectral curvature corrected')
        primaryHDU.header.add_comment('Dark subtracted')
        primaryHDU.header.add_comment('Flat-field corrected')
        primaryHDU.header.add_comment('Processed by RSM_prep')
        primaryHDU.writeto(config.save_dir_path + 'RSM' + file_year + '-' + file_mon + '-' + file_day_seq + '_' + str(
            temp_dict['scan_index']).zfill(4) + '_FE.fits', overwrite=True)
        if_first_print.value = True

    time_end = time.time()
    print('并行进度已完成，所花费时间为：', (time_end - time_start) / 60, 'min(分钟)')
    print('程序结束！')


if __name__ == "__main__":
    main()
