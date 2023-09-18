"""
存储一些程序需要的必要参数
包含TODO标志的参数较为重要 请在运行前检查参数

@author: seu_wxy
"""

# 读入数据存放的文件夹, 相对路径与绝对路径均可
# 文件夹最后的结尾需要带'/'
# TODO: 在运行程序前一定要修改此目录, 此目录应精确到当日的扫描序列文件夹
# 例: /data/chase/Chase/2021/12/16-428/
data_dir_path = "/data/chase/Chase/2023/8/25-010310-3/"
# data_dir_path = "/data/chase/Chase/2023/5/11-8709-2/"
# data_dir_path = "/data/chase/Chase/2023/8/8-010051-2/"


# 数据输出的存放文件夹, 相对路径与绝对路径均可, 需精确到天数
# 请确保此文件夹存在并拥有写入权限ls
# **** 请在程序运行前创建对应月份的文件夹! ****
# 文件夹最后的结尾需要带'/'
# 例: /data/chase/Chase/Lev1/2021/12/16/
# TODO: 在程序运行前需要修改此目录
bias_dir_path = "/data/chase/cyctest/Chase/Lev1/pointing_bias/8/8-010051-2/"



# TODO: 在程序运行前需要修改此目录
# save_dir_path = "/data/chase/cyctest/Chase/Lev1/2023/8/8-010051-2/"
save_dir_path = "/data/chase/cyctest/Chase/Lev1/2022/8/4445/"

# 存储视频的文件夹, 相对路径与绝对路径均可
# 请确保此文件夹存在并拥有写入权限
# 可以不需要手动创建文件夹  会自动创建文件夹
# 文件夹最后的结尾需要带'/'
# TODO: 每个月需要修改此目录
# video_dir_path = "/data/chase/cyctest/Chase/Lev1/2022/8/7-010046-1/video"
# video_dir_path = "/data/chase/cyctest/Chase/Lev1/2023/8/8-010051-2/"
video_dir_path = "/data/chase/cyctest/Chase/Lev1/2022/8/4445/video/"



# 是否bin
# 修改此参数的时候记得修改sun_row_index
# TODO: bin_count务必在程序运行前要确认
bin_count = 2

# 摆扫模式时翻转的序列 序号标号模式为: 0, 1, 2, 3, ...
# 'odd'表示奇数序列翻转 'even'表示偶数序列翻转 'none'表示不翻转
# TODO: 如果运行在摆扫模式 请务必检查此参数
reversal_mode = 'odd'

# 日像汇总结果存放的文件夹, 相对路径与绝对路径均可
# 请确保此文件夹存在并拥有写入权限
# 文件夹最后的结尾需要带'/'
# 目前直接使用与data_dir_path相同的路径
sum_dir_path = save_dir_path

# 读入的暗场fits文件路径, 相对路径与绝对路径均可
# bin = 1的时候使用的暗场
dark_fits_name = 'data/dark.fits'
# bin = 2的时候使用的暗场
dark_fits_name_bin_2 = 'data/dark_binning2.fits'

# 读入的平场参数文件路径, 相对路径与绝对路径均可
# bin = 1的时候
flat_fits_name_bin_1 = 'data/for_flat.fits'
# bin = 2的时候
flat_fits_name_bin_2 = 'data/for_flat_binning2.fits'

# 读入的标准光谱数据文件路径, 相对路径与绝对路径均可
# sun_std_name可以删除
HA_absorption_path = 'data/HA_absorption.txt'
FE_absorption_path = 'data/FE_absorption.txt'

# 读入的色谱(作为日像汇总结果png格式的输入)文件路径, 相对路径与绝对路径均可
color_camp_name = 'data/color_map.txt'

# 读入的星图文件, 需要以url的格式输入, 因此为绝对路径
# 此路径不要随意改变, 若需要改变, 则首先确保绝对路径内拥有此星图文件
de_file_url = 'file:///data/home/wangxinyu/de430.bsp'

# 读入的HA头部与FE头部参数文件 注意是txt格式 一般不做修改
# 如需增添新的头部 可以修改对应的header.txt, 但具体写入何值则需要修改header.py
header_file = 'data/header.txt'

# 多核并行数
# 若为'default' 则程序自动根据获得的cpu核数设定并行数
# 默认使用可支配核数-4的并行数量
multiprocess_count = 'default'

# sit stare模式
# 若为Ture, 则不对序列进行数额判断, 所有序列都会尝试进行处理
sit_stare_mode = False

# bin = 1的时候, 程序将使用的参数
# 一个扫描序列最多包含多少行(一组数据包含多少文件)
sun_row_count_bin_1 = 4625
# 矫正所需的日心文件序号
standard_offset_index_bin_1 = 2313
# 谱线弯曲矫正x0 参数
curve_cor_x0_bin_1 = 2321.26
# 谱线弯曲矫正C 参数
curve_cor_C_bin_1 = 1.92909e-011
# 谱线分辨率
wavelength_resolution_bin_1 = 0.024202301
# 合并HA日像时使用的行数
sum_row_index_HA_bin_1 = 136
# 合并FE日像时使用的行数
sum_row_index_FE_bin_1 = 20

# bin = 2的时候, 程序将使用的参数
# 一个扫描序列最多包含多少行(一组数据包含多少文件)
sun_row_count_bin_2 = 2313
# 矫正所需的日心文件序号
standard_offset_index_bin_2 = 1156
# 谱线弯曲矫正x0 参数
curve_cor_x0_bin_2 = 1161.58
# 谱线弯曲矫正C 参数
curve_cor_C_bin_2 = 7.639e-011
# 谱线分辨率
wavelength_resolution_bin_2 = 0.048404602
# 合并HA日像时使用的行数
sum_row_index_HA_bin_2 = 68
# 合并FE日像时使用的行数
sum_row_index_FE_bin_2 = 10

# 通用参数
# 红蓝移修正内 波长相关数据
HA_start = 6559.58  # HA窗口零像素波长
FE_start = 6567.66  # FE窗口零像素波长
# 写入头文件内的HA与FEI相关数据
HA_lineCore = 6562.82  # HA线心处波长
FE_lineCore = 6569.22  # FE线心处波长
# 中值滤波窗口大小
filter_kernel_size = 3
# 图像窗口数据 包含ha窗口与fe窗口的长度
height_Ha = 260
height_Fe = 116
# 通用像素分辨率
pixel_resolution = 0.52
# 右侧多少行的数据会被置0
pixel_to_zero_right_count = 100
# 左侧多少行的数据会被置0
pixel_to_zero_left_count = 100
# SIT_MODE下 预设数组大小(bin=2的情况下会默认除以2)
sit_stare_array_size = 2400
# 输出图像日期标记的字体大小与字体厚度
# 使用数字输入
date_font_size = 4
date_font_thick = 2
# 插值参数，使用几次插值
interpolation_parameter = 3
# 生成视频的帧率相关
# 注意: 一定为正整数
# 视频帧率(fps)(一秒钟画面显示多少画面)
frame_pre_sec = 10
# 单张图片的写入次数
# 举例: 如果下述参数值为10, 视频帧率也为10
# 实际效果为: 视频每秒10张画面, 且这十张画面均为单张HA图像, 也就是一张图展示1s
write_to_video_count = 1
# 是否需要跳帧
# 举例: 如果需要 [每隔n个] 观测序列写入一次图片 就将此值置为n
# 注意: 一定为整数
# 默认为0(没有间隔, 所有图片都写入视频)
write_to_video_bias = 0

# 摆扫模式参数
# 计算日心位置时的均值阈值
center_mean_threshold = 150
# 使用的线心位置
# 目前使用Fe窗口的第10行作为获取的线心初始 往下4行求平均
center_mean_index = int(height_Ha / 2) + 9
center_mean_count = 4

# 序列扫描时间偏移量(单位: 秒)
# 会在写入头部的时候将对应的STR_TIME/OBS_TIME/END_TIME时间加上时间偏移量
scan_time_offset = 0

# 输出何种类型的图片
# 'fts'输出fts格式灰度文件, 'default'则输出以color map为camp的png图片
save_img_form = 'default'

# 质心计算依赖数值
Ha_lower = 94
Ha_Upper = 172
Fe_lower = 38
Fe_Upper = 84

winsize = 32 # 畸变窗口大小,值越大，对畸变的改正效果越差，但太阳本身的运动的影响越小