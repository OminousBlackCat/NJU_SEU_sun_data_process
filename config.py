# 存储一些程序需要的必要参数

# 读入数据存放的文件夹, 相对路径与绝对路径均可
# 文件夹最后的结尾需要带'/'
data_dir_path = "data/raw_data_test/"

# 数据输出的存放文件夹, 相对路径与绝对路径均可
# 请确保此文件夹存在并拥有写入权限
# 文件夹最后的结尾需要带'/'
save_dir_path = "data/out/"

# 日像汇总结果存放的文件夹, 相对路径与绝对路径均可
# 请确保此文件夹存在并拥有写入权限
# 文件夹最后的结尾需要带'/'
sum_dir_path = 'data/sum/'

# 读入的暗场fits文件路径, 相对路径与绝对路径均可
dark_fits_name = 'data/dark.fits'

# 读入的平场参数文件路径, 相对路径与绝对路径均可
# bin = 1
flat_fits_name_bin_1 = 'data/for_flat.fits'
# bin = 2
flat_fits_name_bin_2 = 'data/for_flat_binning2.fits'

# 读入的标准光谱数据文件路径, 相对路径与绝对路径均可
sun_std_name = 'data/bass2000.txt'
HA_absorption_path = 'data/HA_absorption.txt'
FE_absorption_path = 'data/FE_absorption.txt'

# 读入的色谱(作为日像汇总结果png格式的输入)文件路径, 相对路径与绝对路径均可
color_camp_name = 'data/color_map.txt'

# 读入的星图文件
de_file_url = 'file://'

# 读入的HA头部与FE头部参数文件 注意是txt格式 一般不做修改
header_file = 'data/HA_header.txt'

# 多核并行数
# 若为'default' 则程序自动根据获得的cpu核数设定并行数
# 默认使用可支配核数-4的并行数量
multiprocess_count = 'default'

# 是否bin
# 修改此参数的时候记得修改sun_row_index
# TODO: bin_count务必在程序运行前要确认
bin_count = 1


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
# 合并日像时使用的行数
sum_row_index_bin_1 = 132


# bin = 2的时候, 程序将使用的参数
# 一个扫描序列最多包含多少行(一组数据包含多少文件)
sun_row_count_bin_2 = 2313
# 矫正所需的日心文件序号
standard_offset_index_bin_2 = 1156
# 谱线弯曲矫正x0 参数
curve_cor_x0 = 1161.58
# 谱线弯曲矫正C 参数
curve_cor_C = 7.639e-011
# 谱线分辨率
wavelength_resolution_bin_2 = 0.048404602
# 合并日像时使用的行数
sum_row_index_bin_2 = 66


# 通用参数
# 红蓝移修正内 波长相关数据
HA_start = 6559.58  # HA窗口零像素波长
FE_start = 6567.66  # FE窗口零像素波长
K = 0.024202301  # 分辨率
# 写入头文件内的HA与FEI相关数据
HA_lineCore = 6562.82  # HA线心处波长
FE_lineCore = 6569.22  # FE线心处波长

# TODO: 以下参数仅需修改数值即可, 无需在bin修改时更改
# 读入的日心数据(矫正平场时需要)序号(*不要使用引号引起来*)
standard_offset_index = 2313

# 中值滤波窗口大小
filter_kernel_size = 3

# 图像窗口数据 包含ha窗口与fe窗口的长度
height_Ha = 260
height_Fe = 116

# 输出何种类型的图片
# 'fts'输出fts格式灰度文件, 'default'则输出以color map为camp的png图片
save_img_form = 'fts'
