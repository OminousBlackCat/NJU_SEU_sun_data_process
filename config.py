# 存储一些程序需要的必要参数


# 读入数据存放的文件夹, 相对路径与绝对路径均可
data_dir_path = "data/raw_data_test"

# 数据输出的存放文件夹, 相对路径与绝对路径均可
# 请确保此文件夹存在并拥有写入权限
save_dir_path = "data/out/"

# 日像汇总结果存放的文件夹, 相对路径与绝对路径均可
# 请确保此文件夹存在并拥有写入权限
sum_dir_path = 'data/sum/'

# 读入的暗场fits文件路径, 相对路径与绝对路径均可
dark_fits_name = 'data/dark.fits'

# 读入的平场参数文件路径, 相对路径与绝对路径均可
flat_fits_name = 'data/for_flat_binning2.fits'

# 读入的标准光谱数据文件路径, 相对路径与绝对路径均可
sun_std_name = 'data/bass2000.txt'

# 读入的色谱(作为日像汇总结果png格式的输入)文件路径, 相对路径与绝对路径均可
color_camp_name = 'data/color_map.txt'

# 读入的HA头部与FE头部参数文件 注意是txt格式 一般不做修改
HA_header_file = 'data/HA_header'
FE_header_file = 'data/FE_header'

# 多核并行数
# 若为'default' 则程序自动根据获得的cpu核数设定并行数
# 默认使用可支配核数-4的并行数量
multiprocess_count = 'default'

# 一个扫描序列最多包含多少行(一组数据包含多少文件)(一般不需要改动此参数)
sun_row_count = 4625

# 是否bin
# 修改此参数的时候记得修改sun_row_index
bin_count = 1

# 读入的日心数据(矫正平场时需要)序号(*不要使用引号引起来*)
standard_offset_index = 1256

# 谱线弯曲矫正参数
# 谱线弯曲矫正x0 参数
curve_cor_x0 = 2321.26#1161.58#2321.26

# 谱线弯曲矫正C 参数
curve_cor_C = 1.92909e-011#7.639e-011#1.92909e-011

# 中值滤波窗口大小
filter_kernel_size = 3

# 图像窗口数据 包含ha窗口与fe窗口的长度
height_Ha = 260
height_Fe = 116

# 红蓝移修正内 波长相关数据
HA = 6559.5804
FE = 6569.22
K = 0.024202301

# 合并图像时所拼的行数
# 默认为132行（光谱强度最低点）
# 此参数要伴随bin一起更改
sum_row_index = 132

# 输出何种类型的图片
# 'fts'输出fts格式灰度文件, 'default'则输出以color map为camp的png图片
save_img_form = 'default'
