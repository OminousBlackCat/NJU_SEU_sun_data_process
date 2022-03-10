# 存储一些程序需要的必要参数


# 读入数据存放的文件夹
data_dir_path = "data/raw_data"

# 数据输出的存放文件夹
save_dir_path = "data/out/"

# 汇总结果存放的文件夹
sum_dir_path = 'data/sum/'

# 读入的暗场fits文件名, 和py文件需在同级文件夹下
dark_fits_name = 'data/dark.fits'

# 读入的平场参数文件名, 和py文件需在同级文件夹下
flat_fits_name = 'data/for_flat.fits'

# 读入的标准光谱数据文件名, 和py文件需在同级文件夹下
sun_std_name = 'data/bass2000.txt'

# 读入的色谱文件名, 和py文件需在同级文件夹下
color_camp_name = 'data/color_map.txt'

# 读入的日心数据文件名，作为矫正平场的基准文件名
standard_offset_name = 'data/raw_data/RSM20211222T060129-0008-2313.fts'

# 多核并行数
# 若为default(string) 则程序自动根据获得的cpu核数设定并行数
multiprocess_count = 8

# 谱线弯曲矫正参数
# ******此下项目参数在运行程序前应确认是否无误******
# 谱线弯曲矫正x0 参数
curve_cor_x0 = 2321.26

# 谱线弯曲矫正C 参数
curve_cor_C = 1.92909e-011

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
sum_row_index = 132

# 输出何种类型的图片
# 'fts'输出fts格式文件, 'default'则输出以color map为camp的png图片
save_img_form = 'default'
