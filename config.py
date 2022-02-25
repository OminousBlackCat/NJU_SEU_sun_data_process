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

# 多核并行数
# 若为default(string) 则程序自动根据获得的cpu核数设定并行数
multiprocess_count = 14

# 谱线弯曲矫正参数
# ******此下项目参数在运行程序前应确认是否无误******
# 谱线弯曲矫正x0 参数
curve_cor_x0 = 2321.26
# 谱线弯曲矫正C 参数
curve_cor_C = 1.92909e-011

# 滤波窗口大小
filter_kernel_size = 3
