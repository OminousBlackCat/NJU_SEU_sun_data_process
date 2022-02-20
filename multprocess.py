import multiprocessing as mp
import os
import suntools
import time
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits

# cpu核心数量 确定了并行数量
cpu_core_nums = mp.cpu_count()
# 数据文件夹
data_dir_path = "C:\\Users\\seu-wxy\\Desktop\\太阳数据\\6-803\\6-803"
# 存储文件夹
save_dir_path = "C:\\Users\\seu-wxy\\Desktop\\太阳数据\\save\\"


# 读取数据文件夹所有文件
def read_fits_directory():
    arr = os.listdir(data_dir_path)
    return arr


# 定义target task
# 传入一个文件名，读取此文件名对应的fits文件并对其做曲线矫正
def target_task(filename):
    filePath = data_dir_path + "\\" + filename
    image_data = fits.getdata(get_pkg_data_filename(filePath))
    image_data = suntools.Curve_correction(image_data, 2225, 0.06 / (2225 - 770) / (2225 - 770))
    plt.imsave(save_dir_path + filename + "result.jpg", image_data)


def main():
    time_start = time.time()
    lst = read_fits_directory()
    pool = mp.Pool(processes=12)
    pool.map(target_task, lst)
    time_end = time.time()
    print('time cost', time_end - time_start, 's')
    # print('预计时间', (time_end - time_start) / now * N, 's')


if __name__ == "__main__":
    main()
