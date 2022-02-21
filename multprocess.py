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


# 主函数流程：
# i.读取参数文件 得到曲线矫正参数
# ii.读取标准太阳光谱 获得标准太阳光谱数据
# iii.读取日心中间文件 计算平场
# iv.从头文件开始对每个文件进行去暗场再进行曲线矫正
# v.对上一步得到的数据根据标准太阳光谱数据进行红蓝移矫正
# vi.将得到的数据存入新的fits文件中
# vii.必要时再对每个新文件进行读取得到整体太阳数据
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
