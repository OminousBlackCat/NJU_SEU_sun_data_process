## 使用的Python版本：py3.10.1

本程序依赖的Python母程序在用户文件夹内，已经配置好，无需重新配置

## 目前使用的第三方库:

* numpy

* matplotlib

* astropy （用于读入FITS格式的图片）

* scipy（用于中值滤波）

* jplephem (用于读取星表)

## 安装与使用说明书 

0. 打开bash(shell/Terminal), 并将压缩包复制到某文件夹内
1.  使用` tar -xvf seu_solar_data_process`指令将源码解压
2.  使用` cd seu_solar_data_process`切换目录至解压后的源码路径内
3.  使用` vim config.py`指令修改文件运行的参数
4. 使用` source run_solar`指令运行普通模式
5. 使用` source run_scan_solar`运行摆扫模式
   + 注: **一定要确保source指令在程序源码根目录下运行**

## 在运行程序前需检查并更改的参数(config.py)

+ ` data_dir_path`	此参数代表程序本轮运行时处理的数据所在文件夹	
+ `save_dir_path`    此参数代表本轮程序运行后处理结果所在的文件夹
+ `bin_count`            此参数代表本轮程序运行所处理数据的binning模式
+ ` reversal_mode`    此参数代表在摆扫模式运行下的程序翻转的序列模式
+ 其他参数也可以在程序运行之前检查一次，可以查看config内的描述

## 程序运行需要的外部文件(已经附带在程序根目录data文件夹内)

| 文件名                 | 描述                                               | 备注                                                         |
| :--------------------- | -------------------------------------------------- | ------------------------------------------------------------ |
| for_flat.fits          | binning=1对应的平场文件                            |                                                              |
| for_flat_binning2.fits | binning=2对应的平场文件                            |                                                              |
| dark.fits              | 暗场文件                                           |                                                              |
| HA_absorption.txt      | HA窗口吸收系数文件                                 |                                                              |
| FE_absorption.txt      | FE窗口吸收系数文件                                 |                                                              |
| color_map.txt          | 绘制png日像时使用的色表文件                        |                                                              |
| header.txt             | 标准头部参数文件                                   | 修改header内部的键(key)与评论(comment)可以直接体现在最后生成的文件内 |
| de430.bsp              | 计算B0等值时使用的星表文件(此文件在不同的文件夹内) | 此文件夹使用绝对路径索引                                     |

## 程序完整的参数列表
| 参数名                      | 描述                                               | 默认值 |
| :-------------------------- | -------------------------------------------------- | ------ |
| data_dir_path               | binning=1对应的平场文件                            |        |
| save_dir_path               | binning=2对应的平场文件                            |        |
| bin_count                   | 暗场文件                                           |        |
| sum_dir_path                | HA窗口吸收系数文件                                 |        |
| dark_fits_name              | FE窗口吸收系数文件                                 |        |
| flat_fits_name_bin_1        | 绘制png日像时使用的色表文件                        |        |
| flat_fits_name_bin_2        | 标准头部参数文件                                   |        |
| HA_absorption_path          | 计算B0等值时使用的星表文件(此文件在不同的文件夹内) |        |
| FE_absorption_path          |                                                    |        |
| color_camp_name             |                                                    |        |
| de_file_url                 |                                                    |        |
| header_file                 |                                                    |        |
| multiprocess_count          |                                                    |        |
| sun_row_count_bin_1         |                                                    |        |
| standard_offset_index_bin_1 |                                                    |        |
| curve_cor_x0_bin_1          |                                                    |        |
| curve_cor_C_bin_1           |                                                    |        |
| wavelength_resolution_bin_1 |                                                    |        |
| sum_row_index_HA_bin_1      |                                                    |        |
| sum_row_index_FE_bin_1      |                                                    |        |
| sun_row_count_bin_2         |                                                    |        |
| standard_offset_index_bin_2 |                                                    |        |
| curve_cor_x0_bin_2          |                                                    |        |
| curve_cor_C_bin_2           |                                                    |        |
| wavelength_resolution_bin_2 |                                                    |        |
| sum_row_index_HA_bin_2      |                                                    |        |
| sum_row_index_FE_bin_2      |                                                    |        |
| HA_start                    |                                                    |        |
| FE_start                    |                                                    |        |
