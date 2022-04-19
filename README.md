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
4.  使用` source run_solar`指令运行程序
    + 注: **一定要确保source指令在程序源码根目录下运行**

## 在运行程序前需检查并更改的参数(config.py)

+ ` data_dir_path`	此参数代表程序本轮运行时处理的数据所在文件夹	

+ `save_dir_path`    此参数代表本轮程序运行后处理结果所在的文件夹

+ `bin_count`            此参数代表本轮程序运行所处理数据的binning模式

+ 其他参数也可以在程序运行之前检查一次，可以查看config内的描述

## 程序运行需要的外部文件(已经附带在程序根目录data文件夹内)

| 文件名                 | 描述                                               |
| :--------------------- | -------------------------------------------------- |
| for_flat.fits          | binning=1对应的平场文件                            |
| for_flat_binning2.fits | binning=2对应的平场文件                            |
| dark.fits              | 暗场文件                                           |
| HA_absorption.txt      | HA窗口吸收系数文件                                 |
| FE_absorption.txt      | FE窗口吸收系数文件                                 |
| color_map.txt          | 绘制png日像时使用的色表文件                        |
| header.txt             | 标准头部参数文件                                   |
| de430.bsp              | 计算B0等值时使用的星表文件(此文件在不同的文件夹内) |



