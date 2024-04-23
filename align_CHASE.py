# Written by Shihao Rao, Nanjing University


# 1. import modules

import sim
import csv
import copy
import numpy as np
from math import *
import glob, cv2, gc, time # gc这个包大概是想用它来清除缓存（一开始是在个人PC上测试程序所以需要清除），有更好的方法可以更改
from astropy.io import fits
from scipy import interpolate
import matplotlib.pyplot as plt



# 2. global information define

Ddate = '20230803'
Halist = sorted(glob.glob('/data/chase/Chase/Lev1/2022/8/3/*T02*HA.fits'))
Halist2 = sorted(glob.glob('/data/chase/Chase/Lev1/2022/8/3/*T03*HA.fits'))
Halist.extend(Halist2)
Felist = sorted(glob.glob('/data/chase/Chase/Lev1/2022/8/3/*T02*FE.fits'))
Felist2 = sorted(glob.glob('/data/chase/Chase/Lev1/2022/8/3/*T03*FE.fits'))
Felist.extend(Felist2)

se0000_bias = sorted(glob.glob('/data/home/wangxinyu/data/jitter_corr/bias/20220803T021402.csv')) # 读取图像偏差
bin_num = 2 # 可以设置是否为binning模式的扫描
spec_win = 'HA' # 可以设置需要做畸变矫正的窗口，两个窗口（HA和FE）畸变矫正的方法略有区别
cut_range = [550, 850, 1550, 1850] # bottom, top, left, right # 切取小块区域用于检验畸变矫正方法是否可行
L1_5_path = '/data/home/wangxinyu/data/Lev1_5/' # 最终存储数据的路径
tot = len(Halist)
winsize = 32 # 畸变窗口大小,值越大，对畸变的改正效果越差，但太阳本身的运动的影响越小

headerkeys = ['BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2', 'NAXIS3', 'TELESCOP', \
              'BIN', 'DATE_OBS', 'CRPIX1', 'CUNIT1', 'CDELT1', \
              'CRVAL1', 'CRPIX2', 'CUNIT2', 'CDELT2', 'CRVAL2', 'CRPIX3', \
              'CUNIT3', 'CDELT3', 'CRVAL3', 'STEPTIME', \
              'PRODATE', 'STR_TIME', 'END_TIME', 'EXP_TIME', 'ORID', 'FRM_NUM', \
              'SPECLINE', 'WAVE_LEN', 'SAT_POS1', 'SAT_POS2', \
              'SAT_POS3', 'SAT_VEL1', 'SAT_VEL2', 'SAT_VEL3', 'INST_ROT', \
              'B0', 'OBS_MOD', 'LVL_NUM'] # 因为之前是需要重新写一个fits文件，所以手动设置了头文件，这一点需要修改一下



# 3. processing series_0000
se0000_hacore = fits.open(Halist[0])[1].data[68, :, :]
se0000_header = fits.open(Halist[0])[1].header
se0000_Xo0 = se0000_header['CRPIX1'] # 读取日面中心（o代表原始original数据定出来的坐标）
se0000_Yo0 = se0000_header['CRPIX2'] # 读取日面中心
h, w = se0000_hacore.shape # 读取图片高度和宽度

se0000_hadata = fits.open(Halist[0])[1].data
haa1, haa2, haa3 = se0000_hadata.shape[0], se0000_hadata.shape[1], se0000_hadata.shape[2] # 读取Ha图片尺寸
del se0000_hadata
gc.collect()

se0000_hawing = fits.open(Halist[0])[1].data[110, :, :]
se0000_strtime = se0000_header['STR_TIME'] # 读取时间
se0000_rx = se0000_header['SAT_POS1'] # 读取卫星相对地球位置
se0000_ry = se0000_header['SAT_POS2'] # 读取卫星相对地球位置
se0000_rz = se0000_header['SAT_POS3'] # 读取卫星相对地球位置
se0000_satpos = [se0000_rx, se0000_ry, se0000_rz]

se0000_center = sim.circle_center(se0000_hawing) # 重新计算日心坐标（此时的日面还存在畸变，拟合得到的日心坐标不准确，后续会计算更准确的日心坐标）
se0000_RSUN = sim.theory_rsun(se0000_strtime, se0000_satpos, bin_num) # 计算出理论的太阳视半径（前面pipeline不知道有没有用这个方法计算日面半径，新的pipeline日面半径以这个为准）
se0000_centerx, se0000_centerz = se0000_center[0], se0000_center[1]

l_se0000_biasx = [] # 读取卫星指向偏差（此时是角秒单位）
l_se0000_biasz = [] # 读取卫星指向偏差（此时是角秒单位）
with open(se0000_bias[0]) as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        l_se0000_biasx.append(float(row['biasx']) / (0.5218 * 2)) # 从角秒单位变为pixel单位
        l_se0000_biasz.append(float(row['biasz']) / (0.5218 * 2)) # 从角秒单位变为pixel单位

se0000_x, se0000_z = np.zeros([haa2, haa3]), np.zeros([haa2, haa3])
lse0000_za, lse0000_zb, lse0000_zc = [], [], []

for j in range(haa2):
    for k in range(haa3):
        se0000_x[j, k] = k - l_se0000_biasx[j] # 拍摄时，若镜头向右偏移，重映射时，采用相应图像向右平移
        if k == 0:
            lse0000_za.append(j) # 完整的z方向2312帧坐标
            if j + l_se0000_biasz[j] not in lse0000_zb: # 新图像j + l_se0000_biasz[j]位置由j位置替换
                lse0000_zb.append(j + l_se0000_biasz[j]) # 因为不均匀，故不涵盖所有2312个z方向坐标，该位置作为插值的自变量
                lse0000_zc.append(j) # 该位置作为插值的因变量

se0000_maxtestz = max(lse0000_zb)
se0000_mintestz = min(lse0000_zb)
se0000_tck = interpolate.splrep(lse0000_zb, lse0000_zc)
lse0000_z = interpolate.splev(lse0000_za, se0000_tck, der = 0) # 插值
lse0000_z = list(lse0000_z)

se0000_shiftmask = np.ones([haa2, haa3])
for j in range(haa2):
    for k in range(haa3):
        se0000_z[j, k] = lse0000_z[j] # 插值之后，每个点都对应了
        if j > se0000_maxtestz or j < se0000_mintestz:
            se0000_shiftmask[j, k] = 0 # 外插的点不可信，需要擦除

se0000_x = se0000_x.astype('float32') # 改变数据类型，否则remap时候会报错（下同）
se0000_z = se0000_z.astype('float32')

se0000_hadata = fits.open(Halist[0])[1].data
se0000_haheader = fits.open(Halist[0])[1].header
se0000_hacubeout = np.zeros([haa1, haa2, haa3])

for wl in range(haa1):
    se0000_im_wl = se0000_hadata[wl, :, :]
    se0000_imout = cv2.remap(se0000_im_wl, se0000_x, se0000_z, borderMode = cv2.BORDER_CONSTANT, \
                        interpolation = cv2.INTER_LINEAR) # 非刚性位移改正畸变，边界值默认为常数0以免假信号造成误解
    se0000_imout = np.multiply(se0000_imout, se0000_shiftmask)
    se0000_hacubeout[wl, :, :] = se0000_imout

# calculate the disk center more precisely
se0000_hacore0 = se0000_hacubeout[68, :, :]
se0000_hawing0 = se0000_hacubeout[110, :, :]
se0000_center = sim.circle_center(se0000_hawing0)
se0000_centerX0, se0000_centerY0 = se0000_center[0], se0000_center[1]



# 4. processing other images

for i in range(tot):
    if i == 0:
        print('Starting to process the first scan')
        
        if spec_win == 'HA':
            se0000_header = fits.open(Halist[i])[1].header
            haa1, haa2, haa3 = fits.open(Halist[i])[1].data.shape
        
        if spec_win == 'FE':
            se0000_header = fits.open(Felist[i])[1].header
            fea1, fea2, fea3 = fits.open(Felist[i])[1].data.shape
            se0000_fedata = fits.open(Felist[i])[1].data
            fea1, fea2, fea3 = fits.open(Felist[i])[1].data.shape
        
        p0 = se0000_header['INST_ROT']
        se0000_strtime = se0000_header['STR_TIME']
        se0000_rx = se0000_header['SAT_POS1']
        se0000_ry = se0000_header['SAT_POS2']
        se0000_rz = se0000_header['SAT_POS3']
        se0000_satpos = [se0000_rx, se0000_ry, se0000_rz]
        se00xx_RSUN = sim.theory_rsun(se0000_strtime, se0000_satpos, bin_num)

        if spec_win == 'FE': # Fe线窗口的第一幅全日面拼接像，使用上一步的畸变矩阵进行矫正
            se00xx_cubeout = np.zeros([fea1, fea2, fea3])
            for wl in range(fea1):
                se0000_im_wl = se0000_fedata[wl, :, :]
                se0000_feimout = cv2.remap(se0000_im_wl, se0000_x, se0000_z, borderMode = cv2.BORDER_CONSTANT, \
                                    interpolation = cv2.INTER_LINEAR) #非刚性位移改正畸变
                se00xx_feimout = np.multiply(se0000_feimout, se0000_shiftmask)
                se00xx_cubeout[wl, :, :] = se00xx_feimout
            se00xx_imwing = se00xx_cubeout[10, :, :]
            
            gc.collect()

        if spec_win == 'HA': # Ha线窗口的第一幅全日面拼接像，由于已经用畸变矩阵矫正过
            se00xx_cubeout = copy.deepcopy(se0000_hacubeout)
            se00xx_imwing = se00xx_cubeout[110, :, :]
        
        hacore0 = se0000_hacubeout[68, :, :] # 取出矫正过Ha线线心图像
        
    else:
        print('Starting to process scan: ', i)
        
        if spec_win == 'HA': # 读取光球图像用于定日心坐标
            se00xx_wing = fits.open(Halist[i])[1].data[110, :, :]
            se00xx_header = fits.open(Halist[i])[1].header
        if spec_win == 'FE':
            se00xx_wing = fits.open(Felist[i])[1].data[10, :, :]
            se00xx_header = fits.open(Felist[i])[1].header
        
        se00xx_strtime = se00xx_header['STR_TIME']
        se00xx_rx = se00xx_header['SAT_POS1']
        se00xx_ry = se00xx_header['SAT_POS2']
        se00xx_rz = se00xx_header['SAT_POS3']
        se00xx_satpos = [se00xx_rx, se00xx_ry, se00xx_rz]
        p0 = se00xx_header['INST_ROT']

        se00xx_center = sim.circle_center(se00xx_wing)
        se00xx_RSUN = sim.theory_rsun(se00xx_strtime, se00xx_satpos, bin_num)
        se00xx_centerx, se00xx_centery = se00xx_center[0], se00xx_center[1]
        se00xx_hacore = fits.open(Halist[i])[1].data[68, :, :].astype('float32')
        dx, dy = -(se00xx_centerx - se0000_centerX0), -(se00xx_centery - se0000_centerY0)
        print(dx, dy, se00xx_RSUN)
        se00xx_hacore2 = sim.imshift(se00xx_hacore, [int(dy), int(dx)]) #刚性对齐
        
        print('Processing distortion correction...')
        
        se00xx_flow = cv2.calcOpticalFlowFarneback(hacore0, se00xx_hacore2, flow = None, pyr_scale = 0.5, \
                                            levels = 3, winsize = winsize, iterations = 5, poly_n = 5, \
                                            poly_sigma = 1.2, flags = 0) #计算畸变

        se00xx_x1, se00xx_y1 = np.meshgrid(np.arange(w), np.arange(h))
        
        if spec_win == 'HA':
            se00xx_data = fits.open(Halist[i])[1].data
            se00xx_cubeout = np.zeros([haa1, haa2, haa3])
            a1, a2, a3 = haa1, haa2, haa3
        if spec_win == 'FE':
            se00xx_data = fits.open(Felist[i])[1].data
            se00xx_cubeout = np.zeros([fea1, fea2, fea3])
            a1, a2, a3 = fea1, fea2, fea3

        se00xx_x2 = copy.deepcopy(se00xx_x1)
        se00xx_y2 = copy.deepcopy(se00xx_y1)
        se00xx_flow2 = copy.deepcopy(se00xx_flow)
        
        for j in range(a2): # 只选取日面范围的区域（2.1日面半径的正方形）计算每步扫描的指向偏差，或许可以自适应调节区域大小
            if abs(j - se00xx_centery) <= 1.05 * se00xx_RSUN:
                xrange = [se00xx_centerx - sqrt((1.05 * se00xx_RSUN) ** 2 - (se00xx_centery - j) ** 2), \
                          se00xx_centerx + sqrt((1.05 * se00xx_RSUN) ** 2 - (se00xx_centery - j) ** 2)]
                flowx_mean = se00xx_flow[j, int(xrange[0]):int(xrange[1]), 0].mean()
                flowy_mean = se00xx_flow[j, int(xrange[0]):int(xrange[1]), 1].mean()
                for k in range(a3):
                    se00xx_flow2[j, k, 0] = flowx_mean # 同一步的图像需要刚性平移改正畸变
                    se00xx_flow2[j, k, 1] = flowy_mean
        
        se00xx_x2 = se00xx_x2 + se00xx_flow2[:, :, 0]
        se00xx_y2 = se00xx_y2 + se00xx_flow2[:, :, 1]
        se00xx_x2 = se00xx_x2.astype('float32')
        se00xx_y2 = se00xx_y2.astype('float32')

        for j in range(a1):
            if spec_win == 'HA':
                se00xx_im3tmp = se00xx_data[j, :, :]
                se00xx_im3 = sim.imshift(se00xx_im3tmp, [int(dy), int(dx)]).astype('float32')
                if j == 110:
                    se00xx_imwing = se00xx_data[j, :, :].astype('float32')
            if spec_win == 'FE':
                se00xx_im3tmp = se00xx_data[j, :, :]
                se00xx_im3 = sim.imshift(se00xx_im3tmp, [int(dy), int(dx)]).astype('float32')
                if j == 10:
                    se00xx_imwing = se00xx_data[j, :, :].astype('float32')
            se00xx_imout = cv2.remap(se00xx_im3, se00xx_x2, se00xx_y2, borderMode = cv2.BORDER_CONSTANT, \
                            interpolation = cv2.INTER_LINEAR) #非刚性位移改正畸变
            se00xx_cubeout[j, :, :] = se00xx_imout
        
        print('Distortion correction successful')
        
    # 将INST_ROT旋转回来
    se00xx_center = sim.circle_center(se00xx_imwing)
    se00xx_centerx, se00xx_centery = se00xx_center[0], se00xx_center[1]
    print(se00xx_centerx, se00xx_centery, se00xx_RSUN)

    print('Rotating...')
    if spec_win == 'HA':
        a1, a2, a3 = haa1, haa2, haa3
    if spec_win == 'FE':
        a1, a2, a3 = fea1, fea2, fea3
    se00xx_cubeout2 = np.zeros([a1, a2, a3])
    se00xx_rotate1 = np.zeros([a2, a3])
    se00xx_rotate2 = np.zeros([a2, a3])

    for j in range(a2):
        for k in range(a3):
            original_coor = [k - se00xx_centerx, 0, j - se00xx_centery]
            local_rot = np.dot(sim.rotation_matrix3('y', p0 / 180 * pi), \
                            original_coor)
            se00xx_rotate1[j, k] = local_rot[0] + se00xx_centerx
            se00xx_rotate2[j, k] = local_rot[2] + se00xx_centery

    se00xx_rotate1 = se00xx_rotate1.astype('float32')
    se00xx_rotate2 = se00xx_rotate2.astype('float32')

    for j in range(a1):
        se00xx_cubeout2[j, :, :] = cv2.remap(se00xx_cubeout[j, :, :], se00xx_rotate1, se00xx_rotate2, \
                                    borderMode = cv2.BORDER_CONSTANT, \
                                        interpolation = cv2.INTER_LINEAR)

    print('Rotation successful')
    
    se00xx_mask = np.ones([a2, a3])
    for j in range(a2):
        for k in range(a3):
            if (j - se00xx_centery) ** 2 + (k - se00xx_centerx) ** 2 >= (1.1 * se00xx_RSUN) ** 2:
                se00xx_mask[j, k] = 0
    if spec_win == 'HA':
        se00xx_wing2 = copy.deepcopy(np.multiply(se00xx_cubeout2[110, :, :], se00xx_mask))
    if spec_win == 'FE':
        se00xx_wing2 = copy.deepcopy(np.multiply(se00xx_cubeout2[10, :, :], se00xx_mask))
    
    se00xx_center = sim.circle_center(se00xx_wing2)
    se00xx_centerx, se00xx_centery = se00xx_center[0], se00xx_center[1]

    print('Final disk parameter for series {}: '.format(i), se00xx_centerx, se00xx_centery, se00xx_RSUN)

    # 绘制快速预览图
    if spec_win == 'HA':
        core_img = se00xx_cubeout2[68]
        wing_img = se00xx_cubeout2[110]
        filename = Halist[i].split('/')[-1]
    if spec_win == 'FE':
        core_img = se00xx_cubeout2[22]
        wing_img = se00xx_cubeout2[10]
        filename = Felist[i].split('/')[-1]

    # 线心全日面像
    plt.figure(figsize = (20, 20))
    plt.imshow(core_img[int(se00xx_centery) - 1000:int(se00xx_centery) + 1000, int(se00xx_centerx) - 1000:int(se00xx_centerx) + 1000], origin = 'lower', cmap = 'afmhot', vmin = 0, vmax = 4 * core_img.mean())
    plt.title(filename[3:7] + '-' + filename[7:9] + '-' + filename[9:11] + 'UT' + filename[12:14] + ':' + filename[14:16] + ':' + filename[16:18], size = 30)
    if spec_win == 'HA':
        plt.savefig(L1_5_path + 'figure/ha/' + filename.split('.')[0] + '.png', dpi = 300)
    if spec_win == 'FE':
        plt.savefig(L1_5_path + 'figure/fe/' + filename.split('.')[0] + '.png', dpi = 300)
    plt.close()

    # 线心局部区域像
    plt.figure(figsize = (20, 20))
    plt.imshow(core_img[int(se00xx_centery) - 1000:int(se00xx_centery) + 1000, int(se00xx_centerx) - 1000:int(se00xx_centerx) + 1000][cut_range[0]:cut_range[1], cut_range[2]:cut_range[3]], \
    origin = 'lower', cmap = 'afmhot', vmin = 0, vmax = 4 * core_img.mean())
    plt.title(filename[3:7] + '-' + filename[7:9] + '-' + filename[9:11] + 'UT' + filename[12:14] + ':' + filename[14:16] + ':' + filename[16:18], size = 30)
    if spec_win == 'HA':
        plt.savefig(L1_5_path + 'figure/ha/cut_core' + filename[3:].split('.')[0] + '.png', dpi = 300)
    if spec_win == 'FE':
        plt.savefig(L1_5_path + 'figure/fe/cut_core' + filename[3:].split('.')[0] + '.png', dpi = 300)
    plt.close()    

    # 线翼局部区域像
    plt.figure(figsize = (20, 20))
    plt.imshow(wing_img[int(se00xx_centery) - 1000:int(se00xx_centery) + 1000, int(se00xx_centerx) - 1000:int(se00xx_centerx) + 1000][cut_range[0]:cut_range[1], cut_range[2]:cut_range[3]], \
    origin = 'lower', cmap = 'afmhot', vmin = 0, vmax = 4 * wing_img.mean())
    plt.title(filename[3:7] + '-' + filename[7:9] + '-' + filename[9:11] + 'UT' + filename[12:14] + ':' + filename[14:16] + ':' + filename[16:18], size = 30)
    if spec_win == 'HA':
        plt.savefig(L1_5_path + 'figure/ha/cut_wing' + filename[3:].split('.')[0] + '.png', dpi = 300)
    if spec_win == 'FE':
        plt.savefig(L1_5_path + 'figure/fe/cut_wing' + filename[3:].split('.')[0] + '.png', dpi = 300)
    plt.close()       

    # 存为fits文件
    dataf = fits.PrimaryHDU(se00xx_cubeout2)
    dataHDU = fits.HDUList([dataf])
    hdr = dataHDU[0].header
    
    for j in range(39):
        try:
            hdr[headerkeys[j]] = se00xx_header[headerkeys[j]]
        except:
            try:
                hdr[headerkeys[j]] = se0000_header[headerkeys[j]]
            except:
                print(headerkeys[j] + 'not exist.')
    hdr['TELESCOP'] = 'CHASE'
    hdr['SPECLINE'] = 'HA'
    hdr['OBS_MOD'] = 'RSM'
    hdr['CRPIX1'] = round(se00xx_centerx, 6)
    hdr['CRPIX2'] = round(se00xx_centery, 6)
    hdr['R_SUN'] = round(se00xx_RSUN, 6)
    hdr['RSUN_OBS'] = round(se00xx_RSUN * bin_num * 0.5218, 6)

    print('Writing...')
    dataHDU.writeto(L1_5_path + filename, overwrite = True)
    print('Writing successful')














