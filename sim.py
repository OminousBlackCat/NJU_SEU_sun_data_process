# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 16:05:28 2018

@author: jkf
pip install opencv-contrib-python==3.4.2.16
"""

#import cv2
import matplotlib.pyplot as plt
import numpy as np
from math import *
from skimage.transform import warp, EuclideanTransform,SimilarityTransform, rotate,rescale,resize
from skimage import transform
import astropy
import jplephem
from astropy.time import Time
from astropy import coordinates
import scipy.signal as signal
from scipy.optimize import curve_fit
import scipy.fftpack as fft

import config


#import cupy as cp


def rotation_matrix3(xyz, theta):
    # theta in rad
    # vector don't rotate, but axis rotate, right hand rotate
    # np.dot(rotation_matrix3('x', theta), vector) returns the coordinate of the vector in a rotated new coordinate system
    if xyz == 'x':
        R = np.array([[1, 0, 0], [0, cos(theta), sin(theta)], [0, -sin(theta), cos(theta)]])
    if xyz == 'y':
        R = np.array([[cos(theta), 0, -sin(theta)], [0, 1, 0], [sin(theta), 0, cos(theta)]])
    if xyz == 'z':
        R = np.array([[cos(theta), sin(theta), 0], [-sin(theta), cos(theta), 0], [0, 0, 1]])
    return R


def fitswrite(fileout, im, header=None):
    from astropy.io import fits
    import os
    if os.path.exists(fileout):
        os.remove(fileout)
    if header is None:
        fits.writeto(fileout, im, output_verify='fix', overwrite=True, checksum=False)
    else:        
        fits.writeto(fileout, im, header, output_verify='fix', overwrite=True, checksum=False)


def fitsread(filein):
    from astropy.io import fits
    head = '  '
    hdul = fits.open(filein)

    try:
        data0 = hdul[0].data.astype(np.float32)
        head = hdul[0].header
    except:
        hdul.verify('silentfix')
        data0 = hdul[1].data
        head = hdul[1].header

    return data0, head

def removelimb(im, center=None, RSUN=None):
  #  pip install polarTransform
    import polarTransform as pT
    from scipy import signal

    radiusSize, angleSize = 1024, 1800
    im = removenan(im)
    im2=im.copy()
    if center is None:
        T = (im.max() - im.min()) * 0.2 + im.min()
        arr = (im > T)
        import scipy.ndimage.morphology as snm
        arr=snm.binary_fill_holes(arr)
#        im2=(im-T)*arr
        Y, X = np.mgrid[:im.shape[0], :im.shape[1]]
        xc = (X * arr).astype(float).sum() / (arr*1).sum()
        yc = (Y * arr).astype(float).sum() / (arr*1).sum()
        center = (xc, yc)
        RSUN = np.sqrt(arr.sum() / np.pi)

    Disk = np.int8(disk(im.shape[0], im.shape[1], RSUN * 0.95))
    impolar, Ptsetting = pT.convertToPolarImage(im, center, radiusSize=radiusSize, angleSize=angleSize)
    profile = np.median(impolar, axis=0)
    profile = signal.savgol_filter(profile, 11, 3)
    Z = profile.reshape(-1, 1).T.repeat(impolar.shape[0], axis=0)
    limb=Ptsetting.convertToCartesianImage(Z)
#    im2 = removenan(im / Ptsetting.convertToCartesianImage(Z))-1
#    im2 = im2 * Disk
    im = removenan(im /limb)
    im= im*Disk
    return im, center, RSUN, profile,limb

def diskcenter(im):
    from skimage import filters

  #  pip install polarTransform
    # import polarTransform as pT
    # from scipy import signal

    # radiusSize, angleSize = 1024, 1800
    im = removenan(im)
    im2=im.copy()
    T=filters.thresholding.threshold_isodata(im)
    #print(T)
    # T = (im.max() - im.min()) * 0.2 + im.min()
    # T=800
    arr = (im > T)
    import scipy.ndimage.morphology as snm
    arr=snm.binary_fill_holes(arr)
#        im2=(im-T)*arr
    arr2=arr*np.maximum(im-T,0)
    Y, X = np.mgrid[:im.shape[0], :im.shape[1]]
    xc = (X * arr2).astype(float).sum() / (arr2*1).sum()-0.5
    yc = (Y * arr2).astype(float).sum() / (arr2*1).sum()-0.5
    center = (xc, yc)
    RSUN = np.sqrt(arr.sum() / np.pi)

    
    return  center, RSUN, arr

def rebina(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)

def canny_limb(image):
    
    a, b = image.shape
    
    kernel_b = np.zeros([25, 25])
    radius = 12
    fwhm = 5
    for y in range(- radius, radius + 1):
        for x in range(- radius, radius + 1):
            v = 1.0 / (2 * pi * fwhm ** 2) * e ** (-1.0 / (2 * fwhm ** 2) * (x ** 2 + y ** 2))
            kernel_b[y + radius, x + radius] = v
    kernel_b = kernel_b / np.sum(kernel_b)
    convx = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    convy = np.array([[-1, -2, -1], [0, 0, 0], [1, 2, 1]])
    
    img_blur = signal.convolve2d(image, kernel_b, 'same')
    gradientx = signal.convolve2d(img_blur, convx, 'same')
    gradienty = signal.convolve2d(img_blur, convy, 'same')
    gradient_ori = np.abs(gradientx) + np.abs(gradienty) 
    gradient2 = np.zeros([a, b])
    tan22_5 = tan(22.5 / 180 * pi)
    tan67_5 = tan(67.5 / 180 * pi)
    for i in range(1, a - 1):
        for j in range(1, b - 1):
            loc_g = gradient_ori[i, j]
            if gradientx[i, j] != 0:
                yx_ratio = gradienty[i, j] / gradientx[i, j]
            else:
                yx_ratio = 9999
            nei_a, nei_b = 0, 0
            if yx_ratio >= tan67_5 or yx_ratio < - tan67_5:
                nei_a = gradient_ori[i + 1, j]
                nei_b = gradient_ori[i - 1, j]
            elif yx_ratio >= tan22_5 and yx_ratio < tan67_5:
                nei_a = gradient_ori[i + 1, j + 1]
                nei_b = gradient_ori[i - 1, j - 1]                
            elif yx_ratio >= - tan22_5 and yx_ratio < tan22_5:
                nei_a = gradient_ori[i, j + 1]
                nei_b = gradient_ori[i, j - 1]                   
            elif yx_ratio >= - tan67_5 and yx_ratio < - tan22_5:
                nei_a = gradient_ori[i + 1, j - 1]
                nei_b = gradient_ori[i - 1, j + 1]
            if loc_g >= nei_a and loc_g >= nei_b:
                gradient2[i, j] = loc_g
    gradient = np.zeros([a, b])
    thred_l = 0.1 * np.max(gradient2)
    thred_h = 0.3 * np.max(gradient2)
    for i in range(1, a - 1):
        for j in range(1, b - 1):
            if gradient2[i, j] <= thred_l:
                gradient[i, j] = 0
            elif gradient2[i, j] >= thred_h:
                gradient[i, j] = 1
            elif (gradient2[i - 1:i + 1, j - 1:j + 1] >= thred_h).any():
                gradient[i, j] = 1
    
    return gradient

def median_filter(image, k):

    height, width = image.shape

    edge = int((k - 1)/2)
    if height - 1 - edge <= edge or width - 1 - edge <= edge:
        return None
    
    new_arr = np.zeros([height, width])
    for i in range(height):
        for j in range(width):
            if i <= edge - 1 or i >= height - 1 - edge or j <= edge - 1 or j >= height - edge - 1:
                new_arr[i, j] = image[i, j]
            else:
                new_arr[i, j] = np.median(image[i - edge:i + edge + 1, j - edge:j + edge + 1])
    
    return new_arr

def ini_diskpara(xc, yc):
    
    lenth = len(xc)
    lx, lx2, ly, ly2, lxy, lx2y, lxy2, lx3, ly3 = xc, [xc[i] ** 2 for i in range(lenth)], yc, [yc[i] ** 2 for i in range(lenth)], \
    [xc[i] * yc[i] for i in range(lenth)], [xc[i] ** 2 * yc[i] for i in range(lenth)], [xc[i] * yc[i] ** 2 for i in range(lenth)], \
    [xc[i] ** 3 for i in range(lenth)], [yc[i] ** 3 for i in range(lenth)]
    CC = lenth * sum(lx2) - sum(lx) * sum(lx)
    DD = lenth * sum(lxy) - sum(lx) * sum(ly)
    EE = lenth * sum(lx3) + lenth * sum(lxy2) - (sum(lx2) + sum(ly2)) * sum(lx)
    GG = lenth * sum(ly2) - sum(ly) * sum(ly)
    HH = lenth * sum(lx2y) + lenth * sum(ly3) - (sum(lx2) + sum(ly2)) * sum(ly)
    aa, bb = (HH * DD - EE * GG) / (CC * GG - DD ** 2), (HH * CC - EE * DD) / (DD ** 2 - CC * GG)
    cc = -(sum(lx2) + sum(ly2) + aa * sum(lx) + bb * sum(ly)) / lenth
    
    return -aa / 2, -bb / 2

def circle_center(image, shrink_num = 4, scale_num = 100):
    
    filter_num = 100 / shrink_num
    
    height, width = image.shape
    print(image.shape)
    s_height, s_width = floor(height / shrink_num) * shrink_num, floor(width / shrink_num) * shrink_num
    
    im1 = rebina(image[:s_height, :s_width], [int(height / shrink_num), int(width / shrink_num)])
    img_med1 = median_filter(im1, filter_num)
    canny_img1 = canny_limb(img_med1)
    yc = list(np.where(canny_img1 == 1)[0] / scale_num)
    xc = list(np.where(canny_img1 == 1)[1] / scale_num)

    ccx, ccy = ini_diskpara(xc, yc)

    return (ccx * scale_num + 0.5) * shrink_num, (ccy * scale_num + 0.5) * shrink_num

def theory_rsun(strtime, satpos, bin_num = 1):
    t = Time(strtime)
    sat_rx = satpos[0]
    sat_ry = satpos[1]
    sat_rz = satpos[2]
    astropy.coordinates.solar_system_ephemeris.set(config.de_file_url)
    sun_6 = astropy.coordinates.get_body_barycentric_posvel('sun', t)
    earth_6 = astropy.coordinates.get_body_barycentric_posvel('earth', t)

    sun_pos = sun_6[0]
    earth_pos = earth_6[0]

    # sun position is in ICRF (aberration not corrected), satellite position is in J2000.0 coordinates (need to further examine how the two parameters change due to different coordinate systems)
    # sat2sun_pos infers the vector pointing from satellite to the Sun
    sat2sun_rx = sun_pos.x.value - earth_pos.x.value - sat_rx
    sat2sun_ry = sun_pos.y.value - earth_pos.y.value - sat_ry
    sat2sun_rz = sun_pos.z.value - earth_pos.z.value - sat_rz
    sat2sun_pos = np.array([[sat2sun_rx], [sat2sun_ry], [sat2sun_rz]])

    normalize_factor = sqrt(sat2sun_pos[0] ** 2 + sat2sun_pos[1] ** 2 + sat2sun_pos[2] ** 2)
    sat2sun_pos_normalize = sat2sun_pos / normalize_factor # the satellite position to the sun position (point from satellite to sun)

    r = asin((6.963 * 10 ** 5) / normalize_factor) * 180 / pi * 60 * 60 / 0.5218 / bin_num # in pixel

    return r

def imnorm(im, mx=0, mi=0):
    #   图像最大最小归一化 0-1
    if mx != 0 and mi != 0:
        pass
    else:
        mi, mx = np.min(im), np.max(im)

    im2 = removenan((im - mi) / (mx - mi))

    arr1 = (im2 > 1)
    im2[arr1] = 1
    arr0 = (im2 < 0)
    im2[arr0] = 0

    return im2


def removenan(im, key=0):
    """
    remove NAN and INF in an image
    """
    im2 = np.copy(im)
    arr = np.isnan(im2)
    im2[arr] = key
    arr2 = np.isinf(im2)
    im2[arr2] = key

    return im2


def showim(im,k=3,cmap='gray',interpolation='nearest'):
    mi = np.max([im.min(), im.mean() - k * im.std()])
    mx = np.min([im.max(), im.mean() + k * im.std()])
    if len(im.shape) == 3:
        plt.imshow(im, vmin=mi, vmax=mx)
    else:
        plt.imshow(im, vmin=mi, vmax=mx, cmap=cmap,interpolation=interpolation)

    return


def zscore2(im):
    im = (im - np.mean(im)) / im.std()
    return im


def disk(M, N, r0):
    X, Y = np.meshgrid(np.arange(int(-(N / 2)), int(N / 2)), np.linspace(-int(M / 2), int(M / 2) - 1, M))
    r = (X) ** 2 + (Y) ** 2
    r = (r ** 0.5)
    im = r < r0
    return im


#def fgauss(M, N, I, x0, y0, r):
#    # 产生高斯图像
#
#    r = r * r * 2
#    x = np.arange(0, M)
#    x = x - M / 2 + x0 - 1
#    y = np.arange(0, N)
#    y = y - N / 2 + y0 - 1
#    w1 = np.exp(-x ** 2 / r)
#    w2 = np.exp(-y ** 2 / r)
#    w2 = np.reshape(w2, (-1, 1))
#    f = I * w1 * w2
#    return f
#
#
#def showmesh(im):
#    X, Y = np.mgrid[:im.shape[0], :im.shape[1]]
#    from mpl_toolkits.mplot3d import Axes3D
#    figure = plt.figure('mesh')
#    axes = Axes3D(figure)
#
#    axes.plot_surface(X, Y, im, cmap='rainbow')
#    return


def create_gif(images, gif_name, duration=1):
    import imageio
    frames = []
    # Read
    T = images.shape[2]
    for i in range(T):
        frames.append(np.uint8(imnorm(images[:, :, i]) * 255))
    #    # Save
    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)

    return


def immove2(im,dx=0,dy=0):
#    im2,para=array2img(im)
    im2=im.copy()
    tform = SimilarityTransform(translation=(dx,dy))
    im2 = warp(im2, tform.inverse, output_shape=(im2.shape[0], im2.shape[1]),mode='constant',cval=0)
#    im2=img2array(im2,para)
    return im2



def imcenterpix(im):
    X0=(im.shape[0])//2
    Y0=(im.shape[1])//2
    cen=(X0,Y0)
    return cen



def xcorrcenter(standimage, compimage, R0=2, flag=0):
    # if flag==1,standimage 是FFT以后的图像，这是为了简化整数象元迭代的运算量。直接输入FFT以后的结果，不用每次都重复计算
    try:
        M, N = standimage.shape

        standimage = zscore2(standimage)
        s = fft.fft2(standimage)

        compimage = zscore2(compimage)
        c = np.fft.ifft2(compimage)

        sc = s * c
        im = np.abs(fft.fftshift(fft.ifft2(sc)))  # /(M*N-1);%./(1+w1.^2);
        cor = im.max()
        if cor == 0:
            return 0, 0, 0

        M0, N0 = np.where(im == cor)
        m, n = M0[0], N0[0]

        if flag:
            m -= M / 2
            n -= N / 2
            # 判断图像尺寸的奇偶
            if np.mod(M, 2): m += 0.5
            if np.mod(N, 2): n += 0.5

            return m, n, cor
        # 求顶点周围区域的最小值
        immin = im[(m - R0):(m + R0 + 1), (n - R0):(n + R0 + 1)].min()
        # 减去最小值
        im = np.maximum(im - immin, 0)
        # 计算重心
        x, y = np.mgrid[:M, :N]
        area = im.sum()
        m = (np.double(im) * x).sum() / area
        n = (np.double(im) * y).sum() / area
        # 归算到原始图像
        m -= M / 2
        n -= N / 2
        # 判断图像尺寸的奇偶
        if np.mod(M, 2): m += 0.5
        if np.mod(N, 2): n += 0.5
    except:
        print('Err in align_Subpix routine!')
        m, n, cor = 0, 0, 0
    return m, n, cor
def cc(standimage, compimage, flag=0,win=1):
 
        M, N = standimage.shape
        if flag==0:
            standimage = zscore2(standimage)
            s = fft.fft2(standimage)
        else:    
            s=standimage
            
        c = zscore2(compimage)
        c = fft.fft2(c)
    
        sc = s * np.conj(c)*win
        im = np.abs(fft.fftshift(fft.ifft2(sc)))  # /(M*N-1);%./(1+w1.^2);
#        im=im/(im.shape[0]*im.shape[1])
        cor = im.max()
        if cor == 0:
            return 0, 0, 0

        M0, N0 = np.where(im == cor)
        m, n = M0[0], N0[0]

        m -= M / 2
        n -= N / 2
        # 判断图像尺寸的奇偶
        if np.mod(M, 2): m += 0.5
        if np.mod(N, 2): n += 0.5

 #       c=np.abs(c)
        cor/=standimage.size
        return m, n, cor,im
    
        # 求顶点周围区域的最小值

def gc(standimage, compimage, flag=0,win=1):
 
        M, N = standimage.shape
        if flag==0:
            standimage = zscore2(standimage)
            s = fft.fft2(standimage)
        else:    
            s=standimage
            
        c = zscore2(compimage)
        c=compimage.copy()
        c = fft.fft2(c)
    
        sc = s * np.conj(c)*win
        im = np.abs(fft.fftshift(fft.ifft2(sc/np.abs(sc))))  # /(M*N-1);%./(1+w1.^2);
#        im=im/(im.shape[0]*im.shape[1])
        cor = im.max()
        if cor == 0:
            return 0, 0, 0

        M0, N0 = np.where(im == cor)
        m, n = M0[0], N0[0]

        m -= M / 2
        n -= N / 2
        # 判断图像尺寸的奇偶
        if np.mod(M, 2): m += 0.5
        if np.mod(N, 2): n += 0.5

 #       c=np.abs(c)

        return m, n, cor,im
    
        # 求顶点周围区域的最小值


def immove(image, dx, dy):
    """
    image shift by subpix
    """
    # The shift corresponds to the pixel offset relative to the reference image
    from scipy.ndimage import fourier_shift
    if dx == 0 and dy == 0:
        offset_image = image
    else:
        shift = (dx, dy)
        offset_image = fourier_shift(fft.fft2(image), shift)
        offset_image = np.real(fft.ifft2(offset_image))

    return offset_image

def immove_fft(fftim, dx, dy):
    """
    image shift by subpix
    """
    # The shift corresponds to the pixel offset relative to the reference image
    from scipy.ndimage import fourier_shift
    if abs(dx) < 0.1 and abs(dy) < 0.1:
        offset_image = fftim
    else:
        shift = (dx, dy)
        offset_image = fourier_shift(fftim, shift)

#        offset_image = np.real(fft.ifft2(offset_image))

    return offset_image

def combin_img(z,Ncol):
    import numpy as np
    Nrow=z.shape[0]//Ncol

    for i in range(Ncol):
        for j in range(Nrow):
            if j==0:
                row=(z[j+i*Nrow])
            else:
                row=np.hstack((row,z[j+i*Nrow]))
        if i==0:        
            col=row
        else:
            col=np.vstack((col,row))
    return col

def rebin(arr, nbin):

    m=arr.shape[0]//nbin
    n=arr.shape[1]//nbin
    shape = (m, nbin,n, nbin)
    return arr.reshape(shape).sum(-1).sum(1)

def showmesh(im,flag=0):
#    x=np.arange(0,im.shape[0])
#    y=np.arange(0,im.shape[1])
    X,Y=np.mgrid[:im.shape[0],:im.shape[1]]
    from mpl_toolkits.mplot3d import Axes3D
    figure = plt.figure('mesh '+str(flag))
    axes = Axes3D(figure)
#    plt.show()
    axes.plot_surface(X,Y,im,cmap='rainbow')
    return  

def cc_gpu(standimage, compimage, flag=1):
	# if flag==1,standimage 是FFT以后的图像，这是为了简化整数象元迭代的运算量。直接输入FFT以后的结果，不用每次都重复计算
        import cupyx.scipy.fftpack as gfft
        if cp.cuda.Device(0):
            cp.cuda.Device(0).use()
            if flag==0:
#                standimage = zscore2(standimage)

#                s_gpu = cp.ndarray((M,N),dtype=np.complex64)
                s_gpu =standimage
                s_gpu = gfft.fft2(s_gpu)
            else:
                s_gpu=standimage.copy()
            #im = np.ndarray((M,N),dtype=np.complex64)
            M, N = s_gpu.shape
            #prepare 3 arrays on gpu

            # c_gpu = cp.ndarray((M,N),dtype=np.complex64)
            # sc_gpu = cp.ndarray((M,N),dtype=np.complex64)
            #copy images from host to gpu
            
            c_gpu = compimage
            c_gpu = gfft.ifft2(c_gpu)
            
            sc_gpu = cp.multiply(s_gpu,c_gpu)
            sc_gpu = gfft.ifft2(sc_gpu)
            
            im = cp.abs(cp.fft.fftshift(sc_gpu))     
#            im = cp.asnumpy(sc_gpu)
#            im=sc_gpu
        else:

            if flag==0:
                standimage = zscore2(standimage)
                s = fft.fft2(standimage)
            else:    
                s=standimage.copy()
                
                c = zscore2(compimage)
                c = fft.ifft2(c)
        
                sc = s * c
                im = np.abs(fft.fftshift(fft.ifft2(sc)))  # /(M*N-1);%./(1+w1.^2) 
    		#print("GPU ended!")
        cor = im.max()
        if cor == 0:
            return 0, 0, 0


        M0, N0 = np.where(im == cor)
        m, n = M0[0].get(), N0[0].get()

	
        R0=3
        try:
            immin = im[(m - R0):(m + R0 + 1), (n - R0):(n + R0 + 1)].min()
            # 减去最小值
            im = cp.maximum(im - immin, 0)
            # 计算重心
            x, y = cp.mgrid[:M, :N]
            area = im.sum()
            m = (im * x).sum() / area
            n = (im * y).sum() / area
            # 归算到原始图像
            m -= M / 2
            n -= N / 2
            # 判断图像尺寸的奇偶
            if cp.mod(M, 2): m += 0.5
            if cp.mod(N, 2): n += 0.5
        except:
            print('Err in align_Subpix routine!')
            m, n, cor = 0, 0, 0
        return cp.int(m+0.5), cp.int(n+0.5), cor
def cc_gpu_all(s_gpu, c_gpu):
    
    M, N = s_gpu.shape

    # c_gpu = cp.asarray(c_gpu)
    # c_gpu = zscore2_gpu(c_gpu)
    # c_gpu = fft.ifft2(c_gpu)
    
    sc_gpu = cp.multiply(s_gpu,c_gpu.conj())
    sc_gpu = cp.fft.ifft2(sc_gpu)
    
    im = cp.abs(cp.fft.fftshift(sc_gpu))     
    cor = cp.max(im)
#         if cor == 0:
#             return 0, 0, 0

#         # M0, N0 = cp.unravel_index(cp.argmax(im),(M,N))
#         # m, n = M0.get(), N0.get()


    M0, N0 = cp.where(im == cor)
    m, n = M0[0], N0[0]

    m =m- M / 2
    n =n- N / 2
		# 判断图像尺寸的奇偶
    if cp.mod(M, 2): m += 0.5
    if cp.mod(N, 2): n += 0.5
#        showmesh(im.get())
    return m, n, cor/s_gpu.size  



def ccsub_gpu_all(s_gpu, c_gpu):
    
    M, N = s_gpu.shape

    # c_gpu = cp.asarray(c_gpu)
    # c_gpu = zscore2_gpu(c_gpu)
    # c_gpu = fft.ifft2(c_gpu)
    
    sc_gpu = cp.multiply(s_gpu,c_gpu.conj())
    sc_gpu = cp.fft.ifft2(sc_gpu)
    
    im = cp.abs(cp.fft.fftshift(sc_gpu))     
    cor = cp.max(im)
#         if cor == 0:
#             return 0, 0, 0

#         # M0, N0 = cp.unravel_index(cp.argmax(im),(M,N))
#         # m, n = M0.get(), N0.get()


    M0, N0 = cp.where(im == cor)
    m, n = M0[0], N0[0]

#         m =m- M / 2
#         n =n- N / 2
# 		# 判断图像尺寸的奇偶
#         if cp.mod(M, 2): m += 0.5
#         if cp.mod(N, 2): n += 0.5
# #        showmesh(im.get())
#         return m, n, cor  
    R0=3
    try:
        immin = im[(m - R0):(m + R0 + 1), (n - R0):(n + R0 + 1)].min()
        # 减去最小值
        im = cp.maximum(im - immin, 0)
        # 计算重心
        x, y = cp.mgrid[:M, :N]
        area = im.sum()
        m = (im * x).sum() / area
        n = (im * y).sum() / area
        # 归算到原始图像
        m -= M / 2
        n -= N / 2
        # 判断图像尺寸的奇偶
        if cp.mod(M, 2): m += 0.5
        if cp.mod(N, 2): n += 0.5
    except:
        print('Err in align_Subpix routine!')
        m, n, cor = 0, 0, 0
    return m, n, cor
    
def pc_gpu_all(s_gpu, c_gpu):
    
        M, N = s_gpu.shape

        # c_gpu = cp.asarray(c_gpu)
        # c_gpu = zscore2_gpu(c_gpu)
        # c_gpu = fft.ifft2(c_gpu)
        
        sc_gpu = cp.multiply(s_gpu,c_gpu.conj())
        sc_gpu = cp.fft.ifft2(sc_gpu/cp.abs(sc_gpu))
        
        im = cp.abs(cp.fft.fftshift(sc_gpu))     

        cor = cp.max(im)
        if cor == 0:
            return 0, 0, 0

        # M0, N0 = cp.unravel_index(cp.argmax(im),(M,N))
        # m, n = M0.get(), N0.get()


        m, n = cp.where(im == cor)
#        m, n = M0[0], N0[0]

        m =m- M / 2
        n =n- N / 2
		# 判断图像尺寸的奇偶
        if cp.mod(M, 2): m += 0.5
        if cp.mod(N, 2): n += 0.5
#        showmesh(im.get())
        return m, n, cor        
def gc_gpu(standimage, compimage, flag=0):
	# if flag==1,standimage 是FFT以后的图像，这是为了简化整数象元迭代的运算量。直接输入FFT以后的结果，不用每次都重复计算
        M, N = standimage.shape
        if cp.cuda.Device(0):
            cp.cuda.Device(0).use()
            if flag==0:
                standimage = zscore2(standimage)

#                s_gpu = cp.ndarray((M,N),dtype=np.complex64)
                s_gpu = cp.asarray(standimage)
                s_gpu = cp.fft.fft2(s_gpu)
            else:
                s_gpu=standimage.copy()
            #im = np.ndarray((M,N),dtype=np.complex64)

            #prepare 3 arrays on gpu
            compimage = zscore2(compimage)
            # c_gpu = cp.ndarray((M,N),dtype=np.complex64)
            # sc_gpu = cp.ndarray((M,N),dtype=np.complex64)
            #copy images from host to gpu
            
            c_gpu = cp.asarray(compimage)
            c_gpu = cp.fft.ifft2(c_gpu)
            
            sc_gpu = cp.multiply(s_gpu,c_gpu)
            sc_gpu=sc_gpu/cp.abs(sc_gpu)
            sc_gpu = cp.fft.ifft2(sc_gpu)
            
            sc_gpu = cp.abs(cp.fft.fftshift(sc_gpu))     
            im = cp.asnumpy(sc_gpu)
        else:

            if flag==0:
                standimage = zscore2(standimage)
                s = fft.fft2(standimage)
            else:    
                s=standimage.copy()
                
                c = zscore2(compimage)
                c = np.fft.ifft2(c)
        
                sc = s * c
                sc=sc/np.abs(sc)
                im = np.abs(fft.fftshift(fft.ifft2(sc)))  # /(M*N-1);%./(1+w1.^2) 
    		#print("GPU ended!")
        cor = im.max()
        if cor == 0:
            return 0, 0, 0


        M0, N0 = np.where(im == cor)
        m, n = M0[0], N0[0]

	
        m -= M / 2
        n -= N / 2
		# 判断图像尺寸的奇偶
        if np.mod(M, 2): m += 0.5
        if np.mod(N, 2): n += 0.5
        
        
        return m, n, cor  
def makeGaussian(size, sigma = 3, center=None):

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / sigma**2)   

def imgcut(A, X, Y):
    """
    get subimage
    A: image narray
    X,Y: narray [2,3,4,5,.....100]
    """
    try:
        B = A[X[0]:X[-1], Y[0]:Y[-1]]
    except:
        B = A
        print('Warning:ROI is out of Image range. Whole image is selected as ROI')

    return B
def imshift(im,translation=[0,0]):

    # tform = EuclideanTransform(translation=translation)
    # im = warp(im2, tform.inverse, output_shape=(im2.shape[0], im2.shape[1]),mode='reflect')

    """
    shift an image by pixels
    """
    translation=(np.array(translation)).astype('int')
    im1 = im.copy()
    im1 = np.roll(im1, translation[0], axis=0)
    im1 = np.roll(im1, translation[1], axis=1)
    return im1  
def imshift_gpu(im,translation=[0,0]):

    # tform = EuclideanTransform(translation=translation)
    # im = warp(im2, tform.inverse, output_shape=(im2.shape[0], im2.shape[1]),mode='reflect')

    """
    shift an image by pixels
    """
    translation=(cp.int32(translation))
    im1 = im.copy()
    im1 = cp.roll(im1, translation[0], axis=0)
    im1 = cp.roll(im1, translation[1], axis=1)
    return im1  
def imrotate(im,para):
    return rotate(im,para,mode='reflect')

def imresize(im,para):
    return resize(im,para,mode='reflect')

def imrescale(im,para):
    return rescale(im,para,mode='reflect')

def imtransform(im,scale=1,rot=0,translation=[0,0]):
    im2=im.copy()
    tform = SimilarityTransform(translation=translation)
    im2 = warp(im2, tform.inverse, output_shape=(im2.shape[0], im2.shape[1]),mode='reflect')

    im2=rotate(im2,rot,mode='reflect')
    im2=rescale(im2,scale,mode='reflect')

    return im2   

def mkdir(path):    # 引入模块
    import os

    path = path.strip()
    # 去除尾部 \ 符号
    path = path.rstrip("\\")

    # 判断路径是否存在
    isExists = os.path.exists(path)

    # 判断结果
    if isExists:
        return False
    else:
        os.makedirs(path)
        return True
    
def ring_window(N):

    x = np.arange(0, N, 1, float)
    y = x[:,np.newaxis]
    
    window = np.exp(-(np.sqrt((x - N // 2) ** 2 + (y - N // 2) ** 2) - N // 5) ** 2 / (0.02 * N * N))
    return window 

def polpow(im0,order=1,method=0):
  #  pip install polarTransform
    import polarTransform as pT

    m=im0.shape[0]//2
    impolar, Ptsetting = pT.convertToPolarImage(im0,finalRadius=m,radiusSize=m, angleSize=360,order=order)
    if method==0:
        profile = np.median(impolar, axis=0)
    elif method ==2:     
        profile = np.max(impolar, axis=0)
    else:     
        profile = np.mean(impolar, axis=0)
        
    return profile,impolar
   
def R_tukey(width,alpha=0.5,order=1):
  #  pip install polarTransform
    import polarTransform as pT
    from scipy import signal
    profile=signal.windows.tukey(width,alpha)[width//2:]   
    z = profile.reshape(1,-1).repeat(360, axis=0)

    im=pT.convertToCartesianImage(z,order=order)
    return im[0]   

def R_tukey2(width,alpha=0.5):
  #  pip install polarTransform
    import polarTransform as pT
    from scipy import signal
    profile=signal.windows.tukey(width,alpha)  
    
    im=np.dot(profile.reshape(-1,1),profile.reshape(1,-1))

    return im  
def tukey(width,alpha=0.5):

    from scipy import signal

    s=signal.windows.tukey(width,alpha=alpha)
    s=s[np.newaxis]
    z=np.dot(s.T,s)
    return z
def H_tukey(width,alpha=0.5):

    from scipy import signal
    profile=signal.windows.tukey(width,alpha)   
    z = profile.reshape(1,-1).repeat(width, axis=0)

    im=z*rotate(z,60)*rotate(z,120)
    return im   

def gauss_win(size,k=10):
#    std=k
    from scipy import signal
    profile=signal.windows.gaussian(size,std=k)
    win=profile.reshape(-1,1)*profile
    return win

def hamming_win(size):
    std=size/10
    from scipy import signal
    profile=signal.windows.hamming(size+1)
    win=profile.reshape(-1,1)*profile
    return win[:-1,:-1]
def hanning_win(size):
    std=size/10
    from scipy import signal
    profile=signal.windows.hanning(size)
    win=profile.reshape(-1,1)*profile
    return win
def test_gpu():
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()

    print("used bytes: %s"%mempool.used_bytes())
    print("total bytes: %s\n"%mempool.total_bytes())
def cal_spe_gpu(im):
    im=cp.asarray(im)
    spe=cp.abs(fft.fft2(im))
    spe=cp.fft.fftshift(spe)
    spe*=spe
    return spe

def nonrigid(im1,im2,size=21,pyr_scale=0.5, levels=3,  iterations=5, poly_n=5, poly_sigma=1.2, flags=0):
    import cv2
    flow = cv2.calcOpticalFlowFarneback(im1, im2, flow=None, pyr_scale= pyr_scale, levels=levels, winsize=size, iterations=iterations, poly_n=poly_n, poly_sigma=poly_sigma, flags=0)

    h,w=im1.shape
    x1, y1 = np.meshgrid(np.arange(w), np.arange(h))
    x1=x1.astype('float32')+flow[:,:,0]
    y1=y1.astype('float32')+flow[:,:,1]
    
    out=cv2.remap(im2,x1,y1,borderMode=cv2.BORDER_REFLECT_101, interpolation=cv2.INTER_NEAREST)
    return flow[:,:,0],flow[:,:,1],out

def create_gif3(images, gif_name, duration=1,sigma=3):
    import imageio
    frames = []
    # Read
    T = images.shape[2]
    for i in range(T):
        mm=images[:, :, i].mean()
        st=images[:, :, i].std()
        frames.append(np.uint8(imnorm(images[:, :, i],mm+sigma*st,max(mm-sigma*st,0)) * 255))
    #    # Save
    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
def create_gif_maxmin(images, gif_name, duration=1,mx=100,mi=-100):
    import imageio
    frames = []
    # Read
    T = images.shape[2]
    for i in range(T):

        frames.append(np.uint8(imnorm(images[:, :, i],mx=mx,mi=mi) * 255))
    #    # Save
    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)    
    
def create_gif2(images, gif_name, duration=1,sigma=3):
    import imageio
    frames = []
    w0=images.shape[0]//4
    w1=images.shape[1]//4
    # Read
    T = images.shape[2]
    for i in range(T):
        mm=images[w0:-w0, w1:-w1, i].mean()
        st=images[w0:-w0, w1:-w1, i].std()
        frames.append(np.uint8(imnorm(images[:, :, i],mm+sigma*st,mm-sigma*st) * 255))
    #    # Save
    imageio.mimsave(gif_name, frames, 'GIF', duration=duration)
    
def butterworthPassFilter(image,d,n):
    f = np.fft.fft2(image)
    fshift = np.fft.fftshift(f)
    
    def make_transform_matrix(d):
        transfor_matrix = np.zeros(image.shape)
        center_point = tuple(map(lambda x:(x)/2,f.shape))
        for i in range(transfor_matrix.shape[0]):
            for j in range(transfor_matrix.shape[1]):
                def cal_distance(pa,pb):
                    from math import sqrt
                    dis = sqrt((pa[0]-pb[0])**2+(pa[1]-pb[1])**2)
                    return dis
                dis = cal_distance(center_point,(i,j))
                transfor_matrix[i,j] = 1/((1+(dis/d))**n)
        return transfor_matrix
    d_matrix = make_transform_matrix(d)
    new_img = np.abs(np.fft.ifft2(np.fft.ifftshift(fshift*d_matrix)))
    return new_img,d_matrix