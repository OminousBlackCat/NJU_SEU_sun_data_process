from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 13:16:29 2020

@author: jkf
@editor: wxy
"""


def fitsread(filein):
    from astropy.io import fits
    head = '  '
    hdul = fits.open(filein, do_not_scale_image_data=True)

    try:
        data0 = hdul[0].data.astype(np.float32)
        head = hdul[0].header
    except:
        hdul.verify('silentfix')
        data0 = hdul[1].data
        head = hdul[1].header

    return data0, head


def show_SPEC(SPEC, DSPEC):
    X = SPEC.shape[1]
    Y = SPEC.shape[0]
    figSPEC = plt.figure('Ha Image', figsize=(12, 7))
    ax = plt.subplot(121)
    ax.imshow(SPEC, cmap='gray')
    plt.draw()
    ax.set_title('Ha Image')
    bx = plt.subplot(122)
    bx.set_title('Ha Specturm')
    plt.tight_layout()

    def onclick(event):
        ix, iy = event.xdata, event.ydata
        if (ix < X) and (iy < Y):
            get_x(ix, iy)

    def get_x(pX, pY):
        x = np.round(pX).astype('int')
        y = np.round(pY).astype('int')
        print(x, y)

        ax.plot(x, y, '+', ms=10)
        plt.draw()

        bx.plot(DSPEC[:, y, x], '.-')
        plt.tight_layout()

    cid = figSPEC.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    return SPEC


def main():
    def test(array: list):
        array[0] = 1

    def test_numpy(array: np.array):
        array[0] = 1

    test_arr = [0, 0, 0]
    test_np = np.zeros((3,))
    test(test_arr[1: 2])
    test_numpy(test_np[2: 3])
    print(test_arr)
    print(test_np)

    # filename = 'data/time_test/RSM20220913T112208_0000_HA.fits'
    # hdu = fits.open(filename, do_not_scale_image_data=True)
    # time_data = hdu[0].data
    # sun_data = hdu[1].data
    # print(time_data)
    # print(sun_data)
    # img = plt.imshow(time_data)
    # plt.show()
    # print('Reading fits file....')


if __name__ == '__main__':
    main()
