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
    filename = 'data/decompression_test/RSM20221017T182100_0005_FE.fts'
    file = filename
    print('Reading fits file....')
    dat = fitsread(file)[0]

    SPEC = dat[26]
    show_SPEC(SPEC, dat)

if __name__ == '__main__':
    main()
