from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import numpy as np


def main():
    result = None
    filename = 'data/decompression_test/RSM20220514T214718_0021_HA.fits'
    with fits.open(filename) as hdul:
        print(hdul.info())


if __name__ == '__main__':
    main()
