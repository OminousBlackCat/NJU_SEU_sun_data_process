from astropy.io import fits


def main():
    result = None
    filename = 'data/decompression_test/RSM20220326T124304_0016_HA.fits'
    with fits.open(filename,  do_not_scale_image_data=True) as hdul:
        imgData = hdul[1].data
        header = hdul[1].header
        print(repr(header))
        print('Done!')


if __name__ == '__main__':
    main()
