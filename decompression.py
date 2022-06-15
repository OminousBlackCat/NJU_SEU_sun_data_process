from astropy.io import fits


def main():
    result = None
    filename = 'data/decompression_test/RSM20220130T191819_0000_FE.fits'
    with fits.open(filename) as hdul:
        imgData = hdul[1].data
        header = hdul[1].header
        print(repr(header))
        print('Done!')


if __name__ == '__main__':
    main()
