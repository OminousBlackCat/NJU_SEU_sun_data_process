import config
import numpy as np
from astropy.io import fits


def main():
    file_data = fits.open(config.data_dir_path + '/' + 'RSM20220120T062139-0013-1243.fts')
    image_data = np.array(file_data[0].data, dtype=float)
    a = file_data[0].header
    a['SIMPLE']
    print(file_data[0].header)


if __name__ == '__main__':
    main()