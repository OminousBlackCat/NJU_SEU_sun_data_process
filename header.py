import config
import numpy as np
from astropy.io import fits
import os

standard_header_list_HA = [
    {
        'key': 'SIMPLE',
        'value': 'T',
        'comment': ''
    },
    {
        'key': 'BITPIX',
        'value': 16,
        'comment': 'number of bits per data pixel'
    },
    {
        'key': 'NAXIS',
        'value': 3,
        'comment': 'Number of data axes'
    },
    {
        'key': 'NAXIS1',
        'value': 'T',
        'comment': 'Length of data axis 1 (slit dimension)'
    },
    {
        'key': 'NAXIS2',
        'value': 'T',
        'comment': ' Length of data axis 2 (scanning steps)'
    },
    {
        'key': 'NAXIS3',
        'value': 'T',
        'comment': 'Length of data axis 3 (wavelength dimension)'
    },
    {
        'key': 'TELESCOP',
        'value': 'CHASE',
        'comment': 'HIS'
    },
    {
        'key': 'BIN',
        'value': 2,
        'comment': 'Binning mode'
    },
    {
        'key': 'DATE_OBS',
        'value': '2022-01-20T06:08:39',
        'comment': 'Observation time'
    },
    {
        'key': 'CRPIX1',
        'value': 'T',
        'comment': ' X coordinate of solar center in pixel'
    },
    {
        'key': 'CRVAL1',
        'value': 'T',
        'comment': 'E-W arcsec of solar center'
    },
    {
        'key': 'CRPIX2',
        'value': 'T',
        'comment': 'Y coordinate of solar center in pixel'
    },
    {
        'key': 'CUNIT2',
        'value': 'T',
        'comment': 'unit of y-axis'
    },
    {
        'key': 'CDELT2',
        'value': 'T',
        'comment': 'Pixel spatial resolution along y-axis'
    },
    {
        'key': 'CRVAL2',
        'value': 'T',
        'comment': 'N-S arcsec of solar center '
    },
    {
        'key': 'CRPIX3',
        'value': 'T',
        'comment': 'Reference pixel at the wavelength dimension'
    },
    {
        'key': 'CUNIT3',
        'value': 'T',
        'comment': 'unit of z-axis'
    },
    {
        'key': 'CDELT3',
        'value': 'T',
        'comment': 'Angstrom per pixel'
    },
    {
        'key': 'CRVAL3',
        'value': 'T',
        'comment': 'Wavelength at the reference pixel'
    },
    {
        'key': 'BZERO',
        'value': 'T',
        'comment': 'Data is Unsigned Integer'
    },
    {
        'key': 'BSCALE',
        'value': 'T',
        'comment': 'default scaling factor'
    },
    {
        'key': 'STEPTIME',
        'value': 'T',
        'comment': 'One-step scanning duration,unit(s)'
    },
    {
        'key': 'PRODATE',
        'value': '2022-03-31T19:40:12',
        'comment': 'File creation date'
    },
    {
        'key': 'STR_TIME',
        'value': 'T',
        'comment': 'scanning start time'
    },
    {
        'key': 'END_TIME',
        'value': 'T',
        'comment': 'scanning end time'
    },
    {
        'key': 'EXP_TIME',
        'value': 'T',
        'comment': 'Exposure time for raster scanning,unit(s)'
    },
    {
        'key': 'FRM_NUM',
        'value': 'T',
        'comment': 'Frame number of a scanning series'
    },
    {
        'key': 'SPECLINE',
        'value': 'T',
        'comment': 'Spectral line'
    },
    {
        'key': 'LINECORE',
        'value': 'T',
        'comment': 'Wavelength of Ha line core'
    },
    {
        'key': 'SAT_POS1',
        'value': 'T',
        'comment': 'Satellite x-position in J2000 coordinate system'
    },
    {
        'key': 'SAT_POS2',
        'value': 'T',
        'comment': 'Satellite y-position in J2000 coordinate system'
    },
    {
        'key': 'SAT_POS3',
        'value': 'T',
        'comment': 'Satellite z-position in J2000 coordinate system'
    },
    {
        'key': 'SAT_VEL1',
        'value': 'T',
        'comment': 'Satellite x-velocity in J2000 coordinate system'
    },
    {
        'key': 'SAT_VEL2',
        'value': 'T',
        'comment': 'Satellite y-velocity in J2000 coordinate system'
    },
    {
        'key': 'SAT_VEL3',
        'value': 'T',
        'comment': 'Satellite z-velocity in J2000 coordinate system'
    },
    {
        'key': 'INST_ROT',
        'value': 'T',
        'comment': '[deg] Angle of Solar North pole wrt y-axis.'
    },
    {
        'key': 'R_SUN',
        'value': 'T',
        'comment': 'Radius of the Sun in pixels on the CCD'
    },
    {
        'key': 'RSUN_OBS',
        'value': 'T',
        'comment': 'Apparent solar radius seen by CHASE'
    },
    {
        'key': 'B0',
        'value': 'T',
        'comment': '[deg] latitude of point at disk center. '
    },
]


def read_header_from_txt(txtPath):
    tempList = []
    with open(txtPath) as f:
        lines = f.readlines()
    for line in lines:
        tempDic = {
            'key': '',
            'value': '',
            'comment': ''
        }
        l_split_first = line.split('=')
        if len(l_split_first) < 2:
            continue
        l_split_second = l_split_first[1].split('/')
        tempDic['key'] = l_split_first[0].strip()
        tempDic['comment'] = l_split_second[1].strip()
        tempList.append(tempDic)
    return tempList


def main():
    file_data = fits.open(config.data_dir_path + '/' + 'RSM20220120T062539-0017-1256.fts')
    b = fits.header.Header()
    temp = read_header_from_txt('D:\\QQ\\WechatTemp\\WeChat Files\\wxid_vch2mdk54tie22\\FileStorage\\File\\2022-03'
                                '\\HA_header.txt')
    for h in temp:
        b.set(h['key'], value=0, comment=h['comment'])
    primaryHDU = fits.PrimaryHDU(data=file_data[0].data, header=b)
    greyHDU = fits.HDUList([primaryHDU])
    greyHDU.writeto(config.data_dir_path + '/' + 'RSM20220120T062139-0013-1244_test.fts', overwrite=True)
    file_data.close()
    file_data = fits.open(config.data_dir_path + '/' + 'RSM20220120T062139-0013-1244_test.fts')
    a = file_data[0].header
    print(repr(a))


if __name__ == '__main__':
    main()
