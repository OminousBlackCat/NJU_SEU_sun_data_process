import config
from astropy.io import fits
import suntools
from astropy.io.fits import HDUList, Header

# 所有在header中需要制定默认值的项
# 在下述dict中直接编辑value项即可改变输入值
# 更改此文件仅会更改默认值
# 如需更改基本的项名称或备注 请更改对应的HA(FE)_header.txt
static_header_items = [
    {
        'key': 'TELESCOP',
        'value': 'CHASE/HIS'
    },
    {
        'key': 'CUNIT1',
        'value': 'arcsec'
    },
    {
        'key': 'CDELT1',
        'value': 1.04
    },
    {
        'key': 'CRVAL1',
        'value': 0.00000000
    },
    {
        'key': 'CUNIT2',
        'value': 'arcsec'
    },
    {
        'key': 'CDELT2',
        'value': 1.04
    },
    {
        'key': 'CRVAL2',
        'value': 0.00000000
    },
    {
        'key': 'CRPIX3',
        'value': 0.00000000
    },
    {
        'key': 'CUNIT3',
        'value': 'angstrom'
    },
    {
        'key': 'CDELT3',
        'value': 0.0484
    },
    {
        'key': 'STEPTIME',
        'value': 0.01
    },
    {
        'key': 'SPECLINE',
        'value': 'HA'  # FE
    },
    {
        'key': 'BZERO',
        'value': 32768
    },
    {
        'key': 'BSCALE',
        'value': 1
    }
]

# 直接复制粘贴的项
# 一般选取经过日心的文件
copy_header_items = [
    {
        'key': 'EXP_TIME',
        'typeOfValue': 'number'
    },
    {
        'key': 'SAT_POS1',
        'typeOfValue': 'number'
    },
    {
        'key': 'SAT_POS2',
        'typeOfValue': 'number'
    },
    {
        'key': 'SAT_POS3',
        'typeOfValue': 'number'
    },
    {
        'key': 'SAT_VEL1',
        'typeOfValue': 'number'
    },
    {
        'key': 'SAT_VEL2',
        'typeOfValue': 'number'
    },
    {
        'key': 'SAT_VEL3',
        'typeOfValue': 'number'
    }
]

calculate_header_items = [
    {
        'key': 'BIN',
        'typeOfValue': 'number'
    },
    {
        'key': 'DATE_OBS',
        'typeOfValue': 'string'
    },
    {
        'key': 'CRPIX1',  # TODO
        'typeOfValue': 'number'
    },
    {
        'key': 'CRPIX2',  # TODO
        'typeOfValue': 'number'
    },
    {
        'key': 'PRODATE',
        'typeOfValue': 'string'
    },
    {
        'key': 'STR_TIME',
        'typeOfValue': 'string'
    },
    {
        'key': 'END_TIME',
        'typeOfValue': 'string'
    },
    {
        'key': 'INST_ROT',
        'typeOfValue': 'number'
    },
    {
        'key': 'R_SUN',  # TODO
        'typeOfValue': 'number'
    },
    {
        'key': 'RSUN_OBS',  # TODO
        'typeOfValue': 'number'
    },
    {
        'key': 'B0',
        'typeOfValue': 'number'
    },
    {
        'key': 'FRM_NUM',
        'typeOfValue': 'string'
    }
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
        tempDic['value'] = l_split_first[1].strip()
        tempDic['comment'] = l_split_second[1].strip()
        tempList.append(tempDic)
    return tempList


def main():
    file_data = fits.open(config.data_dir_path + '/' + 'RSM20211222T060105-0008-0001.fts')
    h = file_data[0].header
    print(repr(h))
    q0 = h['Q0']
    q1 = h['Q1']
    q2 = h['Q2']
    q3 = h['Q3']
    strtime = h['STR_TIME']
    B0, INST = suntools.getB0P0(q0, q1, q2, q3, strtime)
    print(B0)

    # b = fits.header.Header()
    # temp = read_header_from_txt('D:\\QQ\\WechatTemp\\WeChat Files\\wxid_vch2mdk54tie22\\FileStorage\\File\\2022-03'
    #                             '\\HA_header.txt')
    # for h in temp:
    #     b.set(h['key'], value=0, comment=h['comment'])
    # primaryHDU = fits.PrimaryHDU(data=file_data[0].data, header=b)
    # greyHDU = fits.HDUList([primaryHDU])
    # greyHDU.writeto(config.data_dir_path + '/' + 'RSM20220120T062139-0013-1244_test.fts', overwrite=True)
    # file_data.close()
    # file_data = fits.open(config.data_dir_path + '/' + 'RSM20220120T062139-0013-1244_test.fts')
    # a = file_data[0].header
    # print(repr(a))


if __name__ == '__main__':
    main()
