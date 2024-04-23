"""
此文件代码主要目的是读取data/header.txt中的头部项
并将这些项构造成为一个astropy的header对象
@author: seu_wxy
"""

import config
from astropy.io import fits

'''
read_header_from_text()
此方法用来读取txt 并将其构造成为一个dict并返回
所有的value都先会尝试使用int构造
再使用float构造: 如果是浮点数 请在txt内表明小数点 例如:5.000
如果都失败才会使用String构造: 不要使用纯数字的字符串值
'''
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
    },
    {
         'key': 'ORID',
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
        true_value = None
        try:
            true_value = int(l_split_second[0].strip())
        except ValueError:
            try:
                true_value = float(l_split_second[0].strip())
            except ValueError:
                true_value = l_split_second[0].replace("'", "")
                true_value = true_value.strip()
        tempDic['value'] = true_value
        tempDic['comment'] = l_split_second[1].strip()
        tempList.append(tempDic)
    return tempList


def main():
    file_data = fits.open('C:\\Users\\seu-wxy\\Desktop\\太阳\\data\\raw_data_test\\RSM20220903T093109_0000_FE.fits')
    h = file_data[1].header
    print(repr(h))

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
