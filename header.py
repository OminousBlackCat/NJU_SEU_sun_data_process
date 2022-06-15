import config
from astropy.io import fits
import suntools

'''
此文件代码主要目的是读取data/header.txt中的头部项
并将这些项构造成为一个astropy的header对象
'''

'''
read_header_from_text()
此方法用来读取txt 并将其构造成为一个dict并返回
所有的value都先会尝试使用int构造
再使用float构造: 如果是浮点数 请在txt内表明小数点 例如:5.000
如果都失败才会使用String构造: 不要使用纯数字的字符串值
'''


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
