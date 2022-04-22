import config
import os
if config.bin_count == 1:
    limit_size = config.sun_row_count_bin_1
else:
    limit_size = config.sun_row_count_bin_2
def getTarget(filePath):
    dbtype_list = os.listdir(filePath)
    target_list = []
    out_list = []
    if len(dbtype_list) >= limit_size:
        target_list.append(filePath)
        nums = filePath.split('/')
        out_list .append( nums[-1].split('-')[0])
    err = 0
    for dbtype in dbtype_list:
        try:
            if os.path.isdir(filePath + dbtype):
                filelist = os.listdir(filePath + dbtype)
            else:
                filelist = []
        except IOError:
            err += 1
        else:
            if (len(filelist)) >= limit_size:
                target_list.append(filePath + dbtype)
                nums = dbtype.split('/')
                out_list .append( nums[-1].split('-')[0])

    print(target_list)
    print(out_list)
    return target_list,out_list

if __name__ == "__main__":
    #getTarget('D:/NJU_SEU_sun_data_process/')
    filePath = 'D:/NJU_SEU_sun_data_process/3-23'
    nums = filePath.split('/')
    out_list = nums[-1].split('-')[0]
    print(out_list)

