import config
import os
if config.bin_count == 1:
    limit_size = config.sun_row_count_bin_1
else:
    limit_size = config.sun_row_count_bin_2
def getTarget(filePath):
    dbtype_list = os.listdir(filePath)
    target_list = []
    if len(dbtype_list) >= limit_size:
        target_list.append(filePath)
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

    print(target_list)
    return target_list

if __name__ == "__main__":
    getTarget('D:/NJU_SEU_sun_data_process/')

