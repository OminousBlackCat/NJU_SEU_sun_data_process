# 检查输出文件的文件合法性
# 主要检查相邻序列号是否一致
# def check_outputs(outputs_dir: str):
#     arr = os.listdir(outputs_dir)
#     fits_HA_list = []
#     fits_FE_list = []
#     delete_HA_list = []
#     delete_FE_list = []
#     for filename in arr:
#         if filename.split('.')[-1] == 'fits':
#             if filename.split('_')[-1].split('.')[0] == 'HA':
#                 fits_HA_list.append(filename)
#             if filename.split('_')[-1].split('.')[0] == 'FE':
#                 fits_FE_list.append(filename)
#     fits_HA_list.sort()
#     fits_FE_list.sort()
#     # 遍历fits序列 检查时间差
#     log("正在检查输出文件命名有效性...")
#     # 首先排除相同时间的文件
#     for i in range(len(fits_HA_list)):
#         if i > 0:
#             # 示例文件名:
#             # 012 3456 7890 123 45 67 8 9012 3 4567890
#             # RSM 2023 0411 T04 51 06 _ 0001 _ FE.fits
#             current_datetime = datetime.datetime(year=1, month=1, day=1,
#                                                  hour=int(fits_HA_list[i].split('_')[0][12: 14]),
#                                                  minute=int(fits_HA_list[i].split('_')[0][14: 16]),
#                                                  second=int(fits_HA_list[i].split('_')[0][16: 18]))
#             upper_datetime = datetime.datetime(year=1, month=1, day=1,
#                                                hour=int(fits_HA_list[i - 1].split('_')[0][12: 14]),
#                                                minute=int(fits_HA_list[i - 1].split('_')[0][14: 16]),
#                                                second=int(fits_HA_list[i - 1].split('_')[0][16: 18]))
#             # 如果时间差距为0, 则说明两个扫描序列重复出现, 则此时删除较小的文件, 保留较大的文件
#             if (current_datetime - upper_datetime).seconds == 0:
#                 log("文件:" + fits_HA_list[i - 1] + "与文件" + fits_HA_list[i] + "时间一致, 将剔除大小较小的...")
#                 current_HA_file_size = os.path.getsize(os.path.join(outputs_dir, fits_HA_list[i]))
#                 upper_HA_file_size = os.path.getsize(os.path.join(outputs_dir, fits_HA_list[i - 1]))
#                 # 当前 新轨道的 文件较大, 则让新轨道的文件直接和上一个序列文件名一样
#                 if current_HA_file_size > upper_HA_file_size:
#                     log("删除文件:" + fits_HA_list[i - 1] + "!")
#                     delete_HA_list.append(fits_HA_list[i - 1])
#                     delete_FE_list.append(fits_FE_list[i - 1])
#                     # os.remove(os.path.join(outputs_dir, fits_HA_list[i - 1]))
#                     # os.remove(os.path.join(outputs_dir, fits_FE_list[i - 1]))
#                     # os.rename(os.path.join(outputs_dir, fits_HA_list[i]),
#                     #           os.path.join(outputs_dir, fits_HA_list[i - 1]))
#                     # os.rename(os.path.join(outputs_dir, fits_FE_list[i]),
#                     #           os.path.join(outputs_dir, fits_FE_list[i - 1]))
#                 # 上轨道的文件较大或相等, 则直接删除新轨道文件即可
#                 else:
#                     log("删除文件:" + fits_HA_list[i] + "!")
#                     delete_HA_list.append(fits_HA_list[i])
#                     delete_FE_list.append(fits_FE_list[i])
#                     # os.remove(os.path.join(outputs_dir, fits_HA_list[i]))
#                     # os.remove(os.path.join(outputs_dir, fits_FE_list[i]))
#     # 删除delete list 内的文件
#     for f in delete_HA_list:
#         fits_HA_list.remove(f)
#     for f in delete_FE_list:
#         fits_FE_list.remove(f)
#     # 再开始判断是否有轨道错误
#     for i in range(len(fits_HA_list)):
#         try:
#             if i > 0:
#                 # 开启了一个新的扫描轨道
#                 if int(fits_HA_list[i].split('_')[-2]) - int(fits_HA_list[i - 1].split('_')[-2]) != 1:
#                     # 示例文件名:
#                     # 012 3456 7890 123 45 67 8 9012 3 4567890
#                     # RSM 2023 0411 T04 51 06 _ 0001 _ FE.fits
#                     current_datetime = datetime.datetime(year=1, month=1, day=1,
#                                                          hour=int(fits_HA_list[i].split('_')[0][12: 14]),
#                                                          minute=int(fits_HA_list[i].split('_')[0][14: 16]),
#                                                          second=int(fits_HA_list[i].split('_')[0][16: 18]))
#                     upper_datetime = datetime.datetime(year=1, month=1, day=1,
#                                                        hour=int(fits_HA_list[i - 1].split('_')[0][12: 14]),
#                                                        minute=int(fits_HA_list[i - 1].split('_')[0][14: 16]),
#                                                        second=int(fits_HA_list[i - 1].split('_')[0][16: 18]))
#                     # 如果当前时间差小于规定时间 则应该和上一序列处于同一轨道
#                     if (current_datetime - upper_datetime).seconds < config.validate_time_period:
#                         old_HA_name = fits_HA_list[i]
#                         old_FE_name = fits_FE_list[i]
#                         old_HA_png = old_HA_name.split('.')[0] + '.png'
#                         old_FE_png = old_FE_name.split('.')[0] + '.png'
#                         fits_HA_list[i] = fits_HA_list[i].split('_')[0] + "_" + str(
#                             int(fits_HA_list[i - 1][19: 23]) + 1).zfill(4) + "_" + fits_HA_list[i].split('_')[-1]
#                         fits_FE_list[i] = fits_FE_list[i].split('_')[0] + "_" + str(
#                             int(fits_FE_list[i - 1][19: 23]) + 1).zfill(4) + "_" + fits_FE_list[i].split('_')[-1]
#                         new_HA_png = fits_HA_list[i].split('.')[0] + '.png'
#                         new_FE_png = fits_FE_list[i].split('.')[0] + '.png'
#                         log("文件:" + old_HA_name + "应处于上个扫描轨道, 重命名为:" + fits_HA_list[i])
#                         log("文件:" + old_FE_name + "应处于上个扫描轨道, 重命名为:" + fits_FE_list[i])
#                         # os.rename(os.path.join(outputs_dir, old_HA_name),
#                         #                         #           os.path.join(outputs_dir, fits_HA_list[i]))
#                         #                         # os.rename(os.path.join(outputs_dir, old_FE_name),
#                         #                         #           os.path.join(outputs_dir, fits_FE_list[i]))
#                         #                         # os.rename(os.path.join(outputs_dir, old_HA_png), os.path.join(outputs_dir, new_HA_png))
#                         #                         # os.rename(os.path.join(outputs_dir, old_FE_png), os.path.join(outputs_dir, new_FE_png))
#         except:
#             log(traceback.print_exc())
#             log("修改文件名失败, 已跳过")

def test():
    matplotlib.rcParams['font.sans-serif'] = ['KaiTi']
    filepath_result = "testResult/"
    filepath_test = "testData/"
    filepath_bash = "bass2000.txt"

    # print(base)
    image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    dark_data = np.array(fits.getdata(image_file), dtype=float)
    dark_data = change(dark_data)
    image_file = get_pkg_data_filename(filepath_test + 'for_flat_binning2.fits')
    flat_data = np.array(fits.getdata(image_file), dtype=float)
    # RSM20211222T215254-0010-2313-基准.fts    RSM20211222T215555-0013-2367-测试.fts
    # RSM20220120T062536-0017-1081.fts
    H, W = flat_data.shape
    print(H, W)
    print(bin_count)
    flat_data, b, d = curve_correction(flat_data - dark_data, x0, C)
    flat_data = getFlat(flat_data)
    filename = filepath_test + 'RSM20221002T135645-0000-1855.fts'
    image_file = get_pkg_data_filename(filename)
    imgData = np.array(fits.getdata(image_file), dtype=float)
    print("pic0")
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    imgData = moveImg(imgData, -1)
    imgData, HofHa, HofFe = curve_correction(imgData - dark_data, x0, C)
    print("pic1")
    plt.figure()
    plt.imshow(imgData, cmap="gray", aspect='auto')
    plt.show()
    # print(HofHa, HofFe)
    imgData = DivFlat(imgData, flat_data)
    base = get_Sunstd(filepath_bash)
    filepathHA = "HA_absorption.txt"
    filepathFE = "FE_absorption.txt"
    abortion = get_Absorstd(filepathHA, filepathFE, HofHa, HofFe)
    plt.figure()
    plt.plot(base)
    plt.show()
    plt.figure()
    plt.imshow(flat_data, cmap="gray", aspect='auto')
    plt.show()
    # filelist = os.listdir(filepath_test)
    # image_file, imgData = entireWork(filepath_test + 'RSM20221002T135645-0000-1855.fts', dark_data, flat_data, abortion)
    # #
    # print("OK")
    # plt.figure()
    # plt.imshow(imgData, cmap="gray", aspect='auto')
    # plt.show()

    # grey = fits.PrimaryHDU(image_file)

if __name__ == "__main__":
    test()

    # filepath_test = "testData/"
    # image_file = get_pkg_data_filename(filepath_test + 'dark.fits')
    # dark_data = np.array(fits.getdata(image_file), dtype=float)
    # dark_data = change(dark_data)
    # image_file = get_pkg_data_filename(filepath_test + 'for_flat_binning2.fits')
    # flat_data = np.array(fits.getdata(image_file), dtype=float)
    # print(flat_data.shape)
    # flat_data, b, d = curve_correction(flat_data - dark_data, x0, C)
    # flat_data = getFlat(flat_data)
    # flat_data = FlatNormalization(flat_data)
    # primaryHDU = fits.PrimaryHDU(data=flat_data)
    # greyHDU = fits.HDUList([primaryHDU])
    # greyHDU.writeto('FLAT.fts', overwrite=True)
    # plt.figure()
    # plt.imshow(flat_data, cmap="gray", aspect='auto')
    # plt.show()
    # print(min(min(row) for row in flat_data))
    # print(max(max(row) for row in flat_data))

    # print(FlatNormalization(np.array([[0.66,3],[6,9]])))
    # height_ha = int(height_Ha / bin_count) - int(24 / bin_count)
    # height_fe = int(height_Fe / bin_count) - int(24 / bin_count)
    # print(height_ha)
    # A = config.FE_start - config.HA_start
    # B = A / K
    # C = (5430 - 5096) / bin_count
    # print(A, B, C)
    # A = height_Ha
    # B = 5430 - 5096
    # print(A, B)
    # Q = np.array([1, 2, 3, 4, 5, 6, 7])
    # D = 3
    # print(Q[3:int(D) + 3])
    # print(cal_center_mean(np.zeros((164, 3, 3))))
    # cal_center_mean(np.array([[[1,2],[3,4]],[[5,6],[7,8]]]))
    # test()
    # I = Image.open("123.png")

    # I_array = np.array(I.convert('L'))
    # H, W = I_array.shape
    # text1 = str(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    # I_array = add_time(I_array, text1)
    # plt.figure()
    # plt.imshow(I_array, cmap="gray")
    # plt.show()
    # testPath = "circle/circle/"
    # type = "check"
    # if type == "test":
    #     Filelist = os.listdir(testPath)
    #     if True:
    #         id = 105
    #         Filelist = os.listdir(testPath)
    #         I = Image.open(testPath + Filelist[1 + id])
    #         I_array = np.array(I.convert('L'))
    #
    #         # image_file = get_pkg_data_filename(testPath + 'sum8.fts')
    #         # I_array = np.array(fits.getdata(image_file), dtype=float)
    #         # # print(np.shape(I_array))
    #         rx, ry, r = getCircle(I_array)
    #         print(id, rx, ry, r)
    #         H, W = I_array.shape
    #         for i in range(H):
    #             for j in range(W):
    #                 if abs((i - rx) * (i - rx) + (j - ry) * (j - ry) - r * r) < 10000:
    #                     I_array[i][j] = 240
    #         # for point in points:
    #         #     for i in range(20):
    #         #         for j in range(20):
    #         #             I_array[point[0]+i-8][point[1]+j-8] = 240
    #         # print(rx, ry, r * 0.52)
    #         plt.figure()
    #         plt.imshow(I_array)
    #         plt.show()
    # H,W = I_array.shape
    # print(rx,ry,r)
    # for i in range(H):
    #     for j in range(W):
    #         if abs((i-rx)*(i-rx) + (j-ry)*(j-ry) -r*r) <10000:
    #             I_array[i][j]=240
    # plt.figure()
    # plt.imshow(I_array)
    # plt.show()
    # if_first_print = True
    # for i in range(100):
    #     remaining_count = mp.Value('i', int(i))
    #     if if_first_print:
    #         print('当前进度:' + str(remaining_count.value) + '/' + str(file_count.value), end='')
    #         if_first_print = False
    #     else:
    #         print('\b' * (5 + len(str(remaining_count)) + 1 + len(str(file_count.value))) + '当前进度:' + str(
    #             remaining_count.value) + '/' + str(file_count.value), end='')
    #     time.sleep(0.5)
    # test()
    # if type == "check":
    #     Filelist = os.listdir(testPath)
    #     print(Filelist)
    #     L = len(Filelist) - 1
    #     err = []
    #     id = 0
    #     while id < L:
    #         I = Image.open(testPath + Filelist[id + 1])
    #         I_array = np.array(I.convert('L'))
    #
    #         # image_file = get_pkg_data_filename(testPath + 'sum8.fts')
    #         # I_array = np.array(fits.getdata(image_file), dtype=float)
    #         # # print(np.shape(I_array))
    #         ry, rx, r = getCircle(I_array, id)
    #         print(id, rx, ry, r)
    #         H, W = I_array.shape
    #         for i in range(H):
    #             for j in range(W):
    #                 if abs((i - rx) * (i - rx) + (j - ry) * (j - ry) - r * r) < 10000 / bin_count:
    #                     I_array[i][j] = 240
    #         plt.imsave("Result/result/" + str(id) + ".jpg", I_array)
    #         # for point in points:
    #         #     for i in range(20):
    #         #         for j in range(20):
    #         #             I_array[point[0]+i-8][point[1]+j-8] = 240
    #         # print(rx, ry, r * 0.52)=
    #         id += 1
    # test()