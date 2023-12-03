"""
此文件用来保存与生成视频相关的参数
此文件的部分参数仅对生成视频脚本生效

@auther: seu_wxy
"""

# 目标png文件所在的文件夹 绝对路径
# 结尾请带"/"
# 示例: /data/chase/Chase/Lev1/2023/11/23/
path_to_target_png = "/"

# 存放视频的文件夹 绝对路径
# 结尾请记得带"/"
# 示例: /data/chase/Chase/Lev1/2023/11/movie/
path_to_video_save = "/"

# 帧大小, 元组(width, height)
# 均为int
# 此参数会显著影响生成视频的大小
# 默认为 (1000, 1000)
video_frame_size = (1000, 1000)

# FPS, 视频每秒包含的帧
# int
# 此参数会显著影响生成视频的时长
# 默认为 6
video_frame_per_sec = 6

# png写入帧的频率
# 默认为1 即所有png都写入视频
# 若为n 则所有n的倍数的png会写入视频
target_png_interval = 2