if [[ "$1" == "-h" ]]; then
	echo "此脚本用于激活虚拟环境并执行生成预览太阳像视频"
	echo "执行此脚本请先修改video_config.py "
	exit 0
fi

echo "正在激活Python 3.10.1虚拟环境中..."
CURRENT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
export CURRENT_DIR
source bin/activate
echo "生成视频中..."
python save_png_video.py
deactivate
echo "退出程序中..."

