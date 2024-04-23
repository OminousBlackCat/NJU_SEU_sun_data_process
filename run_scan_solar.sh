echo "正在激活Python 3.10.1虚拟环境中..."
CURRENT_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
echo "当前工作路径为:"
echo "$CURRENT_DIR"
export CURRENT_DIR
source bin/activate
echo "确保config.py已经经过编辑..."
echo "运行程序中..."
python swing_scan_multiprocess_up.py
echo "此轮程序运行完毕!"
deactivate
echo "退出程序中..."
