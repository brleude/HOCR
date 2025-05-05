import json
import os

# 设置工作路径为当前文件所在的目录
current_file_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_file_dir)

# 保存得分
def save_scores(data_dict, file_path):
    """将字典保存到 txt 文件中

    参数：
         data_dict: 要保存的字典
         file_path: 保存的文件路径
    """
    try:
        with open(file_path, 'w') as file:
            json.dump(data_dict, file)
        print(f"已成功保存到 {file_path}")
    except Exception as e:
        print(f"保存时出现错误: {e}")

# 加载得分
def load_scores(file_path):
    """从 txt 文件中加载字典

    参数：
         file_path: 要加载的文件路径
        
    返回：
         要加载的字典，如果出现错误则返回空字典
    """
    try:
        with open(file_path, 'r') as file:
            data_dict = json.load(file)
        return data_dict
    except Exception as e:
        print(f"加载时出现错误: {e}")
        return {}
