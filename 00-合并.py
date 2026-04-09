import os
import pandas as pd

root_dir = "/Users/lemon/Library/Containers/com.tencent.xinWeChat/Data/Documents/xwechat_files/wxid_ont20j4kijbz22_4927/msg/file/2026-04/Origin"

output_A = "Bio/RIPs.csv"
output_B = "Bio/selected_ori_regions.csv"

# 确保输出目录存在
os.makedirs("Bio", exist_ok=True)

written_A = False
written_B = False

for subdir, _, files in os.walk(root_dir):
    for file in files:
        if file.endswith(".csv"):
            file_path = os.path.join(subdir, file)
            print("处理：", file_path)

            try:
                # 👇 分块读取（核心优化）
                for chunk in pd.read_csv(file_path, chunksize=20000):

                    # 可选：记录来源文件（做数据追踪很有用）
                    # chunk["source_file"] = file

                    if "RIP" in file:
                        chunk.to_csv(
                            output_A,
                            mode="a",
                            index=False,
                            header=not written_A
                        )
                        written_A = True

                    elif "selected_ori_regions" in file:
                        chunk.to_csv(
                            output_B,
                            mode="a",
                            index=False,
                            header=not written_B
                        )
                        written_B = True

            except Exception as e:
                print(f"❌ 读取失败：{file_path} -> {e}")

print("✅ 完成！")