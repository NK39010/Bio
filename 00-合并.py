import os
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description="Merge RIP and Ori-region CSV files.")
    parser.add_argument("--root-dir", default="D:/dataprocess/OriVresult", help="Input root directory.")
    parser.add_argument("--output-rip", default="data/RIPs.csv", help="Output RIP CSV path.")
    parser.add_argument(
        "--output-ori",
        default="data/selected_ori_regions.csv",
        help="Output selected_ori_regions CSV path.",
    )
    parser.add_argument("--chunksize", type=int, default=20000, help="CSV chunk size.")
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to existing output files instead of replacing them.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    root_dir = args.root_dir
    output_A = args.output_rip
    output_B = args.output_ori

    os.makedirs(os.path.dirname(output_A) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(output_B) or ".", exist_ok=True)

    if not args.append:
        for output_file in (output_A, output_B):
            if os.path.exists(output_file):
                os.remove(output_file)
                print(f"已删除旧输出：{output_file}")

    written_A = args.append and os.path.exists(output_A) and os.path.getsize(output_A) > 0
    written_B = args.append and os.path.exists(output_B) and os.path.getsize(output_B) > 0
    output_paths = {os.path.abspath(output_A), os.path.abspath(output_B)}

    for subdir, dirs, files in os.walk(root_dir):
        dirs.sort()
        for file in sorted(files):
            if file.endswith(".csv"):
                file_path = os.path.join(subdir, file)
                if os.path.abspath(file_path) in output_paths:
                    continue

                print("处理：", file_path)

                try:
                    for chunk in pd.read_csv(file_path, chunksize=args.chunksize):
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


if __name__ == "__main__":
    main()
