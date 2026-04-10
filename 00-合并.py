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
    return parser.parse_args()


def main():
    args = parse_args()
    root_dir = args.root_dir
    output_A = args.output_rip
    output_B = args.output_ori

    os.makedirs(os.path.dirname(output_A) or ".", exist_ok=True)
    os.makedirs(os.path.dirname(output_B) or ".", exist_ok=True)

    written_A = False
    written_B = False

    for subdir, _, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".csv"):
                file_path = os.path.join(subdir, file)
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