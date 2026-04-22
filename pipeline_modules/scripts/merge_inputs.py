import argparse
import os

import pandas as pd


RIP_COLUMNS = [
    "Accession_Number",
    "Sequence_length",
    "gene_start",
    "gene_end",
    "gene_strand",
    "translation",
    "gene_order",
    "gene_id",
    "gene",
    "product",
    "Gene_ID",
    "mmseqs_hit",
    "mmseqs_Identity",
    "mmseqs_Alignmentlength",
    "mmseqs_Querylength",
    "mmseqs_Subjectlength",
    "mmseqs_Evalue",
    "mmseqs_Bitscore",
    "RIP",
]

ORI_COLUMNS = [
    "Accession_Number",
    "Intergenic_Start",
    "Intergenic_End",
    "Evidence",
    "Intergenic_Sequence",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Merge RIP and Ori-region CSV files.")
    parser.add_argument("--root-dir", default="D:/dataprocess/OriVresult", help="Input root directory.")
    parser.add_argument("--output-rip", default="data/RIPs.csv", help="Output RIP CSV path.")
    parser.add_argument(
        "--output-ori",
        default="data/selected_ori_regions.csv",
        help="Output selected_ori_regions CSV path.",
    )
    parser.add_argument(
        "--chunksize",
        "--chunk-size",
        dest="chunksize",
        type=int,
        default=20000,
        help="CSV chunk size.",
    )
    parser.add_argument(
        "--append",
        action="store_true",
        help="Append to existing output files instead of replacing them.",
    )
    return parser.parse_args()


def normalize_chunk(chunk, columns):
    col_map = {str(col).strip().lstrip("\ufeff"): col for col in chunk.columns}

    if all(col in col_map for col in columns):
        normalized = chunk[[col_map[col] for col in columns]].copy()
        normalized.columns = columns
    else:
        normalized = chunk.iloc[:, : len(columns)].copy()
        normalized.columns = columns[: normalized.shape[1]]
        for col in columns[normalized.shape[1] :]:
            normalized[col] = pd.NA
        normalized = normalized[columns]

    header_rows = normalized["Accession_Number"].astype(str).str.strip().eq("Accession_Number")
    return normalized.loc[~header_rows]


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
                print(f"Removed old output: {output_file}")

    written_A = args.append and os.path.exists(output_A) and os.path.getsize(output_A) > 0
    written_B = args.append and os.path.exists(output_B) and os.path.getsize(output_B) > 0
    output_paths = {os.path.abspath(output_A), os.path.abspath(output_B)}

    rip_rows = 0
    ori_rows = 0

    for subdir, dirs, files in os.walk(root_dir):
        dirs.sort()
        for file in sorted(files):
            if not file.endswith(".csv"):
                continue

            file_path = os.path.join(subdir, file)
            if os.path.abspath(file_path) in output_paths:
                continue

            print("Processing:", file_path)

            try:
                for chunk in pd.read_csv(file_path, chunksize=args.chunksize):
                    if "RIP" in file:
                        chunk = normalize_chunk(chunk, RIP_COLUMNS)
                        if chunk.empty:
                            continue
                        chunk.to_csv(output_A, mode="a", index=False, header=not written_A)
                        written_A = True
                        rip_rows += len(chunk)

                    elif "selected_ori_regions" in file:
                        chunk = normalize_chunk(chunk, ORI_COLUMNS)
                        if chunk.empty:
                            continue
                        chunk.to_csv(output_B, mode="a", index=False, header=not written_B)
                        written_B = True
                        ori_rows += len(chunk)

            except Exception as e:
                print(f"Read failed: {file_path} -> {e}")

    print("Done.")
    print(f"RIP rows written: {rip_rows}")
    print(f"Ori rows written: {ori_rows}")


if __name__ == "__main__":
    main()
