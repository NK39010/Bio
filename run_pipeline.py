from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from pathlib import Path

def run_cmd(cmd: list[str]) -> None:
    print(f"\n>>> Running: {' '.join(shlex.quote(x) for x in cmd)}")
    subprocess.run(cmd, check=True)


def run_python(script: str, args: list[str] | None = None) -> None:
    cmd = [sys.executable, script]
    if args:
        cmd.extend(args)
    run_cmd(cmd)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="End-to-end Ori/Rep/CD-HIT pipeline runner.")
    parser.add_argument("--root-dir", required=True, help="Root folder for step-00 merge input CSVs.")
    parser.add_argument("--cdhit-mode", choices=["mac", "wsl"], default="mac", help="How to run cd-hit-est.")
    parser.add_argument("--cdhit-c", type=float, default=0.9, help="CD-HIT identity threshold.")
    parser.add_argument("--cdhit-n", type=int, default=5, help="CD-HIT word size.")
    parser.add_argument("--cdhit-t", type=int, default=8, help="CD-HIT threads.")
    parser.add_argument("--cdhit-m", type=int, default=0, help="CD-HIT memory limit MB (0 means unlimited).")
    parser.add_argument(
        "--analysis-first",
        action="store_true",
        help="Deprecated. Pre-CD-HIT analysis now runs by default unless skipped.",
    )
    parser.add_argument("--skip-00", action="store_true")
    parser.add_argument("--skip-01", action="store_true")
    parser.add_argument("--skip-02", action="store_true")
    parser.add_argument("--skip-03", action="store_true")
    parser.add_argument("--skip-annotate", action="store_true")
    parser.add_argument("--skip-04", action="store_true")
    parser.add_argument("--skip-05", action="store_true")
    parser.add_argument("--skip-06", action="store_true")
    parser.add_argument("--skip-07", action="store_true")
    parser.add_argument("--skip-08", action="store_true")
    parser.add_argument("--skip-09", action="store_true")
    parser.add_argument("--kmer", type=int, default=6, help="K-mer size for token output.")
    parser.add_argument("--final-csv", default="final_data.csv", help="Final filtered CSV output path.")
    parser.add_argument("--final-parquet", default="final_data.parquet", help="Final parquet path.")
    parser.add_argument("--token-json", default="final_data_tokens.json", help="Token JSON output path.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    merged_csv = "merged_final_optimized.csv"
    annotated_csv = "merged_final_optimized_replication_annotated.csv"
    ori_fasta = "ori.fasta"
    cdhit_fasta = "cd_hit_process.fasta"
    cdhit_clstr = f"{cdhit_fasta}.clstr"

    if not args.skip_00:
        run_python(
            "00-合并.py",
            ["--root-dir", args.root_dir, "--output-rip", "data/RIPs.csv", "--output-ori", "data/selected_ori_regions.csv"],
        )

    if not args.skip_01:
        run_python("01-OriV&RIFs-mix.py")

    if not args.skip_02:
        run_python(
            "02-种属分布.py",
            ["--input", merged_csv, "--output-excel", "species_before.xlsx", "--output-image", "species_before.png"],
        )
    if not args.skip_03:
        run_python("03-长度分布.py", ["--input", merged_csv, "--output-image", "length_before.png"])

    if not args.skip_annotate:
        run_python("02-annotate_replication_mechanism.py", ["--input", merged_csv, "--output", annotated_csv])

    if not args.skip_04:
        run_python("04-cd-hit序列构建.py", ["--input", annotated_csv, "--output", ori_fasta])

    if not args.skip_05:
        cdhit_cmd = [
            "cd-hit-est",
            "-i",
            ori_fasta,
            "-o",
            cdhit_fasta,
            "-c",
            str(args.cdhit_c),
            "-n",
            str(args.cdhit_n),
            "-M",
            str(args.cdhit_m),
            "-T",
            str(args.cdhit_t),
            "-d",
            "0",
        ]
        if args.cdhit_mode == "wsl":
            cdhit_cmd = ["wsl"] + cdhit_cmd
        run_cmd(cdhit_cmd)

    if not args.skip_06:
        run_python("06-cd-hit去冗余结果可视化.py", ["--clstr", cdhit_clstr, "--output", "cdhit_redundancy.png"])

    if not args.skip_07:
        run_python(
            "07-cd-hit序列返还csv构建.py",
            ["--fasta", cdhit_fasta, "--csv", annotated_csv, "--output", args.final_csv, "--id-column", "ori_id"],
        )

    if not args.skip_02:
        run_python(
            "02-种属分布.py",
            ["--input", args.final_csv, "--output-excel", "species_after.xlsx", "--output-image", "species_after.png"],
        )
    if not args.skip_03:
        run_python("03-长度分布.py", ["--input", args.final_csv, "--output-image", "length_after.png"])

    token_base_path = Path(args.token_json).with_name(
        f"{Path(args.token_json).stem}_without_mechanism{Path(args.token_json).suffix}"
    )
    token_with_mech_path = Path(args.token_json).with_name(
        f"{Path(args.token_json).stem}_with_mechanism{Path(args.token_json).suffix}"
    )

    if not args.skip_08:
        run_python("08-构建训练用的parquet.py", ["--input", args.final_csv, "--parquet", args.final_parquet])

    if not args.skip_09:
        run_python(
            "09-构建词元.py",
            [
                "--input",
                args.final_csv,
                "--token-base-json",
                str(token_base_path),
                "--token-with-mechanism-json",
                str(token_with_mech_path),
                "--kmer",
                str(args.kmer),
            ],
        )

    print("\n✅ Pipeline finished.")
    print(f"Final CSV: {args.final_csv}")
    print(f"Final Parquet: {args.final_parquet}")
    print(f"Token JSON (without mechanism): {token_base_path}")
    print(f"Token JSON (with mechanism): {token_with_mech_path}")


if __name__ == "__main__":
    main()
