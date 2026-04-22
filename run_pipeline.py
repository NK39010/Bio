from __future__ import annotations

import argparse
from pathlib import Path, PureWindowsPath

from pipeline_modules.common import run_cmd, run_python_module, run_python_script


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="End-to-end Ori/Rep/CD-HIT pipeline runner.")
    parser.add_argument("--root-dir", required=True, help="Root folder for step-00 merge input CSVs.")
    parser.add_argument("--output-dir", default="output", help="Unified output directory.")
    parser.add_argument("--cache-file", default="taxonomy_cache.pkl", help="Shared cache file in repo root.")

    parser.add_argument(
        "--cdhit-mode",
        choices=["native", "wsl", "mac", "linux"],
        default="native",
        help=(
            "How to run cd-hit-est. "
            "'native' runs directly in the current shell (macOS/Linux). "
            "'wsl' runs via WSL on Windows. "
            "'mac' and 'linux' are backward-compatible aliases of 'native'."
        ),
    )
    parser.add_argument("--cdhit-c", type=float, default=0.9, help="CD-HIT identity threshold.")
    parser.add_argument("--cdhit-n", type=int, default=5, help="CD-HIT word size.")
    parser.add_argument("--cdhit-t", type=int, default=8, help="CD-HIT threads.")
    parser.add_argument("--cdhit-m", type=int, default=0, help="CD-HIT memory limit MB (0 means unlimited).")
    parser.add_argument(
        "--cdhit-dedup-mode",
        choices=["keep-all", "representative-only"],
        default="keep-all",
        help="How step 07 returns CD-HIT results to CSV.",
    )
    parser.add_argument("--validation-ratio", type=float, default=0.25, help="Split fallback ratio in step 08.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for split fallback.")
    parser.add_argument("--species-min-count", type=int, default=6, help="Tokenizer species min count.")
    parser.add_argument(
        "--include-mechanism-column",
        action="store_true",
        help="Include replication mechanism column in generated dataset files.",
    )
    parser.add_argument(
        "--mechanism-column",
        default="replication_mechanism_term",
        help="Source mechanism column name for dataset generation.",
    )
    parser.add_argument(
        "--include-mechanism-tokens",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Whether tokenizer includes mechanism tokens.",
    )

    parser.add_argument(
        "--hmm-db",
        default="pipeline_modules/resources/replication_annotation/rep_plasmid_core_db.hmm",
        help="HMM DB file for annotation.",
    )
    parser.add_argument(
        "--rules",
        default="pipeline_modules/resources/replication_annotation/pfam_replication_mechanism_table.tsv",
        help="Pfam-to-mechanism rules TSV file.",
    )
    parser.add_argument("--cpus", type=int, default=8, help="Threads for pyhmmer.")

    parser.add_argument("--skip-data-processing", action="store_true")
    parser.add_argument("--skip-cdhit", action="store_true")
    parser.add_argument("--skip-dataset-build", action="store_true")
    parser.add_argument("--skip-viz", action="store_true")
    return parser.parse_args()


def wsl_path(path: Path) -> str:
    resolved = PureWindowsPath(path.resolve())
    drive = resolved.drive.rstrip(":").lower()
    parts = "/".join(resolved.parts[1:])
    return f"/mnt/{drive}/{parts}"


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    data_dir = output_dir / "data"
    train_data_dir = output_dir / "train-data"
    tokenizer_dir = output_dir / "tokenizer"
    log_dir = output_dir / "logs"
    viz_dir = output_dir / "visualizations"

    for folder in (data_dir, train_data_dir, tokenizer_dir, log_dir, viz_dir):
        folder.mkdir(parents=True, exist_ok=True)

    merged_csv = data_dir / "merged_final_optimized.csv"
    annotated_csv = data_dir / "merged_final_optimized_replication_annotated.csv"
    final_csv = data_dir / "final_data.csv"
    ori_fasta = data_dir / "ori.fasta"
    cdhit_fasta = data_dir / "cd_hit_process.fasta"
    cdhit_clstr = Path(str(cdhit_fasta) + ".clstr")

    if not args.skip_data_processing:
        run_python_script(
            "pipeline_modules/scripts/merge_inputs.py",
            [
                "--root-dir",
                args.root_dir,
                "--output-rip",
                str(data_dir / "RIPs.csv"),
                "--output-ori",
                str(data_dir / "selected_ori_regions.csv"),
            ],
        )
        run_python_script(
            "pipeline_modules/scripts/pair_oriv_rep.py",
            [
                "--rip-csv",
                str(data_dir / "RIPs.csv"),
                "--ori-csv",
                str(data_dir / "selected_ori_regions.csv"),
                "--output",
                str(merged_csv),
                "--cache-file",
                args.cache_file,
            ],
        )
        run_python_script(
            "pipeline_modules/scripts/annotate_replication_mechanism.py",
            [
                "--input",
                str(merged_csv),
                "--output",
                str(annotated_csv),
                "--hmm-db",
                args.hmm_db,
                "--rules",
                args.rules,
                "--cpus",
                str(args.cpus),
                "--log-file",
                str(log_dir / "02_annotate_replication_mechanism.log"),
            ],
        )
        run_python_script(
            "pipeline_modules/scripts/build_cdhit_fasta.py",
            [
                "--input",
                str(annotated_csv),
                "--output",
                str(ori_fasta),
                "--log-file",
                str(log_dir / "04_build_cdhit_fasta.log"),
            ],
        )

    if not args.skip_viz and merged_csv.exists():
        run_python_module(
            "pipeline_modules.visualization.species_distribution",
            [
                "--input",
                str(merged_csv),
                "--output-excel",
                str(viz_dir / "species_before.xlsx"),
                "--output-image",
                str(viz_dir / "species_before.png"),
            ],
        )
        run_python_module(
            "pipeline_modules.visualization.length_distribution",
            [
                "--input",
                str(merged_csv),
                "--output-image",
                str(viz_dir / "length_before.png"),
            ],
        )

    if not args.skip_cdhit:
        cdhit_input = str(ori_fasta)
        cdhit_output = str(cdhit_fasta)
        if args.cdhit_mode == "wsl":
            cdhit_input = wsl_path(ori_fasta)
            cdhit_output = wsl_path(cdhit_fasta)

        cdhit_cmd = [
            "cd-hit-est",
            "-i",
            cdhit_input,
            "-o",
            cdhit_output,
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

    if not args.skip_viz and cdhit_clstr.exists():
        run_python_module(
            "pipeline_modules.visualization.cdhit_redundancy",
            [
                "--clstr",
                str(cdhit_clstr),
                "--output",
                str(viz_dir / "cdhit_redundancy.png"),
            ],
        )

    if not args.skip_dataset_build:
        run_python_script(
            "pipeline_modules/scripts/restore_cdhit_to_csv.py",
            [
                "--fasta",
                str(cdhit_fasta),
                "--csv",
                str(annotated_csv),
                "--output",
                str(final_csv),
                "--id-column",
                "ori_id",
                "--dedup-mode",
                args.cdhit_dedup_mode,
                "--log-file",
                str(log_dir / "07_return_cdhit_to_csv.log"),
            ],
        )
        dataset_cmd = [
            "--input",
            str(final_csv),
            "--output-dir",
            str(train_data_dir),
            "--validation-ratio",
            str(args.validation_ratio),
            "--seed",
            str(args.seed),
            "--mechanism-column",
            args.mechanism_column,
            "--log-file",
            str(log_dir / "08_build_dataset.log"),
        ]
        if args.include_mechanism_column:
            dataset_cmd.append("--include-mechanism-column")
        run_python_script("pipeline_modules/scripts/build_training_dataset.py", dataset_cmd)
        tokenizer_cmd = [
            "--input",
            str(train_data_dir / "all_data.csv"),
            "--tokenizer-json",
            str(tokenizer_dir / "tokenizer.json"),
            "--species-min-count",
            str(args.species_min_count),
            "--tokenizer-config-json",
            str(tokenizer_dir / "tokenizer_config.json"),
            "--special-tokens-map-json",
            str(tokenizer_dir / "special_tokens_map.json"),
            "--log-file",
            str(log_dir / "09_build_tokenizer.log"),
        ]
        tokenizer_cmd.append("--include-mechanism-tokens" if args.include_mechanism_tokens else "--no-include-mechanism-tokens")
        run_python_script("pipeline_modules/scripts/build_tokenizer.py", tokenizer_cmd)

    if not args.skip_viz and final_csv.exists():
        run_python_module(
            "pipeline_modules.visualization.species_distribution",
            [
                "--input",
                str(final_csv),
                "--output-excel",
                str(viz_dir / "species_after.xlsx"),
                "--output-image",
                str(viz_dir / "species_after.png"),
            ],
        )
        run_python_module(
            "pipeline_modules.visualization.length_distribution",
            [
                "--input",
                str(final_csv),
                "--output-image",
                str(viz_dir / "length_after.png"),
            ],
        )

    print("\nPipeline finished.")
    print(f"Output root: {output_dir}")
    print(f"Data: {data_dir}")
    print(f"Train-data: {train_data_dir}")
    print(f"Tokenizer: {tokenizer_dir}")
    print(f"Logs: {log_dir}")
    print(f"Visualizations: {viz_dir}")
    print(f"Shared cache (root): {args.cache_file}")


if __name__ == "__main__":
    main()
