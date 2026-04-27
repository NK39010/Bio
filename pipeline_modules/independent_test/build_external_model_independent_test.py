from __future__ import annotations

"""Build an independent test set for a model trained on external train/validation.

Training reference:
  - train-00000-of-00001 (1).parquet
  - validation-00000-of-00001 (1).parquet

Independent-test candidates:
  - output-lab/train-data/all_data.csv
  - test-00000-of-00001.parquet

The candidates are filtered against the training reference by exact oriV/plasmid
matches and BLASTn similarity, using the OriGen-style default thresholds:
pident > 95 and qcovs > 95.
"""

import argparse
import logging
import shutil
from pathlib import Path

import pandas as pd

try:
    from .build_independent_test_blastn import (
        STANDARD_COLS,
        apply_internal_dedup,
        build_reference_exact_mask,
        clean_text,
        read_external_tables,
        run_blastn,
        write_parquet,
    )
except ImportError:
    from build_independent_test_blastn import (
        STANDARD_COLS,
        apply_internal_dedup,
        build_reference_exact_mask,
        clean_text,
        read_external_tables,
        run_blastn,
        write_parquet,
    )


logger = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an independent test set for the external train/validation model."
    )
    parser.add_argument(
        "--training-reference-parquet",
        nargs="+",
        default=[
            "train-00000-of-00001 (1).parquet",
            "validation-00000-of-00001 (1).parquet",
        ],
        help="Parquet files used to train/validate the external model; these become the BLASTn reference.",
    )
    parser.add_argument(
        "--candidate-csv",
        nargs="*",
        default=["output-lab/train-data/all_data.csv"],
        help="CSV candidate files to include in the independent-test pool.",
    )
    parser.add_argument(
        "--candidate-parquet",
        nargs="*",
        default=["test-00000-of-00001.parquet"],
        help="Parquet candidate files to include in the independent-test pool.",
    )
    parser.add_argument(
        "--output-dir",
        default="output-lab/external-model-independent-test",
        help="Directory for outputs and audit files.",
    )
    parser.add_argument("--pident-threshold", type=float, default=95.0, help="BLAST pident cutoff.")
    parser.add_argument("--qcovs-threshold", type=float, default=95.0, help="BLAST qcovs cutoff.")
    parser.add_argument(
        "--candidate-internal-dedup-key",
        choices=["none", "row", "oriv", "rep-oriv", "plasmid-oriv"],
        default="row",
        help="How to deduplicate the independent-test candidate pool before filtering.",
    )
    parser.add_argument(
        "--reference-exact-key",
        choices=["oriv", "oriv-plasmid", "strict"],
        default="oriv-plasmid",
        help="Exact filtering against training reference before BLASTn.",
    )
    parser.add_argument("--keep-blast-workdir", action="store_true", help="Keep BLAST intermediate files.")
    parser.add_argument("--log-level", default="INFO", help="Logging level: DEBUG, INFO, WARNING, ERROR.")
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


def load_standard_csv(path: Path) -> pd.DataFrame:
    logger.info("Loading candidate CSV: %s", path)
    df = pd.read_csv(path, low_memory=False)
    missing = [col for col in STANDARD_COLS if col not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing standard columns: {', '.join(missing)}")

    out = df[STANDARD_COLS].copy()
    out["rep_protein"] = out["rep_protein"].apply(lambda x: clean_text(x, upper=True))
    out["oriv_sequence"] = out["oriv_sequence"].apply(lambda x: clean_text(x, lower=True))
    out["host_species"] = out["host_species"].apply(clean_text)
    out["database_source"] = out["database_source"].apply(clean_text)
    out["plasmid_id"] = out["plasmid_id"].apply(clean_text)
    out["split"] = "independent_test"
    out["source_file"] = path.name
    return out


def load_candidate_pool(csv_paths: list[str], parquet_paths: list[str]) -> pd.DataFrame:
    frames = []
    for path_text in csv_paths:
        frames.append(load_standard_csv(Path(path_text)))
    if parquet_paths:
        frames.append(read_external_tables(parquet_paths))
    if not frames:
        raise ValueError("No independent-test candidates were provided.")
    merged = pd.concat(frames, ignore_index=True)
    logger.info("Merged candidate rows: %d", len(merged))
    return merged


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "blast_work"

    logger.info("Loading external train/validation as training reference")
    reference_df = read_external_tables(args.training_reference_parquet)
    reference_df = reference_df[reference_df["oriv_sequence"] != ""].copy()
    raw_reference_rows = len(reference_df)
    reference_df, removed_reference_duplicate_rows = apply_internal_dedup(reference_df, "row")
    logger.info("Training reference rows after exact row dedup: %d", len(reference_df))

    candidates = load_candidate_pool(args.candidate_csv, args.candidate_parquet)
    raw_candidate_rows = len(candidates)
    candidates = candidates[candidates["oriv_sequence"] != ""].copy()
    removed_empty_oriv = raw_candidate_rows - len(candidates)
    candidates, removed_candidate_internal_duplicates = apply_internal_dedup(
        candidates,
        args.candidate_internal_dedup_key,
    )
    logger.info(
        "Candidate rows after internal dedup (%s): %d",
        args.candidate_internal_dedup_key,
        len(candidates),
    )

    exact_mask = build_reference_exact_mask(candidates, reference_df, args.reference_exact_key)
    exact_rejected = candidates[exact_mask].copy()
    filtered_candidates = candidates[~exact_mask].copy().reset_index(drop=True)
    filtered_candidates["candidate_id"] = "cand_" + filtered_candidates.index.astype(str)
    logger.info("Candidate rows after exact filtering against training reference: %d", len(filtered_candidates))

    blast_hits = run_blastn(reference_df, filtered_candidates, work_dir)
    blast_hits.to_csv(output_dir / "blast_hits.tsv", sep="\t", index=False)

    if blast_hits.empty:
        blast_rejected_orivs: set[str] = set()
        high_similarity_hits = blast_hits.copy()
    else:
        high_similarity_hits = blast_hits[
            (blast_hits["pident"] > args.pident_threshold) & (blast_hits["qcovs"] > args.qcovs_threshold)
        ].copy()
        blast_rejected_orivs = set(high_similarity_hits["oriv_sequence"].dropna().astype(str))
    high_similarity_hits.to_csv(output_dir / "blast_rejected_hits.tsv", sep="\t", index=False)

    blast_rejected = filtered_candidates[filtered_candidates["oriv_sequence"].isin(blast_rejected_orivs)].copy()
    independent_df = filtered_candidates[~filtered_candidates["oriv_sequence"].isin(blast_rejected_orivs)].copy()

    audit_cols = STANDARD_COLS + ["source_file"]
    exact_rejected[audit_cols].to_csv(output_dir / "rejected_exact.csv", index=False, encoding="utf-8-sig")
    blast_rejected[audit_cols + ["candidate_id"]].to_csv(
        output_dir / "rejected_blastn.csv",
        index=False,
        encoding="utf-8-sig",
    )

    independent_df = independent_df[audit_cols].copy()
    independent_csv = output_dir / "independent_test.csv"
    independent_parquet = output_dir / "independent_test.parquet"
    independent_df.to_csv(independent_csv, index=False, encoding="utf-8-sig")
    write_parquet(independent_df, independent_parquet)

    report = output_dir / "dedup_report.txt"
    report.write_text(
        "\n".join(
            [
                "External train/validation model independent-test report",
                f"training_reference_parquet: {', '.join(args.training_reference_parquet)}",
                f"candidate_csv: {', '.join(args.candidate_csv)}",
                f"candidate_parquet: {', '.join(args.candidate_parquet)}",
                f"raw_training_reference_rows: {raw_reference_rows}",
                f"removed_training_reference_duplicate_rows: {removed_reference_duplicate_rows}",
                f"raw_candidate_rows: {raw_candidate_rows}",
                f"removed_empty_oriv_candidates: {removed_empty_oriv}",
                f"candidate_internal_dedup_key: {args.candidate_internal_dedup_key}",
                f"removed_candidate_internal_duplicates: {removed_candidate_internal_duplicates}",
                f"reference_exact_key: {args.reference_exact_key}",
                f"removed_exact_against_training_reference: {len(exact_rejected)}",
                f"blastn_filter: pident > {args.pident_threshold} and qcovs > {args.qcovs_threshold}",
                f"removed_blastn_against_training_reference: {len(blast_rejected)}",
                f"final_independent_test_rows: {len(independent_df)}",
                f"independent_csv: {independent_csv}",
                f"independent_parquet: {independent_parquet}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    if not args.keep_blast_workdir and work_dir.exists():
        shutil.rmtree(work_dir)

    logger.info("Independent test parquet: %s", independent_parquet)
    logger.info("Independent test CSV: %s", independent_csv)
    logger.info("Report: %s", report)


if __name__ == "__main__":
    main()
