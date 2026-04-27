from __future__ import annotations

"""Build an independent test set by BLASTn de-duplicating external oriVs.

This follows the OriGen-style held-out test filtering rule: remove external
candidate oriVs that are highly similar to any training oriV
(default: pident > 95 and qcovs > 95).
"""

import argparse
import logging
import shutil
import subprocess
from pathlib import Path

try:
    import pandas as pd
    import pyarrow as pa
    import pyarrow.parquet as pq
except ModuleNotFoundError:
    pd = None
    pa = None
    pq = None


logger = logging.getLogger(__name__)


STANDARD_COLS = [
    "rep_protein",
    "oriv_sequence",
    "host_species",
    "database_source",
    "plasmid_id",
    "split",
]

EXTERNAL_COL_ALIASES = {
    "rep_protein": ["rep_protein", "rep_seq", "rep_sequence", "protein_sequence", "Rep protein"],
    "oriv_sequence": ["oriv_sequence", "oriV_sequence", "ori_sequence", "OriC sequence", "sequence"],
    "host_species": ["host_species", "species", "host", "organism"],
    "database_source": ["database_source", "source", "db_source"],
    "plasmid_id": ["plasmid_id", "plasmid", "accession", "genome_id", "id"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge external parquet files and build a BLASTn-deduplicated independent test set."
    )
    parser.add_argument(
        "--reference-csv",
        default="output-lab/train-data/all_data.csv",
        help="Your full training/validation CSV used as the BLASTn reference.",
    )
    parser.add_argument(
        "--external-parquet",
        nargs="+",
        default=[
            "test-00000-of-00001.parquet",
            "train-00000-of-00001 (1).parquet",
            "validation-00000-of-00001 (1).parquet",
        ],
        help="External parquet files to merge into candidate independent-test data.",
    )
    parser.add_argument(
        "--output-dir",
        default="output-lab/independent-test",
        help="Directory for independent-test outputs and audit files.",
    )
    parser.add_argument("--pident-threshold", type=float, default=95.0, help="BLAST pident cutoff.")
    parser.add_argument("--qcovs-threshold", type=float, default=95.0, help="BLAST qcovs cutoff.")
    parser.add_argument(
        "--internal-dedup-key",
        choices=["none", "row", "oriv", "rep-oriv", "plasmid-oriv"],
        default="row",
        help=(
            "How to deduplicate merged external files before comparing to your data. "
            "'row' only removes exact duplicate standardized records; 'oriv' is stricter."
        ),
    )
    parser.add_argument(
        "--reference-exact-key",
        choices=["oriv", "oriv-plasmid", "strict"],
        default="oriv-plasmid",
        help=(
            "Exact filtering against your reference before BLASTn. "
            "'oriv-plasmid' removes exact oriV matches and same plasmid_id matches; "
            "'strict' also removes exact rep_protein and rep+oriV matches."
        ),
    )
    parser.add_argument(
        "--keep-blast-workdir",
        action="store_true",
        help="Keep intermediate FASTA, BLAST database, and raw BLAST TSV files.",
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level: DEBUG, INFO, WARNING, ERROR.")
    return parser.parse_args()


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        force=True,
    )


def clean_text(value: object, upper: bool = False, lower: bool = False) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if upper:
        return text.upper()
    if lower:
        return text.lower()
    return text


def find_col(df: pd.DataFrame, aliases: list[str]) -> str | None:
    exact = {col: col for col in df.columns}
    lowered = {str(col).strip().lower(): col for col in df.columns}
    for alias in aliases:
        if alias in exact:
            return exact[alias]
        hit = lowered.get(alias.strip().lower())
        if hit is not None:
            return hit
    return None


def series_or_empty(df: pd.DataFrame, aliases: list[str]) -> pd.Series:
    col = find_col(df, aliases)
    if col is None:
        return pd.Series([""] * len(df), index=df.index, dtype="object")
    return df[col]


def standardize_external(df: pd.DataFrame, source_file: str) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    out["rep_protein"] = series_or_empty(df, EXTERNAL_COL_ALIASES["rep_protein"]).apply(
        lambda x: clean_text(x, upper=True)
    )
    out["oriv_sequence"] = series_or_empty(df, EXTERNAL_COL_ALIASES["oriv_sequence"]).apply(
        lambda x: clean_text(x, lower=True)
    )
    out["host_species"] = series_or_empty(df, EXTERNAL_COL_ALIASES["host_species"]).apply(clean_text)
    source = series_or_empty(df, EXTERNAL_COL_ALIASES["database_source"]).apply(clean_text)
    out["database_source"] = source.where(source != "", "external")
    out["plasmid_id"] = series_or_empty(df, EXTERNAL_COL_ALIASES["plasmid_id"]).apply(clean_text)
    out["split"] = "independent_test"
    out["source_file"] = source_file
    return out


def read_external_tables(paths: list[str]) -> pd.DataFrame:
    frames = []
    for path_text in paths:
        path = Path(path_text)
        logger.info("Loading external parquet: %s", path)
        df = pd.read_parquet(path)
        frames.append(standardize_external(df, path.name))
    if not frames:
        raise ValueError("No external parquet files were provided.")
    merged = pd.concat(frames, ignore_index=True)
    logger.info("Merged external rows: %d", len(merged))
    return merged


def write_fasta(records: pd.DataFrame, seq_col: str, id_col: str, path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in records.itertuples(index=False):
            seq = getattr(row, seq_col)
            if not seq:
                continue
            handle.write(f">{getattr(row, id_col)}\n")
            for start in range(0, len(seq), 80):
                handle.write(seq[start : start + 80] + "\n")


def run_command(command: list[str]) -> None:
    logger.info("Running: %s", " ".join(command))
    subprocess.run(command, check=True)


def apply_internal_dedup(df: pd.DataFrame, mode: str) -> tuple[pd.DataFrame, int]:
    if mode == "none":
        return df.copy(), 0

    before = len(df)
    if mode == "row":
        subset = STANDARD_COLS
    elif mode == "oriv":
        subset = ["oriv_sequence"]
    elif mode == "rep-oriv":
        subset = ["rep_protein", "oriv_sequence"]
    elif mode == "plasmid-oriv":
        subset = ["plasmid_id", "oriv_sequence"]
    else:
        raise ValueError(f"Unsupported internal dedup mode: {mode}")

    out = df.drop_duplicates(subset=subset).copy()
    return out, before - len(out)


def build_reference_exact_mask(
    external_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    mode: str,
) -> pd.Series:
    ref_orivs = set(reference_df["oriv_sequence"].dropna().astype(str).str.strip().str.lower())
    mask = external_df["oriv_sequence"].isin(ref_orivs)

    if mode in {"oriv-plasmid", "strict"} and "plasmid_id" in reference_df.columns:
        ref_plasmids = set(reference_df["plasmid_id"].dropna().astype(str).str.strip())
        mask |= external_df["plasmid_id"].isin(ref_plasmids) & (external_df["plasmid_id"] != "")

    if mode == "strict" and "rep_protein" in reference_df.columns:
        ref_reps = set(reference_df["rep_protein"].dropna().astype(str).str.strip().str.upper())
        ref_pairs = set(zip(reference_df["rep_protein"], reference_df["oriv_sequence"], strict=False))
        mask |= external_df["rep_protein"].isin(ref_reps) & (external_df["rep_protein"] != "")
        mask |= pd.Series(
            list(zip(external_df["rep_protein"], external_df["oriv_sequence"], strict=False)),
            index=external_df.index,
        ).isin(ref_pairs)

    return mask


def run_blastn(reference_df: pd.DataFrame, candidate_df: pd.DataFrame, work_dir: Path) -> pd.DataFrame:
    blastn = shutil.which("blastn")
    makeblastdb = shutil.which("makeblastdb")
    if blastn is None or makeblastdb is None:
        raise RuntimeError("blastn and makeblastdb must be available on PATH.")

    work_dir.mkdir(parents=True, exist_ok=True)
    ref_fasta = work_dir / "reference_oriv.fasta"
    query_fasta = work_dir / "candidate_oriv.fasta"
    db_prefix = work_dir / "reference_oriv_db"
    blast_tsv = work_dir / "candidate_vs_reference.blastn.tsv"

    ref_records = (
        reference_df[["oriv_sequence"]]
        .dropna()
        .assign(oriv_sequence=lambda x: x["oriv_sequence"].astype(str).str.strip().str.lower())
    )
    ref_records = ref_records[ref_records["oriv_sequence"] != ""].drop_duplicates("oriv_sequence").reset_index(drop=True)
    ref_records["seq_id"] = "ref_" + ref_records.index.astype(str)

    query_records = candidate_df[["oriv_sequence"]].copy()
    query_records["oriv_sequence"] = query_records["oriv_sequence"].astype(str).str.strip().str.lower()
    query_records = query_records[query_records["oriv_sequence"] != ""]
    query_records = query_records.drop_duplicates("oriv_sequence").reset_index(drop=True)
    query_records["query_id"] = "query_" + query_records.index.astype(str)

    if ref_records.empty:
        raise ValueError("Reference CSV contains no usable oriv_sequence values.")
    if query_records.empty:
        return pd.DataFrame(columns=["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qcovs"])

    write_fasta(ref_records, "oriv_sequence", "seq_id", ref_fasta)
    write_fasta(query_records, "oriv_sequence", "query_id", query_fasta)

    run_command([makeblastdb, "-in", str(ref_fasta), "-dbtype", "nucl", "-out", str(db_prefix)])
    outfmt = "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovs"
    run_command(
        [
            blastn,
            "-query",
            str(query_fasta),
            "-db",
            str(db_prefix),
            "-out",
            str(blast_tsv),
            "-outfmt",
            outfmt,
            "-max_target_seqs",
            "10",
        ]
    )

    cols = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qcovs"]
    if blast_tsv.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    hits = pd.read_csv(blast_tsv, sep="\t", names=cols)
    return hits.merge(query_records[["query_id", "oriv_sequence"]], left_on="qseqid", right_on="query_id", how="left")


def write_parquet(df: pd.DataFrame, path: Path) -> None:
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path, compression="snappy")


def main() -> None:
    args = parse_args()
    configure_logging(args.log_level)
    if pd is None or pa is None or pq is None:
        raise RuntimeError(
            "This script requires pandas and pyarrow. Install them with: "
            "python -m pip install pandas pyarrow"
        )

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "blast_work"

    logger.info("Loading reference CSV: %s", args.reference_csv)
    reference_df = pd.read_csv(args.reference_csv, low_memory=False)
    if "oriv_sequence" not in reference_df.columns:
        raise ValueError("Reference CSV must contain an 'oriv_sequence' column.")

    reference_df["oriv_sequence"] = reference_df["oriv_sequence"].apply(lambda x: clean_text(x, lower=True))
    if "rep_protein" in reference_df.columns:
        reference_df["rep_protein"] = reference_df["rep_protein"].apply(lambda x: clean_text(x, upper=True))
    if "plasmid_id" in reference_df.columns:
        reference_df["plasmid_id"] = reference_df["plasmid_id"].apply(clean_text)

    external_df = read_external_tables(args.external_parquet)
    raw_external_rows = len(external_df)

    external_df = external_df[external_df["oriv_sequence"] != ""].copy()
    removed_empty_oriv = raw_external_rows - len(external_df)
    logger.info("Rows after removing empty oriV sequences: %d", len(external_df))

    external_df, removed_internal_duplicates = apply_internal_dedup(external_df, args.internal_dedup_key)
    logger.info(
        "Rows after external internal dedup (%s): %d",
        args.internal_dedup_key,
        len(external_df),
    )

    exact_mask = build_reference_exact_mask(external_df, reference_df, args.reference_exact_key)
    exact_rejected = external_df[exact_mask].copy()
    candidates = external_df[~exact_mask].copy()
    logger.info("Rows after exact dedup against reference: %d", len(candidates))

    candidates = candidates.reset_index(drop=True)
    candidates["candidate_id"] = "cand_" + candidates.index.astype(str)
    blast_hits = run_blastn(reference_df, candidates, work_dir)
    blast_hits.to_csv(output_dir / "blast_hits.tsv", sep="\t", index=False)

    if blast_hits.empty:
        blast_rejected_orivs: set[str] = set()
    else:
        high_similarity_hits = blast_hits[
            (blast_hits["pident"] > args.pident_threshold) & (blast_hits["qcovs"] > args.qcovs_threshold)
        ].copy()
        high_similarity_hits.to_csv(output_dir / "blast_rejected_hits.tsv", sep="\t", index=False)
        blast_rejected_orivs = set(high_similarity_hits["oriv_sequence"].dropna().astype(str))

    blast_rejected = candidates[candidates["oriv_sequence"].isin(blast_rejected_orivs)].copy()
    independent_df = candidates[~candidates["oriv_sequence"].isin(blast_rejected_orivs)].copy()

    audit_cols = STANDARD_COLS + ["source_file", "candidate_id"]
    exact_rejected.to_csv(output_dir / "rejected_exact.csv", index=False, encoding="utf-8-sig")
    blast_rejected[audit_cols].to_csv(output_dir / "rejected_blastn.csv", index=False, encoding="utf-8-sig")

    independent_df = independent_df[STANDARD_COLS + ["source_file"]].copy()
    independent_csv = output_dir / "independent_test.csv"
    independent_parquet = output_dir / "independent_test.parquet"
    independent_df.to_csv(independent_csv, index=False, encoding="utf-8-sig")
    write_parquet(independent_df, independent_parquet)

    report = output_dir / "dedup_report.txt"
    report.write_text(
        "\n".join(
            [
                "Independent test set BLASTn dedup report",
                f"reference_csv: {args.reference_csv}",
                f"external_parquet: {', '.join(args.external_parquet)}",
                f"raw_external_rows: {raw_external_rows}",
                f"removed_empty_oriv: {removed_empty_oriv}",
                f"internal_dedup_key: {args.internal_dedup_key}",
                f"removed_external_internal_duplicates: {removed_internal_duplicates}",
                f"reference_exact_key: {args.reference_exact_key}",
                f"removed_exact_against_reference: {len(exact_rejected)}",
                f"blastn_filter: pident > {args.pident_threshold} and qcovs > {args.qcovs_threshold}",
                f"removed_blastn_against_reference: {len(blast_rejected)}",
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
