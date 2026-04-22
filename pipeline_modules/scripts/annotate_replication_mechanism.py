from __future__ import annotations

"""Annotate Rep proteins with HMM hits and infer replication mechanism labels."""

import argparse
import logging
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, cast

import pandas as pd
import pyhmmer
from tqdm import tqdm


DEFAULT_INPUT = "merged_final_optimized.csv"
DEFAULT_OUTPUT = "merged_final_optimized_replication_annotated.csv"
DEFAULT_DB = "pipeline_modules/resources/replication_annotation/rep_plasmid_core_db.hmm"
DEFAULT_RULES = "pipeline_modules/resources/replication_annotation/pfam_replication_mechanism_table.tsv"

logger = logging.getLogger(__name__)


PAIR_RULES = [
    ({"Rep3_N", "Rep3_C"}, "theta", "paired_theta_marker:Rep3_N+Rep3_C", "high"),
    (
        {"RepA_N", "Bac_RepA_C"},
        "theta",
        "paired_theta_marker:RepA_N+Bac_RepA_C",
        "high",
    ),
    (
        {"WHD_RepA_N", "Bac_RepA_C"},
        "theta",
        "paired_theta_marker:WHD_RepA_N+Bac_RepA_C",
        "high",
    ),
    ({"Rep_OBD", "RepB_C"}, "theta", "paired_theta_marker:Rep_OBD+RepB_C", "high"),
    (
        {"Replitron_HUH", "Replitron_C"},
        "RCR_like",
        "paired_rcr_like_marker:Replitron_HUH+Replitron_C",
        "high",
    ),
    (
        {"RepD-like_N", "Rep_trans"},
        "RCR",
        "paired_rcr_marker:RepD-like_N+Rep_trans",
        "high",
    ),
]

HIGH_CONFIDENCE_SINGLE = {
    "RCR": {"Rep_1", "RepL", "Rol_Rep_N", "RepD-like_N"},
    "theta": {"Rep3_N", "Rep3_C", "RepA_N", "Bac_RepA_C", "WHD_RepA_N", "RepA_C", "Rep_OBD", "RepC", "TrfA"},
}

LIKE_HINTS = {
    "RCR_like": {"Replitron_HUH", "Replitron_C", "RepB-RCR_reg"},
    "theta_like": {"Rep_assoc_RepA4", "L_lactis_RepB_C", "RepB_primase", "RepB_primase_C", "RepA1_leader"},
}


@dataclass
class HitRecord:
    pfam_name: str
    pfam_acc: str
    score: float
    evalue: float
    included: bool
    reported: bool
    domains: int


def load_rules(path: Path) -> pd.DataFrame:
    rules = pd.read_csv(path, sep="\t")
    rules["pfam_name"] = rules["pfam_name"].astype(str)
    return rules


def build_rule_lookup(rules: pd.DataFrame) -> Dict[str, dict]:
    lookup: Dict[str, dict] = {}
    for row in rules.to_dict(orient="records"):
        lookup[row["pfam_name"]] = row
    return lookup


def clean_sequence(seq: object) -> str:
    if seq is None or (isinstance(seq, float) and math.isnan(seq)):
        return ""
    text = str(seq).strip().upper()
    return "".join(ch for ch in text if "A" <= ch <= "Z")


def sequence_batches(
    unique_sequences: Sequence[str], batch_size: int
) -> Iterable[Sequence[str]]:
    for start in range(0, len(unique_sequences), batch_size):
        yield unique_sequences[start : start + batch_size]


def load_hmms(hmm_path: Path) -> List[pyhmmer.plan7.HMM]:
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as hmm_file:
        return list(hmm_file)


def scan_sequences(
    sequences: Sequence[str],
    hmms: Sequence[pyhmmer.plan7.HMM],
    cpus: int,
    evalue_threshold: float,
) -> Dict[str, List[HitRecord]]:
    alphabet = pyhmmer.easel.Alphabet.amino()
    results: Dict[str, List[HitRecord]] = {}
    batch_size = 2000
    total_batches = max(1, math.ceil(len(sequences) / batch_size))

    for batch in tqdm(
        sequence_batches(sequences, batch_size),
        total=total_batches,
        desc="HMM batches",
        unit="batch",
    ):
        digital_sequences = [
            pyhmmer.easel.TextSequence(
                name=f"seq_{idx}".encode(),
                sequence=sequence,
            ).digitize(alphabet)
            for idx, sequence in enumerate(batch)
        ]

        query_sequences = cast(
            "Iterable[pyhmmer.easel.DigitalSequence[pyhmmer.easel.Alphabet]]",
            digital_sequences,
        )

        top_hits_iter = pyhmmer.hmmer.hmmscan(
            query_sequences,
            hmms,
            cpus=cpus,
            E=evalue_threshold,
            incE=evalue_threshold,
            domE=evalue_threshold,
            incdomE=evalue_threshold,
        )

        for sequence, top_hits in zip(batch, top_hits_iter):
            hits: List[HitRecord] = []
            for hit in top_hits:
                hits.append(
                    HitRecord(
                        pfam_name=str(hit.name),
                        pfam_acc="" if hit.accession is None else str(hit.accession),
                        score=float(hit.score),
                        evalue=float(hit.evalue),
                        included=bool(hit.included),
                        reported=bool(hit.reported),
                        domains=len(hit.domains),
                    )
                )
            results[sequence] = hits

    return results


def rank_hits(hits: Sequence[HitRecord], rule_lookup: Dict[str, dict]) -> List[HitRecord]:
    def rank_key(hit: HitRecord) -> tuple:
        rule = rule_lookup.get(hit.pfam_name, {})
        use_for_call = 1 if str(rule.get("use_for_call", "")).lower() == "yes" else 0
        confidence = {"high": 3, "medium": 2, "low": 1}.get(
            str(rule.get("confidence", "")).lower(),
            0,
        )
        return (use_for_call, confidence, hit.score, -hit.evalue)

    return sorted(hits, key=rank_key, reverse=True)


def infer_mechanism(hits: Sequence[HitRecord], rule_lookup: Dict[str, dict]) -> dict:
    if not hits:
        return {
            "rep_hmm_top_hit": "",
            "rep_hmm_top_acc": "",
            "rep_hmm_top_score": "",
            "rep_hmm_top_evalue": "",
            "rep_hmm_hits": "",
            "replication_mechanism_call": "unresolved",
            "replication_mechanism_basis": "no_significant_hmm_hit",
            "replication_mechanism_confidence": "none",
        }

    ranked_hits = rank_hits(hits, rule_lookup)
    hit_names = {hit.pfam_name for hit in ranked_hits}
    top_hit = ranked_hits[0]

    all_hits_text = ";".join(
        f"{hit.pfam_name}|{hit.pfam_acc}|score={hit.score:.2f}|evalue={hit.evalue:.2e}"
        for hit in ranked_hits
    )

    for required_hits, mechanism, basis, confidence in PAIR_RULES:
        if required_hits.issubset(hit_names):
            return {
                "rep_hmm_top_hit": top_hit.pfam_name,
                "rep_hmm_top_acc": top_hit.pfam_acc,
                "rep_hmm_top_score": round(top_hit.score, 3),
                "rep_hmm_top_evalue": top_hit.evalue,
                "rep_hmm_hits": all_hits_text,
                "replication_mechanism_call": mechanism,
                "replication_mechanism_basis": basis,
                "replication_mechanism_confidence": confidence,
            }

    direct_calls = set()
    direct_bases = []
    for hit in ranked_hits:
        if hit.pfam_name in HIGH_CONFIDENCE_SINGLE["RCR"]:
            direct_calls.add("RCR")
            direct_bases.append(f"single_core_marker:{hit.pfam_name}")
        elif hit.pfam_name in HIGH_CONFIDENCE_SINGLE["theta"]:
            direct_calls.add("theta")
            direct_bases.append(f"single_core_marker:{hit.pfam_name}")

    if len(direct_calls) == 1:
        mechanism = next(iter(direct_calls))
        return {
            "rep_hmm_top_hit": top_hit.pfam_name,
            "rep_hmm_top_acc": top_hit.pfam_acc,
            "rep_hmm_top_score": round(top_hit.score, 3),
            "rep_hmm_top_evalue": top_hit.evalue,
            "rep_hmm_hits": all_hits_text,
            "replication_mechanism_call": mechanism,
            "replication_mechanism_basis": ",".join(direct_bases),
            "replication_mechanism_confidence": "high",
        }

    if len(direct_calls) > 1:
        return {
            "rep_hmm_top_hit": top_hit.pfam_name,
            "rep_hmm_top_acc": top_hit.pfam_acc,
            "rep_hmm_top_score": round(top_hit.score, 3),
            "rep_hmm_top_evalue": top_hit.evalue,
            "rep_hmm_hits": all_hits_text,
            "replication_mechanism_call": "unresolved",
            "replication_mechanism_basis": "conflicting_core_markers",
            "replication_mechanism_confidence": "medium",
        }

    if hit_names & LIKE_HINTS["RCR_like"]:
        basis = ",".join(
            f"support_or_like_marker:{name}"
            for name in ranked_hits_names(ranked_hits)
            if name in LIKE_HINTS["RCR_like"]
        )
        return {
            "rep_hmm_top_hit": top_hit.pfam_name,
            "rep_hmm_top_acc": top_hit.pfam_acc,
            "rep_hmm_top_score": round(top_hit.score, 3),
            "rep_hmm_top_evalue": top_hit.evalue,
            "rep_hmm_hits": all_hits_text,
            "replication_mechanism_call": "RCR_like",
            "replication_mechanism_basis": basis,
            "replication_mechanism_confidence": "medium",
        }

    if hit_names & LIKE_HINTS["theta_like"]:
        basis = ",".join(
            f"support_or_like_marker:{name}"
            for name in ranked_hits_names(ranked_hits)
            if name in LIKE_HINTS["theta_like"]
        )
        return {
            "rep_hmm_top_hit": top_hit.pfam_name,
            "rep_hmm_top_acc": top_hit.pfam_acc,
            "rep_hmm_top_score": round(top_hit.score, 3),
            "rep_hmm_top_evalue": top_hit.evalue,
            "rep_hmm_hits": all_hits_text,
            "replication_mechanism_call": "theta_like",
            "replication_mechanism_basis": basis,
            "replication_mechanism_confidence": "medium",
        }

    return {
        "rep_hmm_top_hit": top_hit.pfam_name,
        "rep_hmm_top_acc": top_hit.pfam_acc,
        "rep_hmm_top_score": round(top_hit.score, 3),
        "rep_hmm_top_evalue": top_hit.evalue,
        "rep_hmm_hits": all_hits_text,
        "replication_mechanism_call": "unresolved",
        "replication_mechanism_basis": "mechanism_unresolved_from_available_hits",
        "replication_mechanism_confidence": "low",
    }


def ranked_hits_names(hits: Sequence[HitRecord]) -> List[str]:
    return [hit.pfam_name for hit in hits]


def annotate_chunk(
    chunk: pd.DataFrame,
    cache: Dict[str, dict],
    hmms: Sequence[pyhmmer.plan7.HMM],
    rule_lookup: Dict[str, dict],
    cpus: int,
    evalue_threshold: float,
) -> pd.DataFrame:
    cleaned = chunk["rep_seq"].map(clean_sequence)
    chunk = chunk.copy()
    chunk["_clean_rep_seq"] = cleaned

    sequences_to_scan = [
        seq for seq in chunk["_clean_rep_seq"].unique().tolist() if seq and seq not in cache
    ]
    if sequences_to_scan:
        scan_results = scan_sequences(
            sequences=sequences_to_scan,
            hmms=hmms,
            cpus=cpus,
            evalue_threshold=evalue_threshold,
        )
        for sequence, hits in scan_results.items():
            cache[sequence] = infer_mechanism(hits, rule_lookup)

    blank_annotation = infer_mechanism([], rule_lookup)
    annotations = []
    for sequence in chunk["_clean_rep_seq"]:
        annotations.append(cache.get(sequence, blank_annotation) if sequence else blank_annotation)

    annotation_df = pd.DataFrame(annotations, index=chunk.index)
    chunk = pd.concat([chunk.drop(columns=["_clean_rep_seq"]), annotation_df], axis=1)
    return chunk


def run(
    input_csv: Path,
    output_csv: Path,
    hmm_path: Path,
    rules_path: Path,
    chunksize: int,
    cpus: int,
    evalue_threshold: float,
) -> None:
    logger.info("Loading rules from %s", rules_path)
    rules = load_rules(rules_path)
    rule_lookup = build_rule_lookup(rules)
    cache: Dict[str, dict] = {}
    logger.info("Loading HMM database from %s", hmm_path)
    hmms = load_hmms(hmm_path)

    with open(input_csv, "r", encoding="utf-8", errors="ignore") as handle:
        total_rows = sum(1 for _ in handle) - 1
    first_chunk = True

    for chunk in tqdm(
        pd.read_csv(input_csv, chunksize=chunksize),
        total=max(1, math.ceil(total_rows / chunksize)),
        desc="CSV chunks",
        unit="chunk",
    ):
        annotated_chunk = annotate_chunk(
            chunk=chunk,
            cache=cache,
            hmms=hmms,
            rule_lookup=rule_lookup,
            cpus=cpus,
            evalue_threshold=evalue_threshold,
        )
        annotated_chunk.to_csv(
            output_csv,
            mode="w" if first_chunk else "a",
            header=first_chunk,
            index=False,
        )
        first_chunk = False
    logger.info("Annotation completed. Output written to %s", output_csv)


def configure_logging(level: str, log_file: str | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_path, encoding="utf-8"))

    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=handlers,
        force=True,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate protein sequences with rep_plasmid_core_db.hmm and infer replication mechanism."
    )
    parser.add_argument("--input", "--input-csv", dest="input", default=DEFAULT_INPUT, help="Input CSV path.")
    parser.add_argument("--output", "--output-csv", dest="output", default=DEFAULT_OUTPUT, help="Output CSV path.")
    parser.add_argument("--hmm-db", default=DEFAULT_DB, help="HMM database path.")
    parser.add_argument("--rules", default=DEFAULT_RULES, help="Pfam-to-mechanism rules TSV.")
    parser.add_argument("--chunk-size", "--chunksize", dest="chunksize", type=int, default=20000, help="CSV rows to process per chunk.")
    parser.add_argument("--cpus", type=int, default=max(1, min(os.cpu_count() or 1, 8)), help="Threads for pyhmmer.")
    parser.add_argument(
        "--evalue-threshold",
        type=float,
        default=1e-3,
        help="Per-hit e-value threshold used in hmmscan.",
    )
    parser.add_argument("--log-level", default="INFO", help="Logging level: DEBUG, INFO, WARNING, ERROR.")
    parser.add_argument("--log-file", default=None, help="Optional log file path.")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    configure_logging(args.log_level, args.log_file)
    run(
        input_csv=Path(args.input),
        output_csv=Path(args.output),
        hmm_path=Path(args.hmm_db),
        rules_path=Path(args.rules),
        chunksize=args.chunksize,
        cpus=args.cpus,
        evalue_threshold=args.evalue_threshold,
    )
