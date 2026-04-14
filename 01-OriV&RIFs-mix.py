import os
import pickle
import re
import time
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio import Entrez


print("Start processing data...")

Entrez.email = "nk3901@foxmail.com"


def normalize_accession(value):
    if pd.isna(value):
        return np.nan

    value = str(value).strip()
    if not value or value.lower() == "nan":
        return np.nan

    return value.split(".")[0]


def is_placeholder_accession(value):
    if pd.isna(value):
        return True

    value = str(value).strip()
    return value in {"", "Accession_Number", "accession", "accession_number"}


def normalize_esummary_records(records: Any) -> list[dict[str, Any]]:
    if isinstance(records, Mapping):
        if "DocumentSummarySet" in records:
            summary_set = records["DocumentSummarySet"]
            if isinstance(summary_set, Mapping):
                document_summary = summary_set.get("DocumentSummary", [])
                if isinstance(document_summary, Iterable) and not isinstance(
                    document_summary, (str, bytes, Mapping)
                ):
                    return [dict(item) for item in document_summary if isinstance(item, Mapping)]
        document_summary = records.get("DocumentSummary", [])
        if isinstance(document_summary, Iterable) and not isinstance(
            document_summary, (str, bytes, Mapping)
        ):
            return [dict(item) for item in document_summary if isinstance(item, Mapping)]
        return []

    if isinstance(records, Iterable) and not isinstance(records, (str, bytes, Mapping)):
        return [dict(item) for item in records if isinstance(item, Mapping)]

    return []


def fetch_accession_taxids(accessions, batch_size=200, max_retries=3):
    taxid_map = {}
    acc_list = list(accessions)

    for i in range(0, len(acc_list), batch_size):
        batch = [acc for acc in acc_list[i:i + batch_size] if pd.notna(acc)]
        if not batch:
            continue

        for attempt in range(max_retries):
            try:
                handle = Entrez.esummary(db="nuccore", id=",".join(batch))
                try:
                    records = normalize_esummary_records(Entrez.read(handle))
                finally:
                    handle.close()

                for rec in records:
                    acc_version = rec.get("AccessionVersion", "")
                    acc = acc_version.split(".")[0] if acc_version else None
                    if not acc:
                        continue

                    taxid_map[acc] = rec.get("TaxId")

                break
            except Exception as e:
                print(f"Batch {i} attempt {attempt + 1} failed: {e}")
                time.sleep(2 + attempt * 2)
        else:
            print(f"Batch {i} failed after {max_retries} retries.")

        time.sleep(0.5)

        if i % 2000 == 0:
            print(f"Processed {i}/{len(acc_list)} accessions")

    return taxid_map


def fetch_single_accession_taxid(accession, max_retries=3):
    accession = normalize_accession(accession)
    if is_placeholder_accession(accession):
        return None

    for attempt in range(max_retries):
        try:
            handle = Entrez.esummary(db="nuccore", id=accession)
            try:
                records = normalize_esummary_records(Entrez.read(handle))
            finally:
                handle.close()

            for rec in records:
                acc_version = rec.get("AccessionVersion", "")
                acc = acc_version.split(".")[0] if acc_version else accession
                if acc:
                    return rec.get("TaxId")
            return None
        except Exception as e:
            print(f"Single accession {accession} attempt {attempt + 1} failed: {e}")
            time.sleep(2 + attempt * 2)

    return None


def extract_species_from_taxonomy_record(record):
    if not isinstance(record, Mapping):
        return "unknown"

    rank = str(record.get("Rank", "")).lower()
    scientific_name = str(record.get("ScientificName", "")).strip()

    if rank == "species" and scientific_name:
        return scientific_name

    lineage = record.get("LineageEx", [])
    if isinstance(lineage, Iterable) and not isinstance(lineage, (str, bytes, Mapping)):
        for node in reversed(list(lineage)):
            if not isinstance(node, Mapping):
                continue
            if str(node.get("Rank", "")).lower() == "species":
                candidate = str(node.get("ScientificName", "")).strip()
                if candidate:
                    return candidate

    return scientific_name or "unknown"


def fetch_species_by_taxid(taxids, batch_size=200, max_retries=3):
    species_by_taxid = {}
    taxid_list = [str(taxid) for taxid in taxids if pd.notna(taxid)]

    for i in range(0, len(taxid_list), batch_size):
        batch = taxid_list[i:i + batch_size]
        if not batch:
            continue

        for attempt in range(max_retries):
            try:
                handle = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
                try:
                    records = Entrez.read(handle)
                finally:
                    handle.close()

                for rec in records:
                    taxid = str(rec.get("TaxId", "")).strip()
                    if taxid:
                        species_by_taxid[taxid] = extract_species_from_taxonomy_record(rec)

                break
            except Exception as e:
                print(f"Taxonomy batch {i} attempt {attempt + 1} failed: {e}")
                time.sleep(2 + attempt * 2)
        else:
            print(f"Taxonomy batch {i} failed after {max_retries} retries.")

        time.sleep(0.5)

    return species_by_taxid


def update_taxonomy_cache(cache_file, taxid_map, species_title_map, species_name_by_taxid):
    with open(cache_file, "wb") as f:
        pickle.dump((taxid_map, species_title_map, species_name_by_taxid), f)


def fetch_missing_accessions_and_update_cache(
    missing_accessions,
    cache_file,
    taxid_map,
    species_title_map,
    species_name_by_taxid,
):
    missing_accessions = [
        acc for acc in (normalize_accession(acc) for acc in missing_accessions)
        if not is_placeholder_accession(acc) and acc not in taxid_map
    ]

    if not missing_accessions:
        return taxid_map, species_name_by_taxid

    print(f"Fetching {len(missing_accessions)} missing accessions...")
    fetched_taxids = fetch_accession_taxids(missing_accessions, batch_size=50)
    taxid_map.update(fetched_taxids)

    still_missing = [acc for acc in missing_accessions if acc not in taxid_map]
    if still_missing:
        print(f"Retrying {len(still_missing)} accessions one by one...")
        for acc in still_missing:
            taxid = fetch_single_accession_taxid(acc)
            if taxid is not None:
                taxid_map[acc] = taxid

    new_taxids = sorted(
        {
            str(taxid)
            for acc in missing_accessions
            for taxid in [taxid_map.get(acc)]
            if pd.notna(taxid) and str(taxid) not in species_name_by_taxid
        }
    )
    if new_taxids:
        species_name_by_taxid.update(fetch_species_by_taxid(new_taxids))

    update_taxonomy_cache(cache_file, taxid_map, species_title_map, species_name_by_taxid)

    unresolved = [acc for acc in missing_accessions if acc not in taxid_map]
    if unresolved:
        print(f"Still unresolved after single-accession retry: {len(unresolved)}")
        print("Unresolved examples:", unresolved[:20])

    return taxid_map, species_name_by_taxid


def clean_species(name):
    if pd.isna(name):
        return "unknown"

    name = str(name).strip()
    name = re.sub(r" plasmid.*", "", name, flags=re.IGNORECASE)
    parts = name.split()

    if not parts:
        return "unknown"

    genus = parts[0].rstrip(":;,")

    if len(parts) == 1:
        return genus

    species = parts[1].rstrip(":;,")

    if species.lower() in {"sp", "sp."}:
        return genus

    if not species.islower():
        return genus

    return f"{genus} {species}"


def valid_replicon_length(value):
    if pd.isna(value):
        return np.nan

    try:
        length = float(value)
    except (TypeError, ValueError):
        return np.nan

    return length if length > 0 else np.nan


def feature_midpoint(start, end, replicon_length=np.nan):
    start = float(start)
    end = float(end)

    if pd.notna(replicon_length) and start > end:
        midpoint = (start + end + float(replicon_length)) / 2.0
        if midpoint > float(replicon_length):
            midpoint -= float(replicon_length)
        return midpoint

    return (start + end) / 2.0


def circular_distance(position_a, position_b, replicon_length=np.nan):
    linear_distance = abs(float(position_a) - float(position_b))
    if pd.isna(replicon_length):
        return linear_distance

    length = float(replicon_length)
    if length <= 0:
        return linear_distance

    wrapped_distance = linear_distance % length
    return min(wrapped_distance, length - wrapped_distance)


def pair_reps_to_nearest_oris(df_ori_small, df_rip_small):
    """
    Pair each rep with the nearest OriC within the same plasmid.

    If Sequence_length is available, distances are measured on circular plasmid
    coordinates so features spanning the coordinate origin are handled correctly.
    This avoids Cartesian expansion and keeps one paired row per rep.
    """
    paired_rows = []

    ori_groups = {
        acc: group.reset_index(drop=True)
        for acc, group in df_ori_small.groupby("Accession_Number", sort=False)
    }
    rep_groups = {
        acc: group.reset_index(drop=True)
        for acc, group in df_rip_small.groupby("Accession_Number", sort=False)
    }

    for acc, rep_group in rep_groups.items():
        ori_group = ori_groups.get(acc)
        if ori_group is None or ori_group.empty or rep_group.empty:
            continue

        ori_group = ori_group.dropna(subset=["ori_start", "ori_end"]).reset_index(drop=True)
        rep_group = rep_group.dropna(subset=["rep_start", "rep_end"]).reset_index(drop=True)
        if ori_group.empty or rep_group.empty:
            continue

        length_values = rep_group["Sequence_length"].dropna().map(valid_replicon_length).dropna()
        replicon_length = length_values.iloc[0] if not length_values.empty else np.nan

        ori_mid = np.array(
            [
                feature_midpoint(row["ori_start"], row["ori_end"], replicon_length)
                for _, row in ori_group.iterrows()
            ],
            dtype=float,
        )

        rep_mid = np.array(
            [
                feature_midpoint(row["rep_start"], row["rep_end"], replicon_length)
                for _, row in rep_group.iterrows()
            ],
            dtype=float,
        )

        for rep_idx, rep_row in rep_group.iterrows():
            rep_mid_value = float(rep_mid[rep_idx])
            distances = np.array(
                [
                    circular_distance(ori_mid_value, rep_mid_value, replicon_length)
                    for ori_mid_value in ori_mid
                ],
                dtype=float,
            )
            best_pos = int(np.argmin(distances))

            ori_row = ori_group.iloc[best_pos]
            distance = float(distances[best_pos])
            pairing_method = (
                "nearest_ori_by_rep_circular"
                if pd.notna(replicon_length)
                else "nearest_ori_by_rep_linear"
            )

            paired_rows.append(
                {
                    "OriC sequence": ori_row["Intergenic_Sequence"],
                    "ori_id": ori_row["ori_id"],
                    "plasmid_id": acc,
                    "species": ori_row["species_clean"],
                    "taxid": ori_row["taxid"],
                    "source": "PLSDB",
                    "pfamid_fast": np.nan,
                    "rep_id": rep_row["rep_id"],
                    "Rep_type_fast": rep_row["product"],
                    "rep_seq": rep_row["translation"],
                    "rep_dna_seq": np.nan,
                    "full_replicon_seq": np.nan,
                    "split": np.nan,
                    "__index_level_0__": np.nan,
                    "ori_start": ori_row["ori_start"],
                    "ori_end": ori_row["ori_end"],
                    "rep_start": rep_row["rep_start"],
                    "rep_end": rep_row["rep_end"],
                    "distance": distance,
                    "pairing_method": pairing_method,
                    "ori_midpoint": float(ori_mid[best_pos]),
                    "rep_midpoint": rep_mid_value,
                    "replicon_length": replicon_length,
                }
            )

    if not paired_rows:
        return pd.DataFrame(
            columns=[
                "OriC sequence",
                "ori_id",
                "plasmid_id",
                "species",
                "taxid",
                "source",
                "pfamid_fast",
                "rep_id",
                "Rep_type_fast",
                "rep_seq",
                "rep_dna_seq",
                "full_replicon_seq",
                "split",
                "__index_level_0__",
                "ori_start",
                "ori_end",
                "rep_start",
                "rep_end",
                "distance",
                "pairing_method",
                "ori_midpoint",
                "rep_midpoint",
                "replicon_length",
            ]
        )

    return pd.DataFrame(paired_rows)


df_rip = pd.read_csv("data/RIPs.csv", low_memory=False)
df_ori = pd.read_csv("data/selected_ori_regions.csv", usecols=range(5), low_memory=False)

print(f"RIP rows: {len(df_rip)}")
print(f"OriC rows: {len(df_ori)}")

print("Fetching taxonomy...")

header_mask = df_ori["Accession_Number"].astype(str).str.strip().eq("Accession_Number")
if header_mask.any():
    print(f"Dropping {header_mask.sum()} embedded header rows from OriC table")
    df_ori = df_ori.loc[~header_mask].copy()

df_ori["Accession_Number"] = df_ori["Accession_Number"].apply(normalize_accession)
df_ori["Intergenic_Start"] = pd.to_numeric(df_ori["Intergenic_Start"], errors="coerce")
df_ori["Intergenic_End"] = pd.to_numeric(df_ori["Intergenic_End"], errors="coerce")
valid_acc_mask = ~df_ori["Accession_Number"].apply(is_placeholder_accession)
if (~valid_acc_mask).any():
    print(f"Dropping {(~valid_acc_mask).sum()} invalid accession rows from OriC table")
    df_ori = df_ori.loc[valid_acc_mask].copy()

unique_acc = df_ori["Accession_Number"].dropna().unique()

print(f"Need to query {len(unique_acc)} accessions")

cache_file = "taxonomy_cache.pkl"

if os.path.exists(cache_file):
    print("Using cached taxonomy results...")
    with open(cache_file, "rb") as f:
        cache_data = pickle.load(f)

    if isinstance(cache_data, tuple) and len(cache_data) == 3:
        taxid_map, species_title_map, species_name_by_taxid = cache_data
    elif isinstance(cache_data, tuple) and len(cache_data) == 2:
        taxid_map, species_title_map = cache_data
        species_name_by_taxid = {}
    else:
        taxid_map = {}
        species_title_map = {}
        species_name_by_taxid = {}
else:
    taxid_map = fetch_accession_taxids(unique_acc)
    unique_taxids = pd.Series(taxid_map.values()).dropna().astype(str).unique()
    species_name_by_taxid = fetch_species_by_taxid(unique_taxids)
    species_title_map = {}
    update_taxonomy_cache(cache_file, taxid_map, species_title_map, species_name_by_taxid)

missing_accessions = [acc for acc in unique_acc if acc not in taxid_map]
if missing_accessions:
    taxid_map, species_name_by_taxid = fetch_missing_accessions_and_update_cache(
        missing_accessions,
        cache_file,
        taxid_map,
        species_title_map,
        species_name_by_taxid,
    )

df_ori["taxid"] = df_ori["Accession_Number"].map(taxid_map)
if not species_name_by_taxid and taxid_map:
    unique_taxids = pd.Series(taxid_map.values()).dropna().astype(str).unique()
    species_name_by_taxid = fetch_species_by_taxid(unique_taxids)
    update_taxonomy_cache(cache_file, taxid_map, species_title_map, species_name_by_taxid)

df_ori["species"] = df_ori["taxid"].astype("Int64").astype(str).map(species_name_by_taxid)
df_ori.loc[df_ori["taxid"].isna(), "species"] = np.nan
df_ori["species_clean"] = df_ori["species"].apply(clean_species)

df_ori["ori_start"] = df_ori["Intergenic_Start"]
df_ori["ori_end"] = df_ori["Intergenic_End"]
df_ori["ori_id"] = (
    df_ori["Accession_Number"].astype(str)
    + "_"
    + df_ori["ori_start"].astype(str)
    + "_"
    + df_ori["ori_end"].astype(str)
)

df_rip["Accession_Number"] = df_rip["Accession_Number"].astype(str).str.split(".").str[0]
df_rip["Sequence_length"] = pd.to_numeric(df_rip["Sequence_length"], errors="coerce")
df_rip["gene_start"] = pd.to_numeric(df_rip["gene_start"], errors="coerce")
df_rip["gene_end"] = pd.to_numeric(df_rip["gene_end"], errors="coerce")
df_rip["rep_start"] = df_rip["gene_start"]
df_rip["rep_end"] = df_rip["gene_end"]
df_rip["rep_id"] = (
    df_rip["Accession_Number"].astype(str)
    + "_"
    + df_rip["rep_start"].astype(str)
    + "_"
    + df_rip["rep_end"].astype(str)
)

df_ori_small = df_ori[
    [
        "Accession_Number",
        "ori_id",
        "ori_start",
        "ori_end",
        "Intergenic_Sequence",
        "species_clean",
        "taxid",
    ]
]

df_rip_small = df_rip[
    [
        "Accession_Number",
        "rep_id",
        "rep_start",
        "rep_end",
        "Sequence_length",
        "gene_id",
        "product",
        "translation",
    ]
]

print("Deduplicating Ori and Rep records by ori_id/rep_id before nearest-position pairing...")
df_ori_small = df_ori_small.drop_duplicates(subset=["ori_id"])
df_rip_small = df_rip_small.drop_duplicates(subset=["rep_id"])
print(f"Unique OriC rows: {len(df_ori_small)}, Unique Rep rows: {len(df_rip_small)}")

print("Pairing each rep with its nearest OriC within the same plasmid...")
df_final = pair_reps_to_nearest_oris(df_ori_small, df_rip_small)

print(f"Paired rows: {len(df_final)}")
output_file = Path("merged_final_optimized.csv")
temp_file = output_file.with_suffix(output_file.suffix + ".tmp")
try:
    df_final.to_csv(temp_file, index=False)
    temp_file.replace(output_file)
except PermissionError:
    print(f"Permission denied writing {output_file}. It may already exist and be open or locked by another process.")
    print("Please close the file if it is open, delete it if needed, and rerun the script.")
    raise
except Exception as exc:
    print(f"Failed to write output file {output_file}: {exc}")
    raise

print("Done.")
print(f"Final row count: {len(df_final)}")
print(f"Output file: {output_file}")
print("\nQuality checks:")
print("Missing taxid ratio:", df_final["taxid"].isna().mean())
print("Unknown species ratio:", (df_final["species"] == "unknown").mean())
