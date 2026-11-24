#!/usr/bin/env python3
# 01_extract_metadata.py
# Extract metadata for 315 SINV sequences.
# Input:  data/SINV/315.gb
# Output: results/tables/metadata.csv
# Dependencies: biopython

from pathlib import Path
from Bio import SeqIO
import csv

ROOT     = Path(__file__).resolve().parent.parent
GB_FILE  = ROOT / "data" / "SINV" / "315.gb"
OUT_FILE = ROOT / "results" / "tables" / "metadata.csv"

OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

rows = []
for record in SeqIO.parse(str(GB_FILE), "genbank"):
    host, country, isolate, strain, date = "", "", "", "", ""
    
    for feature in record.features:
        if feature.type == "source":
            q = feature.qualifiers
            host    = q.get("host",    [""])[0]
            country = q.get("country", q.get("geo_loc_name", [""]))[0]
            isolate = q.get("isolate", [""])[0]
            strain  = q.get("strain",  [""])[0]
            date    = q.get("collection_date", [""])[0]
            break

    acc         = record.id.split(".")[0]
    version_str = record.id

    rows.append({
        "Accession": acc,
        "Version": version_str,
        "Definition": record.description,
        "Organism": record.annotations.get("organism", ""),
        "Host": host,
        "Country": country,
        "Isolate": isolate if isolate else strain,
        "Collection_Date": date
    })

with open(OUT_FILE, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print(f"Extracted {len(rows)} sequences.")
print(f"Saved to {OUT_FILE}")
