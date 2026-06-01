#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import textwrap
import pandas as pd

MORPH_MAP = {
    "podoviridae": "Podo", "podovirus": "Podo", "podo": "Podo",
    "myoviridae": "Myo", "myovirus": "Myo", "myo": "Myo",
    "siphoviridae": "Sipho", "siphiviridae": "Sipho", "siphovirus": "Sipho",
    "sipho": "Sipho", "sypho": "Sipho",
}

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ*-")


def safe_field(value):
    if value is None or pd.isna(value):
        return "NA"
    value = str(value).strip()
    if not value:
        return "NA"
    value = value.replace("|", "_")
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    return value or "NA"


def morph_short(value):
    if value is None or pd.isna(value):
        return "NA"
    key = str(value).strip().lower()
    return MORPH_MAP.get(key, safe_field(value))


def clean_sequence(value):
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"\s+", "", str(value)).upper()


def read_fasta(path):
    records = []
    header = None
    seq = []
    if not path or not Path(path).exists():
        return records
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq)))
                header = line[1:].replace(" ", "_")
                seq = []
            else:
                seq.append(line)
        if header is not None:
            records.append((header, "".join(seq)))
    return records


def write_fasta(records, output):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            out.write("\n".join(textwrap.wrap(seq, 80)) + "\n")


def parse_outgroup_metadata(clean_header):
    """
    Parse fixed outgroup metadata.

    Expected header format:
    OUTGROUP_01|Escherichia_phage_T4|tail_fiber_assembly|NP_049864.1|OUTGROUP

    Metadata convention:
    - host: bacterial host genus
    - phage: bacteriophage name
    - morphotype: Myo/Sipho/Podo
    - role: non-RBP tail protein role
    """
    parts = clean_header.split("|")
    outgroup_id = parts[0] if len(parts) > 0 else clean_header
    accession = parts[3] if len(parts) > 3 else "NA"

    fixed_metadata = {
        "OUTGROUP_01": {
            "host": "Escherichia",
            "phage": "Bacteriophage_T4",
            "morphotype": "Myo",
            "role": "tail_fiber_assembly",
        },
        "OUTGROUP_02": {
            "host": "Escherichia",
            "phage": "Bacteriophage_Lambda",
            "morphotype": "Sipho",
            "role": "tail_terminator",
        },
        "OUTGROUP_03": {
            "host": "Escherichia",
            "phage": "Bacteriophage_T5",
            "morphotype": "Sipho",
            "role": "tail_length_tape_measure_protein",
        },
        "OUTGROUP_04": {
            "host": "Mycobacterium",
            "phage": "Bacteriophage_Che8",
            "morphotype": "Sipho",
            "role": "tail_terminator",
        },
        "OUTGROUP_05": {
            "host": "Pseudomonas",
            "phage": "Bacteriophage_LUZ7",
            "morphotype": "Podo",
            "role": "tail_length_tape_measure_protein",
        },
    }

    meta = fixed_metadata.get(outgroup_id)

    if meta:
        return (
            accession,
            meta["host"],
            meta["phage"],
            meta["morphotype"],
            meta["role"],
        )

    # Fallback for unexpected outgroups.
    raw_phage = parts[1] if len(parts) > 1 else "NA"
    raw_role = parts[2] if len(parts) > 2 else "OUTGROUP"

    return accession, raw_phage, raw_phage, "OUTGROUP", raw_role

def main():
    parser = argparse.ArgumentParser(
        description="Build annotated RBP FASTA and combined BLAST FASTA with outgroups."
    )
    parser.add_argument("--excel", default="data/RBP_validated_final.xlsx")
    parser.add_argument("--sheet", default="Sheet1")
    parser.add_argument("--output", default="data/rbp_db_annotated.fasta")
    parser.add_argument("--blast-output", default="data/blast_db_with_outgroups.fasta")
    parser.add_argument("--outgroup", default="data/outgroup.fasta")
    parser.add_argument("--metadata", default="results/database/rbp_db_metadata.tsv")
    args = parser.parse_args()

    df = pd.read_excel(args.excel, sheet_name=args.sheet)
    required = ["RBP Accession", "Host", "Bacteriophage", "Morphotype", "CleanSeq"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in Excel: {missing}")

    seen = set()
    rbp_records = []
    outgroup_records = []
    metadata_rows = []

    for idx, row in df.iterrows():
        acc = safe_field(row["RBP Accession"])
        host = safe_field(row["Host"])
        phage = safe_field(row["Bacteriophage"])
        morph = morph_short(row["Morphotype"])
        seq = clean_sequence(row["CleanSeq"])
        if acc == "NA" or not seq:
            continue
        bad_chars = sorted(set(seq) - AA_ALLOWED)
        if bad_chars:
            raise ValueError(f"Invalid amino-acid characters in {acc}: {bad_chars}")
        seq_id = f"{acc}|{host}|{phage}|{morph}|RBP"
        if seq_id in seen:
            seq_id = f"{seq_id}_{idx}"
        seen.add(seq_id)
        rbp_records.append((seq_id, seq))
        metadata_rows.append((seq_id, acc, host, phage, morph, "RBP", len(seq)))

    for header, seq in read_fasta(args.outgroup):
        clean_header = header.replace(" ", "_")
        if "OUTGROUP" not in clean_header.upper():
            clean_header = f"{clean_header}|OUTGROUP"
        seq = seq.upper()
        bad_chars = sorted(set(seq) - AA_ALLOWED)
        if bad_chars:
            raise ValueError(f"Invalid amino-acid characters in outgroup {clean_header}: {bad_chars}")
        accession, host, phage_or_role, morphotype, role = parse_outgroup_metadata(clean_header)
        outgroup_records.append((clean_header, seq))
        metadata_rows.append((clean_header, accession, host, phage_or_role, morphotype, role, len(seq)))

    write_fasta(rbp_records, args.output)
    write_fasta(rbp_records + outgroup_records, args.blast_output)

    metadata = Path(args.metadata)
    metadata.parent.mkdir(parents=True, exist_ok=True)
    with open(metadata, "w") as out:
        out.write("seq_id\taccession\thost\tphage\tmorphotype\trole\tlength\n")
        for row in metadata_rows:
            out.write("\t".join(map(str, row)) + "\n")

    print(f"RBP-only FASTA written to: {args.output}")
    print(f"Combined BLAST FASTA written to: {args.blast_output}")
    print(f"Metadata written to: {args.metadata}")
    print(f"RBP records: {len(rbp_records)}")
    print(f"Outgroup records added to BLAST DB: {len(outgroup_records)}")


if __name__ == "__main__":
    main()
