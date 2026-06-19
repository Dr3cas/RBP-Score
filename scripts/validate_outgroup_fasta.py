#!/usr/bin/env python3
from pathlib import Path
import argparse
import re

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ*-")

def main():
    parser = argparse.ArgumentParser(description="Validate outgroup FASTA for global phylogeny.")
    parser.add_argument("input", default="data/outgroup.fasta")
    parser.add_argument("output", default="results/phylogeny/outgroup.validated.txt")
    parser.add_argument("--min-n", type=int, default=3)
    parser.add_argument("--max-n", type=int, default=7)
    args = parser.parse_args()

    inp = Path(args.input)
    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)

    if not inp.exists():
        raise FileNotFoundError(
            "Missing data/outgroup.fasta. Add 4-5 non-RBP tail proteins, "
            "preferably tail chaperone/tail-associated proteins, with ROLE=OUTGROUP in the header."
        )

    ids = []
    header = None
    seq = []
    records = []
    with open(inp) as f:
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

    if not (args.min_n <= len(records) <= args.max_n):
        raise ValueError(f"Expected {args.min_n}-{args.max_n} outgroup sequences, found {len(records)}")

    for h, s in records:
        if "|OUTGROUP" not in h and "ROLE=OUTGROUP" not in h:
            raise ValueError(f"Outgroup header should include OUTGROUP/ROLE=OUTGROUP: {h}")
        seq_id = h.split()[0]
        if seq_id in ids:
            raise ValueError(f"Duplicated outgroup ID: {seq_id}")
        ids.append(seq_id)
        clean_seq = re.sub(r"\s+", "", s).upper()
        bad = sorted(set(clean_seq) - AA_ALLOWED)
        if bad:
            raise ValueError(f"Invalid amino-acid characters in outgroup {seq_id}: {bad}")

    outp.write_text(f"validated_outgroups\t{len(records)}\n")
    print(f"Outgroup FASTA validated: {len(records)} sequences")

if __name__ == "__main__":
    main()
