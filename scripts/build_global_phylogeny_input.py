#!/usr/bin/env python3
from pathlib import Path
import argparse
import textwrap


def read_fasta(path):
    records = []
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header.replace(" ", "_"), "".join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            records.append((header.replace(" ", "_"), "".join(seq)))
    return records


def main():
    parser = argparse.ArgumentParser(description="Build one global FASTA for phylogeny: DB + queries + outgroups.")
    parser.add_argument("--db", default="data/rbp_db_annotated.fasta")
    parser.add_argument("--query", default="data/query_clean.fasta")
    parser.add_argument("--outgroup", default="data/outgroup.fasta")
    parser.add_argument("--output", default="results/phylogeny/global_phylogeny_input.fasta")
    args = parser.parse_args()

    db_records = read_fasta(Path(args.db))
    query_records = read_fasta(Path(args.query))
    outgroup_records = read_fasta(Path(args.outgroup))
    all_records = db_records + query_records + outgroup_records

    ids = [h.split()[0] for h, _ in all_records]
    duplicates = sorted({x for x in ids if ids.count(x) > 1})
    if duplicates:
        raise ValueError(f"Duplicated FASTA IDs in global phylogeny input: {duplicates}")

    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)
    with open(outp, "w") as out:
        for h, s in all_records:
            out.write(f">{h}\n")
            out.write("\n".join(textwrap.wrap(s, 80)) + "\n")

    print(f"Global phylogeny input written to: {outp}")
    print(f"Reference RBPs: {len(db_records)}")
    print(f"Queries: {len(query_records)}")
    print(f"Outgroups: {len(outgroup_records)}")
    print(f"Total: {len(all_records)}")

if __name__ == "__main__":
    main()
