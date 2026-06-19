#!/usr/bin/env python3
from pathlib import Path
import argparse
import re
import textwrap

AA_ALLOWED = set("ACDEFGHIKLMNPQRSTVWYXBZUOJ*-")

def safe_header(header: str) -> str:
    header = header.strip().replace(">", "")
    header = header.replace(" ", "_")
    header = header.replace("unknownn", "unknown")
    header = header.replace("unknow_", "unknown_")
    header = header.replace("|", "|")
    fields = []
    for field in header.split("|"):
        field = re.sub(r"[^A-Za-z0-9_.=-]+", "_", field)
        field = re.sub(r"_+", "_", field).strip("_")
        fields.append(field or "NA")
    return "|".join(fields)


def read_fasta(path: Path):
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, "".join(seq)


def main():
    parser = argparse.ArgumentParser(description="Clean query FASTA headers for reproducible pipeline use.")
    parser.add_argument("input")
    parser.add_argument("output")
    args = parser.parse_args()

    inp = Path(args.input)
    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)

    seen = set()
    records = []
    for header, seq in read_fasta(inp):
        new_header = safe_header(header)
        if new_header in seen:
            raise ValueError(f"Duplicated query header after cleaning: {new_header}")
        seen.add(new_header)
        seq = re.sub(r"\s+", "", seq).upper()
        bad = sorted(set(seq) - AA_ALLOWED)
        if bad:
            raise ValueError(f"Invalid amino-acid characters in query {new_header}: {bad}")
        records.append((new_header, seq))

    with open(outp, "w") as out:
        for h, s in records:
            out.write(f">{h}\n")
            out.write("\n".join(textwrap.wrap(s, 80)) + "\n")

    print(f"Clean query FASTA written to: {outp}")
    print(f"Queries: {len(records)}")

if __name__ == "__main__":
    main()
