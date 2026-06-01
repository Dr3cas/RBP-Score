from pathlib import Path
import re

input_fasta = Path("data/outgroup.fasta")
output_dir = Path("reference_structures/outgroup_fastas")
output_dir.mkdir(parents=True, exist_ok=True)

def safe_name(header):
    header = header.strip().replace(">", "")
    header = header.split()[0]
    header = header.replace("|", "_")
    header = re.sub(r"[^A-Za-z0-9_.-]", "_", header)
    return header

records = []
header = None
seq = []

with open(input_fasta) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq)))
            header = line
            seq = []
        else:
            seq.append(line)

if header:
    records.append((header, "".join(seq)))

for h, s in records:
    name = safe_name(h)
    with open(output_dir / f"{name}.fasta", "w") as out:
        out.write(f"{h}\n{s}\n")

print(f"Outgroup FASTA files written: {len(records)}")
