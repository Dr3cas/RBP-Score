from pathlib import Path
import re

input_fasta = Path("data/rbp_db.fasta")
output_dir = Path("reference_structures/input_fastas")
output_dir.mkdir(parents=True, exist_ok=True)


def safe_name(header):
    header = header.strip().replace(">", "")
    header = header.split()[0]
    header = header.replace("|", "_")
    header = header.replace("/", "_")
    header = header.replace("\\", "_")
    header = re.sub(r"[^A-Za-z0-9_.-]", "_", header)
    return header


current_header = None
current_seq = []

records = []

with open(input_fasta) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                records.append((current_header, "".join(current_seq)))

            current_header = line
            current_seq = []
        else:
            current_seq.append(line)

    if current_header is not None:
        records.append((current_header, "".join(current_seq)))


for header, seq in records:
    name = safe_name(header)
    out = output_dir / f"{name}.fasta"

    with open(out, "w") as f:
        f.write(f"{header}\n")
        f.write(f"{seq}\n")

print(f"Created {len(records)} individual FASTA files in {output_dir}")