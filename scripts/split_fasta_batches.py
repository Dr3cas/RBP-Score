from pathlib import Path

input_fasta = Path("data/rbp_db.fasta")
out_dir = Path("batches")
out_dir.mkdir(exist_ok=True)

batch_size = 10

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

for i in range(0, len(records), batch_size):
    batch = records[i:i+batch_size]
    out = out_dir / f"batch_{i//batch_size + 1}.fasta"
    with open(out, "w") as f:
        for h, s in batch:
            f.write(f"{h}\n{s}\n")

print(f"Total sequences: {len(records)}")
print(f"Batches created: {len(list(out_dir.glob('*.fasta')))}")
