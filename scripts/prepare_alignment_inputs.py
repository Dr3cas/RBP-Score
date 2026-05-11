from pathlib import Path
import sys
from collections import defaultdict

blast_file = Path(sys.argv[1])
query_fasta = Path(sys.argv[2])
db_fasta = Path(sys.argv[3])
output_dir = Path(sys.argv[4])

TOP_N = 10

#filtros
EVALUE_CUTOFF = 1e-5 #reduz hits ao acaso, reduzir para 1e-5 para default
PIDENT_CUTOFF = 30 #garante identidade minima, diminuir se der poucos hits 
COVERAGE_CUTOFF = 0.5 #evita alinhamentos curtos sem valor bio, same...pode se reduzir para o.5 ou 0.4

output_dir.mkdir(parents=True, exist_ok=True)


def read_fasta(path):
    records = {}
    current_header = None
    current_seq = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records[current_header] = "".join(current_seq)

                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            records[current_header] = "".join(current_seq)

    return records


query_records = read_fasta(query_fasta)
db_records = read_fasta(db_fasta)

hits_by_query = defaultdict(list)

with open(blast_file, "r") as f:
    for line in f:
        if not line.strip():
            continue

        cols = line.strip().split("\t")
        
        # ignorar header
        if cols[0] == "qseqid":
            continue

        qseqid = cols[0]
        sseqid = cols[1]
        pident = float(cols[2])
        length = float(cols[3])
        evalue = float(cols[10])
        bitscore = float(cols[11])
        qlen = float(cols[12])

        coverage = length / qlen if qlen > 0 else 0

        if (
            evalue < EVALUE_CUTOFF
            and pident >= PIDENT_CUTOFF
            and coverage >= COVERAGE_CUTOFF
        ):
            hits_by_query[qseqid].append({
                "sseqid": sseqid,
                "pident": pident,
                "coverage": coverage,
                "evalue": evalue,
                "bitscore": bitscore,
            })


summary_lines = []

for qseqid, qseq in query_records.items():
    hits = hits_by_query.get(qseqid, [])

    hits = sorted(hits, key=lambda x: x["bitscore"], reverse=True)

    selected_subjects = []
    seen = set()

    for hit in hits:
        sid = hit["sseqid"]
        if sid not in seen and sid in db_records:
            selected_subjects.append(sid)
            seen.add(sid)

        if len(selected_subjects) >= TOP_N:
            break

    safe_qseqid = qseqid.replace("|", "_").replace("/", "_").replace("\\", "_")
    out_fasta = output_dir / f"{safe_qseqid}.fasta"

    with open(out_fasta, "w") as out:
        out.write(f">{qseqid}\n{qseq}\n")

        for sid in selected_subjects:
            out.write(f">{sid}\n{db_records[sid]}\n")

    summary_lines.append(
        f"{qseqid}\t{len(selected_subjects)}\t{out_fasta}"
    )

summary_file = output_dir / "alignment_inputs_summary.tsv"

with open(summary_file, "w") as out:
    out.write("query_id\tn_hits_used\talignment_input_file\n")
    out.write("\n".join(summary_lines) + "\n")

print("FASTA files for alignment generated.")
print(f"Queries processed: {len(query_records)}")
print(f"Summary: {summary_file}")