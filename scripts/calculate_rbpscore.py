from pathlib import Path
import sys
import math
import csv

try:
    from Bio import Phylo
except ImportError:
    raise ImportError("Instala primeiro: pip install biopython")


blast_file = Path(sys.argv[1])
query_fasta = Path(sys.argv[2])
tree_dir = Path(sys.argv[3])
output_file = Path(sys.argv[4])

output_file.parent.mkdir(parents=True, exist_ok=True)


def safe_id(qid):
    return qid.replace("|", "_").replace("/", "_").replace("\\", "_")


def read_query_ids(fasta):
    ids = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def evalue_to_score(evalue):
    if evalue == 0:
        return 1.0
    score = -math.log10(evalue)
    return min(score / 50, 1.0)


def read_tree(tree_path):
    if not tree_path.exists():
        return None

    text = tree_path.read_text().strip()

    if not text or text == "NO_TREE":
        return None

    try:
        return Phylo.read(str(tree_path), "newick")
    except Exception:
        return None


def get_tree_distance(tree, query_id, subject_id):
    if tree is None:
        return None

    terminals = {term.name: term for term in tree.get_terminals()}

    if query_id not in terminals or subject_id not in terminals:
        return None

    try:
        return tree.distance(terminals[query_id], terminals[subject_id])
    except Exception:
        return None


query_ids = read_query_ids(query_fasta)

hits = []

with open(blast_file) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        qseqid = row["qseqid"]
        sseqid = row["sseqid"]

        pident = float(row["pident"])
        length = float(row["length"])
        qlen = float(row["qlen"])
        evalue = float(row["evalue"])
        bitscore = float(row["bitscore"])

        coverage = length / qlen if qlen > 0 else 0

        hits.append({
            "query_id": qseqid,
            "subject_id": sseqid,
            "pident": pident,
            "coverage": coverage,
            "evalue": evalue,
            "bitscore": bitscore,
        })


# máximo bitscore por query para normalização
max_bitscore_by_query = {}

for h in hits:
    q = h["query_id"]
    max_bitscore_by_query[q] = max(
        max_bitscore_by_query.get(q, 0),
        h["bitscore"]
    )


scored_rows = []

for h in hits:
    q = h["query_id"]
    s = h["subject_id"]

    identity_score = h["pident"] / 100
    coverage_score = h["coverage"]

    max_bitscore = max_bitscore_by_query.get(q, 1)
    bitscore_score = h["bitscore"] / max_bitscore if max_bitscore > 0 else 0

    evalue_score = evalue_to_score(h["evalue"])

#mudar valores aqui se necessário

    blast_score = (
        0.35 * identity_score +
        0.30 * coverage_score +
        0.25 * bitscore_score +
        0.10 * evalue_score
    )

    tree_path = tree_dir / f"{safe_id(q)}.tree.nwk"
    tree = read_tree(tree_path)

    tree_distance = get_tree_distance(tree, q, s)

    if tree_distance is None:
        tree_score = None
        final_score = blast_score
        score_basis = "blast_only"
    else:
        tree_score = 1 / (1 + tree_distance)
        final_score = 0.70 * blast_score + 0.30 * tree_score  #mudar valores aqui tbm
        score_basis = "blast_plus_tree"

    scored_rows.append({
        "query_id": q,
        "subject_id": s,
        "pident": h["pident"],
        "coverage": h["coverage"],
        "evalue": h["evalue"],
        "bitscore": h["bitscore"],
        "blast_score": blast_score,
        "tree_distance": tree_distance if tree_distance is not None else "NA",
        "tree_score": tree_score if tree_score is not None else "NA",
        "final_rbp_score": final_score,
        "score_basis": score_basis
    })


# adicionar queries sem hits
queries_with_hits = {h["query_id"] for h in hits}

for q in query_ids:
    if q not in queries_with_hits:
        scored_rows.append({
            "query_id": q,
            "subject_id": "NA",
            "pident": "NA",
            "coverage": "NA",
            "evalue": "NA",
            "bitscore": "NA",
            "blast_score": 0,
            "tree_distance": "NA",
            "tree_score": "NA",
            "final_rbp_score": 0,
            "score_basis": "no_valid_blast_hits"
        })


scored_rows = sorted(
    scored_rows,
    key=lambda x: (x["query_id"], -float(x["final_rbp_score"]))
)


fieldnames = [
    "query_id",
    "subject_id",
    "pident",
    "coverage",
    "evalue",
    "bitscore",
    "blast_score",
    "tree_distance",
    "tree_score",
    "final_rbp_score",
    "score_basis"
]

with open(output_file, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(scored_rows)

print(f"RBP-Score table written to: {output_file}")