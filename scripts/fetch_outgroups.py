"""
fetch_outgroups.py
------------------
Fetches non-RBP tail chaperone / tail-associated proteins from a curated
list of EXTERNAL phages (NOT from your dataset) to use as phylogenetic
outgroups.

Outgroup phages were chosen to be:
  - Well-characterised, with complete annotated genomes
  - Taxonomically distant from each other
  - Free of the host-range bias of your query dataset

Usage:
    pip install biopython
    python fetch_outgroups.py
"""

from Bio import Entrez, SeqIO
from pathlib import Path
import re
import time

# ── CONFIGURE ────────────────────────────────────────────────────────────────
Entrez.email = "pg59751@uminho.pt"   # <-- altera aqui

OUT_FASTA = Path("outgroup.fasta")
OUT_TSV   = Path("outgroup_candidates.tsv")
N_OUTGROUPS = 5

# ── EXTERNAL REFERENCE PHAGES ────────────────────────────────────────────────
# One outgroup candidate per phage; all well-annotated complete genomes.
# These are deliberately OUTSIDE typical Caudoviricetes RBP study scope.
EXTERNAL_PHAGES = [
    "NC_000866",   # Enterobacteria phage T4        (Myovirus)
    "NC_001416",   # Enterobacteria phage lambda    (Siphovirus)
    "NC_005859",   # Bacteriophage T5               (Siphovirus, distant from T4)
    "NC_004680",   # Bacillus phage SPP1            (Siphovirus, Gram-positive host)
    "NC_013691",   # Pseudomonas phage phiKZ        (Giant phage, very distant)
    "NC_004813",   # Staphylococcus phage 11        (Siphovirus, fallback)
    "NC_007041",   # Mycobacterium phage D29        (Siphovirus, fallback)
]

# ── ANNOTATION FILTERS ───────────────────────────────────────────────────────
INCLUDE_TERMS = [
    "tail assembly chaperone",
    "tail chaperone",
    "tape measure protein",
    "tail completion protein",
    "tail terminator",
    "major tail protein",
    "tail tube protein",
    # T4-specific annotations that are non-RBP tail-associated
    "tail fiber assembly",   # e.g. gp57A in T4 — chaperone, not RBP itself
]

EXCLUDE_TERMS = [
    "receptor",
    "receptor-binding",
    " rbp",
    "tail spike",
    "baseplate",
    "depolymerase",
    "host specificity",
    "adhesin",
    "lysin",
    "holin",
    "capsid",
    "portal",
    "terminase",
    "integrase",
    # exclude the actual fiber structural proteins (RBPs) but keep chaperones
    "long tail fiber",
    "short tail fiber",
    "tail fiber protein",
    "tail fibre protein",
]

# ── HELPERS ──────────────────────────────────────────────────────────────────
def clean(text: str) -> str:
    """Replace spaces and special chars for FASTA header use."""
    text = str(text).strip().replace(" ", "_").replace("|", "_")
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", text)

def is_candidate(product: str) -> bool:
    p = product.lower()
    if not any(t in p for t in INCLUDE_TERMS):
        return False
    if any(t in p for t in EXCLUDE_TERMS):
        return False
    return True

# ── MAIN ─────────────────────────────────────────────────────────────────────
candidates = []

for genome_acc in EXTERNAL_PHAGES:
    print(f"[+] Fetching genome {genome_acc} ...")
    try:
        handle = Entrez.efetch(
            db="nuccore", id=genome_acc, rettype="gb", retmode="text"
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()
        time.sleep(0.4)   # NCBI rate limit courtesy
    except Exception as e:
        print(f"    ERROR fetching {genome_acc}: {e}")
        continue

    # Use organism name from record
    organism = record.annotations.get("organism", genome_acc)

    for feature in record.features:
        if feature.type != "CDS":
            continue

        product    = feature.qualifiers.get("product",    [""])[0]
        protein_id = feature.qualifiers.get("protein_id", ["NA"])[0]
        translation= feature.qualifiers.get("translation",[""])[0]

        if not translation or len(translation) < 80:
            continue
        if not is_candidate(product):
            continue

        candidates.append({
            "genome_accession": genome_acc,
            "organism":         organism,
            "protein_id":       protein_id,
            "product":          product,
            "length":           len(translation),
            "sequence":         translation,
        })
        # Take only the FIRST match per phage to ensure diversity
        break

    if not candidates or candidates[-1]["genome_accession"] != genome_acc:
        print(f"    WARNING: no suitable candidate found in {genome_acc}")

# ── WRITE TSV (all candidates for inspection) ────────────────────────────────
with open(OUT_TSV, "w") as fh:
    fh.write("genome_accession\torganism\tprotein_id\tproduct\tlength\tsequence\n")
    for c in candidates:
        fh.write(
            f"{c['genome_accession']}\t{c['organism']}\t{c['protein_id']}\t"
            f"{c['product']}\t{c['length']}\t{c['sequence']}\n"
        )

# ── WRITE FASTA (template format) ────────────────────────────────────────────
selected = candidates[:N_OUTGROUPS]

with open(OUT_FASTA, "w") as fh:
    for i, c in enumerate(selected, start=1):
        # Header format: >OUTGROUP_XX|Organism|Protein_name|Accession|OUTGROUP
        header = (
            f">OUTGROUP_{i:02d}"
            f"|{clean(c['organism'])}"
            f"|{clean(c['product'])}"
            f"|{c['protein_id']}"
            f"|OUTGROUP"
        )
        # Wrap sequence at 60 chars
        seq = c["sequence"]
        wrapped = "\n".join(seq[j:j+60] for j in range(0, len(seq), 60))
        fh.write(f"{header}\n{wrapped}\n\n")

print(f"\nDone.")
print(f"  Candidates TSV : {OUT_TSV}  ({len(candidates)} entries)")
print(f"  Outgroup FASTA : {OUT_FASTA}  ({len(selected)} sequences)")