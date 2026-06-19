from pathlib import Path
import shutil
import re

raw_dir = Path("reference_structures/colabfold_raw")
out_dir = Path("reference_structures/pdb")
out_dir.mkdir(parents=True, exist_ok=True)

pdb_files = list(raw_dir.rglob("*.pdb"))

selected = []

for pdb in pdb_files:
    name = pdb.name.lower()

    # tentar apanhar o melhor modelo/ranking
    if "rank_001" in name or "rank_1" in name or "rank_001_" in name:
        selected.append(pdb)

# fallback: se não encontrar rank, copiar todos os pdbs únicos
if not selected:
    selected = pdb_files

copied = 0
seen_prefixes = set()

for pdb in selected:
    # simplificar nome
    stem = pdb.stem
    prefix = re.sub(r"_unrelaxed.*", "", stem)
    prefix = re.sub(r"_relaxed.*", "", prefix)
    prefix = re.sub(r"_rank.*", "", prefix)

    if prefix in seen_prefixes:
        continue

    seen_prefixes.add(prefix)
    out_name = f"{prefix}.pdb"
    shutil.copy2(pdb, out_dir / out_name)
    copied += 1

print(f"Found PDB files: {len(pdb_files)}")
print(f"Copied best/reference PDBs: {copied}")
print(f"Output directory: {out_dir}")