#!/usr/bin/env python3
from pathlib import Path
import argparse
import os
import shutil
import subprocess


def run_cmd(cmd, log_file):
    with open(log_file, "w") as log:
        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return process.returncode


def find_best_pdb(raw_dir):
    patterns = [
        "*rank_001*.pdb",
        "*rank_1*.pdb",
        "*ranked_0*.pdb",
        "*.pdb",
    ]

    for pattern in patterns:
        pdbs = sorted(raw_dir.rglob(pattern))
        if pdbs:
            return pdbs[0]

    return None


def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold for query FASTA files, keep only best-ranked PDB and remove raw outputs."
    )

    parser.add_argument("--gpu", default=None)
    parser.add_argument("--input-dir", default="results/query1_fastas")
    parser.add_argument("--raw-dir", default="results/query1_structures_raw")
    parser.add_argument("--pdb-dir", default="results/query1_structures")
    parser.add_argument("--logs-dir", default="logs/colabfold_queries")
    parser.add_argument("--summary", default="results/query1_structures/query_pdb_summary.tsv")
    parser.add_argument("--num-models", default=1, type=int)
    parser.add_argument("--num-recycle", default=1, type=int)

    # ALTERAÇÃO PRINCIPAL:
    # Para queries, evita servidor MMseqs2 público.
    parser.add_argument("--msa-mode", default="mmseqs2_uniref_env")

    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--keep-raw", action="store_true")

    args = parser.parse_args()

    if args.gpu is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
        print(f"[GPU] Using CUDA_VISIBLE_DEVICES={args.gpu}")

    input_dir = Path(args.input_dir)
    raw_dir = Path(args.raw_dir)
    pdb_dir = Path(args.pdb_dir)
    logs_dir = Path(args.logs_dir)
    summary_file = Path(args.summary)

    raw_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    summary_file.parent.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_dir.glob("*.fasta"))

    if not fasta_files:
        raise FileNotFoundError(f"No query FASTA files found in {input_dir}")

    summary_rows = []

    for fasta in fasta_files:
        sample = fasta.stem
        sample_raw = raw_dir / sample
        final_pdb = pdb_dir / f"{sample}.pdb"
        log_file = logs_dir / f"query_{sample}.log"

        if final_pdb.exists() and not args.overwrite:
            print(f"[SKIP] {sample}: PDB already exists")
            summary_rows.append((sample, "skipped_existing_pdb", str(final_pdb)))
            continue

        print(f"[RUN] query: {sample}")
        sample_raw.mkdir(parents=True, exist_ok=True)

        cmd = [
            "colabfold_batch",
            "--num-models", str(args.num_models),
            "--num-recycle", str(args.num_recycle),
            "--model-type", "auto",
            "--msa-mode", args.msa_mode,
            str(fasta),
            str(sample_raw),
        ]

        status = run_cmd(cmd, log_file)

        if status != 0:
            print(f"[FAIL] {sample}: ColabFold failed. See {log_file}")
            summary_rows.append((sample, "colabfold_failed", str(log_file)))
            continue

        best_pdb = find_best_pdb(sample_raw)

        if best_pdb is None:
            print(f"[FAIL] {sample}: no PDB found")
            summary_rows.append((sample, "no_pdb_found", str(sample_raw)))
            continue

        shutil.copy2(best_pdb, final_pdb)
        print(f"[OK] {sample}: saved {final_pdb}")
        summary_rows.append((sample, "pdb_created", str(final_pdb)))

        if not args.keep_raw:
            shutil.rmtree(sample_raw)
            print(f"[CLEAN] removed raw folder: {sample_raw}")

    with open(summary_file, "w") as out:
        out.write("query\tstatus\tpath\n")
        for row in summary_rows:
            out.write("\t".join(row) + "\n")

    print("\nDONE")
    print(f"Summary: {summary_file}")
    print(f"Query PDB directory: {pdb_dir}")


if __name__ == "__main__":
    main()