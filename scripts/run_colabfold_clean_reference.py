#!/usr/bin/env python3

from pathlib import Path
import subprocess
import shutil
import argparse
import os
import sys

def run_cmd(cmd, log_file):
    with open(log_file, "w") as log:
        process = subprocess.run(
            cmd,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )
    return process.returncode

def find_best_pdb(raw_dir):
    patterns = [
        "*rank_001*.pdb",
        "*rank_1*.pdb",
        "*ranked_0*.pdb",
        "*.pdb"
    ]

    for pattern in patterns:
        pdbs = sorted(raw_dir.rglob(pattern))
        if pdbs:
            return pdbs[0]

    return None

def run_colabfold_for_folder(input_dir, raw_dir, pdb_dir, logs_dir, label, args):
    fasta_files = sorted(Path(input_dir).glob("*.fasta"))

    if not fasta_files:
        print(f"[WARNING] No FASTA files found in {input_dir}")
        return []

    summary = []

    for fasta in fasta_files:
        sample = fasta.stem
        sample_raw = raw_dir / sample
        final_pdb = pdb_dir / f"{sample}.pdb"
        log_file = logs_dir / f"{label}_{sample}.log"

        if final_pdb.exists() and not args.overwrite:
            print(f"[SKIP] {sample}: PDB already exists")
            summary.append((sample, label, "skipped_existing_pdb", str(final_pdb)))
            continue

        print(f"[RUN] {label}: {sample}")

        sample_raw.mkdir(parents=True, exist_ok=True)

        cmd = [
            "colabfold_batch",
            "--num-models", str(args.num_models),
            "--num-recycle", str(args.num_recycle),
            "--model-type", "auto",
            "--msa-mode", args.msa_mode,
            str(fasta),
            str(sample_raw)
        ]

        status = run_cmd(cmd, log_file)

        if status != 0:
            print(f"[FAIL] {sample}: ColabFold failed. See {log_file}")
            summary.append((sample, label, "colabfold_failed", str(log_file)))
            continue

        best_pdb = find_best_pdb(sample_raw)

        if best_pdb is None:
            print(f"[FAIL] {sample}: no PDB found")
            summary.append((sample, label, "no_pdb_found", str(sample_raw)))
            continue

        shutil.copy2(best_pdb, final_pdb)
        print(f"[OK] {sample}: saved {final_pdb}")

        summary.append((sample, label, "pdb_created", str(final_pdb)))

        if args.clean_raw:
            shutil.rmtree(sample_raw)
            print(f"[CLEAN] removed raw folder: {sample_raw}")

    return summary

def main():
    parser = argparse.ArgumentParser(
        description="Run ColabFold for reference RBPs and outgroups, keep only best PDB, remove raw files."
    )

    parser.add_argument("--gpu", default=None, help="GPU index to use, e.g. 0, 1, 2, 3")
    parser.add_argument("--reference-dir", default="reference_structures/input_fastas")
    parser.add_argument("--outgroup-dir", default="reference_structures/outgroup_fastas")
    parser.add_argument("--raw-dir", default="reference_structures/colabfold_raw")
    parser.add_argument("--pdb-dir", default="reference_structures/pdb")
    parser.add_argument("--logs-dir", default="logs/colabfold_reference")
    parser.add_argument("--num-models", default=1, type=int)
    parser.add_argument("--num-recycle", default=1, type=int)
    parser.add_argument("--msa-mode", default="mmseqs2_uniref_env")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--clean-raw", action="store_true", default=True)

    args = parser.parse_args()

    if args.gpu is not None:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
        print(f"[GPU] Using CUDA_VISIBLE_DEVICES={args.gpu}")

    raw_dir = Path(args.raw_dir)
    pdb_dir = Path(args.pdb_dir)
    logs_dir = Path(args.logs_dir)

    raw_dir.mkdir(parents=True, exist_ok=True)
    pdb_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    all_summary = []

    all_summary += run_colabfold_for_folder(
        input_dir=args.reference_dir,
        raw_dir=raw_dir / "references",
        pdb_dir=pdb_dir,
        logs_dir=logs_dir,
        label="reference",
        args=args
    )

    all_summary += run_colabfold_for_folder(
        input_dir=args.outgroup_dir,
        raw_dir=raw_dir / "outgroups",
        pdb_dir=pdb_dir,
        logs_dir=logs_dir,
        label="outgroup",
        args=args
    )

    summary_file = pdb_dir / "reference_and_outgroup_pdb_summary.tsv"

    with open(summary_file, "w") as out:
        out.write("sample\ttype\tstatus\tpath\n")
        for row in all_summary:
            out.write("\t".join(row) + "\n")

    print("\nDONE")
    print(f"Summary: {summary_file}")
    print(f"PDB directory: {pdb_dir}")

if __name__ == "__main__":
    main()
