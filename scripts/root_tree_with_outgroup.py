#!/usr/bin/env python3
from pathlib import Path
import argparse
from Bio import Phylo


def read_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().replace(" ", "_").split()[0])
    return ids


def main():
    parser = argparse.ArgumentParser(description="Root a Newick tree using one or more outgroups.")
    parser.add_argument("--tree", required=True)
    parser.add_argument("--outgroup", default="data/outgroup.fasta")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    tree = Phylo.read(args.tree, "newick")
    outgroup_ids = read_fasta_ids(Path(args.outgroup))
    terminals = {t.name: t for t in tree.get_terminals()}
    present = [terminals[x] for x in outgroup_ids if x in terminals]

    if not present:
        raise ValueError("No outgroup sequence was found in the tree. Check outgroup FASTA headers.")

    if len(present) == 1:
        tree.root_with_outgroup(present[0])
    else:
        tree.root_with_outgroup(tree.common_ancestor(present))

    outp = Path(args.output)
    outp.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(outp), "newick")
    print(f"Rooted tree written to: {outp}")
    print(f"Outgroups used: {len(present)}")

if __name__ == "__main__":
    main()
