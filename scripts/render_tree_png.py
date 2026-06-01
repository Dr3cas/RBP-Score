#!/usr/bin/env python3
from pathlib import Path
import argparse
from Bio import Phylo
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Render a Newick tree as PNG.")
    parser.add_argument("tree")
    parser.add_argument("output")
    parser.add_argument("--title", default="Global phylogenetic tree")
    args = parser.parse_args()

    tree_file = Path(args.tree)
    output_png = Path(args.output)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    tree = Phylo.read(str(tree_file), "newick")
    n = len(tree.get_terminals())
    fig_height = max(12, n * 0.23)
    fig_width = 20

    fig = plt.figure(figsize=(fig_width, fig_height))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
    ax.set_title(args.title)
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Tree PNG written to: {output_png}")

if __name__ == "__main__":
    main()
