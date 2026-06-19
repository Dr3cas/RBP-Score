# RBP-Score

**RBP-Score** is a reproducible bioinformatic pipeline for the identification and prioritisation of candidate phage receptor-binding proteins (RBPs).

The pipeline integrates three complementary evidence layers:

1. **Sequence similarity** using BLASTp;
2. **Phylogenetic context** using multiple sequence alignment and tree-based distance;
3. **Structural similarity** using predicted protein structures and Foldseek.

These evidence layers are combined into a weighted prioritisation score, with an additional outgroup-based penalty using non-RBP capsid/major head proteins as negative controls.

> **Important:** RBP-Score is currently a prototype and should be interpreted as a candidate prioritisation framework, not as a fully calibrated RBP classifier.

---

## Background

Bacteriophages infect bacteria with high host specificity. This specificity is largely mediated by receptor-binding proteins (RBPs), which recognise bacterial surface receptors and are therefore central to phage host range, diagnostics, and phage engineering.

Computational identification of RBPs remains difficult because these proteins are often highly diverse, modular and poorly conserved at sequence level. RBP-Score addresses this challenge by combining sequence, phylogenetic and structure-level evidence rather than relying on a single annotation strategy.

---

## Main features

- Uses a curated reference dataset of experimentally validated RBPs.
- Incorporates non-RBP capsid/major head proteins as negative outgroups.
- Builds a BLASTp database containing reference RBPs and outgroups.
- Cleans and standardises query protein identifiers.
- Performs sequence-based screening with BLASTp.
- Builds a global phylogenetic input FASTA containing references, queries and outgroups.
- Aligns sequences and infers a global phylogenetic tree.
- Roots the tree using the outgroup common ancestor.
- Predicts query protein structures with ColabFold.
- Compares predicted query structures against reference and outgroup structures using Foldseek.
- Calculates raw, outgroup-penalised and effective component scores.
- Produces final RBP-Score rankings and detailed per-query outputs.

---

## Repository structure

```text
RBP-Score/
├── config/
│   └── config.yaml
├── data/
│   ├── RBP_validated_final.xlsx
│   ├── blast_db_with_outgroups.fasta
│   ├── outgroup.fasta
│   ├── query.fasta
│   ├── query_clean.fasta
│   └── rbp_db_annotated.fasta
├── db/
│   └── rbp_db.*
├── docs/
│   ├── manuscript/
│   │   ├── main.tex
│   │   ├── references.bib
│   │   ├── global_tree_rooted.png
│   │   └── eval_transform.png
│   ├── RBPScore_FINAL.pptx
│   ├── RBPScore_Final_Submission_19-06.pdf
│   ├── RBP_Score_pipeline_for_phage_RBPs_1.pdf
│   ├── Script_Projeto.odt
│   └── Tema Projeto em Bioinformatica_SBS(2026).pdf
├── envs/
│   └── rbpscore.yaml
├── reference_structures/
│   ├── input_fastas/
│   └── pdb/
├── results/
│   ├── alignments/
│   ├── database/
│   ├── foldseek/
│   ├── phylogeny/
│   ├── query_fastas/
│   ├── query_structures/
│   ├── rbpscore/
│   ├── trees/
│   ├── trees_png/
│   └── blast_results.tsv
├── scripts/
│   ├── build_global_phylogeny_input.py
│   ├── build_rbp_db_with_morphotype.py
│   ├── calculate_rbpscore.py
│   ├── clean_query_fasta.py
│   ├── collect_best_pdbs.py
│   ├── fetch_outgroups.py
│   ├── render_tree_png.py
│   ├── root_tree_with_outgroup.py
│   ├── run_colabfold_clean_queries.py
│   ├── run_colabfold_clean_reference.py
│   ├── split_outgroup_fasta.py
│   ├── split_query_fasta.py
│   ├── split_reference_fasta.py
│   └── validate_outgroup_fasta.py
├── workflow/
│   └── snakefile
└── README.md
```

---

## Input data

The main input files are stored in `data/`.

| File | Description |
|---|---|
| `RBP_validated_final.xlsx` | Curated table of experimentally validated RBP references and metadata. |
| `rbp_db_annotated.fasta` | Annotated reference RBP FASTA generated from the curated database. |
| `outgroup.fasta` | Final non-RBP capsid/major head outgroup set. |
| `blast_db_with_outgroups.fasta` | FASTA containing reference RBPs and outgroups for BLASTp searches. |
| `query.fasta` | Input query protein sequences. |
| `query_clean.fasta` | Cleaned query FASTA with standardised identifiers. |

Additional alternative outgroup FASTA files are included in `data/`, reflecting exploratory outgroup selection during development.

---

## Configuration

Main parameters are defined in:

```text
config/config.yaml
```

The configuration includes:

```yaml
blast:
  evalue_cutoff: 1e-5
  pident_cutoff: 30.0
  coverage_cutoff: 0.50

blast_score_weights:
  identity: 0.35
  coverage: 0.30
  bitscore: 0.25
  evalue: 0.10

structural_score_weights:
  qtmscore: 0.50
  alntmscore: 0.30
  rmsd: 0.20

final_rbpscore_weights:
  sequence: 0.35
  phylogeny: 0.25
  structural: 0.40

outgroup_penalty:
  blast: 0.10
  phylogeny: 0.15
  structural: 0.20

foldseek:
  evalue: 10
  max_seqs: 1000
  alignment_type: 1

colabfold:
  gpu: 0
  num_models: 1
  num_recycle: 1
  msa_mode: mmseqs2_uniref_env

outgroup_filters:
  min_identity: 20.0
  min_length: 10
  max_evalue: 1000.0

normalization:
  evalue_max_log: 150.0
```

These weights are configurable and exploratory in the current prototype.

---

## Conda environment

The repository includes a Conda environment file in:

```text
envs/rbpscore.yaml
```

This environment is intended to make the pipeline easier to reproduce by installing the main command-line and Python dependencies used by the workflow.

To create the environment from the repository root, run:

```bash
conda env create -f envs/rbpscore.yaml
```

Then activate it with:

```bash
conda activate rbpscore
```

The current environment file defines:

```yaml
name: rbpscore

channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.11
  - blast
  - muscle
  - fasttree
  - foldseek
  - hmmer
  - biopython
  - pandas
  - pyyaml
  - matplotlib
  - networkx
  - pip

  - pip:
      - snakemake
```

This environment includes Python 3.11, BLAST, MUSCLE, FastTree, Foldseek, HMMER, Biopython, pandas, PyYAML, matplotlib, NetworkX and Snakemake.

> **Note:** ColabFold is used by the workflow for query structure prediction, but it is not installed by the current `envs/rbpscore.yaml` file. If running the structural prediction step, make sure ColabFold is installed separately or available in the execution environment.

---

## Pipeline overview

The Snakemake workflow is defined in:

```text
workflow/snakefile
```

A simplified representation of the workflow is:

```text
Query protein FASTA
        │
        ▼
Curated RBP references + non-RBP capsid outgroups
        │
        ▼
Evidence layers
  ├── BLASTp sequence similarity
  ├── MUSCLE/FastTree phylogenetic context
  └── ColabFold/Foldseek structural similarity
        │
        ▼
Outgroup penalty + weighted RBP-Score
        │
        ▼
Final ranking + partial score inspection
```

---

## Scoring logic

RBP-Score combines sequence, phylogenetic and structural evidence.

### Sequence score

The sequence score integrates:

- BLASTp percentage identity;
- query coverage;
- bitscore relative to the best RBP hit for the query;
- transformed BLASTp e-value.

### Phylogenetic score

Phylogenetic support is calculated from tree distance:

```text
S_phylo = 1 / (1 + d)
```

where `d` is the tree distance between a query protein and a reference RBP.

### Structural score

The structural score combines:

- query TM-score;
- aligned-region TM-score;
- RMSD transformed into a bounded similarity score.

### Outgroup penalty

For each evidence layer, the best matching non-RBP outgroup is used as negative evidence:

```text
effective score = max(0, raw score - penalty weight × best outgroup score)
```

### Final score

The final score is calculated from the effective component scores:

```text
RBP-Score = 0.35 × sequence + 0.25 × phylogeny + 0.40 × structural
```

The final score should be interpreted as a prioritisation score, not as a probability of RBP function.

---

## Running the workflow

First create and activate the Conda environment:

```bash
conda env create -f envs/rbpscore.yaml
conda activate rbpscore
```

From the repository root, run the workflow with:

```bash
snakemake -s workflow/snakefile --cores 1
```

For a dry run:

```bash
snakemake -s workflow/snakefile --cores 1 -n
```

For a more complete run using all available CPU cores:

```bash
snakemake -s workflow/snakefile --cores all
```

The Conda environment provides the main dependencies required by the workflow, including BLAST, MUSCLE, FastTree, Foldseek, HMMER and Snakemake. ColabFold should be installed or made available separately if the structure prediction rules are executed.

---

## Main outputs

The main RBP-Score outputs are stored in:

```text
results/rbpscore/
```

| Output | Description |
|---|---|
| `final_rbpscore_summary.tsv` | Final per-query ranking with RBP-Score and key support metrics. |
| `detailed_rbpscore.tsv` | Detailed query-reference score table. |
| `details/` | Per-query detailed score outputs. |

Other relevant outputs include:

| Directory / file | Description |
|---|---|
| `results/blast_results.tsv` | BLASTp query-reference and query-outgroup results. |
| `results/alignments/` | Multiple sequence alignments. |
| `results/phylogeny/` | Global phylogeny input and outgroup validation files. |
| `results/trees/` | Newick tree files. |
| `results/trees_png/` | Rendered phylogenetic tree images. |
| `results/foldseek/` | Structural comparison results. |
| `results/query_structures/` | Query structure prediction summaries. |
| `results/database/` | Reference database metadata generated during the workflow. |

---

## Current interpretation

The current repository represents a functional prototype/test run. It demonstrates that the pipeline can:

- integrate curated RBP references and non-RBP outgroups;
- combine sequence, phylogenetic and structural evidence;
- penalise candidates with stronger similarity to non-RBP controls;
- rank query proteins using an interpretable score;
- retain partial scores for manual inspection.

The current implementation is not yet a fully validated classifier. Larger positive and negative validation sets are required to estimate thresholds, tune weights and evaluate sensitivity/specificity.

---

## Documentation

The `docs/` directory contains:

- project reports;
- presentation material;
- manuscript files and supporting documents.

The `docs/manuscript/` directory is intended for LaTeX source files, bibliography files and manuscript figures.

---

## Limitations

- The scoring weights are exploratory and configurable.
- The non-RBP outgroup set is intentionally narrow and based on homologous capsid/major head proteins.
- Structural similarity may be informative for divergent proteins, but can also produce non-specific support when proteins share common architectural features.
- The current score should not be interpreted as a calibrated probability.
- Manual inspection of partial sequence, phylogenetic and structural scores remains essential.

---

## Suggested citation

If using or adapting this repository, cite it as:

```text
Gomes A, Santos S. RBP-Score: An integrative bioinformatic pipeline for the identification and prioritization of phage receptor-binding proteins. University of Minho, 2026.
```

---

## Author

**André Gomes**  
Centre of Biological Engineering, University of Minho  
Braga, Portugal
