# RBP-Score: A Bioinformatics Pipeline for Phage Receptor-Binding Protein Discovery

This repository contains the code, datasets, and documentation developed for a bioinformatics project focused on detecting and ranking receptor-binding proteins (RBPs) in bacteriophages.

---

## Project Description

Bacteriophages (phages) infect bacteria with high specificity, a process largely mediated by receptor-binding proteins (RBPs). These proteins are responsible for recognizing and binding to bacterial surface receptors, making them essential for host specificity and valuable for applications such as diagnostics and phage engineering.

However, identifying RBPs computationally remains challenging due to their high sequence variability and poor conservation across phages.

The **RBP-Score** project addresses this limitation by implementing a reproducible and integrative pipeline that combines multiple sources of evidence into a single scoring system to prioritize candidate RBPs.

---

## Main Goals

- Build a curated dataset of experimentally validated RBPs  
- Detect potential RBPs in phage genomes  
- Evaluate candidates using different analytical approaches  
- Combine multiple evidence layers into a unified scoring system  
- Rank proteins based on their likelihood of being true RBPs  
- Assess the performance of the pipeline using reference data  

---

## Pipeline Overview

The workflow is designed to integrate complementary analyses:

- **Sequence-based search**  
  Identification of similar proteins using BLASTp against a curated RBP dataset  

- **Evolutionary analysis**  
  Multiple sequence alignment and distance calculation to evaluate divergence  

- **Structure-based comparison**  
  Structural similarity assessment using AlphaFold models and tools like DALI  

- **Scoring system**  
  Integration of all evidence into a global **RBP-Score**, allowing ranking of candidates  

---

## Computacional Tools

- Python  
- Snakemake  
- BLASTp  
- MUSCLE  
- ColabFold  
- Foldseek  
- FASTA files and protein datasets  

---

## Repository Organization

```bash
├── config/            # Configuration files
├── data/              # Input datasets (reference RBPs, query genomes)
├── docs/              # Reports and project documentation
│   └── manuscript/     # LaTeX files (.tex, .bib, figures)
├── results/           # Output results and scoring tables
├── scripts/           # Custom scripts for analysis and scoring
├── workflow/          # Snakemake pipeline files
└── README.md
