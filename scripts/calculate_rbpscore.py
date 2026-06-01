#!/usr/bin/env python3
from pathlib import Path
import argparse
import csv
import math
import yaml
from Bio import Phylo

DEFAULT_CONFIG = {
    "blast": {"evalue_cutoff": 1e-5, "pident_cutoff": 30.0, "coverage_cutoff": 0.50},
    "blast_score_weights": {"identity": 0.35, "coverage": 0.30, "bitscore": 0.25, "evalue": 0.10},
    "structural_score_weights": {"qtmscore": 0.50, "alntmscore": 0.30, "rmsd": 0.20},
    "final_rbpscore_weights": {"sequence": 0.35, "phylogeny": 0.25, "structural": 0.40},
    "outgroup_penalty": {"blast": 0.20, "phylogeny": 0.20, "structural": 0.20},
    "outgroup_filters": {"min_identity": 20.0, "min_length": 10, "max_evalue": 1000.0},
}


def clamp_score(value, min_value=0.0, max_value=1.0):
    """
    Keep any score-like value within a closed interval, by default [0, 1].
    This prevents partial scores such as coverage, TM-score-derived values,
    structural scores or final RBP-Scores from exceeding their intended scale.
    """
    try:
        if value in (None, "", "NA", "nan"):
            return min_value
        value = float(value)
    except Exception:
        return min_value

    if math.isnan(value):
        return min_value

    return max(min_value, min(value, max_value))


def load_config(path):
    if path is None or not Path(path).exists():
        return DEFAULT_CONFIG
    with open(path) as f:
        user_cfg = yaml.safe_load(f) or {}
    cfg = {section: dict(values) for section, values in DEFAULT_CONFIG.items()}
    for section, values in user_cfg.items():
        if isinstance(values, dict):
            cfg.setdefault(section, {})
            cfg[section].update(values)
        else:
            cfg[section] = values
    return cfg


def safe_id(x):
    parts = str(x).split("|")
    parts = [p for p in parts if not p.lower().startswith("label=")]
    x = "|".join(parts)

    return (
        x.replace("|", "_")
        .replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
    )



def display_query_id(x):
    """Return a clean query identifier for TSV outputs, without label= metadata."""
    parts = str(x).split("|")
    parts = [p for p in parts if not p.lower().startswith("label=")]
    return "|".join(parts)


def normalize_seq_id(x):
    x = str(x).strip().split()[0]
    x = Path(x).name
    if x.lower().endswith(".pdb"):
        x = x[:-4]
    return safe_id(x)


def normalize_pdb_id(x):
    return normalize_seq_id(x)


def is_outgroup_id(seq_id):
    return "OUTGROUP" in str(seq_id).upper()


def passes_outgroup_blast_filters(pident, length, evalue, filters):
    return (
        pident >= float(filters["min_identity"])
        and length >= float(filters["min_length"])
        and evalue <= float(filters["max_evalue"])
    )


def read_fasta_headers(path):
    headers = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                header = line[1:].strip().replace(" ", "_")
                seq_id = header.split()[0]
                headers[seq_id] = header
    return headers


def query_label(query_id, header):
    h = header.lower()
    q = query_id.lower()
    if "label=positive" in h:
        return "positive_control"
    if "label=negative" in h:
        return "negative_control"
    if "label=random" in h:
        return "random_candidate"
    if "label=pharbp" in h:
        return "pharbp_candidate"
    if "label=outgroup" in h:
        return "outgroup"
    if "random" in h or "phage_test" in h:
        return "random_candidate"
    return "unknown_candidate"


def safe_float(value, default=None):
    try:
        if value in (None, "", "NA", "nan"):
            return default
        return float(value)
    except Exception:
        return default


def evalue_to_score(evalue):
    if evalue == 0:
        return 1.0
    return clamp_score((-math.log10(evalue)) / 50.0)


def phylogeny_score(distance):
    if distance is None:
        return 0.0
    return clamp_score(1.0 / (1.0 + float(distance)))


def penalise(raw_score, outgroup_score, penalty_weight):
    return clamp_score(float(raw_score) - float(penalty_weight) * float(outgroup_score))


def structural_score(qtmscore, alntmscore, rmsd, weights):
    qtm = clamp_score(safe_float(qtmscore, 0.0))
    atm = clamp_score(safe_float(alntmscore, 0.0))
    r = safe_float(rmsd, None)
    rmsd_score = 0.0 if r is None else clamp_score(1.0 / (1.0 + max(r, 0.0)))

    score = (
        float(weights["qtmscore"]) * qtm
        + float(weights["alntmscore"]) * atm
        + float(weights["rmsd"]) * rmsd_score
    )
    return clamp_score(score)


def load_structural_hits(foldseek_file, query_safe_to_id, db_safe_to_id, outgroup_safe_to_id, structural_weights):
    pair_hits = {}
    outgroup_hits_by_query = {}
    if not foldseek_file or not Path(foldseek_file).exists():
        return pair_hits, outgroup_hits_by_query

    expected_cols = ["query", "target", "evalue", "bits", "alntmscore", "qtmscore", "ttmscore", "rmsd"]
    with open(foldseek_file) as f:
        first_line = f.readline().strip()
        f.seek(0)
        if first_line.startswith("query\ttarget"):
            reader = csv.DictReader(f, delimiter="\t")
        else:
            reader = (dict(zip(expected_cols, row)) for row in csv.reader(f, delimiter="\t") if row and len(row) >= len(expected_cols))
        for row in reader:
            q_safe = normalize_pdb_id(row["query"])
            t_safe = normalize_pdb_id(row["target"])
            q = query_safe_to_id.get(q_safe)
            if q is None:
                continue
            score = structural_score(row["qtmscore"], row["alntmscore"], row["rmsd"], structural_weights)
            hit = {
                "foldseek_target": t_safe,
                "foldseek_evalue": safe_float(row["evalue"], "NA"),
                "foldseek_bits": safe_float(row["bits"], "NA"),
                "alntmscore": safe_float(row["alntmscore"], "NA"),
                "qtmscore": safe_float(row["qtmscore"], "NA"),
                "ttmscore": safe_float(row["ttmscore"], "NA"),
                "rmsd": safe_float(row["rmsd"], "NA"),
                "structural_score": score,
            }
            if t_safe in db_safe_to_id:
                s = db_safe_to_id[t_safe]
                key = (q, s)
                if key not in pair_hits or score > pair_hits[key]["structural_score"]:
                    pair_hits[key] = hit
            elif t_safe in outgroup_safe_to_id or is_outgroup_id(t_safe):
                outgroup_id = outgroup_safe_to_id.get(t_safe, t_safe)
                current = outgroup_hits_by_query.get(q)
                if current is None or score > current["structural_score"]:
                    out_hit = dict(hit)
                    out_hit["outgroup_id"] = outgroup_id
                    outgroup_hits_by_query[q] = out_hit
    return pair_hits, outgroup_hits_by_query


def main():
    parser = argparse.ArgumentParser(description="Calculate integrated RBP-Score from BLAST, phylogeny and Foldseek, with outgroup penalty.")
    parser.add_argument("--blast", required=True)
    parser.add_argument("--query", required=True)
    parser.add_argument("--db", required=True)
    parser.add_argument("--tree", required=True)
    parser.add_argument("--foldseek", required=False, default=None)
    parser.add_argument("--outgroup", required=False, default="data/outgroup.fasta")
    parser.add_argument("--summary", required=True)
    parser.add_argument("--detailed", required=True)
    parser.add_argument("--details-dir", required=True)
    parser.add_argument("--config", default="config/config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    blast_cfg = cfg["blast"]
    blast_weights = cfg["blast_score_weights"]
    structural_weights = cfg["structural_score_weights"]
    final_weights = cfg["final_rbpscore_weights"]
    outgroup_penalty = cfg["outgroup_penalty"]
    outgroup_filters = cfg["outgroup_filters"]

    query_headers = read_fasta_headers(Path(args.query))
    db_headers = read_fasta_headers(Path(args.db))
    outgroup_headers = read_fasta_headers(Path(args.outgroup)) if Path(args.outgroup).exists() else {}

    query_ids = list(query_headers.keys())
    db_ids = list(db_headers.keys())
    outgroup_ids = list(outgroup_headers.keys())

    query_safe_to_id = {normalize_seq_id(q): q for q in query_ids}
    db_safe_to_id = {normalize_seq_id(s): s for s in db_ids}
    outgroup_safe_to_id = {normalize_seq_id(o): o for o in outgroup_ids}

    tree = Phylo.read(args.tree, "newick")
    terminals = {t.name: t for t in tree.get_terminals()}
    terminals_safe = {normalize_seq_id(t.name): t for t in tree.get_terminals()}

    def get_distance(a, b):
        if a in terminals and b in terminals:
            return tree.distance(terminals[a], terminals[b])
        a_safe = normalize_seq_id(a)
        b_safe = normalize_seq_id(b)
        if a_safe in terminals_safe and b_safe in terminals_safe:
            return tree.distance(terminals_safe[a_safe], terminals_safe[b_safe])
        return None

    blast_pairs = {}
    with open(args.blast) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            q_norm = normalize_seq_id(row["qseqid"])
            s_norm = normalize_seq_id(row["sseqid"])
            q_original = query_safe_to_id.get(q_norm)
            if q_original is None:
                continue
            pident = float(row["pident"])
            length = float(row["length"])
            qlen = float(row["qlen"])
            evalue = float(row["evalue"])
            bitscore = float(row["bitscore"])
            coverage = clamp_score(length / qlen) if qlen > 0 else 0.0
            is_outgroup_hit = is_outgroup_id(s_norm)

            passes_positive_rbp_filters = (
                evalue <= float(blast_cfg["evalue_cutoff"])
                and pident >= float(blast_cfg["pident_cutoff"])
                and coverage >= float(blast_cfg["coverage_cutoff"])
            )

            passes_negative_outgroup_filters = (
                is_outgroup_hit
                and passes_outgroup_blast_filters(pident, length, evalue, outgroup_filters)
            )

            if passes_positive_rbp_filters or passes_negative_outgroup_filters:
                key = (q_original, s_norm)
                hit = {
                    "original_subject_id": row["sseqid"],
                    "normalised_subject_id": s_norm,
                    "pident": pident,
                    "coverage": coverage,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "is_outgroup_hit": is_outgroup_hit,
                }
                if key not in blast_pairs or bitscore > blast_pairs[key]["bitscore"]:
                    blast_pairs[key] = hit

    max_bitscore_by_query = {}
    for (q, s_norm), hit in blast_pairs.items():
        if s_norm in db_safe_to_id:
            max_bitscore_by_query[q] = max(max_bitscore_by_query.get(q, 0.0), hit["bitscore"])

    structural_pairs, best_outgroup_by_query = load_structural_hits(args.foldseek, query_safe_to_id, db_safe_to_id, outgroup_safe_to_id, structural_weights)

    all_rows = []
    summary_rows = []
    details_dir = Path(args.details_dir)
    details_dir.mkdir(parents=True, exist_ok=True)
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.detailed).parent.mkdir(parents=True, exist_ok=True)

    for q in query_ids:
        per_query = []
        # Best BLAST outgroup hit for this query.
        # Used only as negative evidence / penalty.
        outgroup_blast_hits = []
        for (q_hit, s_norm), hit in blast_pairs.items():
            if q_hit == q and (s_norm in outgroup_safe_to_id or is_outgroup_id(s_norm)):
                outgroup_id = outgroup_safe_to_id.get(s_norm, hit.get("original_subject_id", s_norm))
                outgroup_score = clamp_score(hit["pident"] / 100.0)
                outgroup_blast_hits.append((outgroup_score, outgroup_id))

        if outgroup_blast_hits:
            best_outgroup_blast_score, best_outgroup_blast_id = max(
                outgroup_blast_hits,
                key=lambda x: x[0],
            )
        else:
            best_outgroup_blast_score, best_outgroup_blast_id = 0.0, "NA"

        # Best phylogenetic outgroup for this query.
        # This identifies the non-RBP outgroup closest to the query in the tree.
        outgroup_phylo_hits = []
        for out in outgroup_ids:
            dist_to_out = get_distance(q, out)
            if dist_to_out is not None:
                outgroup_phylo_hits.append((phylogeny_score(dist_to_out), out, dist_to_out))

        if outgroup_phylo_hits:
            best_outgroup_phylogeny_score, best_outgroup_phylogeny_id, best_outgroup_phylogeny_distance = max(
                outgroup_phylo_hits,
                key=lambda x: x[0],
            )
        else:
            best_outgroup_phylogeny_score, best_outgroup_phylogeny_id, best_outgroup_phylogeny_distance = 0.0, "NA", "NA"

        # Best structural outgroup for this query.
        outgroup_hit = best_outgroup_by_query.get(q)
        best_outgroup_structural_score = clamp_score(outgroup_hit["structural_score"]) if outgroup_hit else 0.0
        best_outgroup_structural_id = outgroup_hit["outgroup_id"] if outgroup_hit else "NA"

        for s in db_ids:
            s_norm = normalize_seq_id(s)
            hit = blast_pairs.get((q, s_norm))
            if hit:
                identity_score = clamp_score(hit["pident"] / 100.0)
                coverage_score = hit["coverage"]
                max_bits = max_bitscore_by_query.get(q, 1.0)
                bitscore_score = clamp_score(hit["bitscore"] / max_bits) if max_bits > 0 else 0.0
                ev_score = evalue_to_score(hit["evalue"])
                raw_blast_score = clamp_score(
                    float(blast_weights["identity"]) * clamp_score(identity_score)
                    + float(blast_weights["coverage"]) * clamp_score(coverage_score)
                    + float(blast_weights["bitscore"]) * clamp_score(bitscore_score)
                    + float(blast_weights["evalue"]) * clamp_score(ev_score)
                )
                blast_present = "yes"
                pident, coverage, evalue, bitscore = hit["pident"], hit["coverage"], hit["evalue"], hit["bitscore"]
            else:
                raw_blast_score = 0.0
                blast_present = "no"
                pident = coverage = evalue = bitscore = "NA"

            effective_blast_score = penalise(raw_blast_score, best_outgroup_blast_score, float(outgroup_penalty["blast"]))
            dist = get_distance(q, s)
            raw_phylogeny_score = phylogeny_score(dist)
            effective_phylogeny_score = penalise(raw_phylogeny_score, best_outgroup_phylogeny_score, float(outgroup_penalty["phylogeny"]))

            structural_hit = structural_pairs.get((q, s))
            if structural_hit:
                struct_present = "yes"
                raw_structural_score = clamp_score(structural_hit["structural_score"])
                foldseek_evalue, foldseek_bits = structural_hit["foldseek_evalue"], structural_hit["foldseek_bits"]
                alntmscore, qtmscore, ttmscore, rmsd = structural_hit["alntmscore"], structural_hit["qtmscore"], structural_hit["ttmscore"], structural_hit["rmsd"]
            else:
                struct_present = "no"
                raw_structural_score = 0.0
                foldseek_evalue = foldseek_bits = alntmscore = qtmscore = ttmscore = rmsd = "NA"

            effective_structural_score = penalise(raw_structural_score, best_outgroup_structural_score, float(outgroup_penalty["structural"]))
            rbp_score = clamp_score(
                float(final_weights["sequence"]) * effective_blast_score
                + float(final_weights["phylogeny"]) * effective_phylogeny_score
                + float(final_weights["structural"]) * effective_structural_score
            )

            row = {
                "query_id": display_query_id(q), "reference_rbp_id": s, "blast_present": blast_present,
                "pident": pident, "coverage": coverage, "evalue": evalue, "bitscore": bitscore,
                "raw_blast_score": raw_blast_score, "best_blast_outgroup": best_outgroup_blast_id,
                "outgroup_blast_score": best_outgroup_blast_score, "effective_blast_score": effective_blast_score,
                "tree_distance": dist if dist is not None else "NA", "raw_phylogeny_score": raw_phylogeny_score,
                "best_phylogeny_outgroup": best_outgroup_phylogeny_id,
                "outgroup_phylogeny_distance": best_outgroup_phylogeny_distance,
                "outgroup_phylogeny_score": best_outgroup_phylogeny_score, "effective_phylogeny_score": effective_phylogeny_score,
                "structure_present": struct_present, "foldseek_evalue": foldseek_evalue, "foldseek_bits": foldseek_bits,
                "alntmscore": alntmscore, "qtmscore": qtmscore, "ttmscore": ttmscore, "rmsd": rmsd,
                "raw_structural_score": raw_structural_score, "best_structural_outgroup": best_outgroup_structural_id,
                "outgroup_structural_score": best_outgroup_structural_score,
                "effective_structural_score": effective_structural_score, "rbp_score": rbp_score,
            }
            per_query.append(row)
            all_rows.append(row)

        per_query.sort(key=lambda r: float(r["rbp_score"]), reverse=True)
        detail_path = details_dir / f"{safe_id(q)}.details.tsv"
        with open(detail_path, "w", newline="") as out:
            writer = csv.DictWriter(out, fieldnames=list(per_query[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(per_query)

        blast_rows = [r for r in per_query if r["blast_present"] == "yes"]
        struct_rows = [r for r in per_query if r["structure_present"] == "yes"]
        best_blast = max(blast_rows, key=lambda r: float(r["effective_blast_score"])) if blast_rows else None
        best_phylo = max(per_query, key=lambda r: float(r["effective_phylogeny_score"]))
        best_struct = max(struct_rows, key=lambda r: float(r["effective_structural_score"])) if struct_rows else None
        best_combined = per_query[0]

        summary_rows.append({
            "query_id": display_query_id(q), "query_label": query_label(q, query_headers[q]),
            "best_blast_reference": best_blast["reference_rbp_id"] if best_blast else "NA",
            "raw_blast_score": best_blast["raw_blast_score"] if best_blast else 0.0,
            "best_blast_outgroup": best_outgroup_blast_id,
            "outgroup_blast_score": best_outgroup_blast_score,
            "effective_blast_score": best_blast["effective_blast_score"] if best_blast else 0.0,
            "best_phylogeny_reference": best_phylo["reference_rbp_id"], "tree_distance": best_phylo["tree_distance"],
            "raw_phylogeny_score": best_phylo["raw_phylogeny_score"],
            "best_phylogeny_outgroup": best_outgroup_phylogeny_id,
            "outgroup_phylogeny_distance": best_outgroup_phylogeny_distance,
            "outgroup_phylogeny_score": best_outgroup_phylogeny_score,
            "effective_phylogeny_score": best_phylo["effective_phylogeny_score"],
            "best_structural_reference": best_struct["reference_rbp_id"] if best_struct else "NA",
            "raw_structural_score": best_struct["raw_structural_score"] if best_struct else 0.0,
            "best_structural_outgroup": best_outgroup_structural_id,
            "outgroup_structural_score": best_outgroup_structural_score,
            "effective_structural_score": best_struct["effective_structural_score"] if best_struct else 0.0,
            "best_combined_reference": best_combined["reference_rbp_id"], "rbp_score": best_combined["rbp_score"],
            "n_valid_blast_hits": len(blast_rows), "n_structural_hits": len(struct_rows),
        })

    all_rows.sort(key=lambda r: (r["query_id"], -float(r["rbp_score"])))
    with open(args.detailed, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=list(all_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    summary_rows.sort(key=lambda r: float(r["rbp_score"]), reverse=True)
    for i, row in enumerate(summary_rows, start=1):
        row["rank"] = i

    summary_fields = [
        "rank", "query_id", "query_label",
        "best_blast_reference", "raw_blast_score", "best_blast_outgroup", "outgroup_blast_score",
        "effective_blast_score",
        "best_phylogeny_reference", "tree_distance", "raw_phylogeny_score",
        "best_phylogeny_outgroup", "outgroup_phylogeny_distance", "outgroup_phylogeny_score",
        "effective_phylogeny_score",
        "best_structural_reference", "raw_structural_score",
        "best_structural_outgroup", "outgroup_structural_score", "effective_structural_score",
        "best_combined_reference", "rbp_score", "n_valid_blast_hits", "n_structural_hits",
    ]
    with open(args.summary, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Summary written to: {args.summary}")
    print(f"Detailed global table written to: {args.detailed}")
    print(f"Per-query details written to: {args.details_dir}")


if __name__ == "__main__":
    main()
