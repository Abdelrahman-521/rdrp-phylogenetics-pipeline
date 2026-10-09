#!/usr/bin/env python3
"""Step 7: compare the ClustalW and T-Coffee results.

New after the course. The course version compared the two trees by eye in
FigTree. This step puts numbers on it:

* alignment stats: length, share of gaps, mean pairwise identity
* IQ-TREE stats: best-fit model, log-likelihood, informative sites
* tree distance: Robinson-Foulds distance between the two unrooted trees
* support: mean ultrafast-bootstrap support of the clades both trees agree
  on, versus the clades only one tree has
* a tanglegram (the two trees face to face, matching viruses joined by lines)
"""

from __future__ import annotations

import argparse
import itertools
import re
from io import StringIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import Phylo  # noqa: E402

from phylobot_utils import parse_fasta_text  # noqa: E402


# ---------- alignments ----------

def alignment_stats(path):
    records, _ = parse_fasta_text(Path(path).read_text())
    seqs = [s.upper() for _, s in records]
    length = len(seqs[0])
    if any(len(s) != length for s in seqs):
        raise ValueError(f"{path} is not an alignment (rows differ in length)")
    gaps = sum(s.count("-") for s in seqs) / (len(seqs) * length)
    identities = []
    for a, b in itertools.combinations(seqs, 2):
        both = [(x, y) for x, y in zip(a, b) if x != "-" and y != "-"]
        if both:
            identities.append(sum(x == y for x, y in both) / len(both))
    return {
        "sequences": len(seqs),
        "columns": length,
        "gap_fraction": gaps,
        "mean_pairwise_identity": sum(identities) / len(identities),
    }


# ---------- IQ-TREE report ----------

def iqtree_stats(path):
    text = Path(path).read_text()

    def find(pattern, cast=str):
        m = re.search(pattern, text)
        return cast(m.group(1)) if m else None

    return {
        "best_model": find(r"Best-fit model according to BIC:\s*(\S+)"),
        "log_likelihood": find(r"Log-likelihood of the tree:\s*(-?[\d.]+)", float),
        "informative_sites": find(r"Number of parsimony informative sites:\s*(\d+)", int),
        "iqtree_version": find(r"(?m)^IQ-TREE (?:multicore version )?(\d\S*)"),
    }


# ---------- trees ----------

def read_tree(path):
    return Phylo.read(StringIO(Path(path).read_text()), "newick")


def splits_with_support(tree, taxa):
    """Non-trivial bipartitions of an unrooted tree, keyed on the side without
    the alphabetically first taxon, mapped to bootstrap support."""
    ref = min(taxa)
    out = {}
    for clade in tree.find_clades():
        leaves = frozenset(t.name for t in clade.get_terminals())
        if 1 < len(leaves) < len(taxa) - 1:
            side = leaves if ref not in leaves else frozenset(taxa) - leaves
            support = clade.confidence
            if support is None and clade.name and re.fullmatch(r"[\d.]+", clade.name):
                support = float(clade.name)
            out[side] = support
    return out


def compare(tree_a, tree_b):
    taxa_a = {t.name for t in tree_a.get_terminals()}
    taxa_b = {t.name for t in tree_b.get_terminals()}
    if taxa_a != taxa_b:
        raise ValueError(f"trees have different taxa: {sorted(taxa_a ^ taxa_b)}")
    sa, sb = splits_with_support(tree_a, taxa_a), splits_with_support(tree_b, taxa_a)
    shared = sa.keys() & sb.keys()
    rf = len(sa.keys() ^ sb.keys())
    max_rf = 2 * (len(taxa_a) - 3)

    def mean(values):
        values = [v for v in values if v is not None]
        return sum(values) / len(values) if values else float("nan")

    def smaller_side(split):
        other = frozenset(taxa_a) - split
        return split if len(split) <= len(other) else other

    only_a = sorted(((sa[s], smaller_side(s)) for s in sa.keys() - shared), key=lambda x: -(x[0] or 0))
    only_b = sorted(((sb[s], smaller_side(s)) for s in sb.keys() - shared), key=lambda x: -(x[0] or 0))
    return {
        "taxa": len(taxa_a),
        "rf": rf,
        "rf_normalised": rf / max_rf,
        "shared_splits": len(shared),
        "splits_a": len(sa),
        "splits_b": len(sb),
        "support_shared_a": mean(sa[s] for s in shared),
        "support_shared_b": mean(sb[s] for s in shared),
        "support_only_a": mean(v for v, _ in only_a),
        "support_only_b": mean(v for v, _ in only_b),
        "only_a": only_a,
        "only_b": only_b,
    }


def sister_groups(tree):
    """Closest relatives of each virus in a midpoint-rooted tree."""
    tree.root_at_midpoint()
    out = {}
    for leaf in tree.get_terminals():
        parent = tree.get_path(leaf)[-2] if len(tree.get_path(leaf)) > 1 else tree.root
        out[leaf.name] = frozenset(t.name for t in parent.get_terminals()) - {leaf.name}
    return out


def changed_relatives(tree_a, tree_b):
    sa, sb = sister_groups(tree_a), sister_groups(tree_b)
    return {name: (sa[name], sb[name]) for name in sa if sa[name] != sb[name]}


# ---------- tanglegram ----------

def layout(tree, order_hint=None):
    """x = distance from root, y = leaf slot. Children are rotated so leaves
    line up with `order_hint` (the other tree's leaf order) where possible."""
    def key(clade):
        leaves = clade.get_terminals()
        if order_hint:
            return sum(order_hint.get(t.name, 0) for t in leaves) / len(leaves)
        return len(leaves)

    def sort(clade):
        clade.clades.sort(key=key)
        for child in clade.clades:
            sort(child)

    sort(tree.root)
    depths = tree.depths()
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    xs, ys = {}, {}
    leaves = tree.get_terminals()
    for i, leaf in enumerate(leaves):
        ys[leaf] = i
    for clade in tree.find_clades(order="postorder"):
        xs[clade] = depths[clade]
        if not clade.is_terminal():
            ys[clade] = (ys[clade.clades[0]] + ys[clade.clades[-1]]) / 2
    return xs, ys, [leaf.name for leaf in leaves]


def draw_tree(ax, tree, xs, ys, x0, width, mirror):
    span = max(xs.values()) or 1.0

    def X(clade):
        offset = xs[clade] / span * width
        return x0 - offset if mirror else x0 + offset

    for clade in tree.find_clades():
        for child in clade.clades:
            ax.plot([X(clade), X(clade)], [ys[clade], ys[child]], color="#444", lw=0.8)
            ax.plot([X(clade), X(child)], [ys[child], ys[child]], color="#444", lw=0.8)
            if not child.is_terminal() and child.confidence is not None and child.confidence < 95:
                ax.text(X(child), ys[child] + 0.25, f"{child.confidence:.0f}", fontsize=5.5,
                        color="#b23", ha="right" if not mirror else "left", va="bottom")
    return {leaf.name: (X(leaf), ys[leaf]) for leaf in tree.get_terminals()}


def tanglegram(tree_a, tree_b, label_a, label_b, out_png, conflicted=frozenset()):
    tree_a.root_at_midpoint()
    tree_b.root_at_midpoint()
    tree_a.ladderize()
    xs_a, ys_a, order_a = layout(tree_a)
    xs_b, ys_b, order_b = layout(tree_b, {name: i for i, name in enumerate(order_a)})
    # A couple of alternating passes untangle most of the crossings.
    for _ in range(2):
        xs_a, ys_a, order_a = layout(tree_a, {name: i for i, name in enumerate(order_b)})
        xs_b, ys_b, order_b = layout(tree_b, {name: i for i, name in enumerate(order_a)})

    n = len(order_a)
    fig, ax = plt.subplots(figsize=(11, max(6, n * 0.24)))
    tips_a = draw_tree(ax, tree_a, xs_a, ys_a, 0.0, 1.0, mirror=False)
    tips_b = draw_tree(ax, tree_b, xs_b, ys_b, 3.0, 1.0, mirror=True)
    for name, (x, y) in tips_a.items():
        ax.text(x + 0.03, y, name, fontsize=6.5, va="center")
    for name, (x, y) in tips_b.items():
        ax.text(x - 0.03, y, name, fontsize=6.5, va="center", ha="right")
    for name in order_a:
        ya, yb = tips_a[name][1], tips_b[name][1]
        hit = name in conflicted
        ax.plot([1.38, 1.62], [ya, yb], color="#d55" if hit else "#9bd", lw=0.9 if hit else 0.6)
    ax.set_xlim(-0.05, 3.05)
    ax.set_ylim(-1, n + 1)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title(f"{label_a} (left) vs {label_b} (right): midpoint-rooted ML trees\n"
                 f"red links: viruses whose closest relatives differ between the trees; red numbers: bootstrap support < 95",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ---------- report ----------

def fmt_clade(split, limit=6):
    names = sorted(split)
    text = ", ".join(names[:limit])
    return text + (f" (+{len(names) - limit} more)" if len(names) > limit else "")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clustal-aln", required=True)
    p.add_argument("--tcoffee-aln", required=True)
    p.add_argument("--clustal-tree", required=True)
    p.add_argument("--tcoffee-tree", required=True)
    p.add_argument("--clustal-report", required=True)
    p.add_argument("--tcoffee-report", required=True)
    p.add_argument("--prefix", default="comparison")
    args = p.parse_args()

    aln = {"ClustalW": alignment_stats(args.clustal_aln), "T-Coffee": alignment_stats(args.tcoffee_aln)}
    iq = {"ClustalW": iqtree_stats(args.clustal_report), "T-Coffee": iqtree_stats(args.tcoffee_report)}
    tree_c, tree_t = read_tree(args.clustal_tree), read_tree(args.tcoffee_tree)
    cmp_ = compare(tree_c, tree_t)
    moved = changed_relatives(read_tree(args.clustal_tree), read_tree(args.tcoffee_tree))
    conflicted = frozenset(moved)
    tanglegram(read_tree(args.clustal_tree), read_tree(args.tcoffee_tree),
               "ClustalW", "T-Coffee", f"{args.prefix}_tanglegram.png", conflicted)

    with open(f"{args.prefix}.tsv", "w") as out:
        out.write("metric\tClustalW\tT-Coffee\n")
        for key in ["sequences", "columns", "gap_fraction", "mean_pairwise_identity"]:
            out.write(f"{key}\t{aln['ClustalW'][key]}\t{aln['T-Coffee'][key]}\n")
        for key in ["best_model", "log_likelihood", "informative_sites", "iqtree_version"]:
            out.write(f"{key}\t{iq['ClustalW'][key]}\t{iq['T-Coffee'][key]}\n")
        out.write(f"robinson_foulds\t{cmp_['rf']}\t{cmp_['rf']}\n")
        out.write(f"robinson_foulds_normalised\t{cmp_['rf_normalised']:.4f}\t{cmp_['rf_normalised']:.4f}\n")

    a, t = aln["ClustalW"], aln["T-Coffee"]
    ia, it = iq["ClustalW"], iq["T-Coffee"]
    lines = [
        "# ClustalW vs T-Coffee: what changes in the RdRp tree?",
        "",
        f"{cmp_['taxa']} viral RdRp domains, aligned two ways, each alignment turned into a "
        "maximum-likelihood tree with IQ-TREE (ModelFinder + 1,000 ultrafast bootstraps).",
        "",
        "| | ClustalW | T-Coffee |",
        "|---|---|---|",
        f"| Alignment columns | {a['columns']} | {t['columns']} |",
        f"| Gaps | {a['gap_fraction']:.1%} | {t['gap_fraction']:.1%} |",
        f"| Mean pairwise identity | {a['mean_pairwise_identity']:.1%} | {t['mean_pairwise_identity']:.1%} |",
        f"| Parsimony-informative sites | {ia['informative_sites']} | {it['informative_sites']} |",
        f"| Best-fit model (BIC) | {ia['best_model']} | {it['best_model']} |",
        f"| Tree log-likelihood | {ia['log_likelihood']} | {it['log_likelihood']} |",
        f"| Mean support, clades both trees share | {cmp_['support_shared_a']:.1f} | {cmp_['support_shared_b']:.1f} |",
        f"| Mean support, clades only this tree has | {cmp_['support_only_a']:.1f} | {cmp_['support_only_b']:.1f} |",
        "",
        f"**Robinson-Foulds distance: {cmp_['rf']}** out of a possible {2 * (cmp_['taxa'] - 3)} "
        f"({cmp_['rf_normalised']:.0%}). The trees agree on {cmp_['shared_splits']} of "
        f"{cmp_['splits_a']} clades, and {len(moved)} of {cmp_['taxa']} viruses have a different "
        "closest relative (red links in the tanglegram).",
        "",
        "Log-likelihoods are not comparable across the two alignments (different columns), "
        "so the useful signals are the topology difference and where the support is weak.",
        "",
        "## Clades only one tree has (strongest first)",
        "",
        "| Tree | Support | Clade |",
        "|---|---|---|",
    ]
    for label, items in (("ClustalW", cmp_["only_a"]), ("T-Coffee", cmp_["only_b"])):
        for support, split in items[:5]:
            lines.append(f"| {label} | {support if support is not None else '-'} | {fmt_clade(split)} |")
    lines += ["", "## Viruses whose closest relatives change", "",
              "| Virus | ClustalW sister group | T-Coffee sister group |", "|---|---|---|"]
    for name in sorted(moved):
        left, right = moved[name]
        lines.append(f"| {name} | {fmt_clade(left, 4)} | {fmt_clade(right, 4)} |")
    lines += ["", f"![Tanglegram]({Path(args.prefix).name}_tanglegram.png)", ""]
    Path(f"{args.prefix}.md").write_text("\n".join(lines))
    print("\n".join(lines[:20]))


if __name__ == "__main__":
    main()
