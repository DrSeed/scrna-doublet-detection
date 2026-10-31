# Self-contained demo: 4 cell types, 8% doublets (some same-type), scored and called.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score

from doublets import doublet_scores, call_doublets


def simulate(n_cells=1500, n_genes=800, n_types=4, doublet_rate=0.08, seed=4):
    rng = np.random.default_rng(seed)
    base = rng.lognormal(-1.0, 1.2, n_genes)
    profiles = np.tile(base, (n_types, 1))
    for t in range(n_types):
        profiles[t, rng.choice(n_genes, 60, replace=False)] *= 8
    profiles /= profiles.sum(axis=1, keepdims=True)

    def cells(types):
        lib = rng.lognormal(np.log(3000), 0.35, len(types))
        mu = profiles[types] * lib[:, None]
        return rng.negative_binomial(5, 5 / (5 + mu))

    types = rng.integers(0, n_types, n_cells)
    counts = cells(types)
    kind = np.array(["singlet"] * n_cells, dtype=object)
    dbl = rng.choice(n_cells, int(doublet_rate * n_cells), replace=False)
    partner = rng.integers(0, n_types, len(dbl))
    counts[dbl] = counts[dbl] + cells(partner)
    kind[dbl] = np.where(partner == types[dbl], "homotypic doublet", "heterotypic doublet")
    return counts, types, kind


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    counts, types, kind = simulate()
    scores, pcs = doublet_scores(counts)
    called, threshold = call_doublets(scores, expected_rate=0.08)

    rows = []
    for k in ["heterotypic doublet", "homotypic doublet"]:
        m = (kind == k) | (kind == "singlet")
        rows.append({"class": k, "n": int((kind == k).sum()),
                     "auc_vs_singlets": round(roc_auc_score(kind[m] == k, scores[m]), 3),
                     "pct_detected": round(100 * called[kind == k].mean(), 1)})
    rows.append({"class": "singlet", "n": int((kind == "singlet").sum()), "auc_vs_singlets": None,
                 "pct_detected": round(100 * called[kind == "singlet"].mean(), 1)})
    table = pd.DataFrame(rows)
    table.to_csv("results/summary.csv", index=False)

    colours = {"singlet": "#cccccc", "homotypic doublet": "#1f77b4", "heterotypic doublet": "#d62728"}
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for k, c in colours.items():
        m = kind == k
        axes[0].scatter(pcs[m, 0], pcs[m, 1], s=5 if k == "singlet" else 12, color=c, label=k)
    axes[0].set_xlabel("PC1")
    axes[0].set_ylabel("PC2")
    axes[0].set_title("Heterotypic doublets sit between clusters")
    axes[0].legend(fontsize=8, markerscale=2)

    bins = np.linspace(0, 1, 31)
    for k, c in colours.items():
        axes[1].hist(scores[kind == k], bins=bins, color=c, alpha=0.7, density=True, label=k)
    axes[1].axvline(threshold, color="k", ls="--", label="call threshold")
    axes[1].set_xlabel("doublet score (fraction of simulated neighbours)")
    axes[1].set_ylabel("density")
    axes[1].set_title("Score distributions")
    axes[1].legend(fontsize=8)

    axes[2].bar(table["class"], table["pct_detected"], color=[colours[c] for c in table["class"]])
    axes[2].set_ylabel("% flagged as doublet")
    axes[2].set_title("Same-type doublets are much harder to catch")
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
