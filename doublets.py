# Scrublet-style doublet scoring: simulate doublets, embed, count simulated neighbours.
import numpy as np
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


def log_normalise(counts, target=1e4):
    lib = counts.sum(axis=1, keepdims=True)
    return np.log1p(counts / lib * target)


def simulate_doublets(counts, n, rng):
    i = rng.integers(0, len(counts), n)
    j = rng.integers(0, len(counts), n)
    return counts[i] + counts[j]


def doublet_scores(counts, sim_ratio=2.0, n_pcs=20, k=None, seed=0):
    """Score = fraction of a cell's nearest neighbours that are simulated doublets."""
    rng = np.random.default_rng(seed)
    n = len(counts)
    fake = simulate_doublets(counts, int(sim_ratio * n), rng)
    X = log_normalise(np.vstack([counts, fake]))
    mu, sd = X[:n].mean(axis=0), X[:n].std(axis=0) + 1e-8
    Z = (X - mu) / sd
    pcs = PCA(n_components=n_pcs, random_state=seed).fit(Z[:n]).transform(Z)
    k = k or max(5, int(round(0.5 * np.sqrt(n))))
    idx = NearestNeighbors(n_neighbors=k + 1).fit(pcs).kneighbors(pcs[:n], return_distance=False)[:, 1:]
    return (idx >= n).mean(axis=1), pcs[:n]


def call_doublets(scores, expected_rate=0.08):
    # Flag the top expected_rate fraction of cells, as you would with a known loading rate.
    threshold = np.quantile(scores, 1 - expected_rate)
    return scores > threshold, float(threshold)
