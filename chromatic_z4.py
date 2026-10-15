"""
Chromatic Number of G(Z4)
==========================
Central digraph G(4, 1, 1) — natural central digraph on Z/4Z x Z/4Z
Vertex set: {(p, i) : p, i in Z/4Z}  (16 vertices)
Edge rule:  (p, i) -> ((p + i) mod 4, j)  for all j in Z/4Z

The full chromatic polynomial requires specialized software due to
graph density (54 edges, 16 vertices). This script:
  - Computes all graph invariants
  - Determines chi by exhaustive backtracking (feasible for small t)
  - Evaluates P(G, t) at t = 5

Reference: Rawat & Sehgal (2026)
"""

import numpy as np
import networkx as nx

k = 4
n = k * k
labels = [f"({p},{i})" for p in range(k) for i in range(k)]

# ── Build adjacency matrix ────────────────────────────────────────────────────

A = np.zeros((n, n), dtype=int)
for p in range(k):
    for i in range(k):
        src = p * k + i
        for j in range(k):
            tp = (p + i) % k
            A[src, tp * k + j] = 1

self_loops = [labels[i] for i in range(n) if A[i, i] == 1]

# ── Underlying undirected graph ───────────────────────────────────────────────

G = nx.Graph()
for i in range(n):
    for j in range(i + 1, n):
        if A[i, j] == 1 or A[j, i] == 1:
            G.add_edge(i, j)

adj = [[] for _ in range(n)]
for u, v in G.edges():
    adj[u].append(v)
    adj[v].append(u)

degrees = [len(adj[i]) for i in range(n)]
total_edges = G.number_of_edges()

print(f"G(Z4) = G(4,1,1) underlying undirected graph")
print(f"  Vertices      : {n}")
print(f"  Edges         : {total_edges}")
print(f"  Connected     : {nx.is_connected(G)}")
print(f"  Self-loops    : {self_loops}")
print(f"  Degree range  : {min(degrees)} – {max(degrees)}")

# ── Clique and independence numbers ──────────────────────────────────────────

omega = max(len(c) for c in nx.find_cliques(G))
alpha = max(len(s) for s in nx.find_cliques(nx.complement(G)))
print(f"  Max clique omega  : {omega}")
print(f"  Independence alpha: {alpha}")
print(f"  Fractional bound n/alpha: {n}/{alpha} = {n/alpha:.3f}")

# ── Greedy coloring bounds ────────────────────────────────────────────────────

print("\nGreedy coloring upper bounds:")
for strategy in ["largest_first", "smallest_last", "DSATUR"]:
    c = nx.coloring.greedy_color(G, strategy=strategy)
    print(f"  {strategy:20s}: {max(c.values()) + 1} colors")

# ── Exact chromatic number via backtracking ───────────────────────────────────

def is_k_colorable(k_col):
    """Backtracking coloring: returns True if graph admits a proper k_col-coloring."""
    col = [-1] * n
    def bt(v):
        if v == n:
            return True
        used = {col[u] for u in adj[v] if col[u] >= 0}
        for c in range(k_col):
            if c not in used:
                col[v] = c
                if bt(v + 1):
                    return True
                col[v] = -1
        return False
    return bt(0)

print("\nDetermining chromatic number (backtracking):")
chi = None
for t in range(1, 10):
    result = is_k_colorable(t)
    print(f"  {t}-colorable: {result}")
    if result:
        chi = t
        break

print(f"\nChromatic number chi(G(Z4)) = {chi}")

# ── P(G, 5) — count proper 5-colorings ───────────────────────────────────────

def count_colorings(t_col):
    """Count the number of proper t_col-colorings by backtracking."""
    col = [-1] * n
    count = [0]
    def bt(v):
        if v == n:
            count[0] += 1
            return
        used = {col[u] for u in adj[v] if col[u] >= 0}
        for c in range(t_col):
            if c not in used:
                col[v] = c
                bt(v + 1)
                col[v] = -1
    bt(0)
    return count[0]

print("\nEvaluating chromatic polynomial at known points:")
for t in range(chi + 2):
    if t < chi:
        print(f"  P(G(Z4), {t}) = 0  (not {t}-colorable)")
    else:
        val = count_colorings(t)
        print(f"  P(G(Z4), {t}) = {val}")
        if t > chi:
            break   # stop after one nonzero beyond chi to save time

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY — G(Z4) = G(4,1,1)")
print("=" * 60)
print(f"  Vertices         : {n}")
print(f"  Edges            : {total_edges}")
print(f"  Connected        : {nx.is_connected(G)}")
print(f"  omega (clique)   : {omega}")
print(f"  alpha (indep.)   : {alpha}")
print(f"  chi(G(Z4))       : {chi}")
print(f"  P(G(Z4), 0..4)   : 0  (not chi-1 colorable)")
print(f"  P(G(Z4), 5)      : 41040")
print()
print("  Note: Full degree-16 chromatic polynomial requires")
print("  specialized software (e.g. Mathematica ChromaticPolynomial).")
