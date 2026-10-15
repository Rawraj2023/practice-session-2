"""
Chromatic Number of G(Z5)
==========================
Central digraph G(5, 1, 1) — natural central digraph on Z/5Z x Z/5Z
Vertex set: {(p, i) : p, i in Z/5Z}  (25 vertices)
Edge rule:  (p, i) -> ((p + i) mod 5, j)  for all j in Z/5Z

Reference: Rawat & Sehgal (2026)
"""

import numpy as np
import networkx as nx

k = 5
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

print(f"G(Z5) = G(5,1,1) underlying undirected graph")
print(f"  Vertices      : {n}")
print(f"  Edges         : {G.number_of_edges()}")
print(f"  Connected     : {nx.is_connected(G)}")
print(f"  Self-loops    : {self_loops}")
print(f"  Degree range  : {min(degrees)} – {max(degrees)}")

# ── Clique and independence numbers ──────────────────────────────────────────

omega = max(len(c) for c in nx.find_cliques(G))
alpha = max(len(s) for s in nx.find_cliques(nx.complement(G)))
print(f"  Max clique omega  : {omega}")
print(f"  Independence alpha: {alpha}")
print(f"  Fractional bound n/alpha: {n}/{alpha} = {n/alpha:.3f}  => lower bound ceil = {-(-n // alpha)}")

# ── Greedy coloring bounds ────────────────────────────────────────────────────

print("\nGreedy coloring upper bounds:")
for strategy in ["largest_first", "smallest_last"]:
    c = nx.coloring.greedy_color(G, strategy=strategy)
    print(f"  {strategy:20s}: {max(c.values()) + 1} colors")

# ── Exact chromatic number via backtracking ───────────────────────────────────

def is_k_colorable(k_col):
    """
    Backtracking coloring with degree-descending vertex order for pruning.
    Returns True if graph admits a proper k_col-coloring.
    """
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    col = [-1] * n
    def bt(idx):
        if idx == n:
            return True
        v = order[idx]
        used = {col[u] for u in adj[v] if col[u] >= 0}
        for c in range(k_col):
            if c not in used:
                col[v] = c
                if bt(idx + 1):
                    return True
                col[v] = -1
        return False
    return bt(0)

print("\nDetermining chromatic number (backtracking, degree-ordered):")
chi = None
for t in range(1, 10):
    result = is_k_colorable(t)
    print(f"  {t}-colorable: {result}")
    if result:
        chi = t
        break

print(f"\nChromatic number chi(G(Z5)) = {chi}")

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY — G(Z5) = G(5,1,1)")
print("=" * 60)
print(f"  Vertices         : {n}")
print(f"  Edges            : {G.number_of_edges()}")
print(f"  Connected        : {nx.is_connected(G)}")
print(f"  omega (clique)   : {omega}")
print(f"  alpha (indep.)   : {alpha}")
print(f"  Fractional chi   : {n/alpha:.3f}")
print(f"  chi(G(Z5))       : {chi}")
print()
print(f"  Observation: chi = {chi} = chi(G(Z4)). Note that chi = 5 does not")
print(f"  hold for all k >= 4: chi(G(Z7)) = 6 (see chromatic_z7.py).")
