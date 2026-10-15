"""
Chromatic Number of G(Z8)
==========================
Central digraph G(8, 1, 1) — natural central digraph on Z/8Z x Z/8Z
Vertex set: {(p, i) : p, i in Z/8Z}  (64 vertices)
Edge rule:  (p, i) -> ((p + i) mod 8, j)  for all j in Z/8Z

Uses the same exhaustive DSATUR search as chromatic_z7.py: branch on the
vertex with the most distinct neighbour colours, and let a vertex open only
the next unused colour (colour-symmetry breaking). A "False" result is a
proof of non-colourability.

Result: G(Z8) is NOT 5-colourable, so chi(G(Z8)) = 6 = chi(G(Z7)).
Unlike k = 7, the bound n/alpha = 64/16 = 4 gives only chi >= 4 here.

Reference: Rawat & Sehgal (2026)
"""

import time
import numpy as np
import networkx as nx

k = 8
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
G.add_nodes_from(range(n))
for i in range(n):
    for j in range(i + 1, n):
        if A[i, j] == 1 or A[j, i] == 1:
            G.add_edge(i, j)

adj = [set(G[v]) for v in range(n)]
degrees = [len(adj[i]) for i in range(n)]

print(f"G(Z8) = G(8,1,1) underlying undirected graph")
print(f"  Vertices      : {n}")
print(f"  Edges         : {G.number_of_edges()}")
print(f"  Connected     : {nx.is_connected(G)}")
print(f"  Self-loops    : {self_loops}")
print(f"  Degree range  : {min(degrees)} – {max(degrees)}")

# ── Clique and independence numbers ──────────────────────────────────────────

omega = max(len(c) for c in nx.find_cliques(G))
_, alpha = nx.max_weight_clique(nx.complement(G), weight=None)
print(f"  Max clique omega  : {omega}")
print(f"  Independence alpha: {alpha}")
print(f"  Fractional bound n/alpha: {n}/{alpha} = {n/alpha:.3f}  => lower bound ceil = {-(-n // alpha)}")

# ── Greedy coloring bounds ────────────────────────────────────────────────────

print("\nGreedy coloring upper bounds:")
for strategy in ["largest_first", "smallest_last", "DSATUR"]:
    c = nx.coloring.greedy_color(G, strategy=strategy)
    print(f"  {strategy:20s}: {max(c.values()) + 1} colors")

# ── Exact chromatic number via DSATUR backtracking ────────────────────────────

def find_coloring(k_col):
    """
    Exhaustive DSATUR backtracking with colour-symmetry breaking.
    Returns (coloring, nodes_searched); coloring is None if no proper
    k_col-coloring exists.
    """
    col = [-1] * n
    nodes = [0]

    def pick():
        # Uncoloured vertex with max saturation, ties broken by degree
        best, best_key = None, None
        for v in range(n):
            if col[v] < 0:
                sat = len({col[u] for u in adj[v] if col[u] >= 0})
                key = (sat, degrees[v])
                if best_key is None or key > best_key:
                    best, best_key = v, key
        return best

    def bt(used_colors):
        nodes[0] += 1
        v = pick()
        if v is None:
            return True
        forbidden = {col[u] for u in adj[v] if col[u] >= 0}
        # Colours are interchangeable, so trying one brand-new colour suffices
        for c in range(min(k_col, used_colors + 1)):
            if c not in forbidden:
                col[v] = c
                if bt(max(used_colors, c + 1)):
                    return True
                col[v] = -1
        return False

    found = bt(0)
    return (list(col) if found else None), nodes[0]

print("\nDetermining chromatic number (DSATUR backtracking, exhaustive):")
chi, coloring = None, None
for t in range(1, 10):
    start = time.time()
    coloring, searched = find_coloring(t)
    print(f"  {t}-colorable: {coloring is not None}  ({searched} search nodes, {time.time() - start:.2f}s)")
    if coloring is not None:
        chi = t
        break

print(f"\nChromatic number chi(G(Z8)) = {chi}")

# ── Cross-check: verify the coloring found by the search ─────────────────────

print(f"\nVerifying the {chi}-coloring found by the search:")
valid = all(coloring[u] != coloring[v] for u, v in G.edges())
class_sizes = [coloring.count(c) for c in range(chi)]
print(f"  Colors used          : {len(set(coloring))}")
print(f"  Valid proper coloring: {valid}")
print(f"  Color class sizes    : {class_sizes}  (each <= alpha = {alpha})")

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY — G(Z8) = G(8,1,1)")
print("=" * 60)
print(f"  Vertices         : {n}")
print(f"  Edges            : {G.number_of_edges()}")
print(f"  Connected        : {nx.is_connected(G)}")
print(f"  omega (clique)   : {omega}")
print(f"  alpha (indep.)   : {alpha}")
print(f"  Fractional chi   : {n/alpha:.3f}")
print(f"  chi(G(Z8))       : {chi}")
print()
if chi == 6:
    print(f"  Observation: chi = 6 = chi(G(Z7)), so chi(G(Zk)) = 6 for k = 7, 8")
    print(f"  (vs. 5 for k = 4, 5, 6).")
else:
    print(f"  Observation: chi = {chi}, differs from chi(G(Z7)) = 6.")
