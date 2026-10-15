"""
Chromatic Number of G(Z9)
==========================
Central digraph G(9, 1, 1) — natural central digraph on Z/9Z x Z/9Z
Vertex set: {(p, i) : p, i in Z/9Z}  (81 vertices)
Edge rule:  (p, i) -> ((p + i) mod 9, j)  for all j in Z/9Z

The pure-Python DSATUR search used for k = 7, 8 is too slow at 81 vertices,
so this script encodes "G is K-colourable" as a SAT problem and hands it to
the CaDiCaL solver (python-sat). An UNSAT answer is a proof that no proper
K-colouring exists; a SAT answer comes with a colouring, which is then
checked edge by edge in plain Python.

Encoding (variable x[v][c] = "vertex v gets colour c"):
  - every vertex gets at least one colour
  - every vertex gets at most one colour
  - adjacent vertices never share a colour
  - symmetry breaking: the vertices of one triangle are fixed to colours
    0, 1, 2 (sound because colours are interchangeable)

Result: G(Z9) is NOT 5-colourable, so chi(G(Z9)) = 6 = chi(G(Z7)) = chi(G(Z8)).

Requires: pip install python-sat

Reference: Rawat & Sehgal (2026)
"""

import time
import itertools
import numpy as np
import networkx as nx
from pysat.solvers import Cadical153

k = 9
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

degrees = [G.degree(v) for v in range(n)]

print(f"G(Z9) = G(9,1,1) underlying undirected graph")
print(f"  Vertices      : {n}")
print(f"  Edges         : {G.number_of_edges()}")
print(f"  Connected     : {nx.is_connected(G)}")
print(f"  Self-loops    : {self_loops}")
print(f"  Degree range  : {min(degrees)} – {max(degrees)}")

# ── Clique and independence numbers ──────────────────────────────────────────

cliques = list(nx.find_cliques(G))
omega = max(len(c) for c in cliques)
triangle = next(sorted(c) for c in cliques if len(c) == omega)
_, alpha = nx.max_weight_clique(nx.complement(G), weight=None)
print(f"  Max clique omega  : {omega}")
print(f"  Independence alpha: {alpha}")
print(f"  Fractional bound n/alpha: {n}/{alpha} = {n/alpha:.3f}  => lower bound ceil = {-(-n // alpha)}")

# ── Greedy coloring bounds ────────────────────────────────────────────────────

print("\nGreedy coloring upper bounds:")
for strategy in ["largest_first", "smallest_last", "DSATUR"]:
    c = nx.coloring.greedy_color(G, strategy=strategy)
    print(f"  {strategy:20s}: {max(c.values()) + 1} colors")

# ── Exact chromatic number via SAT ────────────────────────────────────────────

def find_coloring(k_col):
    """
    Decide k_col-colourability with a SAT solver.
    Returns (coloring, clause_count); coloring is None if the formula is UNSAT.
    """
    var = lambda v, c: v * k_col + c + 1
    clauses = []
    for v in range(n):
        clauses.append([var(v, c) for c in range(k_col)])
        for a, b in itertools.combinations(range(k_col), 2):
            clauses.append([-var(v, a), -var(v, b)])
    for u, v in G.edges():
        for c in range(k_col):
            clauses.append([-var(u, c), -var(v, c)])
    # Symmetry breaking: pin the triangle to the first colours
    for c, v in enumerate(triangle[:k_col]):
        clauses.append([var(v, c)])

    with Cadical153(bootstrap_with=clauses) as solver:
        if not solver.solve():
            return None, len(clauses)
        true_vars = {x for x in solver.get_model() if x > 0}
    coloring = [next(c for c in range(k_col) if var(v, c) in true_vars) for v in range(n)]
    return coloring, len(clauses)

print(f"\nDetermining chromatic number (SAT solver: CaDiCaL 1.5.3):")
chi, coloring = None, None
for t in range(omega, 10):     # fewer than omega colours is impossible
    start = time.time()
    coloring, num_clauses = find_coloring(t)
    status = "SAT" if coloring is not None else "UNSAT"
    print(f"  {t}-colorable: {coloring is not None}  ({status}, {num_clauses} clauses, {time.time() - start:.2f}s)")
    if coloring is not None:
        chi = t
        break

print(f"\nChromatic number chi(G(Z9)) = {chi}")

# ── Cross-check: verify the coloring returned by the solver ──────────────────

print(f"\nVerifying the {chi}-coloring returned by the solver:")
valid = all(coloring[u] != coloring[v] for u, v in G.edges())
class_sizes = [coloring.count(c) for c in range(chi)]
print(f"  Colors used          : {len(set(coloring))}")
print(f"  Valid proper coloring: {valid}")
print(f"  Color class sizes    : {class_sizes}  (each <= alpha = {alpha})")

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY — G(Z9) = G(9,1,1)")
print("=" * 60)
print(f"  Vertices         : {n}")
print(f"  Edges            : {G.number_of_edges()}")
print(f"  Connected        : {nx.is_connected(G)}")
print(f"  omega (clique)   : {omega}")
print(f"  alpha (indep.)   : {alpha}")
print(f"  Fractional chi   : {n/alpha:.3f}")
print(f"  chi(G(Z9))       : {chi}")
print()
if chi == 6:
    print(f"  Observation: chi = 6 = chi(G(Z7)) = chi(G(Z8)), so chi(G(Zk)) = 6")
    print(f"  for k = 7, 8, 9 (vs. 5 for k = 4, 5, 6).")
else:
    print(f"  Observation: chi = {chi}, differs from chi(G(Z7)) = chi(G(Z8)) = 6.")
