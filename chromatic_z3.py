"""
Chromatic Polynomial and Chromatic Number of G(Z3)
====================================================
Central digraph G(3, 1, 1) — natural central digraph on Z/3Z x Z/3Z
Vertex set: {(p, i) : p, i in Z/3Z}  (9 vertices)
Edge rule:  (p, i) -> ((p + i) mod 3, j)  for all j in Z/3Z

Reference: Rawat & Sehgal (2026)
"""

import numpy as np

k = 3
n = k * k
labels = [f"({p},{i})" for p in range(k) for i in range(k)]

# ── Build adjacency matrix ────────────────────────────────────────────────────

A = np.zeros((n, n), dtype=int)
for p in range(k):
    for i in range(k):
        src = p * k + i
        for j in range(k):
            tp = (p + i) % k          # w = 1, so w^j = 1 for all j
            A[src, tp * k + j] = 1

print("Adjacency matrix A of G(Z3):")
print("     " + "  ".join(f"{l:>5}" for l in labels))
for r in range(n):
    row = "  ".join(f"{A[r, c]:>5}" for c in range(n))
    print(f"{labels[r]:>5} {row}")

self_loops = [labels[i] for i in range(n) if A[i, i] == 1]
print(f"\nSelf-loops (digraph): {self_loops}")

# ── Underlying undirected graph (no self-loops) ───────────────────────────────

U = np.zeros((n, n), dtype=int)
for i in range(n):
    for j in range(n):
        if i != j and (A[i, j] == 1 or A[j, i] == 1):
            U[i, j] = 1

degrees = list(U.sum(axis=1))
total_edges = int(U.sum()) // 2
print(f"\nUndirected graph: {n} vertices, {total_edges} edges")
print(f"Degree sequence: {degrees}")

# BFS connected components
visited = [False] * n
components = []
for start in range(n):
    if not visited[start]:
        comp, queue = [], [start]
        visited[start] = True
        while queue:
            v = queue.pop(0)
            comp.append(v)
            for u in range(n):
                if U[v, u] == 1 and not visited[u]:
                    visited[u] = True
                    queue.append(u)
        components.append(comp)
print(f"Components: {len(components)}, sizes: {[len(c) for c in components]}")

# ── Chromatic polynomial via deletion-contraction ─────────────────────────────
# Polynomials stored as dict {degree: coefficient}

def poly_sub(p1, p2):
    result = dict(p1)
    for deg, coef in p2.items():
        result[deg] = result.get(deg, 0) - coef
    return {d: c for d, c in result.items() if c != 0}

def poly_mul(p1, p2):
    result = {}
    for d1, c1 in p1.items():
        for d2, c2 in p2.items():
            d = d1 + d2
            result[d] = result.get(d, 0) + c1 * c2
    return result

def poly_eval(p, t):
    return sum(c * (t ** d) for d, c in p.items())

def poly_str(p):
    terms = []
    for deg in sorted(p.keys(), reverse=True):
        c = p[deg]
        if c == 0:
            continue
        sign = "+" if c > 0 else "-"
        ac = abs(c)
        if deg == 0:
            term = str(ac)
        elif deg == 1:
            term = f"{ac}t" if ac != 1 else "t"
        else:
            term = f"{ac}t^{deg}" if ac != 1 else f"t^{deg}"
        terms.append((sign, term))
    if not terms:
        return "0"
    result = ("-" + terms[0][1]) if terms[0][0] == "-" else terms[0][1]
    for sign, term in terms[1:]:
        result += f" {sign} {term}"
    return result

def chromatic_poly(adj, vertices):
    """Deletion-contraction algorithm. Returns polynomial as {degree: coeff}."""
    vertices = list(vertices)
    # Find first edge
    edge = None
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            if adj[vertices[i], vertices[j]] == 1:
                edge = (vertices[i], vertices[j])
                break
        if edge:
            break
    # Base case: no edges
    if edge is None:
        return {len(vertices): 1}
    u, v = edge
    # Deletion: remove edge u-v
    A_del = adj.copy()
    A_del[u, v] = A_del[v, u] = 0
    p_del = chromatic_poly(A_del, vertices)
    # Contraction: merge v into u
    A_con = adj.copy()
    new_verts = [x for x in vertices if x != v]
    for w in new_verts:
        if w != u and (A_con[v, w] == 1 or A_con[u, w] == 1):
            A_con[u, w] = A_con[w, u] = 1
    A_con[v, :] = A_con[:, v] = 0
    p_con = chromatic_poly(A_con, new_verts)
    return poly_sub(p_del, p_con)

import sys
sys.setrecursionlimit(100000)

print("\nComputing chromatic polynomial by deletion-contraction...")
component_polys = []
for comp in components:
    p = chromatic_poly(U.copy(), comp)
    component_polys.append(p)
    print(f"  Component {[labels[v] for v in comp]}:")
    print(f"    P = {poly_str(p)}")

# Multiply polynomials across components
total_poly = {0: 1}
for p in component_polys:
    total_poly = poly_mul(total_poly, p)

print(f"\nChromatic polynomial:")
print(f"  P(G(Z3), t) = {poly_str(total_poly)}")

# ── Evaluations ───────────────────────────────────────────────────────────────

print("\nEvaluations:")
for t in range(8):
    val = poly_eval(total_poly, t)
    marker = "  <-- first nonzero" if val > 0 and poly_eval(total_poly, t - 1) == 0 else ""
    print(f"  P(G(Z3), {t}) = {val}{marker}")

chi = next(t for t in range(1, 20) if poly_eval(total_poly, t) > 0)
print(f"\nChromatic number chi(G(Z3)) = {chi}")

# ── Summary ───────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print("SUMMARY — G(Z3) = G(3,1,1)")
print("=" * 60)
print(f"  Vertices      : {n}")
print(f"  Edges         : {total_edges}")
print(f"  Connected     : {'Yes' if len(components) == 1 else 'No'}")
print(f"  Self-loops    : {self_loops}")
print(f"  Degree seq.   : {sorted(set(degrees))} (counts {[(d, degrees.count(d)) for d in sorted(set(degrees))]})")
print(f"  chi(G(Z3))    : {chi}")
print(f"  P(G(Z3), t)   : {poly_str(total_poly)}")
