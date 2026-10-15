# practice-session-2

Computational study of the **chromatic number of the natural central digraph G(ℤₖ)** for k = 3 to 9.

Reference: Rawat & Sehgal (2026).

## The graph

For each k ≥ 3, the central digraph **G(k, 1, 1)** on ℤ/kℤ × ℤ/kℤ is defined by:

- **Vertices:** all pairs (p, i) with p, i ∈ ℤ/kℤ, so there are n = k² vertices.
- **Arcs:** (p, i) → ((p + i) mod k, j) for every j ∈ ℤ/kℤ.

The vertices (p, 0) have self-loops, because (p, 0) → (p, 0). Every computation below uses the **underlying undirected graph**: arc directions are ignored and self-loops are removed. "Colouring" means a proper vertex colouring of that graph.

## Main result

| k | \|V\| = k² | \|E\| | Degrees | ω (clique) | α (indep.) | n/α | **χ** |
|---|---|---|---|---|---|---|---|
| 3 | 9  | 21  | 4–5   | 3 | 3  | 3.00 | **4** |
| 4 | 16 | 54  | 6–7   | 3 | 4  | 4.00 | **5** |
| 5 | 25 | 110 | 8–9   | 3 | 6  | 4.17 | **5** |
| 6 | 36 | 195 | 10–11 | 3 | 9  | 4.00 | **5** |
| 7 | 49 | 315 | 12–13 | 3 | 12 | 4.08 | **6** |
| 8 | 64 | 476 | 14–15 | 3 | 16 | 4.00 | **6** |
| 9 | 81 | 684 | 16–17 | 3 | 20 | 4.05 | **6** |

**The conjecture χ(G(ℤₖ)) = 5 for all k ≥ 4 is false.** The values for k = 4, 5 and 6 suggested that χ settles at 5, but G(ℤ₇) is not 5-colourable, so χ(G(ℤ₇)) = 6.

Evidence for k = 7:

- An exhaustive DSATUR backtracking search (`chromatic_z7.py`) finds no proper 5-colouring. It visits 94,744 search nodes in about 2 seconds.
- The same search finds a proper 6-colouring, which the script checks edge by edge. Its colour classes have sizes 9, 10, 8, 7, 7, 8.
- An independent check with a SAT solver (CaDiCaL, via `python-sat`) gives the same answer. It also reproduces χ = 4, 5, 5, 5 for k = 3 to 6.

**k = 8 also has χ = 6.** `chromatic_z8.py` proves G(ℤ₈) is not 5-colourable (374,084 search nodes, about 12 s) and finds a valid 6-colouring with class sizes 13, 13, 11, 9, 9, 9. The SAT solver agrees. Here n/α = 64/16 = 4 gives only χ ≥ 4, so the jump to 6 cannot be explained by independence counting alone.

**k = 9 also has χ = 6.** `chromatic_z9.py` uses a SAT solver instead of DSATUR. The solver proves the 5-colouring problem has no solution (UNSAT) and returns a 6-colouring, which the script checks edge by edge. Its class sizes are 9, 16, 13, 15, 13, 15. So far χ = 6 for k = 7, 8 and 9.

## Properties seen for every k from 3 to 9

- **Connected:** the undirected graph is always connected.
- **ω = 3:** it always contains triangles but never a K₄. So the clique number gives only χ ≥ 3, far below the true values.
- **Near-regular:** every vertex has degree 2k − 2 or 2k − 1.
- **Edge count:** |E| = 21, 54, 110, 195, 315, 476, 684 for k = 3 to 9.
- **Self-loops** in the digraph are exactly at (p, 0) for p ∈ ℤ/kℤ.
- **Independence bound:** χ ≥ ⌈n/α⌉. This gives χ ≥ 5 for k = 5, 7 and 9, and χ ≥ 4 for k = 4, 6 and 8.
- **Greedy colouring is a poor guide:** greedy strategies need 5 to 10 colours, well above the true χ.

## Chromatic polynomials

**k = 3.** The full polynomial, computed exactly by deletion–contraction:

```
P(G(ℤ₃), t) = t⁹ − 21t⁸ + 196t⁷ − 1061t⁶ + 3625t⁵ − 7932t⁴ + 10722t³ − 8062t² + 2532t
```

Its values are P(0) through P(3) = 0, P(4) = 48, P(5) = 6,360, P(6) = 120,240 and P(7) = 1,062,600.

**k = 4.** Only some values, found by counting colourings directly: P(0) through P(4) = 0, P(5) = 41,040 and P(6) = 25,536,960. The full polynomial has degree 16 and has not been computed.

**k ≥ 5.** No polynomials have been computed. The graphs are too dense for plain deletion–contraction.

## Files

| File | Method | Run time* |
|---|---|---|
| `chromatic_z3.py` | Adjacency matrix, connected components, exact chromatic polynomial by deletion–contraction | ~1 s |
| `chromatic_z4.py` | Clique and independence numbers, greedy bounds, exact χ by backtracking, counts colourings for P(5) and P(6) | ~50 s |
| `chromatic_z5.py` | As above, with vertices ordered by degree (highest first) to prune the search | ~1 s |
| `chromatic_z6.py` | Backtracking with forward checking, then checks the 5-colouring it found | ~2 s |
| `chromatic_z7.py` | Exhaustive DSATUR search with colour-symmetry breaking, then checks the 6-colouring it found | ~3 s |
| `chromatic_z8.py` | Same search as `chromatic_z7.py`, for k = 8 | ~50 s |
| `chromatic_z9.py` | SAT encoding solved with CaDiCaL (`python-sat`), triangle fixed to colours 0–2 for symmetry breaking, then checks the 6-colouring it returns | ~1 s |
| `Chromatic Invariants of G(Zk).html` | Summary report: table of invariants, the refuted conjecture, polynomials | n/a |

\*Measured on a single CPU core. None of the scripts use the GPU.

### How the k = 7 and k = 8 searches work

The fixed-order backtracking used for k ≤ 6 does not scale to 49 or 64 vertices, so `chromatic_z7.py` and `chromatic_z8.py` add two techniques:

1. **DSATUR branching.** The next vertex to colour is the uncoloured vertex whose neighbours already use the most distinct colours. Ties go to the vertex with the higher degree. This finds contradictions early.
2. **Colour-symmetry breaking.** Colours are interchangeable, so a vertex may use any colour already in play or only the *next* unused colour. This removes up to K! equivalent branches.

Both techniques keep the search exhaustive. A "not K-colourable" result is therefore a proof, not a heuristic.

### How the k = 9 SAT check works

At 81 vertices the pure-Python DSATUR search becomes too slow, so `chromatic_z9.py` writes "G is K-colourable" as a Boolean formula. Variable x(v, c) means vertex v gets colour c. The clauses say every vertex has exactly one colour and adjacent vertices never share one. One triangle is fixed to colours 0, 1, 2, which is safe because colours are interchangeable. If the CaDiCaL solver answers UNSAT, that proves no K-colouring exists; if it answers SAT, it returns a colouring, which the script checks independently.

## Running

Requirements: Python 3, `numpy` and `networkx`. `chromatic_z9.py` also needs `python-sat`.

```bash
pip install numpy networkx python-sat
python chromatic_z3.py   # or chromatic_z4.py ... chromatic_z9.py
```

## Changes in the second revision

- **Added `chromatic_z7.py`.** It settles the k = 7 case, previously listed as "5 or 6", at χ = 6.
- **Fixed `chromatic_z6.py`.** The final check used to test a greedy colouring with 8 colours. It now tests the 5-colouring found by the search.
- **Fixed `chromatic_z4.py`.** The summary used to print a hard-coded P(5) = 41040. It now prints the computed values.
- **Fixed `chromatic_z3.py`.** Degrees now print as plain integers instead of `np.int64(...)`.
- **Fixed the HTML report.** The degree claim was wrong: it said "k−1 or k" where the correct values are "2k−2 or 2k−1". The report also now has the k = 7 result, a conjecture section stating that the conjecture fails at k = 7, and P(6) for ℤ₄.

## Changes in the third revision

- **Added `chromatic_z8.py`.** It confirms χ(G(ℤ₈)) = 6.

## Changes in the fourth revision

- **Added `chromatic_z9.py`.** It uses a SAT solver and confirms χ(G(ℤ₉)) = 6.

## Open questions

- Find χ(G(ℤₖ)) for k ≥ 10, and work out whether χ stays at 6 or keeps growing.
- Explain structurally why G(ℤ₇), G(ℤ₈) and G(ℤ₉) need 6 colours when their largest clique has only 3 vertices.
- Compute the full chromatic polynomial for k = 4.
