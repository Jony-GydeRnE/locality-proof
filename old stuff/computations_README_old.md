# Computations

Symbolic and numerical verification of the hidden-zero theorem for phi^3 amplitudes.

---

## Dependencies

- **SymPy** (`pip install sympy`) — for n=5 symbolic computations
- **NumPy** (`pip install numpy`) — for n=6 numerical computations

---

## Scripts

### `verify_tree_zeros.py`

**What:** Direct symbolic verification that A_5^tree vanishes on each of the 5 cyclic 1-zero loci.

**Relation to paper:** Verifies forward direction of Theorem 4.1 (A_5 satisfies hidden zeros).

**Runtime:** < 1 second.

```bash
python3 verify_tree_zeros.py
```

---

### `full_ansatz_n5.py`

**What:** Proves the converse at n=5: A_5 is the UNIQUE rational function of mass dimension -4 with simple poles that vanishes on all 5 cyclic 1-zero loci. Builds a 45-parameter ansatz, imposes 115 constraint equations, computes nullspace.

**Relation to paper:** Computational verification of Theorem 4.1 (Part I base case).

**Result:** Matrix rank 44/45, nullspace dimension 1, nullspace = A_5^tree.

**Runtime:** ~10-15 seconds.

```bash
python3 full_ansatz_n5.py
```

---

### `full_ansatz_n6.py`

**What:** Numerical verification of the hidden-zero theorem at n=6. Uses the 84-dimensional simple-poles ansatz and the broader 165-dimensional ansatz (with higher-order poles). Samples random kinematic points on each 1-zero locus, builds constraint matrix, computes nullspace via SVD.

**Relation to paper:** Verifies Proposition 5.4 (Part II, n=6 case).

**Result:** Both ansatz spaces yield nullspace dimension 1, spanned by A_6^tree (14 triangulations of the hexagon). Singular value gap exceeds 13 orders of magnitude. Stable across 4 random seeds.

**Runtime:** < 5 seconds.

```bash
python3 full_ansatz_n6.py
```

---

### `independence_check.py`

**What:** Determines how many 1-zero loci are needed to uniquely fix the amplitude. At n=5 (symbolic) and n=6 (numerical), tracks the constraint rank as loci are added progressively.

**Relation to paper:** Verifies Corollary 4.2 (n=5) and Remark 5.7 (independence count). Corrects an earlier conjecture: at n=6, five of six loci are needed (not floor(n/2)=3).

**Key results:**
| n | Ansatz dim | Loci needed | Redundant loci |
|---|---|---|---|
| 5 | 45 (all pairs) | 2 | 3 |
| 5 | 15 (planar only) | 3 | 2 |
| 6 | 84 | 5 | 1 |

**Runtime:** ~20 seconds.

```bash
python3 independence_check.py
```

---

### `verify_factorization.py`

**What:** Verifies BCFW factorization of A_n^tree at n=5,6: at each propagator pole X_{ij}=0, the residue factorizes as A_L * A_R, where the diagonal (i,j) splits the n-gon into two sub-polygons. The number of contributing terms equals C_{|L|-2} * C_{|R|-2} (product of Catalan numbers).

**Relation to paper:** Supports the factorization argument in Theorem 5.6, Step 4. Verifies factorization for A_tree; the remaining gap is whether general B in B_n shares this structure.

**Runtime:** < 2 seconds.

```bash
python3 verify_factorization.py
```

---

## How the scripts relate to the paper

| Paper section | Script | What is verified |
|---|---|---|
| Prop. 1.3 (A_5 vanishes on Z_k) | `verify_tree_zeros.py` | Forward direction |
| Thm. 4.1 (A_5 is unique, n=5) | `full_ansatz_n5.py` | Full theorem (base case) |
| Prop. 5.4 (n=6 verification) | `full_ansatz_n6.py` | Nullspace = 1, = A_6^tree |
| Cor. 4.2 / Remark 5.7 (independence) | `independence_check.py` | How many loci needed |
| Thm. 5.6, Step 4 (factorization) | `verify_factorization.py` | BCFW residue structure |

## Variable naming

The paper uses `X_{ij}` notation (planar Mandelstam variables). The n=5 script
`full_ansatz_n5.py` uses `s_ij` notation (2-particle Mandelstam invariants). The mapping is:

| Paper | Scripts (n=5) | Formula |
|---|---|---|
| X_{13} | s13 | (p1+p2)^2 |
| X_{14} | s14 | (p1+p2+p3)^2 |
| X_{24} | s24 | (p2+p3)^2 |
| X_{25} | s25 | (p2+p3+p4)^2 |
| X_{35} | s35 | (p3+p4)^2 |

The n=6 scripts use `X_{ij}` notation directly, matching the paper.
