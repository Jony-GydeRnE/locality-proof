# Rodina Locality Proof

A purely algebraic proof that Rodina's hidden-zero conditions uniquely determine
the phi^3 tree amplitudes. No Feynman diagrams. No D-subsets. Just linear algebra.

**Paper:** Rodina, arXiv:2406.04234v5 — "Hidden zeros are equivalent to enhanced
ultraviolet scaling and lead to unique amplitudes in Tr(phi^3) theory"

---

## What this proves

**Theorem (general n):** Let B be a rational function with at most simple
poles in the planar Mandelstam variables, homogeneous of mass dimension -(n-3).
If B vanishes on each of the n cyclic 1-zero loci Z_1,...,Z_n, then

  B = c * A_n^tree

where A_n^tree is the phi^3 tree amplitude and c is a scalar.

**In words:** The 1-zero conditions *uniquely fix* the amplitude. Locality (only
Feynman poles survive) is *derived* from the 1-zeros, not assumed.

---

## Structure

```
paper/
  main.tex           Full formal proof, Part I (n=5) + Part II (general n)
  part2_induction.tex  Part II: induction, n=6 verification, independence
  references.bib     Bibliography

computations/
  verify_tree_zeros.py    Forward direction: A_5 = 0 on each Z_k
  full_ansatz_n5.py       Hard direction at n=5: A_5 is UNIQUE (symbolic)
  full_ansatz_n6.py       Hard direction at n=6: A_6 is UNIQUE (numerical)
  independence_check.py   How many loci are needed at n=5, n=6
  verify_factorization.py BCFW residue factorization at n=5, n=6
  README.md               Script documentation

working-notes/
  claude_proof_exploration.md    Exploratory session notes
  rodina_1zero_analysis.md       Computation analysis
  session-results.md             Session log (2026-04-15)
  homogeneity_lemma_for_claude.pdf  Reference PDF
```

---

## Quick start

```bash
pip install sympy numpy

# n=5: verify A_5 vanishes on all loci (< 1 second)
python3 computations/verify_tree_zeros.py

# n=5: verify uniqueness via nullspace computation (~10-15 seconds)
python3 computations/full_ansatz_n5.py

# n=6: verify uniqueness numerically (< 5 seconds)
python3 computations/full_ansatz_n6.py

# Independence: how many loci needed (~20 seconds)
python3 computations/independence_check.py

# BCFW factorization check (< 2 seconds)
python3 computations/verify_factorization.py
```

---

## Key results

### Part I (n=5 base case)
- 15-dimensional ansatz (planar variables) or 45-dimensional (all pairs)
- 5 cyclic 1-zero loci -> rank 44/45, nullspace dimension 1
- Unique solution = A_5^tree. QED.
- Only 3 of 5 loci are independent (legs 1,2,3 suffice)

### Part II (general n)
- Homogeneity Lemma: BCFW weight components inherit 1-zero vanishing
- Dimensional Reduction: restriction to Z_k maps B_n into B_{n-1}
- Induction: base case n=5 + BCFW weight decomposition + reduction

### n=6 verification
- 84-dimensional ansatz (simple poles): rank 83/84, nullspace dim 1
- 165-dimensional ansatz (with higher poles): rank 164/165, nullspace dim 1
- Both give A_6^tree as unique solution
- 5 of 6 loci needed (corrects floor(n/2) conjecture)

### Remaining gap
One flagged gap in the induction: the BCFW residue factorization for
general B in B_n (Step 4 of Theorem 5.6). At n=5 and n=6 this is
bypassed by direct computation. For n>=7, either a direct computation
or a self-contained factorization proof is needed.

---

## Status

- Part I (n=5 base case): **complete** — fully rigorous
- Part II (general n): **complete with one flagged gap** — factorization step
- n=6 verification: **complete** — independent numerical confirmation
- All computational scripts: **complete and passing**

---

## Author

Jonathan Valenzuela, 2026
